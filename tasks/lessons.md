# Lessons learned

## Power / Times radical canonicalisation are coupled (2026-05-24)

### `Power[N, p/q]` splitting and Times "generalised radical fusion" are inverses

Adding a `Power[Integer, Rational]` split that produces
`2^(1/3) * 3^(2/3)` from `18^(1/3)` triggered immediate infinite
recursion: the Times canonicalizer in `src/times.c` has a
"Generalized radical fusion" rule that re-combines
`Power[a, e_i] * Power[b, e_j]` into a single `Power` whenever one
exponent is an integer multiple of the other (k = ±1, ±2, ...). My
split's `Times[Power[2, 1/3], Power[3, 2/3]]` matched `k = 2` and
fused back to `Power[18, 1/3]`, which re-entered my split.

Lesson: when changing the canonical form on the Power side, you must
*also* relax the symmetric simplification on the Times side, or the
two rules pump against each other and recursion limit hits within a
handful of evaluator passes.

The fix in `times.c` was to gate the `|k| > 1` branch on
`gcd(base_i, base_j) > 1`. Same-prime cancellations
(`2^(1/3) * 8^(2/3) -> 4 * 2^(1/3)`,
`12^(1/3) * 2^(-2/3) -> 3^(1/3)`) still fuse because the GCD is
non-trivial; coprime-prime pairs with `|k| > 1` stay split (the new
canonical form).

## zero_test / PossibleZeroQ (2026-05-24)

### Mathilda's `numericalize` keeps `Rational[Real, Real]` un-collapsed

`N[1/10^30]` returns `Rational[1.0, 1.0e+30]` (a function with two Real
args) rather than a single `EXPR_REAL`. Helpers that assume "if it's
numeric, `is_rational(e, &n, &d)` extracts int64 components" silently
fail and treat the value as non-numeric. Always also accept the
`Rational[any-numeric, any-numeric]` shape and divide manually.

### Stage-3 sampling: `evaluate(sub)` collapses cancellation context

For Schwartz–Zippel substitution, *don't* call `evaluate` on the
substituted expression — the evaluator eagerly numericalizes (e.g.
`Sin[Complex[19, -16]]^2 + Cos[...]^2 - 1` → `Complex[0.078, -0.24]`),
discarding the `Plus` structure that the cancellation-aware threshold
needs. Pass the substituted-but-unevaluated form straight to
`decide_numeric`; its first numericalize rung will do the evaluation
while `magnitude_scale_at` still sees the original operand magnitudes.

### MPFR doesn't propagate through every numeric path

Mathilda's `Sin` / `Cos` of `Complex[Real, Real]` evaluates at machine
precision regardless of the requested MPFR precision — the result is a
double-precision number padded out to the printed digit count. A
naive ladder that tightens its threshold by 2^(-p/2) per rung will
spuriously declare "non-zero" because the residual never shrinks.
Detect this by checking whether the magnitude actually drops between
rungs (`m < prev_mag * 0.5`); if not, accept the prior verdict.

### `add_test(...)` in `tests/CMakeLists.txt` is a no-op

`tests/CMakeLists.txt` never calls `enable_testing()`, so `add_test`
lines are silently ignored. The project's test convention is to invoke
each `*_tests` binary directly, not via `ctest`.

## simp_factorial (2026-05-07)

### Mathilda's Factor changes behaviour inside Simplify

`builtin_factor` at `src/facpoly.c:2896` checks
`bool inside_simplify = (factor_memo_top() != NULL);` and uses a
**combined** num/den variable list (vs the **separate** lists used
outside). The combined-scope behaviour can refuse to factor a
denominator like `Factorial[n] + n*Factorial[n]` -> `Factorial[n]*(n+1)`
because the polynomial-in-`n` viewer treats `Factorial[n]` as a
non-factorable coefficient.

Workaround: `factor_memo_push(NULL)` around the Factor call to
force-disable the inside-Simplify branch for that one invocation,
then `factor_memo_pop()`. Direct user `Factor[a]` works because no
memo is active.

This is now documented in the `simp_factorial` source comments.

### Mathilda does NOT auto-coalesce `Power[a,-1] * Power[b,-1]`

Mathematica's evaluator combines `Times[Power[a, -1], Power[b, -1]]`
into `Power[Times[a, b], -1]`; Mathilda's evaluator does not. The
un-coalesced form scores higher under SimplifyCount than the
coalesced form (count 12 vs 9 on `1/(n*(n-1))`), so a factorial
rewrite that lands at a Times-of-inverses can lose the round-loop
tiebreak even though it is the canonical answer.

`simp_fact_combine_inverses` handles this manually: at every Times
node, partition children into "no `-1` power" vs "carries exponent
`-1`", coalesce the latter into a single `Power[Times[...], -1]`.
Conservative -- only exponent `-1` is collapsed, not arbitrary
negative exponents.

### `Together` expands polynomial denominators

Mathilda's `Together` returns `n/(Factorial[n] + n*Factorial[n])`, NOT
`n/(Factorial[n]*(n+1))`. The factored form has to be recovered via
a follow-up `Factor` (with the memo workaround above).

### simp_classify must route factorial inputs to the general pipeline

The rational / polynomial pipelines (`simp_pipeline_rational`,
`simp_pipeline_polynomial`) don't seed `FactorialRules`, so factorial
inputs that classified as `SHAPE_RATIONAL` would silently miss the
factorial rewrite. Added `if (contains_factorial(e)) return
SIMP_SHAPE_GENERAL;` to `simp_classify`.

### Force-take vs SimplifyCount

A factorial-free form often scores **higher** under SimplifyCount than
the factorial-bearing input (the input has fewer leaves). Without a
force-take, the round loop reverts to the input. The fix mirrors
LogExpRules / AssumptionRules: when the candidate strictly reduces the
factorial-atom count, force-take it as the new best regardless of
SimplifyCount tiebreak.

### Forward declarations across the simp.c monolith

`simp.c` is ~10K lines with a single-pass C99 compile, so a helper
defined after `simp_search` cannot be called from `simp_search` /
`transform_can_fire` without a forward decl. Putting the forward
decls right above `SIMP_TRANSFORMS[]` (the earliest place they're
needed) keeps the cluster discoverable.

## ASSERT() elides its argument under Release / NDEBUG

`tests/test_utils.h` has `#define ASSERT(cond) assert(cond)`. The
project's CMake `Release` build (the default) defines `-DNDEBUG`, so
`assert(cond)` becomes `((void)0)` — and the **expression `cond` is
not evaluated**.

Placing a side-effecting call inside `ASSERT` will silently skip the
call in Release builds. Symptom seen 2026-05-07: `test_qaupoly_*`
tests crashed in cleanup because

```c
ASSERT(qaupoly_divrem(a, b, &q, &r));   // evaluated only in Debug
qaupoly_free(q);                         // q is uninitialized in Release
```

The fix is the standard pattern used by `test_zupoly.c`:

```c
bool ok = qaupoly_divrem(a, b, &q, &r);
ASSERT(ok);
qaupoly_free(q);
```

Or override `ASSERT` at the top of the test file to always evaluate
its argument (the trick `test_zupoly.c` uses):

```c
#undef ASSERT
#define ASSERT(cond) do { if (!(cond)) { fprintf(stderr, "FAIL: %s\n", #cond); exit(1); } } while (0)
```

When auditing other test files, look for `ASSERT(funcname(...))`
where `funcname` allocates / mutates / writes via output-pointer.
Pure predicate checks (`ASSERT(qa_eq(a, b))`, `ASSERT(p->deg == 1)`)
are safe to elide in Release because they only weaken the test, not
break setup.

## Numeric helpers must accept BigInt-backed Rationals (2026-05-09)

`add_numbers` (plus.c) and `multiply_numbers` (times.c) had a fast
path for `Rational[Integer, Integer]` but no fallback when the
numerator or denominator was already an `EXPR_BIGINT`.  The helpers
would then hit `return NULL`, and the callers in `builtin_plus` /
`builtin_times` blindly fed that into `is_overflow()` and crashed.

This was latent for years and only surfaced when the rational-
function integration corpus exercised resultant computations
whose intermediate coefficients overflowed 64 bits.

Lesson: every numeric helper that branches on operand type must
treat `Rational[<Integer-or-BigInt>, <Integer-or-BigInt>]` as a
single rational case.  The signature `is_rational(e, &n, &d)` with
int64-out-pointers is too narrow on its own — pair it with a
fall-through GMP path that recognises `Rational[BigInt, ...]` (and
the matching defensive NULL handling at every call site that does
`x = helper(...); is_overflow(x)`).

## Subresultant PRS over Q(α): Power[α, k/m] vs Times[α^q, Sqrt[α]] don't combine via Plus (2026-05-09)

When implementing Bronstein's subresultant PRS for Resultant, naive
`pseudo_rem` over a Q(α)[t] coefficient ring (with α a radical, e.g.
Sqrt[3]) blows up geometrically.  The chain element coefficients
should stay bounded — Bronstein's β-scaling keeps the chain in
D[x] — but our system's algebraic-number canonicalisation has a
subtle issue: `Sqrt[3]^3` auto-simplifies to `Power[3, 3/2]`, NOT
to `Times[3, Sqrt[3]]`, and `Plus` treats those as different terms
because they're structurally distinct.  So the same algebraic
value accumulates in many forms each chain step, doubling the
expression size every iteration.

`Together` recognizes the equivalence (it does its work in a
canonical-fraction representation), but is too slow to call
per-coefficient or even per-chain-step on big polynomial inputs.

Lesson: when an algorithm's correctness depends on
"algebraically-equal sub-expressions canonicalize identically,"
verify that condition holds in the host CAS for the specific
coefficient ring — and gate the algorithm out for rings where it
fails.  A conservative `subres_has_algebraic` predicate (anything
of shape `Power[X, Rational[a, b>1]]`) is enough to keep the new
fast path correct in practice while routing alg-number cases
through the existing matrix path.

The fundamental fix would be a proper Q(α) substrate (qaupoly
or similar) where coefficients are reduced modulo α's minimal
polynomial after every operation, but that's a much larger change.

## Recursion on tree size, not on "variables stripped" (2026-05-09)

`is_zero_poly` recurses by stripping one variable per descent
via `CoefficientList(expanded, vars[0])`.  The induction is
"polynomials in fewer variables are simpler," which holds only
when `vars[0]` is a real polynomial variable.  When
`collect_variables` returns an algebraic constant like `Sqrt[5]`
as `vars[0]` and the polynomial mixes several radicals,
`CoefficientList` does not actually strip anything and the
recursive call sees the same expression — unbounded recursion,
EXC_BAD_ACCESS at the next deep call site.

Lesson: any recursive simplification whose termination relies on
"each call sees a smaller subproblem" needs an explicit depth
bound when the smallness predicate (here: "fewer
non-numeric leaves") can be defeated by the pre-processing.
Pick a bound well above any genuine tree depth and bail out
conservatively (return the safe answer for the caller) when it
is exhausted.

## `exact_poly_div` field-vs-ring soundness (2026-05-10)

The `var_count == 0` base case in `exact_poly_div` (poly.c) used
to fall through to a symbolic `Times[A, Power[B, -1]]` for any
non-bigint pair, on the assumption that the coefficient ring is a
field. That assumption holds for Q and Q[i] but breaks the moment
a non-rational atom appears (Sqrt[2], Sqrt[3], ...): in
`Q[Sqrt[2], Sqrt[3], ...]` the divisor doesn't actually divide the
dividend, and the symbolic Times propagates a `Power[Plus, -1]` up
into intermediate polynomials. The downstream `PolynomialGCD` call
then runs multivariate Euclid on a rational input — that's the
case-13 Together hang.

Lesson: if a function is named "exact" division it must return
NULL on non-exactness. Returning a symbolic `Times[A, B^{-1}]` as
a "fallback" hides correctness bugs from callers and only ever
comes back as a hang or a wrong answer. Restrict the symbolic
fallback to the strict cases where it's actually exact (operands
in Q or Q[i]); for everything else, return NULL and let callers
choose what to do.

## Plus auto-distribute `Times[-1, Plus[…]]` (2026-05-10)

Mathilda's Plus auto-eval groups by `(coeff, base)`. When a Plus
arg is `Times[-1, Plus[A, B]]`, get_coeff_base returns
`(-1, Plus[A, B])` — but the OTHER args have bases `A`, `B`
(distinct from `Plus[A, B]`), so no cancellation fires. The
canonical Mathematica behaviour is to distribute the leading -1
into the inner Plus before grouping, so `a + b - (a + b)` reduces
to 0.

Lesson: when adding distribution rules to Plus, gate them on the
literal `-1` coefficient — distributing arbitrary `c·Plus[…]` into
the outer Plus would expand harmless products like `2 (a + b)`,
diverging from MMA's behaviour. The `-1` case is the cancellation-
enabling step; other coefficients stay as Times factors.

## exact_poly_div NULL contract: callers must check (2026-05-11)

`exact_poly_div` (`src/poly.c`) was tightened on 2026-05-10 to return
NULL when its operands are not in a field (anything beyond `Q` /
`Q[i]`, e.g. `Q[Sqrt[a]]`).  Several callers were not updated to
check for NULL; one in `heuristic_factor_impl` (`src/facpoly.c`)
fed the NULL straight into a recursive `heuristic_factor` call,
which dereferenced `P->type` and crashed `Simplify` on inputs like
`Simplify[(Sqrt[a] - Sqrt[a] x)/(2 Sqrt[a] + 2 Sqrt[a] x)]`.

Lesson: when a helper's contract changes from "always returns a
result" to "may return NULL on non-applicability", every existing
call site needs an audit.  `grep -n exact_poly_div src/` should be a
required step after that kind of soundness tightening.  Adding a
defensive `if (!P) return NULL;` at the top of
`heuristic_factor` is good belt-and-braces for any future callers
that forget.

## intrat: per-summand `c · piece` accumulator needs explicit Expand (2026-05-11)

After `intrat_integrate_summands` builds each Apart piece's integral
as `c_k · piece_int_k` (with `c_k` the constant from
`extractConstants`), the result accumulates as
`Plus[Log[x], Times[2, Plus[1/4 ArcTan[…], -1/8 Log[…]]], …]`.
Mathilda's `Times` has `ATTR_FLAT | ATTR_ORDERLESS | …` but NO
auto-distribution over `Plus`, so the literal `Times[2, Plus[…]]`
survives all the way to print.  Mathematica's
`IntegrateRational.m:99` calls `Collect[intlog // Expand, …, simproot]`
exactly to flatten this — we need the analogous Expand pass.

Lesson: when porting a Mathematica pipeline, every `// Expand` /
`Collect` / `// Distribute` in the source is load-bearing.  Don't
assume Mathilda's evaluator does the equivalent — it doesn't.  The
two cheap post-passes (`expr_expand` to distribute Times-over-Plus,
plus a Log-arg constant-stripper for `Log[c · p] -> Log[p]`) want to
run both BEFORE and AFTER `intrat_log_to_arctanh` so the log-pairing
rule sees fully-distributed sums.

## intrat: sign-pos-assumption can't decide `Sqrt[5] − 5 < 0` (2026-05-11)

`intrat_sign_pos_assumption` treats free symbols as positive reals
and walks Plus arg-by-arg; for `Plus[−5, 2 Sqrt[5]]` (one negative
constant, one positive radical) it bails to sign-unknown.  This
is correct for parametric inputs but blocks the closed-form
palindromic-quartic dispatch in `logtoreal_dispatch` and
`expand_palindromic_quartic_real`, where we routinely need to
decide inequalities like `((1 + Sqrt[5])/2)^2 − 4 < 0`.

Fix: a thin `intrat_numeric_sign(e)` helper that evaluates `N[e]`
and returns ±1 when the result is a definite-sign Real / Integer
(`|v| > 1e-12` deadzone for round-off).  Use it as a fallback
**only** when `intrat_sign_pos_assumption` returns 0; it preserves
parametric correctness (`N[Plus[a, …]]` with a symbolic doesn't
reduce to a number) while resolving radical-only inequalities.

Lesson: numeric evaluation is a legitimate fallback for sign
decisions when symbolic positive-walks bail.  The two-tier pattern
(symbolic-then-numeric) is the right structure: symbolic for
parametric soundness, numeric for closed-radical decidability.

## intrat: scaled-palindromic LRT Q breaks LogToReal substitution chain (2026-05-11)

For `1/(x^5 + 1)`, the LRT producer's Q-in-t is the *scaled*
palindromic `625 t^4 + 125 t^3 + 25 t^2 + 5 t + 1` (palindromic
under `u = 5 t`, not under `u = t`).  A first cut at the
scaled-palindromic case added a detector + factor-via-`u = r t`
substitution to `logtoreal_dispatch`; the math is right, but the
resulting `logtoreal_quadratic` calls hand
`subst_t(s, t, u_root + I v_root)` nested-radical complex values to
`split_re_im` / `LogToAtan`, and LogToAtan's polynomial-GCD
machinery wedges for 2+ minutes.

Workaround taken: restrict the `logtoreal_dispatch` palindromic
branch to the pure-palindromic case `r = 1` (cheap radicals after
substitution), AND in parallel add an `expand_palindromic_quartic_real`
hook inside `intrat_naive_log_part` that builds the real form
directly from `(a(α) / d'(α)) · Log[x − α] + (conj …)` for each
conjugate root pair — bypassing LogToReal / LogToAtan entirely.

Lesson: when porting a Mathematica routine that delegates to
`Solve` / `ToRadicals` / general factor-over-extension, the C port's
hand-coded substitution chain (LogToReal → split_re_im → LogToAtan)
has very different cost characteristics than Mathematica's.  Inputs
where the LRT Q has clean radicals can route through LogToReal;
inputs where the LRT Q's radicals get cubed / squared by the
substitution should be routed around it — handle palindromic /
cyclotomic structure at the NaiveLogPart layer (closer to the
integrand's d(x)) where the radicals stay shallow.

## intrat: narrow per-degree branches in logtoreal_dispatch leak RootSum at adjacent n (2026-05-11)

The first nth-root fix only handled `c_3 t^3 + c_0` and (via the
older Sophie-Germain shortcut) `c_4 t^4 + c_0` with `q < 0`.
User immediately surfaced nearby failures: `1/(b ± a x^n)` for
n ∈ {4-, 5-, 6±, 8+, 9+} all still leaked `RootSum`, because the
LRT-Q polynomials for those inputs were also sparse nth-root
forms with one degree-specific branch each that hadn't been
written.

**Why:** when a CAS-layer fix is required for a family
`P(c_1, …, c_k, n)` (here: `1/(b ± a x^n)` for varying n), the
fix must be parametric over the family member, not a single
hand-coded case.  Adjacent family members will be the next thing
the user tries.  Reinforces the existing "general algorithms"
memory: heuristics keyed on the literal failing input are not
fixes.

**How to apply:** when implementing a degree-n branch in a
dispatch, ask "what's special about this degree?".  If the answer
is "I'm only doing this degree because that's where the user's
report landed", that's the wrong unit of work — find the structure
that lets the routine extend to all n in the same closure step.
For the LRT nth-root case the structure is the standard
cyclotomic decomposition over R: enumerate angles, pair conjugates,
build quadratics with `Cos[k π / n]` coefficients.  One helper
(`logtoreal_nthroot_sparse`) handles all n ≥ 3 in ~140 lines.
Adjacent family members (different n, different sign of `−c_0/c_n`)
fall out for free.

## CRC integral rules: branch-correct forms (2026-05-15)

Trigger: 11 DIFF-NONZERO regressions in the CRC corpus, several of
which were "correct on the principal real branch but wrong as
expressions".  Classical CRC-table forms like
`∫Sqrt[1 + Cos[a x]] dx = (2 Sqrt[2]/a) Sin[a x/2]` reduce the
integrand's radical via `1 + Cos[a x] = 2 Cos²[a x/2]`, but the
resulting closed form's derivative is `Sqrt[2] Cos[a x/2]`, which
matches `|Cos[a x/2]|` (the actual `Sqrt[…]` value) only when
`Cos[a x/2] ≥ 0`.

**Why:** the `D[r, x] - integrand` corpus check evaluates the
diff under ordinary `Simplify`/numeric eval, which respects branch
choices — it does NOT do PowerExpand.  So a rule that
"folds the radical out" silently picks one branch and fails outside
it.  Symbolic-only verification hides this; numeric sampling
across the natural domain exposes it.

**How to apply:**

1. When writing a CRC-style rule that handles `Sqrt[trig identity]`,
   keep the integrand's radical literally in the antiderivative
   (e.g. `(2/a) Tan[a x/2] Sqrt[1 + Cos[a x]]` instead of `Sqrt[2]
   Sin[a x/2]`).  The two forms are equal on the principal branch,
   but only the first has a derivative that the simplifier can
   show equal to the integrand without PowerExpand.
2. Same trick for algebraic integrands: for `∫ Sqrt[(a+x)/(a-x)] dx`,
   the form `(x − a) Sqrt[(a+x)/(a-x)] + a ArcSin[x/a]` keeps the
   integrand intact in the antiderivative.
3. When a corpus runner uses numeric sampling, sample only inside
   the integrand's natural real domain — evaluate `N[integrand]` at
   each candidate point first, and skip the ones that come back
   complex / infinite / indeterminate.  Otherwise you get false
   positives at branch-cut sample points (e.g. `Sqrt[(1+x)/(1-x)]`
   at `x = 1.7`) that are not regressions, just samples outside the
   rule's natural domain.

## Pattern guards: missing FreeQ is silent corruption (2026-05-15)

Trigger: `1/x_ Power[b_. x_^2 + a_, -1/2]` was a top-of-file CRC
rule with NO `FreeQ` guard.  It matched
`1/(x Sqrt[a + b x + c x²])` with `b → c`, `a → (a + b x)`,
producing `-ArcTanh[Sqrt[a + b x + c x²]/Sqrt[a + b x]]/Sqrt[a + b x]`
— wrong on every quadratic-surd integrand, shadowing the correct
Formula 246 below it.

**Why:** unguarded `a_` matches *any* subexpression, including
ones that contain `x_`.  The default-coefficient `b_.` makes the
shadow even worse: it lets the pattern accept a Plus with the
"wrong" number of terms.  Subsequent rules never see the input.

**How to apply:** every `_.`/`_` pattern variable that names a
"constant" (in the sense of the surrounding formula's preconditions)
needs an explicit `FreeQ[{vars}, x]` guard — even rules at the very
top of the table, especially those that match an aggressive
super-pattern (Power with negative exponents, Plus with optional
coefficients, etc.).  When a CRC corpus run shows a closed-form
result with sub-expressions that *contain the integration variable*
in positions meant to be "constants", look for an unguarded
top-level pattern intercepting the dispatch.

## Pattern matcher: trust the repro over the task note (2026-05-16)

Trigger: `tasks/crc_corpus_2026-05-15.md` "Out-of-scope findings"
described an "underlying matcher gap on a + b x + c x² Plus patterns"
affecting ~80 CRC formulas.  When a user asked about the issue, I
quoted the note as authoritative.  The user immediately
disproved it with `MatchQ[2 + 3 x + 5 x^2, a_ + b_. x_ + c_. x_^2] ==
True`.  The actual gap was in the §3.5 DownValue dispatch filter:
held LHSs sit in parser-shape (`Power[Times[x_, Sqrt[...]], -1]`)
while runtime inputs arrive in evaluated-shape (`Times[Power[x,-1],
Power[Sqrt[...],-1]]`), and the filter's pointer-equal head compare
rejects the rule before `match()` runs.  Fixed by canonicalizing the
LHS at insertion time (`pattern_canonicalize` in `src/symtab.c`).

**Why:** task notes capture what the author *thought* the cause was
at the time, not necessarily the right cause.  "Out-of-scope" findings
are by definition unverified — the author didn't trace them to ground.
Quoting them carries that uncertainty forward, and a wrong root-cause
attribution wastes time on the wrong fix.

**How to apply:** when a task note names a root cause and the user
asks about it, run the smallest possible repro first — `MatchQ`,
`FullForm`, a one-line REPL test — and let the actual behaviour
correct or confirm the note.  Then update the note (or the
changelog) with the verified cause.  Especially for matcher /
dispatch / evaluator bugs, where the surface symptom and the
underlying mechanism often have nothing to do with each other.

## Memory ownership patterns (2026-05-16)

### evaluate() and expr_expand() do NOT consume their input

Both `evaluate(e)` and `expr_expand(e)` start with an internal
`expr_copy(e)` (refcount++) and return a fresh result; the caller
still owns `e` and must free it.  The leak-bait pattern is:

```c
Expr* x = evaluate(internal_times(...));      // temp leaks
Expr* y = expr_expand(internal_power(...));   // temp leaks
foo->slot = evaluate(expr_new_function(...)); // temp leaks
```

The fresh `internal_X(...)` / `expr_new_function(...)` result has
nowhere to go — no variable, no free.  Fix either by introducing a
local, calling, then freeing, or by using a wrapper that consumes:
`eval_and_free(e)` from `src/eval.h` does this for `evaluate`; for
`expr_expand`, `intrat.c` has a `expand_and_free` helper.

**Why:** Found while fixing valgrind leaks in Integrate`RischNorman
/ BronsteinRational unit tests (May 2026).  Multiple modules
(`intrischnorman.c`, `intrat.c`, `symtab.c`, `deriv.c`) had this
pattern in wrappers like `eval_expand`/`eval_cancel`/etc. that
themselves called `expr_copy(f)` internally — meaning callers also
had to free `f` separately, which they often forgot for fresh
temps.  The eval_* wrappers in `intrischnorman.c` were converted to
take ownership of their argument; that single contract change
eliminated 70+ leak sites.

**How to apply:** when reviewing a new helper that wraps
`evaluate` or `expr_expand`, decide explicitly whether it consumes
or borrows.  If it borrows, every call site passing a fresh temp is
a latent leak.  Prefer the consuming contract (less code at call
sites, no leak risk).  When a callee both copies internally AND the
caller never references the input again, drop the redundant copy.

### expr_new_function memcpys args but leaves the array to caller

`expr_new_function(head, args, count)` allocates its own backing
store and `memcpy`s the `args[]` slot pointers in.  The new
function "owns" the referenced Expr*s (they're decremented when it
frees), but the **`args` malloc itself** is still the caller's —
must be `free`'d if heap-allocated.  Stack/compound literal `args`
arrays are fine.

**Why:** Found two leak sites with this exact bug:
- `src/context.c:context_path_as_list` — `malloc`'d args, passed to
  `expr_new_function`, never `free`'d.  Leaked once per process.
- `src/intrat.c:intrat_apart_list` — same shape, leaked once per
  call.

**How to apply:** any `Expr** args = malloc(...); ... expr_new_function(h, args, n)`
must have a matching `free(args)` after the call (NOT `expr_free` —
that would double-decrement the slots).

### upoly_div_rem_mod overwrites *out_r unconditionally

`upoly_div_rem_mod(a, b, mod, &q, &r)` writes to `*r` without
freeing whatever it pointed to.  Most call sites pass an
uninitialized local or freshly-freed slot, so this works.  But one
site in `cz_ddf` passed `&x_pow_p` while `x_pow_p` still held a
live UPoly — leaking it.

**Why:** Out-param contracts are easy to get wrong when the slot
is shared across loop iterations.  Found via valgrind during the
May 2026 leak hunt.

**How to apply:** when a helper assigns through an out-param
pointer, ensure the slot is empty (NULL or just-freed) before the
call; or use a fresh local and assign after.

### Module does not substitute locals into HoldAll bodies (Table, etc.)

`Module[{lu = ...}, Table[lu[[i, j]], {i, n}, {j, n}]]` leaves the
inner `lu` references as the literal symbol `lu` — Mathilda's
`Module` does not propagate its renaming into Hold-* arguments of
nested calls.  `Block[{lu = ...}, Table[lu[[i, j]], ...]]` works
correctly because `Block` uses dynamic scoping (the symbol's
existing value is temporarily replaced).

**Why:** Found writing the LUDecomposition unit tests
(2026-05-22); `Module`'s scoping rule did not reach into
`Table[..., {i, n}]`, so the `lu` in the body printed as a free
symbol.

**How to apply:** for tests / scripts that want to bind local
matrix data and then iterate over it with `Table`, `Sum`, `Map`,
etc., use `Block` rather than `Module`.

### Iteration variables in Table can collide with symbolic matrix entries

`Table[If[i > j, lu[[i, j]], 0], {i, n}, {j, n}]` applied to a
matrix `{{a, b, c}, {d, e, f}, {g, h, i}}` silently corrupts the
result: the `i` from the matrix takes precedence in the `If`
condition (or in the indexing expression), giving wrong output
without any error.

**Why:** Found writing the LUDecomposition symbolic test for
`{{a, b, c}, {d, e, f}, {g, h, i}}` — the iteration variable `i`
clashed with the matrix entry `i`, and the residual identity
check failed with `{0, 0, -3 + i}` instead of all zeros.

**How to apply:** when a `Table` (or any iterator) is going to
operate over general symbolic data, use iteration variables that
are unlikely to appear in the data — `ii`, `jj`, `kk`, or
`Module[{i}, Table[..., {i, n}]]`.

## Integration by parts for unknown functions (Roach §1.7) — 2026-06-06

**Lesson 1 — zero-test a rational difference with `Cancel[Together[Expand[...]]]`,
not Expand or Together alone.** The unknown-function integrator's linear check
and residual (`newI`) must collapse mathematically-zero expressions to literal
`0`. `Expand` alone fails on different-denominator fractions (e.g.
`1/(1+(g'/f)^2)` vs `1/(f^2+g'^2)` from an ArcTan derivative); `Together`/`Cancel`
alone fail to combine syntactically-distinct-but-equal *products* like
`g(1+x^2)` vs `(g + g x^2)`. Only all three passes together are reliable. A bug
where `(1+x^2)` coefficients failed while `x^2` worked traced exactly to a
`canon` that omitted `Expand`.

**Why:** `Together`/`Cancel` operate on fraction structure (common denominator,
cancel common factors) but do not distribute a sum factor over a product;
`Expand` does the distribution but cannot merge fractions over a common
denominator. They are complementary, not redundant.

**Lesson 2 — never split an integrable sum into term-by-term integrals.** A
residual like `f'[x]g'[x] + f[x]g''[x]` is the exact derivative `(f g')'` and
integrates cleanly *as a whole*, but its individual terms (`f'g'`) have no
closed form and send integration-by-parts into an infinite cycle
(`∫f'g' = fg' - ∫fg''`, `∫fg'' = fg' - ∫f'g'`, …). Hand the whole residual back
to the integrator; only factor out genuine `x`-free constants from a *single*
term for cosmetics. Splitting caused a segfault (a NULL function-arg built from
the runaway recursion).

**How to apply:** any by-parts / linearity engine that recurses through a global
`Integrate` must (a) keep sums intact across the hand-off, and (b) carry a
canonical-form cycle guard (stack of in-flight integrands compared with
`expr_eq` after `Cancel[Together[Expand[...]]]`) so genuinely non-elementary
inputs terminate unevaluated instead of looping.
