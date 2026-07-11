# Lessons learned

## No NIntegrate crosscheck inside Integrate — verify correct-by-construction (2026-07-08)

The definite-integral methods in `Integrate` must NEVER validate a symbolic
result against `NIntegrate` (or any numeric quadrature). The project philosophy
(see the residue method `integrate_residue.c` and the cascade-ordering lesson) is
**correct-by-construction**: deterministic methods gated by symbolic
convergence/assumption conditions, with symbolic self-verification only
(`PossibleZeroQ`, `Simplify`, exact base values).

For the DiffUnderInt (Feynman) method specifically, this is both a rule and a
gift: the Conrad §12 conditional-convergence trap (differentiating `∫sin(tx)/x`
into the divergent `∫cos(tx)`) is caught automatically because the inner
`Integrate[∂_p f, {x,a,b}]` fails to close / returns divergent under the engine's
own gates — so we skip that parameter. Verification = symbolic derivative check
`PossibleZeroQ[D[I,p] − J]` + an EXACT base value `I(p0)` (zero-integrand or the
engine's exact `Integrate` of `f|_{p→p0}`), never a numeric compare.

Test-side numeric comparison (`N[result - expected]` in `test_*.c`) is fine — the
prohibition is on numerics *inside* the `Integrate` code path.

## Carving a C file into regions: grep ALL return types, not just the obvious one (2026-06-07)

When splitting `arithmetic.c` into `numbertheory.c`, I mapped function
boundaries with `grep '^(static )?Expr\*|^void .*_init|^static Expr'`. That
pattern misses functions returning `int`/`bool`/other types. Two such groups
sat *inside* the span I was moving and silently broke the cut:
- PowerMod's `static int` modular-root helpers (modroot_brute, tonelli_shanks,
  hensel_lift, modular_root) — needed to move WITH PowerMod.
- The `bool`/`int` numeric predicates (is_infinity_sym, expr_numeric_sign,
  is_neg_infinity_form) — core helpers used by plus/times/power that had to
  STAY, even though they were physically interleaved among NT builtins.

Rule: before choosing line ranges, enumerate every top-level definition with a
return-type-agnostic grep like
`grep -nE '^[A-Za-z_].*\b[a-z_]+\s*\(' file.c` (or list all of `Expr*`, `int`,
`bool`, `void`, `static ...`). Then classify each by concern, not by file
position — interleaved code means contiguous line ranges rarely equal one
concern. The main build can link fine while the test build (fixed COMMON_SRC)
exposes the misplacement, so always build BOTH after a move.

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

---

## Polynomial-ideal frameworks: fix the *boundary*, don't change the generator (2026-06-09)

**Context.** `Integrate[Sqrt[Cot[x]], x]` threw `Power::infy` (1/0) where the
symmetric `Sqrt[Tan[x]]` worked. Root cause in `src/simp/trigrat.c`: a radical
`Sqrt[g]` is carried as `l` with `l^2 = g`; for `Cot = Cos/Sin = c·s^(-1)` the
radicand is rational with the *odd* generator `s` in its denominator, so
reducing `l^2` injects `s^(-1)`, and the conjugate `den|_{s->0}` evaluates
`Power[0,-1]`. (`Tan = s/c` escaped: its inverse generator is the *even* `c`,
never substitute-zeroed.)

**Correction from user.** My first fix normalised the radicand globally via
`Sqrt[N/D] = Sqrt[N D]/D` (so `l^2 = N·D` is always polynomial). It removed the
crash but made the case **~20× slower** — `l^2 = c·s` (degree 2) instead of the
natural `c/s` raises every downstream `Together`/`Cancel` degree. User caught the
hang.

**Lesson.** In a polynomial-ideal reduction (Gröbner-like normal forms,
conjugate rationalisation), changing the generator's defining relation to dodge a
degenerate boundary case usually inflates degree and tanks performance across the
*whole* pipeline. Prefer a **local, conditionally-triggered repair at the exact
boundary** that misbehaves. Here: keep the natural low-degree `l^2 = g`, and just
before the conjugate substitute-zero, clear inverse powers of the generator being
zeroed (`tr_has_neg_sgen_power` / `tr_clear_neg_sgen`: recombine → re-split →
re-reduce). Gate the trigger narrowly (odd-generator inverse powers only) so the
already-working path (`Tan`, `c^(-1)`) stays byte-for-byte unchanged and pays
nothing.

**How to apply.** When a symmetric pair (A works, B crashes) diverges inside a
canonical-form engine, find the *single* operation that assumes a precondition B
violates (here: "den is polynomial in the var being cleared"), and restore that
precondition in place — don't re-architect the representation both cases share.
Verify the fast case is untouched (diff its output/timing) and the slow case
matches the fast one's cost.

---

## DerivativeDivides hang: a loop guard must fix *termination* AND *cost* (2026-06-09)

**Symptom.** `Integrate[x Sin[x^2], x]` hung. The reduced sub-integral
`Integrate[Sin[u]/2, u]` re-enters the *full* derivative-divides stage; each
level mints a fresh substitution variable, so overlapping branches regenerate
the same integrand and fan out exponentially with the expensive Eliminate/Solve
search at every node.

**Correction from user (twice).** (1) My first instinct — gate the Eliminate
search off on recursive calls — was rejected: "we still should try
derivative-divides on recursive calls, but guard against infinite loops by
checking the current integrand against previous ones." So I built an **integrand
memo** (canonicalise by renaming the integration var to a fixed sentinel so
gensym'd duplicates compare equal; short-circuit anything seen in this descent).
(2) The memo *terminated* but left a 7–18 s near-hang — it bounds the *count* of
nodes but not the ~0.1–1 s cost of Eliminate at each. I surfaced the data and
the user then chose memo **+** restricting the heavyweight Eliminate/Solve search
to the outermost call (direct derivative-divides still runs recursively).

**Lesson.** "Stop the infinite loop" ≠ "make it usable." A memo/cycle guard on a
recursive symbolic routine restores *termination* but does nothing for a heavy
per-node cost — a bounded-but-exponential-work tree still reads as a hang.
Always **measure wall-clock after the guard**, not just "does it return." When a
correctness guard isn't enough, the lever is usually to confine the *expensive*
strategy to where it pays off (here: the heavyweight Eliminate search only earns
its keep on the original integrand; reduced sub-integrals are finished cheaply by
the rest of the cascade).

**How to apply.** Before declaring a hang fixed: (a) confirm it terminates, then
(b) time the actual user-facing case and compare against a known-fast sibling
(`Integrate[Sqrt[Tan[x]], x]` here). If still slow, separate the two failure
modes — *non-termination* (fix with a memo/cycle guard) vs *too-much-work* (fix
by scoping the costly stage) — and address both. When the user has stated a
preferred mechanism that you find insufficient, implement it faithfully, then
present before/after numbers and let them choose the augmentation rather than
silently overriding.

## Fix the root cause, not a cap (Gamma[201/2] / LogGamma, 2026-06-09)

When an exact path delegates to another builtin (LogGamma's half-integer path
called `Gamma`, which calls `Factorial`), do NOT assume the delegate is correct
across the full input range. `Factorial` of a half-integer built the odd double
factorial and `2^k` denominator in `int64_t` and silently overflowed past
~`37/2`, so `Gamma[201/2]` returned garbage and LogGamma inherited it.

My first instinct was to *cap* LogGamma's exact path to a "safe" magnitude and
fall through above it. The user's correction ("We must fix this bug:
Gamma[201/2]") was right: cap = papering over a real overflow that also breaks
`Gamma`/`Factorial` directly. The fix was to rebuild the coefficient in GMP
(`mpz` + `mpz_pair_to_rational_expr`). **Pattern:** if a "safe range" cap is
hiding a defect in a shared primitive, fix the primitive — a silent-wrong-value
bug in a building block is worse and more widely felt than the symptom you hit.
Verify big exact values with folding identities (functional equation,
reflection) so the test asserts a clean result instead of a giant literal.

## Integrate cascade ordering: deterministic domain-specific methods go FIRST (2026-06-09)

While adding the Jeffrey–Rich Weierstrass integrator I first placed it *after*
`DerivativeDivides` (after the radical-substitution stages, before Risch). The
user corrected: a domain-specific algorithm that is **guaranteed to succeed and
is correct by construction** (no differentiate-back verification needed) should
run *ahead* of the search-and-verify methods (`DerivativeDivides`'
Eliminate/Solve branch search) and ahead of `RischNorman` (whose complex-log
forms for trig rationals are ugly/discontinuous). **Pattern:** order the
`Integrate` cascade by *(a)* confidence — methods that close deterministically
for their domain before heuristic/search methods — then *(b)* output quality
(continuous real forms before complex-log forms). Cheap, precise gating
(`wj_has_kernel_in_denominator`, kernel detection) keeps an early stage from
clobbering integrands other stages handle more cleanly.

## TrigExpand pre-pass for trig/hyperbolic substitution integrators (2026-06-09)

The Weierstrass detector requires every kernel argument to be the bare variable
`x`, so `Cosh[x] Cosh[2 x]` (multiple angle) and `Sin[x + 1]` (sum angle) failed
detection. The user pointed out `TrigExpand` rewrites those into kernels of the
bare `x` (`Cosh[x] Cosh[2 x] // TrigExpand` = `Cosh[x]^3 + Cosh[x] Sinh[x]^2`).
**Pattern:** for any method keyed on "kernel of the bare variable", try the
integrand verbatim first (fast common path), then retry on `TrigExpand[f]` only
on failure. `TrigExpand` leaves single-angle kernels untouched (so the verbatim
path and clean rational denominators are unaffected) and reduces multiple/sum
angles — exactly the normalisation these substitutions need.

## Verify Tan[x/2]-rational antiderivatives with PossibleZeroQ, not Simplify (2026-06-09)

The test predicate `Simplify[D[Integrate[f,x] /. Floor[_]->0, x] - f] === 0`
times out on integrands whose antiderivative carries `Tan[x/2]^3` (e.g.
`1/(1 + Sin[x]^2)`): the residue is a deep nest of half-angle rationals that
defeats `Simplify`. `PossibleZeroQ[...]` (the numeric two-phase sampler) returns
`True` instantly and is the right correctness oracle for these. (Strip the
secular `Floor` term first with `/. Floor[_] -> 0`, since its symbolic `D` is
`Derivative[1][Floor]`, not 0.)

## 2026-06-14 — Rewriting a numerical engine: prefer hybrid over wholesale replacement

Context: fixing 3 NSum deficiencies. Planned to *replace* symbolic Euler–Maclaurin
derivatives with numerical contour derivatives. That broke things the original did
fine (contour is fragile for geometric/oscillatory summands; an over-broad
"black-box → never EM" rule killed valid nested EM and a passing multidim test).

Lessons:
- **Make the new path a SUPERSET, not a replacement.** The robust design kept
  symbolic D as primary (byte-identical to original for simple summands → zero
  regressions) and used the new mechanism only where the old one fails (composite
  summands that balloon). When a rewrite "fixes case A but regresses B", the answer
  is usually a hybrid keyed on the property that distinguishes A from B.
- **Validate a "fix" against the EXISTING test suite early, not just the target
  case.** I confirmed the target (Log WP35) before running `nsum_tests`; the suite
  caught `(-5)^i/i!`, the multidim cases, and a forced-CVZ peaked case I'd have
  missed. Run the affected `*_tests` binary after each behavioural change.
- **Adding per-call work can amplify a PRE-EXISTING leak.** NSum's evaluator leaks
  a GMP rational per summand eval (present on `main`). My far-tail ladder (16
  evals/profile) and an extra oscillatory-probe eval doubled multidim valgrind
  blocks. Fix was to not add evals (skip ladder on monotone heads; read signs from
  existing head terms), not to chase the shared-evaluator leak. Always valgrind a
  representative input against the ORIGINAL binary to separate your delta from
  baseline noise (file-swap `git show HEAD:path`, not `git stash`).
- **MPFR convergence gates must sit ABOVE the roundoff floor.** Setting a DE-quad
  reltol *below* achievable precision means it never trips → refines to a
  catastrophic node count (looked like a hang). Scale reltol to `target-2` digits.

## NIntegrate: a "new method" must be named & listed (2026-06-14)
- When a fix adds a genuinely new quadrature *strategy* (not just a tweak),
  expose it as a named `Method -> "..."` string AND list it in the docstring +
  docs/spec method table — the user explicitly requires this. Wire it in
  `ni_method_from_string` + `ni_method_implemented`, give it a forced path, and
  keep it as an `Automatic` fallback too. (Added `"OscillatorySingularity"`.)
- **A multi-method dispatcher's selection is only as honest as each method's
  conv/err.** The osc finite path returned `conv=true, abserr=0` unconditionally;
  AUTO then preferred its garbage over correct non-converged estimates. Any
  method that competes via `ni_consider` MUST report a real error estimate.
- **An exponential endpoint map x=a±(b−a)e^{−t} samples the singular endpoint by
  rounding**, not just underflow: for end≠0, `1 − tiny == 1` → 1/0. Gate the
  sample on `mapped_abscissa == endpoint`, not on Jacobian underflow alone.
- **The wrong mirror of an endpoint transform converges to a wrong value.** EXP_HI
  on a left-singular integrand samples the (undamped) left singularity at t=0.
  Only transform a *detected-singular* endpoint; never blindly try both.
- **Between-the-zeros marchers need an adaptive step** (sized from the last gap),
  not a fixed half-period — otherwise an accelerating chirp leaps over lobes.

## Post-integration normalisation (intsimp_finalize, 2026-06-15)
- **Result-cleanup passes must respect the verifier's assumptions.** The
  intrat corpus checks `Cancel[Together[Expand[D[result,x]-f]]] == 0` WITHOUT
  assuming free symbols positive. A cleanup that rewrites an x-free constant
  under a positivity assumption (`(1/a)^(1/3) -> a^(-1/3)`, or PowerExpand on a
  negative/complex numeric base like `(-3)^(1/2)`) changes the constant's value
  and the differentiate-back no longer cancels. Scope such rewrites to
  positive bases AND to the result class that actually needs them
  (radical-bearing antiderivatives) — never blanket-apply to all Integrate
  output.
- **Never re-`evaluate()` an integration result that still contains a nested
  `Integrate[...]`.** It re-enters `builtin_integrate` and blows the 1024
  recursion limit → segfault. Guard with a contains-Integrate tree-walk.
- **ArcTan/ArcTanh oddness will undo a hand-rolled sign pull.** Building
  `Times[-1, ArcTan[Times[-1, negPlus]]]` collapses straight back: ArcTan pulls
  the inner `-1` out again. `Expand` the negation into a genuinely positive
  `Plus` before wrapping, so the argument has no leading `-1` to re-pull.
- **The macOS valgrind baseline is 12,800 B / 400 blocks** (dyld/Accelerate);
  grep leak stacks for your own source frames rather than trusting the total.

## Multi-generator radical simplification (simp_radical_rational, 2026-06-15)
- **Reduce before you rationalise.** Combining a multi-radical rational gives a
  big `P/Q` in the generators. Reduce `P` and `Q` modulo the relation ideal
  (`PolynomialRemainder`) FIRST — they shrink dramatically — then rationalise the
  denominator with `PolynomialExtendedGCD`. Rationalising the *un-reduced*
  numerator (`Expand[P*u]`) detonates the multivariate GCD: in dev this left
  three `./Mathilda` processes pegged at 100% CPU for minutes (kill specific PIDs;
  never broad `pkill -f Mathilda` — it also kills the user's REPL). This is the
  "Simplify multi-generator explosion" wall in action.
- **Plain `Cancel` cannot use the generator relation.** After reducing mod
  `{t^3 - (a+b x)}`, `Cancel[P/Q]` still leaves a relation-dependent common factor
  (`s^2 + s t + t^2`) because multivariate GCD treats `s,t,a,b,x` as independent.
  Clear the generators from the denominator via the extended-GCD/norm
  (`PolynomialExtendedGCD[Q, rel, g]`) instead — that's what actually collapses it.
- **FRAMING-2 relations only work when bases nest.** Substituting every base
  (incl. the bare symbol `a -> s^3`) eliminates the symbol; the relation for the
  outer base `t^3 = a+b x` becomes `t^3 = s^3 + b x` (good — `s,b,x` survive). But
  for *independent* bases (e.g. `a^(1/3)` and `(a+b x)^(1/3)` with no other
  occurrence of `a,b,x`), the free vars vanish into the generators and no relation
  can be expressed in the ring. Guard each relation by ring-symbol membership and
  drop the unusable ones; the strict score gate then returns NULL (no regression).
- **Prove a test FAIL is pre-existing by neutralising the new hook, not by
  `git stash`.** New untracked files make `git stash push -- <paths>` abort, and
  `main` already carries a stray stash. Instead replace the new call with
  `Expr* rr = NULL;`, rebuild just the affected `*_tests`, and confirm the same
  FAIL count; then restore. (Did this for the 2 `Sqrt[x^2+6]`/`Sqrt[6]`
  `simplify_tests` soft-asserts — n=1 symbolic base, so the pass is inert.)
- **Corpus `*_tests` binaries load their `.m` via `../`-relative paths.** Run
  `fullsimplify_corpus_tests`/`crc_corpus_tests` from `tests/build/` (so `../` ->
  `tests/`); `intrat_corpus_tests` wants `IntegrateRationalTests.m` in the repo
  root. A "could not load … as a List" failure is a cwd issue, not a result DIFF.

## PossibleZeroQ false-zeros masquerade as integrator/algorithm bugs (2026-06-30)
- **Symptom**: `Integrate[(1-x^3)^(1/3)/x, x, Method->"GoursatAlgebraic"]` (and
  similar cube-root cases) silently declined, looking like a Goursat descent bug.
- **Actual cause**: two bugs in `zero_test.c`'s Schwartz–Zippel sampler made
  `PossibleZeroQ` return `True` for genuine non-zeros with an algebraic constant
  (radical / root of unity) times a free variable. The Goursat descent gates each
  eigenpiece with `is_zero` (`integ_backsub`), so a false-zero collapsed the
  answer to `0` → differentiate-back guard rejected it → decline.
  1. Samples clustered near 0 (`|n|≤64`, `d≤2^16` → `|val|≪1`); a `u^3` term
     vanished below the operand scale. Fix: bound `|sample| ≥ 1`.
  2. `magnitude_scale` scored `Power[denom,-1]` as `denom` (used `|exp|`), wildly
     inflating the cancellation-threshold scale. Fix: use the SIGNED exponent.
- **Lesson**: when a high-level algorithm (Integrate/Simplify/Solve) inexplicably
  declines or returns 0 on an input that is provably correct, suspect
  `PossibleZeroQ` FIRST. Reproduce the exact sub-expression the algorithm tests
  (instrument the gate, print it) and call `PossibleZeroQ` on it directly. A
  *variable-name-dependent* (flaky) verdict is a tell-tale of a sampler/seed bug.
- **Lesson**: a Goursat decline is not proof of non-elementarity. Verify the
  eigenprojection numerically (it's evaluate-only, reliable) before concluding.
  The ω-character criterion `H1==0` (cube) genuinely fails for `F` with a pole at
  a non-ramification point (e.g. `1/(x(1-x^2)^(1/3))`): that integral is
  elementary via Chebyshev's binomial mechanism, NOT the Goursat reduction — the
  WL reference `CubicRootElementaryQ` rejects it too. Hand-derived
  `H1 = -2(-4)^(-1/3) z/(1+z^3) ≠ 0` confirms it's inherent, not a code bug.

## Integrate hang on symbolic-exponent integrands (2026-07-09)
- **Symptom**: `Integrate[x^(k-1)(1-x)^(l-1),{x,0,1}]` (and indefinite
  `Integrate[x^k+x^(k-1),x]`, `Integrate[x^(k-1)(1-x),x]`) hung forever.
- **Root cause**: an integrand with a symbolic-exponent power of x (`x^k`,
  `x^(k-1)`) reaches `Together`/`Cancel` → `PolynomialGCD`, which treats the
  symbolic powers as independent polynomial generators and blows the
  pseudo-remainder sequence up (unbounded). Two entry points hit it: the
  rational-integration classifier `is_rational_in` (its `Together` *probe*, before
  BronsteinRational even runs) and the derivative-divides quotient
  `cancel_together`. `sample <pid>` pinpointed `pseudo_rem`/`poly_gcd_internal`.
- **Fix**: decline these Together-backed paths structurally, up front — a
  `Power[b,e]` with `b` depending on x and non-numeric `e` (`expr_is_numeric_like`
  rejects it) is provably not rational in x and never a productive u-substitution
  kernel. Guards: `has_symbolic_power_in` (integrate.c `try_rational`),
  `has_symbolic_power_of` (integrate_derivdivides.c `dd_core`).
- **Lesson**: when Integrate/Simplify/Together *hangs* (not wrong-answer), it is
  almost always PolynomialGCD/pseudo-remainder on symbolic exponents (see also
  `together_layer4_design.md`, the a^i Together hang). Use macOS `sample <pid> 2`
  on the live process to get the loop; the fix belongs at the *classifier gate*
  (reject before the expensive probe), not inside PolynomialGCD.
- **Process**: `./Mathilda` in pipe mode (non-tty) speaks NDJSON, not bare exprs —
  drive it with `{"id":1,"expr":"..."}` + `{"type":"quit"}`, and beware `| head -1`
  masking the real exit code (use `${PIPESTATUS}` / capture then grep).

## Completeness over aesthetics: don't gate out correct results (2026-07-10)

- **Correction**: For `Integrate\`RischMacsyma`, I added a `FreeQ[_, I]` gate that
  DECLINED `Tan[x]`/`Tanh[x]` because their coupled-hyperexponential answer, via
  the complex substitution `u = I x`, came out I-laden (`I x - Log[1 + E^(2 I x)]`
  `= -Log[Cos[x]]`) and no simplifier reduced it to real form. The user's rule:
  **favour completeness — return the correct antiderivative even when unsimplified,
  and FLAG the case as a `Simplify` improvement opportunity** rather than silently
  declining it.
- **Rule for myself**: a Risch/integration branch that is correct by construction
  must NOT be suppressed just because the output isn't in the prettiest form. If
  the answer is correct (diff-back 0 / SolveAlways-certified), ship it. Record the
  un-simplified shapes as explicit known-gaps (code comment + changelog +
  docs "Simplify improvement opportunity"), so the aesthetic fix can be done in
  `Simplify`, not by dropping capability from the integrator.

## Generalizing an RDE/ansatz solver: gate on genuine rational functions (2026-07-11)

- **Bug I introduced & fixed**: extending the base Risch-DE solver `rm_solve_rde`
  to rational exponents (Phase C, `E^(1/x)`) via a `q = h/Denominator[p]` ansatz, I
  routed every non-polynomial `p` to it. But `rm_exp_poly_case` passes the
  exponential's *coefficient* as `p` — for raw `E^x Sin[x]`, `p = Sin[x]`
  (transcendental). `SolveAlways` then certified a spurious `q = 0`, so
  `Integrate[E^x Sin[x]]` wrongly returned `0` (broke every multi-kernel test).
- **Rule for myself**: when a `SolveAlways`/denominator-theorem ansatz is fed a
  coefficient that may be transcendental, GATE it: require `p` (and `u'`) to be a
  genuine rational function of `x` — `PolynomialQ` on both `Numerator` and
  `Denominator` of `Together[·]`. A transcendental kernel (Sin, Log, another exp)
  must be rejected so the integrand falls through to the case that actually models
  it (expsum / trig front-end), never certified against a truncated ansatz.
- **Build gotcha**: after editing a `src/*.c` that the cmake test binary compiles
  via `COMMON_SRC`, a bare `make` (top-level) and `cmake --build` can run a STALE
  object if a `git stash`/`pop` reordered mtimes — I chased a phantom regression
  from a stale binary. `touch src/<file>.c` (or clean) before trusting test output.

## Eliminate inverse-substitution: don't intercept what the forward pass owns (2026-07-11)
When adding a new pre-pass that keys off a function head, check whether an
existing pass already handles that head *better*. The inverse-function
substitution pass initially included `Log`, which regressed `test_log_power`:
`u==Log[x]` got rewritten to `x->E^M` (a main-variable exponential the Groebner
atomiser can't decompose -> spurious nlin), whereas the pre-existing forward
exp/log algebraisation resolved `Log[x^n] -> n Log[x]` cleanly. Rule: a new
head-triggered transform must be disjoint from existing ones, or measurably
better on their cases. Always run the *existing* target test binary before
declaring done.

## No arbitrary caps / hacks in decision procedures (2026-07-11)
User reaction to shipping a Risch field-RDE degree bound behind an arbitrary
`if(d>5)d=5` cap (and a leftover `nmono<=128` ceiling): "Mathilda should be hack
free! The RDE solver should work for all degrees." A magic-constant degree cap in a
*decision procedure* is not an acceptable "increment" — it silently rejects valid
integrals. Rule: never introduce or leave a magic-number degree cap / resource
ceiling in the Risch (or any decision-procedure) code. Derive the bound from the
problem (here Bronstein RdeBoundDegree: leading-degree balance, monomial-type-aware —
log/x lower degree under D, exp preserves it), shared in a documented helper, with NO
ceiling. Correctness is already guaranteed by SolveAlways-certification + diff-back, so
the bound only affects completeness — which is exactly why it must be principled, not
capped. When a whole family of ansatz sites shares the same hack pattern, say so and
clean them all (or scope explicitly), don't leave siblings capped.

## Pattern rules can't re-use a var nonlinearly (a_^2) — bind linearly, Sqrt on RHS (2026-07-11)
Adding inverse-hyperbolic analogs to the CRC integral table, I first mirrored the
existing trig rules verbatim, e.g. `IntegrateTable[ArcSinh[a_. x_]/Sqrt[1 + a_^2 x_^2]]`.
These NEVER fired — and neither did the trig originals (453–458) they copied. I called it
a "matcher bug"; the user corrected me twice: "That's not a matcher bug, you need to match
to a_ and use Sqrt[a_] on the RHS of the rule," then "All rules in the table should be
updated in this way to fix this issue." Root cause: a pattern that binds `a_` from the
numerator (`ArcSinh[a_. x_]`, a_=2) and then writes `a_^2` in the denominator asks the
matcher to confirm `2^2 == 4` — it does NOT evaluate/invert pattern-var powers, so the
rule silently fails to match (except the trivial a=1 case via the optional default). The
rules "passed" in normal `Integrate` only because the general cascade (DerivativeDivides /
Risch) solved them by another route; through `Method -> "CRCTable"` they were dead.
Fix (apply to EVERY rule with this shape): bind the *quadratic* coefficient linearly as
`a_` using `c_ + a_. x_^2` (the `c_ +` form is required — `1 - a_ x_^2` won't bind the
coefficient because of sign fusion; `c_ + a_. x_^2` binds c=1, a=-4 cleanly), pin the
constant with a `Condition` (`c === 1`), link the numerator coefficient `b_` with
`a === ±b^2`, and recover the linear coefficient as `Sqrt[±a]` on the RHS. For x^n
recurrences use the optional exponent `x_^n_.` so bare `x` (n=1) matches and odd powers
bottom out (the n=1 recurrence term vanishes since its coefficient is n-1=0, and
`0*IntegrateTable[…] -> 0`). Rule: when a table rule "doesn't fire," first check whether the
pattern re-uses a bound variable under a nonlinear op (`a_^2`, `Sqrt[a_]` against a literal)
— the matcher can't invert those. Bind linearly + Condition + reconstruct on the RHS. And
verify rules fire through their ACTUAL dispatch path (`Method -> "CRCTable"` /
`IntegrateTable[...]` directly), not just via top-level `Integrate`, which can mask a dead
rule by solving it another way.

## Squared/cubed pattern constants never match numeric args — sweep the whole table (2026-07-11)
After fixing the inverse-trig/hyperbolic `a_^2 x^2` rules, the user pointed at Formula 488
(`Log[x_^2 + a_^2]`) and said "there are still many cases that will fail for the same reason
... We need to fix every rule in the table that has this issue!" Root cause (general): ANY
IntegrateTable pattern that reuses a constant under a power — `a_^2`, `b_^2`, `c_^3`, `c_^4`
— cannot bind against a numeric argument (the matcher does not invert the power), so the rule
silently never fires for concrete input; it only ever matched symbolic squares. ~169 rule
heads were affected across every family. Fix recipe: bind the constant linearly (`a_^2 -> a_`),
recover the linear part via `Sqrt[a]` on the RHS (even powers halve: a^2->a, a^4->a^2, a^6->a^3;
Abs[a]->Sqrt[a]). Signs: `x_^2 - a_` will NOT bind a_ to a negative literal (verified in Sqrt,
1/(...), (...)^(3/2), Log, and trig contexts), so match BOTH signs with `x_^2 + a_` — even-RHS
pairs merge into one unguarded rule (negative a reproduces the minus form); odd-RHS pairs split
on `Not[TrueQ[a<0]]`/`TrueQ[a<0]` with Sqrt[-a]. `Not[TrueQ[a<0]]` fires for positive-numeric
AND symbolic a but excludes negative, preventing the greedy plus rule from shadowing the minus
sibling. Many squared-constant rules turned out REDUNDANT with linear `a_. + b_. x_^n_` forms
that already fire (c^2 block 43-51, 62, 66, single-power c^3/c^4, 356) — delete those rather than
convert. Two pre-existing CRC transcription bugs (Formulas 181, 216) surfaced only once the rules
began firing — re-derive, don't faithfully copy the bug. Separate pre-existing issues NOT in
scope: `1/(...)^n_` reductions (pattern exponent: `Power[Power[E,n],-1]` != `Power[E,-k]`; explicit
`^2` works) and `1/(x^m Sqrt[...])` form-matching. Verify every rule fires through its real
dispatch (`Method -> "CRCTable"`), not top-level `Integrate` which can mask a dead rule by solving
it another way. General lesson: when a table rule "doesn't fire," first check for a bound variable
reused under a nonlinear op.
