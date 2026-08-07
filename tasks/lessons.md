# Lessons learned

## A coarse global cache key needs a finer *rule* epoch, not a rewrite (2026-08-05)

The evaluator memoizes fixed points against one global clock; every symbol-table
mutation bumps it, so an iterator's OwnValue rebind invalidated *all* cached
values and re-canonicalised loop-invariant data O(n) per step. The instinct is
"add dependency tracking" — a large, risky evaluator change. The elegant fix was
much smaller: a **second, finer epoch** (`g_last_rule_change_clock`) that advances
only on mutations which change how a *head* evaluates, plus a per-node "this value
depends on nothing mutable" flag (GROUND) stolen from a spare bit of an existing
field. A node that is GROUND and stamped after the last rule change is still a
fixed point regardless of the coarse clock.

Two rules this crystallised:
- **Correctness came from a whitelist, not a predicate.** "Any Protected head with
  no rules" *looks* safe but is unsound — `Plus`/`RandomReal` are Protected yet
  their builtins can read mutable global state, so a naive predicate could freeze
  a value that isn't actually constant. Restricting GROUND to six inert structural
  constructors (`List`/`Association`/`Rule`/`RuleDelayed`/`Complex`/`Rational`)
  whose canonical form is a pure function of args is what makes the freeze provably
  correct. When a "clever" predicate almost works, ask what mutable state each
  admitted case can secretly read.
- **Steal a bit before growing the struct.** The flag rode in the top bit of a
  64-bit monotone clock (unreachable 2^63), so `sizeof(Expr)` stayed 48 B on a
  system with millions of nodes. A new `uint8_t` field would have cost 8 B/node to
  alignment. The zero-cost encoding was worth updating the three comparison sites
  to mask.

## "Aware and slow" is a category the audit cannot see (2026-08-01)

`tools/check_packed_aware.py` answers one question: *does every head with an
NDArray fast path opt in?* It reported clean while `Extract`, `MatrixPower`,
`PseudoInverse` and `LeastSquares` were all on the `AWARE` list AND all
throwing the buffer away on their own first line — two of them by 3-4 orders of
magnitude. The opt-in is necessary and not sufficient.

The static check cannot close this: it looks for the *presence* of a dispatch,
and every one of those four had a dispatch (`linalg_delist_and_reeval` is a
dispatch). What sees it is `MATHILDA_PACK_DIAG=gate`, because its report covers
the **post**-gate — the materialisation that happens when a node comes to rest
with a buffer still in it. Run that before believing the audit.

Rule: when a head is on `AWARE` and still measures slow, do not re-read the
allowlist. Run `MATHILDA_PACK_DIAG=gate` and look at whether it appears; if it
does, the builtin is declining, not the gate.

## A fast kernel with a slow twin is where wrong answers hide (2026-08-01)

`sf_machine_productlog` returned **-338.392** for `ProductLog[1.01]` (answer:
0.5707) and nothing caught it, because the interpreter's scalar `ProductLog`
goes to MPFR and was always right. The two paths were never compared, so the
value tests all exercised the correct one. The *only* symptom was a speed row on
a coverage sweep: the array kernel failed on ~1 element in 10^5 and abandoned
the buffer to MPFR, which showed up as 28 us/element.

Two rules out of it:

1. **Whenever a head has both an MPFR path and a machine kernel, the test must
   name the machine one** — `Compile[{x}, f[x]]` or a packed array. An
   `assert_eval_eq("f[1.01]", ...)` tests neither.
2. **A kernel should verify its own answer when the algorithm is iterative.**
   Checking the defining equation costs one `log` against four Halley steps, and
   returning `false` is not an error — it is the kernel contract's "no usable
   value", which degrades to the correct path. "The iteration converges" is a
   claim about the starting guess, and starting guesses have domains.

## Making a head PACK breaks every internal consumer that walks its result (2026-08-01)

Adding `pack_offer` to `RowReduce`'s result was a one-line bonus on top of an
exactness fix. It silently broke `NullSpace`, `MatrixRank`, `PseudoInverse`,
`Inverse`, `Apart` and the eigen solver — six internal call sites that take
`RowReduce`'s output and walk it with `get_tensor_dims` / `flatten_tensor` /
`data.function.arg_count`. `get_tensor_dims` returns 0 for an `EXPR_NDARRAY`, so
`NullSpace` of a machine matrix came back **UNEVALUATED at 16x16 and worked at
14x14** — the break appears only above the 250-element packing threshold.

`pack.h` documents this precisely ("THE ONE GAP THE GATE DOES NOT COVER ... use
`pack_eval_plain` at any internal evaluate() whose result is then walked
structurally") and even lists the heads to watch. `RowReduce` was not on that
list because it had never packed before — which is the point: **the list is of
heads that pack TODAY, so adding a new one makes the list stale in the same
commit.**

Rules:

1. Before adding `pack_offer` to a head, `grep -rn "SYM_<Head>" src/` and fix
   every internal caller with `pack_eval_plain` in the same change.
2. Then add the head to pack.h's list, so the next person greps a current one.
3. **Test on both sides of the 250-element threshold.** A test at one size only
   proves whichever side it happened to land on, and the natural small test case
   lands on the safe side.

Found by asking "what else consumes this?" after the change was already green —
the full suite did not catch it, because no existing test crossed the threshold
for these heads.

## Reasoning by analogy from a fix that just worked is still guessing (2026-08-01)

Having found that `DiagonalMatrix`'s exact zeros were wrong against Mathematica,
I filed `RowReduce` as the same bug — "of a machine-real matrix it should be
uniformly machine-real, exactly as `DiagonalMatrix` should be". It reads like a
deduction. It was an analogy, and it was false: Mathematica's
`RowReduce[{{2., 4.}, {1., 3.}}]` is `{{1, 0.}, {0, 1}}`, heads
`{Integer, Real, Integer, Integer}`, `PackedArrayQ` False. Mathematica's own
RREF of a machine matrix is two-headed.

So within one day the SAME unverified-claim mistake happened twice, the second
time immediately after writing a lesson about the first. The pattern is not
"assume Mathematica agrees" — it is **stating a fact about another system in
prose instead of running it**, which costs ten seconds.

The durable fix is not a resolution, it is a mechanism:
`tools/check_array_exactness.py` will not accept an EXEMPT entry without the
Mathematica output pasted into it. If the claim cannot be quoted, it cannot be
relied on.

(The `RowReduce` change shipped anyway — but as a *stated divergence* from
Mathematica on the project's own rule, which is a different and honest thing
from an unnoticed one.)

## A constraint that excuses a slow path must be tested (2026-08-01)

`DiagonalMatrix` of a `Real` diagonal was 320x behind NumPy, and I wrote down
why: the off-diagonal zeros are exact `Integer`s, so the matrix has two heads
and no uniform buffer holds it — *"that is Mathematica's answer too"*. I put
that in a code comment, in the changelog, in `performance.md`, and filed the row
in the HPC plan against item 10.1 as a known design gap.

It took one `wolframscript` call to disprove. Mathematica gives
`{{1., 0., 0.}, {0., 2., 0.}, {0., 0., 3.}}` — all `Real`, packed array. The
*constraint was the bug*, and it was the only thing keeping the slow path alive.
Correcting the exactness made the matrix one dtype: 70.2 ms → 263 µs.

Rule: when a measurement is bad and the explanation is "we can't, because
correctness requires X", **verify X before writing it down**. A constraint that
conveniently excuses a slow path is exactly the one to check, and the check is
usually cheaper than the sentence justifying it. `wolframscript` is at
`/Applications/Mathematica.app/Contents/MacOS/wolframscript` on this machine.

Corollary: the same session had the same defect in `Subdivide`
(`Subdivide[0, 1., 4]` kept an exact `0`) and I did not notice it until the
Mathematica probe, because I had reasoned from the same wrong premise in both
places.

## An EXEMPT entry needs its reason, and then it pays (2026-08-01)

`ConjugateTranspose` sat in the packed-aware audit's `EXEMPT` table with:
"its NDArray path does not handle rank 1 and comes back UNEVALUATED as
Conjugate[Transpose[v]]". That note turned a 195x-vs-NumPy row into a ten-minute
fix, because it named the precondition instead of just recording a decision.
Keep writing them; the difference between "considered and rejected" and "never
noticed" is worth more than the exemption itself.

## Do not reproduce an Orderless fold's rounding (2026-08-01)

`Subdivide` built `Times[i, span, Power[n, -1]]` per point and let the evaluator
fold it. `Times` is `Orderless`, so its factors are sorted **by value** and
folded left to right — meaning the exact-to-Real transition happens at whichever
factor the `Real` operand happens to sort after, and the last bit of the answer
depends on the magnitudes of the inputs. There is nothing there to preserve.

When replacing such a path, do not try to match it bit-for-bit (I tried; 349719
of 10^6 points agreed). Pick the rule the reference systems use, check
bit-for-bit against **them**, and write down in the code that the change was
deliberate. `numpy.linspace` and Mathematica both compute `min + i*step`, and
Mathilda now agrees with `linspace` on every element.

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

- **Correction**: For `Integrate\`RischTranscendental`, I added a `FreeQ[_, I]` gate that
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

- **Bug I introduced & fixed**: extending the base Risch-DE solver `rt_solve_rde`
  to rational exponents (Phase C, `E^(1/x)`) via a `q = h/Denominator[p]` ansatz, I
  routed every non-polynomial `p` to it. But `rt_exp_poly_case` passes the
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

## Bronstein Risch canonical rep (P0) — 2026-07-12

Building `SplitFactor`/`SplitSquarefreeFactor` (risch_field.c / risch_canonical.c):

- **`Together` does not expand its numerator.** A difference that is algebraically 0
  can come back as an unexpanded, canceling sum. For zero-tests use
  `Expand[Together[a - b]]`, not `Together[a - b]`. (First SplitFactor test "failed"
  on a result that was actually correct.)
- **Field gcd in `k[t]` must special-case a zero operand.** `PolynomialGCD[p, 0]` does
  NOT reliably return `p` here; `gcd(c, 0)` came back as `1`, which made Yun's loop
  never reduce `deg(c)` → infinite loop → `out[]` buffer overflow → heap corruption
  (SIGABRT "free list damaged", detected asynchronously so the backtrace was useless).
  Handle `gcd(a,0)=a`, `gcd(0,b)=b` explicitly (degree `-1` marks the zero poly).
- **Splitting-factorization output is only unique up to units of the field `k=C(x)`.**
  A monic-in-`t` gcd makes the factors monic (`t-1`, not Bronstein's `4x^2(t-1)`), and
  the product of monic factors reconstructs `monic(p)`, not `p` (differ by `lc_t(p)`).
  Test via reconstruction × content and normal/special classification, or pin the
  canonical (monic) special part; don't string-match Bronstein's content-carrying forms.
- **`expr_copy` is a refcount bump, not a deep copy** (`e->refcount++`); `expr_new_function`
  copies the args *array* and *adopts* the element references. So `free(arr)` after
  building a function node is correct and must NOT be paired with freeing the elements.
- New `src/calculus/*.c` is auto-globbed by the makefile but must be added to
  `tests/CMakeLists.txt` `COMMON_SRC` (and the test target added) or the test won't link.

## P1 Phase C — structure-oracle tower replacement is net-negative (2026-07-14)

**Pattern:** Before a "full replacement" of a working heuristic, PROBE whether the
gap it targets is real in the live system. The gaps doc flagged G-A2 (dependent
logs inflate the tower), but the live Risch tower engine already tolerates
dependent generators (each becomes its own tower variable with a valid triangular
derivation, gated by diff-back). A blanket structure-oracle log collapse
REGRESSED the Bronstein `Log[x/Log[x]]` example — a composite log inside a
`1/Log[...]^2` denominator became a coupled two-log denominator the integrator
declines.

**Rule for myself:** (1) A "replacement" that rewrites integrand kernels can turn a
clean single kernel into a harder coupled form — collapsing is only safe for
redundant EXTRA generators, never for kernels appearing in denominators.
(2) Capture a regression baseline (anchor integrals + suite pass/fail) BEFORE
touching a core subsystem, and revert immediately on the first regression rather
than patching around it. (3) Decision routines (LogReducible / ExpReducible /
LogarithmicDerivativeOfRadical) are valuable as reusable builtins even when the
wholesale live-oracle swap is not worthwhile.

## 2026-07-14 — Don't patch Risch deficiencies with Weierstrass
Directive (user): rational trigonometric integrals must be integrated by the
GENUINE Risch machinery — the complex-exponential tower (TrigToExp → integrate over
the E^(ix) monomial → ExpToTrig with REAL reconstruction of the log part) and the
real §5.10 hypertangent case — NOT by the Weierstrass t = Tan[x/2] half-angle
substitution and NOT by routing through Jeffrey–Rich. The half-angle form is ugly
(raw Tan[x/2]) and, more importantly, papers over the real deficiency instead of
fixing it.

What went wrong: I added an rt_weierstrass_case (t=Tan[x/2]) dispatch to close
Sec/Csc/1-over-(2+Cos) and then promoted it ahead of Jeffrey–Rich at the top level.
It regressed Sin[x]^2 (clean (2x-Sin[2x])/4 from the genuine exp route → ugly
Tan[x/2] form) and was a patch, not a fix. Reverted the whole thing (branch reset).

Rules for myself:
1. When an integrand class is handled poorly, FIX the genuine algorithm (real
   reconstruction of the complex log part; complete the rational-of-E^(ix) route;
   relax an over-strict §5.10 gate) — do NOT bolt on a substitution that competes
   with the genuine route.
2. Before "improving" trig integration, capture what the GENUINE Risch route
   (Integrate`RischTranscendental) already produces on clean main — much of it is
   already clean (Sin^2, Sin^3, Tan). Only the true gaps (I-laden Csc, declined
   Sec) are the work.
3. Never promote a new stage ahead of an existing one without checking it does not
   intercept cases the existing stage handles more cleanly (Sin^2 regression).
4. REPL JSON output: parse with python (json), not sed — sed mangles Floor[...] and
   `/` in payloads and falsely reads them as empty/declined.

## Simplify I-laden collapse fallback + masked-failure lesson (2026-07-16)

Task: make Simplify prove `x E^x Sin[x]` Risch diff-backs zero. Root cause: the
exact grid zero-test (`trigexp_rational_is_zero`) declines BY DESIGN on bare
polynomial dependence on the kernel var and on mixed real+imaginary exp kernels
(`E^((1+I)x)=E^x·E^(Ix)` read as independent opaques). Fix: `transform_trigexp_vanish`
now, ONLY when that test returns UNKNOWN (never FALSE) and complexity ≤ 512, tries
`TrigToExp[e]`; a literal `0` (via the evaluator's automatic E^a·E^b→E^(a+b) merge,
NO Together/grid) is a proof. Sound, cheap, leak-free.

Rules for myself:
1. When extending a decision-procedure fast path, gate the new fallback on the
   rigorous test DECLINING — never run it when the rigorous test already proved a
   verdict (wastes work) and never unconditionally (here: unconditional TrigToExp
   reintroduces the exp-form blowup the grid test exists to avoid — the UNKNOWN+cap
   gate is load-bearing).
2. MASKED FAILURES: an assert-active test suite ABORTS at the first failure, hiding
   every later test. Fixing an early failure UNMASKS later ones — that later red is
   NOT necessarily your regression. Prove it: reproduce the newly-visible failure on
   a clean worktree of the true committed HEAD (`git worktree add /tmp/x <sha>`).
   Here Csc[x]^3 returned unevaluated on HEAD 0ee179d too → pre-existing + flaky,
   not caused by my change nor the earlier modularization (which the abort had also
   shielded from ever checking test_circular_trig_integration).
3. Pipe vs C-test divergence: a case can decline deterministically in a cold pipe
   but pass in the warm C-test process (memo/symbol state). Don't conclude "broken"
   from the pipe alone.

## A "clean" fast-path answer can hide a divergence that cancelled (2026-07-28)

`Limit[Cos[1/x] - Cos[1/x + 1], x -> 0]` returned `0`. Substituting `x -> 0`
makes both terms `Cos[ComplexInfinity]`; the `Plus` cancels them; the
`is_divergent(result)` guard sees a clean `0` because by then the divergence is
gone. The function oscillates over `[-2 Sin[1/2], 2 Sin[1/2]]`.

Rule: when a fast path justifies itself with "the substituted answer came out
finite", the guard belongs on the **inputs to the fold**, not the output. Plain
substitution is a limit only where `f` is *continuous*, and a divergent inner
argument IS discontinuity. `has_divergent_inner_arg_at` in `src/calculus/limit.c`
walks every argument of every function node that mentions `x`, substitutes it
alone, and refuses if it diverges — on the **`Together`-cancelled** form, not
the original, or a removable singularity like `(x^2-1)/(x-1)` at 1 gets refused
too (the first cut did exactly that, and `limit_tests` caught it). The same failure shape lurks anywhere the
evaluator normalises infinity-valued subterms — a cancelled `Infinity -
Infinity`, `0 * ComplexInfinity`, `1^ComplexInfinity`.

## A flatness invariant an un-re-evaluated builtin result can violate (2026-07-28)

`TrigToExp[x Sin[x]]` returns `Times[c, Times[x, E^(I x)]]` — a nested `Times`,
even after `Expand`. The evaluator flattens `Times` when *it* evaluates, but a
builtin that hands back a freshly built tree can leave the invariant broken.

Rule: any code that walks `Times`/`Plus` factors structurally must recurse
through nested same-head nodes rather than reading `args` once. In
`limit_osc.c` this bit hard and silently: `Sin[x]` classified correctly while
`x Sin[x]` dropped its `E`-factor and mis-classified the whole limit. The
symptom was "works on the simple case, abstains on the case with one more
factor" — which reads like a missing rule, not a collector bug.

## PossibleZeroQ before a structural check turns a decay into a wrong answer (2026-07-28)

`limit_osc.c` pruned zero-amplitude groups from the normal form *before*
checking that every amplitude was oscillation-free. `PossibleZeroQ` samples
numerically, so `(2 Sin[t])^(t^2)` — which underflows wherever `|2 Sin t| < 1`
— read as zero. The only group got dropped, the normal form came out empty,
and `Limit[(2 Sin[1/x])^(1/x^2), x -> 0]` answered a confident `0` for a
function that is unbounded.

Rule: **structural gates run before numeric ones.** A numeric zero test is a
heuristic that can only be trusted on inputs a structural check has already
vouched for. Where the order is forced the other way, ask what a false "zero"
would do — here a false zero on a genuinely *decaying* amplitude is harmless
(it is exactly what the squeeze rule would have done), so only the impure case
needed protecting.

## `make | grep ... | head -N` can SIGPIPE the build and leave a stale link (2026-07-28)

Filtering a build with `make -j8 2>&1 | grep -E "error" | head -3` tears the
pipeline down early and can kill `make` **after** the objects compile but
**before** the link. Twice in one session that left `./Mathilda` carrying an
older `limit.o` with the new layer absent — and because the binary's mtime then
matched the objects', a follow-up `make` printed "Nothing to be done for
`all'" and changed nothing.

Symptom to recognise: a feature that worked five minutes ago behaves as if it
was never written, and re-running `make` prints nothing. The tell is
`ls -la ./Mathilda src/**/changed.o` showing the **binary older than an
object**.

Rule: never pipe `make` into `head`. Redirect to a log
(`make -j8 > /tmp/build.log 2>&1; echo rc=$?`) and grep the file afterwards.

## A soft-assert `FAIL:` prints at the TOP of the log and still exits 0 (2026-07-28)

`test_limit_assumptions.c` had been failing since the previous commit —
`Limit[x^n, n -> Infinity]` was pinned to the old `E^DirectedInfinity[Log[x]]`,
which the `exp_of_limit` fix correctly turned into an honest unevaluated form.
I ran the suite, saw `All limit_assumptions tests passed!` in a `tail -4`, and
committed.

Two traps compounded:

1. The suites are built with `NDEBUG`, so a failed assertion prints
   `FAIL: <input> / Expected: ... / Actual: ...` and **keeps going**. The exit
   code is 0 and the closing "All ... passed!" banner still prints.
2. The failure was in the *first* test, so it scrolled off the top. A `tail` of
   the log is exactly the wrong window.

Rule: judge a suite by `grep -c 'FAIL'` over its **whole** output, never by the
exit code and never by the tail. Run it as
`./t > /tmp/t.log 2>&1; grep -n FAIL /tmp/t.log` and read the count.

Corollary: when a change alters what a limit/integral *returns* in an edge
case, grep the whole test tree for the old printed form
(`grep -rn "E^DirectedInfinity" tests/`) before assuming nothing pinned it.

## When a fix removes a malformed output, grep for who was consuming it (2026-07-28)

Making `exp_of_limit` refuse a residual infinity was correct — `Limit[E^(I x)/x,
x -> Infinity]` had been answering `E^DirectedInfinity[I]` instead of `0`. But
`E^DirectedInfinity[dir]` was not merely noise: the Frullani pre-pass in
`integrate_ramanujan.c` *pattern-matched on it*, resolving `dir < 0` under the
assumptions to get `f(Infinity) = 0`. Deleting the form silently broke
`Integrate[(E^(-a x) - E^(-b x))/x, {x, 0, Infinity}]` on that route.

A malformed intermediate value is load-bearing more often than it looks — it is
the only channel through which the producer's internal knowledge reached the
consumer.

Rule: when a change stops a function returning some distinctive shape, grep the
tree for that shape (`grep -rn "DirectedInfinity" src/`) before concluding the
fix is local. Then give the consumer a *legitimate* way to get what it needs —
here, the limit of the exponent, which is decidable — rather than restoring the
bad output or leaving the consumer broken.

## Two spellings of one integrand is a routing bug, not a coverage gap (2026-07-28)

`Integrate[E^(1/2 Log[Log x] - 1/Log x)/(x Log x^2), x]` closed to
`-Sqrt[Pi] Erf[1/Sqrt[Log x]]` while the *identical* function written
`E^(-1/Log x)/(x Log x^(3/2))` came back unevaluated. It read like a missing
case in the Knowles erf engine. It was not: the engine was never called. A
scope gate four layers up (`rt_has_algebraic_of_x` in
`builtin_rischtranscendental`) saw `Power[Log[x], Rational[-3,2]]`, classified
it as an algebraic extension, and returned NULL before the dispatcher ran.

Rule: when two spellings of the same expression give different answers, the
defect is almost always in **routing**, not in the engine that should have
answered. Find out whether the engine ran at all before reading a line of its
code — a one-line `getenv`-gated `fprintf` at each entry point localises it in
one build, where reading the engine can burn an hour on the wrong file.

Corollary on relaxing a soundness gate: do not widen the predicate. Admit the
narrow case, transform it into the form the working spelling already takes, and
require an independent certificate (here a diff-back against the *original*
integrand) — then withhold anything the relaxed field cannot justify, which for
a Risch path means never emitting the non-elementary certificate from it.

## A wrong number and a slow number are the same defect (2026-07-29)

`N[Sin[3141592653589793238]]` answered `-0.641653` instead of `-0.446315`.
The root cause was one line — machine-mode `N` converted the exact leaf with
`(double)v` before `Sin` ever saw it — but the sweep written to *test* the fix
found four more defects in the same family, three of which nothing else would
have caught:

- `N[Exp[1000]]` was `inf.0` (machine numbers have a 53-bit mantissa but an
  arbitrary exponent; the codebase already knew this for `N[1001!]`);
- `N[Gamma[2^53+1]]` **aborted the process** inside GMP;
- `ArcSec`/`ArcCsc`/`ArcCoth` took the wrong branch at machine precision;
- `Sqrt[3141592653589793238]` took **24 seconds** (trial division to
  `sqrt(n)`), which the sweep surfaced only as a timeout.

Rules this session earned:

1. **Report an unexplained slowdown as a defect, not as a test-tuning
   problem.** The first move on `Sqrt@pi-digits` timing out was to exclude it
   from the suite so the suite stayed under its alarm. That is backwards: 24 s
   for something `FactorInteger` does in 3.6 ms is the bug. Excluding it hid a
   ~3400x regression that a one-line trial-division bound fixed. If a case has
   to be excluded to keep a suite fast, first prove the case *deserves* to be
   slow.

2. **Give a stress harness a per-case timeout and a fork, from the start.** A
   suite whose whole job is hostile inputs will meet a hang (`N[Erfi[2^53], 30]`
   never terminates) and a crash. Run each case in a forked child with its own
   alarm so both become *recorded outcomes*; a bare `alarm()` in the parent
   kills the run and tells you nothing about which case did it. A one-line
   `current_case` buffer printed from the SIGALRM handler localised the hang in
   a single run.

3. **`expr_to_mpz` initialises its target.** `mpq_init(q)` followed by
   `expr_to_mpz(x, mpq_numref(q))` allocates limbs and immediately leaks them.
   Valgrind only names this if you compare against a *baseline run of the same
   binary* — the macOS startup noise (~13,376 bytes / 418 blocks) is large
   enough to swallow it otherwise. Always diff the leak summary of the workload
   against the leak summary of `1+1`.

4. **The invariant `N[f[x]] == N[f[x], 30]` is worth more than a reference
   table.** It needs no oracle, it applies to every numeric builtin uniformly,
   and it is precisely what the bug violated. Pair it with a gap file that fails
   on both an *unlisted* gap and a *stale* one, so the list can neither grow nor
   rot silently.

## Compiled functional programming (2026-07-29)

- **A compiled fast path must never answer where the interpreter declines.**
  Four heads needed a runtime guard for exactly this: `Fold` over `{}`,
  `Nest`/`FixedPoint` with a negative count, and any unbounded iteration at the
  10^6 cap. In each the interpreter leaves the whole expression unevaluated,
  where a counted loop would silently return the seed. Before compiling a head,
  read its interpreter implementation for the cases that return `NULL`.

- **A branch guard's test must be the FAILING condition, not the good one.**
  `JZ` skips the following instruction when the test is false, so
  `GE_I rc, n, 0; JZ rc -> skip; FAIL` fails on the *good* path. Both guards
  shipped inverted and the disassembler showed it in one line — `CompilePrint`
  before theorising about a wrong answer.

- **Check that a failure signal is actually consumed.** `compiled_eval_real`
  computed `failed` and never read it. Harmless for as long as no opcode an
  all-Real program could contain was able to fail; the first such opcode made it
  a wrong answer. When adding a way to fail, grep every entry point for the flag.

- **Compiling a head can require fixing the INTERPRETER first.** `Fold` left an
  `NDArray` unevaluated (an NDArray is atomic, so the element walk missed it),
  so there was nothing to be parity with. Same for `Select`, `Join`, `First`,
  `Differences`, `RotateLeft/Right`, `Riffle`, `Partition`, `TakeWhile`,
  `AllTrue`/`AnyTrue`/`NoneTrue`.

- **Benchmark the compiled array path with a PACKED argument.** Over a plain
  `List` a compiled `Map` measured 1.0x — both sides dominated by packing 200k
  `Expr` nodes at the boundary, and the "interpreted" side already using the
  legacy `numloop` fast path. Packed in and out, the same body is 277x.

- **A test helper's variable name is part of its contract.** `ref_at` binds
  `xq`; eight `Table` parity cases written with `x` all reported "shape/kind
  mismatch" against a reference that was still symbolic.

## NDArray gaps and the selection heads (2026-07-29)

- **An `NDArray` is ATOMIC, so any builtin that walks
  `arg->data.function.args` silently ignores one** and returns the call
  unevaluated — while the identical `List` call works. Seventeen heads were in
  that state (`Select`, `Join`, `First`, `Differences`, `Riffle`, …). The fix
  that generalises is `ndstruct_delist_repack`: materialise, reuse the ordinary
  List implementation, repack. The ANSWER then comes from the List path and is
  identical to it by construction — only the packing is new.

- **A test helper's comparison decides what it can see.** `arr_cmp` reduces both
  sides to `(Re, Im)`, so it reports a Boolean result as a "shape mismatch". Six
  passing behaviours looked like failures until the assertion changed to
  `expr_eq`, which is also stricter — it compares heads and dtypes, so a Real
  answering where an Integer should now fails rather than compares equal.

- **A legitimate decline is not a test failure.** `parity_arr` counts any
  declined trial as a failure; `Select`/`TakeWhile` decline on an empty result by
  design. Give such tests an explicit `may_decline`, or the harness pushes you
  toward making the code answer where it should not.

## Retargeting an engine: measure the premise first (2026-07-29)

- **"The new engine dominates the old one, so delete the old one" is a claim to
  measure, not to assume.** Mine was false: the Compile[] engine beats
  `src/numloop.c` on `Nest` (2.5x) and `While` (1.6x) but LOSES on `Do` (1.2x)
  and on `Map` over a plain List (1.1x) — the latter because it must pack 200k
  `Expr` nodes at the boundary and unpack them again, which numloop skips by
  walking the List. A wholesale retarget would have slowed the two most common
  loop constructs. Benchmark both engines on identical bodies before planning a
  removal.

- **Read the disassembly before blaming dispatch.** The compiled `Do` was losing
  on INSTRUCTION COUNT, not VM speed: four of eight inner-loop instructions were
  pure loop control, and `OP_LOOP` (increment+test+branch in one) already existed
  and was used by the fused array loops but by no counted loop. One `CompilePrint`
  showed it; no amount of theorising would have.

- **A test that pins codegen shape has to move with the codegen.** An assertion
  that `INC_I` appears in a `Do` disassembly failed when the loop started closing
  with `OP_LOOP`. The right fix was to assert the NEW, more specific opcode (and
  to add a non-unit-step case asserting the old one), not to weaken the check.

- **`#ifdef` on a name a header `#define`s as a VALUE is always true.** The
  overflow build switch was first called `COMPILE_WRAP_INT`, which is also the
  public flag *bit* in `compile.h`. `#ifdef COMPILE_WRAP_INT` in the VM was
  therefore true in every build, the checks were compiled out, and every test
  silently exercised the wrap path — the symptom was "my new code is not in the
  binary", which sent me hunting stale objects for a while. Name a build switch
  after the layer it controls (`VM_NO_INT_CHECK`, like `VM_NO_THREADED`) so it
  cannot collide with a flag value.

- **Before delegating a new element type to an existing library function, check
  what it ACCUMULATES IN.** `NDT_INT64` arrays reached `ndred_total_all` and the
  `flatten_into` packer, both of which route every element through a `double`.
  Nothing errored; `Total[{9007199254740993, 1}]` just came back one short. Exact
  paths (`ndt_get_i`/`ndt_set_i`) had to be added at both ends. A lossy
  delegation is worse than a bail, because a bail is visibly slow and this is
  invisibly wrong.

- **A bail test proves nothing unless you know WHICH bail you got.**
  `must_bail_raw("NestList integer body", …, in2, RI2, 1)` passed for years
  because nargs=1 left `n` a free symbol — it was asserting "a free symbol
  bails", not the integer-history rule its comment claimed, and would have kept
  passing whatever happened to that rule. `CompileDiagnostics` reports the actual
  cause; assert on it, and add the opposite-direction test (`must_compile_raw`)
  when a restriction is lifted.

- **Result HEADS need their own sweep; value tests are structurally blind to
  them.** `35` and `35.` compare equal, so no amount of numeric parity testing
  could see that ten integer-closed heads were coming back as Reals. A 20-line
  sweep over every `NumericFunction` head comparing `Head[h[3]]` against
  `Head[Compile[…][3]]` found all of them in one run, and separated the real bugs
  (integer-closed) from the inherent divergences (symbolic, Rational) — which
  hand-enumeration would not have.

- **When a feature flag turns out to buy nothing, say so.** The
  `"CatchMachineIntegerOverflow" -> False` option was asked for partly for speed;
  measured, it is 0% because the design bakes the choice into each instruction so
  the fast path is byte-identical. Reporting it as a semantics switch (and naming
  the build flag that DOES recover the 4%) is more useful than quietly shipping a
  knob that does nothing for the stated reason.

## Automatic packed arrays (2026-07-30)

- **A representation change is invisible to value tests by construction, so
  differential sweeping is the only thing that finds its bugs.** Switching on
  automatic packing exposed six wrong answers, none of which was visible by
  reading the code. All six came from one mechanical technique: evaluate the same
  expression over a packed value and over the identical plain list
  (`MATHILDA_NO_PACK=1`), for every packed-aware head (134 cases) and every
  registered ndarray kernel (166 cases), and diff. Hand-enumeration found none of
  them; the sweep found all of them in two runs. When adding a second
  representation for an existing value, write the sweep before the producers.

- **An invariant enforced by mutation is not safe under a convergence test that
  cannot see the mutation.** The transparency gate materialises a buffer in place
  and reports `*changed`. But `evaluate`'s fixed point is `expr_eq(current,
  next)`, and `expr_eq` is *deliberately* blind to packing — that blindness is
  what makes Association lookups work. So a step whose only effect was
  materialising looked like no progress and was discarded. Cost: a nested `Table`
  came back as a List of packed rows. If a pass's only effect is invisible to the
  comparison that decides whether to keep it, the pass needs its own signal.

- **An implicit opt-in is the opposite of an audited claim.** `packed_aware` was
  meant to mean "I have been checked against a buffer". But the ndarray kernel
  setters set it as a side effect, so ~80 heads claimed it without anyone
  checking, and six of them (`Floor`, `Ceiling`, `Round`, `IntegerPart`, `Sign`,
  `Im`) answer with different element HEADS than the list does. Having a machine
  kernel is a different claim from matching the list's heads; the registration
  conflated them. Needed a `symtab_clear_packed_aware` to undo.

- **A namespace prefix is not a membership test.** The `$AutoCompilation` /
  `$AutoArrayPacking` assignment hook keyed on the `$` sigil for cheapness, then
  built a probe value — evaluating a delayed right-hand side. The REPL hooks
  (`$Pre`, `$PreRead`, `$Post`, `$PrePrint`, `$Epilog`) share that namespace with
  held right-hand sides, so all of them began evaluating at definition time.
  Check the name against the actual set *before* doing anything with a cost or a
  side effect.

- **A performance gate calibrated on ordinary work is a landmine for any change
  that makes ordinary work faster — and it will be in more than one file.**
  Phase 0 moved `bench_eval`'s calibration off `Total[Range[40000]]` precisely
  because packing would make it ~20x faster. `bench_assoc` had the *same*
  divisor, was missed, and failed all nine of its operations at once with nothing
  actually slower. When defusing a calibration landmine, grep for the pattern
  rather than fixing the file that prompted it.

- **A performance win in one direction can be a regression in the other, and only
  measurement tells you which.** Packing exact-integer lists made `Length` 87000x
  faster and `Total[Range[10^6]]` **1.55x slower**, because the safety gate had to
  materialise every integer buffer before any reduction. The int64 exact paths
  were filed as "pure optimisation, not correctness" — switching on the producers
  promoted them to load-bearing. Benchmark the whole surface after a
  representation change, not just the case that motivated it.

## The Compile[] boundary and packed arrays (2026-07-30)

- **When a fast path is mysteriously not firing, check the GATE before the path.**
  A packed argument to a `CompiledFunction` measured 75x slower than the identical
  visible `NDArray` — two values differing in one enum field. The boundary code was
  fine; the transparency gate was materialising the buffer before
  `compiled_function_apply` ever saw it, because a `CompiledFunction`'s head is an
  `EXPR_COMPILED` with no `SymbolDef` and the gate's allowlist is keyed on the
  symbol table. An allowlist cannot see a head that is not a symbol. Diagnosing
  this by reading the boundary code would have taken arbitrarily long; comparing
  two values that differ only in presentation found it in one measurement.

- **"Derived" and "produced" are different rules and must not share code.** A
  derived array inherits its source's presentation with NO size threshold —
  `Sin[packedList]` is packed however short it is. A producer applies the
  threshold. Using the producer opener for a derived result made
  `Map[#^2 &, ToNDArray[{1., 2., 3., 4.}]]` come back unpacked while `Sin` on the
  same value stayed packed. Two openers (`ndbuild_open` vs `ndbuild_open_like`),
  not one with a flag.

- **One boolean answering two questions is a latent bug waiting for a third
  case.** `packed[i]` at the compile boundary meant both "we own this" and "the
  argument was not an NDArray, so unpack the result". That was sound while a
  borrowed `EXPR_NDARRAY` could only be the visible head. Packing added a third
  input kind and a dtype cast that makes a temporary out of a borrowed argument —
  and the two meanings came apart. Split it into ownership + a three-way kind.

- **Fix the wrong answer in the area before you widen its reach.**
  `Compile[{{u, _Integer, 1}}, u * 2][{1, 2, 3}]` gave `{2., 4., 6.}` — an array
  opcode's scalar operand could only be Real or Complex. Pre-existing, nothing to
  do with packing, and reproducible on a three-element plain List. But returning
  int64 buffers from the boundary would have made it far more visible, so it was a
  prerequisite rather than a separate ticket.

- **Measure the workload you made slower, not only the one you made faster.**
  Automatic packing made `Map[#^2 &, x]` at 10^6 **1.9x slower** and
  `Map[Sin[#] Exp[-#] &, x]` **6.2x slower**, because a packed argument bypassed
  numloop's compiled loop (it tested `EXPR_FUNCTION`) and fell to a per-element
  interpreter walk. Every value test passed. The regression was found only by
  timing the head that matters most on both representations — which is now a habit
  worth keeping for any representation change.

## A missing fast-path opt-in has no symptom — only a timing comparison finds it (2026-07-30)

Writing `docs/design/performance.md` (38 HPC kernels vs Mathematica) found more
defects than any code review of the same subsystem had. Every one was *correct*
and *quiet*: right answer, slow path, no failing test.

The packed-array transparency gate's rule is "materialise for any head that has
not opted in". That makes a missing opt-in fail **safe**, which is the right
design — and it also makes it **invisible**. `Nest` was materialising a 512×512
buffer on every iteration (2.19 s where the same `Do` loop took 0.017 s) and
nothing anywhere said so.

Rules for next time:

1. **When a subsystem is gated by an allowlist, the allowlist is a to-do list.**
   Enumerate the heads that dispatch on the fast path (`grep -l` for the guard)
   and diff that against the allowlist. Three of this round's fixes were exactly
   that diff: 26 linear-algebra heads, the `Nest` family, the structural family.
2. **A missing opt-in can be a crash, not just a slowdown**, when the same test
   selects both the fast path and the correct algorithm. `linalg_call_has_ndarray`
   read false after the gate materialised, so a machine-real solve took the exact
   fraction-free path and overflowed the stack.
3. **Ask the diagnostic tool about the bodies you already believe compile.**
   `CompileDiagnostics` reported `"Compiled" -> False` on a Sieve and a Collatz
   search for one missing form (two-argument `If`). It had existed the whole time;
   nobody asked it.
4. **Benchmark against another system, not against yourself.** A self-comparison
   ("packed vs unpacked") cannot show that both arms are 100× off. Every gap above
   was visible only as a ratio against Wolfram.
5. **A benchmark can be wrong the way a program can.** The Lennard-Jones point set
   built from `Mod[a k, m]` cycles, so it repeated *points*; two coincident bodies
   overflowed the 1/r¹² term, the compiled call correctly bailed on the non-finite
   result, and it read as an 89× performance cliff. Put an answer next to every
   timing — the value column is what caught it.
6. **Measure wall clock.** `Timing[]` is CPU time summed over threads and
   over-reports every threaded/BLAS path by ~the core count. `AbsoluteTiming` had
   to be added before any of this could be measured honestly.

## A fast path that bypasses initialisation is worse than a slow one (2026-07-31)

`RandomInteger[{1,10}, 300]` **hung** on `main`, and `RandomInteger[{1,1024}, 300]`
answered 300 copies of `1`. Both needed the same two conditions: `n` past the
packing threshold, and no earlier random call in the process.

The cause was a fast path added for speed that drew from the xoshiro generator
directly, where every per-element helper called `ensure_rand_init()` for itself.
Against an all-zero state — xoshiro's one fixed point — `xs_next()` returns 0
forever: a power-of-two span silently answers the range minimum, and any other
span spins in the rejection loop.

Rules:

1. **When adding a fast path, list what the slow path does per element that the
   fast path now does once — or not at all.** Seeding, validation, precision
   checks, error reporting. Each is a candidate for exactly this bug. The slow
   path here called `ensure_rand_init()` per element; the fast path called it
   never, and that difference is invisible in a diff that only shows the new loop.
2. **Put the invariant where it cannot be bypassed, and measure the cost rather
   than assuming it.** `ensure_rand_init()` now lives inside `xs_next()` itself.
   The branch was already being paid by `RandomReal` per element, which still
   draws at 2.3 ns — so the argument against it was never true.
3. **Do not make a loud failure quiet in the name of robustness.** The first fix
   also gave the generator state a plausible nonzero default. That removed the
   hang, and it removed the *regression test's ability to fail* — an unseeded
   generator would then produce a valid but identical stream on every run, which
   is a far worse bug than a hang. Reverted: the state is deliberately left zero.
4. **A regression test for a fresh-process bug must run first.** In
   `test_random.c` any earlier test seeds the generator and makes the check
   vacuous. Its position in `main()` *is* the test, and it says so in a comment.
   Verified by removing the fix and confirming it hangs.

## Per-call costs that scale with the data are invisible until you vary the data (2026-07-31)

`Interpolation` evaluation cost time proportional to the size of its own table —
151 µs per point at 10⁴ nodes, 906 µs at 5×10⁴ — so resampling was quadratic and
a 10⁵-node table could not be benchmarked at all (228 s). Four separate causes,
and the interesting part is that a previous round had already found and fixed a
*fifth* one, leaving a comment about it in the file.

Rules:

1. **Hold each input dimension fixed in turn.** One timing tells you nothing.
   Table size fixed / points varied, then points fixed / table varied, separated
   "cost per point" from "cost per call" in two measurements and pointed straight
   at the remaining O(n²) in the grid build.
2. **Profile before the second fix.** The first fix (memoising the value tensor)
   won 27× and the row still looked quadratic; guessing again would have been
   cheap and wrong. `sample` named `interp_apply` and `common_scan_inexact`
   directly.
3. **A memo keyed on `Expr` identity is sound here and O(1) where `expr_eq` is
   O(n)**, because `expr_copy` is a refcount bump — the cache holds a reference,
   so the address cannot be recycled underneath it. The grid cache was already
   testing its own hit with a full structural compare of the table, which cost as
   much per point as it saved.
4. **An insertion sort over already-sorted input is the worst case, not the
   best.** `grid_insert` scanned from the front, and every table built by
   `Table[{x, f[x]}, {x, a, b, dx}]` arrives ascending.

## A representation decision made about one value is paid for by another (2026-07-31)

The third HPC sweep produced seven fixes and six were the same defect: an
operation had a working buffer path and a working List path, and quietly took the
second whenever the two representations met. `PACK_MIN_ELEMENTS` is 250 — a
32-element vector is *correctly* left unpacked — but `X . w` for a 20000×40 `X`
then ran the symbolic loop, 320 ms against 0.31 ms.

Rules:

1. **At a binary operation, pack the small operand up; never materialise the
   large one down.** The threshold judges a value in isolation. A binary
   operation's cost is set by its largest operand, so the isolated judgement is
   the wrong one. Packing is value-preserving by contract, so there is nothing to
   weigh: converting 40 doubles to save 800,000 symbolic multiplies is free.
2. **A fast path that is not opted in does not exist.** Twice in one day I wrote
   a correct fast path, verified it against the List path, and measured no
   change — once because the head was not in `pack.c`'s `AWARE` list so the gate
   had already materialised the arguments, once because the operand needing
   conversion was below the threshold so the function returned NULL before
   reaching the new code. Both times the *values* were right, which is what made
   it silent. **Measure the thing you just changed before believing you changed
   it**, and if the timing does not move, the code is not running.
3. **Check the OUTPUT of a hot function, not just its input.** The single largest
   win was a return statement: a step function took 42 ms on a packed argument
   and 5.75 s on its own output, because building `{a, b, c}` out of packed rows
   materialised every one. The cost appeared in the *next* iteration, so
   profiling the slow call pointed at `Outer` and the arithmetic — all of which
   were innocent and had already been fixed.
4. **Probe the other system before writing a cross-system benchmark.** Two of my
   eight kernels were wrong on the first draft because I assumed
   `matrix - vector` broadcasts. It does not, in either language — `Plus` is
   Listable and threads the *outer* level. One `wolframscript` probe would have
   cost a minute; the wrong assumption cost two rewrites.
5. **A test asserting a fallback is not asserting an invariant.** The `Total[m,
   {k}]` test read `... == Total[list, {k}]` and passed because the spec fell
   through to the List path. When the spec grew a real implementation the test
   failed — correctly — and the right response was to rewrite it to compare
   values across every level range, not to restore the fallback. Comments that
   describe *why the slow path is taken* rot the moment it stops being taken;
   two `packed-arrays.md` claims and two tests had to be rewritten for the same
   reason, and the tests still passed because they used sizes below the
   threshold.
6. **A third measurement separates "slower than a competitor" from "slower than
   the machine".** Adding NumPy — which links the same Accelerate BLAS — moved
   three rows from acceptable to obviously wrong. The sieve is 1.20× *ahead* of
   Mathematica and 26.9× behind NumPy; nothing in a two-system comparison could
   have said so.

## Fourth HPC sweep (2026-07-31)

- **A passing test can be describing the wrong behaviour.** Three ways it
  happened in this codebase, all of which passed for months:
  (a) the assertion compared numeric *distance* where the defect was in element
  *heads* (`Clip` returning `{-1., 0., 1.}` for the List path's `{-1, 0., 1}`);
  (b) the assertion described a *fallback* as though it were an invariant
  (`Total[m, {k}]` "degrades to the List path"); (c) the test data was **below
  `PACK_MIN_ELEMENTS`**, so the path under test never ran. Before trusting a
  test that guards a fast path, check that the path actually executes and that
  the comparison is sensitive to the thing that can break.

- **Compare two spellings of the same operation, not just one against a
  reference.** `FoldList[Plus, 0., NDArray[…]]` and
  `FoldList[Function[{p,q}, p+q], 0., NDArray[…]]` are the same computation;
  after a change they returned different *heads*. No single-path test sees that.

- **An O(1) operation that costs O(n) will not look wrong in a profile** — it
  looks like the function you called. `First[v]` at 123 ms on 10⁶ elements only
  became visible next to `Drop[v, 250]` at 0.88 ms on the same data. When
  probing, always include an operation that *should* cost the same, as a control.

- **An integer division in an inner loop can hide a vectorisable loop behind
  it.** The convolution's affine-stride fast path measured *worse* than the step
  before it, because `(o / stride) % dim` per output was swamping everything. The
  strides are runtime values, so the compiler cannot strength-reduce them; use
  an odometer.

- **A dtype choke point is the right default and the wrong inner loop.**
  `ndt_get`/`ndt_set` is the one place that knows every dtype, which is why it
  should be the fallback — and it is two indirect calls per element where a
  float64 arm is one instruction. `Accumulate` and `Differences` were 4× and 16×
  off NumPy for this reason alone.

- **A loop whose trip count is a runtime 1 still pays a full loop.**
  `Accumulate`'s rank-1 case ran the generic `T`-inner loop with `T == 1`: the
  prologue, the index multiplies and the trip test, all around a single add.
  Specialising it was 3×.

- **Measure the benchmark harness, not just the benchmark.** One row (Jacobi)
  read 223 ms in a full run and 128 ms alone. Bisecting the *prefix* found the
  cause — a preceding `dgemm` leaves Accelerate's threads competing with ours,
  reproducibly, at 1.45×. That is a real finding about the system, and it would
  have been written off as noise.

## Round 8 — two representations, opposite gates (2026-08-01)

- **A safety guard on one representation can be the reason another is
  unguarded.** The packing transparency gate materialises an `int64` buffer for
  any head that has not claimed exactness on one, which is why `Sin` of a
  *packed* integer list has always been right. A **visible** `NDArray` is
  deliberately not gated — naming the head is the request for the buffer — and
  nothing else stood between an integer buffer and a `double` kernel. So every
  `real_closed` kernel truncated: `Sin[NDArray[{1,2,3}, DataType -> "int64"]]`
  was `{0, 0, 0}` and `Exp` was `{2, 7, 20}`. When a guard exists, ask which
  paths *bypass* it, not just whether it works.

- **Timing one operation at a time attributes the cost to the consumer and
  never names the producer.** `Union` at 807 ms, `Tally` at 68 ms and
  `DeleteDuplicates` at 67 ms were all one defect: `Mod` dropped the buffer, so
  they were handed 10⁶ boxed Integers. All three were already on `AWARE` *and*
  `INT64_OK` and were fast the moment the data arrived packed. Packing is a
  chain — audit whether a head **returns** a packed array (`--survival`), not
  only whether it consumes one quickly.

- **`ToNDArray` does not make an `NDArray`.** It returns a *packed List*
  (`Head` is `List`); `NDArray[...]` is the visible constructor. Building the
  audit's third surface out of `ToNDArray` would have made it a force-repack of
  the second and left the surface with no gate in front of it untested — the
  surface where the wrong answers were. Check what a constructor actually
  returns before making it the basis of a comparison.

- **A checksum can disagree with itself across representations.** The sweep's
  `ck` is `Total[Flatten[{x}]]`, and `Flatten` descends only into `List` heads —
  so on a visible `NDArray` it does not reduce, `Total` stays unevaluated, and
  *every correct row* read as a value mismatch. Before believing a wall of
  DISAGREE rows, verify the comparison on a case known to be right.

- **A test that pins the buffer path to the buffer path proves nothing.**
  `test_ndarray_functions.c` asserted
  `Quotient[NDArray[{1.,2.,3.,4.,5.}], 3] == NDArray[{0.0, 0.0, 1.0, 1.0, 1.0}]`
  — Reals. The List path has always given the exact `{0, 0, 1, 1, 1}`, and so
  does Mathematica. The assertion enshrined the defect it was meant to guard,
  because it was written against the implementation rather than against the
  other path. Assert buffer output equals **List** output.

- **`is_packed_list` where `is_ndarray` was meant is a silent 100× cliff.**
  `setop_i64` tested the packed form only, so the whole rank-1 int64 set-op path
  was unreachable from a visible `NDArray` (`Union` 850 ms against 5.85 ms). The
  two predicates read almost identically and mean quite different things; a
  *dispatch* gate almost always wants `is_ndarray`, and only a *presentation*
  decision wants `is_packed_list`.

- **The same asymmetry surfaced three times, in three unrelated subsystems.**
  The kernel engine truncated a visible int64 array; `setop_i64` could not see
  one; and `MemberQ`/`Count`/`Position`/`Cases`/`FreeQ` searched one for
  elements it does not expose and answered `False`/`0`/`{}` — a *confident wrong
  answer*. In all three the packed form was fine, because the transparency gate
  had already materialised it. **Nothing in `patterns.c` mentions arrays at
  all**, so no amount of reading dispatch sites could have found the third. When
  a guard protects one representation, enumerate what reaches the code *without*
  passing through it.

- **"Answered nowhere" is not "answered differently."** The audit's first
  classifier reported 29 rows as `DISAGREE` that were simply undefined
  functions, because `numeric_sweep.agree()` returns False for `"UNEVAL"`
  against itself — correct by its own contract, wrong as a building block. A
  wall of red findings is itself a signal to re-check the classifier before
  chasing any of them.

- **A tool that over-reports is worth fixing before its output is acted on.**
  The survival check flagged `Extract` (3 elements), `TakeLargest` (10) and
  `Cross` (3) as producers that "dropped packing" — but `ToNDArray` deliberately
  ignores `PACK_MIN_ELEMENTS`, so a short result reads as packable when not
  packing it is correct. Six of 44 findings were noise until the threshold was
  applied.

## Ninth round — buffer paths for the 26 heads that declined a visible NDArray (2026-08-02)

- **"Make it evaluate" is not "make it fast", and an audit row cannot tell them
  apart.** The plan inherited from the eighth round was a post-gate in `eval.c`
  that materialises a visible array when a node comes to rest. It would have
  cleared all 26 `ND-UNSUPPORTED` rows and left every one of them materialising
  10⁶ `Expr` nodes — the audit would have gone green while nothing got faster.
  The user's correction ("the fix should be highly efficient, fast paths for
  NDArray objects") is the general rule: **when a metric can be satisfied by
  making the symptom disappear, satisfying it is not evidence the cause is
  gone.** `check_packed_aware.py` and `numeric_coverage.py` are already
  documented as failing this way — registration is not speed — and a post-gate
  is the same mistake with evaluation in place of registration.

- **Count the heads that share a missing capability before writing head-level
  code.** `ArcTan[a, b]` and `Beta[p, q]` looked like two separate gaps. Both
  were `ndarray_map_binary` requiring exactly ONE array operand, and so were
  thirteen other registered-but-unreachable kernels. One function
  (`ndarray_map_binary2`) fixed fifteen heads. The tell was that the kernels
  were already *registered*: a registered kernel that never fires is an engine
  question, not a head question.

- **Marking a head `AWARE` changes what reaches its SLOW path too.** The gate
  stops materialising, so a buffer now arrives at code that tests
  `type == EXPR_FUNCTION` and falls out unevaluated. `Inner` did exactly that
  for any operator pair but `Times`/`Plus`; `Prime`, `PowerMod`,
  `IntegerDigits` and `HypergeometricPFQ` each needed an explicit
  `ndarray_delist_and_reeval` degrade. **Adding a head to `AWARE` is a promise
  about every input it can receive, not only the ones the new fast path
  handles.**

- **A checksum comparison validates the algorithm; only `===` validates the
  kernel.** `ndk_ArcTan2_c` computed `arg` via `csqrt`/`clog` where the scalar
  builtin calls `atan2`, and had been 1–2 ulp out on 68 of 400 elements since it
  was written. The sweep's tolerance is 1e-5 — correctly, since it compares
  three systems that print differently — so it was never going to report it. A
  kernel and its scalar twin are the *same* system and should be held to
  element-wise identity.

- **A negative performance result is worth recording precisely.** A float64 hot
  lane was added to both binary map chunks by analogy with `ndu_hot_chunk`,
  which is most of the unary elementary set's win. Measured: **8% on `ArcTan`,
  nothing on `Beta`.** Two-argument kernels are libm-bound, not marshalling-
  bound. Writing "8% and nothing" into the comment stops the next round
  re-deriving the same disappointment.

- **Checking whether a finding is *right* is worth as much as fixing it.** The
  eighth round listed the 6×6 `Dot`/`Inverse`/`LinearSolve` rows among the 26.
  Re-probed, all three answered correctly and always had — and running down why
  the row existed found `Flatten` treating a visible `NDArray` inside a `List`
  as an atom, a fourth instance of the same asymmetry in a fourth subsystem. The
  bug was in the *checksum*, which nobody would have audited.

- **The RNG is part of the answer.** `RandomSample`/`RandomChoice` gather from
  the buffer through the *same* `fisher_yates_sample` / `random_index` calls the
  List path makes, not a new shuffle. A separate implementation would have made
  `SeedRandom[1]; RandomSample[v]` depend on whether `v` happened to be packed —
  a surface disagreement with no error message and no way to see it coming.

- **`chk_eq` compares against a literal string, not an evaluated expression.**
  Half the first draft of the new tests asserted `f[NDArray[...]]` against
  `"f[{...}]"` and compared an `NDArray[...]` printout to that text. Use
  `Normal[lhs] === rhs` → `"True"`, and remember the printer writes `2.0`, not
  `2.`.

- **A control row is what turns a measurement into a diagnosis.** `GCD` and
  `LCM` measured 2.6×/3.8× slower packed than visible and could have been
  written down as "packed set-op weirdness". Timing `Mod` and `Quotient`
  alongside them — same kernel machinery, same shape, **not** `Orderless` —
  gave 2.19 ms on both surfaces and named the cause in one step:
  `expr_compare` materialises a packed list to compare it against a scalar, and
  `Orderless` is what makes that comparison happen. Pick the control by the
  attribute you suspect, not by size.

- **A timed expression that constructs its own operands measures the
  constructor.** The first measurement pass here built the second `NDArray`
  operand *inside* the timed call for `ArcTan[v, w]` and `Beta[p, p]`:
  `t["ArcTan visible", ArcTan[vn, NDArray[Reverse[v], …]]]`. A ~90 ms
  conversion in front of a ~2.5 ms kernel does not add noise, it **is** the
  measurement. Reported 3.6× and 2.7× where the truth is 112× and 222×.
  Worse than the wrong numbers was the wrong *conclusion* they supported —
  "these two are libm-bound, marshalling is not the cost" — which then made a
  float64 hot lane look worthless (8% and nothing) when a real A/B gives 23%
  and 6%. **Hoist every conversion into a preamble**, which is exactly what
  `nd_surface_audit.py` does and why its numbers disagreed with mine by two
  orders of magnitude.

- **When your harness and the project's tool disagree, the tool is right until
  proven otherwise.** The disagreement was visible for hours in
  `AbsoluteTiming` numbers I had already written into three documents. Treating
  a 30× discrepancy as "different data sizes" instead of stopping to reconcile
  it is how a wrong number reaches a changelog.

- **Prefer the existing tool's output to a hand-rolled table.** The audit
  computes plain/packed/visible, gain and skew per probe with its arrays built
  in a preamble. Every hand-rolled table is a second harness that has to be got
  right again, and this one was not.

- **A safety margin is a claim about a trade, and a trade has two measurable
  sides.** `PACK_MIN_ELEMENTS = 250` carried its own reasoning in a comment —
  chosen for "blast radius, not cost", with break-even already known to be
  around n = 2. The reasoning was sound and the constant still cost **544×** on
  `Det` of a 6×6, because neither side had been measured: not what the margin
  bought (swept: no regression anywhere down to 4) nor what it cost (a complete
  LAPACK path unreachable). When a constant is defended by an argument rather
  than a number, that is the constant to sweep.

- **The fast path for "this comparison cannot depend on the contents" should
  still call the real comparison.** Fixing `expr_compare`'s packed-vs-scalar
  materialise, the short version was `return 1`. The version that shipped builds
  an empty `List` stand-in and calls `expr_compare` on it — same tier, same
  answer, but the ordering rule stays in exactly one place. A hardcoded verdict
  in a fast path is a second copy of a rule that something else owns.

- **Lowering a threshold does not add behaviour, it makes existing behaviour
  reachable — including the parts nobody looked at.** `PACK_MIN_ELEMENTS`
  250 → 4 put a 2×2 on the LAPACK path, and LAPACK's `-0.0` for a
  subtraction-reached zero had been there all along at ≥250 elements. Nobody
  had seen it because a matrix small enough to read had never been packed. The
  changed sizes are exactly the sizes a human inspects, so a latent cosmetic
  wart became a visible one. **Before changing a size gate, ask what code the
  new sizes will now reach, not just what the new sizes will cost.**

- **Five test failures from a threshold change were all test-side assumptions,
  and that was still the useful signal.** Every one asserted
  `type == EXPR_FUNCTION` on a result that had merely started packing — no value
  changed. But it is exactly the "blast radius" the old constant was defending:
  code that assumes a numeric List is an `EXPR_FUNCTION`. Product code is
  protected by the transparency gate; C tests bypass it. Hence `test_delist()`
  in `test_utils.h`, and `pack_set_min_elements()` in the one test whose subject
  IS the threshold — a test about what happens below a boundary must name the
  boundary rather than inherit it.

- **A timing tool sharing a machine with anything else is measuring the other
  thing too.** The final audit ran in the same job as a
  `check-array-exactness` pass and reported four `ND-SLOW` rows. The tell was
  internal inconsistency: `nonnegative` tripped at 0.59× while `positive` and
  `negative` — the same operation on the same data — sat at 0.74/0.73. Idle,
  `nonnegative` costs 24.6 ms against the 65–110 ms the audit recorded, and all
  four probes are level. Same lesson as the hand-rolled harness earlier in the
  round, one level up: **check a surprising row against a control that should
  behave identically before writing it down.**

## 2026-08-02 — the sixth sweep: auditing the Compile surface, and a tool that names nobody

- **An audit that reads a curated list can only find what someone already
  thought of.** Four tools guarded the packed surface and all four were green
  on the day `Commonest` cost 880 ms where `Tally` of the identical buffer cost
  21.5 ms. Two read the *source* (so a head with no dispatch is invisible), one
  read *element heads* (Commonest's were right), and one *measured* but from a
  hand-written probe list with no `commonest` entry. The fix is not a better
  list — it is a tool that enumerates the whole symbol table and **discovers
  each head's call shapes by trial**, so nothing depends on remembering a name.

- **One signal is never enough to separate "slow" from "not on the buffer".**
  Cost per element cannot tell a missing fast path from an expensive function;
  the packed/unpacked ratio cannot tell one from a head that never touches
  elements. Both together are decisive. Design the discriminator before the
  harness — otherwise you get a ranked list you cannot act on.

- **A wrong answer can need the head to be nested before it shows.**
  `Compile[{{v,_Real,1}}, Max[v]]` was *rejected*, so nothing looked wrong;
  `Max[v] + 1.` compiled to `v + 1.` and returned a list where the interpreter
  returns a number. Likewise `Floor[v]` alone was right and `Total[Floor[v]]`
  returned `2.96439*10^-323`. **Probing a head in isolation is not probing the
  head.** Every new coverage probe should have a nested form beside the bare
  one.

- **A declared type is a promise to the CONSUMER, not a description.** The
  narrowing kernels declared Real and wrote `NDT_INT64`. Standalone that is
  invisible, because the caller re-reads the buffer's real dtype. Inside one
  program the next opcode reads the slot as declared and gets the integer's
  bits as a double. Wherever two layers each decide a type, one must be derived
  from the other — `nd_unary_elem` now mirrors `ndarray_map_unary` condition
  for condition.

- **A missing arm is not one fact.** `!k->cplx` meant "degrade sentinel" for
  years. It now means three different things (sentinel / narrowing-only /
  integer-only), and a guard keyed on the wrong field re-breaks a class the
  previous fix just closed: keying on `to_int_r` instead of `to_int` left
  `MoebiusMu` uncompilable *after* the narrowing case was fixed. When a
  sentinel gains meanings, grep every test of it, not just the one you came for.

- **A gate that fails from the day it lands is not a gate.** 52 heads are
  legitimately still unlowered. Shipping `check-compile-coverage` as
  "non-empty ⇒ fail" would have made it noise within a week. It ratchets
  against a checked-in `BASELINE` instead: a NEW gap fails, the backlog is
  reported, and a head that starts compiling is named so the line gets deleted.

- **The sweep's own hazards were the expensive part.** `PadLeft[vi, wi]` reads
  its second argument as a dimension spec — two 200000-element integer vectors
  ask for `prod(wi)` elements — and hung the first full run for half an hour
  behind 90-second single-probe retries. Enumerating dangerous heads by name
  missed `ParametricPlot3D`; matching by suffix does not. **For an exhaustive
  prober, the skip table and the timeout budget are the design, not the
  boilerplate.**

## 2026-08-02 (later) — what the sweep's own bugs taught

- **A `Listable` head accepts everything, so "the head changed" is not a test.**
  The first gate run reported half the symbol table, because `BesselJ[v]`
  threads into a List of 50000 *unevaluated* `BesselJ` calls and the result's
  head is `List`. The test has to be `FreeQ[result, H]` — the head must not
  appear anywhere in the answer. Before trusting a probe's notion of "this
  worked", ask what the answer looks like when it *didn't*.

- **`symtab_add_builtin` is not the symbol table.** 313 live names have no C
  registration — everything in `src/internal/*.m`, plus options and colours. A
  source scan cannot even name them. Ask the binary (`Names["*"]`); the
  discovery step drops the ones that are not functions for free, which is
  cheaper than a hand-maintained exclusion list and cannot go stale.

- **A per-process counter needs a per-process question.** The gate diagnostic
  totals by head over the whole run, so gating every shape at once mixed the
  accepted calls with the unevaluated ones and the element count stopped
  belonging to the call being reported. Discover first, then gate only the
  shapes that mean something.

- **A cap that changes the answer is a wrong answer, not a limit.**
  `Length[Range[2000000]]` was `1000001`, silently. Found by accident while
  sizing a benchmark vector — the only way a silent truncation is ever found is
  downstream of the damage. Two rules follow: put the ceiling on the resource
  that actually runs out (bytes, not elements), and *decline* rather than
  truncate, so the caller has something to test for. Worth grepping for other
  `if (count > CAP) count = CAP;`.

- **Marking a head packed-aware is a commitment to its REWRITE TARGET's shapes
  too.** `Subtract` and `Divide` rewrite to `Plus`/`Times`/`Power`, so making
  them aware would have inherited those heads' handling of a symbolic operand —
  which was to leave `packedList + x` unevaluated. The bug had to be fixed at
  the source rather than propagated to two more heads.

- **A benchmark helper that is not `HoldAll` measures nothing.**
  `best[e_] := ... Do[t = Min[t, First[AbsoluteTiming[r = e]]], {5}]` evaluates
  `e` once on the way in and then times an assignment: every row read 0.001 ms.
  Generate the timing code textually, with a DIFFERENT constant per repetition
  so no cache can answer it.

- **Measure the "before", do not infer it.** The gate diagnostic said `Subtract`
  materialised 200000 elements, and the obvious conclusion — `v - 1.` is slow —
  was wrong: the parser desugars infix `-` to `Plus`, so only the *explicit*
  `Subtract[v, 1.]` was ever slow. A `git stash` of `src/`, a rebuild and one
  measurement is four minutes and turns a plausible number into a real one.

- **A diagnostic can name something you never probed.** The first
  baseline-validation run of `nd_fastpath_sweep.py` died with
  `KeyError: 'List'`: an internal `List[...]` assembling a result materialises
  a buffer too, and `MATHILDA_PACK_DIAG=gate` records it by name like any other
  head — so the verify phase went looking for a probe expression that never
  existed. When you join a diagnostic's output back against your own input,
  the diagnostic's key set is not a subset of yours. Filter, do not index.

- **A zero-test that is sound in one direction is not a pivot test.**
  `is_zero_poly` proves zero by polynomial identity: when it says *zero* it is
  right, when it says *nonzero* it may just have failed. `RowReduce` consumes
  precisely the unreliable answer — it needs "this entry is nonzero, pivot on
  it" — so for a *casus irreducibilis* eigenvalue it row-reduced a singular
  matrix to the **identity** and `Eigenvectors` returned three zero vectors.
  Before wiring a partial decision procedure into a caller, ask which of its
  two answers the caller actually acts on.

- **Prefer removing the need for a hard test to strengthening it.** The
  tempting fix was a stronger algebraic zero test (qqbar) inside `RowReduce` —
  broad blast radius, and slow. The adjugate identity
  `M · adj(M) = det(M) · I ≡ 0 (mod q)` yields the same eigenvector out of
  cofactor determinants: no division, no pivoting, so nothing has to be decided
  zero *while it is computed*. The one test left runs after reduction mod the
  minimal polynomial, on a univariate polynomial over `Q`, where the same
  `is_zero_poly` is exact and complete. Move the test to where it is decidable.

- **A test the defining identity cannot fail is not a test.** `m.v == λv` holds
  trivially for `v = 0`, which was exactly the bug — so the regression test has
  to assert *nonzero* first. Same shape as the SVD lesson that reconstruction
  cannot see a wrong basis for a zero singular value: when the wrong answer is a
  zero, every identity that multiplies by it passes.

- **Check the make target exists before believing "tests pass".**
  `make -j8 test_eigen` in `tests/build` prints `No rule to make target` — which
  contains neither "error" nor "Error", so a `grep -E "error|Error"` filter
  reports a clean build and the stale binaries from the last session run
  happily and pass. The targets are `eigen_tests`, not `test_eigen`. Grep the
  build for `No rule` too, or confirm the new test's name appears in the run.

- **A constant that routes around a blow-up is a symptom, not the bug.** Issue
  #41: `Limit[c ArcTan[Sqrt[-1+x]/Sqrt[2]], x->Infinity]` hung while the bare
  `ArcTan[...]` worked. Tempting fix: pull the constant out of the limit. But
  the bare form only worked because `compose_at_infinity` caught it *before* the
  Series layer; `Series[ArcTan[Sqrt[-1+x]/Sqrt[2]], {x,oo,2}]` hangs on its own,
  no constant. The constant merely denied the fast path. Always test whether the
  "working" sibling works for the right reason before declaring the difference
  the cause. The real bug was general (all bounded kernels over shifted radicals
  at infinity) and lived in the shared series machinery, not the Limit layer.

- **Horner series composition nests radical coefficients because Times does not
  distribute over Plus.** `so_compose_scalar_kernel` builds `Σ aₖ uᵏ` by
  repeated `result = result*u + aₖ`; each step's coefficient becomes
  `scalar + u_coef*(prev)`, and `simp`/evaluate leaves `Sqrt[2]*(a + b Sqrt[2])`
  un-multiplied, so depth-N composition nests exponentially. `Expand` performs
  exactly the missing distribution + like-radical collection
  (`Sqrt[2]*(a+b Sqrt[2]) -> a Sqrt[2] + 2b`). Fix: normalize radical-bearing
  coefficients between Horner steps. General across every at-zero kernel.

- **Extracting a square factor from `c^2 * r` must not trial-divide `c`.** The
  squarefree part is invariant under multiplication by a perfect square, so
  `squarefree(c^2 * r) == squarefree(r)`: factor only the small radicand, fold
  `c` into the extracted root, and recover the coefficient<->radicand
  cross-cancellation with gcds. Folding `c^2` in and re-factoring is O(c) and
  hangs on large series coefficients. Prove equivalence to the old form
  exhaustively (a Python model over millions of small inputs) before porting.

- **A node-cached memoization index survives only if the node is stable, and the
  eval clock decides that.** Caching a key→position index on an Association node
  (in `EXPR_FUNCTION` union slack, so `sizeof(Expr)` is unchanged) makes the
  lookup PRIMITIVE O(1), but the interpreter re-*builds* function nodes from
  evaluated args every `evaluate_step` (`eval.c:1197`), so an eagerly-attached
  index lands on the node the fixed-point logic then DISCARDS (it keeps the
  original `current`, frees the structurally-equal rebuilt `next`). Two things
  fix it: (1) attach the index LAZILY in the reader, on the node that actually
  survives; (2) add an in-loop timestamp short-circuit in `evaluate()` (mirroring
  the entry check at `eval.c:1939`) so a value already stamped under the current
  clock is not re-canonicalised each time it is reached mid-loop. Verify O(1) on
  the PRIMITIVE (direct C calls), never through a loop — see next.

- **`Do`/`Table`/`Fold` bump the global eval clock every iteration; `Map` does
  not.** Iterator binding goes through `symtab_add_own_value` (`iter.c:417,442`),
  a symbol-table mutation that bumps the clock, invalidating ALL `last_evaluated_at`
  memoization each step — so a large loop-invariant value (an association) is
  re-evaluated O(n) per iteration regardless of any node-cached index. `Map` with
  a pure function binds no named symbol, so the clock is stable across its
  elements and the memoization/index holds. Consequence: benchmark an "O(1)
  lookup" claim on the primitive or on `Map`, not on `Do[Lookup[a,k],{...}]`
  (which measures O(n) re-canonicalisation, not the lookup). A `Do`-loop timing
  that looks O(n) after an "O(1)" change is this, not a broken index.

- **Immutable-by-convention is the precondition for any node-attached cache.**
  Before caching anything on a shared expression node, audit that nothing mutates
  it in place after construction (`grep` every `data.function.args[i] =` and
  `expr_unshare` site). Associations passed the audit — every update rebuilds via
  `assoc_from_rules`/`assoc_entry_with_value` — except one latent aliasing bug
  (`part.c` `delete_path` wrote a refcount-shared entry's value in place). The
  cache was safe; the bug was orthogonal but real. `expr_copy` is a refcount
  bump, so a "copy" shares the node and its cache; only a physical copy
  (`expr_unshare`) must null the cache pointer to avoid a double free.

## Verify that your verification tool actually verifies (2026-08-05)

`site/verify_tutorial.py` piped a tutorial's `In[]` lines into `./Mathilda` and
parsed the stdout for `Out[]=` lines. But Mathilda serves a **non-tty stdin over
the NDJSON pipe protocol** (`src/repl.c pipe_mode_loop`), not the interactive
`In[]/Out[]` transcript — so the tool matched **zero** `Out[]` lines,
`zip(pairs, [])` iterated nothing, and it printed `OK` for **every** tutorial,
forever. "Verified against the binary" was aspirational, not enforced.

- **Rule**: a checker that can pass with an empty result set is not a checker.
  Before trusting one, feed it a deliberately-wrong expectation and confirm it
  FAILS. (One `MISMATCH` from the fixed tool is worth more than a hundred green
  runs from the broken one.)
- **How to drive Mathilda headlessly**: send `{"id": k, "expr": "..."}` lines on
  stdin, read `{"id": k, "type": "expr", "payload": ...}` back; the payload is
  the exact `OutputForm` a REPL user sees as `Out[k]=`. Match by id, not by
  position, so a `;`-suppressed line (payload `"Null"`) and raw `Print` /
  `CompilePrint` stdout can't shift the alignment.
- **Tutorial `In[k]:=` expressions must be single-line.** The pipe protocol (and
  the verifier) send one line per expression; a wrapped `Compile[...]` across two
  markdown lines truncates to the first line, and the definition silently never
  binds. Caught only because the fixed verifier flagged `logistic[0.5, 20]`
  coming back unevaluated.

## A packed-int64 buffer materialises for any aware head lacking `packed_int64_ok` (2026-08-05)

`N[Range[10^6]]` unpacked to a list of boxed reals while `N[Range[1., 10^6]]`
stayed packed. Root cause is NOT in `N`: `eval.c`'s transparency gate
materialises an `int64` packed buffer to a nested `List` for any packed-aware
head that has not claimed `packed_int64_ok` (the `Sin[int64]={0,0,0}` guard
class), so `N` received a plain list and `numericalize`'s element rebuild dropped
packing. Real buffers skip the gate branch and pass through untouched — hence the
asymmetry between the same operation on int vs real data. Fix = claim
`packed_int64_ok` (so the buffer reaches the builtin) AND make the builtin
actually handle the buffer (`numericalize_rec`'s `EXPR_NDARRAY` case widens
int64→float64 via `expr_new_ndarray_like`, inheriting presentation). One without
the other is incomplete: the claim alone would hand `N` an int64 array it copies
verbatim (wrong: int, not real); the widening alone is never reached.

## A reported build warning can be the visible tip of a chronically-broken degrade config (2026-08-06)

A handful of `-Wunused-function` warnings in `rat.c`/`rootreduce.c` (plus a
failing Linux CI email) turned out to sit on top of a no-MPFR/no-FLINT config
that had been broken far past the reported symptoms: `groebner.h` +
`nsolve_system.c` included `<mpfr.h>` unguarded (the actual CI compile failure),
`nd_ac_prec_free` was defined inside a file-wide `#ifdef USE_MPFR` yet called
unconditionally (a *link* error), and **17** static functions were dead code in
the degrade config. The CI `build-no-mpfr` job had only ever reached the compile
stage — it aborted early on the mpfr.h include — so its green-until-now history
never covered the link or the tail of the file list. Lessons: (1) when a
degrade-config warning appears, fix the whole *class*, not the two files that
happened to be reported, and add a gate (`-Werror=unused-function`) so it can't
silently return; (2) macOS masks this entire class because `mpfr.h` is on the
include path — reproduce the Linux build locally with a **poison stub `mpfr.h`**
(`#error`) injected via `CC="gcc-NN -I<dir>"` ahead of `-I/usr/local/include`,
and `make -k` to enumerate the whole backlog in one pass rather than iterating
through CI one failure at a time. See [[project_use_mpfr_zero_build]].

## A failing test can be right and a deliberate code change wrong (2026-08-07)

`test_quotient` had failed since 2026-07-27: it expected
`Quotient[17.5 + 6 I, 1 + 2 I] == Complex[6, -6]` but the code returned
`Complex[5, -6]`. The tempting read is "stale test, update the expectation to
match the code" — especially because the code change that broke it (commit
`c0c8dcb`) came with a *detailed justification comment*: complex `Quotient` was
switched from `round` to `floor` to "agree with the real branch", on the premise
that `Quotient` is `Floor[(m-d)/n]` by definition.

The premise was false. That is the definition for **real** arguments only; for
complex arguments `Quotient` is Gaussian-integer division, which rounds each
part of the ratio to the nearest integer. One `wolframscript` batch settled it —
and *every* example the commit cited as support was wrong against Mathematica
(`Quotient[5.5 + 1. I, 3.]` is `2` not `1`, `Quotient[10 + 2 I, 3 + I]` is `3`
not `3 - I`). The test had been correct all along.

Rules: (1) when a test and the code disagree, **the test is a hypothesis about
ground truth, not automatically the stale side** — check the authority (here the
local Mathematica kernel) before deciding which to change; the direction of the
fix is the whole question. (2) A confident justification comment is not
evidence; it is exactly what a `wolframscript` probe is for — the more detailed
the rationale, the more it is worth the one call to check. Same
`/Applications/Mathematica.app/Contents/MacOS/wolframscript` pattern as the
`DiagonalMatrix` lesson above. (3) The Mathematica-faithful "nearest" is
round-half-to-**even** (`Quotient[5 + 3 I, 2] == 2 + 2 I`, ratio `2.5 + 1.5 I`),
so the fix reused a `round_half_even` helper, not C `round()` (ties away from
zero).
