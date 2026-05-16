# 6. Performance Notes

## 6.1. Polynomial Algebra (`poly.c`)

The polynomial subsystem -- `Coefficient`, `CoefficientList`, `PolynomialGCD`,
`PolynomialQuotient`, `PolynomialRemainder`, `Resultant`, `HornerForm`,
`Discriminant`, `Decompose`, `Collect`, etc. -- shares a small set of inner
loops. Two optimisations dominate:

1. **Direct coefficient extraction.** The internal helper
   `get_coeff_expanded(expanded, var, n)` walks an already-expanded
   polynomial once per query, summing the contributions of each summand
   directly. This avoids the full evaluator pipeline (which would
   construct `Coefficient[expanded, var, n]`, look up its symbol,
   apply `Listable`/`Flat`/`Orderless`, then re-expand inside
   `builtin_coefficient`). A bulk variant `get_all_coeffs_expanded`
   produces every coefficient in a single pass and is used by
   `poly_content`, `coeff_list_rec`, `poly_derivative` and the
   resultant routine.

2. **Skipping `Together` + `Cancel` in pure-integer division.** The
   univariate division loops in `exact_poly_div` and `poly_div_rem`
   used to call `internal_together` and `internal_cancel` on every
   iteration to unify denominators introduced by symbolic leading
   coefficients. We now track whether the quotient coefficient was
   produced by an exact integer / bigint division (in which case no
   new fractions can appear) and only invoke the heavy unification
   path for the genuinely symbolic case.

Combined effect on the test bench (compiled at `-O3`):

| Test                             | Before | After  | Speedup |
|----------------------------------|-------:|-------:|--------:|
| `tests/poly_tests`               | 1.86s  | 0.53s  |   3.5x  |
| `tests/facpoly_tests`            | 4.24s  | 1.60s  |   2.7x  |
| `PolynomialGCD` (50 large iter.) | 10.07s | 2.40s  |   4.2x  |
| `CoefficientList` (30 iter.)     | 3.81s  | 0.92s  |   4.1x  |
| `HornerForm` (200 iter.)         | 16.76s | 5.13s  |   3.3x  |
| `PolynomialQuotient` (500 iter.) | 2.55s  | 1.15s  |   2.2x  |

A small clean-up was made alongside the perf work:
* `poly.h` had stray declarations after `#endif`; they are now inside
  the include guard.
* `PolynomialQ` was registered twice in `poly_init`; the duplicate
  registration has been removed.

## 6.2. Differentiation (`deriv.c`)

`D`, `Dt` and `Derivative` used to be defined by ~60 DownValues in
`src/internal/deriv.m`. Each call walked the whole rule list linearly,
re-ran the pattern matcher (including side-conditions such as
`FreeQ[c, x]`), and re-entered the evaluator for the replacement. The
native C implementation (added in this release) replaces that with:

1. **Direct head-symbol dispatch.** A single `strcmp` cascade in
   `compute_deriv` picks the right closed-form rule for `Plus`,
   `Times`, `Power`, `Sqrt`, `Exp`, `Log`, the trig/hyperbolic
   families and their inverses.
2. **Allocation-free constant detection.** A tailored `expr_free_of`
   walk replaces the generic `FreeQ` builtin on the hot path; it
   short-circuits on the first match against the differentiation
   variable without allocating a `MatchEnv`.
3. **Direct `Derivative[...][f][g]` construction.** The chain rule
   for unknown single- and multi-argument heads builds the
   `Derivative` head, the intermediate `Derivative[...][f]` and the
   final application directly, avoiding the N^2 DownValue-search
   behaviour the matcher produced for nested heads.

Measured on the stock interpreter at `-O3`:

| Workload                                               | Rule-based | C      | Speedup |
|--------------------------------------------------------|-----------:|-------:|--------:|
| 1000 x simple mixed (`D[Sin[x^2]+Cos[a x]Exp[x]+x^x+Log[x], x]`) | 1.24s      | 0.54s  |  2.3x   |
| 1000 x deep chain (`D[Sin[Cos[Tan[ArcSin[x^2+1]]]], x]`)         | 3.37s      | 2.95s  |  1.1x   |
| 200 x higher-order (`D[Sin[a x] + x^3 Cos[b x], {x, 6}]`)        | 33.26s     | 8.60s  |  3.9x   |
| 500 x mixed partials (`D[f[x,y,z]*Sin[x+y]*Cos[x y], x, y, z]`)  | 17.78s     | 2.27s  |  7.8x   |

The largest wins come from higher-order and mixed-partial calls,
where the rule-based implementation had to re-traverse the entire
DownValue list at every inner step; the native dispatch visits the
relevant head exactly once per sub-expression.

## 6.3. Limits (`limit.c`)

Mathilda now has a native `Limit` built-in implemented in C in
`src/limit.c`, registered by `limit_init()` in the standard
`core_init()` chain. The design follows the layered dispatch outlined
in `limit_candidate_spec.md`, with each layer either resolving the
limit and short-circuiting or passing the problem down:

```
Layer 0 -- Interface normalization (calling forms + Direction option)
Layer 1 -- Structural fast paths (free-of-var, identity, continuous subst.)
Layer 3 -- Rational function P(x)/Q(x) short-cut
Layer 5.3 -- Logarithmic reduction for f^g indeterminate forms
Layer 2 -- Series-based evaluation (the workhorse)
Layer 5.1 -- L'Hospital with complexity-growth guardrail
Layer 6 -- Bounded-oscillation Interval returns
```

**Calling forms** (Mathematica-compatible):

| Form                                            | Meaning                            |
|-------------------------------------------------|------------------------------------|
| `Limit[f, x -> a]`                              | single-variable limit              |
| `Limit[f, {x1 -> a1, ..., xn -> an}]`           | iterated limit (rightmost first)   |
| `Limit[f, {x1, ..., xn} -> {a1, ..., an}]`      | multivariate joint limit           |
| `Limit[f, x -> a, Direction -> spec]`           | directional approach               |

**Direction settings:**

| Setting                        | Internal meaning                  |
|--------------------------------|-----------------------------------|
| `Reals`, `"TwoSided"`          | two-sided (default on the reals)  |
| `"FromAbove"` or `-1`          | approach from above (x -> a^+)    |
| `"FromBelow"` or `+1`          | approach from below (x -> a^-)    |
| `Complexes`                    | complex / radial-in-all-directions|

The Mathematica sign convention (`-1 == FromAbove`) is applied once in
the option parser so that the math layers only ever see internal
`+1 / -1 / 0` tags.

**Return values:** finite expression, `Infinity`, `-Infinity`,
`ComplexInfinity`, `Indeterminate`, `Interval[{lo, hi}]`, or the
original unevaluated `Limit[...]` when the system cannot determine a
value.

**Integration:**
* Registered under `Limit` with attributes
  `Protected | ReadProtected | HoldAll`. `HoldAll` prevents the
  second-argument rule from being evaluated prematurely against any
  OwnValue of the limit variable.
* Documented via `symtab_set_docstring` for `?Limit` in the REPL.
* Covered by `tests/test_limit.c`, which exercises each layer and
  every Direction setting.

**Layer highlights:**
1. *Continuous substitution* (Layer 1) does not just evaluate
   `f /. x -> a`; it Together-normalizes first and checks that the
   denominator does not vanish, so `Sin[x]/x` at x = 0 correctly skips
   the fast path instead of silently returning 0 (as Mathilda's
   arithmetic would fold `Sin[0] * 0^(-1)`). Expressions of the form
   `Power[_, expr_with_x]` are likewise refused here -- they are
   handled by Layer 5.3 instead.
2. *Series-based evaluation* (Layer 2) calls `Series[f, {x, a, k}]`
   with increasing orders (4, 8, 16, 32) until a nonzero leading term
   is found, then reads the limit off the leading exponent. Limits at
   `-Infinity` substitute x -> -y internally and recurse at +Infinity.
3. *Log reduction* (Layer 5.3) fires for any `Power[base, exp]` whose
   exponent depends on the limit variable, rewriting as
   `Exp[Limit[exp * Log[base]]]` and recursing. It sits above Series in
   the dispatch order because Series has no kernel for `b^g` with
   g a non-trivial function of x.
4. *L'Hospital* (Layer 5.1) only fires when pointwise evaluation is
   strictly `0/0` or `inf/inf`, and aborts if the leaf count of the
   quotient grows for three consecutive derivations -- this avoids
   spinning on Sin/Cos expansions where the series layer is the right
   tool.
5. *Bound analysis* (Layer 6) handles the canonical
   `Limit[Sin[1/x], x -> 0] = Interval[{-1, 1}]` shape by recognising
   a bounded trigonometric head wrapped around a divergent argument.

**Worked examples (cross-checked against Mathematica):**

| Input                                                   | Result        |
|---------------------------------------------------------|---------------|
| `Limit[Sin[x]/x, x -> 0]`                               | `1`           |
| `Limit[(Cos[x]-1)/(Exp[x^2]-1), x -> 0]`                | `-1/2`        |
| `Limit[(Tan[x]-x)/x^3, x -> 0]`                         | `1/3`         |
| `Limit[Sin[2 x]/Sin[x], x -> Pi]`                       | `-2`          |
| `Limit[(x^a - 1)/a, a -> 0]`                            | `Log[x]`      |
| `Limit[(1 + A x)^(1/x), x -> 0]`                        | `E^A`         |
| `Limit[(1 + a/x)^(b x), x -> Infinity]`                 | `E^(a b)`     |
| `Limit[x/(x+1), x -> Infinity]`                         | `1`           |
| `Limit[3 x^2/(x^2 - 2), x -> -Infinity]`                | `3`           |
| `Limit[Sqrt[x^2+a x]-Sqrt[x^2+b x], x -> Infinity]`     | `(a-b)/2`     |
| `Limit[1/x, x -> 0, Direction -> "FromAbove"]`          | `Infinity`    |
| `Limit[1/x, x -> 0, Direction -> "FromBelow"]`          | `-Infinity`   |
| `Limit[Sin[1/x], x -> 0]`                               | `Interval[{-1, 1}]` |
| `Limit[x/(x + y), {x -> 0, y -> 0}]`                    | `1` (iterated)|

**Regression-suite additions (2026-04-20):** four new layers were
bolted on in response to a batch of REPL-driven test cases. Each is
composable with the original dispatcher and sits at the order indicated
below:

1. *Reciprocal-trig normalization* (runs at the top of
   `compute_limit`, before any other layer). Rewrites
   `Csc[z] -> 1/Sin[z]`, `Sec[z] -> 1/Cos[z]`, `Cot[z] -> Cos[z]/Sin[z]`,
   plus hyperbolic twins `Csch/Sech/Coth`, and also `Tan[z] -> Sin[z]/Cos[z]`,
   `Tanh[z] -> Sinh[z]/Cosh[z]`. The rewrite is purely structural (tree
   walk) followed by a single `evaluate()` pass. It turns the
   `0 * ComplexInfinity` folding trap (e.g. `x Csc[x] -> 0`) into a
   structural `0/0` that Series and L'Hospital can resolve. This one
   change resolves roughly 15 REPL-failing cases on its own.
2. *ArcTan / ArcCot at divergent inner argument* (runs before the
   rational layer). Computes `Limit[inner, ...]`; on
   `+Infinity / -Infinity` maps `ArcTan` to `+/- Pi/2` and `ArcCot` to
   `0 / Pi`. Previously `ArcTan[x^2 - x^4]` at `Infinity` sat on a
   Series that could not expand the inner polynomial at infinity and
   emitted `Power::infy` warnings.
3. *Bounded envelope (squeeze theorem)* (gated to `+Infinity` limits,
   runs before the Series layer). Uses a structural magnitude-bound
   walk to produce `|f| <= g(x)`: `|Sin/Cos/Tanh[anything]| <= 1`,
   `|ArcTan/ArcCot[anything]| <= Pi/2`, triangle inequality on Plus,
   multiplicativity on Times, and `|a^n| = |a|^n` for non-negative
   integer `n`. If `Limit[g(x), x -> Infinity] = 0`, the original limit
   is zero. Covers `Sin[t^2]/t^2`, `(1 +/- Cos[x])/x`,
   `(x Sin[x])/(5 + x^2)`, and friends at infinity, plus `x^2 Sin[1/x]`
   at 0 (which already worked through a different path but now routes
   here without Series warnings).
4. *Generalized log-reduction* (replaces the old single-`Power` form).
   Fires for a top-level `Power[b, g]` with x in `g` OR a top-level
   `Times[P1, P2, ...]` whose non-constant factors are all of that
   shape. Rewrites as `Exp[Limit[Sum[g_i * Log[b_i]]]]`. A
   post-processor maps a `+Infinity / -Infinity` log-limit to
   `Infinity / 0`, and refuses to produce an answer for
   `ComplexInfinity` / `Indeterminate`. Covers
   `((-3+2x)/(5+2x))^(1+2x)` at `Infinity`, `x^Tan[Pi x/2]` at 1,
   `(1+Sin[ax])^Cot[bx]` at 0, and the `E^(-x) x^20` shape (via
   log-reduction of the implicit `Power[E, -x] * Power[x, 20]`).

A smaller tweak: continuous substitution now folds `Power[0, positive]`
and `Sqrt[0]` to `0` locally, so `Limit[Sqrt[x-1]/x, x -> 1]` reports
`0` instead of the un-folded `Sqrt[0]` that Mathilda's Power evaluator
leaves in place.

**Regression-suite additions (2026-04-20, batch 2):** three further
layers and one post-processing pass were added after a second round
of REPL-driven cases.

1. *Log of a finite-limit inner* (runs before Series). If `f = Log[g]`
   and `Limit[g]` is finite, returns `Log[c]`. Enables shapes such as
   `Log[1 + 2 Exp[-x]]` at `Infinity`.
2. *Log / polynomial merge at infinity* (runs after the Log-of-finite
   layer). Rewrites `Sum(Log[g_i]) + Sum(h_j)` into a single
   `Log[(Product g_i) * Exp[Sum h_j]]` whose inner product is expanded
   to expose cancellations (e.g. `Exp[x] * Exp[-x] -> 1`). Resolves
   `-x + Log[2 + E^x]` at `Infinity` to `0`.
3. *Term-wise Plus at infinity*. If every summand of a `Plus` at
   `+/- Infinity` has an individually finite limit, return the sum of
   those limits. Refuses as soon as any term is divergent or
   unresolved. Keeps the outer log-merge sound, because individually
   divergent shapes still bail and reach the other layers.
4. *Radical fusion in `Times`* (global, not Limit-specific). Moved
   out of `limit.c` and into `builtin_times` so that
   `Sqrt[6]/Sqrt[2]` becomes `Sqrt[3]` system-wide, not only as a
   Limit post-pass. The rule fuses `Power[a, q] * Power[b, -q]` into
   `Power[a/b, q]` whenever `a` and `b` are both positive numeric
   (integer, bigint, rational, or real) -- `a > 0` and `b > 0`
   places us on the principal branch. The ratio `a/b` may itself
   be an integer (`Sqrt[6]/Sqrt[2] -> Sqrt[3]`) or a rational
   (`9^(1/3)/3^(1/3) -> 3^(1/3)`, `2^(1/3)/3^(1/3) -> (2/3)^(1/3)`,
   `Sqrt[3]/Sqrt[6] -> Sqrt[1/2]`). Applies after same-base
   grouping, restarts on each fusion so chained reductions like
   `Sqrt[210]/Sqrt[6]/Sqrt[5] -> Sqrt[7]` converge.
5. *`Power[0, b]` folding* in `builtin_power`. `0^b` now evaluates to
   `0` for any positive `b` (integer, rational, or real) and to
   `ComplexInfinity` for any negative `b`. This eliminates the
   `Sqrt[0]` leak that previously bubbled up from continuous
   substitution in Limit; the auxiliary `reduces_to_zero` predicate
   has been removed from `limit.c`.

Also fixed: `heuristic_factor` in `facpoly.c` recursed indefinitely
when factoring `Power[a, rational]` where `a` was a constant (no
variables, but non-unit content). The guard is a single early return
for `v_count == 0` (no variables -> return the constant); segfault in
`Limit[1/(t Sqrt[t+1]) - 1/t, t -> 0]` was caused by this path and is
now fixed.

**Regression-suite additions (2026-04-21, batch 3):** follow-up batch
from a REPL-driven failure report. The fixes are small, targeted, and
preserve the layered dispatcher rather than adding new layers.

1. *Unary-minus extraction in the `Times` printer* (`src/print.c`). A
   leading negative literal (integer, real, or negative-numerator
   `Rational`) is now emitted as a prefixed `-` sign. This normalises
   `Times[-1, Infinity]` to `-Infinity`, `Times[-1, a, b]` to `-a b`,
   `Times[-3, x]` to `-3 x`, and `Times[Rational[-2, 3], x]` to
   `-2/3 x`. Previously these printed with a literal `-1` factor
   (`-1 Infinity`), which both looked wrong and tripped equality
   checks. Existing tests whose fixtures asserted the old spelling
   (`test_logexp`, `test_list`, `test_limit`) were updated.
2. *Opaque-head guard* (`src/limit.c`). Before running any analytic
   layer, `compute_limit` now refuses the problem if the expression
   applies an undefined or known-discontinuous head to an argument
   containing the limit variable. Unknown heads are detected via the
   symbol table (no `builtin_func`, no `down_values`, no `own_values`).
   A curated whitelist of discontinuous heads -- `Floor`, `Ceiling`,
   `Round`, `FractionalPart`, `IntegerPart`, `Sign`, `UnitStep`,
   `HeavisideTheta`, `KroneckerDelta`, `DiscreteDelta`, `Piecewise`,
   `Boole`, `Mod`, `Quotient` -- is treated as opaque even when a
   builtin exists, because plain substitution would silently pick one
   side's value at a jump. Fixes `Limit[f[x], x -> a]`,
   `Limit[Ceiling[5 - x^2], x -> 1]`, and similar cases that previously
   returned a misleading numeric answer.
3. *Even-order pole handling in the rational layer*. Layer 3 used to
   short-cut `Limit[N(x)/D(x), x -> a]` to `ComplexInfinity` whenever
   `D(a) = 0` and the direction was two-sided. This lost the parity
   information, so `1/(x-2)^2` at x = 2 came back as `ComplexInfinity`
   instead of `+Infinity`. The layer now defers unconditionally to the
   Series layer in that branch, which inspects the leading exponent:
   even-order poles return `+/-Infinity` (signed by the leading
   coefficient), odd-order poles return `ComplexInfinity`.
4. *Oscillatory limits return `Indeterminate`*. Layer 6 previously
   returned `Interval[{-1, 1}]` for shapes like `Sin[1/x]` at 0 and
   `Sin[x]` at `Infinity`. That interval is a correct bound but not a
   limit value; Mathematica returns `Indeterminate`, which we now match.
   The layer still only fires after the squeeze / substitution / Series
   paths have failed, so bounded shapes that actually converge to a
   fixed value are unaffected.
5. *Leaked limit variable in directional-infinity coefficient*. When
   the Series leading coefficient is a non-constant expression (e.g.
   `Log[x]` for `Log[x]/x` at 0), emitting
   `DirectedInfinity[Log[x]]` gave back a pseudo-answer that still
   mentioned the limit variable. `read_leading_term_limit` now bails
   in this case so the other layers (or the unevaluated fall-through)
   can run.
6. *`Direction -> I` and other purely-imaginary directions*. The
   option parser now accepts `I`, `k*I`, `-I`, and `Complex[0, k]` as
   valid directions, routing them through `LIMIT_DIR_COMPLEX`. Full
   branch-cut analysis isn't implemented yet, so such limits typically
   still return unevaluated -- but they no longer reject the option
   outright.

**Known limitations carried forward from earlier batches:**
* Joint multivariate limits with path dependence
  (`x y / (x^2 + y^2)` at `{0,0}`, `ArcTan[y/x]` at `{Infinity,
  Infinity}`, `y/(x+y)` at `{0,0}`, `(x^3+y^3)/(x^2+y^2)` at `{0,0}`):
  the polar-substitution heuristic is still future work. A handful of
  cases that happen to be continuous at the substitution point do
  resolve via `run_multivariate`.
* `(1 + Sinh[x])/Exp[x]` at `Infinity`: Series at Infinity does not
  recognise the natural `Exp[-x]` factorisation and L'Hospital stalls
  on the Cosh/Sinh derivative cycle.
* `(25^x - 5^x)/(4^x - 2^x)` at 0: Series infinite-escalates on the
  log-expansion of `k^x`; currently hangs if invoked (the REPL must be
  interrupted). Avoid until a dedicated `b^x` expansion kernel is
  added.
* `Log[1 - (Log[Exp[z]/z - 1] + Log[z])/z]/z` at a finite non-singular
  point (e.g. z = 100): one of the sub-expressions triggers repeated
  series retries inside `evaluate()`. Regression to be triaged.
* `Tan[x]` at `Pi/2, Direction -> Reals`: should return `Indeterminate`
  (the two sides diverge with opposite signs) but currently yields
  `ComplexInfinity`. A pole-with-sign-disagreement classifier is
  future work.
* `Sin[n Pi]` at `n -> Infinity`: should be `Indeterminate`; currently
  the log-reduction pathway picks an incorrect sign via the bounded
  envelope and returns `-Infinity`-like output in some code paths.
* Complex `Exp[Exp[-x/(1+Exp[-x])]]`-style stacked-exponential limits
  at `Infinity`: the bounded envelope does not see through nested
  `Exp` applied to bounded-but-decaying arguments.
* `Exp[Log[Log[x + ...]]/Log[Log[Log[Exp[x] + x + Log[x]]]]]` at
  `Infinity`: this requires comparative asymptotic expansion of
  iterated logarithms, which Series does not currently model.
* ~~The `$SeriesInvVar$` reverse-substitution symbol can leak back into
  user output for deeply nested `Log[Log[...]]` chains.~~ **Fixed
  2026-04-21**: `do_series_single` now substitutes the internal inverse
  variable `u -> 1/x` in every coefficient after the at-Infinity
  expansion, not only in `SeriesObj::x`. Previously a coefficient that
  was treated as "free of u" (e.g. `Log[1/u] -> -Log[u]` folded as a
  constant-in-u term) would retain the placeholder. Regression test
  `test_series_infinity_no_inv_var_leak` grep-checks the `$SeriesInvVar$`
  literal out of both InputForm and FullForm output for four
  representative shapes.

---

**Work-package-driven additions (2026-04-21, batches 4–12):** a
coordinated series of work packages closed out most of the residual
Limit cases from the REPL report. Each landed as its own bounded
regression in `tests/test_limit.c` (`test_wp8_*`, `test_wp2_*`,
`test_wp4_*`, `test_wp1_*`, `test_wp3_*`, `test_wp6_*`, `test_wp5_*`,
`test_wp9_*`, `test_wp7_*`).

1. **WP-8 (numeric-point fast path).** `layer1_fast_paths` gains a
   `try_numeric_point_substitution` that substitutes a plain numeric
   limit point into the Together-normalised form and returns the
   result when clean, skipping the Series / L'Hospital pipeline. Gate
   checks: Together's denominator does not vanish at the point, and no
   `Power[base, x-dep-exponent]` with a divergent exponent remains
   (1^inf / 0^0 / inf^0 indeterminate seeds). Fixes
   `Limit[Log[1 - (Log[Exp[z]/z - 1] + Log[z])/z]/z, z -> 100]` hang.
2. **WP-2 (b^x series kernel).** `so_inv` caps its iteration count
   based on input-coefficient leaf size; previously 1/(2^x - 3^x) at
   x = 0 spun in the O(N^2) simp-per-iteration loop because Mathilda
   doesn't canonicalise polynomials in symbolic `Log[2], Log[3]`. The
   cap is sized so numeric/trivial inputs retain the full order and
   heavy symbolic inputs still produce a valid leading-term Laurent.
   Also: the single-simp-per-iteration rewrite in `so_inv` cuts the
   per-step cost by roughly 3x on expression growth.
3. **WP-4 (pole sign-disagreement classifier).** New
   `LIMIT_DIR_REALS` tag (distinct from the implicit `LIMIT_DIR_TWOSIDED`
   default). At an odd-order pole with two sides disagreeing in sign,
   `Direction -> Reals` returns `Indeterminate` (the real-line answer
   matches Mathematica's behaviour for `Tan[x] at Pi/2, Reals`), while
   the default two-sided direction keeps the old `ComplexInfinity`
   fall-back for unqualified rational-function limits.
   `Direction -> Complexes` returns `ComplexInfinity` for any pole
   regardless of parity (radial interpretation).
4. **WP-1 (multivariate path-dependence).** `run_multivariate` now
   performs polar/spherical substitution at the origin (2D/3D) and
   cross-checks with direction sampling along axes and diagonals.
   Returns `Indeterminate` when sampled paths disagree and the common
   value otherwise. Joint-at-Infinity works for all-positive-orthant
   via straight-line paths. Covered by `test_wp1_multivariate`
   (`Tan[x y]/(x y)` → 1, `y/(x+y)` → Indeterminate,
   `(x^3 + y^3)/(x^2 + y^2)` → 0, `ArcTan[y^2/(x^2+x^3)]` →
   Indeterminate, `x z/(x^2+y^2+z^2)` → Indeterminate,
   `ArcTan[y/x] at {Infinity, Infinity}` → Indeterminate). The
   origin-fast-path gains an inner-divide-by-zero scanner that catches
   0/0 shapes buried inside ArcTan/Sin/etc. that Mathilda's arithmetic
   would silently fold to 0.
5. **WP-3 (Series of x^a at nonzero point).** Two-part fix. First,
   `do_series_single`'s `try_factor_power_prefactor` is now gated to
   x0 = 0 or Infinity -- at any other x0 the Power[x, alpha] is
   expanded via the (1+u)^alpha binomial kernel rather than pulled
   out as a symbolic prefactor (which would leave the other factors
   to be expanded in isolation, producing 0). Second,
   `try_apart_preprocess` refuses to Apart-preprocess expressions
   containing `Power[base, non-rational]` because Mathilda's Apart
   collapses such inputs to 0. Together these resolve
   `Limit[3 (x^a - a x + a - 1)/(x-1)^2, x -> 1]` to
   `(3/2) a (a - 1)`.
6. **WP-6 (Sinh/Cosh exponentialisation at Infinity).** A new
   pre-Series tree walk (`rewrite_hyperbolic_to_exp`) replaces
   `Sinh[z] -> (E^z - E^-z)/2`, `Cosh[z] -> (E^z + E^-z)/2`,
   `Tanh[z] -> (E^z - E^-z)/(E^z + E^-z)` whenever the limit point is
   ±Infinity. The result is `Expand[]`-ed so the term-wise Plus layer
   can fold cancelling Exp[kx] summands. Fixes
   `(1 + Sinh[x])/Exp[x] at Infinity -> 1/2` and the Cosh twin.
   Limits at finite points are unchanged (Taylor via D suffices).
7. **WP-5 (sign-aware envelope / dominant term).** `layer_plus_termwise`
   gained a `growth_exponent_upper` helper: polynomial degree upper
   bound that treats Sin/Cos/Tanh/ArcTan/ArcCot as 0, multiplies
   through Times, takes max over Plus, multiplies by exponent in
   `Power[base, nonnegative_int]`, and treats
   `Power[base, negative_int]` with growing base as 0 (reciprocal of
   something diverging is bounded by 0). When one Plus summand has
   strictly larger growth than every other and its limit is ±Infinity,
   that summand wins and we return its limit. Covers
   `x^2 + x Sin[x^2]` → Infinity, `x + Sin[x]` → Infinity (bounded
   oscillator absorbed by dominant polynomial).
8. **WP-9 (complex-direction branch cuts).** New dir tag
   `LIMIT_DIR_IMAGINARY` (for numeric-imaginary directions: `I`, `k I`,
   `Complex[0, k > 0]`). The analytic layers compute the principal-
   branch result as usual; a post-pass in `builtin_limit` then
   conjugates the imaginary part via `ReplaceAll[I -> -I]` when the
   direction was `I` (landing on the branch below the cut), and
   returns `Indeterminate` for `Direction -> Complexes` when the
   principal-branch result picked up an imaginary part (a branch-point
   value). Any signed-infinity pole under `Complexes` collapses to
   `ComplexInfinity`. Also: `Limit[{f1, f2, ...}, ...]` now threads
   over a top-level List in its first argument.
9. **WP-7 (Gruntz-lite for Log[sum]).** A partial Gruntz-style
   rewrite for the iterated-log family: `Log[dom + rest...]` at
   +Infinity with a unique dominant summand gets rewritten as
   `Log[dom] + Log[1 + rest/dom]`, then recursed. Handles
   `Log[x + Log[x]] - Log[x] -> 0` and
   `Log[x^2 + x] - 2 Log[x] -> 0` directly. **Full Hardy-field
   comparative asymptotics** (stacked `Exp[Exp[-x/(1+Exp[-x])]]` etc.,
   multi-level log-exp dominance, `Sin[x] + Log[x-a]/Log[E^x-E^a]` as
   x → a) remain future work -- a fully compliant Gruntz algorithm
   is ~600-800 LOC on top of a proper MRV (most-rapidly-varying)
   rewriter, which we have not yet built. The groundwork here
   (growth_exponent_upper, Log[sum] rewrite, hyperbolic
   exponentialisation) is the platform for that next step.

**Known limitations carried forward (still unevaluated):**
* Stacked-exponential limits like `Exp[Exp[-x/(1+Exp[-x])]] ...` at
  Infinity -- requires full Gruntz MRV.
* `Limit[Exp[x] (Exp[1/x - Exp[-x]] - Exp[1/x]), x -> Infinity]`
  (should be -1): first-order asymptotic cancellation between two
  close exponentials, out of reach without Gruntz.
* `Limit[Sin[x] + Log[x-a]/Log[E^x - E^a], x -> a]` (should be
  `1 + Sin[a]`): the Log/Log ratio needs two iterations of
  L'Hospital with analytic simplification in between; the current
  L'Hospital guardrail trips before it gets there.
* `Limit[Exp[Log[Log[x + Exp[Log[x] Log[Log[x]]]]]/
                Log[Log[Log[Exp[x] + x + Log[x]]]]], x -> Infinity]`
  (should be `E`): classic Gruntz case with three levels of
  logarithmic dominance ranking.
* Iterated limits involving `ArcTan[y/x]` where the inner limit
  variable's sign must be inferred from the outer context -- e.g.
  `Limit[ArcTan[y/x], {x -> Infinity, y -> Infinity}]`: we'd need
  `Assumptions -> y > 0` style inference.

**Known limitations** (remaining; these return the original
`Limit[...]` unevaluated):
* Multivariate limits that are truly path-dependent (e.g.
  `x y / (x^2 + y^2)` at `{0,0}`) are left unevaluated rather than
  returning `Indeterminate`; the polar-substitution heuristic from
  Layer 5a of the spec is future work.
* `(1 + Sinh[x])/Exp[x]` at `Infinity`: Series at Infinity does not
  recognise the natural `Exp[-x]` factorisation, and L'Hospital stalls
  on the Cosh/Sinh derivative cycle.
* `Csc[x]/E^x` at `Infinity` and the `(x Sin[x])/(x + Sin[x])` /
  `(x^2(1+Sin[x]^2))/(x+Sin[x])^2` shapes: these are genuinely
  `Indeterminate` (bounded but sign-oscillating numerator divided by a
  factor that does not dominate). A conservative "Indeterminate
  classifier" is future work.
* `(1 - E^(-x))^E^x` at `Infinity`: should reduce to `1/E` via
  Series-of-log but the inner limit `Log[1 - E^(-x)] * E^x` sits
  outside the current log-expansion capabilities.
* `Limit[(1/(1-x))^(-1/x^2), x -> 0]`: two-sided, the from-above and
  from-below limits are respectively `0` and `Infinity`, so the
  two-sided answer is ambiguous. Our dispatcher returns `1` via a
  fall-through path (strictly speaking both `Indeterminate` and `1`
  are defensible depending on convention; Mathematica returns
  `Indeterminate`).


---

