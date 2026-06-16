# Mathilda bugs & discrepancies

Observations surfaced while authoring documentation examples (trivial +
non-trivial) across all builtins. Recorded for follow-up, not fixed here.
Each "wrong result / broken" item below was reproduced directly against the
current `./Mathilda` build.

## Wrong results / broken functionality

- **`Coefficient[poly, x, 0]` leaks negative-power terms on a Laurent
  polynomial.** `Coefficient[Expand[(x + 1/x)^6], x, 0]` returns
  `20 + 1/x^6 + 6/x^4 + 15/x^2` instead of the constant coefficient `20`.
  The `x^0` extraction collects every non-positive-power term, not just the
  constant.
- **`EvenQ`/`OddQ` of a symbol returns `False` instead of staying
  unevaluated, corrupting sums.** `Sum[k Boole[EvenQ[k]], {k, 1, 10}]`
  returns `0` (should be `30`): the summand is evaluated once with symbolic
  `k`, `EvenQ[k] -> False`, `Boole -> 0`, so every term is zero before the
  iterator substitutes values. Root cause: `EvenQ[k]`/`OddQ[k]` should remain
  unevaluated (or the summand should not be pre-evaluated for predicates).
- **Matrix `Norm` is non-functional.** `Norm[{{1,2},{3,4}}]` and the
  `p = 1, 2, Infinity` forms all return unevaluated, although the docstring
  promises Frobenius and induced operator norms. Scalar and vector norms work.
- **`Flatten[expr, Infinity]` returns unevaluated.** The symbolic `Infinity`
  level spec is unhandled; integer levels and the no-level form work.
- **`N[Root[Function[t, poly], k]]` does not numericalize.** The documented
  primary `Function[t, p[t]]` input form stays a held `Root[...]`; only the
  `Slot` form `N[Root[#^2 - 2 &, k]]` evaluates numerically.
- **`Variables[lhs == rhs]` wraps the whole equation.**
  `Variables[2 x + 3 y == 5]` returns `{2 x + 3 y == 5}` instead of `{x, y}`
  (standard Mathematica returns the variable list). Polynomial inputs work.
- **Definite integration mis-threads the bound list.**
  `Integrate[f, {x, a, b}]` threads `{x, a, b}` as a Listable list argument
  and returns garbage (e.g. `{F[x], Integrate[f, a], Integrate[f, b]}`)
  rather than the definite integral. (Already documented in the Integrate
  overlay; only the indefinite single-variable form is supported.)
- **`FixedPoint` returns unevaluated on a non-converging (oscillating)
  iteration** instead of stopping at the default `MaxIterations`. Repro:
  `FixedPoint[# - (#^2 - 2)/(2 #) &, 1.0]` stalls oscillating between adjacent
  machine doubles and yields the unevaluated expression.

## Missing simplifications / discrepancies vs standard Mathematica

These return unevaluated or unsimplified rather than wrong — gaps, not defects.

- `Erf[x] + Erfc[x]` does not simplify to `1`.
- `Element[Sqrt[2], Algebraics]` stays unevaluated (algebraic-number
  recognition of radicals not implemented).
- `Im`/`Re` do not distribute symbolically: `Im[Sin[x + I y]]`,
  `Im[E^(I Pi/3)]`, `Im[Gamma[1 + I]]` stay unevaluated (only explicit
  `Complex`/numeric forms resolve).
- `Log[I]` and `Log[E^(I Pi/3)]` stay unevaluated (imaginary-axis / complex
  log); real-axis cases like `Log[-1] -> I Pi` work.
- Trig on degree/period arguments not reduced: `Sin[30 Degree]` is not exact
  `1/2` symbolically (only under `N`); `Tan[x + Pi]` periodicity not applied.
- `Limit[Tanh[x], x -> Infinity]` returns unevaluated (should be `1`).
- `TrigFactor` sum-to-product not implemented:
  `TrigFactor[Sin[x] + Sin[y]]` stays unevaluated.
- `Simplify[GoldenRatio^2 - GoldenRatio - 1]` does not reduce to `0`
  (GoldenRatio is an opaque protected constant; only collapses under `N`).
- `Quartics -> True` is a no-op for the general irreducible quartic
  (`Solve[x^4 + x + 1 == 0, x, Quartics -> True]` returns the same `Root[]`
  objects as the default); biquadratic/special quartics do reduce to radicals.
- `Sum[1/n^2, {n, 1, Infinity}]` and other zeta-type infinite series return
  unevaluated (no closed form); geometric and exponential infinite sums work.
- `PrimePi[Pi]` (real argument) stays unevaluated.
- `Coefficient` / `CoefficientList` do not accept a `SeriesData` directly
  (must `Normal` first).
- `ChebyshevT[n, x]` / `Cyclotomic[n, x]` are not expanded by `Expand`
  (stay in inert head form), so `CoefficientList` of them yields one element.
- Nested-radical denesting (`Sqrt[3 + 2 Sqrt[2]]`) is not automatic; needs
  `FullSimplify`.

## Unimplemented functions encountered (return unevaluated)

`Divisors`, `DivisorSigma`, `Divisible`, `CoprimeQ`, `Prime`, `Reduce`,
`Refine`, `FunctionExpand`, `ComplexExpand`, `MapIndexed`, `StringCases`,
`MatrixForm`. Note: `Map[f, Divisors[n]]` silently descends into the
unevaluated `Divisors[n]` (e.g. `Map[EulerPhi, Divisors[30]] -> Divisors[8]`),
which looks like a wrong answer rather than an inert passthrough.
