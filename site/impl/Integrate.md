---
references:
  - "M. Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (Springer, 2005)."
  - "B. M. Trager, \"Algebraic Factoring and Rational Function Integration\", ACM SYMSAC 1976."
  - "K. O. Geddes, S. R. Czapor, G. Labahn, *Algorithms for Computer Algebra* (Kluwer, 1992)."
  - "D. J. Jeffrey, A. D. Rich, \"The Evaluation of Trigonometric Integrals Avoiding Spurious Discontinuities\", ACM TOMS 20(1), 1994, 124–135."
  - "G. W. Cherry, \"Integration in Finite Terms with Special Functions: The Logarithmic Integral\", SIAM J. Comput. 1986."
  - "G. W. Cherry, \"An Analysis of the Rational Exponential Integral\", SIAM J. Comput. 18(5), 1989."
  - "S. Blake, \"Parallel Integration over Simple Radical Extensions II\" (ParallelMixedTower / ParallelMixedSpecial)."
  - "S. Blake, \"A Generalisation of Goursat's Algorithm for Integration in Finite Terms\" and \"A Fourth-Root Generalisation of Goursat's Algorithm\"."
  - "K. Roach, \"Symbolic Integration\", 1992 (undefined-function integrands, §1.7)."
  - "Boulnois, arXiv:2308.09619 — differentiation under the integral sign framed as a first-order ODE in the parameter."
  - "CRC Standard Mathematical Tables, 31st ed. — reference integral tables."
source:
  - src/calculus/integrate.c
  - src/calculus/intrat.c
  - src/calculus/integrate_risch_transcendental.c
  - src/calculus/integrate_newton_leibniz.c
  - src/calculus/integrate_residue.c
  - src/calculus/integrate_ramanujan.c
  - src/calculus/integrate_piecewise.c
  - src/calculus/integrate_king.c
  - src/calculus/integrate_series.c
  - src/calculus/integrate_feynparam.c
  - src/calculus/integrate_intrep.c
  - src/calculus/integrate_diffunderint.c
  - src/calculus/integrate_line.c
  - src/internal/PolyLogLadder.m
  - src/internal/ParallelMixedSpecial.m
---
**Algorithm.** `builtin_integrate` (src/calculus/integrate.c) first threads over a `List`
integrand, then splits on the second argument. A range or contour spec goes to
`integrate_definite`; a symbol goes to the indefinite cascade. Before that cascade runs, an
applied `InterpolatingFunction` gives a new antiderivative `InterpolatingFunction`, and a
`SeriesData` is integrated term by term (`series_integrate`).

*Indefinite cascade.* Under `Method -> Automatic` each stage runs only if every earlier stage
declined (`IP_STAGE`). The order is:

1. Undefined-function integrands (`try_undefined`, Roach §1.7).
2. `BronsteinRational`: Hermite reduction, then the Lazard–Rioboo–Trager log part with
   Rioboo's `LogToReal`.
3. Radical substitutions: `LinearRadicals`, `QuadraticRadicals` (Euler), `LinearRatioRadicals`
   (Möbius), and `ChebychevAlgebraic` (binomial differentials).
4. Recognisers: `Fresnel` (Sin/Cos of a quadratic) and `GammaPower` (`x^p E^(a x^m)` with a
   symbolic `p`, giving an incomplete Gamma).
5. Linearity over a `Plus`, then linear-argument substitution.
6. `Weierstrass`: the Jeffrey–Rich continuous half-angle form, with a secular
   `K Floor[(x-b)/p]` term.
7. `DerivativeDivides`.
8. `RischTranscendental`: recursive transcendental Risch with Cherry's special-function
   extensions (`Erf`, `ExpIntegralEi`, `LogIntegral`, `PolyLog`).
9. `CRCTable`, loaded lazily.
10. `PolyLogLadder`, loaded lazily from src/internal/PolyLogLadder.m. It covers `R(x) Log[a x + b]^k`
    over partial fractions, and a polynomial times a trig rational via `TrigToExp`.
11. `LogByParts`: `Log[g]` times a trig kernel, integrated by parts once.
12. `ParallelMixedTower`: the parallel Risch–Norman method over a mixed radical tower, from
    src/internal/mixed/ParallelMixed.m. It is skipped for pseudo-elliptic `F/R^p`.
13. `GeneralisedGoursat`: pseudo-elliptic reductions, and their cube- and fourth-root versions,
    to genus-0 curves.
14. `ParallelMixedSpecial`: the only stage allowed to answer with non-elementary kernels
    (`Erf`, `Gamma[s, z]`, `SinIntegral`/`CosIntegral`, `EllipticF`/`E`/`Pi`). It runs last
    because an elementary answer is always preferred.

The heuristic stages are accepted only after a differentiate-back check.
`Integrate::nonelem` is held until the cascade ends, and is emitted only if no later stage
answers.

*Definite driver.* `integrate_definite` folds the specs from the innermost outward. A complex
or polyline spec goes to `integrate_line_from_spec`. A real interval runs the following chain;
apart from step 1, every step is limited to `Automatic` or its own pinned method:

1. Exact identities: an empty interval gives 0, and `DiracDelta` sifting gives a point value.
2. Whole-line `Log[Abs[affine]]/quadratic`.
3. `Piecewise`: split at the breakpoints of `Abs`, `Sign`, `UnitStep`, `Min`/`Max` and
   `Piecewise`. Each branch's sign is proved.
4. `Residue`: upper half-plane, unit circle and Fourier/Jordan families. A pole on the axis
   gives `Integrate::idiv`.
5. `Symmetry`: odd integrands give 0, and even ones halve the interval.
6. The Mellin `Log`-power family, then `RationalLog` on [0, 1], `SeriesExpansion` (term by
   term, resummed to Catalan / `DirichletBeta` / `Zeta`) and `FeynmanParameter` (Frullani-type
   `(Σ c_i x^e_i)/Log[x]`).
7. `Beta`: Euler Beta, Beta·₂F₁ and half-Beta.
8. `TrigPower`, log-sine and Poisson log integrals.
9. King reflection: `x g(x)`, and the general `f(x) + f(a+b-x) = const` case.
10. Raabe's integral.
11. `Mellin`: the Ramanujan Master Theorem table, with each convergence strip proved under
    `Assumptions`.
12. The principal-value `x^(c-1)/(1-x^d)` case.
13. `IntegralRepresentation`: Laplace–Bessel, Bessel-K and Airy.
    `MeijerConvolution` follows it: `C x^alpha K1(x) K2(x)` on (0, ∞) for two kernels from
    `BesselJ`/`Y`/`K`, `E^(-p x)`, `E^(-p x^2)`, `Cos`, `Sin` and `(a^2 + x^2)^-r`. Each is written
    as a Meijer G of x², and the convolution theorem (PBM 2.24.1.1) gives one `MeijerG`.
    Convergence is proved under the assumptions; an unproved condition returns a
    `ConditionalExpression` (`src/internal/MeijerConvolution.m`, v0.460).
14. `NewtonLeibniz`: an indefinite antiderivative, split at the interior poles, with the
    endpoints and poles evaluated by one-sided `Limit`. With `PrincipalValue -> True`, interior
    poles are combined symmetrically.
15. `DiffUnderInt`: differentiation under the integral sign as a first-order ODE in the
    parameter. It is skipped under `PrincipalValue`.

A result that still contains the integration variable is discarded. With an explicit
`Assumptions`, the folded value is `Refine`d, and the refinement is rejected if it makes a
real value complex.

**Data structures.** Everything is an `Expr*` tree. Stages talk to each other through
context-qualified heads (`call_stage` evaluating ``Integrate`Stage[f, x]``): a stage that returns its
own head unevaluated, or `$Failed`, has declined. `IntegrateMethod` holds the parsed
`Method`: `parse_method_option` handles the indefinite form, including
`{"DerivativeDivides", "Substitution" -> u}`, and `definite_parse_method` handles the
definite form. There are two kinds of per-evaluation state, both keyed on
`eval_toplevel_id()` so they reset with each user command:
- a 32-slot fail memo (`intg_fail_*`), so a fixed-point re-entry does not search the cascade
  a second time;
- one absolute `TimeConstraint` deadline, checked between stages at every depth.

Inexact integrands are rationalised by `common_rationalize_input`, integrated exactly, and
numericalised back at the input's minimum precision. Large `.m` stages (the CRC table,
PolyLogLadder, ParallelMixed*) load on first use.

**Complexity / limits.**
- The rational stage is complete: Hermite reduction plus LRT.
- Transcendental coverage is a decision procedure for the towers that
  `RischTranscendental` models, and a heuristic beyond them.
- The automatic indefinite cascade is bounded by `TimeConstraint` (default 3 s).
  `DerivativeDivides` is measured as nearly all of the cost on integrands that fail. An
  explicit `Method` is unbounded unless a `TimeConstraint` is given.
- Definite stages are narrow, shape-gated recognisers. They either hold by construction (the
  residue theorem, Beta and Mellin strips) or are checked with `NIntegrate` (King,
  series, Feynman). Anything undecidable is left unevaluated.
- `Integrate` is symbolic only. It has no packed, NDArray or `Compile[]` surface; the
  numerical counterpart is `NIntegrate`.
- Multivariate integration is supported only as iterated definite integrals.
