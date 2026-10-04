# NIntegrate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NIntegrate[f, {x, xmin, xmax}]`**

gives a numerical approximation to the integral of f with respect to x from xmin to xmax.

**`NIntegrate[f, {x, xmin, xmax}, {y, ymin, ymax}, ...] evaluates a multidimensional integral by adaptive cubature over a constant box, or iterated 1D quadrature when an inner bound depends on an outer variable. The variable is localised (HoldAll). xmin/xmax may be Infinity, -Infinity, or complex (a straight-line contour); extra nodes {x, x0, x1, ..., xk} give a piecewise-linear contour or mark interior singularities. Method -> Automatic chooses globally-adaptive Gauss-Kronrod for smooth finite integrands, double-exponential (tanh-sinh / sinh-sinh / exp-sinh) for endpoint singularities and infinite ranges and high precision, a Levin/zeros scheme for oscillatory integrands, an exponential endpoint map plus integration-between-the-zeros for an oscillatory endpoint singularity, and Monte-Carlo for high dimensions and region (Boole) integrands. Machine or arbitrary precision via WorkingPrecision.`**

<details>
<summary>Notes</summary>

Options: Method (Automatic | GlobalAdaptive | GaussKronrodRule | DoubleExponential | TrapezoidalRule | LevinRule | OscillatorySingularity | MonteCarlo | QuasiMonteCarlo | AdaptiveMonteCarlo | PrincipalValue), WorkingPrecision (default MachinePrecision), PrecisionGoal, AccuracyGoal, MaxRecursion, MinRecursion, MaxPoints, Exclusions.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= NIntegrate[Cos[x], {x, 0, Pi/2}]
Out[1]= 1.0

In[2]:= NIntegrate[Exp[-x^2], {x, 0, Infinity}]
Out[2]= 0.886227

In[3]:= NIntegrate[1/Sqrt[x], {x, 0, 1}]
Out[3]= 2.0

In[4]:= NIntegrate[Sin[x]/x, {x, 0, Infinity}]
Out[4]= 1.5708

In[5]:= NIntegrate[Exp[-x^2 - y^2], {x, -Infinity, Infinity}, {y, -Infinity, Infinity}]
Out[5]= 3.14159
```

### Applications (5)

```mathematica
In[6]:= NIntegrate[Sin[x], {x, 0, Pi}]
Out[6]= 2.0

In[7]:= NIntegrate[Exp[-x^2], {x, -Infinity, Infinity}]
Out[7]= 1.77245

In[8]:= NIntegrate[Exp[-x^2], {x, -Infinity, Infinity}, WorkingPrecision -> 30]
Out[8]= 1.772453850905516027298167483341

In[9]:= NIntegrate[Sin[x]/x, {x, 0, Infinity}]
Out[9]= 1.5708

In[10]:= NIntegrate[Log[x] Log[1 - x], {x, 0, 1}]
Out[10]= 0.355066
```

## Algorithm

```text
nint.c — NIntegrate[f, {x, xmin, xmax}, opts]   (see nint.h)
```

Phase 1: one-dimensional integrals over a finite real interval at machine

```text
precision, via globally-adaptive Gauss-Kronrod (gkadapt).  HoldAll: the
```

integrand and bounds are held, the bounds are evaluated to numbers, then the integration variable is Block-localised and the integrand is evaluated /

```text
numericalised at each sample point.  Subsequent phases layer endpoint
```

singularities, infinite ranges, complex contours, arbitrary precision, multidimensional iteration, oscillatory and Monte-Carlo methods, Exclusions and principal values on top of this same sampling machinery.

Memory contract: never frees `res`; returns a fresh Expr* or NULL; restores the variable binding on every return path.

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| NI 50-digit Gaussian | 6.87 s | 2.11 s | 1.1 s |
| NI 2-D ridge | 2.54 s | 43.5 s | 0.685 s |
| NI oscillatory k=40 | 0.515 s | 4.05 s | 0.002 s |
| NI oscillatory k=200 | 0.477 s | 21.5 s | 0.002 s |
| NI oscillatory k=1000 | 0.389 s | 140 s | 0.002 s |
| NI oscillatory k=1001 nonzero | 0.382 s | 1.03 s | 0.626 s |

## Implementation notes

**Algorithm.** `builtin_nintegrate` peels options into an `NiOpts`, validates the
`{x, a, b}` specs, and dispatches. A finite real interval under `Method ->
Automatic` runs a best-of cascade in `ni_core_finite` — adaptive Gauss–Kronrod,
then tanh–sinh, then Levin collocation (if an `f·{cos|sin|e^{i}}` kernel is
detected), then an exponential endpoint-singularity map, then
integration-between-the-zeros — each tried only while no earlier engine has
converged (`ni_consider` keeps the converged result, else the smallest error
estimate). The engines:

- **Gauss–Kronrod** (`gkadapt.c`): the 7-point Gauss rule embedded in the
  15-point Kronrod rule (G7–K15, QUADPACK abscissae), with subintervals in a
  binary max-heap keyed by local error so the worst panel bisects first (QAG);
  running totals feed Wynn's ε (shared `seqaccel`) for endpoint singularities
  (QAGS). Error = `|resk − resg|` refined by the QUADPACK `resasc` heuristic.
- **Double-exponential** (`denint.c`, `dequad.c`): tanh–sinh on a finite
  interval, exp–sinh on a half line (`x = a + exp((π/2)sinh t)`), sinh–sinh on
  the whole line; trapezoidal on the transformed variable, step-halved level by
  level, converging double-exponentially.
- **Oscillatory**: Levin collocation (`levincoll.c`) solves the Levin ODE
  `p' + i g' p = f` in a Chebyshev basis at Chebyshev–Gauss–Lobatto nodes via a
  LAPACK complex LU (accuracy *improves* with oscillation rate; the factored
  matrix is reused across right-hand sides); `oscint.c` locates half-periods by
  sign changes and sums/Wynn-extrapolates one lobe per panel; `oscde.c` is the
  MPFR Ooura–Mori Fourier double-exponential for an aligned `amp·{sin|cos}(ωx)`
  half-line integrand.
- **Multidimensional**: adaptive Genz–Malik degree-7 (embedded degree-5 error)
  box cubature (`cubature.c`, region max-heap), or iterated 1-D quadrature
  (dependent inner bounds via `HoldAll` + `Block` localisation), or a 2-D Levin
  reduction, or quasi-Monte-Carlo (Halton radical-inverse, `mcint.c`) chosen
  automatically for a `Boole`/region integrand or ≥ 6 specs.

Infinite and complex endpoints give exp-sinh rays and piecewise-linear contours
(`ni_run_contour`, each segment `z0 + t(z1−z0)`); `PrincipalValue` uses a
symmetric mirror sample so a simple pole cancels; `Exclusions` are solved with
`Solve` and split the interval so the tanh–sinh rule straddles each singularity.

**Data structures.** Machine path is `double _Complex` throughout, with a
Block-style snapshot/restore of the variable's OwnValues and
`arith_warnings_mute` around each sample. Arbitrary precision is **MPFR**
(`mpfr_t`, run at `target + 64` guard bits and rounded back), with complex
values carried as hand-rolled `mpfr_t` pairs (`qcx`/`ncpx`) — there is no
Arb/MPC dependency; every engine has an MPFR twin. The integrand is
**auto-compiled to bytecode** lazily on the first sample (`autocompile`), with
the interpreter as the exact fallback at any point the compiled real program
reports non-real/singular. Panels (GK) and regions (cubature) live in binary
max-heaps; Levin caches its LU factors.

**Complexity / limits.** GK is adaptive QAG/QAGS; DE converges
double-exponentially (hundreds of digits in a few thousand smooth samples);
Levin is `O(n³)` per solve (`O(n²)` prepared) and cheaper the faster the
oscillation; MC is `O(1/√N)`. Options: `WorkingPrecision`, `AccuracyGoal`
(default `MachinePrecision`), `PrecisionGoal` (default `WorkingPrecision/2`),
`MaxRecursion`, `MinRecursion`, `MaxPoints`, `Exclusions`. An explicitly
requested but unimplemented method (`ClenshawCurtisRule`, `LobattoKronrodRule`,
…) warns `NIntegrate::method` and stays unevaluated rather than silently
approximating; all diagnostics route through `mth_message` (gated on the
arithmetic-warning mute so an inner iterated integral stays quiet).

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Block](../../scoping-constructs/Block/), [Boole](../../control-flow/Boole/), [UnitStep](../../elementary-functions/UnitStep/), [Cos](../../elementary-functions/Cos/), [Sin](../../elementary-functions/Sin/), [D](../../calculus/D/), [PrecisionGoal](../../other-advanced/PrecisionGoal/), [AccuracyGoal](../../other-advanced/AccuracyGoal/)

- R. Piessens, E. de Doncker-Kapenga, C. W. Überhuber and D. K. Kahaner, *QUADPACK: A Subroutine Package for Automatic Integration* (Springer, 1983) — the QAG/QAGS adaptive Gauss–Kronrod strategy.
- A. C. Genz and A. A. Malik, *An adaptive algorithm for numerical integration over an N-dimensional rectangular region*, J. Comput. Appl. Math. **6** (1980) 295–302.
- H. Takahasi and M. Mori, *Double exponential formulas for numerical integration*, Publ. RIMS Kyoto **9** (1974) 721–741.
- T. Ooura and M. Mori, *A robust double exponential formula for Fourier-type integrals*, J. Comput. Appl. Math. **112** (1999) 229–241.
- D. Levin, *Fast integration of rapidly oscillatory functions*, J. Comput. Appl. Math. **67** (1996) 95–101.
- Source: [`src/numerical_calculus/nint.c`](https://github.com/stblake/mathilda/blob/main/src/numerical_calculus/nint.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_integrate_newton_leibniz.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_newton_leibniz.c)
- Tests: [`tests/test_nint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nint.c)

## Notes & additional examples

### Notes

`NIntegrate[f, {x, a, b}]` approximates a definite integral. The Gaussian
example reproduces `Sqrt[Pi] = 1.77245...`, computed to 30 digits with
`WorkingPrecision -> 30` via the double-exponential rule on the infinite range.
The Dirichlet integral `Sin[x]/x` over `[0, Infinity]` is the oscillatory case,
returning `Pi/2`. The final integral has the closed form `2 - Pi^2/6 =
0.355066...`. `Method -> Automatic` selects globally-adaptive Gauss-Kronrod,
double-exponential, Levin oscillatory, or Monte-Carlo schemes per region.
Endpoints may be infinite or complex (a contour).
