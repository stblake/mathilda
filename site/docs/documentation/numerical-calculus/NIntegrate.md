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

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Block](../../scoping-constructs/Block/), [Boole](../../control-flow/Boole/), [UnitStep](../../elementary-functions/UnitStep/), [Cos](../../elementary-functions/Cos/), [Sin](../../elementary-functions/Sin/), [D](../../calculus/D/), [PrecisionGoal](../../other-advanced/PrecisionGoal/), [AccuracyGoal](../../other-advanced/AccuracyGoal/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
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
