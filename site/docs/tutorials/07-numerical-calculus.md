# Numerical calculus

The [calculus tutorial](06-calculus.md) worked symbolically: `Integrate`, `Sum`,
`Limit`, and `Series` returned exact, closed-form answers. But not every problem
has a closed form — many integrals are non-elementary, many sums and products
have no formula, and some series live around an essential singularity that the
symbolic engine cannot expand. For those, Mathilda offers a family of
**numerical** routines, all named with a leading `N`, that return a
machine-precision (or, on request, arbitrary-precision) number.

This tutorial tours all seven: `ND`, `NIntegrate`, `NSum`, `NProduct`, `NLimit`,
`NSeries`, and `NResidue`. (The numerical *equation* solvers `FindRoot`,
`FindMinimum`, and `FindMaximum` are covered at the end of the
[calculus tutorial](06-calculus.md#numerical-calculus).)

Every transcript below was produced by the actual Mathilda binary. The numbers
you see are real output — `NSum[1/n^2, {n, 1, Infinity}]` really does print
`1.64493`, which you can recognise as π²/6.

!!! note "Numerical answers are approximations"
    A numerical routine returns the best estimate it can compute, not an exact
    value. Several of these methods sample a function in the complex plane and so
    leave tiny residuals (around `10^-10`) in coefficients that *should* be zero.
    When that happens, wrap the result in `Chop` to clean it up — you will see
    this done below. See the
    [precision tutorial](04-machine-and-arbitrary-precision-arithmetic.md) for
    what the printed decimals mean.

## Numerical integration with NIntegrate

`NIntegrate[f, {x, xmin, xmax}]` computes a **definite** integral numerically.
Unlike symbolic `Integrate` (which in this build does indefinite integrals only),
`NIntegrate` takes explicit bounds and hands back a number:

```mathematica
In[1]:= NIntegrate[Cos[x], {x, 0, Pi/2}]
Out[1]= 1.0

In[2]:= NIntegrate[Exp[-x^2], {x, 0, Infinity}]
Out[2]= 0.886227
```

`In[1]` is the area under one quarter-wave of cosine, exactly `1`. `In[2]` is the
Gaussian integral `∫₀^∞ e^(-x²) dx = √π/2 ≈ 0.886227` — a non-elementary
integrand (it has no antiderivative in closed form) over an *infinite* range,
both of which `NIntegrate` handles automatically.

The bounds can be infinite, and the integrand may have an integrable singularity
at an endpoint. Mathilda detects each case and switches to an appropriate rule
(a double-exponential transform for endpoint singularities and infinite ranges,
an oscillatory scheme for slowly-decaying waves):

```mathematica
In[1]:= NIntegrate[1/Sqrt[x], {x, 0, 1}]
Out[1]= 2.0

In[2]:= NIntegrate[Sin[x]/x, {x, 0, Infinity}]
Out[2]= 1.5708
```

`In[1]` integrates `1/√x`, which blows up at the left endpoint `x = 0` yet has
finite area `2`. `In[2]` is the famous Dirichlet integral `∫₀^∞ sin(x)/x dx =
π/2 ≈ 1.5708`, an oscillating integrand whose tails decay only like `1/x`.

Several integrals can be nested for a multidimensional integral. Here is the
two-dimensional Gaussian over the whole plane, which equals `π`:

```mathematica
In[1]:= NIntegrate[Exp[-x^2 - y^2], {x, -Infinity, Infinity}, {y, -Infinity, Infinity}]
Out[1]= 3.14159
```

You can ask for more digits with `WorkingPrecision`, and steer the method by hand
with the `Method` option (`GaussKronrodRule`, `DoubleExponential`, `LevinRule`,
`MonteCarlo`, …) when you know the structure of your integrand. The default,
`Method -> Automatic`, picks for you.

## Numerical differentiation with ND

`ND[expr, x, x0]` estimates the derivative of `expr` with respect to `x` at the
point `x0`. It is the numerical companion to symbolic `D`, useful when `expr` is
awkward to differentiate by rule or you only need the value at one point:

```mathematica
In[1]:= ND[Sin[x], x, 1]
Out[1]= 0.540302

In[2]:= ND[Exp[x], {x, 2}, 0]
Out[2]= 1.0

In[3]:= ND[Gamma[x], x, 1]
Out[3]= -0.577216
```

`In[1]` is `cos(1) ≈ 0.540302`, the derivative of sine at `x = 1`. The
`{x, n}` form in `In[2]` takes the `n`-th derivative — here the second derivative
of `eˣ` at the origin, which is `1`. `In[3]` is `Γ′(1) = -γ`, the negative of
Euler's constant `EulerGamma ≈ 0.577216` — a derivative you would not want to
work out by hand.

By default `ND` uses Richardson extrapolation of finite differences, which works
even for non-analytic functions. Switching to `Method -> NIntegrate` uses
Cauchy's integral formula instead, which additionally allows *fractional* and
*complex* derivative orders.

## Numerical summation with NSum

`NSum[f, {i, imin, imax}]` adds up a series numerically. Crucially, `imax` may be
`Infinity` — so `NSum` evaluates **infinite** sums that symbolic `Sum` (finite
only, in this build) leaves alone:

```mathematica
In[1]:= NSum[1/n^2, {n, 1, Infinity}]
Out[1]= 1.64493

In[2]:= NSum[(-1)^(n+1)/n, {n, 1, Infinity}]
Out[2]= 0.693147
```

`In[1]` is the Basel sum `Σ 1/n² = ζ(2) = π²/6 ≈ 1.64493`. `In[2]` is the
alternating harmonic series `Σ (-1)^(n+1)/n = log 2 ≈ 0.693147`, which converges
only conditionally — `NSum` recognises the alternating pattern and accelerates it
with the Cohen–Villegas–Zagier method rather than summing terms naively.

Ask for more precision with `WorkingPrecision`:

```mathematica
In[1]:= NSum[1/n^2, {n, 1, Infinity}, WorkingPrecision -> 30]
Out[1]= 1.644934066848226436472415166649
```

That is `π²/6` to 30 digits. If you feed `NSum` a series that does not converge,
it tells you so by returning `ComplexInfinity` (you can switch this convergence
check off with `VerifyConvergence -> False`).

## Numerical products with NProduct

`NProduct[f, {i, imin, imax}]` is the multiplicative analogue of `NSum`.
Internally it computes `Exp[NSum[Log[f], …]]`, so it inherits the same
acceleration and convergence machinery:

```mathematica
In[1]:= NProduct[1 - 1/n^2, {n, 2, Infinity}]
Out[1]= 0.5
```

This is the telescoping product `∏_{n≥2} (1 - 1/n²) = 1/2`. As with `NSum`, the
upper bound may be `Infinity`, the index is local (you needn't clear `n`
afterwards), and a divergent product returns `ComplexInfinity`.

## Numerical limits with NLimit

`NLimit[expr, z -> z0]` estimates a limit by sampling `expr` along a sequence of
points marching toward `z0` and accelerating that sequence. It is the numerical
counterpart to symbolic `Limit`:

```mathematica
In[1]:= NLimit[Sin[x]/x, x -> 0]
Out[1]= 1.0

In[2]:= NLimit[(1 + 1/n)^n, n -> Infinity]
Out[2]= 2.71828
```

`In[1]` is the indeterminate `0/0` form `sin(x)/x → 1`. `In[2]` is the limit that
*defines* Euler's number, `(1 + 1/n)ⁿ → e ≈ 2.71828`; note the target `z0` can be
`Infinity` (or even a complex infinite point like `I Infinity`). The approach
direction and method (`EulerSum` Richardson extrapolation, or `SequenceLimit`
Wynn's epsilon) are both options.

## Numerical series with NSeries

`NSeries[f, {x, x0, n}]` builds a power-series approximation by sampling `f` on a
circle in the complex plane and taking a discrete Fourier transform of the
samples — Cauchy's integral formula in numerical clothing. It returns an ordinary
`SeriesData` object, just like symbolic `Series`. Because the method is numerical
it leaves small residuals in the zero coefficients, so we `Chop` the result:

```mathematica
In[1]:= Chop[NSeries[Exp[x], {x, 0, 4}], 10^-6]
Out[1]= 1.0 + 1.0 x + 0.5 x^2 + 0.166667 x^3 + 0.0416667 x^4 + O[x]^5
```

Those coefficients are the Taylor coefficients `1/k!` of `eˣ`, in decimal:
`1, 1, 0.5, 0.1667, 0.04167`. Unlike symbolic `Series`, `NSeries` can produce a
**Laurent** series — with negative powers — for a function that has a pole at the
expansion point:

```mathematica
In[1]:= Chop[NSeries[Exp[x]/x^2, {x, 0, 3}], 10^-6]
Out[1]= 1.0/x^2 + 1.0/x + 0.5 + 0.166667 x + 0.0416667 x^2 + 0.00833333 x^3 + O[x]^4
```

The `1.0/x^2 + 1.0/x` head captures the double pole at the origin. `NSeries` even
works around *essential* singularities (such as `Sin[x + 1/x]`) where the
symbolic `Series` cannot — provided the sampling circle avoids any branch cut.

## Numerical residues with NResidue

A **residue** is the coefficient of `(z - z0)^-1` in the Laurent expansion of a
function about a singularity — the quantity at the heart of the residue theorem.
`NResidue[expr, {z, z0}]` computes it by integrating `expr` around a small circle
enclosing `z0`. Because the contour integral is computed numerically it leaves a
tiny imaginary residual, so — as with `NSeries` — we `Chop` the result:

```mathematica
In[1]:= Chop[NResidue[1/z, {z, 0}]]
Out[1]= 1.0

In[2]:= Chop[NResidue[Cot[z], {z, 0}]]
Out[2]= 1.0
```

`In[1]` is the residue of the simple pole of `1/z` at the origin, `1`. `In[2]`
is the residue of `cot z = cos z / sin z` at `z = 0`, also `1`. Like `NSeries`,
`NResidue` succeeds at essential singularities the symbolic `Residue` cannot
touch — here `e^(1/z)`, whose residue is the `1/1!` coefficient of its Laurent
series. The default contour radius (`1/100`) makes `e^(1/z)` enormous, so we
widen it with `Radius -> 1` and `Chop` the tiny imaginary residual:

```mathematica
In[1]:= Chop[NResidue[Exp[1/z], {z, 0}, Radius -> 1]]
Out[1]= 1.0
```

Choosing the contour `Radius` is the main control here: it must enclose the
target singularity and no other, and stay on one side of any branch cut.

## Where to next

You have now met Mathilda's full numerical-calculus toolkit: integration
(`NIntegrate`), differentiation (`ND`), summation and products (`NSum`,
`NProduct`), limits (`NLimit`), series (`NSeries`), and residues (`NResidue`).
Together with the symbolic [calculus tutorial](06-calculus.md) — and its
`FindRoot`/`FindMinimum`/`FindMaximum` section — this is the complete calculus
story.

- For the full options of each routine, see the
  [numerical calculus](../documentation/numerical-calculus/index.md) section of
  the function documentation.
- Many of the example answers (`π²/6`, `√π/2`, `log 2`, `e`) are *special*
  values — the [special functions tutorial](08-special-functions.md) introduces
  the `Gamma`, `Zeta`, `PolyLog`, and `Erf` families those constants come from.
- To revisit the guided path, return to the [tutorials index](index.md).

As always, type `?Name` at the prompt — for example `?NIntegrate` — to read a
routine's built-in help, including its full list of options.
