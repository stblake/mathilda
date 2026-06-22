# Numerical Calculus

Numerical (as opposed to symbolic) calculus routines: the numerical derivative
`ND`, numerical integration `NIntegrate`, numerical summation `NSum`, numerical
products `NProduct`, numerical limits `NLimit`, the numerical series expansion
`NSeries`, and the numerical residue `NResidue`. They return machine-precision
or arbitrary-precision (MPFR) numbers, and several succeed on inputs that the
symbolic engine cannot close — essential singularities, non-elementary
integrands, slowly-convergent or alternating sums, and limits with no closed
form. Each routine is built on a shared kernel library (adaptive Gauss–Kronrod
quadrature, double-exponential quadrature, Levin/oscillatory schemes, cubature,
Monte-Carlo, sequence acceleration, and contour integration) under
`src/numerical_calculus/`.

These functions complement the symbolic [`calculus`](calculus.md) routines
(`D`, `Integrate`, `Limit`, `Sum`) and [`power series`](power-series.md)
(`Series`): use the symbolic form when a closed form is wanted and the numerical
form when only a number is needed or the symbolic engine gives up.

## NResidue

Numerical residue by contour integration.  `NResidue[expr, {z, z0}]`
estimates the residue of `expr` at `z = z0` -- the coefficient of
`(z - z0)^-1` in the Laurent expansion -- by integrating around a small
circle in the complex plane.  Implemented natively in C in
`src/numerical_calculus/`: the reusable periodic-trapezoidal contour core
lives in `quadrature.{c,h}`, the builtin in `nresidue.{c,h}`.  Attribute:
`Protected`.

Unlike the symbolic `Residue` (which needs a power series at `z0`),
`NResidue` works for **essential singularities** such as `Exp[1/x]` and
`Sin[1/x]`.  It cannot tell a tiny spurious residual from a true zero --
`Chop` the result when needed -- and returns an incorrect value if the
contour encloses another singularity or crosses a branch cut.

### Method

The residue equals the Cauchy integral over the circle of radius `r`:

```
Res(f, z0) = (1/2 pi i) oint f dz = (r/N) sum_{k=0}^{N-1} f(z0 + r e^{i th_k}) e^{i th_k},  th_k = 2 pi k / N.
```

The integrand is 2*pi*-periodic and analytic in *theta*, so the periodic
trapezoidal rule converges geometrically (Trefethen & Weideman, *SIAM
Review* 2014).  The engine doubles `N` (reusing samples) until
`|S_{2N} - S_N|` meets the precision goal, applies Aitken/Shanks
extrapolation to the doubling sequence, and runs at machine precision or,
when `WorkingPrecision` requests it, in MPFR complex arithmetic.

### Forms

- `NResidue[expr, {z, z0}]` -- residue of `expr` near `z = z0`.
- `NResidue[{e1, e2, ...}, {z, z0}]` -- threads element-wise over the
  first argument (manual threading; `NResidue` is deliberately **not**
  `Listable`, so the `{z, z0}` spec is never split).

### Options

| Option | Default | Meaning |
|--------|---------|---------|
| `Radius`           | `1/100`         | Contour radius, or `Automatic` for an adaptive (Fornberg/Bornemann-style) search that favours fast-converging radii. |
| `WorkingPrecision` | `MachinePrecision` | Machine doubles, or MPFR at the requested decimal precision. |
| `PrecisionGoal`    | `WorkingPrecision - 2` | Target accuracy (decimal digits). |
| `MaxRecursion`     | `10`            | Maximum number of `N`-doublings. |
| `Method`           | `Trapezoidal`   | Only the trapezoidal rule is implemented. |

### Beyond Mathematica's NResidue

- `Radius -> Automatic` removes the manual radius tuning the reference
  implementation requires (its own `Exp[1/x]` example fails until you
  guess `Radius -> 1`).
- A relative-jump / decay-rate diagnostic flags a branch cut that
  **crosses** the contour (`NResidue::bcut`) instead of silently returning
  garbage.  A cut lying entirely *inside* the disk (e.g.
  `Sqrt[x-1] Sqrt[x+1]` on `Radius -> 2`, where the integrand is analytic
  on the circle) is undetectable on the contour and returns the same value
  as Mathematica.
- Aitken/Shanks extrapolation plus a reported error estimate.

### Diagnostics (stderr)

| Tag                | Triggered when |
|--------------------|----------------|
| `NResidue::ivar`   | Second argument is not a `{z, z0}` list, or the variable is not a symbol. |
| `NResidue::nnum`   | `z0` is not numeric, or `expr` did not evaluate to a number on the contour. |
| `NResidue::ncvi`   | Did not converge to the precision goal within `MaxRecursion` (the best estimate is still returned). |
| `NResidue::bcut`   | The integrand appears non-analytic on the contour (a branch-cut crossing); the result is unreliable. |
| `NResidue::badopt` / `NResidue::badmeth` | Invalid option value / unsupported `Method`. |

### Examples

```mathematica
In[1]:= NResidue[1/x, {x, 0}]
Out[1]= 1.

In[2]:= NResidue[Sin[1/(10 x)], {x, 0}] // Chop
Out[2]= 0.1

In[3]:= NResidue[1/(1.7 - 2.7 z + z^2), {z, 1.}] // Chop
Out[3]= -1.42857

In[4]:= NResidue[Exp[1/x], {x, 0}, Radius -> 1] // Chop
Out[4]= 1.

In[5]:= NResidue[{Exp[1/x], Sin[1/x], Cos[1/x]}, {x, 0}, Radius -> 1] // Chop
Out[5]= {1., 1., 0}

In[6]:= NResidue[1/x + 1/(x + 0.005), {x, 0}, Radius -> 0.001] // Chop
Out[6]= 1.

In[7]:= NResidue[Exp[1/x], {x, 0}, Radius -> Automatic] // Chop
Out[7]= 1.

In[8]:= 10! NResidue[Zeta[x]/x^11, {x, 0}, Radius -> 1/2, WorkingPrecision -> 30]
Out[8]= -3.6287999994567658842202915*10^6   (= Derivative[10][Zeta][0])
```


## ND

Numerical derivative.  `ND[expr, x, x0]` gives a numerical approximation to
the derivative of `expr` with respect to `x` at `x = x0`; `ND[expr, {x, n},
x0]` gives the `n`-th derivative.  Implemented natively in C in
`src/numerical_calculus/nderiv.{c,h}`.  Attribute: `Protected`.  Like
`NResidue`, `ND` cannot tell a tiny spurious residual from a true zero --
`Chop` when needed.

### Methods

**`Method -> EulerSum`** (default).  Richardson (Romberg/Neville)
extrapolation of the `n`-th **forward** finite difference taken along the
complex direction `Scale = s`:

```
D(h) = (1/(s h)^n) sum_{k=0}^n (-1)^{n-k} C(n,k) f(x0 + k s h),   h_i = 2^-i,
T(i,0) = D(h_i),   T(i,j) = T(i,j-1) + (T(i,j-1) - T(i-1,j-1)) / (2^j - 1),
result = T(Terms-1, Terms-1).
```

The forward (one-sided) stencil along `s` is what produces directional and
one-sided derivatives -- the left/right derivatives of `Abs`, and complex
directions such as `Scale -> 1 + I`.  Because a forward difference has an
error expansion in *all* powers of `h`, the Richardson denominator is
`2^j - 1` (not the `4^j - 1` of a central stencil).  EulerSum samples only
along `s`, so it works for **non-analytic** `expr` (e.g. `Re[Cos[I y]]`).  It
requires an integer order `n >= 0` and, for high-order derivatives, fights
subtractive cancellation via higher `WorkingPrecision` and more `Terms`.

**`Method -> NIntegrate`**.  Cauchy's integral formula, evaluated by reusing
`NResidue`:

```
f^(n)(x0) = n! Res_{z=x0} f(z)/(z - x0)^(n+1)
          = Gamma(n+1) * NResidue[expr/(x-x0)^(n+1), {x, x0}, Radius -> Scale].
```

`Scale` is the contour radius (default `1`).  `Gamma(n+1)` in place of `n!`
lets the order be **fractional or complex** (e.g. `ND[x, {x, -1/2}, 1]` =
`4/(3 Sqrt[Pi])`).  This method requires `expr` to be **analytic** near `x0`
and silently returns the derivative of the analytic continuation otherwise
(so `Re[Cos[I y]]` gives a wrong answer -- use EulerSum there).

### Forms and options

- `ND[expr, x, x0]` -- first derivative at `x0` (equivalent to `{x, 1}`).
- `ND[expr, {x, n}, x0]` -- `n`-th derivative.
- `ND[{e1, e2, ...}, x, x0]` -- threads element-wise over the first argument
  (manual threading; `ND` is deliberately **not** `Listable`, which would
  split the `{x, n}` spec).

| Option | Default | Meaning |
|--------|---------|---------|
| `Method` | `EulerSum` | `EulerSum` or `NIntegrate`. |
| `Scale` | `1` | EulerSum: step size / complex direction.  NIntegrate: contour radius. |
| `Terms` | `7` | EulerSum extrapolation depth. |
| `WorkingPrecision` | `MachinePrecision` | machine `double` or, for a digit count, MPFR. |
| `PrecisionGoal` | `Automatic` | target accuracy (NIntegrate). |
| `MaxRecursion` | `10` | max contour refinements (NIntegrate). |

| Message | When |
|---------|------|
| `ND::ivar` | The second argument is not `x` or `{x, n}` with `x` a symbol. |
| `ND::nnum` | `x0`/`Scale` is not numeric, or `expr` did not evaluate to a number at a sample point. |
| `ND::ord`  | `Method -> EulerSum` was given a non-integer order (use `NIntegrate`). |
| `ND::badscl` | `Scale` did not evaluate to a nonzero number. |
| `ND::badopt` / `ND::badmeth` | Invalid option value / unsupported `Method`. |

### Examples

```mathematica
In[1]:= ND[Exp[x], x, 1]
Out[1]= 2.71828

In[2]:= ND[Cos[x]^3, {x, 2}, 0]
Out[2]= -3.

In[3]:= ND[Sin[x], x, Pi I]
Out[3]= 11.592 + 1.32527*10^-10 I

In[4]:= ND[{Exp[x], Sin[x]}, x, 1]
Out[4]= {2.71828, 0.540302}

In[5]:= ND[Re[Cos[I y]], y, 1]          (* non-analytic: use EulerSum *)
Out[5]= 1.1752

In[6]:= ND[Abs[x], {x, 1}, 0, Scale -> 1 + I]
Out[6]= 0.707107 - 0.707107 I

In[7]:= ND[Sin[100 x], x, 0, Scale -> 1/100]
Out[7]= 100.

In[8]:= ND[Exp[x^2], {x, 4}, 0, Method -> NIntegrate]
Out[8]= 12. - 3.3723*10^-15 I

In[9]:= ND[x, {x, -1/2}, 1, Method -> NIntegrate]   (* = 4/(3 Sqrt[Pi]) *)
Out[9]= 0.752253

In[10]:= ND[Sin[x^2], {x, 3}, 1, Terms -> 20, WorkingPrecision -> 40]
Out[10]= -14.420070264639875819037588981065446865125
```


## NSeries

Numerical Taylor/Laurent series.  `NSeries[f, {x, x0, n}]` gives a numerical
approximation to the series expansion of `f` about `x = x0`, including the
terms `(x - x0)^-n` through `(x - x0)^n`, as a `SeriesData` object.
Implemented natively in C in `src/numerical_calculus/nseries.{c,h}`.
Attribute: `Protected`.

Unlike the symbolic `Series` (which needs a power series at `x0`), `NSeries`
needs only to **sample `f` numerically**, so it works for functions whose
coefficients have no closed form and for **Laurent expansions about essential
singularities** such as `Sin[x + 1/x]`.  It cannot tell a tiny spurious
residual from a true zero -- `Chop` the result when needed -- and returns an
incorrect value if the disk centred at `x0` contains a branch cut of `f`.
For a Laurent result, the `SeriesData` neglects higher-order poles.

### Method

`f` is sampled at `N` equispaced points on a circle of radius `r` centred at
`x0`, and a discrete Fourier transform of the samples recovers the Laurent
coefficients via Cauchy's integral formula (Lyness & Sande, 1971;
Bornemann, *FoCM* 2011):

```
z_j = x0 + r e^{2 pi i j / N},  j = 0 .. N-1,
c_k = (1/N) sum_j f(z_j) e^{-2 pi i j k / N},     (DFT of the samples)
a_e = c_{e mod N} * r^{-e}        for e = -n .. n  (coeff of (x-x0)^e).
```

The upper-half DFT bins (`k = N - m`) supply the **negative**-power
coefficients, so one transform yields both the principal part and the
analytic part; this is exact when `f` is analytic on an annulus containing
the circle.  The sample count is a power of two with an oversampling margin,
`N = 2^{ceil(log2 n) + 2}`, which pushes the leading aliased term below the
round-off floor.  A direct `O(N^2)` DFT is used (no FFT dependency): `N` is
small and each sample requires a full symbolic evaluation of `f`, which
dominates the runtime.  The same path serves machine (`double _Complex`) and
arbitrary-precision (MPFR) computations.

### Result

`SeriesData[x, x0, {a_-n, ..., a_n}, -n, n+1, 1]` -- the coefficient list runs
from exponent `-n` upward, with an `O[(x-x0)^{n+1}]` term.  Coefficients are
real (`Real`/MPFR) when their imaginary part is exactly zero, else
`Complex[re, im]`.

### Options

| Option | Default | Meaning |
|--------|---------|---------|
| `Radius`           | `1`             | Radius of the sampled circle; picks the annulus within which a Laurent series converges. |
| `WorkingPrecision` | `MachinePrecision` | Machine doubles, or MPFR at the requested decimal precision (also shrinks spurious imaginary residuals). |

### Diagnostics (stderr)

| Tag              | Triggered when |
|------------------|----------------|
| `NSeries::ivar`  | Second argument is not a `{x, x0, n}` list, the variable is not a symbol, or `n` is not a non-negative integer. |
| `NSeries::nnum`  | `x0` is not numeric, or `f` did not evaluate to a number on the contour. |
| `NSeries::badopt`| Invalid option value or unrecognised option. |

### Examples

```mathematica
In[1]:= NSeries[Exp[x], {x, 0, 5}] // Chop
Out[1]= 1. + x + 0.5 x^2 + 0.166667 x^3 + 0.0416667 x^4 + 0.00833333 x^5 + O[x]^6

In[2]:= NSeries[Exp[x], {x, I, 5}] // Chop
Out[2]= (0.540302 + 0.841471 I) + (0.540302 + 0.841471 I) (x - I) + ... + O[x - I]^6

In[3]:= NSeries[Sin[x + 1/x], {x, 0, 10}] // Chop
Out[3]= 2.49234*10^-6/x^9 - 0.000174944/x^7 + ... + 0.576725/x + 0.576725 x - ... + O[x]^11

In[4]:= NSeries[1/((1 + x) (3 + x)), {x, 0, 10}, Radius -> 5] // Chop
Out[4]= 9841./x^10 - 3280./x^9 + 1093./x^8 - 364./x^7 + 121./x^6 - 40./x^5 + 13./x^4 - 4./x^3 + 1/x^2 + O[x]^11

In[5]:= NSeries[Exp[x], {x, 0, 5}, WorkingPrecision -> 30] // Chop
Out[5]= 1. + x + 0.5 x^2 + 0.16666666666666666... x^3 + ... + O[x]^6
```

### Notes

- Unlike Mathematica's `NSeries` (whose documented count is
  `2^{ceil(log2 n) + 1}`), Mathilda oversamples by `+2` for a wider
  anti-aliasing margin.
- `Series[Sin[x + 1/x], {x, 0, 10}]` returns unevaluated (no power series at
  the essential singularity); `NSeries` recovers the Laurent expansion.


## NLimit

Numerical limit.  `NLimit[expr, z -> z0]` numerically finds the limiting value
of `expr` as `z` approaches `z0`.  Implemented natively in C in
`src/numerical_calculus/nlimit.{c,h}`.  Attribute: `Protected`.  Like the other
numerical-calculus builtins, `NLimit` cannot tell a tiny spurious residual from
a true zero -- `Chop` the result when needed.

`NLimit` constructs a geometric sequence of sample points approaching `z0` and
recovers the limit by sequence acceleration:

- **Finite `z0`:** samples `z_k = z0 - d * Scale * 2^-k`, where `d` is the
  (unit) `Direction` vector.  The points lie on the `-d` side of `z0`, so one
  moves *along* `d` to reach it.
- **Infinite `z0`** (`Infinity`, `-Infinity`, `I Infinity`,
  `DirectedInfinity[d]`, `ComplexInfinity`): samples march outward on the
  point's ray from the origin, `z_k = u * Scale * 2^k`.

### Methods

- **`EulerSum`** (default) -- Richardson / Romberg extrapolation of the sample
  sequence, using the all-powers denominator `2^j - 1` (the same convention as
  `ND`'s `EulerSum`).  Best for smooth power-series approaches; depth is set by
  `Terms`.
- **`SequenceLimit`** -- Wynn's epsilon algorithm (iterated Shanks transform).
  Exact in one step for a geometric / exponential tail; the number of
  iterations is `WynnDegree` (which needs at least `2(WynnDegree + 1)` terms).
  The estimate is read from the `ε_{2·WynnDegree}` column at the entry that best
  agrees with its neighbour (avoiding the roundoff-amplified bottom corner).

### Options

| Option | Default | Meaning |
|--------|---------|---------|
| `Method` | `EulerSum` | `EulerSum` or `SequenceLimit`. |
| `WorkingPrecision` | `MachinePrecision` | `MachinePrecision`, or digits → MPFR. |
| `Direction` | `Automatic` (≡ `-1`) | complex approach vector for finite `z0`. |
| `Scale` | `1` | initial step (finite) / distance from origin (infinite). |
| `Terms` | `7` | number of sample points / extrapolation depth. |
| `WynnDegree` | `1` | `SequenceLimit` iterations. |

`Direction -> Automatic` (`-1`) approaches a finite point from larger values;
`Direction -> 1` from smaller; complex rays such as `-Exp[225 Degree I]` select
an arbitrary direction in the complex plane (essential for path-dependent
limits and branch cuts).

### Robustness

The last two extrapolates are compared against the *sample scale* (the largest
`|S_k|`).  A divergent or non-settling sequence -- e.g. a power-law approach to
infinity like `NLimit[1/x, x -> 0]` -- yields `NLimit::noise` and the form is
returned unevaluated.  A bounded but only roughly-resolved limit is returned
rather than refused, matching Mathematica.

### Beyond / unlike Mathematica's NLimit

- `EulerSum` is implemented as Richardson/Romberg extrapolation (consistent
  with `ND`), not Euler series summation; results agree on smooth cases.
- Some sequences that make Mathematica's `EulerSum` fail (e.g.
  `NLimit[Tanh[x], x -> Infinity]` at the default `Terms`) succeed here because
  Richardson does not divide by a vanishing Euler ratio.

### Messages

| Message | When |
|---------|------|
| `NLimit::noise` | Cannot recognise a limiting value (divergent / noisy sequence). |
| `NLimit::notnum` | `expr` not numerical at a sample point, or the point/Direction/Scale is not numerical. |
| `NLimit::ndterm` | Not enough `Terms` for the chosen `Method`. |
| `NLimit::badopt` / `NLimit::badmeth` | Invalid option value / unsupported `Method`. |

### Examples

```mathematica
In[1]:= NLimit[Sin[x]/x, x -> 0]
Out[1]= 1.

In[2]:= NLimit[(1 + 1/n)^n, n -> Infinity]
Out[2]= 2.71828

In[3]:= NLimit[(1 + I/x)^x, x -> Infinity]
Out[3]= 0.540302 + 0.841471 I

In[4]:= NLimit[Tanh[Pi x]/(1 + x^2), x -> I] // Chop
Out[4]= -1.5708 I

In[5]:= NLimit[(10^x - 1)/x, x -> 0, Terms -> 10, Method -> SequenceLimit]
Out[5]= 2.30259                                   (= Log[10])

In[6]:= NLimit[z + Conjugate[z]/z, z -> 0, Direction -> -I] // Chop
Out[6]= -1.

In[7]:= NLimit[Tan[z], z -> Infinity I, Method -> SequenceLimit] // Chop
Out[7]= 1. I

In[8]:= NLimit[(2^x - 1)/x, x -> 0, WorkingPrecision -> 30, Terms -> 14]
Out[8]= 0.693147180559945309417232121458   (= Log[2])

In[9]:= NLimit[1/x, x -> 0]
        NLimit::noise: Cannot recognize a limiting value. ...
Out[9]= NLimit[1/x, x -> 0]
```


## NSum

Numerical summation.  `NSum[f, {i, imin, imax}]` gives a numerical
approximation to the sum of `f` for `i` running from `imin` to `imax` (which
may be `Infinity`).  Implemented natively in C in
`src/numerical_calculus/nsum.{c,h}`.  Attributes: `HoldAll, Protected` -- the
summand and the iterator bounds are held, and the index is `Block`-localised.
Like the other numerical-calculus builtins, `NSum` cannot tell a tiny spurious
residual from a true zero -- `Chop` when needed.

### Forms

- `NSum[f, {i, imin, imax}]` -- sum with unit step (`{i, imax}` means
  `imin = 1`).
- `NSum[f, {i, imin, imax, di}]` -- step `di`; terms are reindexed to
  `x_k = imin + k di`.
- `NSum[f, {i, ...}, {j, ...}, ...]` -- multidimensional sum; an inner bound
  may depend on an outer index (e.g. `{k, 1, n}`), handled by making the outer
  summand an inner `NSum`.

### Methods

The terms are reindexed to `k = 0, 1, 2, …`; the head terms (`NSumTerms`, default
15) are always summed explicitly, after which the tail is approximated.

- **`EulerMaclaurin`** (alias `Integrate`) -- explicit head terms plus the
  Euler–Maclaurin tail
  `(1/di) ∫_N^∞ f dx + f(N)/2 − Σ_{j≥1} B_{2j}/(2j)! · di^{2j-1} f^(2j-1)(N)`,
  with `N = imin + max(NSumTerms, settle)·di` (the head extends through any late
  peak so the tail integral starts in the monotone region).  The tail integral
  uses a self-contained **double-exponential (exp-sinh)** quadrature
  (`dequad.{c,h}`) whose tolerance and depth scale with `WorkingPrecision`.  The
  derivative corrections are **hybrid**: symbolic `D` + `BernoulliB` while the
  derivative tree stays small (robust, and the original path for simple
  summands), switching to **numerical contour (circle-DFT) derivatives** once it
  balloons — so composite summands like `Log[1 + 1/n^2]` no longer truncate
  early.  An oscillatory (sign-alternating) extension stays on the symbolic path
  (its contour DFT is ill-conditioned).  At arbitrary precision the summand is
  evaluated with guard digits so near-1 cancellation (e.g. `Log[1 + 1/x^2]`)
  does not eat into the result.  Best for monotone, slowly-converging series.
- **`AlternatingSigns`** -- the Cohen–Villegas–Zagier (2000) algorithm: a single
  pass over `n` terms with Chebyshev weights `d_n = ((3+√8)^n + (3+√8)^{-n})/2`
  delivering ≈ `2.54 n` bits.  The state of the art for alternating series.
- **`WynnEpsilon`** (alias `SequenceLimit`) -- Wynn's epsilon algorithm applied
  to the partial sums (shared with `NLimit` via `seqaccel.{c,h}`).  General
  fallback; excellent for alternating / geometric tails, weak on monotone ones.
- **`Automatic`** (default) -- probes the first terms, and (when the head is not
  already monotone) a geometric far-tail ladder that locates a late peak or
  sustained growth far beyond the head window.  Chooses `AlternatingSigns` for a
  strictly alternating decreasing summand, `EulerMaclaurin` for a monotone /
  late-settling tail, else `WynnEpsilon`.  The far-tail ladder also drives
  convergence verification, so a summand that merely peaks late (e.g.
  `1/(1 + (k-20)^2)`) is no longer mistaken for divergent.

A **large finite** sum is evaluated as the difference of two infinite tails,
`Σ_{imin}^∞ − Σ_{imax+di}^∞`, when the summand decays.

### Options

| Option | Default | Meaning |
|--------|---------|---------|
| `Method` | `Automatic` | `Automatic`, `EulerMaclaurin`, `AlternatingSigns`, `WynnEpsilon`. |
| `WorkingPrecision` | `MachinePrecision` | `MachinePrecision`, or digits → MPFR. |
| `NSumTerms` | `15` | head terms summed explicitly before extrapolation. |
| `NSumExtraTerms` | auto | length of the Wynn partial-sum sequence. |
| `WynnDegree` | `1` | `WynnEpsilon` iterations. |
| `VerifyConvergence` | `True` | ratio-test divergence check (infinite sums). |
| `AccuracyGoal` / `PrecisionGoal` | `Infinity` / `Automatic` | target tolerances. |

### Convergence

For infinite sums `VerifyConvergence -> True` runs a tail ratio test; a clearly
divergent sum (`|a_{k+1}/a_k| > 1`) yields `NSum::div` and `ComplexInfinity`.
The test is deliberately blind to ratios → 1, so (like Mathematica) it does not
detect the divergence of `Σ 1/k`.  `VerifyConvergence -> False` skips the test
and returns the formal accelerated value (e.g. `NSum[2^i, {i,0,Infinity}]`
gives the Shanks value `-1`).

### Messages

| Message | When |
|---------|------|
| `NSum::div`  | The sum does not appear to converge (returns `ComplexInfinity`). |
| `NSum::ncvg` | The extrapolation did not converge (try more `NSumExtraTerms` or higher `WorkingPrecision`). |
| `NSum::nnum` | The summand did not evaluate to a number at a term. |
| `NSum::badopt` | Invalid option value. |

### Examples

```mathematica
In[1]:= NSum[(-5)^i/i!, {i, 0, Infinity}, NSumTerms -> 25] - Exp[-5]
Out[1]= 1.4*10^-15

In[2]:= NSum[1/i^2, {i, 1, Infinity}] - Pi^2/6 // N
Out[2]= 2.2*10^-16

In[3]:= NSum[1/n^(11/10), {n, 1, Infinity}, WorkingPrecision -> 40] - Zeta[11/10]
Out[3]= -2.78*10^-27

In[4]:= NSum[(-1)^x/(1 + (x - 12)^2), {x, 0, Infinity}, Method -> "AlternatingSigns", WorkingPrecision -> 30]
Out[4]= 0.275193859413953039568971561592

In[5]:= NSum[1/2^i, {i, 0, Infinity, 2}]
Out[5]= 1.33333                                    (= 4/3)

In[6]:= NSum[Log[x]/x^(2 + 2 I), {x, 1, Infinity}]
Out[6]= -0.182175 - 0.136618 I

In[7]:= NSum[1/i^2, {i, 100, 10^6}]
Out[7]= 0.0100492

In[8]:= NSum[(-1)^n (2/n)^k/k^2, {n, 2, Infinity}, {k, 1, n}]
Out[8]= 0.770188

In[9]:= NSum[2^i, {i, 0, Infinity}]
        NSum::div: the sum does not appear to converge
Out[9]= ComplexInfinity
```

### Resolved limitations

Three `NSum` deficiencies surfaced while validating `NProduct` have been fixed
(2026-06-14):

1. **Peaked / late-settling summands** (was false divergence).  A geometric
   far-tail ladder now sees structure beyond the head probe, so a summand that
   peaks late (e.g. `NSum[1/(1+(k-20)^2), {k, 0, Infinity}]` ≈ 3.10462) is summed
   correctly and never mistaken for divergent — while genuinely divergent series
   (`2^k`) are still flagged.
2. **Accuracy at high `WorkingPrecision`** (was an EM ceiling).  Scaling the
   tail-integral tolerance to `WorkingPrecision`, computing the correction series
   from numerical contour derivatives once symbolic `D` balloons, and evaluating
   the summand with guard digits now reach the requested precision on composite
   summands: `NSum[Log[1+1/n^2], …, WP -> 35]` and
   `NSum[Log[E^(-1/(2n))(1+1/(2n))], …, WP -> 35]` give ~33–34 correct digits
   (were ~8–12).

### Remaining limitation

**Nested mixed alternating + smooth summands.**  An infinite-*outer*
multidimensional **product** whose log-summand is alternating with a smooth
(non-alternating) tail — `NProduct[1+(-1)^n (2/n)^k/k^2, {n,2,Infinity},
{k,1,n}]` ≈ 0.607 vs the true ≈ 0.564 — is not summed to high accuracy by any
single classical accelerator (Euler–Maclaurin, Wynn, or Cohen–Villegas–Zagier
each miss it).  The routing is correct (no invalid EM, no hang); the residual
~few-percent error is inherent to the extrapolation.  Infinite-outer
multidimensional **sums** (`NSum[(-1)^n (2/n)^k/k^2, {n,2,Infinity},
{k,1,Infinity}]` ≈ 1.14434, `{k,1,n}` ≈ 0.770188) and all finite ranges are
accurate.


## NProduct

Numerical multiplication.  `NProduct[f, {i, imin, imax}]` gives a numerical
approximation to the product of `f` for `i` from `imin` to `imax` (which may be
`Infinity`); `NProduct[f, {i, imin, imax, di}]` uses step `di`, and
`NProduct[f, {i,…}, {j,…}, …]` is multidimensional (an inner bound may depend on
an outer index).  Implemented in `src/numerical_calculus/nprod.{c,h}`.
Attributes: `HoldAll, Protected`.

Per Keiper (1992) the product is evaluated as `Exp[NSum[Log[f], …]]`, reusing
the full NSum engine: Euler-Maclaurin (`Method -> "EulerMaclaurin"`, default for
monotone factors), Wynn epsilon (`Method -> "WynnEpsilon"`), automatic method
selection, MPFR working precision, large-finite tail differences, and the
convergence test (factors are checked for `-> 1`; a divergent product such as
`∏(1+2^i)` returns `ComplexInfinity`).  Options: `Method`, `WorkingPrecision`,
`NProductFactors` (leading factors taken explicitly, default 15),
`NProductExtraFactors`, `WynnDegree`, `VerifyConvergence`, `AccuracyGoal`,
`PrecisionGoal`.  On the arbitrary-precision path NProduct carries guard digits
because `Exp` turns the exponent's absolute error into the product's relative
error.

Like Mathematica, NProduct can miss the divergence of slowly diverging products
(e.g. `∏(1+1/k)`, whose log-sum is the harmonic series) and may leave a harmless
`+0. I` residue on products of real negative factors.

### Examples

```mathematica
In[1]:= NProduct[1 - 1/n^2, {n, 2, Infinity}]
Out[1]= 0.5

In[2]:= NProduct[(n^2)/(n^2 - 1), {n, 2, Infinity}]
Out[2]= 2.0

In[3]:= NProduct[1 + 1/n^2, {n, 1, Infinity}]
Out[3]= 3.67608
```

`In[1]` is the telescoping product `∏(1−1/n²) = 1/2`; `In[2]` is its reciprocal,
`2`; `In[3]` is `Sinh[Pi]/Pi ≈ 3.67608`.


## NIntegrate

Numerical integration.  `NIntegrate[f, {x, xmin, xmax}]` approximates
∫ f dx; `NIntegrate[f, {x,…}, {y,…}, …]` is multidimensional — adaptive
Genz-Malik cubature over a constant rectangular box, or iterated 1D quadrature
when an inner bound depends on an outer variable (or is infinite/complex).
Implemented in
`src/numerical_calculus/nint.{c,h}` with the rule kernels `gkadapt`
(Gauss-Kronrod), `denint` (tanh-sinh / sinh-sinh), `dequad` (exp-sinh, shared
with NSum), `oscint` (oscillatory), `cubature` (adaptive Genz-Malik
multidimensional cubature), and `mcint` (Monte-Carlo).  Attributes:
`HoldAll, Protected`.  The integration variable is `Block`-localised and the
integrand evaluated/numericalised at each sample point.

### Method selection (`Method -> Automatic`)

| Region / integrand | Engine |
|--------------------|--------|
| finite real, smooth | globally-adaptive Gauss-Kronrod (G7-K15) with QAGS Wynn extrapolation |
| finite real, endpoint singularity | tanh-sinh double-exponential |
| arbitrary `WorkingPrecision` (> machine) | double-exponential at the target precision + guard bits (MPFR) |
| semi-infinite / doubly-infinite | exp-sinh / sinh-sinh |
| oscillatory (many periods / slow tail) | integrate between zeros + Wynn epsilon (finite: half-period panels) |
| doubly-infinite, non-decaying oscillation (e.g. `Exp[I x^2]`, `Cos[x^2]`) | split at 0, each half integrated between the zeros + Wynn epsilon (sinh-sinh fallback) |
| oscillatory endpoint singularity (e.g. `Cos[Log[x]/x]/x` at 0) | exponential endpoint map `x = a + (b−a)e^{−t}` onto a half line, then integrate between the (accelerating) zeros + Wynn epsilon |
| multidimensional, constant rectangular box (2 ≤ d ≤ 5) | adaptive Genz-Malik cubature (degree-7 / degree-5 error) |
| multidimensional, variable-dependent / infinite bounds | iterated 1D quadrature |
| region (`Boole`/`UnitStep`) or dimension ≥ 6 | (quasi-)Monte-Carlo |
| complex `xmin`/`xmax` or extra nodes | straight-line / piecewise-linear contour; complex ∞ gives a ray |

Named methods are accepted as strings: `"GlobalAdaptive"`, `"GaussKronrodRule"`,
`"DoubleExponential"`, `"TrapezoidalRule"`, `"RiemannRule"`, `"NewtonCotesRule"`,
`"LevinRule"`, `"OscillatorySingularity"`, `"MonteCarlo"`, `"QuasiMonteCarlo"`,
`"AdaptiveMonteCarlo"`, `"PrincipalValue"`.  `"OscillatorySingularity"`
forces the exponential endpoint-map + integrate-between-the-zeros engine on the
detected singular endpoint(s) (both are tried when neither is detected singular).
A recognised but **not-yet-implemented**
method (e.g. `"ClenshawCurtisRule"`, `"LobattoKronrodRule"`) emits
`NIntegrate::method` and returns unevaluated rather than silently approximating.

#### Levin collocation (`"LevinRule"`)

For a highly oscillatory integrand of the form `f(x)·{Cos[g]|Sin[g]|Exp[I g]}`
(amplitude `f` slowly varying, kernel oscillating rapidly), `"LevinRule"` uses
Levin's collocation method: it solves the ODE `p'(x) + i g'(x) p(x) = f(x)` in a
Chebyshev basis at Chebyshev–Gauss–Lobatto nodes — a small dense complex linear
system — so the integral becomes the boundary term `p(b)e^{i g(b)} − p(a)e^{i g(a)}`
(`Cos`/`Sin` kernels take the real/imaginary part).  Accuracy **improves** with
the oscillation rate, the opposite of ordinary quadrature: e.g.
`NIntegrate[Cos[100000 x], {x,0,1}]` is resolved exactly where Gauss–Kronrod
cannot.  The phase derivative `g'` is obtained symbolically (`D`); the order is
doubled until two estimates agree.  It runs at machine and arbitrary
`WorkingPrecision` (an in-house complex-MPFR collocation solve).  When the
integrand is not of Levin form, the phase is singular at an endpoint (e.g.
`Sin[1/x]` at 0), or the oscillation is too weak (an ill-conditioned collocation
matrix), it falls back to the ordinary oscillatory cascade.  Automatic selects
it for a detected oscillatory kernel that the smooth rules fail to resolve.

**Multivariate Levin.**  Over a constant rectangular box, an oscillatory axis is
reduced by collocation so the inner integral becomes a single linear solve
rather than a nested adaptive quadrature, and the remaining axis is integrated
by the ordinary 1-D cascade — each of its samples costing one inner solve.  This
makes a separable oscillatory product such as
`NIntegrate[Sin[1/x] Cos[1000 y], {x,0,1}, {y,0,1}, Method -> "LevinRule"]`
tractable.  When the reduction axis's phase derivative is independent of the
outer variable the collocation matrix is factored once and re-used across outer
samples.  (Two-dimensional; higher dimensions fall through to cubature /
iterated quadrature.)

#### Equally-spaced composite rules

`"RiemannRule"`, `"TrapezoidalRule"` and `"NewtonCotesRule"` are the fixed
uniform-sampling Newton–Cotes family (these sample the interval *endpoints*,
unlike Gauss–Kronrod, so a non-numeric endpoint aborts the rule).  Each refines
by panel doubling until `PrecisionGoal` / `AccuracyGoal` is met (bounded by
`MaxRecursion` / `MaxPoints`), at machine **and** arbitrary `WorkingPrecision`.
They take Method sub-options in the `Method -> {"rule", "opt" -> val, …}` form:

| Rule | Sub-option | Values (default) |
|------|-----------|------------------|
| `"RiemannRule"` | `"Type"` | `"Left"` (default), `"Right"`, `"Midpoint"` |
| `"TrapezoidalRule"` | `"RombergQuadrature"` | `True` (default; Richardson/Romberg extrapolation) or `False` (plain piecewise-linear) |
| `"NewtonCotesRule"` | `"Points"` | `2` trapezoid, `3` Simpson (default), `4` Simpson-3/8, `5` Boole |

The Riemann rectangle rules are first order (midpoint second order) and are not
extrapolated; Trapezoidal and Newton–Cotes drive a Romberg table by default.
(Wolfram's same-named methods plug these as *local* rules into the adaptive
strategy, so its crude-`PrecisionGoal` values come from non-uniform subdivision;
the values here are the well-defined uniform-composite estimates.)

### Forms and features

- `NIntegrate[f, {x, x0, x1, …, xk}]` — splits at the interior nodes (handles
  singularities there) or, with complex nodes, integrates a piecewise-linear
  contour in the complex plane.
- `NIntegrate[f, {x, a, b}, Exclusions -> {p1, …}]` or `Exclusions -> (g==0)`
  — splits the domain at the given points / equation roots (best-effort via
  `Solve` for the equation form).
- `Method -> "PrincipalValue"` with `Exclusions` — Cauchy principal value via the
  symmetric mirror rule about each simple pole.
- A `List` (vector/matrix) integrand threads element-wise.

Options: `WorkingPrecision` (default `MachinePrecision`), `PrecisionGoal`,
`AccuracyGoal`, `MaxRecursion`, `MinRecursion`, `MaxPoints`, `Method`,
`Exclusions` (`EvaluationMonitor` accepted and ignored).  Diagnostics:
`NIntegrate::ncvb` (did not converge), `NIntegrate::maxp` (Monte-Carlo did not
reach the goal), `NIntegrate::method`, `NIntegrate::badopt`.

### Limitations

Monte-Carlo and complex-contour paths are machine precision (high
`WorkingPrecision` applies to real ranges).  Finite-range MPFR bounds are taken
to machine precision (exact for rational bounds).  The oscillatory zero finder
tracks a drifting local period (the lobe step adapts to the last zero-to-zero
gap), so an *endpoint* chirp is handled via the `"OscillatorySingularity"`
transform; an *interior* chirp such as `Sin[x^2]` over a wide finite range is not
specially transformed.

### Examples

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

`In[2]` is the Gaussian integral `√π/2`; `In[3]` has an inverse-square-root
endpoint singularity; `In[4]` is the Dirichlet integral `π/2`; `In[5]` is the
two-dimensional Gaussian, `π`.

## NRoots

Numerical roots of a univariate polynomial equation.
`NRoots[lhs == rhs, var]` returns a disjunction of equations
`var==r1 || var==r2 || …` approximating **every** root of `lhs - rhs` in `var`.
A root of multiplicity `k` appears as `k` identical equations; a degree-1
polynomial yields a single bare equation (no `Or`). A numerically vacuous
equation collapses first: `NRoots[1==0, x]` → `False`, `NRoots[1==1, x]` → `True`.

Implemented in `src/numerical_roots/nroots.{c,h}` (orchestration +
CompanionMatrix), `nroots_aberth.{c,h}` (Aberth–Ehrlich + Bini initialization +
shared `ncpx` polynomial helpers), and `nroots_jt.{c,h}` (Jenkins–Traub).
Attribute: `Protected`. Real and complex coefficients are handled at machine and
arbitrary precision; all numeric work is MPFR complex (the `ncpx` toolkit). The
coefficients are extracted symbolically (`Expand` + `get_all_coeffs_expanded`),
numericalized to a complex MPFR coefficient array, a trailing `x^m` factor is
deflated to `m` exact zero roots, and the reduced polynomial is dispatched to the
selected engine. **Exact integer-coefficient polynomials are first squarefree-
decomposed** (Yun's algorithm on the integer polynomial): each squarefree factor
is solved separately — well-conditioned — and its roots emitted with the factor's
multiplicity. This keeps high powers such as `(x^2-2)^30` exact and fast, where
solving the expanded degree-60 polynomial directly would be catastrophically
ill-conditioned. Results are then noise-chopped, conjugate-symmetrized (for real
polynomials), clustered so multiple roots print identically, canonically ordered
(reals ascending, then complex by Re, then |Im|, negative Im first), and rounded
to the target precision.

### Methods (`Method -> …`)

| Method | Algorithm |
|--------|-----------|
| `Automatic` (= `"Aberth"`) | Aberth–Ehrlich simultaneous iteration; cubic convergence, all roots at once, Bini circle initialization from the Newton polygon of `(k, log\|a_k\|)`. |
| `"Aberth"` | as above. |
| `"CompanionMatrix"` | Eigenvalues of the Frobenius companion matrix. Real coefficients use the existing real MPFR QR (`eigen_all_eigenvalues_real_mpfr`) directly; complex coefficients use the `C^{n×n}→R^{2n×2n}` real embedding, with genuine roots selected by residual and multiplicities read from `p` by repeated synthetic division. |
| `"JenkinsTraub"` | The three-stage shifted-deflation algorithm (CPOLY, ACM TOMS 419), used for real and complex coefficients alike; one root per deflation, polished against the original polynomial. The real-arithmetic RPOLY (TOMS 493) is a speed/storage optimization that yields identical roots and is folded into the complex path. |

All three methods agree to tolerance on the same polynomials. Aberth is the
default: it is precision- and complex-agnostic and returns clustered multiple
roots directly.

Options: `Method`, `PrecisionGoal` (`Automatic` ⇒ machine precision; a digit
count selects arbitrary precision via MPFR), `MaxIterations` (caps the Aberth
sweep; `Automatic` ⇒ `100 + 20·degree`), `StepMonitor` (accepted for
compatibility). Diagnostics: `NRoots::neqn` (not an equation),
`NRoots::ivar` (variable not a symbol), `NRoots::npoly` (not polynomial in the
variable), `NRoots::nnum` (non-numeric coefficient), `NRoots::bdmtd` (unknown
method), `NRoots::conv` (a method did not converge).

### Limitations

The complex-coefficient `"CompanionMatrix"` path uses the real `2n` embedding,
which cannot distinguish a complex-coefficient polynomial having a conjugate pair
of roots with *unequal* multiplicities (a measure-zero case); the root values are
still correct. For such inputs use `"Aberth"` or `"JenkinsTraub"`. `StepMonitor`
is accepted but not invoked per step.

### Examples

```mathematica
In[1]:= NRoots[1 + 2 x + 3 x^2 + 4 x^3 == 0, x]
Out[1]= x == -0.60583 || x == -0.0720852 - 0.638327 I || x == -0.0720852 + 0.638327 I

In[2]:= NRoots[x^2 - 2 == 0, x]
Out[2]= x == -1.41421 || x == 1.41421

In[3]:= NRoots[x^2 + 1 == 0, x]
Out[3]= x == -I || x == I

In[4]:= NRoots[(x - 1)^3 == 0, x]
Out[4]= x == 1. || x == 1. || x == 1.

In[5]:= NRoots[x^2 - (3 + 4 I) == 0, x]
Out[5]= x == -2. - I || x == 2. + I

In[6]:= NRoots[x^2 - 2 == 0, x, PrecisionGoal -> 30]
Out[6]= x == -1.414213562373095048801688724209 || x == 1.414213562373095048801688724209
```

`In[1]` is the documentation example (one real root and a conjugate pair);
`In[4]` shows multiplicity as repeated equations; `In[5]` solves a
complex-coefficient equation (`x = ±(2 + I)`); `In[6]` returns 30-digit roots.

## NSolve

Numerical solutions of an equation or system of equations.
`NSolve[expr, vars]` returns approximate solutions as a list of
replacement-rule lists; `NSolve[expr, vars, Reals]` restricts to real
solutions (the default domain is the complexes). `vars` may be a single
symbol or a list, and `NSolve[{e1, e2, …}, vars]` is the conjunction
`e1 && e2 && …`. A working precision may be given via `WorkingPrecision` or as
a trailing positional argument (`NSolve[poly, x, Reals, 30]`); with no
variable list, `NSolve[expr]` collects the variable. Results: `{}` no
solutions, `{{x->s,…},…}` the solutions, `{{}}` the universal solution. For a
single variable, roots are repeated by multiplicity.

Implemented in `src/numerical_roots/nsolve.{c}` (dispatcher) and
`nsolve_system.{c,h}` (the polynomial-system engine). Attribute: `Protected`.
The dispatcher routes each input to the most specific method:

| Input | Method |
|-------|--------|
| Univariate polynomial | `NRoots` (multiplicity repeated; `Reals` discards complex roots). |
| Square / zero-dimensional polynomial system | **Eigenvalue / multiplication-matrix** (Möller–Stetter): a Gröbner basis over `Q` gives the quotient ring `Q[x]/I`; eigenvalues/eigenvectors of multiplication by a generic linear form (via the MPFR real-matrix backend `eigen_all_eigenvectors_real_mpfr`) yield every solution, with coordinates read off as `(M_{x_i} v)[j]/v[j]`. |
| Square polynomial system, `Method -> "Symbolic"` (and the eigenvalue fallback) | **Elimination**: a lexicographic Gröbner basis is solved triangularly by NRoots back-substitution, verifying each completed tuple. |
| Linear / underdetermined / radical / inverse-function | `Solve`, numericalized to the working precision (underdetermined systems give a parametric family; inconsistent give `{}`). |
| Univariate non-polynomial Solve cannot reduce | **FindRoot grid-seeding** (best-effort, a finite sample). |

Solutions are found at the working precision (machine by default, MPFR
otherwise); integer, real, and complex coefficients are all supported, with
solutions over the complexes. Candidate solutions are verified by residual
(controlled by `VerifySolutions`).

For the `Solve` fallback specifically, every numericalized root is re-checked
against the original equation: radical substitution (`t = x^(1/q)`,
`x = t^q`) can introduce extraneous roots that `Solve` does not always discharge
symbolically, so a candidate whose residual numericalizes clearly away from
zero is dropped, and under `Reals` a manifestly complex root is dropped too
(`NSolve[Sqrt[x] + 3 x^(1/3) == 5, x]` → `{{x -> 1.80863}}`, not the two
extraneous complex roots `Solve` reports). Solutions that do not numericalize
(symbolic, or `ConditionalExpression` families with free parameters) are kept.

Guards prevent run-away inputs: a univariate polynomial whose literal degree
exceeds **10000** (e.g. `x^1000000 - 2 x + 3`) is left unevaluated with an
`NSolve::deg` message rather than handed to the (simultaneous) `NRoots` engine
or to `Solve`; and a radical exponent whose denominator-LCM exceeds the radical
solver's cap (e.g. `x^(123451/67890)`) makes the `Solve` radical path bail
instead of building an astronomically high-degree polynomial — `NSolve` then
falls back to grid-seeding for any real roots.

Options: `MaxRoots` (cap on the count), `Method`
(`Automatic | "EndomorphismMatrix" | "Homotopy" | "Symbolic"`),
`WorkingPrecision`, `VerifySolutions`, `RandomSeeding` (seeds the generic
linear form). The four methods agree on the same solution set.

### Beyond / unlike Mathematica's NSolve

The supported envelope is polynomial and linear systems plus the
`Solve`/`FindRoot` fallbacks. Not yet handled (the input is left unevaluated
or routed to a fallback rather than guessed): real-algebraic solving of strict
inequalities `>`, `>=` (no CAD); `MaxRoots -> Infinity` symbolic root families;
quantified systems (`Exists`); bounded-region holomorphic / special-function
solving; and full transcendental completeness (grid-seeding samples a bounded
region only). `!=` is honoured only as a post-filter. `"Homotopy"` is routed to
the eigenvalue engine (no continuation tracker).

### Examples

```mathematica
In[1]:= NSolve[x^5 - 2 x + 3 == 0, x, Reals]
Out[1]= {{x -> -1.42361}}

In[2]:= NSolve[{x^2 + y^2 == 1, x^3 - y^3 == 2}, {x, y}]
Out[2]= {{x -> -1.09791 + 0.839887 I, y -> 1.09791 + 0.839887 I}, … (6 solutions)}

In[3]:= NSolve[{x^2 + y^3 == 1, 2 x + 3 y == 4}, {x, y}, Reals]
Out[3]= {{x -> 7.93641, y -> -3.95761}}

In[4]:= NSolve[x + 2 y + 3 z == 4 && 3 x + 4 y + 5 z == 6 && 6 x + 7 y + 8 z == 0, {x, y, z}]
Out[4]= {}

In[5]:= NSolve[E^x - x == 7, x, Reals]
Out[5]= {{x -> -6.99909}, {x -> 2.22154}}

In[6]:= NSolve[{x^2 + y^2 == 1, x^3 - y^3 == 2}, {x, y}, WorkingPrecision -> 25]
Out[6]= {{x -> -1.097911672722823576416400 + 0.839886921615659203622803 I, …}, …}
```

`In[2]` is the global solver returning all six complex solutions;
`In[4]` is an inconsistent linear system; `In[5]` uses FindRoot grid-seeding
for a transcendental equation Solve cannot reduce.
