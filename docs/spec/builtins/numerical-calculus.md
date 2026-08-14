# Numerical Calculus

Numerical (as opposed to symbolic) calculus routines: the numerical derivative
`ND`, numerical integration `NIntegrate`, numerical summation `NSum`, numerical
products `NProduct`, numerical limits `NLimit`, the numerical series expansion
`NSeries`, the numerical residue `NResidue`, and the numerical differential-
equation solver `NDSolve`. They return machine-precision
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

### Accuracy control (shared)

Most of these routines share a common accuracy contract set by two options.
`AccuracyGoal` (default `MachinePrecision`) is the absolute-error target in
decimal digits and `PrecisionGoal` (default `Automatic`) the relative-error
target; a result is accepted once the estimated error meets the combined
tolerance `10^-AccuracyGoal + |x|·10^-PrecisionGoal`. `MachinePrecision` and
`Automatic` both track `WorkingPrecision` — about two digits below it, so a high
`WorkingPrecision` reaches full precision rather than a fixed 16 digits — and
`Infinity` disables the corresponding term. Refinement grows adaptively up to
each routine's resource cap; if the goal still cannot be met a `Head::accgl`
warning is written to stderr and the best approximation is returned (never
`$Failed` or an unevaluated form). `NDSolve` uses the same two options but
resolves `Automatic` to `WorkingPrecision/2` (see its section).

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
| `AccuracyGoal`     | `MachinePrecision` | Absolute-error target (digits); combines with PrecisionGoal as 10^-a + \|x\|10^-p.  Automatic/MachinePrecision track WorkingPrecision; Infinity disables. |
| `PrecisionGoal`    | `Automatic`     | Relative-error target (digits) in the same combined tolerance. |
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
| `AccuracyGoal` | `MachinePrecision` | Absolute-error target (digits); combines with PrecisionGoal as 10^-a + \|x\|10^-p.  Automatic/MachinePrecision track WorkingPrecision; Infinity disables. |
| `PrecisionGoal` | `Automatic` | Relative-error target (digits) in the same combined tolerance. |
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
| `AccuracyGoal`     | `MachinePrecision` | Absolute-error target (digits); combines with PrecisionGoal as 10^-a + \|x\|10^-p.  Automatic/MachinePrecision track WorkingPrecision; Infinity disables. |
| `PrecisionGoal`    | `Automatic`     | Relative-error target (digits) in the same combined tolerance. |

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

- **`Automatic`** (default) -- runs Richardson/Romberg, Wynn's epsilon
  (at every admissible degree) and Levin's u-transform, and returns the estimate
  whose internal convergence residual is smallest.  Richardson's fixed `2^j - 1`
  denominators can only annihilate an integer-power (analytic) error tail; a
  geometric or fractional-power tail -- e.g. the `sqrt(step)` imaginary part of
  `2 ArcTan[Sqrt[(1+x)/(1-x)]]` as `x -> 1` from larger values -- defeats it but
  is captured by Wynn's epsilon.  Selecting by best self-consistency keeps
  Richardson's accuracy on smooth limits while gaining Wynn's on branch-point /
  algebraic approaches, so this case now returns a real `Pi` (imaginary residual
  `~3e-14` at the default `Terms -> 13`) instead of the spurious `Pi + 0.08 I`
  that plain Richardson produced.  Levin's u-transform is
  admitted only when the sample increments are contracting, so a divergent
  sequence cannot let it collapse to a spurious value.
- **`EulerSum`** -- Richardson / Romberg extrapolation of the sample sequence,
  using the all-powers denominator `2^j - 1` (the same convention as `ND`'s
  `EulerSum`).  Best for smooth power-series approaches; depth is set by `Terms`.
- **`SequenceLimit`** -- Wynn's epsilon algorithm (iterated Shanks transform).
  Exact in one step for a geometric / exponential tail; the number of
  iterations is `WynnDegree` (which needs at least `2(WynnDegree + 1)` terms).
  The estimate is read from the `ε_{2·WynnDegree}` column at the entry that best
  agrees with its neighbour (avoiding the roundoff-amplified bottom corner).
- **`"Levin"`** -- Levin's nonlinear transformation of the sample sequence,
  driven by remainder estimates `ω_i` built from the sample increments
  `a_i = S_i − S_{i-1}`.  `"Levin"`/`"LevinU"` use `ω_i = (β+i) a_i` (the
  u-transform, `β = 1`), `"LevinT"` uses `ω_i = a_i`, `"LevinV"` uses
  `ω_i = a_i a_{i+1}/(a_i − a_{i+1})`.  Strong on logarithmically /
  algebraically convergent approaches.

Every setting is also callable as a head of its own — `NLimit`EulerSum[f,
z -> z0]` is exactly `NLimit[f, z -> z0, Method -> "EulerSum"]` — so a method
can be named without threading an option through:

| Head | Method setting |
|------|----------------|
| `NLimit`Automatic` | `Automatic` |
| `NLimit`EulerSum` | `EulerSum` |
| `NLimit`SequenceLimit` | `SequenceLimit` |
| `NLimit`Levin` | `"Levin"` (≡ `"LevinU"`) |
| `NLimit`LevinU` / `NLimit`LevinT` / `NLimit`LevinV` | the u / t / v remainder estimate |

Each head has its own `Information` string.  The two positional arguments and
every other option (`Direction`, `Scale`, `Terms`, `WynnDegree`,
`WorkingPrecision`) are forwarded untouched; a `Method` option is dropped,
since the head already names the method.  Attributes are
`{Protected, ReadProtected}`.

### Options

| Option | Default | Meaning |
|--------|---------|---------|
| `Method` | `Automatic` | `Automatic` (best of Richardson / Wynn / Levin-u), `EulerSum`, `SequenceLimit`, or `"Levin"` / `"LevinU"` / `"LevinT"` / `"LevinV"`. |
| `WorkingPrecision` | `MachinePrecision` | `MachinePrecision`, or digits → MPFR. |
| `Direction` | `Automatic` (≡ `-1`) | complex approach vector for finite `z0`. |
| `Scale` | `1` | initial step (finite) / distance from origin (infinite). |
| `Terms` | `13` | *starting* number of sample points / extrapolation depth, grown adaptively to meet `AccuracyGoal`.  This, not `WorkingPrecision`, sets the accuracy on a branch-point / fractional-power approach: 13 samples resolve such a tail to ~12 machine digits, where the historical `7` reached only ~3. |
| `WynnDegree` | `1` | `SequenceLimit` iterations. |
| `AccuracyGoal` | `MachinePrecision` | Absolute-error target (digits); combines with PrecisionGoal as 10^-a + \|x\|10^-p.  Automatic/MachinePrecision track WorkingPrecision; Infinity disables. |
| `PrecisionGoal` | `Automatic` | Relative-error target (digits) in the same combined tolerance. |

`Direction -> Automatic` (`-1`) approaches a finite point from larger values;
`Direction -> 1` from smaller; complex rays such as `-Exp[225 Degree I]` select
an arbitrary direction in the complex plane (essential for path-dependent
limits and branch cuts).

### Robustness

Two independent gates, each returning the form unevaluated when it fires.

**Oscillatory divergence (`NLimit::osc`).**  A sequence acceleration returns a
number for *any* input, including a sequence with no limit at all, so before
trusting one `NLimit` checks that the oscillation envelope decays -- the
property that distinguishes `Sin[x]/x` (limit `0`) from `x Sin[x]` (no limit).
The check runs in two stages:

1. *Screen.*  Count the direction reversals of the sample increments
   `S_k − S_{k-1}`.  Monotone and smoothly-converging samples score zero and
   skip stage 2, so an ordinary limit costs exactly what it did before.
2. *Envelope.*  Over the default `Terms -> 13` sampling window a decaying and a
   non-decaying envelope are still indistinguishable, so the diagnosis uses its own
   much wider ladder: 20 octaves sampled twice, at `Scale 2^k` and
   `Scale φ 2^k` with `φ = (Sqrt[5] - 1)/2`.  The φ offset both doubles the
   resolution and breaks the power-of-two aliasing that would otherwise hide a
   function like `Sin[Pi x]`.  Writing `env` for the ratio of `max |f|` over the
   half nearest the limit point to `max |f|` over the half furthest from it, a
   verdict of `env > 0.6` refuses the limit.

The verdict is therefore a property of the expression, not of `Terms`, `Scale`
or `Method`.  It catches growing envelopes (`x Sin[x]`, `Sqrt[x] Sin[x]`,
`Log[x] Sin[x]`), bounded ones (`Sin[x]`, `Cos[1/x]` as `x -> 0`,
`Sin[x] Sin[x^2]`), and oscillations too slow for the sampling window to see at
all (`Sin[Log[x]]`).  Oscillations whose envelope *does* decay keep their limit:
`Sin[x]/x`, `Sin[x]/Sqrt[x]`, `Sin[x]/x^2`, `Exp[-x] Sin[x]`.

Two known blind spots.  An envelope decaying slower than about `x^(-1/10)`
across the diagnostic window cannot be told from a non-decaying one, and is
refused -- the safer of the two errors.  And a function sampled exactly on its
own zeros (`Sin[Pi x]`, whose samples at `x = 2^k` all vanish) presents no
reversal for the screen to see.

**Noise (`NLimit::noise`).**  The last two extrapolates are compared against the
*sample scale* (the largest `|S_k|`).  A divergent or non-settling sequence --
e.g. a power-law approach to infinity like `NLimit[1/x, x -> 0]` -- is returned
unevaluated.  A bounded but only roughly-resolved limit is returned rather than
refused, matching Mathematica.

### Beyond / unlike Mathematica's NLimit

- `EulerSum` is implemented as Richardson/Romberg extrapolation (consistent
  with `ND`), not Euler series summation; results agree on smooth cases.
- Some sequences that make Mathematica's `EulerSum` fail (e.g.
  `NLimit[Tanh[x], x -> Infinity]` at the default `Terms`) succeed here because
  Richardson does not divide by a vanishing Euler ratio.

### Messages

| Message | When |
|---------|------|
| `NLimit::osc` | The sampled values oscillate with a non-decaying envelope, so no limit exists. |
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

In[10]:= NLimit[Sin[x]/x, x -> 0, Method -> "Levin"]
Out[10]= 1.                                (Levin's u-transform; "LevinT"/"LevinV" select the t/v variants)

In[11]:= NLimit[x Sin[x], x -> Infinity]
         NLimit::osc: The sampled values oscillate with a non-decaying envelope; ...
Out[11]= NLimit[x Sin[x], x -> Infinity]

In[12]:= NLimit[Sin[x]/x, x -> Infinity]  (* same oscillation, decaying envelope *)
Out[12]= 0.0142131
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
- **`"Levin"`** (`"LevinU"` / `"LevinT"` / `"LevinV"`) -- Levin's nonlinear
  transformation of the partial sums (shared kernel in `seqaccel.{c,h}`), with
  remainder estimates `ω_k` built from the terms `a_k`: u uses `(β+k) a_k`
  (`β = 1`), t uses `a_k`, v uses `a_k a_{k+1}/(a_k − a_{k+1})`.  Reaches full
  `WorkingPrecision` on smoothly-convergent series (e.g. `Σ 1/n²`).
- **`Automatic`** (default) -- probes the first terms, and (when the head is not
  already monotone) a geometric far-tail ladder that locates a late peak or
  sustained growth far beyond the head window.  Chooses `AlternatingSigns` for a
  strictly alternating decreasing summand, `EulerMaclaurin` for a monotone /
  late-settling tail, else `WynnEpsilon` -- with Levin's u-transform as a purely
  additive last resort when Wynn does not converge (an existing Wynn result is
  never replaced).  The far-tail ladder also drives convergence verification, so
  a summand that merely peaks late (e.g. `1/(1 + (k-20)^2)`) is no longer
  mistaken for divergent.

A **large finite** sum is evaluated as the difference of two infinite tails,
`Σ_{imin}^∞ − Σ_{imax+di}^∞`, when the summand decays.

### Options

| Option | Default | Meaning |
|--------|---------|---------|
| `Method` | `Automatic` | `Automatic`, `EulerMaclaurin`, `AlternatingSigns`, `WynnEpsilon`, `"Levin"` / `"LevinU"` / `"LevinT"` / `"LevinV"`. |
| `WorkingPrecision` | `MachinePrecision` | `MachinePrecision`, or digits → MPFR. |
| `NSumTerms` | `15` | head terms summed explicitly before extrapolation. |
| `NSumExtraTerms` | auto | length of the Wynn partial-sum sequence. |
| `WynnDegree` | `1` | `WynnEpsilon` iterations. |
| `VerifyConvergence` | `True` | ratio-test divergence check (infinite sums). |
| `AccuracyGoal` / `PrecisionGoal` | `MachinePrecision` / `Automatic` | absolute / relative error targets in the combined tolerance (see intro). |

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

In[10]:= NSum[1/n^2, {n, 1, Infinity}, Method -> "Levin", WorkingPrecision -> 30]
Out[10]= 1.644934066848226436472415166646   (= Pi^2/6; Levin reaches full precision on smooth series)
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
with NSum), `oscint` (oscillatory: integrate-between-the-zeros + Wynn epsilon,
machine and MPFR), `oscde` (Ooura–Mori double-exponential rule for semi-infinite
Fourier integrals), `cubature` (adaptive Genz-Malik multidimensional cubature),
and `mcint` (Monte-Carlo).  Attributes:
`HoldAll, Protected`.  The integration variable is `Block`-localised and the
integrand evaluated/numericalised at each sample point.

### Method selection (`Method -> Automatic`)

| Region / integrand | Engine |
|--------------------|--------|
| finite real, smooth | globally-adaptive Gauss-Kronrod (G7-K15) with QAGS Wynn extrapolation |
| finite real, endpoint singularity | tanh-sinh double-exponential |
| arbitrary `WorkingPrecision` (> machine) | double-exponential at the target precision + guard bits (MPFR) |
| semi-infinite / doubly-infinite | exp-sinh / sinh-sinh |
| semi-infinite Fourier `∫_0^∞ amp·{Sin\|Cos}(ω x) dx` | Ooura–Mori double-exponential rule (nodes land on the oscillation's zeros) — the fast, high-precision path (hundreds of digits in a few thousand samples) |
| oscillatory (many periods / slow tail), general | integrate between zeros + Wynn epsilon (finite: half-period panels; MPFR-accelerated at high `WorkingPrecision`) |
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

#### Semi-infinite Fourier integrals (Ooura–Mori DE)

A semi-infinite oscillatory integral `∫_0^∞ amp(x)·{Sin[ω x]|Cos[ω x]} dx` with a
slowly-decaying amplitude (e.g. `x Sin[x]/(x²+4)`, `Sin[x]/x`) defeats exp-sinh —
its tail never decays monotonically — and integrate-between-the-zeros needs one
quadrature panel *per half-period*, which is ruinous at high `WorkingPrecision`
(each panel demands full precision, so hundreds of digits cost hundreds of
thousands of samples).  Automatic instead recognises the aligned Fourier form and
applies the **Ooura–Mori double-exponential rule** (`oscde`): the substitution
`x = (M/ω) φ(t)`, `φ(t) = t/(1 − exp(−2t − α(1−e^{−t}) − β(e^t−1)))`, `M = π/h`,
sends the abscissae `ω x(kh) → kπ` onto the oscillation's zeros double
exponentially, so the offset trapezoidal rule reaches hundreds of digits in a few
thousand samples (reference: Ooura & Mori, *J. Comput. Appl. Math.* 112, 1999).
The abscissa scale `M = π/h` is formed in MPFR — a double-rounded `M` would floor
`sin(x_k)` at ~1e-16 in the tail and cap accuracy at machine precision.  The rule
fires only for a zero-offset linear phase from the lower limit 0 (sine nodes on
`kπ/ω`, cosine on `(k+½)π/ω`); a misaligned phase, a complex-exponential kernel,
or a non-zero lower limit falls back to the integrate-between-the-zeros engine.

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

Options: `Method`, `AccuracyGoal` (default `MachinePrecision`; the absolute
residual a returned root must meet, combined with `PrecisionGoal` as
10^-a + |x|10^-p — a root exceeding the goal triggers an `NRoots::accgl`
warning), `PrecisionGoal` (`Automatic` ⇒ machine precision; a digit
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
`WorkingPrecision`, `AccuracyGoal` (default `MachinePrecision`) and
`PrecisionGoal` (default `Automatic`) — the absolute / relative residual a
solution must meet in the combined tolerance — `VerifySolutions`,
`RandomSeeding` (seeds the generic linear form). The four methods agree on the
same solution set.

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

## NDSolve

`NDSolve[eqns, u, {x, xmin, xmax}]` finds a numerical solution to the ordinary
differential equations `eqns` for the function `u` with independent variable `x`
on `xmin ≤ x ≤ xmax`, returning a list of rules
`{{u -> InterpolatingFunction[...]}}`. The `InterpolatingFunction` can be applied
(`u[x0]`) and differentiated (`u'[x0]`).

- `NDSolve[eqns, {u1, u2, ...}, {x, xmin, xmax}]` — a coupled system.
- `NDSolve[eqns, u[x], {x, xmin, xmax}]` — gives `u[x] -> InterpolatingFunction[...][x]`.
- `NDSolve[eqns, u, {t, tmin, tmax}, {x, xmin, xmax}]` — a partial differential
  equation over a rectangular region (method of lines); see below.
- `NDSolve[eqns, u, {t, tmin, tmax}, {x, xmin, xmax}, {y, ymin, ymax}]` — a PDE
  in two spatial dimensions, giving a 3-D `InterpolatingFunction` `u[t, x, y]`.

Equations are stated with derivatives (`u'[x]`, `u''[x]`, i.e. `Derivative`, not
`Dt`). Higher-order equations are automatically reduced to a first-order system
(state `[u, u', …, u^(n-1)]`); the ODEs must be solvable (linearly) for the
highest derivative. Initial conditions `u[x0] == c`, `u'[x0] == c`, … supply the
starting state; there must be enough to determine the solution.

### Architecture

Two-way modularity under `src/numerical_calculus/ndsolve*`:

1. **Problem class** — a user problem is compiled into a first-order reduced
   system `dY/dt = f(t, Y)` (`NdProblem`). ODE IVPs and 1-D evolution PDEs (by
   the method of lines, `ndsolve_mol.c`) are supported; multi-D PDEs, BVP and DAE
   are the deferred extension seam.
2. **Method** — each integrator is a self-contained module exposing an
   `NdStepper` single-step vtable; a shared adaptive driver owns the time loop,
   error-controlled step-size selection (WRMS norm, Hairer starting step),
   rejection, monitors, and dense output. Every method is reachable both as
   `Method -> "Name"` and as the standalone builtin `NDSolve`Name[...]`.

The symbolic right-hand side is evaluated at numeric points by Block-localizing
the reduced-state symbols and the independent variable (the NIntegrate binding
pattern); implicit steppers reuse a symbolic/finite-difference Jacobian and a
dense Gaussian-elimination solve (the FindRoot pattern). Output is built as
cubic-Hermite `{{t_i}, y_i, y'_i}` triples fed to `Interpolation`.

### Methods

| Method | Order | Kind |
|---|---|---|
| `"ExplicitEuler"` | 1 | explicit, fixed |
| `"ExplicitMidpoint"` | 2 | explicit, fixed |
| `"RK4"` | 4 | classical Runge–Kutta, fixed |
| `"ExplicitRungeKutta"` / `"DOPRI5"` | 5(4) | **adaptive embedded (Automatic default)** |
| `"BackwardEuler"` | 1 | implicit (Newton), A-stable |
| `"ImplicitTrapezoid"` | 2 | implicit (Newton), A-stable |
| `"BDF"` | 1–5 | **adaptive, variable-order** backward-differentiation multistep (stiff) |
| `"Adams"` | 2 | **adaptive** Adams predictor–corrector (PECE) multistep |

The fixed and implicit one-step methods are driven with step-doubling error
control so accuracy goals are honoured. The two multistep methods (**BDF**,
**Adams**) are **adaptive variable-step**: each step forms an explicit predictor
of the same order as the corrector, and the predictor–corrector difference is a
free Milne local-error estimate that drives step accept/reject (WRMS norm ≤ 1),
just like the one-step adaptive pair. **BDF is additionally variable-order
(1–5)** — the stiff workhorse (the LSODE/CVODE family). Its coefficients are the
exact derivatives of the Lagrange interpolant on the actual nonuniform mesh
(so any step ratio and order is exact); the order ramps up while the higher
differences keep decaying and is pulled back where they do not (advection /
near-imaginary spectra, where high-order BDF is unstable); and between changes
the step and order are held constant for a run of steps so the mesh stays uniform
(BDF is zero-stable at every order ≤ 6 on a uniform mesh). It also recovers from a
diverging Newton iteration by halving the step and dropping order toward the
L-stable backward Euler — the key to solving stiff problems with **incompatible
initial/boundary corners** (e.g. the discontinuous-corner heat equation), which
fixed-step BDF could not. High order means tight tolerances are reached in far
fewer steps (a decaying test at `PrecisionGoal`/`AccuracyGoal` 10 uses ~300 steps
where an order-2 method needs thousands). The default `Automatic` method is the
adaptive Dormand–Prince 5(4) pair — the same algorithm as MATLAB's **ode45** —
and it uses the FSAL (first-same-as-last) property so the last stage of an
accepted step becomes both the node slope and the next step's first stage: **~6
function evaluations per accepted step**, matching ode45's cost. (ode23 is the
Bogacki–Shampine 3(2) pair; comparable low-order adaptive behaviour is obtained
here at cruder accuracy goals.)

### Options

`Method`, `WorkingPrecision` (machine default; higher selects the MPFR
integrator), `AccuracyGoal`/`PrecisionGoal` (Automatic = WorkingPrecision/2),
`MaxSteps` (default 10000), `MaxStepSize`, `MaxStepFraction` (1/10),
`StartingStepSize`, `InterpolationOrder` (Automatic = 3, cubic Hermite),
`StepMonitor`, `EvaluationMonitor`.

`Options[NDSolve]` returns the full default list; the `Automatic` sentinels
above are what the solver resolves internally when an option is omitted.

### Arbitrary precision

`WorkingPrecision -> p` (with `p > $MachinePrecision`) runs a dedicated MPFR
integrator whose state, independent variable, and step size are all carried at a
guard-padded precision, so non-autonomous right-hand sides are evaluated at full
precision. Explicit methods use the adaptive DOPRI5 pair (or fixed RK4);
**stiff (implicit/multistep) methods run the MPFR variable-order BDF** — the
state, node times, coefficients, Newton residual and linear solve are all MPFR
(the Jacobian stays double, since its accuracy only affects Newton's convergence
rate, not the root), so a stiff problem — where the explicit integrator would
need an impractically tiny step — is solved to full working precision (e.g. a
`1000`-stiff forced oscillator at `WorkingPrecision -> 30` reaches ~1e-19). Note
that accuracy at the achieved node values is bounded by `PrecisionGoal` (default
`p/2`), and interior interpolated values are bounded by the cubic-Hermite
interpolation order — query at an MPFR abscissa (`u[N[t, p]]`) to read a
high-precision node. High precision needs a larger `MaxSteps`.

### Partial differential equations (method of lines)

`NDSolve[eqns, u, {t, tmin, tmax}, {x, xmin, xmax}]` solves a PDE over a
rectangular region by the **method of lines**: the spatial operator is
discretized on a uniform grid (Fornberg finite-difference weights, `ndsolve_stencil.c`),
turning the PDE into the large first-order ODE system the time integrator above
already solves. The result is a **2-D `InterpolatingFunction`** over `(t, x)`,
applied as `u[t, x]`. Stencils are central where they fit and one-sided of the
same order near a boundary; `DifferenceOrder` (default 4) selects the accuracy —
error is `O(h^DifferenceOrder)`. Any spatial derivative order is handled.

The front-end lives in `ndsolve_mol.c` (registered as the controller
`` NDSolve`MethodOfLines ``, reachable via `Method -> "MethodOfLines"`). It
solves the PDE symbolically for its highest temporal derivative, then per grid
node substitutes the finite-difference stencils and the node coordinate into that
expression — so variable-coefficient and **nonlinear** PDEs in one spatial
dimension work through the same Block-localized sampler, with no extra machinery.
Because every accepted time node carries the whole spatial vector, the solution
is a complete tensor grid handed to the multidimensional `Interpolation` builtin.

**Compiled operator (efficiency).** When the discretized system is linear —
`dU/dt = A·U + s(t)` with a constant matrix — the front-end (`ndsolve_operator.c`)
compiles the numeric banded `A` and forcing `s(t)`, so the RHS is a
matrix–vector product (no per-call symbolic evaluation) and the Jacobian is
exactly `A` (free). These kernels use **BLAS/LAPACK** (Accelerate on macOS,
system BLAS/LAPACK elsewhere; scalar fallback without `USE_LAPACK`): the matvec
is `cblas_dgbmv`/`cblas_dgemv`, and the implicit iteration matrix `I − hθA` is
**LU-factored once per Newton solve** (constant Jacobian) and back-substituted
per iteration — a pivoted **banded** factor (`dgbtrf`/`dgbtrs`, `O(d·bw)`) for
the narrow FD structure, a pivoted dense factor for wide bands. Nonlinear
(symbolic-Jacobian) stiff problems detect the band and use the same banded
factor per iteration. Measured ~10× on stiff coupled grids over the previous
per-iteration dense refactorization; `Compiled -> False` forces the symbolic
sampler (identical results).

**Nonlinear RHS compiler.** When the operator does not apply (a nonlinear or
non-affine RHS), each reduced component is compiled once into a numeric
stack-machine program (`ndsolve_compile.c`) so the stepper evaluates the RHS as
bytecode over the state vector — no symbol binding, no `Expr` copy, no
`numericalize` per call. The Jacobian is then formed by Curtis–Powell–Reid
**colored finite differences** over the bytecode (`O(bandwidth)` evaluations for
a banded discretization). Any construct the VM does not support makes the
compile bail and the symbolic sampler takes over (an `EvaluationMonitor` also
forces the symbolic path). Measured ~8× on an explicit shallow-water dam break
and >100× on a stiff nonlinear-diffusion grid (the colored-FD Jacobian replaces
the per-entry symbolic Jacobian). Nonlinear ODEs share this path.

**Stiffness auto-selection:** parabolic problems (a diffusion term with
first-order time evolution) default to `"BDF"` when no time-integration method
is given.

Currently supported: one spatial dimension; one dependent function; temporal
order 1 (parabolic — heat, advection–diffusion, reaction–diffusion, Burgers) and
2 (hyperbolic — wave); arbitrary spatial derivative order; boundary conditions —
**Dirichlet** `u[t,x0]==g`, **Neumann** `Derivative[0,1][u][t,x0]==g`, **Robin**
`a u[t,x0] + b Derivative[0,1][u][t,x0]==g`, and **Periodic** `u[t,xmin]==u[t,xmax]`
(all constant or time-dependent); arbitrary-order (Fornberg) stencils via
`DifferenceOrder`; machine precision. Neumann/Robin edges are eliminated with a
one-sided finite-difference of the same order; periodic domains use cyclic
stencils. The grid and stencil are set with
`Method -> {"MethodOfLines", "SpatialDiscretization" -> {"TensorProductGrid", "MinPoints" -> n, "DifferenceOrder" -> q}}`.
Diffusion-dominated (stiff) problems should use `Method -> "BDF"` (auto-selected
for parabolic problems). **Two spatial dimensions** are also supported on a
rectangle — `NDSolve[eqns, u, {t,..}, {x,..}, {y,..}]` — with general linear
boundary conditions (**Dirichlet** `u[t,x0,y]==g`, **Neumann**
`Derivative[0,1,0][u][t,x0,y]==g`, **Robin** `a u + b u_x + r == 0`, and their
`y`-edge analogues) on each of the four edges and unmixed spatial derivatives
(2-D heat, wave, advection–diffusion); the result is a 3-D
`InterpolatingFunction` `u[t,x,y]`. As in 1-D, each Neumann/Robin boundary node
is eliminated into the interior stencils via a one-sided first-derivative
formula; corners resolve through the transverse edge.

**Coupled systems of 1-D PDEs.** `NDSolve[eqns, {u1, u2, ...}, {t,..}, {x,..}]`
solves a coupled system in one spatial dimension: each dependent function is
discretized on its own block of interior unknowns within one global first-order
ODE vector, and each function's per-node right-hand side substitutes the
stencils/node-values of *every* function, so coupling terms — `(h u)_x`,
`u u_x`, `g h_x`, cross-diffusion/reaction — resolve. The result is one
`InterpolatingFunction` per function: `{{u1 -> if1, u2 -> if2, ...}}`. Each
function may be temporal order 1 or 2 (mixed within a system), and each evolution
equation must be solvable for exactly one function's highest temporal derivative
(no coupled mass matrix). Per-function Dirichlet/Neumann/Robin or periodic BCs
plus full initial data are required. Parabolic systems auto-select `BDF`;
first-order hyperbolic systems (e.g. **shallow water** `{h_t+(h u)_x==0,
u_t+u u_x+g h_x==0}`) keep the explicit default. A system missing a condition, a
coupled mass matrix, a complex-valued system, or a 2-D system stays unevaluated
with a diagnostic. Navier–Stokes (incompressible) is **not** supported — it needs
a pressure/divergence-free constraint (a DAE) the method-of-lines driver lacks.

**Upwind schemes for advective terms** (opt-in; centered Fornberg is the
default). Two first-order schemes stabilize hyperbolic/advection-dominated
problems that centered high-order stencils would ring on:
- **Donor-cell upwind** (scalar) — `"DifferenceOrder" -> 1` or `"Upwind" -> True`
  in the `SpatialDiscretization` list. The first-derivative stencil is biased by
  the sign of the local advection speed (`wind = -∂G/∂u_x`), evaluated per node
  so the upwind direction follows a state-dependent (nonlinear) wind at runtime.
  Monotone and non-oscillatory; first-order accurate and convergent.
- **Lax–Friedrichs viscosity** (systems and scalar) — `"LaxFriedrichs" -> True`
  (or `"DifferenceOrder" -> 1` on a system, where donor-cell's per-characteristic
  wind is ill-defined). Adds grid-scaled artificial dissipation to purely
  hyperbolic equations; robust for nonlinear systems such as a dam break.

**Complex-valued PDEs (Schrödinger).** When the solved RHS carries the imaginary
unit (e.g. `I D[u[t,x],t] == -D[u[t,x],{x,2}]`), the front-end **realifies** the
system — each complex unknown is split into interleaved (Re, Im) real unknowns
via `ComplexExpand` — and returns the solution as
`u -> Function[{t,x}, ifRe[t,x] + I ifIm[t,x]]` (two real `InterpolatingFunction`s),
so `u[t,x]` evaluates to a complex number. The realified system is real and
linear, so the compiled operator still applies. Verified against the exact
semi-discrete Schrödinger eigenmode and by norm conservation (`Σ|ψ|²` drifts
~1e-10 under a potential well).

**Arbitrary precision.** `WorkingPrecision -> p` (p > machine) runs the MPFR
integrator on the discretized system, giving an MPFR-valued 2-D
`InterpolatingFunction` (1-D PDEs). Explicit methods handle non-stiff PDEs (wave,
advection); **stiff diffusion uses the MPFR variable-order BDF** (`Method -> "BDF"`,
the auto-selection for parabolic problems), so stiff PDEs now run at arbitrary
precision too. Achievable precision at interior query points is bounded by the
cubic-Hermite interpolation (query at MPFR node abscissae for full precision).

Adaptive-implicit stepping (variable-step variable-order BDF with Newton-failure
recovery) handles the incompatible IC/BC corners that previously diverged.
Remaining (deferred): periodic BCs and mixed spatial derivatives in 2-D (a
periodic coupling `u[t,xmin,y]==u[t,xmax,y]` is detected and reported rather than
mis-solved).

### Beyond / unlike Mathematica's NDSolve

Supported: ODE initial-value problems — scalar, systems, and higher-order,
including **complex-valued** ODEs (via realification: each complex unknown is
split into interleaved Re/Im real unknowns and the solution recombines to
`u -> Function[{x}, ifRe[x] + I ifIm[x]]`, integrated at machine precision) and
initial conditions posed at an **interior** point `x0 != xmin` (the driver
integrates both forward to `xmax` and backward to `xmin`) — plus 1-D and 2-D
evolution PDEs by the method of lines (above), at machine (PDEs) or arbitrary
(ODEs) precision. Not yet handled (deferred, with the `NdProblem`/`NdStepper`
seams in place): DAEs, boundary-value problems, event location
(`"EventLocator"`), and the controller methods (`"Projection"`, `"Splitting"`,
`"Composition"`, `"Extrapolation"`, symplectic integrators).
`"StiffnessSwitching"` currently maps to `"BDF"`. BDF is adaptive variable-step
variable-order (1–5); Adams is adaptive variable-step but capped at order 2. When
`MaxSteps` is exhausted the partial `InterpolatingFunction` is returned with an
`NDSolve::mxst` warning; querying it beyond the reached domain additionally emits
`InterpolatingFunction::dmval` (extrapolation), as in Mathematica.

### Examples

```mathematica
In[1]:= sol = NDSolve[{y'[x] == -y[x], y[0] == 1}, y, {x, 0, 5}];
        y[1] /. sol
Out[1]= {0.367879}                          (* = E^-1 *)

In[2]:= NDSolve[{y''[x] + y[x] == 0, y[0] == 1, y'[0] == 0}, y, {x, 0, 6}];
        y[3.0] /. %                          (* = Cos[3] *)
Out[2]= {-0.989992}

In[3]:= NDSolve[{x'[t] == y[t], y'[t] == -x[t], x[0] == 1, y[0] == 0},
                {x, y}, {t, 0, 6}]           (* circle: {Cos t, -Sin t} *)

In[4]:= NDSolve[{y'[x] == -1000 (y[x] - Cos[x]) - Sin[x], y[0] == 1},
                y, {x, 0, 3}, Method -> "BackwardEuler"]     (* stiff *)

In[5]:= NDSolve[{y'[x] == y[x], y[0] == 1}, y, {x, 0, 1},
                WorkingPrecision -> 30, PrecisionGoal -> 22, MaxSteps -> 200000]

In[6]:= (* heat equation u_t = u_xx, Dirichlet, method of lines *)
        sol = NDSolve[{D[u[t, x], t] == D[u[t, x], {x, 2}], u[0, x] == Sin[Pi x],
                       u[t, 0] == 0, u[t, 1] == 0}, u, {t, 0, 0.05}, {x, 0, 1},
                      Method -> "BDF"];
        u[0.05, 0.5] /. sol                 (* ~ E^(-Pi^2 0.05) Sin[Pi/2] *)
Out[6]= {0.612973}

In[7]:= (* wave equation u_tt = u_xx (default adaptive DOPRI5) *)
        NDSolve[{D[u[t, x], {t, 2}] == D[u[t, x], {x, 2}], u[0, x] == Sin[Pi x],
                 Derivative[1, 0][u][0, x] == 0, u[t, 0] == 0, u[t, 1] == 0},
                u, {t, 0, 0.5}, {x, 0, 1}]
```

## FindMinimum / FindMaximum

Iterative local optimisation.  Implemented natively in C in
`src/numerical_calculus/findmin.c`.  Both have `HoldAll, Protected` attributes and use
`Block`-style local binding of the search variables.  `FindMaximum[f, ...]`
is a thin wrapper around `FindMinimum[-f, ...]` that negates the
objective value before returning.

### Forms

- `FindMinimum[f, {x, x0}]` -- 1D from a single start (default Brent).
- `FindMinimum[f, {x, x0, x1}]` -- bracket Brent on `[x0, x1]`.
- `FindMinimum[f, {x, xstart, xmin, xmax}]` -- bracket Brent on `[xmin, xmax]`.
- `FindMinimum[f, {{x, x0}, {y, y0}, ...}]` -- n-D from a user start.
- `FindMinimum[f, {x, y, ...}]` -- n-D auto-start at 1 for each variable
  (matches Mathematica; avoids the common saddle-at-origin trap for
  oscillatory objectives like `Sin[x] Sin[2 y]` whose gradient vanishes
  at the origin).
- `FindMinimum[{f, cons}, vars]` -- constrained minimisation.
- `FindMaximum[...]` -- same forms; maximises `f` (equivalent to negating
  the objective and the f-value of the result).

### Output

`{f_min, {var1 -> v1, ..., varN -> vN}}` -- a 2-element list whose first
element is the function value and whose second is the rule list for the
optimising variable assignments.

### Method dispatch

| Spec                  | Default method |
|-----------------------|----------------|
| n = 1                 | Brent          |
| n >= 2                | QuasiNewton (BFGS) |
| `{x, x0, x1}` (1D)    | Brent (bracket) |

Methods overridable via `Method -> "Brent" | "Newton" | "QuasiNewton"
| "ConjugateGradient"`.  Brent is golden-section search with parabolic
interpolation (derivative-free), QuasiNewton is BFGS with Armijo
backtracking line search, ConjugateGradient is Polak-Ribière+ with
restart, Newton uses the symbolic Hessian (via repeated `D[]`) with
modified-Cholesky safeguarding.  Gradients default to a symbolic
gradient (`D[f, x_i]` per variable) with a central-difference fallback;
override via `Gradient -> {dfdx1, dfdx2, ...}`.

### Constraints

Inside the `{f, cons}` form, `cons` is a boolean tree of comparisons:

- `Less / LessEqual / Greater / GreaterEqual` between a bare iteration
  variable and a numeric constant become **box constraints** (enforced
  by per-step projection).
- Other inequalities (`g(x) <= 0`) and equalities (`h(x) == 0`) feed a
  **quadratic-penalty** wrapper around the inner solver.  The outer μ
  schedule starts at 1 and multiplies by 10 each round until feasible
  (max 9 rounds, μ up to 10^8).  The inner BFGS/CG/Newton iterations
  drive the *augmented* objective `f + μ·Σ_k max(0, g_k)^2 + μ·Σ_j h_j^2`
  using a matching *augmented* gradient `∇f + 2μ·Σ_k (active) g_k ∇g_k +
  2μ·Σ_j h_j ∇h_j` — the gradient of each constraint expression is
  computed symbolically at setup and falls back to central differences
  per-constraint when symbolic differentiation fails.
- `Or[...]`, `Element[...]`, `x ∈ Integers` and the rest of the
  Mathematica constraint surface are not yet implemented -- they emit
  `FindMinimum::nimpl`.

The penalty wrapper converges from feasible *or* infeasible starting
points on smooth nonlinear constraints (linear/quadratic inequalities,
linear/quadratic equalities, intersections thereof).  At very high μ
the inner solver may exit early on line-search exhaustion; the outer
loop's feasibility check is authoritative, and only emits a diagnostic
(`FindMinimum::infeas`) when the final iterate genuinely fails the
constraint tolerance (1e-12).

### Options

| Option              | Default        | Effect |
|---------------------|----------------|--------|
| `Method`            | `Automatic`    | `"Brent"`, `"QuasiNewton"`, `"ConjugateGradient"`, `"Newton"`, or `Automatic`. |
| `WorkingPrecision`  | `MachinePrecision` | `MachinePrecision`, or a digit count (>= ~16 routes through MPFR).  Lifts the 1D `Brent` and n-D `QuasiNewton` iterations into MPFR at the requested precision so the result `{f_min, {x -> ...}}` carries MPFR leaves with that many digits.  Explicit `Method -> "Newton"` or `"ConjugateGradient"` at MPFR currently falls back to `QuasiNewton` with a `FindMinimum::nimpl` diagnostic; general (non-box) constraints at MPFR fall back to machine precision similarly. |
| `MaxIterations`     | `500`          | Iteration limit on the inner loop. |
| `AccuracyGoal`      | `Automatic`    | Digit count `n` ⇒ stop when `\|grad\| < 10^{-n}`. `Infinity` disables. `Automatic` resolves to `WorkingPrecision/2`. |
| `PrecisionGoal`     | `Automatic`    | Digit count `n` ⇒ stop when `\|step\| < \|x\| * 10^{-n}`. |
| `Gradient`          | `Automatic`    | Explicit `{dfdx1, ..., dfdxN}` overrides the symbolic gradient. |
| `StepMonitor`       | `None`         | A held expression evaluated after each step. |
| `EvaluationMonitor` | `None`         | A held expression evaluated each time `f` (or any partial) is evaluated. |

### Diagnostics (stderr)

| Tag                  | Triggered when |
|----------------------|----------------|
| `FindMinimum::argt`    | Wrong arg count. |
| `FindMinimum::ivar`    | Malformed variable spec. |
| `FindMinimum::vecvar`  | Vector-valued variable (deferred). |
| `FindMinimum::badmeth` | Unknown `Method` value. |
| `FindMinimum::badopt`  | Unknown option name or invalid value. |
| `FindMinimum::nimpl`   | Method, constraint shape, or domain restriction not yet supported. |
| `FindMinimum::nlnum`   | `f`, gradient, or constraint did not evaluate to a number. |
| `FindMinimum::cvmit`   | Inner-loop `MaxIterations` exhausted. |
| `FindMinimum::lstol`   | Line search could not find a sufficient decrease. |
| `FindMinimum::dsing`   | Hessian non-positive-definite (Newton). |
| `FindMinimum::infeas`  | Penalty outer loop could not satisfy all constraints. |

### Examples

```mathematica
In[1]:= FindMinimum[(x - 3)^2 + 1, {x, 0}]
Out[1]= {1.0, {x -> 3.0}}

In[2]:= FindMinimum[x Cos[x], {x, 2}]
Out[2]= {-3.28837, {x -> 3.42562}}

In[3]:= FindMinimum[x Cos[x], {x, 7, 1, 15}]
Out[3]= {-9.47729, {x -> 9.52933}}

In[4]:= FindMinimum[Sin[x] Sin[2 y], {{x, 2}, {y, 2}}]
Out[4]= {-1.0, {x -> 1.5708, y -> 2.35619}}

In[5]:= FindMinimum[{x Cos[x], 1 <= x && x <= 15}, {x, 7}]
Out[5]= {-9.47729, {x -> 9.52933}}

(* Chained `lo <= x <= hi` parses as an Inequality[...] node; FindMinimum
   walks its (value, op, value) triples and registers each as a box bound. *)
In[5b]:= FindMaximum[{x Cos[x], 1 <= x <= 15}, {x, 7}]
Out[5b]= {6.36096, {x -> 6.4373}}

In[6]:= FindMinimum[(1-x)^2 + 100 (y-x^2)^2, {{x, 0}, {y, 0}}]
Out[6]= {0.0, {x -> 1.0, y -> 1.0}}

In[7]:= FindMaximum[Cos[x], {x, 0}]
Out[7]= {1.0, {x -> 0.0}}

In[8]:= FindMinimum[(x - 3)^2, {x, 0}, Method -> "ConjugateGradient"]
Out[8]= {0.0, {x -> 3.0}}

(* Arbitrary precision via WorkingPrecision: the 1D Brent and n-D BFGS
   iterations both run in MPFR at the requested precision and the
   returned `{f_min, {x -> ...}}` carries MPFR leaves with that many
   digits. *)
In[9]:= FindMinimum[(x - Pi)^2, {x, 0}, WorkingPrecision -> 50]
Out[9]= {0.0, {x -> 3.1415926535897932384626433832795028841971693993751}}

In[10]:= FindMinimum[x Cos[x], {x, 2}, WorkingPrecision -> 80]
Out[10]= {-3.2883713955908964865125964571235283975158511553846230554230811211040811736596049,
         {x -> 3.42561845948172814647771386218545617769640923939753965919739613085112431446169}}
```

## NMinimize / NMaximize

Numerical **global** optimisation.  Implemented natively in C in
`src/numerical_calculus/findmin.c`, layered on the `FindMinimum` machinery (it reuses the same
`Block`-style variable binding, `{f, cons}` constraint parsing, penalty/BFGS
local solvers, MPFR path, and result construction).  Both are `Protected` but
**not** `HoldAll` (matching Mathematica's `Attributes[NMinimize] == {Protected}`);
their variables must be unbound symbols, and an assigned optimization variable
makes the call return unevaluated with `NMinimize::ivar`.  `NMaximize[f, ...]` is a thin wrapper around
`NMinimize[-f, ...]` that negates the objective value before returning.

Where `FindMinimum` descends from a single start point to the nearest local
minimum, `NMinimize` runs a stochastic global search over a bounded region and
then polishes the best point with the exact local solver.

### Forms

- `NMinimize[f, x]` -- global minimum with respect to one variable.
- `NMinimize[f, {x, y, ...}]` -- global minimum over several variables.
- `NMinimize[{f, cons}, vars]` -- constrained global minimum.  `cons` may be
  a single constraint, an `And[...]` of constraints, or — as Mathematica
  allows — additional list elements `{f, c1, c2, ...}` that are implicitly
  `And`-ed.
- `NMaximize[...]` -- same forms; maximises `f`.

Variables may be given bare (`x`), as a starting interval (`{x, xmin, xmax}`,
used to seed the search region), with an integer domain
(`Element[x, Integers]`, in the variable list or the constraints), or as
**indexed variables** `x[i]`.  Because both arguments are held (`HoldAll`), a
generator that produces the variable list — `Table[x[i], {i, 1, n}]`,
`Array[x, n]` — is evaluated once (with the variable head localized) so it
expands to the concrete list `{x[1], ..., x[n]}`; a held `Table[...]`
constraint list is expanded the same way and treated as an implicit `And`.
Indexed variables are internally rewritten to fresh scalar symbols for the
search and mapped back in the result, so `Rule`s come back keyed by the
original `x[i]`.  This makes the standard `n`-dimensional benchmarks (e.g. the
Rosenbrock valley `Sum[100 (x[i+1]-x[i]^2)^2 + (1-x[i])^2, {i, 1, n-1}]` over
`Table[x[i], {i, 1, n}]`) express directly.

The problem and the variable spec may also be supplied through a bound symbol —
`prob = {f, cons}; vars = Table[x[i], {i, 1, n}]; NMinimize[prob, vars]`. Since
both arguments are held, such a symbol is resolved (with the variable heads
localized) so its `{f, cons}` list or variable list is exposed and then handled
exactly as if written inline; a genuinely unbound symbol stays a single
optimization variable (`NMinimize[f, x]`).

### Output

`{f_min, {var1 -> v1, ...}}`.  Integer-domain variables come back as exact
integers.  If the feasible set is empty, the result is
`{Infinity, {var1 -> Indeterminate, ...}}`.

### Constraints and domains

Constraints are the same relational/boolean forms `FindMinimum` accepts:
`==`, `<`, `<=`, `>`, `>=`, chained inequalities, and their `And`
combinations — plus disjunctive `Or` combinations, which `FindMinimum` does not
accept (see below).  Bare-variable inequalities against a constant become **box
constraints** (used both to bound the search region and to project during the
local polish); other inequalities and equalities are handled with **Deb's
feasibility rules** during the global search (a feasible point always beats an
infeasible one; among feasible points the smaller objective wins; among
infeasible points the smaller total violation wins) and with the
quadratic-penalty local solver during the polish.  Integer variables
(`Element[x, Integers]`) are searched on the integer lattice and refined by
integer coordinate descent.  A domain declaration may name **several variables
at once**, written either as `Element[x | y, Integers]` (Alternatives) or
`Element[{x, y}, Integers]` (List); each named variable gets the domain, exactly
as Mathematica does.  The declaration is accepted in the variable list or the
constraints, and every member must be one of the optimization variables (a
declaration on any other symbol is left in place and rejected as unenforceable).

**Disjunctive (`Or`) constraints.**  A constraint `c1 || c2 || ...`, at the top
level or nested inside an `And`, is feasible when **at least one** branch holds.
Each branch may itself be a single comparison, a chained inequality, or an `And`
of those.  A disjunction is scored by its **minimum-branch penalty** — the
smallest of the branches' squared-violation penalties, which is `0` exactly when
some branch is satisfied — so Deb's feasibility rules select it during the global
search with no extra weighting.  The derivative-free global engines
(DifferentialEvolution / SimulatedAnnealing / NelderMead / RandomSearch) consume
the non-smooth `min` directly; the smooth local polish folds in each
disjunction's currently-active (minimum-penalty) branch so it refines *within*
the feasible region the point already occupies, and the post-polish feasibility
gate keeps the reported point feasible.  For example,
`NMinimize[{(x^2+y-11)^2+(x+y^2-7)^2, (x-3)^2+(y-2)^2<=0.1 || (x+2.8)^2+(y+3.1)^2<=0.1}, {x, y}]`
returns the Himmelblau minimum `0` at `(3, 2)`, the branch that contains it, and
`NMinimize[{x^2, x<=-2 || x>=2}, x]` returns `4` at `x = ±2` — the branch
boundary, not the infeasible `x = 0` in the gap between branches.  (`FindMinimum`,
whose gradient penalty method requires a smooth penalty, still rejects `Or` with
`FindMinimum::nimpl`.)

For a **mixed-integer** problem whose feasible region lies outside the default
`±10` search span, the polish adds a continuous-relaxation recovery step: it
solves the continuous relaxation (the integer variables relaxed to reals, every
coordinate free of the sampling region) to locate the basin, rounds the integer
coordinates, then refines the continuous coordinates with the integers pinned.
The relaxed point is adopted only when it is a Deb-improvement, so it can never
worsen the region-confined result.  This lets e.g.
`NMinimize[{(x-15)^2 + (y-3)^2, Element[y, Integers]}, {x, y}]` recover the
optimum at `(15, 3)` even though `x = 15` is far outside the sampling box.

**Adaptive region expansion.**  When the default `±10` sampling region contains
no feasible point at all — common when the feasible set's location is implied by
nonlinear constraints rather than stated as variable bounds — the fully-unbounded
coordinates are grown by successive powers of ten (up to `±10^5`) and the search
is retried, so the feasible region is found instead of reporting
`{Infinity, ...}`.  Only fully-unbounded coordinates grow; a coordinate carrying
a box bound or a starting-interval hint keeps its resolved region.  The first
attempt is the base region with the base seed, so a problem already feasible
there is solved identically to before; expansion triggers only to rescue
infeasibility and stops at the smallest region that yields feasibility (so it
does not drift into far basins).  Thus
`NMinimize[{(x-50)^2 + (y-40)^2, x + y >= 80}, {x, y}]`, whose feasible set lies
wholly outside `±10`, returns the optimum `(50, 40)`, while a genuinely
infeasible problem (e.g. `x^2 + 1 <= 0`) still returns `{Infinity, ...}`.

A problem that is genuinely **unbounded below** returns whatever feasible point
the (bounded) search settles on; supplying explicit variable bounds (as the
standard engineering MINLP benchmarks do) both bounds the objective and pins the
intended optimum.  For instance the pressure-vessel MINLP written without its
non-negativity bounds is unbounded below (its true global minimum is a large-
negative non-physical point), so it returns a feasible point rather than the
textbook physical optimum until the standard bounds
`1 <= x1, x2 <= 99 && 10 <= x3, x4 <= 240` are added.

Not yet supported (emit `NMinimize::nimpl` and abstain / fall back): vector and
matrix variables (`Vectors[n, dom]`, `Matrices`), geometric-region domains,
`VectorGreaterEqual`, `Or[...]` (disjunctive) constraints, domains other than
`Integers`/`Reals`, and general (non-box) constraints at
`WorkingPrecision > MachinePrecision`.

### Methods

| `Method ->` | Engine |
|-------------|--------|
| `Automatic` | `"DifferentialEvolution"` with a dimension-scaled budget (see below) |
| `"DifferentialEvolution"` | DE/rand/1/bin with Deb feasibility selection (default) |
| `"NelderMead"` | downhill-simplex; each restart's vertex polished into its basin minimum, deepest kept |
| `"RandomSearch"` | multiple random starts, each refined by the local solver, best local minimum kept |
| `"SimulatedAnnealing"` | Metropolis search with geometric cooling |

`"DifferentialEvolution"` keeps every trial point inside the box by *bounce-back*
reinitialisation — a mutant that overshoots a bound is redrawn at random between
its base vector and the violated bound — rather than clamping to the bound.
Clamping would pile members onto the wall, collapse that coordinate's mutation
differentials to zero, and strand the search there whenever a box-constrained
optimum lies in the interior (e.g. the Schwefel function, optimum at
`x_i = 420.97` inside `[-500, 500]`).

On a continuous problem the local polish is a **multi-start** over the final
population, not a single refinement of the global best: the best `Min[2·d, 50]`
*distinct* members (deduplicated by basin) are each polished into their basin
minimum and the deepest is kept. Polishing only the single best strands DE in
whichever basin that one member occupied, so a larger population or a shorter run
— both of which leave the population less converged — could report a *worse*
optimum; ranking basins by their minima instead makes the result improve, not
degrade, with `"SearchPoints"`. On Griewank-10 over `[-600, 600]` the plain
`Method -> "DifferentialEvolution"` reaches the global `0` (Mathematica reports
`~0.175`). `"PostProcess" -> False` disables it and returns the raw global best;
mixed-integer problems keep the single driver polish so their integer-descent cost
is unchanged.

All four multi-start engines share this principle: **polish each independent
start into its basin minimum and keep the deepest, ranking by local minima rather
than by the raw search point that found the basin.** `"RandomSearch"` always did
(one local solve per start). `"SimulatedAnnealing"` polishes each of its
`Min[2·d, 50]` chains, `"NelderMead"` each of its `Min[2·d, 20]` restarts, and
`"DifferentialEvolution"` the distinct members of its final population. The one
engine this cannot rescue is `"RandomSearch"` itself: with no global move, pure
multi-start local search only finds basins a random start lands near, so on a box
far wider than the optimum's basin (Griewank over `[-600, 600]`, where the central
bowl is `~10⁻¹¹` of the 10-D volume) no attainable start count reaches it —
`"DifferentialEvolution"` / `"SimulatedAnnealing"` are the engines for that shape.

**Automatic vs. explicit `"DifferentialEvolution"` — the default budget.**
`Method -> Automatic` (i.e. no `Method`) runs DE with a *dimension-scaled* budget:
population `Clip[10·d, {15, 200}]` and `150·d` generations, rather than the flat
`Clip[10·d, {15, 40}]` / `100` used when `"DifferentialEvolution"` is named
explicitly. Deceptive multimodal landscapes (e.g. the 10-D Michalewicz function,
whose `Sin[...]^20` ridges are near-flat elsewhere so the local post-polish cannot
rescue a weak basin) need far more search than the flat budget provides; the larger
default lets the plain call reach a good basin. This costs nothing on easy problems
— the convergence early-break stops each run as soon as the feasible sub-population
collapses to tolerance, so low-dimensional and convex problems finish in the same
time as before. Naming `"DifferentialEvolution"` explicitly (or supplying your own
`SearchPoints` / `MaxIterations`) selects the flat budget, so a seeded run stays
bit-for-bit reproducible.

A method may carry sub-options as a list, e.g.
`Method -> {"DifferentialEvolution", "SearchPoints" -> 30,
"ScalingFactor" -> 0.6, "CrossProbability" -> 0.9, "RandomSeed" -> 7}`.
The search uses a deterministic PRNG with a fixed default seed, so results are
reproducible; `"RandomSeed"` overrides it.

Recognised sub-options:

| Sub-option | Applies to | Meaning |
|------------|-----------|---------|
| `"SearchPoints" -> n` | DE, NelderMead (restarts), RandomSearch (starts), SimulatedAnnealing (restarts) | population / restart / start count, honored verbatim (NelderMead and RandomSearch were silently capped at 20 / 40 before; an explicit value is now always run). Automatic defaults: DE population `Clip[10·d, {15, 40}]` under explicit `"DifferentialEvolution"` and `Clip[10·d, {15, 200}]` under `Method -> Automatic`; SimulatedAnnealing `Min[Max[2·d, 12], 50]` annealing chains; NelderMead `Min[2·d, 20]` simplex restarts; RandomSearch `Clip[8·d, {4, 40}]` starts |
| `"ScalingFactor" -> F` | DifferentialEvolution | DE differential weight (default 0.6) |
| `"CrossProbability" -> cr` | DifferentialEvolution | DE crossover probability (default 0.9) |
| `"PerturbationScale" -> s` | SimulatedAnnealing | multiplies the trial-step size (default 1.0); a positive real, an invalid value warns (`NMinimize::sopt`) and falls back to 1.0 |
| `"LevelIterations" -> L` | SimulatedAnnealing | trial moves at each temperature level, so the per-chain budget is `MaxIterations · L` (default 50, the previous fixed multiplier). An explicit value is honored verbatim — past the automatic aggregate cap, like `"SearchPoints"`. A positive integer, or `Automatic`/`None` for the default; an invalid value warns (`NMinimize::sopt`) and falls back to Automatic |
| `"BoltzmannExponent" -> f` | SimulatedAnnealing | acceptance-probability exponent for an uphill move: `f[i, df, f0]` (1-based iteration `i`, objective increase `df`, current value `f0`) sets the probability `Exp[f[i, df, f0]]`. `Automatic`/`None` keep the built-in cooling exponent `-df/(T·s)`, where `s` is a per-chain estimate of the objective scale (mean `|Δφ|` of a short start-point probe) so the acceptance temperature tracks the objective's magnitude; an invalid value warns (`NMinimize::bexp`) and falls back to Automatic |
| `"ReflectRatio" -> r` | NelderMead | simplex reflection coefficient (default 1) |
| `"ExpandRatio" -> e` | NelderMead | simplex expansion coefficient (default 2) |
| `"ContractRatio" -> c` | NelderMead | simplex contraction coefficient (default 0.5) |
| `"ShrinkRatio" -> s` | NelderMead | simplex shrink coefficient (default 0.5) |
| `"Tolerance" -> t` | NelderMead | simplex objective-spread convergence threshold |
| `"InitialPoints" -> {{x1,…},…}` | NelderMead | seed the initial simplex (extra points ignored, fewer are filled by perturbation; malformed → random start) |
| `"PostProcess" -> v` | all | final exact local polish. `True`/`Automatic`/a named method (`"InteriorPoint"`, `"FindMinimum"`, `"KKT"`, `{"…", opts}`) → on; `False`/`None` → off |
| `"PenaltyFunction" -> f` | all | function applied to each constraint's violation when scoring infeasible points. `Automatic`/`None` keep the built-in squared penalty; a pure function or function symbol (`#^2 &`, `(10 #) &`, `Sqrt`, …) replaces it |
| `"RandomSeed" -> s` | all | override the default PRNG seed |

Unrecognised sub-options are ignored (matching Mathematica). `"PostProcess"`
defaults to on: the global best is refined by the exact local optimizer (and,
for continuous box/unconstrained problems at `WorkingPrecision > MachinePrecision`,
by an MPFR BFGS step). It accepts the full Mathematica value set — `True`,
`Automatic`, or a named local method as a string (`"InteriorPoint"`,
`"FindMinimum"`, `"KKT"`, …) or `{"method", opts…}` turn the polish on; `False`
or `None` turn it off. Mathilda has a single `FindMinimum`-style local polish
(BFGS for continuous/box problems, a quadratic-penalty solver for general
constraints) that already picks the right inner solver, so a named method
enables post-processing rather than selecting a distinct algorithm; an
unrecognised value warns (`NMinimize::pmeth`) and falls back to `Automatic`.

`"InitialPoints"` seeds the NelderMead simplex, which matters for functions with
a narrow basin in a broad flat region — the Easom function
`-Cos[x] Cos[y] Exp[-((x-π)² + (y-π)²)]` is ≈0 everywhere except a spike of depth
−1 at `(π, π)`, so an undirected search over a large box rarely samples it.
NelderMead's convergence test additionally requires the simplex to be
geometrically small (not merely flat in objective value) before it stops, so a
simplex sitting on the plateau keeps contracting toward its centroid instead of
declaring victory immediately — enough to slide into the spike when the seed
points bracket it. With `"InitialPoints" -> {{-50,-50},{50,50},{10,10}}` (centroid
≈ `(π,π)`) NMinimize returns `{-1., {x -> π, y -> π}}`.

`"PenaltyFunction"` sets how infeasible points are scored during the global
search. By default (`Automatic`) a *general* (non-box) constraint contributes the
square of its violation — `Max[0, g]²` for an inequality `g ≤ 0`, `h²` for an
equality `h == 0` — and box bounds are handled by projection, not penalty. A
supplied function `f` replaces the square: a violated constraint contributes
`f[violation]`, so `Automatic` is exactly `#^2 &`. This affects only the
feasibility scoring used by the global search (Deb's rules, and the NelderMead /
SimulatedAnnealing penalized objective); the final local polish keeps the
differentiable squared penalty its analytic gradient assumes. A satisfied
inequality always contributes 0, so the feasibility test is unchanged; an invalid
value (a number, a string) warns (`NMinimize::penf`) and falls back to
`Automatic`. On a box-only problem the option is inert (there are no general
constraints to penalize).

`"SimulatedAnnealing"` honors four of its own sub-options.
`"SearchPoints" -> K` runs `K` independent annealing chains from random starts
and keeps the best local minimum (default `Automatic = Min[Max[2·d, 12], 50]`).
The `2·d` scaling follows Mathematica, but a rugged, many-basin landscape in low
dimension (Eggholder, Schwefel, Griewank, ...) has far more basins than `2·d = 4`
starts can cover, so a floor of 12 independent starts is what makes the reported
optimum reliably the global one there rather than luck-of-the-trajectory; a
bounded aggregate iteration budget is shared across the chains so a large `K`
stays fast. Each chain's best raw point is polished into its basin minimum before
the chains are ranked — as `"RandomSearch"` does per restart — so the reported
optimum *improves* monotonically with `"SearchPoints"` (ranking by random-walk
lows, which need not sit in the deepest basin, does not).
`"LevelIterations" -> L` sets the number of trial moves at each temperature
level, so a chain's total iteration budget is `MaxIterations · L` (default 50 —
the value the fixed multiplier used to hard-code, so an unset option reproduces
the previous schedule exactly). A larger `L` anneals each chain longer; it is
honored verbatim, past the automatic aggregate cap, since the caller has asked
for that budget. `"PerturbationScale" -> s` multiplies the size of the random
trial step (default 1.0) — a small scale keeps the walk local, a large one
explores more widely. `"BoltzmannExponent" -> f` replaces the Metropolis
acceptance exponent
for an uphill move: the point is accepted with probability `Exp[f[i, df, f0]]`,
where `i` is the (1-based) iteration, `df ≥ 0` the objective increase, and `f0`
the current objective; `Automatic`/`None` keep the built-in cooling exponent
`-df/(T·s)`. The temperature scale `s` is measured once per chain from a short
probe of trial steps at the start point (its mean `|Δφ|`): without it the fixed
`T ∈ [1e-4, 1]` sits far below an objective ranging over hundreds, every uphill
move is rejected, and the "anneal" degenerates into greedy descent from the start
point — so only many restarts, never the walk itself, ever crossed a ridge. The
probe draws from its own PRNG, so the annealing walk's stream is untouched and a
seeded `"SearchPoints" -> 1` chain follows exactly the trajectory it would
without the probe. A downhill move is always accepted, as usual. `"PostProcess"
-> False` skips the per-chain polish and returns the global raw best. On the
strongly-multimodal Griewank function `Sum[x[i]²/4000] - Product[Cos[x[i]/√i]] +
1` over `[-600, 600]^d`, the default multi-chain SimulatedAnnealing reaches
`~0.015` for `d = 10`, below Mathematica's `~0.175` on the same problem; on the
2-D Eggholder over `[-512, 512]²` it reaches the global `-959.64`, below
Mathematica's default `-935.34`.

### Options

| Option | Default | Meaning |
|--------|---------|---------|
| `Method` | `Automatic` | global-search selector (see above) |
| `WorkingPrecision` | `MachinePrecision` | MPFR refinement of continuous, general-constraint-free problems |
| `MaxIterations` | `Automatic` | DE generation cap: `150·d` under `Method -> Automatic`, `100` under explicit `"DifferentialEvolution"`; an explicit value overrides both |
| `AccuracyGoal` / `PrecisionGoal` | `Automatic` | convergence tolerances |
| `EvaluationMonitor` | `None` | `:> body` run on every objective evaluation |
| `StepMonitor` | `None` | `:> body` |

### Evaluation (auto-compilation and message quieting)

At `MachinePrecision` (the default) the global search evaluates the objective —
and each general constraint — at hundreds to thousands of trial points, so
`NMinimize` **auto-compiles** them to bytecode over the search variables once
(`compile_expr`, the same lowering `Compile[]` and NDSolve's RHS use) and runs
the register machine per point instead of the interpreter (`expr_copy` +
`evaluate` + `numericalize`). A body with a construct that cannot be compiled
stays on the interpreter, and every per-point call falls back to the interpreter
on a domain or non-finite result, so the compiled path is a pure speedup with no
change in answer. Example: the 10-D Rosenbrock over `Table[x[i], {i, 1, 10}]`
runs ~11× faster than the interpreter path. At `WorkingPrecision >
MachinePrecision` the exact MPFR interpreter path is used.

The same compiled objective also drives the **local polish** — the exact
solver that refines each candidate — not just the global-search scoring. The
polish's objective values, line searches, and gradient (taken by finite
differences off the compiled objective) all run on the register machine, with
the interpreter as fallback. This matters most for `"RandomSearch"`, which runs
one local solve per `SearchPoints`: a 20-variable `Sum[Abs[x[i]], {i, 1, 20}]`
with `"SearchPoints" -> 1000` went from ~25 s to ~0.5 s (same optimum). Every
method benefits (`DifferentialEvolution` and `NelderMead` also polish), so the
compiled program is now the single evaluation substrate across the whole
machine-precision optimizer.

Expected numeric-domain diagnostics raised while probing the function — e.g.
`Power::infy` from a `1/0` in a gradient term on a non-differentiable ridge
(Bukin N.6, `100 Sqrt[Abs[x2 - 0.01 x1^2]] + ...`) — are **quieted** during the
search: a non-finite value is treated as a bad point and steered away from, so
the messages are noise. This matches Mathematica, which quiets `NMinimize`'s
internal evaluation. (The same applies to `FindMinimum` / `FindMaximum`, which
share the point-evaluation path.)

### Examples

```
In[1]:= NMinimize[x^4 - 3 x^2 - x, x]
Out[1]= {-3.51391, {x -> 1.30084}}

In[2]:= NMinimize[{x + y, x^2 + y^2 <= 9}, {x, y}]
Out[2]= {-4.24264, {x -> -2.12132, y -> -2.12132}}

In[3]:= NMinimize[{x + 2 y, x^2 + 2 y^2 <= 3, x + y == 2, x >= 1}, {x, y}]
Out[3]= {2.33333, {x -> 1.66667, y -> 0.33333}}

In[4]:= NMinimize[{x + y, x + 2 y >= 3, x >= -2},
          {Element[x, Integers], Element[y, Integers]}]
Out[4]= {1, {x -> -1, y -> 2}}

In[5]:= NMinimize[{x, x > 2 && x < 1}, x]
Out[5]= {Infinity, {x -> Indeterminate}}

In[6]:= NMaximize[{x + y, x^2 + y^2 <= 1}, {x, y}]
Out[6]= {1.41421, {x -> 0.70710, y -> 0.70711}}
```

