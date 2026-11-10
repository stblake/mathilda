# Special Functions

Higher transcendental functions: the gamma function `Gamma`, the error function `Erf`, its complement `Erfc` and the imaginary error function `Erfi`, the digamma/polygamma family `PolyGamma`, the log-gamma function `LogGamma`, the harmonic numbers `HarmonicNumber`, the Pochhammer symbol (rising factorial) `Pochhammer`, the Riemann/Hurwitz zeta function `Zeta` (with the inert Stieltjes constants `StieltjesGamma`), the Hurwitz zeta function `HurwitzZeta`, the Bernoulli numbers and polynomials `BernoulliB`, the Euler numbers and polynomials `EulerE`, the polylogarithm `PolyLog`, the Lerch transcendent `LerchPhi`, the Dirichlet L-function `DirichletL`, the hypergeometric family `Hypergeometric0F1`, `Hypergeometric1F1`, `Hypergeometric2F1`, the generalized `HypergeometricPFQ`, and the Tricomi confluent function `HypergeometricU`, the Airy functions `AiryAi` and `AiryBi`, the Lambert W function `ProductLog`, the Legendre functions of the first and second kind (and their associated forms) `LegendreP` and `LegendreQ`, the Jacobi polynomials `JacobiP`, the Riemann–Siegel functions `RiemannSiegelTheta` and `RiemannSiegelZ`, the Chebyshev polynomials of the first and second kind `ChebyshevT` and `ChebyshevU`, the Laguerre polynomials `LaguerreL`, the physicists' Hermite polynomials `HermiteH`, the Whittaker functions `WhittakerM` and `WhittakerW`, the parabolic cylinder function `ParabolicCylinderD`, the Gegenbauer polynomials `GegenbauerC`, the Jacobi theta functions `EllipticTheta` (with `EllipticThetaPrime`), the regularized incomplete beta function `BetaRegularized`, the Appell hypergeometric functions of two variables `AppellF1`, `AppellF2`, `AppellF3`, `AppellF4`, the Weierstrass elliptic functions `WeierstrassP`, `WeierstrassPPrime`, `WeierstrassZeta`, `WeierstrassSigma` (with `WeierstrassHalfPeriods` and `WeierstrassInvariants`), the Hankel functions `HankelH1`, `HankelH2` and the spherical Bessel functions `SphericalBesselJ`, `SphericalBesselY`, `SphericalHankelH1`, `SphericalHankelH2`, the Jacobi elliptic functions `JacobiSN` ... `JacobiSD` (with `JacobiAmplitude` and `InverseJacobiSN`, `InverseJacobiCN`, `InverseJacobiDN`), the even and odd Mathieu functions `MathieuC` and `MathieuS` (with `MathieuCPrime`, `MathieuSPrime`), and the Meijer G-function `MeijerG` (with a minimal `FunctionExpand`, and `MeijerGReduce`, which rewrites an expression in inactive G functions).

**Machine overflow/underflow** on a machine-precision argument never leaks a libm
`inf.0`/`nan`: `Gamma[256.] → 3.35…×10^504`, `Erfi[30.] → 1.38…×10^389`,
`Gamma[900., 1.] → 7.50…×10^2266`, `BesselK[0, 800.] → 1.63…×10^-349`,
`BernoulliB[300.]`, `Hyperfactorial[50.]`, `Hypergeometric1F1[1., 2., 900.]`, … all
promote to an extended-exponent machine number (matching `N[f[exact]]`), and a
magnitude past MPFR's range returns `Overflow[]`/`Indeterminate`. See
[`elementary-functions.md`](elementary-functions.md) and the 2026-09-21 changelog.

**Mixed exact and inexact arguments** (v0.477). A numeric builtin given an
inexact argument numericalizes its exact numeric ones to the lowest precision
among the inexact arguments, as `Plus` and `Times` do: `Zeta[2., Sqrt[2]]`,
`BesselJ[Sqrt[2], 2.]`, `Hypergeometric2F1[1., 2, 3, Sqrt[2]/3]` and
`Zeta[2.`30, Sqrt[2]]` (30 digits) all evaluate. The rule lives once, in the
evaluator (`numeric_mixed_args_call`, `src/numeric.c`): when a `NumericFunction`
builtin declines a call whose arguments are numbers, numeric quantities or Lists
of them, with at least one inexact and at least one irrational exact argument,
the builtin is called once more with the irrational ones numericalized. Integers
and rationals are left as they are (they are often a structural index, such as
`BesselJ`'s order or `ProductLog`'s branch), positions held by `NHoldFirst` /
`NHoldRest` are untouched, `Flat` heads keep their own contagion pass, and a
decline that emitted a message is not retried. A symbolic argument still blocks
it (`Zeta[2., x]`), and so does a user-defined `NumericFunction` head, as in the
reference. Deliberate difference: `HarmonicNumber[Pi/3, 1.36]` evaluates; the
reference leaves that mixed call unevaluated while evaluating the all-machine
one.

## Gamma


- `Gamma[z]` — the Euler gamma function Γ(z) = ∫₀^∞ tᶻ⁻¹ e⁻ᵗ dt.
- `Gamma[a, z]` — the upper incomplete gamma function Γ(a, z) = ∫_z^∞ tᵃ⁻¹ e⁻ᵗ dt.
- `Gamma[a, z0, z1]` — the generalized incomplete gamma Γ(a, z0) − Γ(a, z1).

**Attributes**: `Listable`, `NumericFunction`, `Protected`.

**Features**:
- Exact reductions for `Gamma[z]`:
  - Positive integers: `Gamma[n] = (n-1)!` (exact, with GMP BigInt for large `n`).
  - Non-positive integers are poles: `Gamma[0]`, `Gamma[-n]` → `ComplexInfinity`.
  - Half-integers reduce to rational multiples of `Sqrt[Pi]`, e.g.
    `Gamma[1/2] = Sqrt[Pi]`, `Gamma[5/2] = 3/4 Sqrt[Pi]`,
    `Gamma[-1/2] = -2 Sqrt[Pi]` (via the Factorial functional equation).
  - `Gamma[Infinity]` → `Infinity`, `Gamma[-Infinity]` → `Indeterminate`,
    `Gamma[ComplexInfinity]` → `ComplexInfinity`.
- Exact reductions for the incomplete form:
  - `Gamma[a, 0] = Gamma[a]`, `Gamma[a, Infinity] = 0`.
  - **Positive integer first argument** reduces to its finite closed form
    `Gamma[n, z] = (n-1)! e^-z Σ_{k=0}^{n-1} z^k/k!` for symbolic or exact `z`,
    e.g. `Gamma[1, z] = E^-z`, `Gamma[2, x] = (1 + x) E^-x`,
    `Gamma[3, x] = (2 + 2 x + x^2) E^-x`, and exact `Gamma[2, 3] = 4/E^3`.
- Numeric evaluation:
  - Machine-precision real → libm `tgamma`. Machine-precision complex →
    the Gamma-family double kernel `sfg_gamma_c`
    (`src/special_functions/sf_machine_gamma.c`), shared bit-for-bit by the
    scalar call, packed / `NDArray` maps and `Compile[]`, e.g.
    `Gamma[2.3 + I] = 0.719141 + 0.540614 I`. For `-5 ≤ Re z < 10`,
    `|Im z| < 19/4` the exact recurrence moves z into the strip `1 ≤ Re < 5/2`
    and a Taylor series of lnΓ about the nearest of 30 precomputed centres is
    summed (the constant term carried as a double-double); elsewhere shifted
    Stirling (Re w ≥ 7) with the exponent formed in double-double (Dekker
    products, an extra-precise log|w|), then ONE complex exponential. Accuracy:
    median 0.8 / max 3.5 ulp on `0.5 ≤ Re z ≤ 4, |Im z| ≤ 3`; median 3.6 ulp
    over `-10..30 × ±20i`, where the conditioning |z ψ(z)| itself reaches 10².
    (The former Lanczos g = 7 series was median 9 / max 53 ulp on the same box.)
    `Re z < -30` and overflow fall back to the Lanczos/reflection path.
  - Arbitrary precision (MPFR) real → `mpfr_gamma`, output precision tracking
    the input, e.g. `N[Gamma[22/10], 50]` and `Gamma[2.2`200]`.
  - **Arbitrary-precision complex** `Gamma[z]` → Spouge's approximation,
    whose coefficients are computed at runtime to the requested precision
    (reflection for `Re(z) < 1/2`), e.g.
    `N[Gamma[I], 50] = -0.15494982830181068512… − 0.49801566811835604271… I`.
  - Incomplete real (machine or MPFR) → `mpfr_gamma_inc`, e.g.
    `Gamma[1.5, 7.5] = 0.00160996`, `Gamma[1, 1.1, 2.2] = 0.222068`.
  - **Incomplete complex** `Gamma[a, z]` (machine or arbitrary precision) →
    a lower-incomplete series (`Re(z) < Re(a)+1`) or a Lentz continued
    fraction otherwise, e.g. `Gamma[2.0, 1 + I] = (2 + I) e^-(1+I)` and
    `N[Gamma[3/2, 2 + I], 30] = 0.160487401929263240… − 0.176588715957602346… I`.
- Derivatives: `D[Gamma[a, z], z] = -z^(a-1) E^-z` (chain rule on both
  arguments; the `a`-derivative is the generic `Derivative[1,0][Gamma][a,z]`).
  `D[Gamma[z], z] = Gamma[z] PolyGamma[0, z]`, so higher derivatives compose
  through `PolyGamma`, e.g.
  `D[Gamma[z], {z, 2}] = Gamma[z] PolyGamma[1, z] + Gamma[z] PolyGamma[0, z]^2`.
- All other arguments (e.g. `Gamma[1/3]`, `Gamma[x]`, `Gamma[a, z]`,
  exact non-integer `Gamma[3/2, z]`, exact complex `Gamma[3/2, I]`) stay
  unevaluated.

```mathematica
In[1]:= Gamma[7/2]
Out[1]= 15/8 Sqrt[Pi]

In[2]:= Gamma[2.3 + I]
Out[2]= 0.719141 + 0.540614*I

In[3]:= Gamma[3, x]
Out[3]= E^(-x) (2 + 2 x + x^2)
```

## Beta

- `Beta[a, b]` — the Euler beta function B(a, b) = Γ(a)Γ(b)/Γ(a+b)
  = ∫₀¹ tᵃ⁻¹(1−t)ᵇ⁻¹ dt.
- `Beta[z, a, b]` — the incomplete beta function B_z(a, b)
  = ∫₀ᶻ tᵃ⁻¹(1−t)ᵇ⁻¹ dt (branch cut along z < 0).
- `Beta[z0, z1, a, b]` — the generalized incomplete beta
  B(z0, z1; a, b) = Beta[z1, a, b] − Beta[z0, a, b].

Attributes: `Listable`, `NumericFunction`, `Protected`.
Note that, unlike `Gamma`'s incomplete form, the integration variable `z`
comes **first** in `Beta[z, a, b]`.

- Exact reductions for `Beta[a, b]` (reduced through the `Gamma` machinery):
  - When one argument is a positive integer `n`, the gamma ratio collapses to
    `Beta[n, b] = (n−1)!/Pochhammer[b, n]`, an exact rational for rational `b`,
    e.g. `Beta[2, 3] = 1/12`, `Beta[5, 4] = 1/280`, `Beta[3, 1/3] = 27/14`,
    `Beta[7/2, 2] = 4/63`.
  - Half-integer arguments fold to rational multiples of `Pi`,
    e.g. `Beta[5/2, 7/2] = 3/256 Pi`, `Beta[1/2, 1/2] = Pi`.
  - Other exact arguments keep the residual gamma form,
    e.g. `Beta[1/3, 1/3] = Gamma[1/3]^2/Gamma[2/3]`,
    `Beta[2/3, 4/3] = Gamma[2/3] Gamma[4/3]`.
  - Symmetric: `Beta[a, b] = Beta[b, a]`.
- Poles: when `a`, `b`, or `a+b` lands on a gamma pole (a non-positive
  integer), a surviving pole gives `ComplexInfinity`
  (`Beta[0, b]`, `Beta[a, 0]`, `Beta[0, n]`, `Beta[-2, 5]`,
  `Beta[Infinity, 0]`); a cancelling pair reduces by the finite limit of the
  gamma ratio (`Beta[-2, 1] = -1/2`, `Beta[2, -5] = 1/20`).
- **Machine precision** (real or complex) is computed directly in double by
  `sfg_beta` / `sfg_beta_c` (`sf_machine_gamma.c`), the kernel the packed /
  `NDArray` map and `Compile[]` share — no Gamma-ratio expression is rebuilt.
  Real: `tgamma(a) / tgamma(a+b) · tgamma(b)` with the rounding of `a + b`
  undone to first order through ψ(a+b) (`Γ(s+e) = Γ(s)(1 + ψ(s) e)`; left
  alone it costs up to ~10 ulp at a+b ≈ 10). Complex: the three Gammas
  combined under one exponential (their double-double exponents summed), the
  same ψ correction. Accuracy median 0.8 / max 4.1 ulp real on `[0.3, 5]²`,
  median 1.3 / max 6.3 ulp complex on `[0.5,3] × ±2i`. Poles, `|a|, |b|, |a+b|
  > 170` and underflow fall through to the layered path below.
- Numeric `Beta[a, b]` otherwise: arbitrary-precision (MPFR) real and complex
  arguments evaluate numerically with precision tracking,
  e.g. `Beta[2.3, 3.2] = 0.0540298`, `Beta[2.5 + I, 1 − I] = 0.0831078 +
  0.142164 I`, `N[Beta[22/10, 33/10], 50] =
  0.056485691373282566807051754004491429369537777015241`.
- Incomplete `Beta[z, a, b]` reduces through `Hypergeometric2F1` via
  `B_z(a, b) = zᵃ/a · ₂F₁(a, 1−b; a+1; z)`: `Beta[0, a, b] = 0`,
  `Beta[1, a, b] = Beta[a, b]`, `Beta[z, a, 1] = z^a/a`; a positive-integer `b`
  terminates the series to an exact closed form (`Beta[1/2, 3, 4] = 7/640`),
  and a fully numeric call evaluates to a number. The four-argument form
  evaluates as the difference when all arguments are numeric.
- Derivatives:
  - `D[Beta[a, b], a] = Beta[a, b] (PolyGamma[0, a] − PolyGamma[0, a+b])`
    (symmetric in `b`); higher derivatives compose through `PolyGamma`.
  - `D[Beta[z, a, b], z] = z^(a−1) (1−z)^(b−1)`; the derivatives with respect
    to `a`, `b` are the inert `Derivative[0,1,0]/[0,0,1][Beta][z, a, b]`.
  - `D[Beta[z0, z1, a, b], z1] = z1^(a−1)(1−z1)^(b−1)`,
    `D[…, z0] = −z0^(a−1)(1−z0)^(b−1)`.
- Symbolic arguments (`Beta[a, b]`, `Beta[z, a, b]`, `Beta[z0, z1, a, b]`)
  stay unevaluated.

> **Packed arrays.** `Beta[a, b]` runs elementwise over buffers, including
> the case where **both** arguments are arrays. The machine kernel sits about
> one ulp from the arbitrary-precision path, as every machine special-function
> kernel does.

```mathematica
In[1]:= Beta[3, 5]
Out[1]= 1/105

In[2]:= Beta[-2, 1]
Out[2]= -1/2

In[3]:= Beta[2.5 + I, 1 - I]
Out[3]= 0.0831078 + 0.142164*I
```

## Erf

- `Erf[z]` — the error function erf(z) = (2/√π) ∫₀^z e^(−t²) dt.
- `Erf[z0, z1]` — the generalized error function erf(z1) − erf(z0).

**Attributes**: `Listable`, `NumericFunction`, `Protected`.

`Erf` is an entire function (no branch cuts) and odd in z.

**Features**:
- Exact special values: `Erf[0] = 0`, `Erf[Infinity] = 1`,
  `Erf[-Infinity] = -1`, `Erf[I Infinity] = DirectedInfinity[I]`,
  `Erf[-I Infinity] = DirectedInfinity[-I]`, plus `ComplexInfinity` and
  `Indeterminate` pass through.
- Odd symmetry for symbolic arguments: `Erf[-x] = -Erf[x]`,
  `Erf[-2 x] = -Erf[2 x]`.
- Numeric evaluation:
  - **Machine precision** (real and complex) → the raw-double kernels of
    `src/special_functions/sf_machine_error.c`, shared bit-for-bit by the
    scalar builtin, the packed/`NDArray` element-wise map and `Compile[]`
    (scalar and rank-1, real and complex). The machine gate runs before any
    special-value ladder.
    - Real: |x| < 1/2 the Maclaurin series; otherwise `erf = 1 − e^(−x²)·erfcx(x)`
      with erfcx from 0.5-wide polynomial tables (6 ≤ x ≤ 27.4: x·erfcx as a
      polynomial in 1/x²), e^(−x²) from the exact (Dekker) square — ≤ 1.6 ulp,
      slightly faster than libm. E.g. `Erf[0.95] = 0.820891`, `Erf[1.5] = 0.966105`.
    - Complex: |z| ≤ 0.8 the Maclaurin series (19 terms, Estrin); next to an
      axis (2|xy| < 1/2) the Taylor series about that axis (Hermite recurrence),
      so the SMALL component keeps full relative accuracy
      (`Erf[1.*^-300 + 2. I] = 6.16074*10^-299 + 18.5648 I`, not `1.1*10^-16 + …`);
      |Re z|, |Im z| < 4 `erf = 1 − e^(−z²) w(iz)` with the Faddeeva function w
      from a 0.5×0.5 grid of degree-20 Taylor expansions (generated at 60 digits)
      and e^(−z²) from an exact exponent and exactly reduced phase; elsewhere
      Salzer's series (A&S 7.1.29, an exponentially convergent trapezoidal sum,
      error 2e^(−4π²) relative), with `cosh(nY)` generated outward from the
      Gaussian peak and compensated sums. Odd/conjugate symmetry fold the plane.
      Measured against mpmath: median 0.4 ulp, ≤ 3.5 ulp away from the complex
      zeros (absolute error scaled by |z erf′(z)| ≤ 2.7 ulp everywhere).
      `|Im z| > 26` (e^(y²) overflow) declines to the MPFR path.
  - Arbitrary precision (MPFR) real → `mpfr_erf`, output precision tracking the
    input, e.g.
    `N[Erf[3/2], 50] = 0.96610514647531072706697626164594785868141047925764`
    and `Erf[0.95`100]`.
  - **Complex, arbitrary precision** → the cancellation-aware
    Maclaurin series erf(z) = (2/√π) e^(−z²) Σ 2ⁿ z^(2n+1)/(1·3···(2n+1))
    (DLMF 7.6.2), evaluated in MPFR with `|z|²/ln2` guard bits, e.g.
    `Erf[1.5 - I] = 1.0784 + 0.0279637 I`,
    `N[Erf[1/2 + I], 30] = 1.20484755831421800270211268210 + 1.02440088160844588172486045441 I`.
    A `double complex` series is the fallback for `USE_MPFR=0` builds.
  - **Large complex |z| (> 20)** → arb's `acb_hypgeom_erf` (likewise
    `Erfc`/`Erfi`): the series would need |z|²/ln2 guard bits and ~e|z|² terms
    (16 s at |z| = 100). A machine call whose value leaves the double range
    returns the extended-exponent value, as `Erfc`/`Erfi` already did:
    `Erf[3. + 30. I] = -1.248977832667308×10³⁸⁵ - 1.143536811453673×10³⁸⁵ I`
    (formerly `Overflow[]`), `Erf[0.5 + 100. I] = -1.976×10⁴³⁴⁰ + 3.327×10⁴³⁴⁰ I`
    (formerly unevaluated: the series' term cap was 2|z|², now 3|z|²).
- The two-argument form `Erf[z0, z1]` reduces to `erf(z1) − erf(z0)` only when
  both reduce to something concrete, e.g. `Erf[1.5, 2] = 0.0292171`,
  `Erf[-Infinity, Infinity] = 2`; exact/symbolic pairs such as `Erf[2, 3]` and
  `Erf[a, b]` stay unevaluated.
- Derivative: `D[Erf[z], z] = (2/Sqrt[Pi]) E^(−z²)` (chain rule applies), so
  the origin Taylor series follows from the generic `D`-based fallback, e.g.
  `Series[Erf[x], {x, 0, 5}]` begins `2/Sqrt[Pi] x − 2/(3 Sqrt[Pi]) x^3 + …`.
- All other arguments (symbolic `Erf[x]`, exact `Erf[2]`) stay unevaluated.

```mathematica
In[1]:= Erf[0.95]
Out[1]= 0.820891

In[2]:= Erf[1.*^-300 + 2. I]
Out[2]= 6.16074e-299 + 18.5648*I

In[3]:= Erf[1.5, 2]
Out[3]= 0.0292171
```

## Erfc

- `Erfc[z]` — the complementary error function erfc(z) = 1 − erf(z).

**Attributes**: `Listable`, `NumericFunction`, `Protected`.

`Erfc` is an entire function. Unlike `Erf` it is not odd, so `Erfc[-x]` is left
unexpanded (mathematically erfc(−x) = 2 − erfc(x)).

**Features**:
- Exact special values: `Erfc[0] = 1`, `Erfc[Infinity] = 0`,
  `Erfc[-Infinity] = 2`, `Erfc[I Infinity] = DirectedInfinity[-I]`,
  `Erfc[-I Infinity] = DirectedInfinity[I]` (negated relative to `Erf`), plus
  `ComplexInfinity` and `Indeterminate` pass through.
- Numeric evaluation:
  - **Machine precision** (real and complex) → the shared kernels of
    `sf_machine_error.c` (see `Erf`): real `erfc = e^(−x²)·erfcx(x)` for x ≥ 1/2
    (no cancellation; ≤ 3.1 ulp where the macOS libm measured 5.1), `2 − erfc(−x)`
    for x < 0; complex `erfc = e^(−z²) w(iz)` for Re z ≥ 0 inside the table box,
    `erfc(X) − (erf(z) − erf(X))` from Salzer's series beyond it (Re z > 26
    declines), `1 − erf(z)` for Re z < 0. E.g. `Erfc[0.95] = 0.179109`,
    `Erfc[1.5] = 0.0338949`, `Erfc[1.5 - I] = -0.0783992 - 0.0279637 I`.
  - Arbitrary precision (MPFR) real → `mpfr_erfc`, which is cancellation-free
    even for large positive z (where `1 − erf(z)` would lose all significance),
    output precision tracking the input, e.g.
    `N[Erfc[3/2], 50] = 0.033894853524689272933023738354052141318589520742363`.
  - **Complex, arbitrary precision** → `1 − erf(z)`, with erf(z)
    from the cancellation-aware DLMF 7.6.2 series evaluated in MPFR; the
    complement is formed at working precision (with `|z|²/ln2` guard bits) before
    rounding. A `double complex` series is the fallback for `USE_MPFR=0` builds.
- Derivative: `D[Erfc[z], z] = -(2/Sqrt[Pi]) E^(−z²)` (chain rule applies), so
  the origin Taylor series follows from the generic `D`-based fallback, e.g.
  `Series[Erfc[x], {x, 0, 3}]` begins `1 − 2/Sqrt[Pi] x + …`.
- All other arguments (symbolic `Erfc[x]`, exact `Erfc[2]`) stay unevaluated.

```mathematica
In[1]:= Erfc[10.]
Out[1]= 2.08849e-45

In[2]:= Erfc[1.5 - I]
Out[2]= -0.0783992 - 0.0279637*I

In[3]:= Erfc[-x]
Out[3]= Erfc[-x]
```

## Erfi

- `Erfi[z]` — the imaginary error function erfi(z) = erf(i z)/i
  = (2/√π) ∫₀^z e^(t²) dt.

**Attributes**: `Listable`, `NumericFunction`, `Protected`.

`Erfi` is an entire function (no branch cuts) and odd in z.

**Features**:
- Exact special values: `Erfi[0] = 0`, `Erfi[Infinity] = Infinity`,
  `Erfi[-Infinity] = -Infinity`. The imaginary-axis limits are **finite** (unlike
  `Erf`): `Erfi[I Infinity] = I`, `Erfi[-I Infinity] = -I`, since
  erfi(i y) = -i·erf(-y) → i as y → ∞. `ComplexInfinity` and `Indeterminate`
  pass through.
- Odd symmetry for symbolic arguments: `Erfi[-x] = -Erfi[x]`,
  `Erfi[-2 x] = -Erfi[2 x]`.
- Numeric evaluation (there is no libm/MPFR `erfi`, so the kernels are
  hand-rolled):
  - **Machine precision** → the shared kernels of `sf_machine_error.c`. Real:
    |x| < 1/2 the positive Maclaurin series; otherwise
    `erfi(x) = e^(x²)·g(x)` with g = e^(−x²) erfi(x) = (2/√π)·Dawson(x) from
    0.5-wide tables (x < 10) or its 14-term asymptotic series, e^(x²) from the
    exact square (scaled by e^(−64) where e^(x²) alone would overflow,
    x² ∈ (709, 713.6)); ≤ 2.6 ulp (the old summed series reached 30 ulps at
    |x| ≈ 5, and cost ~80 ns). Complex: `−i erf(iz)` through the complex erf
    kernel. `|x| > 26.72` declines to MPFR.
  - Real, arbitrary precision → the all-positive Maclaurin series
    erfi(x) = (2/√π) Σ x^(2n+1)/(n!(2n+1)). For real x every term shares x's
    sign, so the partial sums climb monotonically to the result — no
    cancellation, no e^(x²) prefactor — evaluated in MPFR with output precision
    tracking the input, e.g. `Erfi[2.5] = 130.396`, `Erfi[0.5] = 0.614952`,
    `N[Erfi[1/2], 50] = 0.61495209469651098083968118562364139305134561789540`.
  - **Complex, arbitrary precision** → erfi(z) = -i erf(i z),
    reusing the cancellation-aware erf series (DLMF 7.6.2) in MPFR with
    `|z|²/ln2` guard bits, e.g. `Erfi[1.5 - I] = -0.70136 - 1.84683 I`,
    `N[Erfi[1/2 + I], 30] = 0.187973467223383313628263810077 + 0.950709728318957173804611826379 I`.
    A `double complex` series is the fallback for `USE_MPFR=0` builds.
- Derivative: `D[Erfi[z], z] = (2/Sqrt[Pi]) E^(z^2)` (positive exponent, vs
  `Erf`'s E^(−z²); chain rule applies), so the origin Taylor series follows from
  the generic `D`-based fallback, e.g. `Series[Erfi[x], {x, 0, 7}]` begins
  `2/Sqrt[Pi] x + 2/(3 Sqrt[Pi]) x^3 + …`.
- All other arguments (symbolic `Erfi[x]`, exact `Erfi[2]`) stay unevaluated.

```mathematica
In[1]:= Erfi[I Infinity]
Out[1]= I

In[2]:= Erfi[2.5]
Out[2]= 130.396

In[3]:= Erfi[1.5 - I]
Out[3]= -0.70136 - 1.84683*I
```

## ExpIntegralEi

- `ExpIntegralEi[z]` — the exponential integral Ei(z), the principal value of
  −∫₋ᵤ^∞ e^(−t)/t dt (z = −u).

**Attributes**: `Listable`, `NumericFunction`, `Protected`.

`ExpIntegralEi` has a branch cut along the negative real axis (−∞, 0); on the cut
the principal value (the average of the two sides) is returned.

**Features**:
- Exact special values: `ExpIntegralEi[0] = -Infinity`,
  `ExpIntegralEi[Infinity] = Infinity`, `ExpIntegralEi[-Infinity] = 0`,
  `ExpIntegralEi[I Infinity] = I Pi`, `ExpIntegralEi[-I Infinity] = -I Pi`;
  `ComplexInfinity` and `Indeterminate` map to `Indeterminate`.
- Exact non-zero arguments stay symbolic (`ExpIntegralEi[2]`, `ExpIntegralEi[1/2]`);
  numeric values follow from a `Real`/MPFR argument or from `N`.
- Numeric evaluation (machine *and* arbitrary precision):
  - Real x > 0 → MPFR `mpfr_eint` (correctly rounded, fast even at very high
    precision: `ExpIntegralEi[1.8] = 4.24987`, `ExpIntegralEi[2.] = 4.95423`,
    `N[ExpIntegralEi[2], 50] = 4.9542343560018901633795051302270352755180535624200`).
  - Real x < 0 → the on-cut convergent series Ei(x) = γ + ln|x| + Σ xᵏ/(k·k!),
    in MPFR with `|x|/ln2` guard bits to absorb the partial-sum cancellation, and
    returns a **real** principal value (`ExpIntegralEi[-1.] = -0.219384`,
    `ExpIntegralEi[-5.] = -0.00114830`).
  - **Complex** → the same series with the principal `Log(z)`, evaluated in MPFR
    with `(|z| + |Re z|)/ln2` guard bits, so machine-precision complex results are
    fully accurate, e.g. `ExpIntegralEi[2. + I] = 4.06998 + 3.40094 I`,
    `N[ExpIntegralEi[2 + I], 30] = 4.06998094789392774228769025521 + 3.40094396980012162163040462603 I`.
    Approaching the cut from above gives +I Pi, from below −I Pi
    (`ExpIntegralEi[-1. + 10^-10 I] ≈ -0.219384 + 3.14159 I`). A `double complex`
    series is the `USE_MPFR=0` fallback.
  - **Large |z|** (real or complex) → the convergent series is infeasible
    (it would need ~`|z|/ln2` guard bits and ~`2|z|` terms), so once `|z|`
    exceeds roughly `prec·ln2` the **asymptotic expansion** takes over:
    `Ei(z) ~ (e^z/z) Σ_{k≥0} k!/z^k + i π sign(Im z)`, summed to its smallest
    term (DLMF 6.12.2; the `i π sign(Im z)` constant is the branch jump, so
    `ExpIntegralEi[±I Infinity] = ±I Pi` is recovered in the limit). The two
    regimes overlap and agree (`ExpIntegralEi[-50 I]` from either route gives
    `-0.00562839 − 3.12241 I`). This also fixes an `mpfr_init2` abort: the old
    guard-bit count `(long)(|z|/ln2)` overflowed for astronomically large `|z|`
    (e.g. `N[ExpIntegralEi[-10^60 I] + I Pi, 20]`).
- **Machine precision** (real and complex) runs on double kernels shared bit-for-bit
  with the packed/NDArray kernels and `Compile[]`
  (`src/special_functions/sf_machine_expint.c`; the scalar builtin calls the same
  function, so scalar == List == packed == NDArray == Compile). A decline (pole,
  overflow, an underflow to a subnormal) falls through to the MPFR path unchanged,
  and arbitrary precision is untouched. Real: Ei = γ + log x + x·R(x) for 0 < x ≤ 16
  (R on six pieces, each a positive-coefficient polynomial, so Horner never
  cancels), e^x/x·(1 + H(1/x)) rational on [16, 64], the asymptotic series beyond;
  Ei(−x) = −E1(x) with E1 from a polynomial for x ≤ 0.6 and e^−x/x·(1 − J(1/x))
  rationals beyond (Cody & Thacher's scheme, coefficients refitted with mpmath by
  `tools/sf_accuracy/fit_expint.py`). Complex: the power series by Horner where its
  terms do not cancel (|z| − Re z ≤ 2), else −E1(−z) ± iπ with E1 from the continued
  fraction evaluated backward from a calibrated depth, and the asymptotic series for
  |z| > 45. Accuracy vs mpmath: median 0.4 ulp (real) / 0.6 (complex), max ≤ 2.6 ulp
  off the zero x₀ = 0.3725 (scaled by |x f′| there).
- Derivative: `D[ExpIntegralEi[z], z] = E^z/z` (chain rule applies, e.g.
  `D[ExpIntegralEi[x^2], x] = (2 E^x^2)/x`); the origin Taylor series at a regular
  point follows from the generic `D`-based fallback.
- Wrong arity emits `ExpIntegralEi::argx` and stays unevaluated.

```mathematica
In[1]:= ExpIntegralEi[2.]
Out[1]= 4.95423

In[2]:= ExpIntegralEi[-1.]
Out[2]= -0.219384

In[3]:= ExpIntegralEi[2. + I]
Out[3]= 4.06998 + 3.40094*I
```

## ExpIntegralE

- `ExpIntegralE[n, z]` — the generalized exponential integral
  E_n(z) = ∫₁^∞ e^(−z t)/t^n dt (Re z > 0), with E_1(z) = Γ(0, z).

**Attributes**: `Listable`, `NumericFunction`, `Protected`.

**Features**:
- Exact values: `ExpIntegralE[0, z] = E^(−z)/z` (the derivative base case),
  `ExpIntegralE[n, 0] = 1/(n−1)` for integer n ≥ 2.
- Numeric: an integer order n ≥ 1 with an inexact real z > 0 evaluates by the
  standard continued-fraction (z > 1) / power-series (0 < z ≤ 1) split
  (`ExpIntegralE[1, 2.] = 0.0489005`, `ExpIntegralE[2, 1.5] = 0.0731008`).
  Exact arguments stay symbolic; `N` forces evaluation.
- **Machine precision** (integer order n ≥ 0): double kernels shared bit-for-bit
  with the packed/NDArray kernels and `Compile[]`
  (`src/special_functions/sf_machine_expint.c`; the scalar builtin calls the same
  function, so scalar == List == packed == NDArray == Compile). A decline (a pole,
  a non-integer order) falls through unchanged. Real x > 0: the upward recurrence E_(k+1) =
  (e^−x − x E_k)/k from E1 for x ≤ 1 (the power series for n > 24), and for x > 1
  the continued fraction (DLMF 8.19.17) evaluated backward from depth 105/x + 12
  (forward Lentz and Wallis forms lose 1e-14 to accumulated rounding). Machine
  **complex** z — previously unevaluated — and a real x < 0 (complex-valued, the
  upper side of the cut, e.g. `ExpIntegralE[2, -1.] = 0.823164 - 3.14159 I`) use the
  complex continued fraction for |z| + Re z > 1 and the series elsewhere. A
  machine-real underflow returns 0., as before (there is no arbitrary-precision path).
  `ExpIntegralE[n, packed]` now has a binary NDArray kernel (it threaded element by
  element before). Accuracy vs mpmath: median 0.4 ulp, max 4.3 ulp (real); max 3.8
  ulp (complex).
- Derivative: `D[ExpIntegralE[n, z], z] = -ExpIntegralE[n-1, z]`.
- Mellin transform: `∫₀^∞ x^(s-1) ExpIntegralE[n, a x] dx = a^(-s) Γ(s)/(s+n-1)`
  (`Re s > 0`), closed by the operational-calculus layer of `Method -> "Mellin"`.

```mathematica
In[1]:= ExpIntegralE[0, z]
Out[1]= E^(-z)/z

In[2]:= {ExpIntegralE[2, 3/2], N[ExpIntegralE[2, 3/2]]}
Out[2]= {ExpIntegralE[2, 3/2], 0.0731008}

In[3]:= ExpIntegralE[3, -1.]
Out[3]= 1.77072 - 1.5708*I

In[4]:= D[ExpIntegralE[n, z], z]
Out[4]= -ExpIntegralE[-1 + n, z]
```

- **Any order** (v0.477): a numeric call with an inexact argument that the
  integer-order kernels do not cover -- a real or complex order, or arbitrary
  precision -- evaluates through arb's `acb_hypgeom_expint` on the principal
  branch (`ExpIntegralE[Pi/3, 0.5] = 0.54325…`, `N[ExpIntegralE[1/3, 2], 25]`).
  `E_s(0.)` for `s ≤ 1` is `ComplexInfinity`.

## StruveH

- `StruveH[nu, z]` — the Struve function H_ν(z), the particular solution of the
  inhomogeneous Bessel equation with forcing `(z/2)^(ν+1)`, with series
  H_ν(z) = (z/2)^(ν+1) (2/(√π Γ(ν+3/2))) ₁F₂(1; 3/2, ν+3/2; −z²/4).

**Attributes**: `Listable`, `NumericFunction`, `Protected`.

**Features**:
- Exact value: `StruveH[nu, 0] = 0` for a real ν > −1.
- **Half-odd-integer order → elementary** (any exact or symbolic `z`), as the
  reference reduces it: `H_{−(n+½)}(z) = (−1)ⁿ J_{n+½}(z)` and
  `H_{n+½}(z) = Y_{n+½}(z) + (1/π) Σ_{j=0}^n Γ(j+½)(z/2)^(n−2j−½)/(n−j)!`
  (DLMF 11.4.4/11.4.5), with the half-integer Bessel functions written out in
  their spherical (Hankel-sum) form, so the cost is linear in n and no Bessel
  rewrite is involved: `StruveH[1/2, x] = Sqrt[2]/(Sqrt[Pi] Sqrt[x]) -
  (Sqrt[2] Cos[x])/(Sqrt[Pi] Sqrt[x])`, `StruveH[1/2, 2] = 1/Sqrt[Pi] -
  Cos[2]/Sqrt[Pi]`; `StruveH[101/2, x]` in 6 ms. Other rational orders stay
  symbolic.
- **Derivative.** `D[StruveH[ν, z], z] = ((z/2)^ν/(√π Γ(ν+3/2)) + H_{ν−1}(z) −
  H_{ν+1}(z))/2`; the ν-derivative is the inert `Derivative[1, 0][StruveH]`.
  `D[StruveL[ν, z], z]` gets the analogue `((z/2)^ν/(√π Γ(ν+3/2)) + L_{ν−1}(z) +
  L_{ν+1}(z))/2` (see StruveL below).
- **Arbitrary precision** (`N[StruveH[ν, z], p]` with an exact or MPFR argument):
  the rigorous FLINT/arb kernel `H_ν(z) = (z/2)^(ν+1) · ₁F̃₂(1; 3/2, ν+3/2; −z²/4)`
  (the *regularized* ₁F₂, so it is entire in ν — no pole to dodge — and valid for
  all ν over the whole complex plane). The working precision carries a
  cancellation allowance (~1.44·|z| bits) so a large |z| is resolved without
  blind doubling; the ball ladder renders only correct digits. Example:
  `N[StruveH[1, 13/10], 50] = 0.32012314597642023952394039617728160520311127053816`
  (matches mpmath).
- **Machine precision**: a fast real `double` kernel — the convergent ₁F₂-form
  series for small/moderate |z|, the `Y_ν + K_ν` asymptotic (DLMF 11.6.1 +
  10.17.4) for large |z| — shared bit-for-bit with the packed/NDArray kernel and
  `Compile[]` (`StruveH[0, 1.5] = 0.736723`, matching the integral representation
  (2/π)∫₀^(π/2) sin(z cos θ) dθ). Accurate to the double floor everywhere; in the
  narrow transition band z ∈ [~17, ~24] (small ν) both elementary series are
  intrinsically ill-conditioned and the result degrades to ~1e-9 absolute (as it
  does for every double implementation, scipy included) — the arb path above
  gives full precision there. A half-integer order makes both Hankel-type
  series finite, so they are summed whole: `ν = −(n+½)` is exact to rounding at
  every `z`, where `H_{−(n+½)}(z) = (−1)ⁿ J_{n+½}(z)` (DLMF 11.4.4) — the
  convergent series drops the terms killed by the `Γ(ν+3/2)` poles and
  `K_ν ≡ 0` (`StruveH[-3/2, 2.] = -0.491294 = −J_{3/2}(2)`; this returned
  `0.0 HypergeometricPFQ[…]` previously). An argument the kernel declines —
  a complex `z`, a real `z < 0` with a complex value — goes to the arb kernel at
  53 bits and comes back as machine numbers
  (`StruveH[-3/2, 2. + I] = -0.646752 - 0.172081 I`).
- Mellin transform: `∫₀^∞ x^(s-1) StruveH[0, a x] dx =
  π 2^(s-1) a^(-s)/(Cos(πs/2) Γ(1−s/2)²)`, closed by `Method -> "Mellin"` through
  the ₁F₂ reduction.

```mathematica
In[1]:= StruveH[0, 1.5]
Out[1]= 0.736723

In[2]:= N[StruveH[1, 13/10], 30]
Out[2]= 0.3201231459764202395239403961771

In[3]:= StruveH[1, {0.5, 1., 2.}]
Out[3]= {0.0521737, 0.198457, 0.646764}

In[4]:= {StruveH[-3/2, 2.], StruveH[-5/2, 2.]}
Out[4]= {-0.491294, 0.223925}
```

## StruveL

- `StruveL[nu, z]` — the modified Struve function L_ν(z), the particular
  solution of `z² y'' + z y' − (z² + ν²) y = 4 (z/2)^(ν+1)/(√π Γ(ν+½))`, with
  series L_ν(z) = (z/2)^(ν+1) Σ_k (z/2)^(2k)/(Γ(k+3/2) Γ(k+ν+3/2)) (DLMF 11.2.2).

**Attributes**: `Listable`, `NumericFunction`, `Protected`.

**Features**:
- Exact value: `StruveL[nu, 0] = 0` for a real ν > −1 (a symbolic ν stays).
- **Half-odd-integer order → elementary**, in `Sinh`/`Cosh`:
  `L_{−(n+½)}(z) = I_{n+½}(z)` and `L_{n+½}(z) = I_{−(n+½)}(z) − (1/π)
  Σ_{j=0}^n (−1)^j Γ(j+½)(z/2)^(n−2j−½)/(n−j)!` (DLMF 11.4.4 and the
  terminating 11.6.2), with the half-integer modified Bessel functions in their
  spherical (Hankel-sum) form: `StruveL[1/2, z] = -Sqrt[2]/(Sqrt[Pi] Sqrt[z]) +
  (Sqrt[2] Cosh[z])/(Sqrt[Pi] Sqrt[z])`, `StruveL[-1/2, z] = Sqrt[2] Sinh[z]/(Sqrt[Pi] Sqrt[z])`.
- **Derivative**: `D[StruveL[ν, z], z] = ((z/2)^ν/(√π Γ(ν+3/2)) + L_{ν−1}(z) +
  L_{ν+1}(z))/2`; the ν-derivative is the inert `Derivative[1, 0][StruveL]`.
- **Machine precision**: a raw-`double` kernel shared bit for bit by the scalar
  path, packed arrays / `NDArray` and `Compile[]` (scalar and rank-1, and so
  auto-compiled `Table`/`Sum`). For ν > −3/2 every series term has one sign, so
  the convergent series is cancellation-free at every real z (measured
  ≤ 1.2e-15 relative against 40-digit arb over ν ∈ [−2.3, 6], z ∈ (0, 100]);
  for z > 30 it switches to `I_{−ν}(z)` plus the DLMF 11.6.2 asymptotic. It
  declines where the value is complex (z < 0 with a non-integer ν) or the
  series cancels (ν < −3/2), and the arb kernel answers at 53 bits
  (`StruveL[0.5, -2.] = 0. − 1.5584 I`).
- **Arbitrary precision / complex**: arb's regularized
  `(z/2)^(ν+1) ₁F̃₂(1; 3/2, ν+3/2; z²/4)`, entire in ν, whole complex plane, with
  StruveH's cancellation allowance and ball ladder.

```mathematica
In[1]:= StruveL[1, 2.]
Out[1]= 1.10276

In[2]:= N[StruveL[1/3, 7], 30]
Out[2]= 166.928466768913900120706714408
```

## AngerJ / WeberE

- `AngerJ[nu, z]` — the Anger function **J**_ν(z) = (1/π)∫₀^π cos(νθ − z sin θ) dθ.
- `WeberE[nu, z]` — the Weber function **E**_ν(z) = (1/π)∫₀^π sin(νθ − z sin θ) dθ.

(The three-argument associated forms are not implemented; they stay unevaluated.)

**Attributes**: `Listable`, `NumericFunction`, `Protected`.

**Features**:
- **Exact reductions** (each checked against the integral definition at 30
  digits in `tests/test_struvel_anger_kelvin.c`):
  - `AngerJ[n, z] = BesselJ[n, z]` for an integer n;
  - `AngerJ[nu, 0] = Sinc[nu Pi]` (evaluated for an exact numeric ν:
    `AngerJ[1/2, 0] = 2/Pi`);
  - `WeberE[n, z] = (1/π) Σ_{k=0}^{⌊(n−1)/2⌋} Γ(k+½)(z/2)^(n−2k−1)/Γ(n+½−k) −
    StruveH[n, z]` for n ≥ 0 (DLMF 11.10.15): `WeberE[0, z] = -StruveH[0, z]`,
    `WeberE[3, z] = 2/(3 Pi) + (2 z^2)/(15 Pi) - StruveH[3, z]`;
  - `WeberE[-n, z] = (-1)^n WeberE[n, z]` (θ → π − θ in the integral);
  - `WeberE[nu, 0] = nu Pi Sinc[nu Pi/2]^2/2` (= (1 − cos νπ)/(νπ), evaluated
    for an exact numeric ν).
- **Derivative**: `D[AngerJ[ν, z], z] = (AngerJ[ν−1, z] − AngerJ[ν+1, z])/2`,
  likewise for `WeberE` (DLMF 11.10.32-33); the ν-derivative is inert.
- **Machine precision** (raw-`double`, shared bit for bit with packed /
  `NDArray` / `Compile[]`): the S1/S2 power series of DLMF 11.10.8-9 for small
  |z| (accepted when its cancellation estimate is within 4 ulp), otherwise
  `J_ν(z) + (sin νπ/π) A(ν)` and `−Y_ν(z) − (A(−ν) + cos νπ A(ν))/π`, with
  `A(μ) = ∫₀^∞ e^(−μt − z sinh t) dt` (DLMF 11.10.4-6) integrated by a 24-point
  Gauss–Legendre rule on [0, T] — the integrand is entire and decays
  double-exponentially, so the rule converges geometrically. Negative z uses
  `J_ν(−z) = J_{−ν}(z)`, `E_ν(−z) = −E_{−ν}(z)`. Errors are judged against the
  local envelope (both functions oscillate through zeros); measured ≤ 2.1e-15 of
  the envelope against 40-digit arb. An integer-order `AngerJ` is the BesselJ
  kernel itself.
- **Arbitrary precision / complex**: two regularized ₁F₂'s in arb,
  `S1 = ₁F̃₂(1; 1+ν/2, 1−ν/2; −z²/4)` and `S2 = (z/2) ₁F̃₂(1; 3/2+ν/2, 3/2−ν/2;
  −z²/4)`, combined with `sin(νπ/2)`, `cos(νπ/2)`; entire in ν and z.

```mathematica
In[1]:= {AngerJ[2, z], WeberE[1, z], AngerJ[nu, 0]}
Out[1]= {BesselJ[2, z], 2/Pi - StruveH[1, z], Sinc[nu Pi]}

In[2]:= {AngerJ[1/3, 2.], WeberE[1/3, 2.]}
Out[2]= {0.55058, -0.565462}
```

## KelvinBer / KelvinBei / KelvinKer / KelvinKei

- `KelvinBer[nu, z]`, `KelvinBei[nu, z]`, `KelvinKer[nu, z]`, `KelvinKei[nu, z]`
  — the Kelvin functions, `ber_ν(x) + i bei_ν(x) = J_ν(x e^(3πi/4))` and
  `ker_ν(x) + i kei_ν(x) = e^(−νπi/2) K_ν(x e^(πi/4))` for x > 0 (DLMF 10.61),
  continued analytically with the only branch cut on the negative real axis
  (values on the cut are taken from above; Schwarz-reflection symmetric).
- `KelvinBer[z]` etc. are the order-0 forms and evaluate to `KelvinBer[0, z]`.

**Attributes**: `Listable`, `NumericFunction`, `Protected`.

**Features**:
- **Exact values**: `head[-n, z] = (-1)^n head[n, z]` for an integer n
  (DLMF 10.61.6); `KelvinBer[0, 0] = 1`, `KelvinBer[nu, 0] = 0` (ν > 0),
  `KelvinBei[nu, 0] = 0` (ν ≥ 0), `KelvinKer[0, 0] = Infinity`,
  `KelvinKei[0, 0] = -Pi/4`, `KelvinKer[1, 0] = KelvinKei[1, 0] = ComplexInfinity`.
- **Derivative** (DLMF 10.63.1, halved and recombined): with R = ber|ker and
  I = bei|kei, `R'_ν = (R_{ν+1} − R_{ν−1} + I_{ν+1} − I_{ν−1})/(2√2)` and
  `I'_ν = (I_{ν+1} − I_{ν−1} − R_{ν+1} + R_{ν−1})/(2√2)`; so
  `D[KelvinBer[z], z] = (2 KelvinBei[1, z] + 2 KelvinBer[1, z])/(2 Sqrt[2])`.
- **Machine precision** (raw-`double`, shared bit for bit with packed /
  `NDArray` / `Compile[]`, including the unary order-0 kernels behind
  `KelvinBer[x]`): the complex Bessel J/K kernels of `sf_machine_bessel.c` at the
  rotated argument. Measured ≤ 4.4e-15 of the envelope `|ber + i bei|` (resp.
  `|ker + i kei|`) against 40-digit arb for x ≤ 100 — the residual is the
  conditioning of the rotated argument, ~x/√2 ulp. A negative x declines (except
  ber/bei at an integer order: `ber_n(−x) = (−1)^n ber_n(x)`) and is answered by
  arb (`KelvinBer[1/2, -2.] = 0. − 0.307238 I`).
- **Arbitrary precision / complex**: ber/bei from two regularized ₀F₁'s at
  ±i z²/4 times `(z/2)^ν`; ker/kei from `K_ν(z e^(±πi/4))` with the
  DLMF 10.34.2 continuation `K_ν(w e^(±πi)) = e^(∓νπi) K_ν(w) ∓ πi I_ν(w)` where
  a rotated argument crosses K's cut. A real result (real ν, z > 0) is
  returned as a Real.

```mathematica
In[1]:= {KelvinBer[1.5], KelvinKei[0, 1.5]}
Out[1]= {0.921072, -0.331396}

In[2]:= N[KelvinKer[1/3, -2 + 3 I], 20]
Out[2]= -10.5906978812710808984 - 4.89677861808647898171*I
```

## LogIntegral

- `LogIntegral[z]` — the logarithmic integral li(z), the principal value of
  ∫₀^z dt/ln t.

**Attributes**: `Listable`, `NumericFunction`, `Protected`.

`LogIntegral` is computed from the identity **li(z) = Ei(Log z)**, reusing the
`ExpIntegralEi` numeric kernels; the principal `Log` supplies the ±i π jump, so
the branch cut runs along (−∞, +1).

**Features**:
- Exact special values: `LogIntegral[0] = 0`, `LogIntegral[1] = -Infinity`,
  `LogIntegral[Infinity] = Infinity`, `LogIntegral[Indeterminate] =
  Indeterminate`; `LogIntegral[ComplexInfinity]` stays unevaluated (the limit
  depends on the direction). An inexact zero is li(0) = 0 (v0.476:
  `LogIntegral[0.] = LogIntegral[-0.] = 0.`, `LogIntegral[0.`30] = 0.`, also in
  packed, `NDArray` and compiled code; it was `Indeterminate`).
- Exact non-special arguments stay symbolic (`LogIntegral[2]`,
  `LogIntegral[1/2]`); numeric values follow from a `Real`/MPFR argument or from
  `N`.
- Numeric evaluation (machine *and* arbitrary precision) routes through
  `ExpIntegralEi[Log[z]]`:
  - Real z > 1 (`Log z > 0`) → MPFR `mpfr_eint`, correctly rounded and fast even
    at very high precision: `LogIntegral[20.] = 9.9053`, `LogIntegral[2.] = 1.04516`,
    `N[LogIntegral[2], 50] = 1.0451637801174927848445888891946131365226155781512`.
  - Real 0 < z < 1 (`Log z < 0`) → the on-cut convergent series, returning a
    **real** principal value: `LogIntegral[0.5] = -0.378671`,
    `LogIntegral[1.2] = -0.933787`.
  - **Complex** (and real z < 0, whose principal `Log` is complex) → the complex
    series with guard bits, so machine-precision complex results are fully
    accurate, e.g. `LogIntegral[2. + I] = 1.41126 + 1.22471 I`,
    `N[Re[LogIntegral[2 + I]], 30] = 1.41125904201780100568439320706`.
- **Machine precision** (real and complex) no longer builds and evaluates
  `ExpIntegralEi[Log[z]]`: double kernels shared bit-for-bit with the packed/NDArray kernels and `Compile[]` (`src/special_functions/sf_machine_expint.c`; the scalar builtin calls the same function, so scalar == List == packed == NDArray == Compile). A decline (pole, overflow, an underflow to a subnormal) falls through to the MPFR path unchanged, and arbitrary precision is untouched. The real kernel evaluates
  Ei(log x) with e^(log x) replaced by x itself in the large-argument forms, so
  the rounding of log x is not amplified by an exponential; a real x < 0 returns
  the complex value directly. The rounding of log x would otherwise cost ~log x
  ulps (9 at x = 1e4) through the e^y in Ei(y); it is recovered as
  dy = (x − e^y)/x and added back as (x/y)·dy. Accuracy vs mpmath: median 0.4 ulp,
  max 1.9 ulp (real, off the zero μ = 1.4513); max 2.7 ulp (complex).
- Derivative: `D[LogIntegral[z], z] = 1/Log[z]` (chain rule applies, e.g.
  `D[LogIntegral[x^2], x] = (2 x)/Log[x^2]`); the Taylor series at a regular
  point follows from the generic `D`-based fallback.
- Wrong arity emits `LogIntegral::argx` and stays unevaluated.

## SinIntegral

- `SinIntegral[z]` — the sine integral Si(z) = ∫₀ᶻ Sin[t]/t dt.

`SinIntegral` is entire and odd, with no branch cuts. The convergent Maclaurin
series `Si(z) = Σ (−1)ᵏ z^(2k+1)/((2k+1)(2k+1)!)` drives moderate arguments; for
large `|z|` the asymptotic form `Si(z) = π/2 − cos(z) f(z) − sin(z) g(z)` is used
instead. Complex arguments run through the shared `ncpx` MPFR-complex toolkit,
with odd symmetry folding the left half-plane onto the right.

**Machine precision** (real and complex) runs on double kernels shared bit-for-bit
with the packed/NDArray kernels and `Compile[]`
(`src/special_functions/sf_machine_expint.c`; the scalar builtin calls the same
function, so scalar == List == packed == NDArray == Compile). A decline (pole,
overflow, an underflow to a subnormal) falls through to the MPFR path unchanged, and
arbitrary precision is untouched. Real (Cephes `sici` scheme, coefficients refitted
with mpmath): x·S(x²) for |x| ≤ 4, and π/2 − f cos x − g sin x beyond with x f = 1 −
w RF(w), x² g = 1 − w RG(w), w = 1/x², rational on [4, 8] and [8, ∞). Complex: the
Maclaurin series by Horner where its terms do not cancel (|z| − |Im z| ≤ 2 or |z| ≤
2.5), else Si = π/2 + (E1(iz) − E1(−iz))/(2i) from the backward continued fraction,
and the f/g asymptotic series for |z| > 45. Accuracy vs mpmath: median 0.3 ulp, max
1.8 ulp (real); max 1.5 ulp scaled by |z f′| near the complex zeros.

- Exact special values: `SinIntegral[0] = 0`, `SinIntegral[Infinity] = Pi/2`,
  `SinIntegral[-Infinity] = -Pi/2`, `SinIntegral[±I Infinity] = ±I Infinity`;
  `ComplexInfinity` and `Indeterminate` map to `Indeterminate`.
- Exact non-special arguments stay symbolic (`SinIntegral[2]`, `SinIntegral[x]`);
  odd symmetry pulls a leading negative out (`SinIntegral[-x] = -SinIntegral[x]`).
- Numeric evaluation at machine or arbitrary (MPFR) precision, tracking the
  input precision: `SinIntegral[2.8] = 1.8321`,
  `N[SinIntegral[2], 50] = 1.6054129768026948485767201481985889408485834223285`.
- Complex arguments are fully accurate, e.g. `SinIntegral[2.5 + I] = 1.99549 +
  0.222995 I`, and `SinIntegral[3. I] = 4.97344 I` (= I Shi(3)).
- Derivative: `D[SinIntegral[z], z] = Sinc[z]` (chain rule applies, e.g.
  `D[SinIntegral[x^2], x] = 2 x Sinc[x^2]`).
- Series: Taylor at the origin (`Series[SinIntegral[x], {x, 0, 7}] =
  x - x^3/18 + x^5/600 - x^7/35280 + O[x]^8`) and the trig-prefactored asymptotic
  expansion at Infinity (`Normal[Series[SinIntegral[x], {x, Infinity, 3}]] =
  Pi/2 - Sin[x]/x^2 + Cos[x] (-1/x + 2/x^3)`). The symbolic
  `DirectedInfinity[z]`-direction form is not yet produced.
- `Listable`, `NumericFunction`, `Protected`. Wrong arity emits
  `SinIntegral::argx` and stays unevaluated.

```mathematica
In[1]:= SinIntegral[2.8]
Out[1]= 1.8321

In[2]:= SinIntegral[2.5 + I]
Out[2]= 1.99549 + 0.222995*I

In[3]:= SinIntegral[-x]
Out[3]= -SinIntegral[x]
```

## CosIntegral

- `CosIntegral[z]` — the cosine integral Ci(z) = −∫_z^∞ Cos[t]/t dt.

`CosIntegral` has a logarithmic singularity at `0` and a branch cut along the
negative real axis `(−∞, 0]`; it is neither entire nor odd. The convergent series
`Ci(z) = EulerGamma + Log(z) + Σ_{k≥1} (−1)ᵏ z^(2k)/(2k (2k)!)` drives moderate
arguments (the principal `Log` supplies the ±i π jump across the cut, so no folding
is needed); for large `|z|` the asymptotic form `Ci(z) = sin(z) f(z) − cos(z) g(z)`
is used, with a piecewise-constant Stokes / reflection term restoring the principal
branch on and left of the imaginary axis. Complex arguments run through the shared
`ncpx` MPFR-complex toolkit.

**Machine precision** (real and complex) runs on double kernels shared bit-for-bit
with the packed/NDArray kernels and `Compile[]`
(`src/special_functions/sf_machine_expint.c`; the scalar builtin calls the same
function, so scalar == List == packed == NDArray == Compile). A decline (pole,
overflow, an underflow to a subnormal) falls through to the MPFR path unchanged, and
arbitrary precision is untouched. Real: γ + log x + x² C(x²) for x ≤ 2, (x −
x₀)·Q(x) through the first zero x₀ = 3.3842 on [2, 4] (the series form cancels
15-fold by x = 4), and the f/g rationals beyond; a real x < 0 returns the complex
value directly. Complex: as for `SinIntegral`, with Ci = −(E1(iz) + E1(−iz))/2 and
Ci(z) = Ci(−z) ± iπ on the left half-plane. Accuracy vs mpmath: median 0.45 ulp, max
3.5 ulp scaled by |x Ci′| (the relative error is unbounded only at the zeros);
complex: series only for |z| ≤ 1 or |z| − |Im z| ≤ 1 (the log cancels against
the sum near the real axis), max 1.8 ulp scaled.

- Exact special values: `CosIntegral[0] = -Infinity`, `CosIntegral[Infinity] = 0`,
  `CosIntegral[-Infinity] = I Pi`, `CosIntegral[±I Infinity] = Infinity`;
  `ComplexInfinity` and `Indeterminate` map to `Indeterminate`.
- Exact non-special arguments stay symbolic (`CosIntegral[2]`, `CosIntegral[x]`,
  `CosIntegral[-x]`); there is no odd/even fold.
- Numeric evaluation at machine or arbitrary (MPFR) precision, tracking the input
  precision: `CosIntegral[2.8] = 0.186488`,
  `N[CosIntegral[2], 50] = 0.42298082877486499569856515319825589413573775630619`.
- Negative real arguments lie on the branch cut and return the from-above value
  `Ci(|x|) + I Pi`: `CosIntegral[-2.] = 0.422981 + 3.14159 I`.
- Complex arguments are fully accurate, e.g. `CosIntegral[2.5 + I] = 0.331465 −
  0.388237 I`, and on the imaginary axis `CosIntegral[3. I] = 4.96039 + 1.5708 I`
  (= Chi(3) + I Pi/2).
- Derivative: `D[CosIntegral[z], z] = Cos[z]/z` (chain rule applies, e.g.
  `D[CosIntegral[x^2], x] = (2 Cos[x^2])/x`).
- Series: the logarithmic Taylor series at the origin
  (`Series[CosIntegral[x], {x, 0, 6}] = EulerGamma + Log[x] - x^2/4 + x^4/96 -
  x^6/4320 + O[x]^7`) and the trig-prefactored asymptotic expansion at Infinity
  (`Normal[Series[CosIntegral[x], {x, Infinity, 3}]] = -Cos[x]/x^2 + Sin[x] (1/x -
  2/x^3)`). Generic-point Taylor expansions follow from the derivative rule.
- `Listable`, `NumericFunction`, `Protected`. Wrong arity emits
  `CosIntegral::argx` and stays unevaluated.

```mathematica
In[1]:= CosIntegral[2.8]
Out[1]= 0.186488

In[2]:= CosIntegral[-2.]
Out[2]= 0.422981 + 3.14159*I

In[3]:= CosIntegral[3. I]
Out[3]= 4.96039 + 1.5708*I
```

## SinhIntegral

- `SinhIntegral[z]` — the hyperbolic sine integral Shi(z) = ∫₀ᶻ Sinh[t]/t dt.

`SinhIntegral` is entire and odd, with no branch cuts — the imaginary-axis
sibling of `SinIntegral` (Shi(z) = −i Si(i z)). The convergent Maclaurin series
is `SinIntegral`'s without the alternating sign (all terms positive), so on the
real axis there is no catastrophic cancellation; large |z| uses the
optimally-truncated asymptotic `cosh(z) F(z) + sinh(z) G(z)` plus a Stokes
constant off the real axis.

**Machine precision** (real and complex) runs on double kernels shared bit-for-bit
with the packed/NDArray kernels and `Compile[]`
(`src/special_functions/sf_machine_expint.c`; the scalar builtin calls the same
function, so scalar == List == packed == NDArray == Compile). A decline (pole,
overflow, an underflow to a subnormal) falls through to the MPFR path unchanged, and
arbitrary precision is untouched. Real (Cephes `shichi` scheme): x·Sh(x²) with
positive coefficients for |x| ≤ 4, (Ei(x) + E1(x))/2 beyond. Complex: Shi(z) = −i
Si(iz), an exact rotation of the `SinIntegral` kernel. Accuracy vs mpmath: median
0.4 ulp, max 2.3 ulp (real); max 2.4 ulp (complex).

- Exact special values: `SinhIntegral[0] = 0`, `SinhIntegral[±Infinity] =
  ±Infinity`, `SinhIntegral[±I Infinity] = ±I Pi/2`; `ComplexInfinity` and
  `Indeterminate` give `Indeterminate`.
- Exact non-special arguments stay symbolic (`SinhIntegral[2]`,
  `SinhIntegral[x]`); odd symmetry pulls a leading negative out
  (`SinhIntegral[-x] = -SinhIntegral[x]`).
- Machine and arbitrary-precision reals evaluate numerically, tracking input
  precision: `SinhIntegral[2.8] = 4.34808`,
  `N[SinhIntegral[2], 50] = 2.5015674333549756414733724827275423989162728736915`.
  Results that overflow a machine double (`SinhIntegral[10.^6] ≈ 1.5166×10^434288`)
  come back as an extended-exponent (MPFR) real rather than `Infinity`.
- Complex arguments are fully accurate, e.g. `SinhIntegral[2.5 + I] = 2.84649 +
  2.22177 I`, and `SinhIntegral[3. I] = 1.84865 I` (= I Si(3)).
- Derivative: `D[SinhIntegral[z], z] = Sinh[z]/z` (chain rule applies, e.g.
  `D[SinhIntegral[x^2], x] = (2 Sinh[x^2])/x`).
- Series: Taylor at the origin (`Series[SinhIntegral[x], {x, 0, 7}] =
  x + x^3/18 + x^5/600 + x^7/35280 + O[x]^8`) and the hyperbolic-prefactored
  asymptotic expansion at Infinity
  (`Normal[Series[SinhIntegral[x], {x, Infinity, 3}]] =
  -I Pi/2 + Sinh[x]/x^2 + Cosh[x] (1/x + 2/x^3)`).
- Integrate: `Integrate[Sinh[a x]/x, x] = SinhIntegral[a x]`.
- `Listable`, `NumericFunction`, `Protected`. Wrong arity emits
  `SinhIntegral::argx` and stays unevaluated.

```mathematica
In[1]:= SinhIntegral[2.8]
Out[1]= 4.34808

In[2]:= SinhIntegral[2.5 + I]
Out[2]= 2.84649 + 2.22177*I

In[3]:= SinhIntegral[-2 x]
Out[3]= -SinhIntegral[2 x]
```

## CoshIntegral

- `CoshIntegral[z]` — the hyperbolic cosine integral
  Chi(z) = EulerGamma + Log[z] + ∫₀ᶻ (Cosh[t] − 1)/t dt.

`CoshIntegral` has a logarithmic singularity at `0` and a branch cut along the
negative real axis `(−∞, 0]` — the imaginary-axis sibling of `CosIntegral`
(Chi(z) = Ci(i z) − i Pi/2). The convergent series is `CosIntegral`'s without
the alternating sign; large |z| uses the asymptotic `sinh(z) F(z) + cosh(z) G(z)`
plus a piecewise Stokes constant, with the principal `Log` supplying the branch.

**Machine precision** (real and complex) runs on double kernels shared bit-for-bit
with the packed/NDArray kernels and `Compile[]`
(`src/special_functions/sf_machine_expint.c`; the scalar builtin calls the same
function, so scalar == List == packed == NDArray == Compile). A decline (pole,
overflow, an underflow to a subnormal) falls through to the MPFR path unchanged, and
arbitrary precision is untouched. Real: γ + log x + x² Ch(x²) for x ≤ 4, (Ei(x) −
E1(x))/2 beyond; a real x < 0 returns the complex value directly. Complex: Chi(z) =
Ci(iz) + (log z − log iz), the constant (−iπ/2, or 3iπ/2 past arg z = π/2) taken
exactly rather than by subtraction. Accuracy vs mpmath: median 0.45 ulp, max 1.6 ulp
scaled (the zero 0.5238); complex max 2.6 ulp.

- Exact special values: `CoshIntegral[0] = -Infinity`, `CoshIntegral[Infinity] =
  Infinity`, `CoshIntegral[-Infinity] = Infinity`, `CoshIntegral[±I Infinity] =
  ±I Pi/2`; `ComplexInfinity` and `Indeterminate` give `Indeterminate`.
- Exact non-special arguments stay symbolic (`CoshIntegral[2]`, `CoshIntegral[x]`,
  `CoshIntegral[-x]`); there is no odd/even fold.
- Machine and arbitrary-precision reals evaluate numerically, tracking input
  precision: `CoshIntegral[2.8] = 4.33122`,
  `N[CoshIntegral[2], 50] = 2.4526669226469145219061326474994928766017806887285`.
  A negative real returns the from-above branch value `Chi(|x|) + I Pi`:
  `CoshIntegral[-2.] = 2.45267 + 3.14159 I`. Overflowing results
  (`CoshIntegral[10.^6]`) come back as an extended-exponent (MPFR) real.
- Complex arguments are fully accurate, e.g. `CoshIntegral[3.5 + I] = 5.36668 +
  4.28423 I`, and on the imaginary axis `CoshIntegral[3. I] = 0.11963 + 1.5708 I`
  (= Ci(3) + I Pi/2).
- Derivative: `D[CoshIntegral[z], z] = Cosh[z]/z` (chain rule applies, e.g.
  `D[CoshIntegral[x^2], x] = (2 Cosh[x^2])/x`).
- Series: the logarithmic Taylor series at the origin
  (`Series[CoshIntegral[x], {x, 0, 6}] = EulerGamma + Log[x] + x^2/4 + x^4/96 +
  x^6/4320 + O[x]^7`) and the hyperbolic-prefactored asymptotic expansion at
  Infinity (`Normal[Series[CoshIntegral[x], {x, Infinity, 3}]] =
  -I Pi/2 + Cosh[x]/x^2 + Sinh[x] (1/x + 2/x^3)`).
- Integrate: `Integrate[Cosh[a x]/x, x] = CoshIntegral[a x]`.
- `Listable`, `NumericFunction`, `Protected`. Wrong arity emits
  `CoshIntegral::argx` and stays unevaluated.

```mathematica
In[1]:= CoshIntegral[2.8]
Out[1]= 4.33122

In[2]:= CoshIntegral[-2.]
Out[2]= 2.45267 + 3.14159*I

In[3]:= CoshIntegral[3.5 + I]
Out[3]= 5.36668 + 4.28423*I
```

## FresnelC and FresnelS

- `FresnelC[z]` — the Fresnel integral C(z) = ∫₀ᶻ Cos[π t²/2] dt.
- `FresnelS[z]` — the Fresnel integral S(z) = ∫₀ᶻ Sin[π t²/2] dt.

Both are entire and odd, with no branch cuts. One numeric kernel computes the
pair (C, S) together: real inputs sum the two convergent Maclaurin series
directly (switching to the DLMF 7.12 asymptotic form for large `|x|`), while
complex inputs sum `A(z) = C + iS = Σ (iπ/2)ᵏ z^(2k+1)/(k!(2k+1))` and its
`i → −i` partner `B` through the shared `ncpx` MPFR-complex toolkit, then recover
`C = (A+B)/2`, `S = (A−B)/(2i)`. Guard bits (`≈ (π/2)|z|²/ln2`) absorb the
partial-sum cancellation exactly; odd symmetry folds negative/left-half-plane
inputs.

- Exact special values: `FresnelC[0] = FresnelS[0] = 0`;
  `FresnelC[±Infinity] = ±1/2`, `FresnelS[±Infinity] = ±1/2`;
  `FresnelC[±I Infinity] = ±I/2`, `FresnelS[±I Infinity] = ∓I/2` (note the sign
  flip on the imaginary directions); `ComplexInfinity`/`Indeterminate` map to
  `Indeterminate`.
- Exact non-special arguments stay symbolic; odd symmetry pulls a leading
  negative out (`FresnelC[-x] = -FresnelC[x]`).
- Numeric evaluation at machine or arbitrary (MPFR) precision, tracking the input
  precision: `FresnelC[1.8] = 0.333633`, `FresnelS[1.8] = 0.450939`,
  `N[FresnelC[2], 50] = 0.48825340607534075450022350335726103768836715450922`.
  **Machine precision** (real and complex; scalar, packed/`NDArray` and
  `Compile[]` share the kernels of `sf_machine_error.c` bit for bit):
  - Real: |x| ≤ 1 the Maclaurin series (x³ formed exactly); 1 < |x| < 8
    `C = 1/2 + f sin(πx²/2) − g cos(πx²/2)`, `S = 1/2 − f cos − g sin` (DLMF
    7.5.3–4) with the auxiliary f, g from 0.5-wide polynomial tables (the Cephes
    `fresnl` structure); |x| ≥ 8 the asymptotic series of f, g (DLMF 7.12.2–3).
    The phase πx²/2 is reduced **exactly** (x² as a Dekker double-double, mod 4),
    so even `FresnelC[123456.789]` is correctly rounded — libm's `sin(π x x/2)`
    carries ~x² ulps of phase error. ≤ 2.5 ulp, median 0.3.
  - Complex: |z| ≤ 1.3 the Maclaurin series in z⁴; beyond,
    `C ± iS = ((1±i)/2) erf((√π/2)(1∓i) z)` through the complex erf kernel, but
    with every exponentiated quantity (`X² = (π/4)(x+y)²`, `XY = (π/4)(y²−x²)`,
    `Y² − X² = ∓πxy`) formed from the exact x, y — rounding the rotated argument
    first costs ~π|z|² ulps (25 at |z| = 2.8). Next to an axis (π|xy| < 1/2) the
    Taylor series about the real axis (via C(iz) = iC(z), S(iz) = −iS(z) for the
    imaginary one), so the small component is accurate relative to itself.
    ≤ 3.5 ulp away from the zeros.
- Complex arguments are fully accurate, e.g. `FresnelC[2.5 + I] = 116.648 −
  105.229 I`, and on the imaginary axis `FresnelC[2. I] = 0. + 0.488253 I` (a
  machine-complex argument gives a machine-complex result, as in the reference).
- Derivatives: `D[FresnelC[z], z] = Cos[π z²/2]`, `D[FresnelS[z], z] =
  Sin[π z²/2]` (chain rule applies, e.g. `D[FresnelC[x^2], x] = 2 x Cos[π x⁴/2]`).
- Series: Taylor at the origin (`Series[FresnelC[x], {x, 0, 9}] = x − π² x⁵/40 +
  π⁴ x⁹/3456 + O[x]^10`; `Series[FresnelS[x], {x, 0, 7}] = π x³/6 − π³ x⁷/336 +
  O[x]^8`) and the trig-prefactored asymptotic expansion at Infinity
  (`Normal[Series[FresnelC[x], {x, Infinity, 3}]] = 1/2 − Cos[π x²/2]/(π² x³) +
  Sin[π x²/2]/(π x)`). The symbolic-index general term is produced as a
  `Piecewise` closed form, e.g. `SeriesCoefficient[FresnelC[x], {x, 0, n}] =
  Piecewise[{{(-1)^((n-1)/4) 2^(1-n) π^(n/2)/(n Γ[(1+n)/4] Γ[(3+n)/4]),
  Mod[n-1, 4] == 0 && n ≥ 1}}, 0]`.
- `Listable`, `NumericFunction`, `Protected`. Wrong arity emits
  `FresnelC::argx` / `FresnelS::argx` and stays unevaluated.

```mathematica
In[1]:= FresnelC[1.8]
Out[1]= 0.333633

In[2]:= FresnelS[1.8]
Out[2]= 0.450939

In[3]:= FresnelC[2.5 + I]
Out[3]= 116.648 - 105.229*I

In[4]:= FresnelS[-I Infinity]
Out[4]= 1/2*I
```

## Sinc

- `Sinc[z]` — the cardinal sine Sin[z]/z, with the removable singularity filled
  in as `Sinc[0] = 1`.

`Sinc` is entire and even. Machine-precision arguments use the kernels of
`sf_machine_error.c` (shared with the packed/`NDArray` map and `Compile[]`):
real `sin(x)/x` (≤ 1.9 ulp), complex `|z| ≤ 1` the even Maclaurin series and
`(sin x cosh y + i cos x sinh y)/z` beyond (≤ 3.7 ulp; `|Im z| > 709` declines).
Arbitrary-precision real arguments use `mpfr_sin(z)/z`; complex ones `sin(z)/z`
via the shared `ncpx` toolkit. Provided as a first-class head so `SinIntegral`
differentiates to it.

- Exact special values: `Sinc[0] = 1`, `Sinc[±Infinity] = 0`;
  `ComplexInfinity` maps to `Indeterminate`.
- Numeric evaluation at machine or arbitrary (MPFR) precision, tracking input
  precision: `Sinc[2.] = 0.454649`,
  `N[Sinc[2], 45] = 0.454648713412840847698009932955872421351127485`. Complex:
  `Sinc[1. + I] = 0.966711 - 0.331747 I`.
- Symbolic arguments stay unevaluated (`Sinc[x]`, `Sinc[-x]`; the reference does not
  auto-fold the even symmetry).
- Derivative: `D[Sinc[z], z] = Cos[z]/z - Sin[z]/z^2`.
- Series at the origin: `Series[Sinc[x], {x, 0, 6}] =
  1 - x^2/6 + x^4/120 - x^6/5040 + O[x]^7`.
- `Listable`, `NumericFunction`, `Protected`. Wrong arity emits `Sinc::argx`
  and stays unevaluated.

## InverseErf

- `InverseErf[s]` — the inverse error function: the z solving s = erf(z).
- `InverseErf[z0, s]` — the inverse of the generalized error function
  `Erf[z0, z]`, i.e. `InverseErf[s + Erf[z0]]`.

**Attributes**: `Listable`, `NumericFunction`, `Protected`.

`InverseErf` is odd in s and (per the reference) produces explicit numerical
values only for **real** s in `[-1, 1]`; complex and out-of-domain (`|s| > 1`)
inputs stay symbolic.

**Features**:
- Exact special values: `InverseErf[0] = 0`, `InverseErf[1] = Infinity`,
  `InverseErf[-1] = -Infinity` (likewise `InverseErf[1.] = Infinity`), plus
  `Indeterminate` passes through.
- Odd symmetry for symbolic arguments: `InverseErf[-x] = -InverseErf[x]`.
- Numeric evaluation (neither C99 libm nor MPFR ships an inverse-erf, so the
  kernel is Newton's iteration on `f(z) = erf(z) − s`,
  `z ← z − (erf(z) − s)(Sqrt[Pi]/2) e^{z²}`):
  - Machine-precision real → `sf_inverse_erf_r` (`sf_machine_error.c`, shared
    with the packed/`NDArray` map and `Compile[]`): |s| < 0.875 a degree-14
    polynomial for erfinv(s)/s on 1/16-wide pieces (no logarithm); beyond, in
    Giles' variable w = −log(1 − s²) (erfinv(s)/s is analytic in the strip
    |Im w| < π), unit-width degree-11 pieces for w < 6 and, past that, a 1e-7 seed
    in 1/√w plus one Halley step on erfc written with the relative residual. No
    erf evaluation on the main path; ≤ 1.5 ulp (the old Winitzki + up-to-four
    Newton steps cost an erf and an exp each). E.g. `InverseErf[0.6] = 0.595116`,
    `InverseErf[1/{2.,3.,4.,5.}] = {0.476936, 0.30457, 0.225312, 0.179143}`.
  - Arbitrary precision (MPFR) real → Newton with precision doubling on
    `mpfr_erf`, output precision tracking the input, e.g.
    `N[InverseErf[33/100], 50] = 0.30133214613370582612850271815839477396582428282853`.
- Two-argument form: `InverseErf[0.4, 0.2] = 0.631776`; if `Erf[z0]` does not
  reduce the call stays in two-argument form, while a reducible `z0` collapses to
  the one-argument inverse, e.g. `InverseErf[0, 1.3] = InverseErf[1.3]`.
- Derivative: `D[InverseErf[z], z] = (Sqrt[Pi]/2) E^(InverseErf[z]²)` (chain rule
  applies), so the origin Taylor series follows from the generic `D`-based
  fallback, e.g. `Series[InverseErf[x], {x, 0, 8}] = (Sqrt[Pi]/2) x +
  (Pi^(3/2)/24) x³ + (7 Pi^(5/2)/960) x⁵ + (127 Pi^(7/2)/80640) x⁷ + O[x]^9`.
- All other arguments (symbolic `InverseErf[x]`, exact `InverseErf[1/2]`,
  out-of-domain `InverseErf[2]`) stay unevaluated.

```mathematica
In[1]:= InverseErf[0.6]
Out[1]= 0.595116

In[2]:= InverseErf[0.4, 0.2]
Out[2]= 0.631776

In[3]:= InverseErf[1]
Out[3]= Infinity
```

## InverseErfc

- `InverseErfc[s]` — the inverse complementary error function: the z solving
  s = erfc(z).

**Attributes**: `Listable`, `NumericFunction`, `Protected`.

Since erfc decreases from 2 (at −∞) to 0 (at +∞), `InverseErfc` produces
explicit numerical values only for **real** s in `[0, 2]`; complex and
out-of-domain (`s < 0` or `s > 2`) inputs stay symbolic.

**Features**:
- Exact special values: `InverseErfc[0] = Infinity`, `InverseErfc[1] = 0`,
  `InverseErfc[2] = -Infinity` (likewise the real `0.`/`2.`), plus
  `Indeterminate` passes through.
- Numeric evaluation. Mathematically `InverseErfc[s] = InverseErf[1 − s]`, but
  for small s (large z) that subtraction loses all significance to
  cancellation, so the kernel instead Newton-iterates directly on the
  cancellation-free `erfc`: `f(z) = erfc(z) − s`,
  `z ← z + (erfc(z) − s)(Sqrt[Pi]/2) e^{z²}`:
  - Machine-precision real → `sf_inverse_erfc_r`, the `InverseErf` kernel fed
    the exact complement q = min(s, 2 − s), so the tail (q down to the least
    normal double, z ≈ 26.5) is resolved from q itself, never from 1 − s;
    ≤ 2.2 ulp. A subnormal s (whose few bits the Halley residual cannot resolve)
    is solved at 53 bits in MPFR — it used to answer `Indeterminate`. E.g.
    `InverseErfc[0.6] = 0.370807`, `InverseErfc[1/{2.,3.,4.,5.}] =
    {0.476936, 0.68407, 0.81342, 0.906194}`, `InverseErfc[1.*^-320] = 27.0732`.
  - Arbitrary precision (MPFR) real → Newton with precision doubling on
    `mpfr_erfc`, output precision tracking the input, e.g.
    `N[InverseErfc[33/100], 50] = 0.68880252811655645040250472890525783544948992349371`
    (≈0.002 s even at 500-digit precision).
- erfc is **not** odd, so unlike `InverseErf` there is no auto-applied symmetry
  rewrite; the reflection `InverseErfc[1.5] = -InverseErfc[0.5]` is computed
  numerically rather than rewritten.
- Derivative: `D[InverseErfc[z], z] = -(Sqrt[Pi]/2) E^(InverseErfc[z]²)` (chain
  rule applies), so higher derivatives and the Taylor series follow from the
  generic `D`-based fallback.
- All other arguments (symbolic `InverseErfc[x]`, exact `InverseErfc[1/2]`,
  out-of-domain `InverseErfc[2.3]`) stay unevaluated.

## HarmonicNumber

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

- `HarmonicNumber[n]` — the n-th harmonic number `H_n = Sum_{i=1}^n 1/i`.
- `HarmonicNumber[n, r]` — the order-r harmonic number
  `H_n^(r) = Sum_{i=1}^n 1/i^r`.

Behaviour:

- **Exact integer n.** A non-negative integer `n` expands to the explicit finite
  sum `Sum_{i=1}^n i^-r`: an exact rational for integer `r`
  (`HarmonicNumber[10] = 7381/2520`, `HarmonicNumber[5, 2] = 5269/3600`), and an
  explicit sum for symbolic order (`HarmonicNumber[4, r] = 1 + 2^-r + 3^-r +
  4^-r`). `HarmonicNumber[0, r] = 0`.
- **At infinity.** `HarmonicNumber[Infinity, r] = Zeta[r]`, so
  `HarmonicNumber[Infinity, 2] = Pi^2/6`, `HarmonicNumber[Infinity, 4] = Pi^4/90`,
  `HarmonicNumber[Infinity, 3] = Zeta[3]`; `HarmonicNumber[Infinity]` →
  `ComplexInfinity`.
- **Non-positive integer order.** `HarmonicNumber[n, -m]` is the Faulhaber
  polynomial `Sum_{i=1}^n i^m`, returned as a polynomial in `n` (built from
  `BernoulliB`): `HarmonicNumber[n, 0] = n`, `HarmonicNumber[n, -1] = n/2 + n^2/2`,
  `HarmonicNumber[z, -4] = -z/30 + z^3/3 + z^4/2 + z^5/5`.
- **Numeric.** With an inexact argument and a numericizable `n`, the value comes
  from the analytic identity `H_n^(r) = Zeta[r] - Zeta[r, n+1]` (and the digamma
  form `EulerGamma + PolyGamma[0, n+1]` for `r = 1`), evaluated at machine or
  arbitrary (`MPFR`) precision with precision tracked from the input:
  `HarmonicNumber[.8, 3] = 0.940124`, `HarmonicNumber[E, 1.] = 1.75002`,
  `N[HarmonicNumber[1/17, 5], 50] =
  0.25327615206118707521034626118754228313433140885976`. Complex order is
  supported: `N[HarmonicNumber[27, 5 - I]] = 1.02598 + 0.0251513 I`.
- Everything else (symbolic `n` with positive order, generic free symbols such as
  `HarmonicNumber[x, 2.5]`, exact non-integer `n` such as `HarmonicNumber[1/17,
  5]`, negative-integer `n`) stays unevaluated.
- **Derivatives** (w.r.t. the continuous index):
  `D[HarmonicNumber[n], n] = Pi^2/6 - HarmonicNumber[n, 2]` and
  `D[HarmonicNumber[n, r], n] = r (Zeta[r+1] - HarmonicNumber[n, r+1])`; the
  `r`-derivative is the generic `Derivative[0,1][HarmonicNumber][n, r]`.

## PolyGamma

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

- `PolyGamma[z]` — the digamma function ψ(z) = Γ′(z)/Γ(z). Rewrites to the
  two-argument form `PolyGamma[0, z]`, except for a machine-precision `z`,
  which is computed directly (`PolyGamma[2.5] === PolyGamma[0, 2.5]`).
- `PolyGamma[n, z]` — the n-th polygamma ψ⁽ⁿ⁾(z) = dⁿ/dzⁿ ψ(z).

Behaviour:

- **Special values.** Non-positive integer arguments are poles:
  `PolyGamma[0]`, `PolyGamma[n, -k]` → `ComplexInfinity`. At infinity,
  `PolyGamma[0, Infinity]` → `Infinity` and `PolyGamma[n, Infinity]` → `0` for
  `n ≥ 1`. `PolyGamma[ComplexInfinity]` and `PolyGamma[n, Indeterminate]`
  → `Indeterminate`.
- **Exact at positive integers.** `PolyGamma[0, m] = H_{m-1} − EulerGamma`
  (a rational minus Euler's constant), e.g. `PolyGamma[5] = 25/12 - EulerGamma`.
  For **odd** order `n ≥ 1` the value closes via ζ(n+1) and is written, as in the
  reference, with `n!` outside the sum: `PolyGamma[3, 2] = 6 (-1 + Pi^4/90)`,
  `PolyGamma[1, 4] = -49/36 + Pi^2/6`, `PolyGamma[3, 1] = Pi^4/15`. For **even**
  order `n ≥ 1` the value involves ζ(odd) and stays symbolic (`PolyGamma[2, 3]`);
  `FunctionExpand` gives `2 (9/8 - Zeta[3])`.
- **Half-integers** (v0.6131, matching the reference). Odd order `n ≥ 1` closes
  through ζ(even): `PolyGamma[1, 1/2] = Pi^2/2`, `PolyGamma[1, 3/2] = -4 + Pi^2/2`,
  `PolyGamma[3, 3/2] = 6 (-16 + Pi^4/6)`. Digamma at any non-integer rational,
  half-integers included, **stays unevaluated** (`PolyGamma[0, 3/2]`), as do even
  orders. `FunctionExpand` gives the closed forms: the Gauss digamma theorem
  (`FunctionExpand[PolyGamma[0, 5/2]] = 8/3 - EulerGamma - Log[4]`,
  `FunctionExpand[PolyGamma[0, 1/3]] = -EulerGamma - Pi/(2 Sqrt[3]) - 3/2 Log[3]`,
  its log terms merged into one log as in the reference) and the ζ(odd) forms
  (`FunctionExpand[PolyGamma[2, 5/2]] = 448/27 - 14 Zeta[3]`). Until v0.6130
  digamma at half-integers evaluated on its own, so a rule written for the
  unevaluated `PolyGamma[0, 3/2]` (FeynCalc's `SimplifyPolyLog`) never fired.
- **Negative order.** `PolyGamma[-1, z] = LogGamma[z]` (the log-gamma function;
  see below). Orders ≤ −2 stay unevaluated.
- **Machine precision** (integer order `0 ≤ n ≤ 40`, real or complex `z`) →
  `sfg_polygamma` / `sfg_polygamma_c` (`sf_machine_gamma.c`), shared with the
  packed / `NDArray` maps and `Compile[]`:
  - ψ, real: a degree-22 Chebyshev interpolant of ψ(x)/(x − x₀) on [1, 2]
    times (x − x₀), x₀ = 1.4616… carried as a double-double, so the value is
    relatively accurate through the zero; integer shifts (exact in double) to
    reach [1, 2] from (0, 10); the asymptotic series for x ≥ 10; reflection
    with an exactly reduced cotangent for x < 0. Median 0.5 / max 3.5 ulp.
  - ψ⁽ⁿ⁾, n ≥ 1, real x > 0: upward recurrence to x ≥ 10 + n/2 and the
    Hurwitz-zeta asymptotic tail (DLMF 25.11.43); every term has one sign,
    each carries a first-order correction for the rounding of x + j, and the
    sum is compensated. `PolyGamma[2, x]` median 0.5 / max 3.3 ulp.
  - ψ, complex: Taylor tiles (as for `Gamma`) on `-5 ≤ Re z < 10`,
    `|Im z| < 19/4`, a Taylor series about x₀ for |z − x₀| ≤ 0.35, the
    asymptotic series after shifting Re z ≥ 10 elsewhere, reflection for
    Re z < -5. ψ⁽ⁿ⁾, n ≥ 1, complex: the recurrence + tail, Re z ≥ 0. Median
    0.7 / max ≈ 4–7 ulp on `0.5 ≤ Re z ≤ 4, |Im z| ≤ 3` (the maxima sit beside
    ψ's real zero, where the error scaled by |z ψ′| is ≤ 2.4 ulp).
  - n ≥ 1 with Re z < 0 declines to the layered path below.
- **Numeric, otherwise.** Inexact real arguments evaluate at arbitrary
  precision (`MPFR`, precision tracked from the input):
  `PolyGamma[100.5] = 4.605174…`, `PolyGamma[1, 3.5] = 0.330358…`,
  `N[PolyGamma[22/10], 50] = 0.5442934367411450377861253708833812285077…`.
  `n = 0` uses MPFR's digamma; higher orders use a recurrence-shift plus the
  Bernoulli asymptotic expansion.
- **Complex.** Inexact complex arguments evaluate (machine or arbitrary
  precision) by the same asymptotic in complex arithmetic, e.g.
  `PolyGamma[2.5 + 3 I] = 1.28127 + 0.979805 I`.
- **Derivatives.** `D[PolyGamma[n, z], z] = PolyGamma[n+1, z]`, e.g.
  `D[PolyGamma[x], x] = PolyGamma[1, x]`.
- Everything else (symbolic `z`, exact non-integer arguments such as
  `PolyGamma[0, 5/2]`, non-integer or complex *order*) stays unevaluated.

*Out of scope:* complex / fractional **order** `n` (so `PolyGamma[6 + I, z]`
stays symbolic), `PolyGamma[n, z]` for `n ≤ −2`, and dedicated `Series`
expansions (only the generic `D`-based Taylor fallback applies).

```mathematica
In[1]:= PolyGamma[5]
Out[1]= 25/12 - EulerGamma

In[2]:= PolyGamma[0, 1/2]
Out[2]= -EulerGamma - Log[4]

In[3]:= PolyGamma[3, 2.5 + I]
Out[3]= 0.0415431 - 0.163057*I
```

- **Non-integer order** (v0.477), in the reference's convention: the
  Riemann–Liouville derivative of order ν, from 0, of ψ(z) + 1/z, with the pole
  term continued in ν,
  ψ^(ν)(z) = z^(-ν)/Γ(1-ν) [-γ + Σ_k (1 - ₂F₁(1,1;1-ν;-z/k))/k] - z^(-1-ν)/Γ(-ν) (Log z - ψ(-ν) - γ)
  (`flint_num_polygamma_frac`, arb; the k-sum runs directly to 2|z| + 2 and its
  tail through Hurwitz zeta values). This is not arb's own `acb_polygamma`
  generalisation, which differs. It agrees with the reference at ν = ±1/2, 3/2,
  z = 1, 2, -3/2, 2 + I to the digits it prints, and tends to the ordinary
  polygamma as ν approaches an integer. Machine and arbitrary precision, real
  and complex z; |z| > 5000 stays unevaluated.

## LogGamma

`LogGamma[z]` is the log-gamma function log Γ(z) — the analytic continuation of
log(Γ(z)) with a single branch cut on the negative real axis (not identical to
`Log[Gamma[z]]`). It stays finite where `Gamma` overflows, so it is the right
primitive for factorial ratios and asymptotics.

**Attributes**: `Listable`, `NumericFunction`, `Protected`.

**Features**:
- **Exact closed forms.** Integers reduce as `LogGamma[n] = Log[(n-1)!]`
  (`LogGamma[5] = Log[24]`); positive half-integers give `Log` of the exact
  `Sqrt[Pi]` form; negative half-integers carry the branch term
  `-Ceiling[-z] Pi I`, e.g. `LogGamma[-3/2] = -2 I Pi + Log[(4 Sqrt[Pi])/3]`.
- **Poles.** Non-positive integers diverge: `LogGamma[0] = LogGamma[-1] = … =
  Infinity`.
- **Symbolic infinities.** `LogGamma[Infinity] = Infinity`,
  `LogGamma[-Infinity] = Indeterminate`, `LogGamma[I Infinity] =
  LogGamma[ComplexInfinity] = ComplexInfinity`.
- **Branch continuity across the reflection.** For `Re[z] < 1/2` the
  arbitrary-precision path (and the machine fallback, `Re z < -30`) use
  `LogGamma[z] = Log[Pi] - Log[Sin[Pi z]] - LogGamma[1-z]`, where the
  log of the sine must be the *continued* one. Factoring
  `Sin[Pi z] = E^(-I Pi z) (E^(2 I Pi z) - 1)/(2 I)` puts the whole winding in
  the exact `-I Pi z` term and leaves a factor near `-1` whose principal log is
  safe. At `Im[z] = 0` this is the limit from above, matching the
  `Im = -Ceiling[-z] Pi` convention used on the real axis. Before 2026-07-27 a
  principal `Log[Pi/Sin[Pi z]]` was used instead: correct only in the strip
  `-1 < Re[z] < 0`, and short by a multiple of `2 Pi I` outside it, so
  `LogGamma[-4.5 + 3. I]` equalled `Log[Gamma[-4.5 + 3. I]]` rather than the
  continuation.
- **Numerics.** Machine real via `lgamma`; arbitrary-precision real via MPFR
  `lgamma`; machine complex via `sfg_loggamma_c` (`sf_machine_gamma.c`, shared
  with the packed / `NDArray` maps and `Compile[]`): Taylor tiles about 30
  centres in the strip `1 ≤ Re < 5/2` reached by the exact recurrence
  (`-5 ≤ Re z < 10`, `|Im z| < 19/4`; on the real-axis row the series is for
  lnΓ/((z−1)(z−2)), so both zeros are carried exactly), shifted Stirling
  elsewhere. The recurrence's N logs are one log of the product, with the
  2πi windings counted as the running product crosses the negative real
  axis — the same continuous branch (Hare 1997). Median 0.5 / max 3.3 ulp on
  `0.5 ≤ Re z ≤ 4, |Im z| ≤ 3` (the Lanczos series it replaces: median 7 /
  max 689 ulp, worst near the zeros at 1 and 2). `Re z < -30` falls back to the
  Lanczos/reflection path described above; arbitrary-precision
  complex via a Stirling series with argument reduction (the continuous
  branch — its imaginary part grows past π where `Log[Gamma]` would wrap).
  Negative real arguments return the complex value `log|Γ(z)| - Pi Ceiling[-z] I`.
- **Derivative.** `D[LogGamma[z], z] = PolyGamma[0, z]`; higher derivatives
  raise the `PolyGamma` order. Produced by `PolyGamma[-1, z]`.

```mathematica
In[1]:= LogGamma[-3/2]
Out[1]= -(2*I) Pi + Log[4/3 Sqrt[Pi]]

In[2]:= LogGamma[-2.5]
Out[2]= -0.0562437 - 9.42478*I

In[3]:= LogGamma[2. + 30. I]
Out[3]= -41.1026 + 74.356*I
```

## Pochhammer


- `Pochhammer[a, n]` — the Pochhammer symbol (rising factorial)
  (a)ₙ = a (a+1) … (a+n-1) = Γ(a+n)/Γ(a).

**Attributes**: `Listable`, `NumericFunction`, `Protected`.

**Features**:
- `Pochhammer[a, 0] = 1` for any `a` (including symbolic and `Infinity`).
- Exact integer order `n` (|n| ≤ 1000) expands to the explicit product of
  `n` linear factors:
  - Symbolic base → a polynomial product, e.g. `Pochhammer[n, 5] =
    n (1 + n) (2 + n) (3 + n) (4 + n)` and `Pochhammer[x, 4] =
    x (1 + x) (2 + x) (3 + x)`.
  - Numeric base → an exact value, e.g. `Pochhammer[10, 6] = 3603600`,
    `Pochhammer[1, 25] = 25!` (GMP BigInt), `Pochhammer[1/2, 3] = 15/8`.
  - Negative `n` gives the reciprocal product, e.g. `Pochhammer[n, -5] =
    1/((-5 + n) (-4 + n) (-3 + n) (-2 + n) (-1 + n))`,
    `Pochhammer[10, -3] = 1/504`.
  - `Pochhammer[0, n] = 0` for positive integer `n` (short-circuited, so
    even `Pochhammer[0, 1285] = 0`); `Pochhammer[Infinity, n] = Infinity`.
- **Machine precision** (real or complex) → `sfg_pochhammer` /
  `sfg_pochhammer_c` (`sf_machine_gamma.c`), shared with the packed /
  `NDArray` maps and `Compile[]`: an integer `n` (|n| ≤ 256) is the product
  `a (a+1) … (a+n−1)` (or `1/((a−1)…(a−|n|))`) in double, any other `n` the
  Gamma ratio with the rounding of `a + n` corrected through ψ, e.g.
  `Pochhammer[2.5, 3] = 39.375`. `Pochhammer[a, 0]` keeps its exact `1`.
- Other numeric arguments evaluate via the Gamma ratio `Γ(a+n)/Γ(a)`,
  reusing the `Gamma` builtin:
  - Exact half-integers reduce to rational multiples of `Sqrt[Pi]`, e.g.
    `Pochhammer[3/2, 1/2] = 2/Sqrt[Pi]`, `Pochhammer[1/2, 1/2] = 1/Sqrt[Pi]`.
  - Machine-precision real → e.g. `Pochhammer[2.4, 8.5] = 2.31022×10⁶`.
  - Arbitrary precision (MPFR) tracks the input precision, e.g.
    `N[Pochhammer[1/3, 7], 50]` and
    `Pochhammer[1.011111111111000000000000000, 8] = 41552.275849087780380888…`.
  - Machine-precision complex → e.g.
    `Pochhammer[2. + 5 I, 8 I] = 2.13868×10⁻⁶ − 1.42187×10⁻⁵ I`.
- Threads over lists (Listable), e.g.
  `Pochhammer[{2, 3, 5, 7, 11}, 3] = {24, 60, 210, 504, 1716}`.
- All other arguments (e.g. `Pochhammer[a, n]`, `Pochhammer[a, 1/2]`,
  `Pochhammer[1/2, 1/3]`) stay unevaluated. (Derivatives and series, which
  the reference expresses through `PolyGamma`, are not yet implemented.)

## HypergeometricPFQ


`HypergeometricPFQ[{a1, …, ap}, {b1, …, bq}, z]` is the generalized
hypergeometric function ₚFₚ(a; b; z), defined by the series

```
HypergeometricPFQ[{a1,…,ap}, {b1,…,bq}, z]
  = Sum_{k>=0} (Pochhammer[a1,k] … Pochhammer[ap,k])
             / (Pochhammer[b1,k] … Pochhammer[bq,k]) * z^k / k!
```

**Features**:
- Attributes `NumericFunction`, `Protected`.
- `HypergeometricPFQ[a, b, 0]` is `1`; threads over a `List` third argument.
- Parameters common to the upper and lower lists cancel; a non-positive
  integer upper parameter terminates the series to a polynomial (valid for
  symbolic `z`).
- Reduces to elementary functions for simple parameters: `0F0 -> E^z`,
  `1F0(a) -> (1-z)^(-a)`, `0F1(1/2) -> Cosh[2 Sqrt[z]]`,
  `0F1(3/2) -> Sinh[2 Sqrt[z]]/(2 Sqrt[z])`, `1F1(1;2) -> (E^z-1)/z`,
  `2F1(1,1;2) -> -Log[1-z]/z`. The central-binomial / arcsin family also closes:
  `2F1(1,1;1/2;z) -> 1/(1-z) + Sqrt[z] ArcSin[Sqrt[z]]/(1-z)^(3/2)`,
  `2F1(1,1;3/2;z) -> ArcSin[Sqrt[z]]/(Sqrt[z] Sqrt[1-z])`,
  `2F1(2,1;3/2;z) -> 1/(2(1-z)) + ArcSin[Sqrt[z]]/(2 Sqrt[z] (1-z)^(3/2))`.
  The error-function family (as the reference evaluates it):
  `1F1(1/2;3/2;-z) -> Sqrt[Pi] Erf[Sqrt[z]]/(2 Sqrt[z])` (`Erfi` for `+z`), and
  every `1F1(a;a+1;z)` with `a` a half-integer or a positive integer through
  `F_(a+1) = (a+1)(E^z - F_a)/z` (`1F1(-1/2;1/2;z) -> E^z - Sqrt[Pi] Sqrt[z]
  Erfi[Sqrt[z]]`, `1F1(2;3;z)` elementary); Kummer's `1F1(1;c;z) = E^z
  1F1(c-1;c;-z)` onto that family. Other `a` stay — the incomplete-gamma form
  `a (-z)^-a Gamma[a, 0, -z]` would have no numerics for real `z > 0`.
  `2F1(1/2,1;3/2;z) -> ArcTanh[Sqrt[z]]/Sqrt[z]` (`ArcTan[Sqrt[-z]]/Sqrt[-z]`
  for a negated `z`), `2F1(1/2,1/2;3/2;z) -> ArcSin[Sqrt[z]]/Sqrt[z]`
  (`ArcSinh` for a negated `z`). Gauss: `2F1(a,b;c;1) = Γ(c)Γ(c−a−b)/(Γ(c−a)Γ(c−b))`
  for rational parameters with `c − a − b > 0` (symbolic ones stay, as in the
  reference). Each identity is verified to 30 digits against the series in
  `tests/test_hypergeopfq.c`.
  A table-backed very-well-poised reduction closes the classic
  central-binomial-cubed Ramanujan `1/Pi` series:
  `4F3({1/2,1/2,1/2,5/4}, {1/4,1,1}, -1) -> 2/Pi` (general `1/Pi` summation is an
  open problem, so this is a recognizer for the known class, not universal).
- Numeric evaluation at machine, arbitrary (MPFR), and complex precision, with
  output precision tracking the input. The general series is summed where it
  converges — `p <= q` (entire) and `p == q+1` with `|z| < 1`; for `p > q+1`,
  and `p == q+1` with `|z| >= 1` other than `2F1`, the call stays unevaluated.
  `2F1` (`{a, b}, {c}`) is analytically continued to the whole plane (see
  `Hypergeometric2F1`) — at arbitrary precision too (arb's `acb_hypgeom_2f1`)
  everywhere off the cut `[1, ∞)`: `N[HypergeometricPFQ[{2/3,11/12},{5/3},-1], 20]`
  is a number.
- **Machine precision: one raw-`double` kernel** (`sf_machine_hyper.c`), real and
  complex (raw re/im), answering every machine-precision call — exact
  (Integer / Rational) parameters included — *before* the exact rewrites, so the
  scalar result is bit-identical to the packed / NDArray buffer and `Compile[]`
  results. `{}/{b}`, `{a}/{b}` and `{a, b}/{c}` dispatch to the
  `Hypergeometric0F1` / `1F1` / `2F1` kernels, so every spelling of the same
  function returns the same bits.
- **Cancellation and rounding are measured, not ignored.** Each series is
  summed once in double with compensated (TwoSum) summation while tracking its
  largest term and where it occurs. If the largest term exceeds the sum by more
  than 2× (cancellation — e.g. the alternating series at negative `z`, whose
  terms reach `e^|z|` times the result) or lies deeper than the 6th term (each
  term then carries the rounding of the ratio products before it), the same
  series is re-summed in double-double arithmetic (~106 bits); beyond 2^45 of
  cancellation the kernel declines and the MPFR series answers.
  `1F1(1;2;-40)` agrees with `(E^z-1)/z` to one ulp; `2F1(1,2;3;0.89)` (380
  terms) is within 1.4 ulp where the plain sum was 13 ulp off.
- Compiles (`Compile`, `Plot`, `NIntegrate`, `NDSolve`, …) at real and complex
  scalar shapes; a complex argument lowers to the kernel's complex arm
  (`OP_KERNN_C`). The rank-1 array shape of the n-ary heads is still open
  (`COMPILE_MISSING.md`).
- `D[HypergeometricPFQ[{a},{b},x], x]
   = (prod a_i / prod b_j) HypergeometricPFQ[{a_i+1},{b_j+1},x]`.

```mathematica
In[1]:= HypergeometricPFQ[{a1, a2, a3}, {b1, b2, b3}, 0]
Out[1]= 1

In[2]:= HypergeometricPFQ[{}, {}, z]
Out[2]= E^z

In[3]:= HypergeometricPFQ[{a, b, c}, {a, d, e}, z]
Out[3]= HypergeometricPFQ[{b, c}, {d, e}, z]

In[4]:= HypergeometricPFQ[{1, 1}, {3, 3, 3}, 2.]
Out[4]= 1.07893

In[5]:= HypergeometricPFQ[{1, 2, 3, 4}, {5, 6, 7}, {0.1, 0.3, 0.5}]
Out[5]= {1.01164, 1.03627, 1.06296}

In[6]:= N[HypergeometricPFQ[{1, 1, 1}, {3/2, 3/2, 3/2}, 10], 50]
Out[6]= 530.19188827362590438855961685444087792733053398358

In[7]:= D[HypergeometricPFQ[{a1, a2}, {b1, b2, b3}, x], x]
Out[7]= (a1 a2 HypergeometricPFQ[{1 + a1, 1 + a2}, {1 + b1, 1 + b2, 1 + b3}, x])/(b1 b2 b3)
```

> **Packed arrays.** With a real or complex machine buffer as `z` and every
> parameter a machine Real or exact Integer / Rational, the machine kernel is
> mapped over the buffer — the same kernel the scalar call uses, so the results
> agree bit-for-bit. For real parameters over a float64 buffer a vector driver
> hoists the parameter work out of the loop (for `2F1`, the series coefficients
> are tabulated once per buffer; the scalar path forms the identical
> coefficients on the fly). A `2F1` element on the cut widens the result to
> complex; an int64 buffer, an MPFR parameter, or a declined element takes the
> element-wise route. For `2F1` series that are one-signed and never grow
> (decided from the parameters), the buffer sums four elements in lock-step with
> the scalar loop's exact arithmetic. 10^6 elements: `1F1(1;2;x)` 11 ms (was
> 1.3 s; SciPy 27 ms), `2F1(1,2;3;x)` 0.14 s (was 2.6 s; SciPy 0.16 s),
> `2F2(1,2;3,4;x)` 0.19 s (was 1.3 s).

- **Integer-shifted parameter pairs** (v0.6490). An upper parameter that exceeds a
  lower one by a positive integer n (n ≤ 12) reduces through
  `pFq[U ∪ {b+n}, L ∪ {b}, z] = P(θ) pFq[U, L, z]`, θ = z d/dz,
  P(k) = (b+n)ₖ/(b)ₖ, when the reduced series closes:
  `HypergeometricPFQ[{2, 2}, {1, 1}, z] -> E^z (1 + 3 z + z^2)` (the Touchard /
  Dobinski sums), `HypergeometricPFQ[{a, c+1}, {c}, z] -> (1-z)^(-a-1) (1 - z + a z/c)`.
  A lower parameter that is a non-positive integer is a pole and is never paired.

**Classical summation theorems (v0.6491).** At z = 1: Gauss 2F1, Pfaff–Saalschütz
(a terminating, balanced 3F2, including a symbolic degree), Dixon 3F2, and
Dougall's very-well-poised 5F4. At z = −1: Kummer 2F1. Parameters are matched in
any order by exact comparison of the linear relations. A numeric parameter set
outside the region of convergence (for example Re(c − a − b) ≤ 0 for Gauss) is
declined. An upper parameter b + n over a lower b (n = 1…12) reduces through
P(z d/dz) when the reduced series closes.

## Hypergeometric0F1


`Hypergeometric0F1[b, z]` is the confluent hypergeometric function ₀F₁,
equal to `HypergeometricPFQ[{}, {b}, z]`.

- Attributes `Listable`, `NumericFunction`, `Protected`.
- Machine kernel (real and complex `b`, `z`): the series where it loses at most
  a bit; for `z < 0` the double-double re-sum up to `|z| ~ 250` and beyond that
  Miller's backward recurrence for `J_{b-1}(2 Sqrt[-z])` written in `-z` itself
  (no rounded `2 Sqrt[-z]`), normalised by the Neumann series (Gautschi 1967) —
  accurate to the envelope of the oscillation; for `z > 900`, the Hankel
  expansion of `I_{b-1}` with `2 Sqrt[z]` formed in double-double. Scalar 0.9 µs.

```mathematica
In[1]:= Hypergeometric0F1[1/2, z]
Out[1]= Cosh[2 Sqrt[z]]
```

## Hypergeometric1F1


`Hypergeometric1F1[a, b, z]` is Kummer's confluent hypergeometric function
₁F₁, equal to `HypergeometricPFQ[{a}, {b}, z]`.

- Attributes `Listable`, `NumericFunction`, `Protected`.
- Machine kernel (real and complex): closed forms for `a = 0`, `a = b` and
  `(a, b) = (1, 2)` (`expm1(z)/z`, complex form free of cancellation at small
  `|z|`); a terminating polynomial for non-positive integer `a` (also past a
  lower pole it never reaches); for `|z| > 20` the DLMF 13.7.2 expansion with
  both exponential sectors (the subdominant one with the Stokes-line `cos`
  factor, exact in the terminating cases), accepted only when its least term is
  below half an ulp; for `Re z < 0` Kummer's transformation
  `e^z 1F1(b-a; b; -z)` (one-signed terms, `b - a` carried in double-double);
  otherwise the series. Scalar 0.15 µs (real), 0.35 µs (complex).

```mathematica
In[1]:= Hypergeometric1F1[a, b, 0]
Out[1]= 1

In[2]:= Hypergeometric1F1[-2, 1, z]
Out[2]= 1 - 2 z + 1/2 z^2

In[3]:= Hypergeometric1F1[1., 2., -40.]
Out[3]= 0.025

In[4]:= Hypergeometric1F1[1 + I, 2, 0.5]
Out[4]= 1.27058 + 0.322122*I
```

## Hypergeometric2F1


`Hypergeometric2F1[a, b, c, z]` is the Gauss hypergeometric function ₂F₁,
equal to `HypergeometricPFQ[{a, b}, {c}, z]`.

- Attributes `Listable`, `NumericFunction`, `Protected`.
- Machine kernel, real `z`: direct series on `[0, 0.9]`, Pfaff's
  `(1-z)^-a F(a, c-b; c; z/(z-1))` on `[-9, 0)`, Gauss's sum at `z = 1`, and the
  connection formulas beyond. Complex `z` (real parameters): the cheapest of the
  direct / Pfaff series (|w| ≤ 0.75), the `1-z` and `1/(1-z)` connection
  formulas (Abramowitz & Stegun 15.3.6, with the degenerate logarithmic cases
  15.3.10–15.3.12 when `c-a-b` or `b-a` is an integer — the digamma differences
  carried by their own recurrences), accepted only while the combination's
  condition number is ≤ 2; otherwise — around `z = e^{±iπ/3}`, where every
  transformation has `|w| ≈ 1`, or for ill-conditioned parameters — analytic
  continuation of the hypergeometric ODE by successive Taylor steps (Bühring
  1987) from `z0 = 1/2 ± 0.45 i`. Complex parameters: series for `|z| ≤ 0.9`.
- **On the cut** `z > 1` a machine-real `z` returns the complex limit from below
  (the reference / mpmath convention): `Hypergeometric2F1[1, 2, 3, 2.5]` is
  `-0.929749 - 1.00531 I` (it used to stay unevaluated, as did every real
  `z < -1` and every complex `|z| ≥ 1`).
- Scalar 0.3 µs (real), 1.7 µs (complex); 10^6 packed elements 0.14 s.
- **Arbitrary precision** (`|z| < 1`) goes through arb's transformation-based
  `acb_hypgeom_2f1` rather than the direct series, whose term count grows like
  `1/(1-|z|)`: `N[Hypergeometric2F1[1/3, 1/2, 1/5, 99/100], 50]` dropped from
  ~60 ms to under 1 ms.

```mathematica
In[1]:= Hypergeometric2F1[1, 1, 2, z]
Out[1]= -Log[1 - z]/z

In[2]:= Hypergeometric2F1[1, 1, 3/2, z]
Out[2]= ArcSin[Sqrt[z]]/(Sqrt[z] Sqrt[1 - z])

In[3]:= Hypergeometric2F1[1, 2, 3, 2.5]
Out[3]= -0.929749 - 1.00531*I

In[4]:= N[Hypergeometric2F1[1/3, 1/2, 1/5, 99/100], 50]
Out[4]= 24.8565393133858972090022500136712230259988212774851
```

## HypergeometricU

`HypergeometricU[a, b, z]` is the Tricomi confluent hypergeometric function
`U(a, b, z)`, the second (recessive) solution of the confluent equation
`z w'' + (b - z) w' - a w = 0`. It has a branch cut in `z` running from `-∞` to
`0` (from the `z^(1-b)` term in its connection to `Hypergeometric1F1`).

- **Exact reductions.** `HypergeometricU[0, b, z] = 1` and
  `HypergeometricU[a, a+1, z] = z^(-a)` (so `HypergeometricU[1, 2, z] = 1/z`,
  `HypergeometricU[2, 3, 4] = 1/16`). A non-positive integer `a = -m`
  terminates (DLMF 13.2.7) and returns the expanded polynomial
  `(-1)^m Σ_{s=0}^m Binomial[m, s] (b+s)_(m-s) (-z)^s` for any `b`, `z`:
  `HypergeometricU[-2, b, z] = b + b^2 - 2 z - 2 b z + z^2`,
  `HypergeometricU[-2, b, 2.5] = 1.25 - 4.0 b + b^2`. (An all-machine call
  still takes the double kernel first.)
- **Numerics.** With at least one inexact argument the value is computed at
  machine or arbitrary (MPFR) precision, real or complex:
  `HypergeometricU[1.5, 7, 4.5] = 0.554931`,
  `N[HypergeometricU[3/2, 7, 9/2], 25]` to 25 digits,
  `HypergeometricU[3 - I, 2 I, 3.0]` complex. A machine-precision call (real
  `a`, `b` exact or machine; `z` machine real or complex) uses a raw-`double`
  kernel, checked before the exact reductions: the asymptotic series for large
  `|z|` (compensated, accepted only when its least term is below half an ulp —
  it used to accept 1e-13), the DLMF 13.2.9 confluent-limit log series for
  integer `b`, the Kummer two-`₁F₁` connection for non-integer `b` (each
  accepted only when its combination amplifies rounding by at most 4× — these
  guards used to allow 1e13 and 1e6, and e.g. `HypergeometricU[5.5, 0.3, 3.]`
  came back 1.5e-10 off), and between them, where both series cancel, the
  integral representation by double-exponential quadrature (exp-sinh rule, 50–90
  nodes, ~4 µs; FLINT took 280 µs there). Real `z < 0` returns the principal
  complex value; a complex machine `z` uses the complex twin of the same
  regimes. A terminating case — `a` or `a-b+1` a non-positive integer `-m` —
  is a finite sum, `U = z^(-a) Σ_{k≤m} (a)_k (a-b+1)_k/k! (-1/z)^k` (for
  `a = -m` the polynomial `(-1)^m (b)_m ₁F₁(-m; b; z)`), summed whole rather than
  cut at a growing term, and declined to FLINT only on measured cancellation:
  `HypergeometricU[-1, 2, 1.] = -1.0` (`z - 2`), `HypergeometricU[-2, 1, 1.37] =
  -1.6031` (`z² - 4z + 2`); these previously returned `0.0`. Complex parameters
  or arbitrary precision use FLINT's rigorous `acb_hypgeom_u`.
- **Derivative.** `D[HypergeometricU[a, b, z], z] = -a HypergeometricU[a+1, b+1, z]`
  (DLMF 13.3.22); `d/da` and `d/db` stay as inert `Derivative[...]`.
- **Performance (machine reals).** A real machine argument with `z > 0` shares a
  hardware-`double` kernel bit-for-bit with the NDArray / packed fast path and
  `Compile[]`, at the scalar and at the rank-1 array shape
  (`Compile[{{v, _Real, 1}}, HypergeometricU[a, b, v]]`; the n-ary kernel threads
  over its last argument, as for `WhittakerM`/`WhittakerW`, `Hypergeometric1F1`,
  `Hypergeometric2F1`, `LerchPhi`).
- **Listable.** `HypergeometricU[1, 2, {2., 4., 5.}] = {0.5, 0.25, 0.2}`.

The general (two-branch) Series at the origin and `FunctionExpand` are not
produced (`U` is not analytic at `0`), matching the confluent siblings.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= HypergeometricU[1, 2, z]
Out[1]= 1/z

In[2]:= HypergeometricU[1.5, 7, 4.5]
Out[2]= 0.554931

In[3]:= HypergeometricU[1/2, 1, 2.]
Out[3]= 0.645694

In[4]:= D[HypergeometricU[a, b, z], z]
Out[4]= -a HypergeometricU[1 + a, 1 + b, z]

In[5]:= {HypergeometricU[-1, 2, 1.], HypergeometricU[-2, 1, 1.37]}
Out[5]= {-1.0, -1.6031}
```

## Zeta

- `Zeta[s]` — the Riemann zeta function ζ(s) = Σ_{k≥1} k⁻ˢ (Re s > 1; analytic
  continuation elsewhere).
- `Zeta[s, a]` — the Hurwitz (generalized) zeta function ζ(s, a) = Σ_{k≥0} (k+a)⁻ˢ,
  summed as ((k+a)²)^(−s/2) with the term k + a = 0 dropped.
- `Zeta[s, a, IncludeSingularTerm -> True]` — keeps that term.

Attributes: `Listable`, `NumericFunction`, `Protected`. Options:
`IncludeSingularTerm -> False`.

- **Exact integer arguments** (via Bernoulli numbers):
  - Even positive integers are rational multiples of πⁿ:
    `Zeta[2] = Pi^2/6`, `Zeta[4] = Pi^4/90`, `Zeta[6] = Pi^6/945`.
  - Negative integers are rational: `Zeta[-1] = -1/12`, `Zeta[-3] = 1/120`;
    negative even integers are the trivial zeros, `Zeta[-2] = 0`.
  - `Zeta[0] = -1/2`, and the pole `Zeta[1] = ComplexInfinity`.
  - Odd positive integers have no closed form and stay symbolic (`Zeta[3]`,
    `Zeta[5]`). `Zeta[Infinity] = 1`.
- **Exact Hurwitz at a positive integer `a`.** `Zeta[s, m]` reduces to
  `Zeta[s] - Σ_{k=1}^{m-1} k⁻ˢ`, valid for symbolic / exact `s`, e.g.
  `Zeta[3, 2] = -1 + Zeta[3]` and `Zeta[4, 5] = Pi^4/90 - 22369/20736`.
  `Zeta[s, 1] = Zeta[s]` and `Zeta[s, 1/2] = (2^s - 1) Zeta[s]`
  (e.g. `Zeta[2, 1/2] = Pi^2/2`, `Zeta[3, 1/2] = 7 Zeta[3]`). For a symbolic or
  non-integer `s` only `a = 2` expands (its one extra term is `1^-s`):
  `Zeta[x, 2] = -1 + Zeta[x]`, while `Zeta[x, 3]` and `Zeta[3/2, 3]` stay; past
  `a = 2`, an integer `s` must also have a closed form `Zeta[s]`, so `Zeta[3, 4]`
  stays too.
- **Exact Hurwitz at `a ≤ 0`** (v0.476). With m = ⌊−a⌋ + 1 shifts,
  `Zeta[s, a] = Σ_{k<m, k+a≠0} |k+a|⁻ˢ + Zeta[s, a+m]`; the dropped term is 0⁻ˢ,
  except at `s = 0`, where 0⁰ = 1 counts. Expanded on the same rules as `a > 0`:
  `a = 0` and `a = −1` for any `s` (`Zeta[x, 0] = Zeta[x]`,
  `Zeta[x, -1] = 1 + Zeta[x]`, `Zeta[2 + I, -1] = 1 + Zeta[2 + I]`); a more
  negative integer `a` when `Zeta[s]` has a closed form
  (`Zeta[2, -10] = 1968329/1270080 + Pi^2/6`, `Zeta[-1, -2] = 35/12`,
  `Zeta[0, -1] = 3/2`, `Zeta[3, -2]` stays); a non-integer rational `a`, either
  sign, when `s` is an exact integer and the shifted `Zeta[s, r]`, r ∈ (0, 1),
  evaluates (`Zeta[2, -3/2] = 40/9 + Pi^2/2`, `Zeta[3, -1/2] = 8 + 7 Zeta[3]`,
  `Zeta[3, 5/2] = -224/27 + 7 Zeta[3]`, `Zeta[-1, -1/2] = 13/24`; `Zeta[2, -2/3]`
  stays). The finite part is one exact GMP rational, so `Zeta[2, -200]` costs a
  sum and a two-term `Plus`. `Zeta[0, a] = 1/2 - a` for a rational `a > 0`.
  Shifts stop at 10⁵ terms (`Zeta[2, -100001]` stays; the reference expands it).
- **`IncludeSingularTerm -> True`.** At an integer `a ≤ 0` (exact or inexact)
  the kept term is 0⁻ˢ: `ComplexInfinity` for `Re s > 0`
  (`Zeta[2, -1, IncludeSingularTerm -> True]`, `Zeta[2., -1., …]`); no change at
  `s = 0` or a negative integer `s` (`Zeta[-1, -1, …] = 11/12`); an exact
  non-integer `s < 0` or a symbolic `s` stays unevaluated (with the option). At a
  non-integer `a` the option changes nothing. `Rule` and `RuleDelayed`, symbol
  and string keys are accepted, and the first setting wins. An unknown option is
  `Zeta::optx`, a non-Boolean value `Zeta::opttf`, a third positional argument
  `Zeta::nonopt` (naming the last one); no arguments is `Zeta::argt`.
  `Options[Zeta]`/`SetOptions` hold the default.
- **Machine numerics.** A machine-precision call (real or complex `s`; for
  `Zeta[s, a]` with `Re a > 0` or a real `a ≤ 0`) is answered by the raw-`double` kernels of
  `src/special_functions/sf_machine_zeta.c` — the same helpers the packed /
  `NDArray` kernels and `Compile[]` use, so all surfaces agree bit for bit:
  - integer `s`: tabulated ζ(2…64) and ζ(0…−60) (`Zeta[2.]` is π²/6 to the
    last bit);
  - `s ≥ 20`: the Dirichlet series summed directly (eight terms);
  - `Re s ≥ 1/2`: Euler–Maclaurin (DLMF 25.2.9) with the direct sum run to
    w ≥ max(10, |s|/2 + 4) and the Bernoulli correction cut at 2⁻⁵⁶;
  - real `s < 1/2`: the functional equation with sin(πs/2) reduced exactly,
    (2π)ˢ from an exact two-product, and the rounding of 1−s corrected to
    first order in Γ and ζ.
  Complex `s` with `Re s < 1/2`, the pole and anything the kernels decline
  (cancellation, overflow) fall through to the MPFR path below. Accuracy
  against mpmath (2000 points): real s ∈ [1.5, 9] median 0.3 / max 1.0 ulp;
  complex s ∈ [1.5, 5] × [−3, 3] median 0.3 / max 1.5 ulp; real s ∈ [−40, ½]
  max ≈ 4.6 ulp. Net latency ≈ 0.35 µs real, 0.85 µs complex (from 75 / 403 µs).
  - real `a ≤ 0` (v0.476, `sf_zeta2_r/_c`): the head Σ_{k<m} |k+a|⁻ˢ
    (TwoSum-formed |k+a|, first-order power correction, Neumaier sum, up to
    m = 10⁷) plus the Hurwitz kernel at a + m ∈ (0, 1]. Within 1.6 ulp of a
    40-digit reference on the probe set (`Zeta[2., -100000.5]` 0.6 ulp);
    `Zeta[2.5, -1.5 - k/10⁴]` costs 1.1 µs (0.34 ms before, through MPFR; the
    reference 13 µs). `HurwitzZeta` keeps the principal powers and calls the
    unshifted kernels. `Zeta[s, a]` with an array in either slot has a binary
    `NDArray` kernel and compiles at scalar and rank-1 shapes.
  - exact even `s`: B₂ₙ from FLINT's `arith_bernoulli_number` (cold `Zeta[400]`
    0.5 ms, from 31 ms through the O(n²) recurrence), and the finished
    `rational · Pi^(2n)` memoised per n.
- **Arbitrary precision.** One-argument zeta routes through arb's
  `acb_dirichlet_zeta` (the `flint_num_zeta` bridge), ~4× faster than `mpfr_zeta`
  and faster than the reference at 50 digits. Two-argument (Hurwitz) inputs, and
  machine inputs the kernels decline, use the Euler–Maclaurin complex-MPFR kernel
  (Riemann is the `a = 1` case). Precision tracks the input, e.g.
  `N[Zeta[3], 50] = 1.2020569031595942853997381615114499907649862923405…`,
  `N[Zeta[5/4], 50] = 4.5951118258429433806853780396946256522810297806048…`,
  `N[Zeta[1/2 + I/2]] = -0.459303 - 0.961254 I`,
  `Zeta[-1.5 + I, 2.5 - I] = 0.0184868 + 1.67553 I`. The two-argument kernel uses
  the symmetric power `((a+k)^2)^(-s/2)` (= `|a+k|^-s` on the real axis), so for
  `Re(a) < 0` it differs from the principal-branch `HurwitzZeta`:
  `Zeta[3, -1.5] = 16.7107` while `HurwitzZeta[3, -1.5] = 0.118102`.
- **Derivatives.** `D[Zeta[s, a], a] = -s Zeta[1+s, a]` (the `s`-derivative is the
  generic `Derivative[1,0][Zeta][s, a]`); `D[Zeta[s], s] = Derivative[1][Zeta][s]`.
  That inert derivative, and `Derivative[n][Zeta]` for `1 <= n <= 64`, has a
  value at an inexact real or complex point (arb's Taylor jet of ζ, a system
  rule through the derivative registry, not a `SubValue`): `Zeta'[2.]` is
  `-0.937548`, a machine argument answers a machine number and an MPFR one a
  number at its precision. An exact point and the pole `s = 1` stay inert.
  `-Zeta'[s]/Zeta[s]` (the Dirichlet series of `MangoldtLambda`) and a Newton
  step on `Zeta` are therefore numbers.
- **Series.** `Series[Zeta[x], {x, 1, n}]` gives the Laurent expansion that
  defines the Stieltjes constants,
  `1/(x-1) + EulerGamma - StieltjesGamma[1](x-1) + ½ StieltjesGamma[2](x-1)² + …`,
  and `Series[Zeta[x], {x, 0, n}]` (n ≤ 3) the Taylor expansion
  `-1/2 - ½ Log[2 Pi] x + …`. Other expansion points fall back to the generic
  differentiation path.

```mathematica
In[1]:= Zeta[2]
Out[1]= 1/6 Pi^2

In[2]:= Series[Zeta[x], {x, 1, 2}] // Normal
Out[2]= EulerGamma + 1/(-1 + x) - StieltjesGamma[1] (-1 + x) + 1/2 StieltjesGamma[2] (-1 + x)^2

In[3]:= N[Zeta[2 + 3 I], 30]
Out[3]= 0.7980219851462757206222945007251 - 0.1137443080529385002159133658573*I

In[4]:= Zeta[-2.5, -0.5]
Out[4]= 0.169765
```

## HurwitzZeta

- `HurwitzZeta[s, a]` — the Hurwitz zeta function ζ(s, a) = Σ_{k≥0} (k+a)⁻ˢ
  (Re s > 1; analytic continuation elsewhere).

Attributes: `Listable`, `NumericFunction`, `Protected`.

`HurwitzZeta` is identical to `Zeta[s, a]` for Re a > 0, but it sums the
**principal-branch** power (k+a)⁻ˢ rather than ((k+a)²)^(−s/2). Two consequences
follow: the functions disagree for non-positive real `a` (e.g.
`HurwitzZeta[3, -3.5] = 0.0307784` versus `Zeta[3, -3.5] = 16.798`), and
`HurwitzZeta` keeps the singular summands `Zeta` discards, so it has poles at
`a = 0, -1, -2, …`.

- **Exact reductions** (when neither argument is inexact):
  - `HurwitzZeta[s, 1] = Zeta[s]`, inheriting all the Riemann closed forms;
    `Table[HurwitzZeta[s, 1], {s, -2, 2}] = {0, -1/12, -1/2, ComplexInfinity, Pi^2/6}`.
  - `HurwitzZeta[s, 1/2] = (2^s - 1) Zeta[s]`.
  - A positive integer `a = m` reduces to `Zeta[s] − Σ_{k=1}^{m-1} k⁻ˢ`, e.g.
    `HurwitzZeta[3, 2] = -1 + Zeta[3]` and `HurwitzZeta[2, 5] = -205/144 + Pi^2/6`.
- **Poles and special points.** `HurwitzZeta[1, a] = ComplexInfinity` for any `a`.
  At a non-positive integer `a`: a positive integer `s` gives `ComplexInfinity`
  (`HurwitzZeta[2, -2] = ComplexInfinity`), while a non-positive integer `s` gives
  the finite Bernoulli-polynomial value `−BernoulliB[1-s, a]/(1-s)`, e.g.
  `HurwitzZeta[0, 0] = 1/2` and `HurwitzZeta[-1, 0] = -1/12`.
- **Machine numerics.** Machine-precision `s` and `a` with `s ≠ 1`,
  `Re s ≥ 1/2` and `Re a > 0` (real or complex) use the raw-`double`
  Euler–Maclaurin kernels `sf_hurwitz_r` / `sf_hurwitz_c`
  (`sf_machine_zeta.c`), shared bit for bit with the packed / `NDArray`
  kernels and `Compile[]`. Each k + a is formed by an exact TwoSum and its
  power corrected to first order; a complex sum whose terms exceed 16× its
  value (large `Im s` with `Im a ≠ 0`) declines rather than return the noise.
  `a = 1` is the Riemann kernel. Accuracy against mpmath (2000 points): real
  s ∈ [1.5, 6], a ∈ [0.5, 4] median 0.4 / max 1.6 ulp; complex
  s ∈ [1.5, 4] × [−2, 2], a ∈ [0.5, 3] × [−1, 1] median 0.8 ulp, max ≈ 8 ulp
  where the sum cancels — within the function's own condition number. Net
  latency ≈ 0.27 µs real, 1.0 µs complex (from 275 / 709 µs).
- **Arbitrary precision and declines.** Everything else evaluates through the
  Euler–Maclaurin complex-MPFR kernel (arb's `acb_dirichlet_hurwitz` above
  machine precision); a numeric pole gives `ComplexInfinity`. Precision tracks
  the input, e.g.
  `HurwitzZeta[3, .2] = 125.739`, `HurwitzZeta[7., 5] = 1.84949×10⁻⁵`,
  `HurwitzZeta[2.3, 8 + I] = 0.0544701 - 0.00944852 I`,
  `N[HurwitzZeta[1/3, 8/7], 50] = -1.1389367444490991746548674334535727810961919460755…`.
- **Re s < 0.** The Euler–Maclaurin head terms grow like (N+|a|)^(−Re s) while
  the value is far smaller, so the kernel adds ~(−Re s)·log₂(N+|a|) guard bits;
  at machine precision such a call goes to arb (rigorous, ~5× faster for large
  −s) and is lowered to machine Reals, or kept as an extended-exponent value
  outside the double range: `HurwitzZeta[-38.5, 0.5] = -1.36481×10¹⁴`,
  `HurwitzZeta[-300.5, 0.5] = 1.675678219245527×10³⁷⁵` (formerly `1.2×10⁴²`
  and `Overflow[]`). For non-integer `s` and `a < 0` the principal-branch
  summands with k + a < 0 are complex, so the value is too:
  `HurwitzZeta[-2.5, -0.5] = -0.00701133 + 0.176777 I` (the imaginary part used
  to be dropped). The same guard applies to `Zeta[s, a]`
  (`Zeta[-38.5, 2.3] = -1.71976×10¹⁴`, formerly `-1.45×10²²`).
- **Derivatives.** `D[HurwitzZeta[s, a], a] = -s HurwitzZeta[1+s, a]` (the
  `s`-derivative is the generic `Derivative[1,0][HurwitzZeta][s, a]`); higher
  `a`-derivatives follow the rising-factorial pattern
  `(-1)^k Pochhammer[s, k] HurwitzZeta[k+s, a]`, including a **symbolic order**
  `D[HurwitzZeta[s, a], {a, n}] = (-1)^n Pochhammer[s, n] HurwitzZeta[n+s, a]`.

```mathematica
In[1]:= HurwitzZeta[s, 1/2]
Out[1]= (-1 + 2^s) Zeta[s]

In[2]:= HurwitzZeta[3, -3.5]
Out[2]= 0.0307784

In[3]:= HurwitzZeta[-3, -2]
Out[3]= -1079/120
```

## StieltjesGamma

`StieltjesGamma[n]` is the n-th Stieltjes constant γₙ, the coefficients of the
Laurent expansion of `Zeta` about s = 1:
ζ(s) = 1/(s-1) + Σ_{n≥0} ((-1)ⁿ/n!) γₙ (s-1)ⁿ. `StieltjesGamma[n, a]` is the
generalized constant γₙ(a), the same coefficients for ζ(s, a).

Attributes: `Listable`, `NHoldFirst`, `Protected` (`NHoldFirst` since v0.477).
An exact call is inert apart from `StieltjesGamma[0] = EulerGamma` and
`StieltjesGamma[n, 1] = StieltjesGamma[n]` (they have no elementary closed form);
it appears in the `Series` expansions of `Zeta`.

- **Numeric** (v0.477): an inexact `a` evaluates through arb's
  `acb_dirichlet_stieltjes`, at `a`'s precision (`StieltjesGamma[1, 2.] =
  -0.0728158…`, complex for `a < 0`), and `N[StieltjesGamma[3]]`, `N[…, 30]`
  evaluate (`numeric_stieltjes_n`; the order stays exact under `N`). `a = 0` or a
  negative integer is `ComplexInfinity`.
- The order must be a non-negative machine integer: `StieltjesGamma[-1]`,
  `StieltjesGamma[1.]` give `StieltjesGamma::intnm` and stay unevaluated.

```mathematica
In[1]:= StieltjesGamma[0]
Out[1]= EulerGamma

In[2]:= N[StieltjesGamma[3]]
Out[2]= 0.00205383
```

## BernoulliB

- `BernoulliB[n]` — the Bernoulli number Bₙ.
- `BernoulliB[n, x]` — the Bernoulli polynomial Bₙ(x).

The Bernoulli polynomials are defined by the generating function
t·eˣᵗ/(eᵗ−1) = Σ_{n≥0} Bₙ(x) tⁿ/n!, and the numbers are Bₙ = Bₙ(0).

- **Exact numbers.** A non-negative integer `n` gives the exact rational Bₙ
  from the recurrence B₀ = 1, Bₘ = −1/(m+1) Σ_{k=0}^{m−1} C(m+1, k) Bₖ
  (lazily cached with GMP, so large indices stay exact). Odd `n > 1` give `0`,
  `BernoulliB[0] = 1`, `BernoulliB[1] = -1/2`; e.g. `BernoulliB[12] = -691/2730`.
- **Polynomials.** For a non-negative integer `n`, `BernoulliB[n, x]` expands
  to Σ_{j=0}^{n} C(n, j) B_{n−j} xʲ with exact rational coefficients, staying
  symbolic in `x` (e.g. `BernoulliB[2, z] = 1/6 - z + z^2`) and evaluating
  numerically when `x` is inexact.
- **Numerics.** An inexact integer-valued order evaluates the rational at
  machine (`BernoulliB[2.] = 0.166667`) or arbitrary precision; `BernoulliB`
  can be evaluated to arbitrary numerical precision via `N`, e.g.
  `N[BernoulliB[12], 30] = -0.253113553113553113553113553113`.
- **Listable.** `BernoulliB[{2, 4, 6}] = {1/6, -1/30, 1/42}`; `BernoulliB[{}] = {}`.
- Negative, non-integer, and symbolic orders stay unevaluated.

Attributes: `Listable`, `Protected`.

```mathematica
In[1]:= Table[BernoulliB[k], {k, 0, 10}]
Out[1]= {1, -1/2, 1/6, 0, -1/30, 0, 1/42, 0, -1/30, 0, 5/66}

In[2]:= BernoulliB[3, z]
Out[2]= 1/2 z - 3/2 z^2 + z^3
```

## EulerE

- `EulerE[n]` — the Euler number Eₙ.
- `EulerE[n, x]` — the Euler polynomial Eₙ(x).

The Euler polynomials are defined by the generating function
2·eˣᵗ/(eᵗ+1) = Σ_{n≥0} Eₙ(x) tⁿ/n!, and the numbers are Eₙ = 2ⁿ Eₙ(1/2).

- **Exact numbers.** A non-negative integer `n` gives the exact integer Eₙ
  from the recurrence E₀ = 1, E_{2m} = −Σ_{k=0}^{m−1} C(2m, 2k) E_{2k}
  (lazily cached with GMP, so large indices stay exact). Odd `n` give `0`,
  `EulerE[0] = 1`, `EulerE[2] = -1`, `EulerE[4] = 5`; e.g. `EulerE[10] = -50521`.
- **Polynomials.** For a non-negative integer `n`, `EulerE[n, x]` expands to the
  degree-`n` polynomial in monomial form with exact rational coefficients,
  staying symbolic in `x` (e.g. `EulerE[2, z] = -z + z^2`) and evaluating
  numerically when `x` is inexact. `EulerE[n, 1/2]` with symbolic `n` folds to
  `2^-n EulerE[n]`.
- **Numerics.** An inexact integer-valued order evaluates the integer at
  machine (`EulerE[2.] = -1.`) or arbitrary precision; an Euler polynomial at an
  exact argument can be evaluated to arbitrary numerical precision via `N`, e.g.
  `N[EulerE[6, 1/3], 30] = -0.825788751714677640603566529492`.
- **Listable.** `EulerE[{2, 4, 6}] = {-1, 5, -61}`; `EulerE[{}] = {}`.
- Negative, non-integer, and symbolic orders stay unevaluated.

Attributes: `Listable`, `Protected`.

```mathematica
In[1]:= Table[EulerE[k], {k, 0, 10}]
Out[1]= {1, 0, -1, 0, 5, 0, -61, 0, 1385, 0, -50521}

In[2]:= EulerE[3, z]
Out[2]= 1/4 - 3/2 z^2 + z^3
```

## PolyLog

- `PolyLog[n, z]` — the polylogarithm Liₙ(z) = Σ_{k≥1} zᵏ/kⁿ (|z| < 1; analytic
  continuation elsewhere, with a branch cut from 1 to ∞).
- `PolyLog[n, p, z]` — the Nielsen generalized polylogarithm S_{n,p}(z):
  evaluated numerically where the reference evaluates it (see the Nielsen
  bullet below); exact calls stay symbolic.

Implemented in `src/special_functions/polylog.c`, modelled on `Gamma` / `Zeta`.

- **Exact closed forms.**
  - `PolyLog[n, 0] = 0`.
  - `PolyLog[1, z] = -Log[1 - z]` and `PolyLog[0, z] = z/(1 - z)`.
  - Negative integer orders are Eulerian-number rational functions:
    `PolyLog[-1, z] = z/(1-z)^2`, `PolyLog[-2, z] = (z+z²)/(1-z)^3`,
    `PolyLog[-3, z] = (z+4z²+z³)/(1-z)^4`. (Built from the Eulerian triangle
    with GMP, so exact `z` reduces to an exact rational and symbolic `z` gives a
    rational function.)
  - For integer `n ≥ 2`: `PolyLog[n, 1] = Zeta[n]` (e.g. `PolyLog[2, 1] = Pi²/6`)
    and `PolyLog[n, -1] = (2^(1-n) - 1) Zeta[n]` (e.g. `PolyLog[2, -1] = -Pi²/12`).
  - `PolyLog[2, I] = I Catalan - Pi²/48` and `PolyLog[2, -I] = -I Catalan - Pi²/48`, as the reference evaluates them (v0.458).
  - Other roots of unity stay automatically (as in the reference) and expand
    under `FunctionExpand` (see there): `FunctionExpand[PolyLog[3, I]]` =
    `I Pi³/32 - 3 Zeta[3]/32`.
  - The dilogarithm and trilogarithm at ½:
    `PolyLog[2, 1/2] = Pi²/12 - Log[2]²/2`,
    `PolyLog[3, 1/2] = Log[2]³/6 - Pi² Log[2]/12 + 7 Zeta[3]/8`.
- **Machine numerics.** A machine-precision call with a real order is answered
  by the raw-`double` kernels of `sf_machine_zeta.c` (`sf_machine_polylog` for
  real z, `sf_polylog_c` for complex z), shared bit for bit with the packed /
  `NDArray` kernels and `Compile[]`:
  - the **dilogarithm** `PolyLog[2, z]` (every real and complex z): the
    Bernoulli series Li₂ = Σ Bₙ uⁿ⁺¹/(n+1)! in u = −log(1−z) ('t Hooft–Veltman;
    Crandall), after the maps z → 1/z (|z| > 1) and z → 1−z (Re z > ½), which
    leave |u| ≤ 1.25 and at most ~14 terms; π²/6 is carried as a two-double
    constant. Real x > 1 answers from below the cut, real x < −1 is real;
  - other real orders: closed forms for s = 0, 1; the power series (double-
    double running power, compensated sum) for |z| ≤ ½; the Jonquière expansion
    in μ = log z (complex μ for negative x) for |μ| < 3.6; the duplication
    formula for real x ∈ [−1, −½) and s ≥ 1; and, for an integer order and
    larger |z|, the inversion formula with the Bernoulli polynomial Bₙ. Near an
    integer order the Jonquière head and ζ(s−k) have cancelling poles: a sum
    whose terms exceed 4× its value (16× for integer s) declines to MPFR.
  Accuracy against mpmath (2000 points): Li₂ real x ∈ [−0.9, 0.9] median 0.4 /
  max 2.3 ulp, complex box |Re|, |Im| ≤ 0.8 median 0.4 / max 3.9 ulp, the wide
  box |Re|, |Im| ≤ 20 max 2.6 ulp; Li₃ real median 0.4 / max 4 ulp. Net latency
  ≈ 0.15 µs real, 0.41 µs complex (from 494 / 6560 µs); 10⁶-element
  `PolyLog[2, v]` ≈ 5 ms (SciPy `spence` 27 ms).
- **Arbitrary precision and declines.** With at least one inexact operand,
  arbitrary-precision (MPFR) real and complex inputs — and machine inputs the
  kernels decline (complex order, |z| out of reach, cancellation) — evaluate
  numerically:
  a direct power series for `|z| ≤ 1/2` (a fast pure-real path for real order and
  real −1 < z < 1), and the Jonquière / zeta expansion
  Liₛ(z) = Γ(1−s)(−ln z)^{s−1} + Σ_{k≥0} ζ(s−k) (ln z)ᵏ/k! (with the
  integer-order variant carrying the Hₙ₋₁ − ln(−ln z) term) for
  1/2 < |z| with |ln z| < 2π. At **arbitrary precision** with a real order and real
  −1 < z < 1 (no branch cut), evaluation routes through arb's `acb_polylog`
  (`flint_num_polylog`), ~45× faster than the MPFR series and faster than
  the reference. So does every other z off the cut at arbitrary precision (non-real
  z, or real z ≤ −1): the cut is [1, ∞) only, so the principal value is unique
  there and arb's agrees; `N[PolyLog[2, (1 + I)/2], 40]` takes ~1 ms where the
  MPFR zeta expansion (one `mpfr_zeta` per term) took 0.3–0.6 s. Machine
  precision and the real cut [1, ∞) keep the branch-cut-careful MPFR paths. The required ζ values for a complex order use
  the ζ functional equation in the left half-plane to stay well-conditioned. For the
  common **integer order** on the real machine/packed path, the Jonquière series'
  ζ(s−k) are constants, read from `static const` tables (ζ at integers ≥2, and the
  Bernoulli-rational ζ at integers ≤0) rather than recomputed per element, and the
  defining series uses an integer power for `kˢ` — together ~6× faster on a 10⁶
  buffer (`PolyLog[2, v]` 453→76 ms), ahead of the reference and within a small
  factor of SciPy's specialised `spence` dilog.
  Examples: `PolyLog[2, 0.9] = 1.29971`, `PolyLog[0, 5.0] = -1.25`,
  `N[PolyLog[1, 1/3], 50] = 0.40546510810816438197801311546434913657199042346249`,
  `PolyLog[.2 + I, .5 - I] = -0.0898526 - 0.595865 I`. The branch cut runs along
  the real axis from 1 to ∞; values on it are taken continuous from below,
  matching the reference (e.g. Im Li₂(x) = −π Log[x] for real x > 1, so
  `PolyLog[2, 2.0] = 2.4674 - 2.17759 I`), and the branch point `PolyLog[2, 1.0]`
  is the finite real `Pi²/6`.
- **Derivatives.** `D[PolyLog[n, z], z] = PolyLog[n-1, z]/z`; the order
  derivative is the generic `Derivative[1, 0][PolyLog][n, z]`.
- **Listable.** `PolyLog[{2, 4}, -1] = {-Pi²/12, -7/720 Pi⁴}`;
  `PolyLog[{}, x] = {}`.
- Other exact/symbolic arguments (e.g. `PolyLog[2, 1/3]`, `PolyLog[5, 1/2]`)
  stay unevaluated.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= PolyLog[3, 1/2]
Out[1]= 1/6 Log[2]^3 - 1/12 Log[2] Pi^2 + 7/8 Zeta[3]

In[2]:= PolyLog[2, 0.9]
Out[2]= 1.29971

In[3]:= PolyLog[2, 3.0]
Out[3]= 2.32018 - 3.45139*I
```

- **Nielsen's `PolyLog[n, p, z]`** (v0.477, extended in v0.478), previously
  always symbolic. `PolyLog[n, 1, z] = PolyLog[n + 1, z]`. An inexact call with
  Re p > 0 evaluates exactly where the reference evaluates it, at the call's
  precision (principal logarithms and powers throughout):
  - integer `n ≤ 0`, any `z`: the closed form (z d/dz)^(-n) (-Log[1 - z])^p / Gamma[p + 1];
  - |z| ≤ 1, `z ≠ 1`, any `n` and `p`: the series S_{n,p}(z) = Σ_{k≥p}
    (c(k,p)/k!) zᵏ/kⁿ (Stirling numbers, integer `p`, |z| ≤ 3/4), the
    generalized-power series (1/Γ(p+1)) Σ c_k z^(p+k)/(p+k)ⁿ (Miller's
    recurrence for the coefficients of (-Log(1-w)/w)^p; |z| ≤ 1/2 or Re n ≤ 0),
    or tanh-sinh quadrature of 1/(Γ(n)Γ(p+1)) ∫₀¹ (-Log t)^(n-1) (-Log(1 - z t))^p dt/t;
  - `z = 1`, Re n > 0: the same quadrature;
  - |z| > 1 with integer `n ≥ 1` and a positive integer `p`: the quadrature,
    split at t = 1/z for a real z > 1 (the reference's side of the cut).
  A non-integer `n` or `p` off the unit disk, and `p ≤ 0`, stay unevaluated, as
  in the reference. A real value comes back Real: integer powers of the complex
  logarithm in the `n ≤ 0` closed form are taken by binary powering, not `cpow`
  (`PolyLog[0, 2, 2.] = -4.9348`, = −π²/2, was `-4.9348 - 6.0e-16 I`), on the
  scalar, packed, `NDArray` and `Compile[]` surfaces alike.
  `D[PolyLog[n, p, z], z] = PolyLog[n - 1, p, z]/z`; the `n`- and
  `p`-derivatives are the inert `Derivative[1, 0, 0]` / `Derivative[0, 1, 0]`. Over a grid of 637 machine calls the evaluated set is the
  reference's exactly, and all 479 values are within 6e-15 of 30-digit
  quadrature; the reference's own machine values are off by up to 2e-6 there
  (`PolyLog[3, 0.5, -1.]`). Deliberate difference: `PolyLog[1, 2, 1.]` is
  `Zeta[3]` = 1.2020569…, where the reference says `Indeterminate`.
- **Machine kernel, arrays and `Compile`** (v0.479): `sf_nielsen_c`
  (`src/special_functions/sf_machine_nielsen.c`) is the double-precision twin of
  the arb engine, with the same domain gate. It uses the closed form for
  `n ≤ 0`, the Stirling and Miller series, and tanh-sinh quadrature, and declines
  to the arb path wherever its own error estimate exceeds 2e-14. The scalar
  builtin, the packed / `NDArray` map and `Compile` (an array in `n`, `p` or `z`)
  all call it, so they agree bit for bit. It is within 6.1e-15 of 30-digit values
  over the 479-value grid. About 17 µs a call; a complex value (real `z > 1`)
  hands a compiled real body back to the interpreter.

## LerchPhi

- `LerchPhi[z, s, a]` — the Lerch transcendent Φ(z, s, a) = Σ_{k≥0} zᵏ/(k+a)ˢ
  (|z| < 1; analytic continuation elsewhere, branch cut z ∈ [1, ∞)). For Re a < 0
  the principal value uses the symmetric power ((k+a)²)^(−s/2) and excludes any
  term with k + a = 0.

Options: `DoublyInfinite -> False`, `IncludeSingularTerm -> False`. Both are
validated (v0.476): an unknown name is `LerchPhi::optx`, a non-Boolean value
`LerchPhi::opttf`, a non-option beyond position 3 `LerchPhi::nonopt`; the first
setting of a name wins.

`LerchPhi` is the common generalization of `Zeta`, `HurwitzZeta` and `PolyLog`:
`LerchPhi[1, s, a] = Zeta[s, a]` and `z LerchPhi[z, s, 1] = PolyLog[s, z]`.
Implemented in `src/special_functions/lerchphi.c`, registered via
`lerchphi_init()` in `core_init()`. `Listable`. On the machine path (machine
arguments, default options, `Re a > 0`) it uses the shared `double` kernels of
`sf_machine_zeta.c` — `sf_machine_lerchphi` for real arguments, `sf_lerch_c` for
complex z, s or a — the same helpers the packed / `NDArray` kernels and the
Compile VM use (a complex argument compiles to the complex n-ary kernel opcode
`KERNN_C`). The defining series carries zᵏ in double-double (an FMA residual per
product), forms each k + a by TwoSum with a first-order power correction, sums
with Neumaier compensation, and declines when the terms exceed 16× the sum.
Accuracy against mpmath: real z ∈ [0, 0.9] median 0.3 / max 1.0 ulp, complex
z in the box ±0.6 median 0.2 / max 0.9 ulp; net latency ≈ 1.3 µs real, 1.6 µs
complex (the reference 149 µs). A single `NDArray`/packed argument is threaded over
the real kernel and returns a packed result (`LerchPhi` is on `pack.c`'s `AWARE`
list), and anything the kernels decline (`|z| ≥ 1`, `a ≤ 0`, cancellation, high
precision) falls through to the accurate MPFR path. A visible `NDArray` is never
left unevaluated. At **arbitrary precision** in the convergent region `|z| ≤ 1`
(default convention) evaluation routes through arb's `acb_dirichlet_lerch_phi`
(`flint_num_lerchphi`), ~15× faster than the MPFR series and ~12× faster than
the reference; `|z| > 1` keeps the Erdélyi continuation and `IncludeSingularTerm`
keeps the series.

- **Exact reductions** (built symbolically, so they also numericalize through the
  reused kernels):
  - `LerchPhi[0, s, a] = a^(-s)` (only the k = 0 term survives).
  - `LerchPhi[z, 0, a] = 1/(1 - z)` (the geometric sum, independent of `a`).
  - `LerchPhi[1, s, a] = Zeta[s, a]`; `LerchPhi[z, s, 1] = PolyLog[s, z]/z`.
  - `LerchPhi[-1, 1, a] = (1/2)(PolyGamma[0, (a+1)/2] - PolyGamma[0, a/2])` (the
    finite digamma reduction; the two-`Zeta` form is a `0/0` at `s = 1`). With the
    Gauss digamma theorem this is elementary at rational `a`, e.g.
    `LerchPhi[-1, 1, 1/2] = Pi/2`, `LerchPhi[-1, 1, 3/2] = 2 - Pi/2`.
  - `LerchPhi[-1, s, a]` at half-integer `a` and integer `s ≥ 2` closes via the
    Dirichlet beta function `LerchPhi[-1, s, 1/2] = 2^s β(s)` plus the recurrence
    `Φ(-1,s,b+1) = b^(-s) - Φ(-1,s,b)`: `LerchPhi[-1, 2, 1/2] = 4 Catalan`,
    `LerchPhi[-1, 3, 1/2] = Pi^3/4`, `LerchPhi[-1, 2, 3/2] = 4 - 4 Catalan`. Even
    `s ≥ 4` (no elementary β) falls back to
    `2^(-s) (Zeta[s, 1/4] - Zeta[s, 3/4])`.
  - A positive integer `a` shifts down to the `PolyLog` series, e.g.
    `LerchPhi[z, s, 2] = (PolyLog[s, z] - z)/z^2`.  This reduction is also taken
    for `z = -1` with integer `a`, so `LerchPhi[-1, 1, 1] = Log[2]` (the
    two-`Zeta` `z = -1` form is indeterminate at integer `s`).
  - A non-positive integer `a = -m` drops the singular `k = m` term (the reference's
    default) and shifts onto `PolyLog`:
    `LerchPhi[z, s, 0] = PolyLog[s, z]` (so `LerchPhi[z, s, 0] == PolyLog[s, z]`,
    `LerchPhi[0, s, 0] = 0`, `LerchPhi[0, 0, 0] = 0`), and in general
    `LerchPhi[z, s, -m] = z^m PolyLog[s, z] + Σ_{j=1}^{m} z^(m-j) (-j)^(-s)`.
  - A negative integer `s = -n` gives a rational function of `z`,
    `(z d/dz + a)^n [1/(1-z)]` (e.g. `LerchPhi[2, -1, a] = 2 - a`,
    `LerchPhi[z, -2, a] = (a² + (1+2a+a²) z + (1-2a-2a²) z² + a² z³)/(1-z)³`).
- **Options.** `DoublyInfinite -> True` sums k from −∞ to ∞, realised as
  `Φ(z,s,a) + z⁻¹ Φ(1/z, s, 1-a)` (in a symmetric case it just doubles, e.g.
  `LerchPhi[1, 2, 0.5, DoublyInfinite -> True] = Pi² = 9.8696`).
  `IncludeSingularTerm -> True` keeps the `k + a = 0` term, giving
  `ComplexInfinity` at a non-positive integer `a`.
- **Numerics.** With an inexact operand and `|z| < 1`, a complex-MPFR power
  series evaluates the value at machine or arbitrary precision (real and
  complex), e.g. `LerchPhi[0.5, 3, 2.5] = 0.0794983`,
  `LerchPhi[0.3 + 0.2 I, 2, 1.5] = 0.495505 + 0.0444653 I`,
  `N[LerchPhi[1/2, 2, 5/2], 30] = 0.219693113434910235649039949138…`.
  For `|z| > 1` off the branch cut `[1, ∞)` the analytic continuation comes
  from arb's rigorous `acb_dirichlet_lerch_phi` (any `s`, any `|Log z|`), with
  the Lerch/Erdélyi expansion
  `Φ(z,s,a) = z^(-a)[Γ(1-s)(-Log z)^(s-1) + Σ_{n≥0} HurwitzZeta(s-n,a)(Log z)^n/n!]`
  as the fallback where it is valid (`s` not a positive integer, `|Log z| < 2π`).
  Both are principal-branch, so for `Re a < 0` the symmetric head
  `Σ_{k<m} zᵏ ((k+a)²)^(-s/2)` is summed first and the continuation applied to
  `zᵐ Φ(z, s, a+m)`, `m = ⌈-Re a⌉`. On the negative real axis the value is real:
  `LerchPhi[-3.0, 1.5, 0.5] = 2.06336`,
  `N[LerchPhi[-3, 3/2, 1/2], 25] = 2.063359848206451659906464…`,
  `LerchPhi[-1.5, 1.5, -0.5] = -0.653307`, `LerchPhi[5. + I, I, I + 2] =
  -0.581502 + 0.384767 I`. (The real-axis values used to be wrong —
  `6.35096 + 1.44603 I` — because the Erdélyi terms `HurwitzZeta[s-n, a]` had
  lost their digits to cancellation, and `Re a < 0` returned the principal
  value.) On the cut `(1, ∞)` (v0.476) arb's continuation gives the limit
  from below, Φ(x − i0, s, a), the convention `PolyLog` uses
  (`LerchPhi[2., 2., 1.] = Pi^2/8 - I Pi Log[2]/2`, `LerchPhi[3., 2., 0.5] =
  4.30831 - 1.99266 I`, `LerchPhi[2., 3, 0.5] = 8.92139 - 0.533649 I`). These
  agree with the reference within 0.46 ulp of a 40-digit value, where the
  reference's own machine values are off by up to 2130 ulp (complex `a`).
  Machine cost ≈ 1.1 ms (arb's |z| > 1 algorithm; 0.75 ms is its floor at 60
  bits), about 3× the reference; `DoublyInfinite` at |z| < 1 adds one such call.
  Without a head shift the operands go to arb at their own precision.
  Exact `z > 1` stays symbolic (`LerchPhi[2, 2, 1/2]`), as in the reference.
- **`a ≤ 0` on the machine path** (v0.476): the same shifted head as `Zeta`,
  Σ_{k<m} zᵏ |k+a|⁻ˢ + zᵐ Φ(z, s, a+m), in the real kernel (m ≤ 10⁵), so
  `LerchPhi[0.5, 2., -1.5]` no longer goes through MPFR, and `z = 1` is
  `Zeta[s, a]` including `a ≤ 0` (`LerchPhi[1., 2., -1.5] = 9.37925`).
- **Series at `z = 0`.** `Series[LerchPhi[z, s, a], {z, 0, n}]` returns the
  defining power series, coefficient of `z^k` being `((k+a)^2)^(-s/2)`:
  `Series[LerchPhi[z, s, a], {z, 0, 5}] // Normal =
  (a^2)^(-s/2) + ((1+a)^2)^(-s/2) z + …`.
- **Derivatives.** `D[LerchPhi[z, s, a], z] = (LerchPhi[z, -1+s, a]
  - a LerchPhi[z, s, a])/z` and `D[LerchPhi[z, s, a], a] = -s LerchPhi[z, 1+s, a]`;
  the `s`-derivative is the generic `Derivative[0, 1, 0][LerchPhi][z, s, a]`.
- **Diagnostics.** Fewer than three arguments emit `LerchPhi::argrx`; a
  non-option beyond position 3 emits `LerchPhi::nonopt`, an unknown option
  `LerchPhi::optx`, a non-Boolean option value `LerchPhi::opttf`; all leave the
  call unevaluated.
- **Deliberate differences.** An exact `DoublyInfinite` call with a closed form
  evaluates (`LerchPhi[1, 2, 1/2, DoublyInfinite -> True] = Pi^2`) where the
  reference leaves it unevaluated; the value agrees with its `N`. An exact
  irrational operand next to an inexact one (`LerchPhi[Exp[I], 2., 1.5]`) stays
  unevaluated, as for every special function here (no mixed-argument
  numericization yet). The compiled form is scalar-only, like every N-ary
  kernel (`WhittakerM`, `HypergeometricU`, …).

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= LerchPhi[z, s, 1]
Out[1]= PolyLog[s, z]/z

In[2]:= LerchPhi[0.5, 3, 2.5]
Out[2]= 0.0794983

In[3]:= LerchPhi[-1, 2, 1/2]
Out[3]= 4 Catalan
```

## DirichletL

- `DirichletL[k, j, s]` — the Dirichlet L-function L(χ, s) for the character
  χ = `DirichletCharacter[k, j, ·]` (see `number-theory.md`); for Re s > 1 it is
  Σ_{n≥1} χ(n) n⁻ˢ, continued to the whole s-plane (entire unless χ is
  principal, which has a pole at s = 1).

Implemented in `src/special_functions/dirichlet_l.c` on the character engine
`dirichlet.{c,h}`, registered via `dirichletl_init()` in `core_init()`.
Attributes `Listable`, `NHoldAll`, `Protected` (v0.476; it was `NumericFunction`, so `NumericQ[DirichletL[4, 2, 3]]` is now `False` and `1. + DirichletL[4, 2, 3]` stays unevaluated, as in the reference). `N` keeps the labels k, j
exact (`numericalize_rec`, as the reference's `NHoldAll` + N-rule does), so
`N[DirichletL[3, 2, 1/2], 50]` evaluates at s = 0.5`50.

Machine accuracy against a 40-digit reference is within 3.4 ulp
(`DirichletL[5, 2, 1.5]` 3.0 ulp; the reference itself is about 11 ulp off
there). Two deliberate differences: `DirichletL[4, 2, 1.]` gives `0.785398`
(π/4) and `DirichletL[4, 2, 1.000001]` is correct to the last digit where the
reference leaves the first unevaluated and is off by 2·10⁻¹¹ on the second;
`DirichletL[1, 1, 1.]` is `ComplexInfinity` (as `Zeta[1.]` is) where the
reference stays unevaluated.

- **Messages.** `argrx` (arity), `intp` (a number that is not a positive integer
  at position 1 or 2), `invl` (j > EulerPhi[k]); all through `mth_message`.
- **k = 1** is `Zeta[s]` on every surface (bit-identical machine values).
- **Exact values** (computed in the cyclotomic field Q(ζ_ord χ) as coefficient
  vectors reduced modulo the cyclotomic polynomial, rendered as rational
  combinations of roots of unity in `DirichletCharacter`'s normal form):
  - s = 1 − m ≤ 0, every character: L(χ, 1−m) = −B_{m,χ*}/m · Π_{p|k, p∤f}
    (1 − χ*(p) p^{m−1}), with the generalized Bernoulli numbers
    B_{m,χ} = Σ C(m,i) B_i f^{i−1} Σ_a χ(a) a^{m−i} of the primitive χ* mod the
    conductor f. E.g. `DirichletL[3, 2, -4] = 2/3`,
    `DirichletL[5, 2, -4] = 148/25 + 86/25 I`, `DirichletL[12, 3, -4] = 205`.
  - s = 1, non-principal χ: the classical closed forms of the primitive χ*
    (odd: π i τ(χ*)/f² Σ χ̄(a) a; even: −2 τ(χ*)/f Σ_{a<f/2} χ̄(a) Log[Sin[π a/f]];
    τ = √f or i√f for a real χ*) times Π (1 − χ*(p)/p). `DirichletL[3, 2, 1] =
    Pi/(3 Sqrt[3])`, `DirichletL[4, 2, 1] = Pi/4`. The reference leaves the
    imprimitive s = 1 values unevaluated; Mathilda gives them.
  - Everything else (s ≥ 2 for k > 1, non-integer exact s, the principal
    character at s = 1 for k > 1) stays unevaluated, as in the reference.
  - Matches the reference on all 1152 exact values it produces for
    k ∈ {2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 16, 20, 21, 24, 25, 27, 32, 35, 40},
    every j, s ∈ −5..1.
- **Machine precision** (`src/special_functions/sf_machine_dirichlet.c`, the one
  raw-double kernel behind the scalar, packed/NDArray and Compile surfaces).
  Reduce χ to χ* and Euler factors; conductor 1 is `sf_zeta_*`. For Re s ≥ 1/2:
  Euler–Maclaurin on L itself — the direct sum to n = M f plus, per residue a,
  the E-M tail of f⁻ˢ ζ(s, M + a/f) with J = 10 Bernoulli terms and M the
  smallest integer meeting the remainder bound. The integral terms are summed
  as −M^{1−s} Σ_a χ(a) ℓ_a φ((1−s)ℓ_a), ℓ_a = log1p(a/(fM)), φ(u) = expm1(u)/u,
  exact because Σ χ(a) = 0: no cancellation as s → 1, and s = 1 is φ(0) = 1. A
  real character at real s runs in real arithmetic. For Re s < 1/2: the
  functional equation with the root number ε = τ(χ*)/(i^δ √f) from the
  character table (ε = 1 for real χ*), exact 0 at the trivial zeros.
  Accuracy: max 6e-14 relative over 3000 random (k ≤ 300, s real in ±12,
  complex |Im s| ≤ 40) against arb; the reference's own machine values are off by
  up to 1e-9 relative at some of these points (e.g. (7, 2, −7/2 + 10 I)).
  Shootout (`benchmarks/97-special-functions`, v0.450): 1.5 µs (real χ mod 5) /
  2.6 µs (complex χ mod 7) per machine call vs the reference 288 / 462 µs and
  mpmath 12 / 18 ms; complex s 3.0 / 4.2 µs vs 425 / 660 µs; a 10⁶ buffer 1.0 /
  1.8 s vs 171 / 289 s; ~35 µs at k = 1000.
- **Arbitrary precision.** arb's `acb_dirichlet_l` on the primitive χ* (set on
  FLINT's own generators from our values, never via Conrey labels) times the
  Euler factors, with the bridge's precision-doubling ladder
  (`flint_num_dirichlet_l`); the result carries the precision of s. 0.73 / 1.09
  ms at 50 digits (the reference 8.9 / 11.2 ms, mpmath 28 / 58 ms). Without FLINT, the HurwitzZeta template
  f⁻ˢ Σ χ*(a) ζ(s, a/f) × Euler factors.
- **Interval** (real characters). A rigorous arb enclosure
  (`flint_dirichlet_l_interval`): pieces whose derivative ball excludes 0 are
  monotone and bounded by their endpoint values; others are bisected (depth
  12) and bounded at the leaves by L(mid) ± r sup|L'|. The jets come from
  deflated Hurwitz series (L = q⁻ˢ Σ χ(a)(ζ(s, a/q) − 1/(s−1)), entire for
  non-principal χ) or, below s = −1/4, the functional equation. Tighter than
  the reference: `DirichletL[5, 1, Interval[{1.23, 1.24}]] =
  Interval[{4.11396, 4.25899}]` (the endpoint values; the reference returns the
  looser {4.09145, 4.28087}), and non-principal characters are enclosed where
  the reference returns unevaluated. A range containing s = 1 is left unevaluated
  (arb's sums have no finite enclosure on a ball through the cancelled poles).
  Exact endpoints (Integer, Rational, or MPFR) are first rounded outward to the
  enclosing doubles (GMP comparison plus `nextafter`), so
  `DirichletL[5, 3, Interval[{1/2, 3/4}]] = Interval[{0.231751, 0.336106}]`
  (previously unevaluated; the reference returns the much looser
  `{-12.92, 13.54}` with dyadic-rational endpoints).
  `CenteredInterval` is not implemented in Mathilda.
- **Surfaces.** A packed / `NDArray` s with scalar k, j threads the machine
  kernel over the buffer with the character table fetched once (`AWARE`; a
  real character over a real buffer stays real, else complex128; an int64
  buffer is exact s and goes to the scalar path). `Compile[]` lowers through the
  n-ary kernel `ndk_DirichletL_n` (scalar `KERNN`/`KERNN_C`, a tile in the s
  slot `VKERNN_R`/`_C`), typed Real only for literal labels of a real character
  at real s. Scalar, unpacked List, packed, `NDArray` and both Compile shapes
  agree bit-for-bit.

```mathematica
In[1]:= DirichletL[4, 2, 1]
Out[1]= 1/4 Pi

In[2]:= DirichletL[12, 3, 1]
Out[2]= 1/3 Pi

In[3]:= {DirichletL[4, 2, 0], DirichletL[4, 2, -2], DirichletL[4, 2, -3]}
Out[3]= {1/2, -1/2, 0}

In[4]:= N[DirichletL[4, 2, 2], 30]
Out[4]= 0.9159655941772190150546035149325
```

## AiryAi

- `AiryAi[z]` — the Airy function Ai(z), the solution of y″ = z y that decays as
  z → +∞. An entire function of z (no branch cuts).

Implemented in `src/special_functions/airyai.c`, modelled on `Erf` (an entire
function with a file-local complex-MPFR toolkit `acx` and a unified core that
returns Ai(z) and Ai′(z) together).

- **Exact values.** `AiryAi[0] = 1/(3^(2/3) Gamma[2/3])`;
  `AiryAi[Infinity] = AiryAi[-Infinity] = 0`. Other exact arguments (e.g.
  `AiryAi[2]`) stay unevaluated and numericalize only under `N`.
- **Machine precision.** Machine real and complex arguments are answered by the
  raw-double kernels in `src/special_functions/sf_machine_bessel.c`, shared bit
  for bit with the NDArray/packed path and `Compile[]`: on the real line,
  piecewise Chebyshev polynomials on [−10, 10] (the functions themselves on
  [−10, 2], the exponentially scaled ones on [2, 10]; generated by
  `tools/sf_accuracy/gen_airy_tables.py`) and the Poincaré expansions beyond,
  with ζ = (2/3)x^{3/2} in double-double; for complex z, the Maclaurin series
  for |z| ≤ 1 and otherwise DLMF 9.6.1–9.6.3 through K and I of order ±1/3, 2/3
  (one core call), continued across Re ζ < 0. Accuracy: ≤ 4 ulp of the local
  magnitude on the real line, median ≈ 1 and max ≤ 10 for complex arguments.
  Values beyond the double range decline to the arbitrary-precision path
  (`AiryAi[120.] = 2.1575×10^-382`, previously flushed to `0.`).
- **Numerics.** With an inexact argument, machine- and arbitrary-precision
  (MPFR) real and complex inputs evaluate numerically. At **arbitrary precision**
  all four Airy heads route through arb's `acb_hypgeom_airy` (`flint_num_airy_*`),
  ~9× faster than the MPFR core and faster than the reference (Airy is entire, so
  arb's values match directly); machine `double` keeps the unified `acx` core,
  which routes on `r = |z|` and the requested precision `P`:
  - a **Maclaurin series** for small/moderate `r` (recurrence
    `a_n = a_{n-3}/(n(n-1))` from Ai″ = z Ai, with `(2/3) r^{3/2}/ln2` guard bits
    to absorb the partial-sum cancellation), and
  - the **asymptotic series** (DLMF 9.7.5) for large `r`, with
    ζ = (2/3) z^{3/2} and `u_k = (6k-5)(6k-3)(6k-1)/((2k-1) 216 k) u_{k-1}`,
    summed to optimal truncation; near the negative real axis
    (2π/3 < |arg z| ≤ π) it applies the DLMF 9.2.12 connection relation
    `Ai(z) = -[ω Ai(ωz) + ω̄ Ai(ω̄z)]`, ω = e^{2πi/3}, so the oscillation falls
    out of two well-conditioned in-sector evaluations.

  Examples: `AiryAi[1.8] = 0.0470362`, `AiryAi[2.5 + I] = -0.00191209 -
  0.0180329 I`,
  `N[AiryAi[2], 50] = 0.034924130423274379135322080791807609761060213897583`.
  Precision tracks the input (`AiryAi[4.8\`40]` is accurate to 40 digits).
  Large imaginary arguments whose magnitude overflows the machine range are
  returned as arbitrary-precision numbers (e.g. `AiryAi[150. I] ≈ 8.1951×10^374
  + 6.3801×10^374 I`).
- **Derivative.** `D[AiryAi[z], z] = AiryAiPrime[z]`; higher derivatives reduce
  via Ai″ = z Ai (e.g. `D[AiryAi[x], {x, 2}] = x AiryAi[x]`).
- **Series at 0.** `Series[AiryAi[x], {x, 0, n}]` gives the closed-form Taylor
  series, e.g. `1/(3^(2/3) Gamma[2/3]) - x/(3^(1/3) Gamma[1/3]) + …`.
- **Series at ∞.** `Series[AiryAi[x], {x, Infinity, n}]` gives the asymptotic
  expansion with the essential-singularity prefactor kept symbolic:
  `E^(-2/3 x^(3/2)) ((1/x)^(1/4)/(2 Sqrt[Pi]) - 5 (1/x)^(7/4)/(96 Sqrt[Pi]) + …)`.
- **Listable.** `AiryAi[{1.2, 1.5, 1.8}] = {0.106126, 0.0717495, 0.0470362}`;
  `AiryAi[{}] = {}`.

`AiryAiPrime[z]` = Ai′(z) is a full numeric evaluator in its own right. Because
the unified core returns Ai(z) and Ai′(z) together, `AiryAiPrime` reuses the very
same Maclaurin / asymptotic / connection machinery and simply selects the
derivative component:

- **Exact values.** `AiryAiPrime[0] = -1/(3^(1/3) Gamma[1/3])` and
  `AiryAiPrime[Infinity] = 0` (Ai′ decays). At −∞ Ai′ oscillates with *growing*
  (~|z|^(1/4)) amplitude and has no limit, so `AiryAiPrime[-Infinity]` is left
  unevaluated. Other exact arguments stay symbolic and numericalize under `N`.
- **Numerics.** Real and complex inputs evaluate at machine and arbitrary (MPFR)
  precision, with the same overflow→MPFR promotion as `AiryAi`. Examples:
  `AiryAiPrime[0.5] = -0.224911`, `AiryAiPrime[2.5 + I] = -0.00187921 +
  0.0310276 I`,
  `N[AiryAiPrime[5/2], 50] = -0.026250881035903230364895496297232509446317838135771`.
- **Derivative.** `D[AiryAiPrime[z], z] = z AiryAi[z]` (from Ai″ = z Ai).
- **Series at 0.** `Series[AiryAiPrime[x], {x, 0, n}]` via the closed-form Taylor
  path: `-1/(3^(1/3) Gamma[1/3]) + x^2/(2·3^(2/3) Gamma[2/3]) - …`.
- **Series at ∞.** `Series[AiryAiPrime[x], {x, Infinity, n}]` gives the DLMF 9.7.6
  expansion `E^(-2/3 x^(3/2)) (-x^(1/4)/(2 Sqrt[Pi]) - 7 (1/x)^(5/4)/(96 Sqrt[Pi])
  + …)`.
- **Listable.** `AiryAiPrime[{1.2, 1.5, 1.8}] = {-0.132785, -0.097382, -0.0685248}`.

Attributes (both heads): `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= AiryAi[0]
Out[1]= 1/(3^(2/3) Gamma[2/3])

In[2]:= AiryAi[1.8]
Out[2]= 0.0470362

In[3]:= AiryAi[2.5 + I]
Out[3]= -0.00191209 - 0.0180329*I
```

## AiryBi

- `AiryBi[z]` — the Airy function Bi(z), the second solution of y″ = z y; it
  grows exponentially as z → +∞ and oscillates with decaying amplitude as
  z → −∞. An entire function of z (no branch cuts), the companion of `AiryAi`.

Implemented in `src/special_functions/airybi.c`, modelled on `AiryAi` (a
file-local complex-MPFR toolkit `acx` and a unified core that returns Bi(z) and
Bi′(z) together).

- **Exact values.** `AiryBi[0] = 1/(3^(1/6) Gamma[2/3])`;
  `AiryBi[Infinity] = Infinity`, `AiryBi[-Infinity] = 0`. Other exact arguments
  (e.g. `AiryBi[2]`) stay unevaluated and numericalize only under `N`.
- **Machine precision.** Machine real and complex arguments are answered by the
  raw-double kernels in `src/special_functions/sf_machine_bessel.c`, shared bit
  for bit with the NDArray/packed path and `Compile[]`: on the real line,
  piecewise Chebyshev polynomials on [−10, 10] (the functions themselves on
  [−10, 2], the exponentially scaled ones on [2, 10]; generated by
  `tools/sf_accuracy/gen_airy_tables.py`) and the Poincaré expansions beyond,
  with ζ = (2/3)x^{3/2} in double-double; for complex z, the Maclaurin series
  for |z| ≤ 1 and otherwise DLMF 9.6.1–9.6.3 through K and I of order ±1/3, 2/3
  (one core call), continued across Re ζ < 0. Accuracy: ≤ 4 ulp of the local
  magnitude on the real line, median ≈ 1 and max ≤ 10 for complex arguments.
  Values beyond the double range decline to the arbitrary-precision path
  (`AiryAi[120.] = 2.1575×10^-382`, previously flushed to `0.`).
- **Numerics.** With an inexact argument, machine- and arbitrary-precision
  (MPFR) real and complex inputs evaluate numerically. The unified core routes
  on `r = |z|`, `θ = arg z`, and the requested precision `P`:
  - a **Maclaurin series** for small/moderate `r` (recurrence
    `b_n = b_{n-3}/(n(n-1))` from Bi″ = z Bi, seeded by `Bi(0)`/`Bi′(0)`, with
    `(2/3) r^{3/2}/ln2` guard bits to absorb the partial-sum cancellation),
  - the **dominant asymptotic series** (DLMF 9.7.7) for large `r` in the central
    sector, with ζ = (2/3) z^{3/2}, the same `u_k`/`v_k` as Ai (no `(-1)^k`,
    prefactor `1/Sqrt[Pi]`), used only where the recessive companion
    `~exp(-2 Re ζ)` is below 2^−P (`Re ζ = (2/3) r^{3/2} cos(3θ/2) > (P ln2)/2`),
    which keeps the exponentially large positive axis at full precision, and
  - near/left of the anti-Stokes line `|arg z| = π/3` (incl. the oscillatory
    negative real axis) the **DLMF 9.2.10 connection to Ai**,
    `Bi(z) = e^{iπ/6} Ai(z e^{2πi/3}) + e^{-iπ/6} Ai(z e^{-2πi/3})`, so the
    oscillation falls out of two well-conditioned in-sector Ai evaluations.

  Examples: `AiryBi[1.8] = 2.59587`, `AiryBi[2.5 + I] = 0.512544 + 5.335 I`,
  `N[AiryBi[2], 50] = 3.2980949999782147102806044252234524220039759634036`.
  Precision tracks the input (`AiryBi[1.8\`30]` is accurate to 30 digits).
  Machine inputs whose magnitude overflows the machine range are returned as
  arbitrary-precision numbers (e.g. `AiryBi[1000.] ≈ 5.4077×10^9154`).
- **Derivative.** `D[AiryBi[z], z] = AiryBiPrime[z]`; higher derivatives reduce
  via Bi″ = z Bi (e.g. `D[AiryBi[x], {x, 2}] = x AiryBi[x]`).
- **Series at 0.** `Series[AiryBi[x], {x, 0, n}]` gives the closed-form Taylor
  series, e.g. `1/(3^(1/6) Gamma[2/3]) + 3^(1/6) x/Gamma[1/3] + …`.
- **Series at ∞.** `Series[AiryBi[x], {x, Infinity, n}]` gives the asymptotic
  expansion with the essential-singularity prefactor kept symbolic:
  `E^(2/3 x^(3/2)) ((1/x)^(1/4)/Sqrt[Pi] + 5 (1/x)^(7/4)/(48 Sqrt[Pi]) + …)`.
- **Listable.** `AiryBi[{1.2, 1.5, 1.8}] = {1.42113, 1.87894, 2.59587}`;
  `AiryBi[{}] = {}`.

`AiryBiPrime[z]` = Bi′(z) is a full numeric evaluator in its own right. Because
the unified core returns Bi(z) and Bi′(z) together, `AiryBiPrime` reuses the very
same Maclaurin / asymptotic / connection machinery and simply selects the
derivative component:

- **Exact values.** `AiryBiPrime[0] = 3^(1/6)/Gamma[1/3]` and
  `AiryBiPrime[Infinity] = Infinity` (Bi′ is the dominant, growing solution). At
  −∞ Bi′ oscillates with *growing* (~|z|^(1/4)) amplitude and has no limit, so
  `AiryBiPrime[-Infinity]` is left unevaluated. Other exact arguments stay
  symbolic and numericalize under `N`.
- **Numerics.** Real and complex inputs evaluate at machine and arbitrary (MPFR)
  precision, with the same overflow→MPFR promotion as `AiryBi`. Examples:
  `AiryBiPrime[1.8] = 2.98554`, `AiryBiPrime[2.5 + I] = -1.20505 + 8.29097 I`,
  `N[AiryBiPrime[5/2], 50] = 9.4214233173343017555823088857282415621646345227564`.
- **Derivative.** `D[AiryBiPrime[z], z] = z AiryBi[z]` (from Bi″ = z Bi).
- **Series at 0.** `Series[AiryBiPrime[x], {x, 0, n}]` via the closed-form Taylor
  path: `3^(1/6)/Gamma[1/3] + x^2/(2·3^(1/6) Gamma[2/3]) + …`.
- **Series at ∞.** `Series[AiryBiPrime[x], {x, Infinity, n}]` gives the DLMF 9.7.8
  expansion `E^(2/3 x^(3/2)) (x^(1/4)/Sqrt[Pi] - 7 (1/x)^(5/4)/(48 Sqrt[Pi]) + …)`.
- **Listable.** `AiryBiPrime[{1.2, 1.5, 1.8}] = {1.22123, 1.88621, 2.98554}`.

Attributes (both heads): `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= AiryBi[0]
Out[1]= 1/(3^(1/6) Gamma[2/3])

In[2]:= AiryBi[1.8]
Out[2]= 2.59587

In[3]:= AiryBi[2.5 + I]
Out[3]= 0.512544 + 5.335*I
```

## BesselJ

- `BesselJ[n, z]` — the Bessel function of the first kind Jₙ(z), the solution of
  z² y″ + z y′ + (z² − n²) y = 0 that is regular at the origin. For non-integer
  order it has a branch cut along the negative real z axis (from the (z/2)ⁿ
  factor).

Implemented in `src/special_functions/bessel.c`. The complex-MPFR arithmetic
is the reusable `ncpx` toolkit (`src/numeric_complex.{c,h}`) — the shared
successor to the file-local `acx`/`ecx`/`gcx` toolkits, so future
`BesselY`/`BesselI`/`BesselK` reuse it.

- **Exact values at the origin.** From the leading `(z/2)^ν/Γ(ν+1)` behaviour:
  `BesselJ[0, 0] = 1`; `0` for integer n ≠ 0 (either sign) or non-integer Re(ν) > 0
  (e.g. `BesselJ[1/2, 0] = 0`); `ComplexInfinity` for non-integer Re(ν) < 0 (e.g.
  `BesselJ[-1/2, 0] = ComplexInfinity`); `Indeterminate` for pure-imaginary order.
  A symbolic order at the origin stays unevaluated. Other exact / symbolic
  arguments stay unevaluated and numericalize only under `N` (e.g.
  `N[BesselJ[0, 4], 50]`).
- **Machine precision.** A machine-number call with a **real order** (machine or
  exact) is answered by the raw-double kernels in
  `src/special_functions/sf_machine_bessel.c` — the same helper the
  NDArray/packed element-wise path and `Compile[]` call, so the scalar, List,
  packed, `NDArray` and compiled answers agree bit for bit. A real argument
  gives a `Real` where the value is real (and a `Complex` on the cut, e.g.
  `BesselJ[1/2, -2.]`); a complex argument gives a `Complex`. Regions:
  - real axis: integer-order J/Y from libm `jn`/`yn`; I₀, I₁, K₀, K₁ from
    Cephes-style regions (series / piecewise Chebyshev fits of the
    exponentially scaled function / Hankel), Kₙ by upward recurrence;
  - everything else by Amos's scheme (ACM TOMS 644): J, Y through I and K of
    u = −iz in the right half plane, with K from Temme's series (|u| ≤ 1.25),
    Taylor continuation from tabulated nodes (1.25 ≤ |u| < 21, orders 0 and
    ±1/3 mod 1; Steed's CF2 for other fractional orders) or Hankel's expansion;
    I from Miller's backward recurrence normalised by a Gegenbauer/Neumann sum.
  Accuracy against mpmath (`tools/sf_accuracy/bessel.py`, ulps of the local
  magnitude): ≤ 3.6 on the real axis, median ≈ 1 and max ≤ 8 for complex
  arguments on the benchmark domain (≤ 12 in a wide band), and median ≈ 1.5,
  max ≤ ~20 for non-integer orders. A kernel
  decline (|ν| > 1000, overflow/underflow, a complex order) falls through to
  the MPFR path below, unchanged.
- **Numerics.** When some argument is inexact, machine- and arbitrary-precision
  (MPFR) real and complex order and argument evaluate numerically. The core
  routes on `r = |z|` and the requested precision `P`:
  - integer order with a real argument takes MPFR's correctly-rounded
    `mpfr_jn` fast path (J₋ₙ = (−1)ⁿ Jₙ);
  - a **power series** for small/moderate `r` (DLMF 10.2.2,
    `t_k = t_{k-1}·(−(z/2)²)/(k(ν+k))`, one `Gamma(ν+1)` via `mpfr_gamma` for
    real order or the `Gamma` builtin for complex order, with ~`2r/ln2` guard
    bits to absorb the alternating-series cancellation);
  - an **asymptotic (Hankel) series** for large `r` away from the negative-real
    axis (DLMF 10.17.3), used only where a double-precision walk of its term
    sizes shows that no term exceeds the leading one and some term reaches
    2^-(P+16) (roughly ν² ≪ r; a half-integer order terminates and is exact);
  - for a **large order at large `r`** (the Hankel walk rejects ν but accepts
    its fractional part μ = ν − ⌊Re ν⌋), the **three-term recurrence in the
    order** from the Hankel values at μ and μ+1: Miller's backward recurrence
    for J and I (start chosen by Olver's forward test, normalised to the Hankel
    pair by least squares), upward for Y and K; a large negative order goes
    through the reflection formulas with a measured cancellation retry. This is
    O(max(ν, r)) steps where the power series needed ~2r/ln2 guard bits.
    Integer-order Y with a real argument and n > 32 uses the upward recurrence
    from `mpfr_y0`/`mpfr_y1` instead of `mpfr_yn` (`BesselY[1500., 5000.]` took
    27 s). The same routing applies to `BesselY`, `BesselI` and `BesselK`.
    Through v0.460 the Hankel series was chosen on `r` alone: `N[BesselJ[51/2,
    100], 20]` gave −0.2001 (true 0.0745832058174051566889) and
    `BesselJ[500.5, 1000.]`, which the machine kernel hands to this path, gave
    1.798 (true −0.00683360793869777).

  Examples: `BesselJ[0, 5.2] = -0.11029`, `BesselJ[7/3 + I, 4.5 - I] =
  1.18908 + 0.715653 I`, `N[BesselJ[0, 4], 50] =
  -0.39714980986384737228659076845169804197561868528939`. Precision tracks the
  input (`BesselJ[0, 4.000000000000000000000000]` carries 25 digits).
- **Half-integer order.** With a **symbolic** argument, `BesselJ[(2k+1)/2, x]`
  rewrites to the elementary spherical-Bessel closed form (e.g.
  `BesselJ[1/2, x] = Sqrt[2/(Pi x)] Sin[x]`); with a numeric argument the call
  uses the accurate general numeric path instead (e.g.
  `BesselJ[35/2, 1.] = 3.55153×10⁻²¹`), and an exact numeric argument such as
  `BesselJ[11/2, 1]` is left unevaluated. These rewrites and the
  negative-integer reflection `BesselJ[-n, x] = (-1)ⁿ BesselJ[n, x]` live in
  the C builtin (`bessel_symbolic_rules` in `bessel.c`), tried only after the
  numeric paths decline; they were formerly DownValues in
  `src/internal/bessel.m`, which every numeric call had to match first.
- **Argument parity.** For **integer** order, Jₙ is entire and obeys
  `BesselJ[n, -z] = (-1)ⁿ BesselJ[n, z]` (even in z for even n, odd for odd n),
  applied whenever the argument is superficially negative — e.g.
  `BesselJ[0, -z] = BesselJ[0, z]`, `BesselJ[1, -z] = -BesselJ[1, z]`,
  `BesselJ[1, -2 x] = -BesselJ[1, 2 x]`. Folded in the C builtin (reusing
  `expr_is_superficially_negative`) only for symbolic / exact arguments; a
  concrete inexact argument is evaluated numerically instead. Non-integer and
  symbolic orders are left unevaluated (no parity across the branch cut).
- **Derivative.** `D[BesselJ[n, x], x] = (BesselJ[n-1, x] - BesselJ[n+1, x])/2`
  (DLMF 10.6.1).
- **Series at 0.** `Series[BesselJ[0, x], {x, 0, 10}] =
  1 - x²/4 + x⁴/64 - x⁶/2304 + …` (Puiseux for half-integer order).
- **Series at ∞.** `Series[BesselJ[n, x], {x, Infinity, k}]` gives the
  asymptotic expansion `Cos[(2n+1)π/4 - x] (Sqrt[2/Pi]/Sqrt[x] + …) +
  Sin[(2n+1)π/4 - x] (…)`, valid for symbolic order too.
- **Listable.** `BesselJ[1, {0.5, 1.0, 1.5}] = {0.242268, 0.440051, 0.557937}`.

Attributes: `Listable`, `NumericFunction`, `Protected`.

Not yet implemented: `Limit[BesselJ[n, x], x -> Infinity]` is left unevaluated
(the limit engine does not yet consume the oscillatory asymptotic leading term);
the symbolic-index general term `SeriesCoefficient[BesselJ[0, x], {x, 0, n}]`
(the Piecewise closed form) — concrete-index `SeriesCoefficient` works.

```mathematica
In[1]:= BesselJ[0, 5.2]
Out[1]= -0.11029

In[2]:= D[BesselJ[n, x], x]
Out[2]= 1/2 (BesselJ[-1 + n, x] - BesselJ[1 + n, x])

In[3]:= BesselJ[1, -2 x]
Out[3]= -BesselJ[1, 2 x]
```

## BesselK

- `BesselK[n, z]` — the modified Bessel function of the second kind Kₙ(z), a
  solution of z² y″ + z y′ − (z² + n²) y = 0. It is even in the order
  (K₋ₙ = Kₙ), decays like e⁻ᶻ as z → ∞, and has a branch cut along the negative
  real z axis.

Implemented in `src/special_functions/bessel.c` (the consolidated Bessel module,
which also houses `BesselJ`), reusing the `ncpx` complex-MPFR toolkit
(`src/numeric_complex.{c,h}`). Unlike `BesselJ` there is **no** MPFR library
routine for the modified Bessel functions, so every value is summed from scratch.

- **Exact values at the origin.** Kₙ diverges for every order with nonzero real
  part: `BesselK[0, 0] = Infinity`; `BesselK[n, 0] = ComplexInfinity` otherwise
  (integer, rational, or symbolic). A **pure-imaginary** order is the lone
  exception — `K_{iτ}` oscillates with bounded amplitude as z → 0, so
  `BesselK[I, 0] = Indeterminate`. Other exact / symbolic arguments stay
  unevaluated and numericalize only under `N`.
- **Machine precision.** A machine-number call with a **real order** (machine or
  exact) is answered by the raw-double kernels in
  `src/special_functions/sf_machine_bessel.c` — the same helper the
  NDArray/packed element-wise path and `Compile[]` call, so the scalar, List,
  packed, `NDArray` and compiled answers agree bit for bit. A real argument
  gives a `Real` where the value is real (and a `Complex` on the cut, e.g.
  `BesselJ[1/2, -2.]`); a complex argument gives a `Complex`. Regions:
  - real axis: integer-order J/Y from libm `jn`/`yn`; I₀, I₁, K₀, K₁ from
    Cephes-style regions (series / piecewise Chebyshev fits of the
    exponentially scaled function / Hankel), Kₙ by upward recurrence;
  - everything else by Amos's scheme (ACM TOMS 644): J, Y through I and K of
    u = −iz in the right half plane, with K from Temme's series (|u| ≤ 1.25),
    Taylor continuation from tabulated nodes (1.25 ≤ |u| < 21, orders 0 and
    ±1/3 mod 1; Steed's CF2 for other fractional orders) or Hankel's expansion;
    I from Miller's backward recurrence normalised by a Gegenbauer/Neumann sum.
  Accuracy against mpmath (`tools/sf_accuracy/bessel.py`, ulps of the local
  magnitude): ≤ 3.6 on the real axis, median ≈ 1 and max ≤ 8 for complex
  arguments on the benchmark domain (≤ 12 in a wide band), and median ≈ 1.5,
  max ≤ ~20 for non-integer orders. A kernel
  decline (|ν| > 1000, overflow/underflow, a complex order) falls through to
  the MPFR path below, unchanged.
- **Numerics.** When some argument is inexact, machine- and arbitrary-precision
  (MPFR) real and complex order and argument evaluate numerically. The core
  routes on `r = |z|` and the requested precision `P`:
  - an **asymptotic series** for large `r`, valid on the whole principal sheet
    (DLMF 10.40.2, `K_ν(z) ~ sqrt(π/(2z)) e⁻ᶻ Σ a_k/zᵏ`, summed to the optimal
    smallest-term truncation);
  - for small/moderate `r` and **non-integer order**, the connection formula
    `K_ν = (π/2)(I₋ᵥ − Iᵥ)/sin(νπ)` (DLMF 10.27.4) with the modified-Bessel
    series `Iᵤ`, plus `−log₂|sin(νπ)|` guard bits to absorb the cancellation as
    ν → integer;
  - for small/moderate `r` and **integer order**, the logarithmic series
    (DLMF 10.31.1), with `ψ(k+1) = −γ + H_k` accumulated incrementally.

  Examples: `BesselK[0, 0.53] = 0.87656`, `BesselK[0, 4.0] = 0.0111597`,
  `BesselK[1 + I, 3.0 - 2 I] = -0.0225108 + 0.0169607 I`,
  `N[BesselK[0, 4], 50] = 0.011159676085853024269745195979833489225009023888474`.
  Precision tracks the input.
- **Half-integer order.** With a **symbolic** argument, `BesselK[(2k+1)/2, x]`
  rewrites to the elementary closed form (e.g.
  `BesselK[1/2, x] = Sqrt[Pi/(2 x)] Exp[-x]`,
  `BesselK[3/2, x] = Sqrt[Pi/(2 x)] Exp[-x] (1 + 1/x)`); with a numeric argument
  the accurate general numeric path is used instead, and an exact numeric
  argument such as `BesselK[11/2, 1]` is left unevaluated. These rewrites and the
  even-order reflection `BesselK[-n, x] = BesselK[n, x]` live in the C builtin
  (`bessel_symbolic_rules`; formerly `src/internal/bessel.m`).
- **Derivative.** `D[BesselK[n, x], x] = -(BesselK[n-1, x] + BesselK[n+1, x])/2`
  (DLMF 10.29.5; note the sign differs from `BesselJ`).
- **Series at 0.** A logarithmic expansion (the coefficients carry `Log[x]`):
  `Series[BesselK[0, x], {x, 0, 3}] =
  (-EulerGamma + Log[2] - Log[x]) + (1 - EulerGamma + Log[2] - Log[x]) x²/4 +
  O[x]⁴`. Integer orders ≥ 1 also carry a 1/xⁿ principal part.
- **Series at ∞.** `Series[BesselK[n, x], {x, Infinity, k}]` gives the
  asymptotic expansion `E⁻ˣ (Sqrt[π/2]/Sqrt[x] − Sqrt[π/2]/(8 x^(3/2)) + …)`,
  valid for symbolic order too. `Limit[BesselK[n, x], x -> Infinity] = 0`.
- **Listable.** `BesselK[{1, 2, 3}, 1.0] = {0.601907, 1.62484, 7.10126}`.

Attributes: `Listable`, `NumericFunction`, `Protected`.

On the branch cut (negative real axis) the principal-branch value is returned
consistently at all precisions; this convention can differ from other systems'
large-order machine-precision heuristics.

```mathematica
In[1]:= BesselK[0, 4.0]
Out[1]= 0.0111597

In[2]:= D[BesselK[n, x], x]
Out[2]= -1/2 (BesselK[-1 + n, x] + BesselK[1 + n, x])

In[3]:= BesselK[0, -1.0]
Out[3]= 0.421024 - 3.97746*I
```

## BesselI

- `BesselI[n, z]` — the modified Bessel function of the first kind Iₙ(z), the
  solution of z² y″ + z y′ − (z² + n²) y = 0 that is regular at the origin. It
  grows like eᶻ as z → ∞, is even in the order for integer index (I₋ₙ = Iₙ), and
  has a branch cut along the negative real z axis for non-integer order.

Implemented in `src/special_functions/bessel.c` (the consolidated Bessel module,
alongside `BesselJ`/`BesselK`), reusing the `ncpx` complex-MPFR toolkit
(`src/numeric_complex.{c,h}`). As with `BesselK` there is **no** MPFR library
routine for the modified Bessel functions, so every value is summed from scratch;
the small-|z| power series is the same `I_μ(z)` series `BesselK` already uses.

- **Exact values at the origin.** Same leading `(z/2)^ν/Γ(ν+1)` law as `BesselJ`:
  `BesselI[0, 0] = 1`; `0` for integer n ≠ 0 or non-integer Re(ν) > 0 (e.g.
  `BesselI[1/2, 0] = 0`); `ComplexInfinity` for non-integer Re(ν) < 0 (e.g.
  `BesselI[-1/2, 0] = ComplexInfinity`); `Indeterminate` for pure-imaginary order;
  symbolic order at the origin stays unevaluated. Other exact / symbolic arguments
  stay unevaluated and numericalize only under `N`.
- **Machine precision.** A machine-number call with a **real order** (machine or
  exact) is answered by the raw-double kernels in
  `src/special_functions/sf_machine_bessel.c` — the same helper the
  NDArray/packed element-wise path and `Compile[]` call, so the scalar, List,
  packed, `NDArray` and compiled answers agree bit for bit. A real argument
  gives a `Real` where the value is real (and a `Complex` on the cut, e.g.
  `BesselJ[1/2, -2.]`); a complex argument gives a `Complex`. Regions:
  - real axis: integer-order J/Y from libm `jn`/`yn`; I₀, I₁, K₀, K₁ from
    Cephes-style regions (series / piecewise Chebyshev fits of the
    exponentially scaled function / Hankel), Kₙ by upward recurrence;
  - everything else by Amos's scheme (ACM TOMS 644): J, Y through I and K of
    u = −iz in the right half plane, with K from Temme's series (|u| ≤ 1.25),
    Taylor continuation from tabulated nodes (1.25 ≤ |u| < 21, orders 0 and
    ±1/3 mod 1; Steed's CF2 for other fractional orders) or Hankel's expansion;
    I from Miller's backward recurrence normalised by a Gegenbauer/Neumann sum.
  Accuracy against mpmath (`tools/sf_accuracy/bessel.py`, ulps of the local
  magnitude): ≤ 3.6 on the real axis, median ≈ 1 and max ≤ 8 for complex
  arguments on the benchmark domain (≤ 12 in a wide band), and median ≈ 1.5,
  max ≤ ~20 for non-integer orders. A kernel
  decline (|ν| > 1000, overflow/underflow, a complex order) falls through to
  the MPFR path below, unchanged.
- **Numerics.** When some argument is inexact, machine- and arbitrary-precision
  (MPFR) real and complex order and argument evaluate numerically. The core
  routes on `r = |z|` and the requested precision `P`:
  - a **power series** for small/moderate `r` (DLMF 10.25.2,
    `t_k = t_{k-1}·(z/2)²/(k(ν+k))`, one `Gamma(ν+1)`; for a non-negative
    integer order `(z/2)ⁿ` is formed by repeated multiplication so a real `z`
    stays exactly real, with ~`2r/ln2` guard bits for complex-argument
    cancellation);
  - an **asymptotic series** for large `r`, valid on the whole principal sheet
    (DLMF 10.40.5, the two-exponential form `eᶻ Σ(−1)ᵏaₖ/zᵏ +
    e^{i(ν+½)π}e⁻ᶻ Σaₖ/zᵏ`, summed to the optimal smallest-term truncation; the
    e⁻ᶻ term dominates on the left half-plane). When the order is real, `z` is
    real, and the result is mathematically real (integer order, or `z ≥ 0`), the
    sub-precision imaginary noise the complex phase would leave is dropped.

  Examples: `BesselI[0, 2.0] = 2.27959`, `BesselI[3 + I, 1.5 - I] =
  −0.25665 + 0.0492771 I`, `N[BesselI[0, 1], 50] =
  1.2660658777520083355982446252147175376076703113550`. Precision tracks the
  input.
- **Half-integer order.** With a **symbolic** argument, `BesselI[(2k+1)/2, x]`
  rewrites to the elementary closed form — unlike `BesselK`, I is **not** even
  for half-integer order, so both base cases are needed
  (`BesselI[1/2, x] = Sqrt[2/(Pi x)] Sinh[x]`,
  `BesselI[-1/2, x] = Sqrt[2/(Pi x)] Cosh[x]`) with an up/down recurrence; with a
  numeric argument the accurate general numeric path is used instead, and an
  exact numeric argument such as `BesselI[11/2, 1]` is left unevaluated. These
  rewrites and the even-order reflection `BesselI[-n, x] = BesselI[n, x]` live in
  the C builtin (`bessel_symbolic_rules`; formerly `src/internal/bessel.m`).
- **Argument parity.** For **integer** order, Iₙ is entire and obeys
  `BesselI[n, -z] = (-1)ⁿ BesselI[n, z]` (e.g. `BesselI[0, -z] = BesselI[0, z]`,
  `BesselI[1, -z] = -BesselI[1, z]`), folded in the C builtin for symbolic /
  exact arguments exactly as for `BesselJ`. `BesselK` carries **no** such fold —
  its negative-real axis is a branch cut, so `BesselK[n, -z]` stays unevaluated.
- **Derivative.** `D[BesselI[n, x], x] = (BesselI[n-1, x] + BesselI[n+1, x])/2`
  (DLMF 10.29.5; the two-term sum carries a `+` like `BesselK`, but the overall
  coefficient is `+1/2`).
- **Series at 0.** `Series[BesselI[0, x], {x, 0, 10}] =
  1 + x²/4 + x⁴/64 + x⁶/2304 + x⁸/147456 + x¹⁰/14745600 + O[x]¹¹` (via the
  generic Taylor path; Puiseux for half-integer order).
- **Series at ∞.** `Series[BesselI[n, x], {x, Infinity, k}]` gives the full
  two-exponential asymptotic (DLMF 10.40.5), e.g. for n = 0
  `E^x (1/(Sqrt[2π] Sqrt[x]) + 1/(8 Sqrt[2π] x^(3/2)) + …) +
  I E^(-x) (1/(Sqrt[2π] Sqrt[x]) − 1/(8 Sqrt[2π] x^(3/2)) + …)`, matching
  the reference's `Normal` form, valid for symbolic order too.
- **Indefinite integrals** (recurrence forms, `src/internal/CRCMathTablesIntegrals.m`):
  `∫ xᵖ⁺¹ Iₚ(x) dx = xᵖ⁺¹ Iₚ₊₁(x)`, `∫ x⁻ᵖ⁺¹ Iₚ(x) dx = x⁻ᵖ⁺¹ Iₚ₋₁(x)`,
  `∫ I₁(x) dx = I₀(x)`.
- **Listable.** `BesselI[{0, 1, 2}, 1.] = {1.26607, 0.565159, 0.135748}`.

Attributes: `Listable`, `NumericFunction`, `Protected`.

Not yet implemented: `Limit[BesselI[n, x], x -> Infinity]` (= ∞) and at
`I Infinity` (= 0) are left unevaluated — the Limit engine does not distribute
over the `Plus` of the two-exponential asymptotic at infinity (the `Series`
itself is exact).

```mathematica
In[1]:= BesselI[0, 2.0]
Out[1]= 2.27959

In[2]:= D[BesselI[n, x], x]
Out[2]= 1/2 (BesselI[-1 + n, x] + BesselI[1 + n, x])

In[3]:= BesselI[1, -z]
Out[3]= -BesselI[1, z]
```

## BesselY

- `BesselY[n, z]` — the Bessel function of the second kind Yₙ(z), the solution of
  `z² y″ + z y′ + (z² − n²) y = 0` that is singular (logarithmically, for integer
  order) at the origin. Second linearly-independent companion of `BesselJ`.

Implemented in `src/special_functions/bessel.c` (the consolidated Bessel module,
alongside `BesselJ`/`BesselK`/`BesselI`), reusing the `ncpx` complex-MPFR toolkit
(`src/numeric_complex.{c,h}`) and the `BesselJ` power-series kernel directly.

- **Exact values at the origin.** `BesselY[0, 0] = -Infinity`. Yₙ diverges for
  every other order **except the negative half-odd-integers**, where `cos(νπ)`
  cancels the divergent term (`Y_{-1/2}(z) = J_{1/2}(z) → 0`,
  `Y_{-3/2} = -J_{3/2} → 0`, …) — so `BesselY[-1/2, 0] = 0`,
  `BesselY[-3/2, 0] = 0`, etc. Everything else gives `ComplexInfinity` (integer
  n ≠ 0, non-integer Re(ν) > 0, negative non-half-odd-integer, and symbolic
  order), with `Indeterminate` for pure-imaginary order. Other exact arguments
  stay symbolic.
- **Machine precision.** A machine-number call with a **real order** (machine or
  exact) is answered by the raw-double kernels in
  `src/special_functions/sf_machine_bessel.c` — the same helper the
  NDArray/packed element-wise path and `Compile[]` call, so the scalar, List,
  packed, `NDArray` and compiled answers agree bit for bit. A real argument
  gives a `Real` where the value is real (and a `Complex` on the cut, e.g.
  `BesselJ[1/2, -2.]`); a complex argument gives a `Complex`. Regions:
  - real axis: integer-order J/Y from libm `jn`/`yn`; I₀, I₁, K₀, K₁ from
    Cephes-style regions (series / piecewise Chebyshev fits of the
    exponentially scaled function / Hankel), Kₙ by upward recurrence;
  - everything else by Amos's scheme (ACM TOMS 644): J, Y through I and K of
    u = −iz in the right half plane, with K from Temme's series (|u| ≤ 1.25),
    Taylor continuation from tabulated nodes (1.25 ≤ |u| < 21, orders 0 and
    ±1/3 mod 1; Steed's CF2 for other fractional orders) or Hankel's expansion;
    I from Miller's backward recurrence normalised by a Gegenbauer/Neumann sum.
  Accuracy against mpmath (`tools/sf_accuracy/bessel.py`, ulps of the local
  magnitude): ≤ 3.6 on the real axis, median ≈ 1 and max ≤ 8 for complex
  arguments on the benchmark domain (≤ 12 in a wide band), and median ≈ 1.5,
  max ≤ ~20 for non-integer orders. A kernel
  decline (|ν| > 1000, overflow/underflow, a complex order) falls through to
  the MPFR path below, unchanged.
- **Numeric evaluation** at machine and arbitrary precision, real and complex,
  routed by `|z|` and order:
  - integer order, real **z > 0** → MPFR-native `mpfr_yn` (correctly rounded);
  - large `|z|` → trig-prefactored asymptotic series (DLMF 10.17.4);
  - small `|z|`, non-integer order → connection formula
    `Y_ν = (J_ν cos(νπ) − J_{−ν}) / sin(νπ)` (DLMF 10.2.3);
  - small `|z|`, integer order → logarithmic series (DLMF 10.8.1).
  Precision tracks the least-precise inexact argument. Examples:
  `BesselY[0, 2.5] = 0.49807`, `BesselY[0.5 I, 3 − I] = 1.04686 + 0.884784 I`,
  `N[BesselY[0, 1], 50] = 0.088256964215676957982926766023515162827817523090676`.
- **Branch cut.** Yₙ(z) has a logarithmic branch point at 0 and a branch cut
  along the negative real axis, so `BesselY[0, -1.0] = 0.0882570 + 1.530395 I` is
  genuinely complex. Unlike `BesselJ`/`BesselI` there is **no** argument-parity
  fold — `BesselY[n, -z]` stays unevaluated (as for `BesselK`).
- **Half-integer order.** With a **symbolic** argument, `BesselY[(2k+1)/2, x]`
  rewrites to the elementary spherical-Bessel closed form
  (`BesselY[1/2, x] = -Sqrt[2/(Pi x)] Cos[x]`,
  `BesselY[-1/2, x] = Sqrt[2/(Pi x)] Sin[x]`) with an up/down recurrence; with a
  numeric argument the call uses the accurate C path, and an exact numeric
  argument such as `BesselY[11/2, 1]` is left unevaluated. These rewrites and the
  negative-integer reflection `BesselY[-n, x] = (-1)ⁿ BesselY[n, x]` live in
  the C builtin (`bessel_symbolic_rules`; formerly `src/internal/bessel.m`).
- **Derivative.** `D[BesselY[n, x], x] = (BesselY[n-1, x] - BesselY[n+1, x])/2`
  (DLMF 10.6.1; same shape as `BesselJ`).
- **Series at 0.** Integer order gives the logarithmic expansion, e.g.
  `Series[BesselY[0, x], {x, 0, 3}] =
  (2 (EulerGamma - Log[2] + Log[x]))/Pi +
  ((1 - EulerGamma + Log[2] - Log[x]) x²)/(2 Pi) + O[x]^4` (DLMF 10.8.1).
- **Series at ∞.** `Series[BesselY[n, x], {x, Infinity, k}]` gives the
  trig-prefactored asymptotic expansion (DLMF 10.17.4), valid for symbolic order
  too.
- **Listable.** `BesselY[0, {1.0, 2.0, 3.0}] = {0.088257, 0.510376, 0.37685}`.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= BesselY[0, 2.5]
Out[1]= 0.49807

In[2]:= D[BesselY[n, x], x]
Out[2]= 1/2 (BesselY[-1 + n, x] - BesselY[1 + n, x])

In[3]:= BesselY[0, -1.0]
Out[3]= 0.088257 + 1.5304*I
```

## Hyperfactorial

`Hyperfactorial[n]` gives the hyperfactorial `prod_{k=1}^{n} k^k`
(`H(0) = H(1) = 1`).  Exact via GMP for a non-negative integer `n`.  A
non-integer **numeric** order (under `N`) evaluates via
`Hyperfactorial[z] = Gamma[z+1]^z / BarnesG[z+1]`, reusing the Barnes-G
continuation (real, complex, arbitrary precision); symbolic orders are left
unevaluated.  Used by `Product` to recognise `Product[k^k, {k, 1, n}]`.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= Hyperfactorial[4]
Out[1]= 27648

In[2]:= N[Hyperfactorial[5.0]]
Out[2]= 8.64*10^7

In[3]:= N[Hyperfactorial[7/2], 30]
Out[3]= 1282.122099453457459415422713168
```

## BarnesG

`BarnesG[z]`, the Barnes G-function, with `G(1) = G(2) = 1` and
`G(z+1) = Gamma[z] G(z)`.  For a positive integer `n`,
`G(n+1) = prod_{k=1}^{n-1} k!` (the superfactorial, exact via GMP), and
`G(m) = 0` for non-positive integer `m`.  Used by `Product` to recognise
`Product[Gamma[i], {i, 1, n-1}]` → `BarnesG[n]`.

Exact non-integer orders are left symbolic, but a non-integer **numeric**
order (i.e. under `N`) is evaluated from the Barnes asymptotic expansion
combined with the recurrence `G(w+1) = Gamma[w] G(w)`, assembled as an
expression so it inherits `Gamma`/`Log`/`Exp`'s complex and MPFR kernels and
the Glaisher–Kinkelin constant at full precision. The series carries 12
Bernoulli terms with coefficient `B_{2k+2}/((2k)(2k+2))`, and the recurrence
shift scales with the requested precision. Real, complex, and arbitrary
precision are all supported.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= BarnesG[5]
Out[1]= 12

In[2]:= N[BarnesG[6.0]]
Out[2]= 288.

In[3]:= N[BarnesG[13/2], 30]
Out[3]= 2548.745769568498989735906104648

In[4]:= N[BarnesG[2.5 + 1.0 I]]
Out[4]= 0.743798 - 0.0953168 I
```

## QPochhammer

`QPochhammer[a, q, n]` gives the q-Pochhammer symbol (q-shifted
factorial) `prod_{k=0}^{n-1} (1 - a q^k)`; `QPochhammer[a, q]` gives the
infinite product `(a;q)_Inf` for `|q| < 1`.  The finite form is
exact/symbolic for a non-negative integer `n` (a symbolic `n` is left
unevaluated, which is what `Product` emits as a closed form); the
infinite form evaluates for machine-real `a`, `q`.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= QPochhammer[a, q, 3]
Out[1]= (1 - a) (1 - a q) (1 - a q^2)
```

- **v0.477**: the machine-real `(a; q)_∞` uses the shared kernel
  `sf_machine_qpochhammer` (also the NDArray / `Compile` one); a complex or
  arbitrary-precision `a, q` evaluates through arb at the call's precision (it
  used to be rounded to a double). A non-integer `n` with an inexact argument is
  `(a; q)_∞ / (a qⁿ; q)_∞` (`QPochhammer[0.5, 0.3, 2.5]`, complex `n` too), and a
  negative integer `n` is `1/Π_{k=1}^{|n|} (1 - a q^-k)` (`QPochhammer[1/2, 1/3,
  -2] = 4/7`). A vanishing factor gives `0.`, not `-0.`.

## ProductLog

- `ProductLog[z]` — the principal branch W₀(z) of the Lambert W function, the
  solution `w` of `z = w e^w` that is real for `z ≥ -1/e`.
- `ProductLog[k, z]` — the k-th branch Wₖ(z) (`k` any integer, `k = 0`
  principal). Branches are ordered by imaginary part.

Implemented in `src/special_functions/productlog.c`. Numeric evaluation runs
through a single complex-MPFR Halley core built on the shared `ncpx` toolkit
(`numeric_complex.h`).

- **Exact values.** `ProductLog[0] = 0`, `ProductLog[E] = 1`,
  `ProductLog[-1/E] = -1`, `ProductLog[-Pi/2] = I Pi/2`,
  `ProductLog[Infinity] = ProductLog[ComplexInfinity] = Infinity`, and
  `ProductLog[k, 0] = -Infinity` for `k ≠ 0`. Other exact arguments (e.g.
  `ProductLog[1/3]`) stay unevaluated and numericalize only under `N`.
- **Numerics.** With an inexact argument, machine- and arbitrary-precision
  (MPFR) real and complex inputs evaluate numerically on any branch. The core
  seeds an initial approximation by region — a branch-point series in
  `p = Sqrt[2(e z + 1)]` near `z = -1/e` (branches 0 and -1), the Maclaurin seed
  `z(1 - z + 3/2 z²)` for the principal branch near 0, and otherwise the
  asymptotic `L1 - L2 + L2/L1` with `L1 = Log z + 2πi k`, `L2 = Log L1` — then
  refines with Halley's cubically-convergent iteration. A real seed keeps the
  iteration exactly real, so real-valued branches return a real leaf without
  imaginary noise; otherwise the result is `Complex[…]`. Examples:
  `ProductLog[2.5] = 0.958586`, `ProductLog[-1.5] = -0.0327837 + 1.54964 I`,
  `ProductLog[1 + 3.5 I] = 1.0546 + 0.703928 I`,
  `Table[ProductLog[k, 2.3], {k, -2, 2}]` spans the five lowest branches.
  Precision tracks the input: `N[ProductLog[1/3], 100]` and `ProductLog[7/3\`100]`
  are accurate to 100 digits.
- **The machine kernel.** A machine-precision argument (`ProductLog[x]`,
  `ProductLog[0, x]` with a machine real or complex `x`) never reaches the MPFR
  core: the scalar builtin, the packed/`NDArray` path and `Compile[]` (real and
  complex) all call the same double kernels in
  `src/special_functions/sf_machine_plog.c`, so they agree bit-for-bit
  (~4 ulp max, ~0.3 ulp median against mpmath). `sf_machine_productlog` (real
  `x ≥ -1/e`) seeds Halley by region — `x(1 - x + 3/2 x²)` below 1, `Log[1 + x]`
  on `[1, e]`, `L1 - L2 + L2/L1` above `e` (only where `Log x > 1`; `L2/L1`
  diverges as `L1 → 0`) — and **verifies** `w + Log w = Log x` before answering.
  `sf_machine_productlog_c` (complex, principal branch) seeds with the Padé
  `z(3 + 6z + z²)/(3 + 9z + 5z²)` near the origin or the asymptotic form
  elsewhere, iterates Halley in re/im doubles, then checks `w + Log w = Log z`
  (mod 2πi) and that `w` lies in W₀'s range (`|Im w| < π`, right of the curve
  `-v cot v + i v`). **Near the branch point** (`|e z + 1| < 1/2`) both kernels
  run Halley in `u = w + 1` on `g(u) - e(z + 1/e)` with
  `g(u) = (u - 1)e^u + 1` summed from its Maclaurin series, and form
  `z + 1/e` with `1/e` split into a double head and tail: the textbook residual
  `w e^w - z` is the difference of two numbers ≈ `-1/e` and is pure rounding
  there. The split also decides the real/complex boundary correctly — the
  double nearest `-1/e` lies *below* it, so `ProductLog[-0.36787944117144233]`
  is `-1. + 8.22×10⁻⁹ I`, as the MPFR core says. A real `x < -1/e` returns the
  complex principal value (the upper side of the cut); a branch `k ≠ 0`, a
  failed verification or a non-finite value falls through to the MPFR core.
  Exact special values (`-1/E`, `-Pi/2`) are matched only for a `Times[...]`
  argument, so a numeric call builds no candidate expressions.
- **Derivative.** `D[ProductLog[z], z] = ProductLog[z]/(z (1 + ProductLog[z]))`;
  the two-argument form differentiates the same way in the argument `z`.
- **Series at 0.** `Series[ProductLog[x], {x, 0, n}]` gives the closed-form
  Taylor series `x - x^2 + 3/2 x^3 - 8/3 x^4 + 125/24 x^5 + …`
  (coefficient `(-k)^(k-1)/k!`).
- **Series at the branch point -1/E.** A Puiseux series in `Sqrt[x + 1/E]`:
  `-1 + Sqrt[2 E] Sqrt[x + 1/E] - 2/3 E (x + 1/E) + …`.
- **Series at ∞.** The nested-logarithm asymptotic expansion (the x⁰
  coefficient): `Log[x] - Log[Log[x]] + Log[Log[x]]/Log[x] - … + O[1/x]`.
- **SeriesCoefficient.** Numeric indices reduce via the series; the symbolic
  general term is `SeriesCoefficient[ProductLog[x], {x, 0, n}] =
  Piecewise[{{(-n)^(n-1)/n!, n ≥ 1}}, 0]`.
- **Listable.** `ProductLog[{1.5, 3.75, 5.5, 7.25}] = {0.725861, 1.16717,
  1.38155, 1.54559}`.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= ProductLog[1.0]
Out[1]= 0.567143

In[2]:= ProductLog[-1/E]
Out[2]= -1

In[3]:= ProductLog[-1, -0.2]
Out[3]= -2.54264

In[4]:= ProductLog[-0.5]
Out[4]= -0.794024 + 0.770112*I

In[5]:= Series[ProductLog[x], {x, 0, 4}]
Out[5]= x - x^2 + 3/2 x^3 - 8/3 x^4 + O[x]^5
```

- **Inexact branch index** (v0.6490). `N` maps over both arguments, so
  `N[ProductLog[-1, 2]]` arrives as `ProductLog[-1., 2.]`; an integer-valued machine
  or MPFR index names the same branch and evaluates (a fractional index does not).

## ChebyshevT

`ChebyshevT[n, x]` gives the Chebyshev polynomial of the first kind `T_n(x)`.

- **Exact polynomial.** An exact integer order `n` produces the explicit
  degree-`|n|` polynomial in `x` with exact integer coefficients (from FLINT's
  `fmpz_poly_chebyshev_t` when available, otherwise the three-term recurrence
  `T_0 = 1`, `T_1 = x`, `T_k = 2x T_{k-1} - T_{k-2}`). `T` is even in its order
  index, so a negative order reflects via `T_{-n} = T_n` (no sign, and
  `T_{-1} = x`). For example `ChebyshevT[10, x] = -1 + 50 x² - 400 x⁴ + 1120 x⁶
  - 1280 x⁸ + 512 x¹⁰`, and `ChebyshevT[2, 3] = 17`.
- **Numerics.** For an integer order with an inexact argument the value comes
  from the *stable* three-term recurrence, not from substituting into the
  explicit polynomial — so `ChebyshevT[60, N[2/3, 16]] = 0.980336` rather than a
  catastrophically-cancelled result. A non-integer or complex order with an
  inexact argument evaluates through `T_v(x) = 2F1(-v, v; 1/2; (1-x)/2)` where
  the Gauss series converges (`|(1-x)/2| < 1`), and otherwise through the
  whole-plane closed form `T_v(x) = Cos[v ArcCos x]`, which takes the reference's
  branch cut from `-∞` to `-1`. All paths run at machine or arbitrary (MPFR)
  precision, real or complex: `ChebyshevT[12, 2.5] = 7.31808× 10⁷`,
  `N[ChebyshevT[1/7, 1/3], 50]` is accurate to 50 digits, `ChebyshevT[2 - 4 I,
  3`40] = 12.23044… - 11.78701… I`, and `ChebyshevT[-1.5, -2.4`100] ≈ -4.85263
  I`.
- **Special values.** `ChebyshevT[0, x] = 1`, `ChebyshevT[n, 0] = Cos[n π/2]`,
  and `ChebyshevT[n, 1] = 1` (the last two hold for symbolic order).
- **Derivative.** `D[ChebyshevT[n, x], x] = n U_{n-1}(x)` (the first kind
  differentiates to the second kind), which also lets the generic
  `NumericFunction` interval threading handle `Interval[]` arguments.
- **Listable.** `ChebyshevT[2, {1, 2, 3}] = {1, 7, 17}`.
- **Performance (machine reals and complexes).** For an integer order with a
  machine-precision argument the recurrence runs in hardware `double` /
  `double`-complex, the same kernel used by the NDArray fast path and
  `Compile[]`, so a scalar `ChebyshevT[n, x]` is bit-identical to its packed and
  compiled forms. A packed **complex** buffer stays on the buffer (the ND kernel
  has a complex arm). An arbitrary-precision argument, a complex argument of
  degree above 1000, or an overflow falls back to the stable MPFR recurrence
  unchanged, so accuracy and the overflow→MPFR behaviour are preserved.
- **High degree (real argument).** Above degree 1000 (machine) and degree 64
  (arbitrary precision) a real `x` takes the closed form `Cos[n ArcCos[x]]` on
  `[-1, 1]` and `Sign[x]^n Cosh[n ArcCosh[Abs[x]]]` off it: `O(1)` operations
  instead of `O(n)`, and free of the recurrence's `~n²/2`-ulp error growth near
  `x = ±1`. The closed form's error is within a factor `2 θ sin θ ≤ 3.7` of the
  error a rounding of `x` already causes (`|T_n'(x)| ε`), so it is as accurate as
  the argument allows. `ChebyshevT[2001, 0.9999997]` takes ~1 µs (was ~400 µs
  through MPFR), the same value on every surface; the machine helper still uses
  the recurrence up to degree 1000, where it is exact on small integers.

An exact non-integer order with an exact/symbolic argument (e.g.
`ChebyshevT[1/2, x]`) stays symbolic; symbolic `Series` / `SeriesCoefficient`
and the half-integer closed forms are deferred.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= ChebyshevT[7, x]
Out[1]= -7 x + 56 x^3 - 112 x^5 + 64 x^7

In[2]:= Table[ChebyshevT[n, Pi/2], {n, 0, 3}]
Out[2]= {1, 1/2 Pi, -1 + 1/2 Pi^2, -3/2 Pi + 1/2 Pi^3}

In[3]:= {ChebyshevT[60, 0.7], Cos[60 ArcCos[0.7]]}
Out[3]= {-0.825313, -0.825313}
```

## ChebyshevU

`ChebyshevU[n, x]` gives the Chebyshev polynomial of the second kind `U_n(x)`.

- **Exact polynomial.** An exact integer order `n` produces the explicit
  degree-`|n'|` polynomial in `x` with exact integer coefficients (from FLINT's
  `fmpz_poly_chebyshev_u` when available, otherwise the three-term recurrence
  `U_0 = 1`, `U_1 = 2x`, `U_k = 2x U_{k-1} - U_{k-2}`). Negative orders reflect
  via `U_{-1} = 0` and `U_{-m} = -U_{m-2}`. For example `ChebyshevU[10, x] = -1
  + 60 x² - 560 x⁴ + 1792 x⁶ - 2304 x⁸ + 1024 x¹⁰`, and `ChebyshevU[2, 3] = 35`.
- **Numerics.** For an integer order with an inexact argument the value comes
  from the *stable* three-term recurrence, not from substituting into the
  explicit polynomial — so `ChebyshevU[60, N[2/3, 16]] = 1.15684` rather than the
  catastrophically-cancelled `0.`. A non-integer or complex order with an
  inexact argument evaluates through `U_v(x) = (v+1) 2F1(-v, v+2; 3/2; (1-x)/2)`
  where the Gauss series converges (`|(1-x)/2| < 1`), and otherwise through the
  whole-plane closed form `U_v(x) = Sin[(v+1) ArcCos x]/Sin[ArcCos x]`, which
  takes the reference's branch cut from `-∞` to `-1`. All paths run at machine or
  arbitrary (MPFR) precision, real or complex: `ChebyshevU[12, 2.5] = 1.53028×
  10⁸`, `N[ChebyshevU[1/7, 1/3], 50]` is accurate to 50 digits, `ChebyshevU[2 -
  4 I, 3`40] = 25.18033… - 24.31072… I`, and `ChebyshevU[-1.5, -2.4`100] ≈
  0.597614 I`.
- **Special values.** `ChebyshevU[0, x] = 1`, `ChebyshevU[n, 0] = Cos[n π/2]`,
  and `ChebyshevU[n, 1] = n + 1` (the last two hold for symbolic order).
- **Derivative.** `D[ChebyshevU[n, x], x] = ((-1-n) U_{n-1}(x) + n x U_n(x)) /
  (x² - 1)`, which also lets the generic `NumericFunction` interval threading
  handle `Interval[]` arguments.
- **Listable.** `ChebyshevU[2, {1, 2, 3}] = {3, 15, 35}`.
- **Performance (machine reals and complexes).** For an integer order with a
  machine-precision argument the recurrence runs in hardware `double` /
  `double`-complex, the same kernel used by the NDArray fast path and
  `Compile[]`, so a scalar `ChebyshevU[n, x]` is bit-identical to its packed and
  compiled forms. A packed **complex** buffer stays on the buffer (the ND kernel
  has a complex arm — `ChebyshevU[10, cvec]` over `10^6` complexes is ~1200×
  faster than materialising to the scalar path, and faster than SciPy's
  `eval_chebyu`). An arbitrary-precision argument, a degree above 1000, or an
  overflow falls back to the stable MPFR recurrence unchanged, so accuracy and
  the overflow→MPFR behaviour are preserved.

An exact non-integer order with an exact/symbolic argument (e.g.
`ChebyshevU[1/2, x]`) stays symbolic; symbolic `Series` / `SeriesCoefficient`
and the half-integer closed forms are deferred.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= ChebyshevU[7, x]
Out[1]= -8 x + 80 x^3 - 192 x^5 + 128 x^7

In[2]:= Table[ChebyshevU[n, Pi/2], {n, 0, 3}]
Out[2]= {1, Pi, -1 + Pi^2, -2 Pi + Pi^3}

In[3]:= ChebyshevU[3, 2. + I]
Out[3]= 8.0 + 84.0*I
```

## LaguerreL

`LaguerreL[n, x]` gives the Laguerre polynomial `L_n(x)`; `LaguerreL[n, a, x]`
gives the generalized Laguerre polynomial `L_n^(a)(x)`, with `LaguerreL[n, x] ==
LaguerreL[n, 0, x]`. They are orthogonal with weight `x^a e^-x` and satisfy
`x y'' + (a + 1 - x) y' + n y = 0`.

- **Exact polynomial.** A non-negative integer order `n` produces the explicit
  degree-`n` polynomial in `x`. For a numeric (or absent) parameter `a` the
  coefficients are exact rationals, built from the three-term recurrence `L_0 =
  1`, `L_1 = 1 + a - x`, `(k+1) L_{k+1} = (2k+1+a-x) L_k - (k+a) L_{k-1}`; for a
  symbolic `a` the polynomial in `a` and `x` is built from the sum form and
  expanded. For example `LaguerreL[5, x] = 1 - 5 x + 5 x² - 5/3 x³ + 5/24 x⁴ -
  1/120 x⁵`, `LaguerreL[2, a, x] = 1 + 3/2 a + 1/2 a² - 2 x - a x + 1/2 x²`, and
  `LaguerreL[2, 3] = -1/2`.
- **Numerics.** For an integer order with an inexact argument the value comes
  from the *stable* recurrence, not from substituting into the explicit
  polynomial — so `LaguerreL[100, 0, N[5, 20]] = 1.4555271625328812607` rather
  than the catastrophically-cancelled `1.4557…` that the expanded polynomial
  gives. A non-integer or complex order with an inexact argument evaluates
  through the entire confluent-hypergeometric form `LaguerreL[n, a, x] =
  Binomial[n+a, n] Hypergeometric1F1[-n, a+1, x]` (a single expression that
  covers the whole plane, since `L_n` is entire in `x`). All paths run at
  machine or arbitrary (MPFR) precision, real or complex: `LaguerreL[6., 5] =
  -2.09028`, `LaguerreL[6. - I, 5] = -2.76039 - 2.46853 I`, `LaguerreL[3, .51,
  .87] = -0.396576`, `LaguerreL[2.5, 1.3, 8 + I] = 8.50391 + 0.576867 I`, and
  `N[LaguerreL[2/5, 1/3, 8/7], 50]` is accurate to 50 digits.
- **Special values.** `LaguerreL[0, x] = 1` and `LaguerreL[n, 0] = Binomial[n+a,
  n]` (so `LaguerreL[n, 0] = 1` for the ordinary polynomial).
- **Derivative.** `D[LaguerreL[n, x], x] = -LaguerreL[n-1, 1, x]` and
  `D[LaguerreL[n, a, x], x] = -LaguerreL[n-1, a+1, x]`, which also lets the
  generic `NumericFunction` interval threading handle `Interval[]` arguments.
- **Listable.** `LaguerreL[3, {1, 2, 3}] = {-2/3, -1/3, 1}`.
- **Performance (machine reals and complexes).** For the two-argument
  `LaguerreL[n, x]` at an integer order with a machine-precision argument the
  recurrence runs in hardware `double` / `double`-complex — the same kernel used
  by the NDArray fast path and `Compile[]`, so a scalar call is bit-identical to
  its packed, NDArray and compiled forms. A packed **complex** buffer stays on
  the buffer (the ND kernel has a complex arm). The generalized three-argument
  `LaguerreL[n, a, x]` threads the **same** machine helper element-wise over a
  packed-List or visible `NDArray` argument `x` (real and complex buffers),
  holding the integer order `n` and machine parameter `a` constant — so its
  scalar, List, packed, NDArray and `Compile[]`-scalar forms are all
  bit-identical, and a million-point array runs at hardware speed rather than
  re-expanding the explicit polynomial per element (a ~20–40× speedup that put
  the form ahead of SciPy's `eval_genlaguerre` and ~50–500× ahead of
  the reference over a packed array). An exact (`"int64"`) buffer, an
  arbitrary-precision argument, a large degree (`n > 1000`), or an overflow
  falls back to the stable exact / MPFR recurrence unchanged.

An exact non-integer order with an exact/symbolic argument (e.g. `LaguerreL[1/2,
x]`) stays symbolic; symbolic `Series` / `SeriesCoefficient`, the Series
expansion at infinity, and the three-argument array `Compile[]` lowering are
deferred.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= LaguerreL[7, x]
Out[1]= 1 - 7 x + 21/2 x^2 - 35/6 x^3 + 35/24 x^4 - 7/40 x^5 + 7/720 x^6 - 1/5040 x^7

In[2]:= Table[LaguerreL[n, a, x], {n, 0, 3}, {a, 0, 1}]
Out[2]= {{1, 1}, {1 - x, 2 - x}, {1 - 2 x + 1/2 x^2, 3 - 3 x + 1/2 x^2}, {1 - 3 x + 3/2 x^2 - 1/6 x^3, 4 - 6 x + 2 x^2 - 1/6 x^3}}

In[3]:= LaguerreL[100, 0, N[5, 20]]
Out[3]= 1.45552716253288126069
```

## HermiteH

`HermiteH[n, x]` gives the physicists' Hermite polynomial `H_n(x)`. They are
orthogonal on `(-∞, ∞)` with weight `e^(-x²)` and satisfy `y'' - 2x y' + 2n y =
0`.

- **Exact polynomial.** A non-negative integer order `n` produces the explicit
  degree-`n` polynomial in `x`, with exact integer coefficients from FLINT's
  `fmpz_poly_hermite_h` (or, without FLINT, the three-term recurrence `H_0 = 1`,
  `H_1 = 2x`, `H_k = 2x H_{k-1} - 2(k-1) H_{k-2}`). For example `HermiteH[10, x] =
  -30240 + 302400 x² - 403200 x⁴ + 161280 x⁶ - 23040 x⁸ + 1024 x¹⁰` and
  `HermiteH[7, x] = -1680 x + 3360 x³ - 1344 x⁵ + 128 x⁷`.
- **Numerics.** For an integer order with an inexact argument the value comes
  from the *stable* recurrence, not from substituting into the explicit
  polynomial — so `HermiteH[100, N[5, 16]] = 8.19319371305915×10⁹⁸` rather than
  the catastrophically-cancelled `0.×10¹⁰⁴` the expanded polynomial gives. A
  non-integer, complex, or negative order with an inexact argument evaluates
  through the entire confluent-hypergeometric form `H_n(z) = 2^n Sqrt[Pi]
  (Hypergeometric1F1[-n/2, 1/2, z²]/Gamma[(1-n)/2] - 2z
  Hypergeometric1F1[(1-n)/2, 3/2, z²]/Gamma[-n/2])`, evaluated with a
  cancellation guard scaled to `Re(z²)` so a large real argument keeps its
  digits. All paths run at machine or arbitrary (MPFR) precision, real or
  complex: `HermiteH[3., 5] = 940.`, `HermiteH[3. - 4.583 I, 5] = -731.082 +
  940.003 I`, `HermiteH[0.71, 0.87] = 1.55872`, `HermiteH[5.3 + I, .8 + I] =
  4.81841 + 157.49 I`, and `N[HermiteH[2/3, 8/7], 50]` is accurate to 50 digits.
- **Special values.** `HermiteH[0, x] = 1` and `HermiteH[n, 0] = 2^n
  Sqrt[Pi]/Gamma[(1-n)/2]` (so `HermiteH[1, 0] = 0`, `HermiteH[10, 0] = -30240`).
- **Derivative.** `D[HermiteH[n, z], z] = 2 n HermiteH[n-1, z]`, which also lets
  the generic `NumericFunction` interval threading handle `Interval[]` arguments.
- **Series.** The Taylor expansion at a finite point follows from the derivative
  rule, e.g. `Series[HermiteH[n, x], {x, x0, 2}] = HermiteH[n, x0] + 2 n
  HermiteH[-1+n, x0] (x-x0) + 2 n (-1+n) HermiteH[-2+n, x0] (x-x0)² + O[x-x0]³`.
  At Infinity the algebraic asymptotic series is produced, `Series[HermiteH[n,
  x], {x, Infinity, 2}] = (2x)^n (1 + (n (1-n))/(4 x²) + O[1/x]⁴)` (the
  subdominant `e^(x²)` branch is not included). `SeriesCoefficient` gives the
  general `x⁰`-expansion term as a `Piecewise` on the parity of the index.
- **Listable.** `HermiteH[3, {1, 2, 3}] = {-4, 40, 180}`.
- **Performance (machine reals and complexes).** At an integer order with a
  machine-precision argument the recurrence runs in hardware `double` /
  `double`-complex — the same kernel used by the NDArray fast path and
  `Compile[]`, so a scalar call is bit-identical to its packed, NDArray and
  compiled forms, and a packed **complex** buffer stays on the buffer (the ND
  kernel has a complex arm). An exact (`"int64"`) buffer, an arbitrary-precision
  argument, a large degree (`n > 1000`), or an overflow falls back to the stable
  exact / MPFR recurrence unchanged.

A non-integer order with an exact/symbolic argument (e.g. `HermiteH[1/2, x]`),
and a negative-integer order with an exact argument, stay symbolic. The
symbolic-index `SeriesCoefficient` is given at `x0 = 0` (the `x0 = 1` holonomic
form in the reference's docs is not produced), and the subdominant exponential
branch of the Infinity asymptotics is not included.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= HermiteH[10, x]
Out[1]= -30240 + 302400 x^2 - 403200 x^4 + 161280 x^6 - 23040 x^8 + 1024 x^10

In[2]:= Table[HermiteH[10, x], {x, 1, 5}]
Out[2]= {8224, 200416, -3093984, 125984224, 3275529760}

In[3]:= HermiteH[100, N[5, 16]]
Out[3]= 8.1931937130591502e+98

In[4]:= D[HermiteH[n, z], z]
Out[4]= 2 n HermiteH[-1 + n, z]
```

## WhittakerM

`WhittakerM[k, m, z]` gives the Whittaker function `M_{k,m}(z)`, a solution of
Whittaker's equation `y'' + (-1/4 + k/z + (1/4 - m²)/z²) y = 0`. It is related to
the Kummer confluent hypergeometric function by `M_{k,m}(z) = e^(-z/2) z^(m+1/2)
Hypergeometric1F1[m-k+1/2, 2m+1, z]`, and is evaluated through that form.

- **Numerics.** With at least one inexact argument the value is computed through
  the confluent form at machine or arbitrary (MPFR) precision, real or complex:
  `WhittakerM[2, 3, 4.5] = 68.967`, `N[WhittakerM[2, 1/3, 1/7], 100]` is accurate
  to 100 digits, `N[WhittakerM[4/5, 2 - I, 2]] = 3.27958 - 3.25602 I`, and
  `WhittakerM[2, 3 - I, 4.5] = -18.2886 - 74.6279 I`. A negative or complex `z`
  returns the complex value off the `z^(m+1/2)` branch cut.
- **Special value.** `M_{k,m}(0) = 0` for `Re(m) > -1/2` (the leading power is
  `z^(m+1/2)`), so `WhittakerM[1, 1, 0] = 0` and `WhittakerM[0, 0, 0] = 0`.
- **Derivative.** `D[WhittakerM[k, m, z], z] = (1/2 - k/z) WhittakerM[k, m, z] +
  ((1/2 + k + m)/z) WhittakerM[k+1, m, z]` (DLMF 13.4.27); `d/dk` and `d/dm` have
  no elementary form and stay as inert `Derivative[1,0,0][WhittakerM][...]` /
  `Derivative[0,1,0][...]`.
- **Series.** At the origin the `z^(1/2+m)` Puiseux form is produced for both a
  concrete and a symbolic `m`, e.g. `Series[WhittakerM[2, 1/2, x], {x, 0, 6}] = x
  - x² + 3/8 x³ - 1/12 x⁴ + 5/384 x⁵ - 1/640 x⁶ + O[x]⁷` and
  `Series[WhittakerM[k, m, x], {x, 0, 4}] = x^(1/2+m) (1 + ... + O[x]⁵)`. At
  Infinity only the *terminating* (elementary) case — `m-k+1/2` a non-positive
  integer — is produced, as the exact closed form, e.g.
  `Series[WhittakerM[2, 1/2, x], {x, Infinity, 5}] = x E^(-x/2) (1 - x/2)`;
  `SeriesCoefficient` at the origin follows for a concrete `m`.
- **Listable.** `WhittakerM[2, 3, {1.5, 2.5, 3.5}] = {2.76149, 13.0096, 34.2394}`.
- **Performance (machine precision).** A machine-precision call — `k`, `m` exact
  or machine real, `z` machine real or complex — uses a raw-`double` kernel
  shared bit-for-bit with the NDArray / packed fast path and `Compile[]` (real
  and complex: `CompileDiagnostics[{{x,_Complex}}, WhittakerM[1, 1/2, x]]`
  reports `Compiled -> True`). Real `z < 0` returns the principal complex value
  (an integer `m + 1/2` gives a real one). Scalar 0.2 µs (real), 0.4 µs
  (complex); exact rational `m` used to fall to the expression template (4.7 µs)
  and its buffers to the element-wise route (5.5 s per 10^6; now 0.08 s).
  Arbitrary precision falls back to the confluent-form evaluator.

The general (non-terminating) asymptotic at Infinity — a two-scale Stokes
transseries — is not produced, `FunctionExpand` does not yet expand it, and the
generic `Interval[]` threading declines (the derivative is self-referential), so
`WhittakerM` over an `Interval[]` stays symbolic.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= WhittakerM[2, 3, 4.5]
Out[1]= 68.967

In[2]:= D[WhittakerM[k, m, z], z]
Out[2]= ((1/2 + k + m) WhittakerM[1 + k, m, z])/z + (1/2 - k/z) WhittakerM[k, m, z]

In[3]:= Series[WhittakerM[2, 1/2, x], {x, 0, 6}]
Out[3]= x - x^2 + 3/8 x^3 - 1/12 x^4 + 5/384 x^5 - 1/640 x^6 + O[x]^7
```

## WhittakerW

`WhittakerW[k, m, z]` gives the Whittaker function `W_{k,m}(z)`, the second
(recessive) solution of Whittaker's equation
`y'' + (-1/4 + k/z + (1/4 - m²)/z²) y = 0`. It is related to the Tricomi confluent
hypergeometric function by `W_{k,m}(z) = e^(-z/2) z^(m+1/2)
HypergeometricU[m-k+1/2, 2m+1, z]`, and is evaluated through that form.

- **Numerics.** With at least one inexact argument the value is computed at
  machine or arbitrary (MPFR) precision, real or complex:
  `WhittakerW[2, 3, 4.5] = 11.3063`, `WhittakerW[6, 4, 1.7] = 1740.55`,
  `N[WhittakerW[3/5, 1/3, 1/7], 30]` to 30 digits,
  `N[WhittakerW[1/5, 2 - I, 2]] = 0.314759 - 1.27345 I`, and the 45-digit complex
  `WhittakerW[2 - 5 I, 3, 4.5\`45]`. A machine-precision call (`k`, `m` exact
  or machine real; `z` machine real or complex) uses the raw-`double`
  `HypergeometricU` kernel (real or complex), shared bit-for-bit with the
  NDArray / packed fast path and `Compile[]`; real `z < 0` gives the principal
  complex value. Arbitrary precision uses the `HypergeometricU` evaluator (FLINT
  `acb_hypgeom_u`). In the common
  integer/half-integer `m` case (`b = 2m+1` an integer) the `double` kernel's
  log series advances its three digammas by the recurrence `ψ(x+1) = ψ(x) + 1/x`
  — three seed evaluations instead of three per term — so the packed path runs
  ~1.8× faster than the SciPy `exp·z^(m+½)·hyperu` composition and ~1000× faster
  than the reference, which has no packed Whittaker path.
- **Values at z = 0.** By the leading exponent `1/2 - |Re m|`: `0` for
  `|Re m| < 1/2` (`WhittakerW[0, 0, 0] = 0`), `ComplexInfinity` for `|Re m| > 1/2`
  (`WhittakerW[0, 1, 0]`), and the finite `1/Gamma[1-k]` for `m = ±1/2`
  (`WhittakerW[k, 1/2, 0] = 1/Gamma[1-k]`; `WhittakerW[2, 0.5, 0] = 0`).
- **Derivative.** `D[WhittakerW[k, m, z], z] = (1/2 - k/z) WhittakerW[k, m, z] -
  WhittakerW[k+1, m, z]/z` (DLMF 13.14.31); `d/dk` and `d/dm` stay inert.
- **Series.** The terminating/elementary case (`a = m-k+1/2` a non-positive
  integer `-n`), via the exact rewrite `HypergeometricU[-n, 2m+1, z] = (-1)^n n!
  LaguerreL[n, 2m, z]`, both at the origin
  (`Series[WhittakerW[2, 1/2, x], {x, 0, 6}] = -2 x + 2 x^2 - 3/4 x^3 + x^4/6 -
  5/192 x^5 + x^6/320 + O[x]^7`) and at Infinity
  (`Series[WhittakerW[2, 1/2, x], {x, Infinity, 5}] = E^(-x/2) (x^2 - 2 x)`). The
  general two-branch origin series and the general asymptotic are not produced,
  as for `WhittakerM`.
- **Listable.** `WhittakerW[2, 3, {1.5, 2.5, 3.5}] = {86.3829, 32.8767, 17.992}`.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= WhittakerW[2, 3, 4.5]
Out[1]= 11.3063

In[2]:= D[WhittakerW[k, m, z], z]
Out[2]= -WhittakerW[1 + k, m, z]/z + (1/2 - k/z) WhittakerW[k, m, z]

In[3]:= WhittakerW[k, 1/2, 0]
Out[3]= 1/Gamma[1 - k]

In[4]:= WhittakerW[2, 1/2, 1.]  (* -E^(-1/2): the U polynomial case *)
Out[4]= -0.606531
```

## ParabolicCylinderD

`ParabolicCylinderD[nu, z]` gives Whittaker's parabolic cylinder function
`D_nu(z)`, the solution of `w'' + (nu + 1/2 - z²/4) w = 0` that decays as
`z -> +Infinity`. It is entire in both `nu` and `z`, and is related to the
confluent functions by `D_nu(z) = 2^(nu/2) e^(-z²/4) HypergeometricU[-nu/2, 1/2,
z²/2]` (`Re z > 0`, DLMF 12.7.14) and by the two-term `Hypergeometric1F1` form
valid everywhere.

- **Exact values.** `ParabolicCylinderD[nu, 0] = 2^(nu/2) Sqrt[Pi]/Gamma[(1-nu)/2]`
  for every order, symbolic included (`ParabolicCylinderD[1, 0] = 0`,
  `ParabolicCylinderD[2, 0] = -1`). Any other exact call stays unevaluated, as in
  the reference system — `ParabolicCylinderD[2, z]` and `ParabolicCylinderD[1, 2]`
  are returned as they are; `N` evaluates them, and `NumericQ` is `True` for
  numeric arguments.
- **FunctionExpand** rewrites an exact integer order `|n| <= 200`: `n >= 0` to the
  Hermite form `2^(-n/2) E^(-z²/4) HermiteH[n, z/Sqrt[2]]` (DLMF 12.7.2), `n < 0`
  through `D_-1 = E^(z²/4) Sqrt[Pi/2] Erfc[z/Sqrt[2]]` and the recurrence
  `D_{k-1} = (z D_k - D_{k+1})/k`.
- **Numerics.** With an inexact argument, at machine or arbitrary precision, real
  or complex order and argument: `ParabolicCylinderD[2.5, 3.] = 1.29848`,
  `ParabolicCylinderD[2.5, -3.] = -0.942235`, `ParabolicCylinderD[2.5, 3. + 2. I]
  = 1.21695 - 6.63911 I`, `N[ParabolicCylinderD[1/3, 2], 30] =
  0.474092755602535213964624995771`. Values outside the double range come back as
  extended-range reals, never `0.` or `Overflow[]`:
  `ParabolicCylinderD[2.5, 60.] = 3.80293*10^-387`.
  - *Machine kernel* (real order `|nu| <= 40`, real `z`; shared bit-for-bit by the
    scalar call, packed / `NDArray` buffers and `Compile[]`): the Poincaré
    expansion for `z >= 9`; for `1.5 < z < 9` the recessive solution integrated
    backward by Taylor steps of the ODE from a point where the expansion holds
    (stable, since the recessive solution grows towards the origin; steps land on
    exactly representable points and `z²/4 - nu - 1/2` is carried in
    double-double); for `-12 <= z <= 1.5` the Maclaurin form summed in double, or
    double-double when the series cancels; for `z < -12` the connection formula
    `D_nu(-t) = cos(pi nu) D_nu(t) + Sqrt[2 Pi]/Gamma[-nu] t^(-nu-1) e^(t²/4) Σ
    (nu+1)_{2s}/(s! (2t²)^s)`. Every exponential is applied through an exactly
    reduced exponent. Complex `z` uses the machine `HypergeometricU` /
    `Hypergeometric1F1` kernels. Whatever the double algorithms decline (a point
    near a zero on the negative axis, `|nu| > 40`, a complex order) is answered
    inside the kernel by arb at 53 bits, so a buffer never leaves its fast path.
    Against mpmath over 2300 random points: real median 0.7 ulp, 99th percentile
    13 ulp; complex median 2.8 ulp, max 31 ulp.
  - *Arbitrary precision:* FLINT/arb, `acb_hypgeom_u` for `Re z > 0` and the
    `acb_hypgeom_m` / `acb_rgamma` form elsewhere, inside a precision ladder that
    also climbs past a confluent kernel that cannot bound its expansion at low
    precision.
- **Derivative.** `D[ParabolicCylinderD[nu, z], z] = z/2 ParabolicCylinderD[nu, z]
  - ParabolicCylinderD[nu + 1, z]` (DLMF 12.8.2); the order derivative stays the
  inert `Derivative[1, 0][ParabolicCylinderD][nu, z]`. `Series` about `z = 0`
  follows from it.
- **Listable**, threads over packed arrays and visible `NDArray`s (real and
  complex buffers, buffers in either slot), and compiles at scalar and rank-1
  shapes (`CompileDiagnostics` reports `Compiled -> True`).
- **Performance** (2.2 GHz Intel, 2026-11-09). Machine scalar ~1.7 µs per call
  in a `Do` loop (wolframscript 1.6 ms; SciPy `pbdv` 12.6 µs from Python); a
  10^6-element packed buffer `ParabolicCylinderD[2.5, x]`, `x` uniform on
  [-5, 5], 1.1 s wall (threaded) where SciPy `pbdv` takes 8.4 s; at 50 digits
  0.14 ms per value (mpmath 0.78 ms); machine complex order 0.12 ms (mpmath
  1.6 ms). The reference system's complex-order values are wrong in the 13th
  digit (`ParabolicCylinderD[2.5 + I, -3. + 2 I]`); Mathilda agrees with mpmath.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= ParabolicCylinderD[2.5, 3.]
Out[1]= 1.29848

In[2]:= ParabolicCylinderD[nu, 0]
Out[2]= (2^(nu/2) Sqrt[Pi])/Gamma[1/2 - nu/2]

In[3]:= D[ParabolicCylinderD[nu, z], z]
Out[3]= 1/2 z ParabolicCylinderD[nu, z] - ParabolicCylinderD[1 + nu, z]

In[4]:= FunctionExpand[ParabolicCylinderD[-2, z]]
Out[4]= E^(-1/4 z^2) - E^(1/4 z^2) Sqrt[1/2 Pi] z Erfc[z/Sqrt[2]]
```

## GegenbauerC

`GegenbauerC[n, lambda, z]` gives the Gegenbauer (ultraspherical) polynomial
`C_n^lambda(z)` and, for any `n`, its continuation
`Gamma[n+2lambda]/(Gamma[2lambda] Gamma[n+1]) Hypergeometric2F1[-n, n+2lambda,
lambda+1/2, (1-z)/2]` (DLMF 18.5.9). `GegenbauerC[n, z]` is `2 ChebyshevT[n,
z]/n`, the `lambda -> 0` limit of `C_n^lambda(z)/lambda`.

- **Exact rules.** An integer `0 <= n <= 400` expands to the explicit polynomial
  `Σ_k (-1)^k (lambda)_{n-k} (2z)^{n-2k}/(k! (n-2k)!)` for every `lambda`
  (`GegenbauerC[3, a, x] = -2 a (1 + a) x + 4/3 a (1 + a) (2 + a) x^3`,
  `GegenbauerC[2, -1/2, x] = 1/2 - x^2/2`, `GegenbauerC[5, -3, x] = -6 x`).
  `GegenbauerC[n, 0, z] = 0` (also for `n = 0`), `GegenbauerC[0, lambda, z] = 1`,
  `GegenbauerC[n, 1, z] = ChebyshevU[n, z]`, `GegenbauerC[n, 1/2, z] =
  LegendreP[n, z]`; a negative integer or non-integer exact degree stays
  unevaluated, as in the reference system. The two-argument form always rewrites
  (`GegenbauerC[n, x] = 2 ChebyshevT[n, x]/n`; `GegenbauerC[0, x]` is
  `ComplexInfinity` with `Power::infy`).
- **Numerics** (machine and arbitrary precision, real and complex in every
  slot). For real `z < -1` and a non-integer degree the value is complex, taken
  above the cut (`GegenbauerC[2.5, 1.5, -3.] = 0.00031844 + 173.394 I`).
  - *Machine kernel* (`sf_machine_gegenbauer_c`, shared bit-for-bit by the
    scalar call, packed / `NDArray` buffers and `Compile[]` at both arities):
    integer degree -- the three-term recurrence in double when no step cancels,
    else in double-double (max 6 ulp against mpmath over 620 points, n up to
    3000); non-integer degree -- the `2F1` series with double-double parameters
    for `z >= -0.2`, Taylor stepping of the Gegenbauer ODE from `z = 0` for
    `-1 < z < -0.2` (where the `2F1` converges slowly), the general machine `2F1`
    otherwise; arb at 53 bits answers whatever these decline (complex
    parameters, Gamma poles, `lambda + 1/2` a non-positive integer). Two-argument
    form: double-double Chebyshev recurrence (integer `n`), `cos(n acos z)` /
    `cosh(n acosh z)` / the exactly reduced `cos(n pi) cosh(n t) + i sin(n pi)
    sinh(n t)` otherwise. Accuracy against mpmath (2000 points): integer degree
    median 0.2 ulp; non-integer degree median 1.9 ulp, p99 110 ulp (the ODE
    stretch near -1); complex median 2 ulp.
  - *Arbitrary precision:* arb -- the recurrence for an integer degree, else
    the regularized `2F1` with `Gamma(lambda+1/2)/Gamma(2lambda)` folded by the
    duplication formula, so `lambda = -1/2, -3/2, ...` has its finite limit
    (`GegenbauerC[2.5, -1.5, 0.3] = 0.666651...`, which `acb_hypgeom_gegenbauer_c`
    and mpmath at the exact parameter both miss).
- **Derivative.** `D[GegenbauerC[n, l, z], z] = 2 l GegenbauerC[n - 1, l + 1, z]`
  (DLMF 18.9.19); `D[GegenbauerC[n, z], z] = 2 ChebyshevU[n - 1, z]` through the
  rewrite; the `n` and `l` derivatives stay inert.
- **Performance** (2.2 GHz Intel). Machine scalar 6.7 µs (wolframscript 52 µs,
  SciPy `eval_gegenbauer` 2.9 µs from Python); 10^6 packed, `n = 20`: 0.27 s
  (wolframscript ~14 s, SciPy 0.11 s); `n = 7.3`: 3.7 s (wolframscript ~146 s,
  SciPy 1.4 s -- SciPy's double recurrences are faster but lose up to ~10^3 ulp
  near zeros); 50 digits 0.33 ms (wolframscript 0.88 ms, mpmath 0.37 ms).

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= GegenbauerC[3, a, x]
Out[1]= -2 a (1 + a) x + 4/3 a (1 + a) (2 + a) x^3

In[2]:= GegenbauerC[2.5, 1.5, 0.3]
Out[2]= -1.58217

In[3]:= D[GegenbauerC[n, a, x], x]
Out[3]= 2 a GegenbauerC[-1 + n, 1 + a, x]
```

## JacobiSN, JacobiCN, JacobiDN and the Jacobi ratios

`JacobiSN[u, m]`, `JacobiCN[u, m]`, `JacobiDN[u, m]` give the Jacobi elliptic
functions `sn(u|m)`, `cn(u|m)`, `dn(u|m)` in the **parameter** convention
`m = k^2` (as `EllipticK`, `EllipticF`). The nine ratios `JacobiCD`, `JacobiCS`,
`JacobiDC`, `JacobiDS`, `JacobiNC`, `JacobiND`, `JacobiNS`, `JacobiSC`,
`JacobiSD` are `pq = p/q` with `n = 1` (`JacobiCD = cn/dn`, `JacobiNS = 1/sn`, ...).
`JacobiAmplitude[u, m]` gives `am(u|m)`, with `Sin[am] = sn`, `Cos[am] = cn`,
continuous in `u`; `InverseJacobiSN[v, m] = EllipticF[ArcSin[v], m]`,
`InverseJacobiCN[v, m] = EllipticF[ArcCos[v], m]`, `InverseJacobiDN[v, m]` the
inverse of `dn`.

- **Exact rules.** `u = 0` (`sn = 0`, `cn = dn = 1`, a ratio over `sn` is
  `ComplexInfinity`), `m = 0` (circular: `JacobiCS[u, 0] = Cot[u]`), `m = 1`
  (hyperbolic: `JacobiCS[u, 1] = Csch[u]`), parity (`JacobiSN[-u, m] =
  -JacobiSN[u, m]`; a ratio with one `s` is odd, the rest even);
  `JacobiAmplitude[u, 0] = u`, `JacobiAmplitude[u, 1] = 2 ArcTan[Tanh[u/2]]` (the
  Gudermannian), `JacobiAmplitude[0, m] = 0`; `InverseJacobiSN[1, m] =
  InverseJacobiCN[0, m] = EllipticK[m]`, `InverseJacobiSN[v, 0] = ArcSin[v]`,
  `InverseJacobiSN[v, 1] = ArcTanh[v]`, `InverseJacobiCN[v, 0] = ArcCos[v]`,
  `InverseJacobiDN[v, 1] = ArcSech[v]`.
- **Branches** (the reference system's). The amplitude winds by `2 Pi` per
  `4K` for every real `m`: `JacobiAmplitude[10., 1.7] = 11.7593`,
  `JacobiAmplitude[10., -2.3] = 13.6764` -- the principal `-I Log[cn + I sn]` plus
  `2 Pi Round[(Pi Re[u]/(2 Re[K]) - Re[phi])/(2 Pi)]`. `InverseJacobiDN[v, m]` for
  `0 <= v < Sqrt[1 - m]` is `K + I y` (`InverseJacobiDN[0.4, 0.3] = 1.71389 +
  1.56202 I`; the plain `EllipticF` form gives the conjugate there).
- **Numerics.** Machine kernel `sf_machine_jacobi`, shared by the scalar call,
  packed / `NDArray` buffers (any slot) and `Compile[]`: for real `u` and
  `0 <= m <= 1`, `u` is reduced modulo `4K` and then to the nearest quarter
  period with `4K = 2 Pi/AGM(1, Sqrt[1-m])` in double-double, and the
  arithmetic-geometric-mean method (DLMF 22.20) runs on the small residual (so
  values near the zeros of `sn`, `cn` keep full relative accuracy). Complex `u`
  uses the addition formulas over `(x|m)` and `(y|1-m)`; real `m > 1` and
  `m < 0` the reciprocal and imaginary-modulus transformations, carried in
  double-double. Complex `m` and anything declined is answered by arb at 53 bits
  inside the kernel. Arbitrary precision: arb theta quotients at
  `tau = I K(1-m)/K(m)`, `z = u/(2K)`. The inverses use `EllipticF`'s machine
  Carlson kernel where real, and `EllipticF[ArcSin[...]]` elsewhere. Accuracy
  against mpmath (2550 points; |u| up to 40, every m regime): real median 0.4
  ulp, max 19 ulp; complex median 0.7 ulp, max 7 ulp; amplitude max 2.6 ulp;
  inverses max 3.6 ulp.
- **Derivatives** in `u` (DLMF 22.13): `D[JacobiSN[u, m], u] = JacobiCN[u, m]
  JacobiDN[u, m]`, ..., `D[JacobiAmplitude[u, m], u] = JacobiDN[u, m]`,
  `D[InverseJacobiSN[v, m], v] = JacobiCD[InverseJacobiSN[v, m], m]/(1 - v^2)`;
  the `m` derivative stays inert.
- **Performance** (2.2 GHz Intel): machine scalar 1.1 µs (wolframscript 6.9 µs,
  SciPy `ellipj` 3.9 µs from Python); 10^6 packed `JacobiSN[x, 0.3]` 0.13 s
  (wolframscript ~13 s, SciPy `ellipj` 0.34 s); 50 digits 56 µs (wolframscript
  74 µs, mpmath 336 µs); machine complex `m` 53 µs (mpmath 1.0 ms).

Attributes: `Listable`, `NumericFunction`, `Protected` (all sixteen).

```mathematica
In[1]:= JacobiSN[0.7, 0.3]
Out[1]= 0.632305

In[2]:= D[JacobiCD[u, m], u]
Out[2]= (-1 + m) JacobiND[u, m] JacobiSD[u, m]

In[3]:= JacobiAmplitude[10., 1.7]
Out[3]= 11.7593
```

## HankelH1, HankelH2 and the spherical Bessel functions

`HankelH1[nu, z] = BesselJ[nu, z] + I BesselY[nu, z]` and `HankelH2[nu, z] =
BesselJ[nu, z] - I BesselY[nu, z]`; `SphericalBesselJ[n, z]`,
`SphericalBesselY[n, z]` are `Sqrt[Pi]/Sqrt[2 z]` times `BesselJ[n + 1/2, z]`,
`BesselY[n + 1/2, z]` (principal `Sqrt[2 z]`, the reference's form), and
`SphericalHankelH1/H2[n, z] = SphericalBesselJ[n, z] +- I SphericalBesselY[n, z]`.
Principal branches, cut along the negative real axis.

- **Exact rules** (as the reference system): exact calls stay unevaluated except
  at `z = 0` -- `HankelH1[0, 0] = ComplexInfinity`, `SphericalBesselJ[0, 0] = 1`,
  `SphericalBesselJ[n, 0] = 0` (integer `n > 0`), `SphericalBesselY[n, 0] =
  ComplexInfinity`. Integer order on the negative axis uses the parity relations,
  so `SphericalBesselY[3, -2.] = -1.48437` is real.
- **Numerics.** Machine kernel `sf_machine_hankel_sph` (scalar, packed /
  `NDArray`, `Compile[]`) over the shared machine Bessel kernels: where `J + I Y`
  would cancel (H1 in the upper half plane, H2 in the lower) the Hankel function
  is `(2/Pi) e^(-+ I Pi (nu+1)/2) K_nu(-+ I z)` (DLMF 10.27.8), elsewhere `J +- I Y`
  (on the real axis just the real and imaginary parts). Spherical functions of
  integer order `n <= 200` and real `z > 0` use the upward recurrence from
  `j_0, j_1, y_0, y_1` (always for `y`; for `j` where `z > n`). Complex orders and
  declined points: arb at 53 bits in the kernel. Arbitrary precision: arb, with
  the same `K` route. Accuracy against mpmath (2400 points): real `z` median
  0.2-1.5 ulp, max 101 ulp (near a zero of `j`); complex `z` median 1.7-2.5 ulp,
  p99 51 ulp, max 1483 ulp (the shared complex Bessel kernel, rare).
- **Derivatives.** `D[HankelH1[nu, z], z] = (HankelH1[nu - 1, z] - HankelH1[nu + 1,
  z])/2`; `D[SphericalBesselJ[n, z], z] = -SphericalBesselJ[n, z]/(2 z) +
  (SphericalBesselJ[n - 1, z] - SphericalBesselJ[n + 1, z])/2` (likewise for the
  other three); the order derivative stays inert.
- **Performance** (2.2 GHz Intel; the machine was shared during these runs):
  `SphericalBesselJ[3, x]` over 10^6 packed reals 0.038 s (SciPy `spherical_jn`
  0.056 s, wolframscript ~42 s); `HankelH1[1.5, x]` 0.32 s (SciPy `hankel1`
  0.15 s -- both evaluate J and Y; ours run the general half-integer Bessel
  kernels at full accuracy, wolframscript ~22 s); scalar `HankelH1[2.5, 3.]`
  1.3-2.7 µs (SciPy 1.6 µs, wolframscript 26 µs), `SphericalBesselJ[3, 2.5]`
  0.7-1.4 µs (SciPy 45 µs, wolframscript 32 µs); 50 digits 55-140 µs (mpmath
  281 µs, wolframscript 14 µs).

Attributes: `Listable`, `NumericFunction`, `Protected` (all six).

```mathematica
In[1]:= HankelH1[2.5, 3.]
Out[1]= 0.41271 - 0.369041 I

In[2]:= SphericalBesselJ[2.5, -2.]
Out[2]= 0. + 0.114273 I

In[3]:= D[SphericalBesselJ[n, z], z]
Out[3]= -SphericalBesselJ[n, z]/(2 z) + 1/2 (SphericalBesselJ[-1 + n, z] - SphericalBesselJ[1 + n, z])
```

## BetaRegularized

`BetaRegularized[z, a, b]` gives the regularized incomplete beta function
`I_z(a, b) = Beta[z, a, b]/Beta[a, b]`, with `Beta[z, a, b] = z^a/a
Hypergeometric2F1[a, 1 - b, a + 1, z]` (its analytic continuation to complex
arguments and `z` outside `[0, 1]`). `BetaRegularized[z0, z1, a, b]` gives
`I_z1(a, b) - I_z0(a, b)`.

- **Exact rules** (the reference system's): a non-positive integer `a` gives `1`, a
  non-positive integer `b` gives `0` (also for machine `0.`); `a = 1` gives
  `1 - (1 - z)^b`, `b = 1` gives `z^a - 0^a`; an exact numeric `z` with positive
  integer `a`, `b` the terminating binomial sum (`BetaRegularized[1/2, 2, 3] =
  11/16`, `BetaRegularized[1/5, 1/2, 2, 3] = 5067/10000`); other exact calls
  stay unevaluated. At `z = 1.` the value is the machine `1.` (the reference
  prints an exact `1`, which would make a packed result two-headed).
- **Machine kernel** `sf_machine_betareg` (scalar, packed / `NDArray` through the
  n-ary map -- one kernel for both arities, `nargs_alt` -- and `Compile[]` at scalar
  and rank 1 for both forms): for `0 <= z <= 1`, `a, b > 0` the continued fraction
  (modified Lentz, DLMF 8.17.22) on `I_z(a, b)` or `1 - I_(1-z)(b, a)`, times the
  prefactor `z^a (1-z)^b/(a Beta[a, b])` -- `pow` and the shared machine `Beta` for
  small parameters, the Stirling form with a `log1p(u) - u` exponent (whose linear
  parts cancel exactly) for `a, b >= 10` -- with `1 - z` as a double-double. A
  complement below `1/8`, an exponent beyond 64, a non-converging fraction, and
  every complex, out-of-range or negative-parameter call go to arb
  (`acb_hypgeom_beta_lower`, regularized) inside the kernel.
- **Arbitrary precision**: `arb_hypgeom_beta_lower` (real `z` in `[0, 1]`) /
  `acb_hypgeom_beta_lower` on the shared precision ladder; the four-argument form
  on its own ladder, so a cancelling difference climbs.
- **Derivatives.** `D` in `z`: `(1 - z)^(b - 1) z^(a - 1)/Beta[a, b]`; in `a`:
  `-z^a HypergeometricPFQ[{a, a, 1 - b}, {a + 1, a + 1}, z]/(a^2 Beta[a, b]) +
  BetaRegularized[z, a, b] (Log[z] - PolyGamma[0, a] + PolyGamma[0, a + b])`; in
  `b` the same through `I_z(a, b) = 1 - I_(1-z)(b, a)` (checked against central
  differences at 30 digits); the four-argument form as the difference.
- **Accuracy** against mpmath at 40 digits (600 points): `0 <= z <= 1`, `a, b` in
  `[0.01, 100]` median 0.9 ulp, max 59 ulp; complex / out-of-range / negative
  parameters (arb) max 0.9 ulp. The reference system prints
  `BetaRegularized[0.3, 200., 300.] = 1.0496985243291939*10^-6` (true
  `1.04969852432929467*10^-6`).
- **Performance** (2.2 GHz Intel): scalar 0.44 µs (wolframscript 1.5 µs; SciPy
  `betainc` 2.2 µs per Python call, mpmath 156 µs); 10^6 packed 0.049 s (SciPy
  0.36 s, wolframscript 1.8 s); `a, b = 200, 300` over 10^5 points 0.058 s (SciPy
  0.037 s); complex 68 µs (wolframscript 95 µs, mpmath 247 µs); 50 digits 130 µs
  (wolframscript 15 µs, mpmath 235 µs). Open: the 50-digit case is arb-bound, and
  far tails at large parameters (`a, b ~ 10^3`, values below `e^-64`) fall back to
  arb at milliseconds per point, where SciPy answers in double.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= BetaRegularized[0.3, 2.5, 3.7]
Out[1]= 0.319003

In[2]:= BetaRegularized[1/2, 2, 3]
Out[2]= 11/16

In[3]:= BetaRegularized[z, a, 1]
Out[3]= -0^a + z^a
```

## AppellF1, AppellF2, AppellF3 and AppellF4

The Appell hypergeometric functions of two variables (DLMF 16.13), in the
reference system's argument order:

- `AppellF1[a, b1, b2, c, x, y]` = Σ (a)_(m+n) (b1)_m (b2)_n / ((c)_(m+n) m! n!) x^m y^n (|x|, |y| < 1)
- `AppellF2[a, b1, b2, c1, c2, x, y]` = Σ (a)_(m+n) (b1)_m (b2)_n / ((c1)_m (c2)_n m! n!) x^m y^n (|x| + |y| < 1)
- `AppellF3[a1, a2, b1, b2, c, x, y]` = Σ (a1)_m (a2)_n (b1)_m (b2)_n / ((c)_(m+n) m! n!) x^m y^n (|x|, |y| < 1)
- `AppellF4[a, b, c1, c2, x, y]` = Σ (a)_(m+n) (b)_(m+n) / ((c1)_m (c2)_n m! n!) x^m y^n (√|x| + √|y| < 1)

each continued analytically; the cuts `x, y` real `> 1` are taken from below
(`x - i0`, the convention of `Hypergeometric2F1`).

- **Exact rules** (the reference system's): `x = 0`, `y = 0` and a vanishing
  numerator parameter reduce to `Hypergeometric2F1`, `a = 0` (F1, F2, F4) gives
  `1`; F1 also `x = y` -> `2F1(a, b1 + b2; c; x)`, `a = c` -> `(1-x)^-b1
  (1-y)^-b2`, `c = b1 + b2` -> `(1-y)^-a 2F1(a, b1; b1 + b2; (x-y)/(1-y))`, `x = 1`
  -> `2F1(a, b1; c; 1) 2F1(a, b2; c - b1; y)` (and `y = 1`), and a negative integer
  `a` the terminating polynomial (`AppellF1[-1, 1, 2, 3, 1/2, 1/3] = 11/18`).
  Other exact calls stay unevaluated; a pole (`c` a non-positive integer) stays
  unevaluated at every precision. Not done: the reference system's elementary
  closed forms for integer parameters (`AppellF1[1, 2, 3, 4, x, y]`).
- **Method.** One outer series Σ_m T_m G_m whose terms carry a Gauss function of
  the other variable, continued analytically: F1 `x^m 2F1(a+m, b2; c+m; y)`, F2
  `2F1(a+m, b2; c2; y)`, F3 `2F1(a2, b2; c+m; y)`, F4 `2F1(a+m, b+m; c2; y)`. The
  outer series converges for `|x| < 1` (F1, F3), `|x| < min(1, |1-y|)` (F2) and
  `|x| < min |1 ∓ √y|^2` (F4) for any `y` off its cut, so either variable can be
  outer. F1 adds the five Euler-integral substitutions mapping `[0, 1]` onto itself
  (DLMF 16.16; two of them valid only on one side of an extra ray, tested), F2 the
  three Pfaff-type ones; the planner takes the smallest outer ratio. The `G_m`
  come from two Gauss values and a three-term recurrence in `m` run in its stable
  direction (backward from directly summed start values when `G_m` is the minimal
  solution: F1 with `|y| < 1`, F3 with `|1 - 1/y| > 1`), or one by one when the
  recurrence's roots are within a factor 1.5. Where no expansion converges (F1
  with both variables large and `|1 - y/x|, |1 - x/y| >= 1`), F1 is its Euler
  integral (`Re c > Re a > 0`): series at the singular endpoints (`0`, `1`, `1/x`
  for an `x` on its cut), `acb_calc_integrate` between. F2, F3, F4 there stay
  unevaluated (open).
- **Machine kernel** `sf_machine_appell` (scalar, packed / `NDArray` through the
  n-ary map, `Compile[]` at scalar and rank 1, real and complex; bit-identical on
  every surface): the plan in double complex. A sum with cancellation (Σ|terms| >
  4 |sum|), an outer ratio above 0.97, a denominator parameter left of 1/2, F3 with
  an inner variable beyond 0.6, complex parameters beyond the machine 2F1, and the
  Euler region go to arb at 53 bits inside the kernel.
- **Arbitrary precision**: the same plan in ball arithmetic on the shared
  precision ladder -- the outer coefficients and the direct Gauss series on
  midpoints with an explicit error bound (a complex ball product overestimates by
  up to √2 a step), the recurrence on midpoints with `8 log2 M` guard bits, the
  start values' radii carried; the truncation bound is asymptotic in the growth of
  `G_m`. Inside `|z| < 0.99` the start values are the direct Gauss series (FLINT
  3.6's `acb_hypgeom_2f1` returns an indeterminate ball for, e.g., `2F1(-0.0721,
  1.15; 102.33; 0.514 - 0.614 i)`), beyond it Pfaff when `|z/(z-1)| < 0.99`.
- **Arbitrary-precision Compile[]**: one row per head in the special-function
  registry (`compile_mgd_special.c`), real and complex.
- **Derivatives.** `D` in `x` and `y` by the contiguous-parameter formulas:
  `∂x F1 = a b1/c F1(a+1, b1+1, b2, c+1)`, `∂y F1 = a b2/c F1(a+1, b1, b2+1, c+1)`;
  `∂x F2 = a b1/c1 F2(a+1, b1+1, b2, c1+1, c2)`; `∂x F3 = a1 b1/c F3(a1+1, a2, b1+1, b2, c+1)`;
  `∂x F4 = a b/c1 F4(a+1, b+1, c1+1, c2)` (and the `y` twins). A parameter derivative
  is the inert `Derivative[...][AppellFk][...]`.
- **Accuracy** against mpmath (100 random points per function, parameters in
  `[-1.5, 3]`, real and complex variables inside the double-series regions, mpmath
  at 40 and 80 digits): machine median / max — F1 1.0 / 24 ulp, F2 0.5 / 36 ulp,
  F3 0.5 / 45 ulp, F4 0.4 / 5 ulp; 30 digits: all within `1.5*10^-35` relative.
  mpmath itself is unreliable near the region boundaries (`appellf4` at `(0.2, 0.3)`
  returns 1.1075 at 60 digits; the value is 1.0248673980162395192, confirmed by the
  direct outer sum), and has no continuation for F1 with both variables large. The
  reference system prints `AppellF1[0.3, 0.2, 0.5, 1.7, 10., -20.]` as `0.661148 -
  0.097611 I` (and `0.6611582714...` at 20 digits); the value is
  `0.66115898412418623440 - 0.09761940125683091707 I`. At `(0.5 + 2 I, 2.)`,
  `(2., 0.5 + 2 I)` and `(2., 3.)` its values disagree with the Euler integral
  entirely (a different continuation).
- **Performance** (2.2 GHz Intel): scalar F1 5 µs, complex 10-300 µs (one-by-one
  Gauss values), outside the series region (`(3., -2.)`) 35 µs, the Euler region
  1.8 ms; F2 2.6-46 µs, F3 2.9 µs (arb route 0.12-8 ms), F4 3.2 µs (9 ms at
  `√x + √y = 0.995`). wolframscript: F1 6.2 ms, complex 10 ms, `(3., -2.)` 75 ms,
  F2 6.3 ms, F3 0.88 ms, F4 6.8 ms (at `(0.2, 0.3)` no answer within an hour);
  mpmath 1.5-6 ms. 10^6 packed: F1 5.4 s, F2 2.1 s, F4 1.3 s, F3 38 s (the arb
  route for `|x| > 0.6`). 50 digits: F1 0.78 ms, F2 0.74 ms, F3 0.53 ms, F4 1.1 ms
  (wolframscript 4.4 / 3.4 / 1.5 / 5.2 ms; mpmath 17 / 18 / 7.7 / 23 ms). Open:
  the hard F3 arb route (8 ms against wolframscript's 0.9 ms) and the F4 boundary.

Attributes: `Listable`, `NumericFunction`, `Protected` (the reference system
gives these heads neither `Listable` nor `NumericFunction`; both are needed for
the packed and `Compile[]` surfaces).

```mathematica
In[1]:= AppellF1[0.3, 0.2, 0.5, 1.7, 0.4, 0.3]
Out[1]= 1.04745

In[2]:= AppellF1[a, b1, b2, c, x, x]
Out[2]= Hypergeometric2F1[a, b1 + b2, c, x]

In[3]:= D[AppellF4[a, b, c1, c2, x, y], x]
Out[3]= (a b AppellF4[1 + a, 1 + b, 1 + c1, c2, x, y])/c1

In[4]:= N[AppellF1[3/10, 1/5, 1/2, 17/10, 10, -20], 20]
Out[4]= 0.66115898412418623440 - 0.097619401256830917067 I
```

## EllipticTheta and EllipticThetaPrime

`EllipticTheta[a, u, q]` gives the Jacobi theta function `theta_a(u, q)`,
`a = 1, 2, 3, 4`, of **nome** `q` (`|q| < 1`): `theta_1 = 2 Σ (-1)^n
q^((n+1/2)^2) Sin[(2n+1) u]`, `theta_2` the same without `(-1)^n` and with `Cos`,
`theta_3 = 1 + 2 Σ q^(n^2) Cos[2 n u]`, `theta_4` with `(-1)^n`; `q^(1/4)` is the
principal root. `EllipticThetaPrime[a, u, q]` is the `u`-derivative;
`EllipticTheta[a, q]` / `EllipticThetaPrime[a, q]` are the values at `u = 0`.

- **Exact rules** (the reference system's): `q = 0` gives `0, 0, 1, 1` (and `0`
  for every derivative); the vanishing values `theta_1(0, q)`, `theta'_2(0, q)`,
  `theta'_3(0, q)`, `theta'_4(0, q)` are `ConditionalExpression[0, Abs[q] < 1]`;
  other exact calls stay unevaluated (`EllipticTheta[3, q]`); a first argument
  outside `1..4` gives `EllipticTheta::etype`.
- **Numerics.** Machine kernel `sf_machine_theta` (scalar, packed / `NDArray`,
  `Compile[]`): for `|q| <= 0.9` the q-series in double complex arithmetic, `Re u`
  first reduced modulo `Pi` in double-double, declining when the series cancels by
  more than `2^3`; otherwise arb (`acb_modular_theta[_jet]` at `tau =
  Log[q]/(I Pi)`, `z = u/Pi`, which reduces `tau` to the fundamental domain), in
  the kernel. The arb ladder requires **relative** accuracy even for a ball that
  contains zero: theta values far below `2^-53` (`EllipticThetaPrime[3, 2.045,
  0.99752] = 1.01*10^-206`) are computed, not guessed. Accuracy against mpmath
  at 500 digits (2000 points; mpmath at 130 digits is itself wrong for `|q| >
  0.99`): median 0.2-1 ulp, max 95 ulp (complex `q`); `|q| > 0.95` max 0.7 ulp.
  The reference system loses this in machine precision: it prints
  `EllipticThetaPrime[4, 0.7, 0.99] = 6.7*10^-16` for the true `5.24*10^-30`.
- **Derivatives.** `D[EllipticTheta[a, u, q], u] = EllipticThetaPrime[a, u, q]`;
  the other partials stay inert `Derivative` forms (as the reference).
- **Performance** (2.2 GHz Intel): scalar 0.5 µs (wolframscript 4.0 µs, mpmath
  37 µs); `|q| = 0.97` 13 µs (arb); 10^6 packed `EllipticTheta[3, x, 0.3]`
  0.075 s (wolframscript ~5.6 s); 50 digits 35 µs (wolframscript 22 µs, mpmath
  55 µs). SciPy has no theta functions.

Attributes: `Listable`, `NHoldFirst`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= EllipticTheta[1, 0.7, 0.3]
Out[1]= 0.838179

In[2]:= EllipticTheta[1, 0, q]
Out[2]= ConditionalExpression[0, Abs[q] < 1]

In[3]:= D[EllipticTheta[3, u, q], u]
Out[3]= EllipticThetaPrime[3, u, q]
```

## WeierstrassP, WeierstrassPPrime, WeierstrassZeta, WeierstrassSigma, WeierstrassHalfPeriods, WeierstrassInvariants

`WeierstrassP[u, {g2, g3}]` gives the Weierstrass elliptic function `wp(u)` with
invariants `g2`, `g3`: the solution of `wp'^2 = 4 wp^3 - g2 wp - g3` with a double
pole at `u = 0`. `WeierstrassPPrime` is its `u`-derivative, `WeierstrassZeta` the
zeta function (`zeta' = -wp`, `zeta ~ 1/u`) and `WeierstrassSigma` the sigma
function (`sigma'/sigma = zeta`, `sigma ~ u`). `WeierstrassHalfPeriods[{g2, g3}]`
gives the half-periods `{w1, w3}` of the period lattice, and
`WeierstrassInvariants[{w1, w3}]` the invariants of the lattice with those
half-periods.

- **Half-period convention** (the reference system's, recovered from its outputs
  on 380 sampled invariants and reproduced on all of them): the basis with
  `tau = w3/w1` in the closed fundamental domain `|Re tau| <= 1/2`, `|tau| >= 1`;
  on its boundary `Re tau = +1/2` (not `-1/2`) and `Re tau >= 0` on `|tau| = 1`;
  then `Re w1 > 0` (`Im w1 < 0` when `w1` is imaginary). Real invariants with
  `g2^3 - 27 g3^2 < 0` therefore give `{w1, w1/2 + I y}` (`w1 > 0`), and with a
  positive discriminant a real and an imaginary half-period. The square lattice
  (`g3 = 0`) takes `arg w1` in `[-Pi/4, Pi/4)` with `w3 = I w1`; the hexagonal one
  (`g2 = 0`) `arg w1` in `[-Pi/6, Pi/6)` and `arg w3` in `[Pi/3, 2 Pi/3)`. Function
  values do not depend on the choice.
- **Exact rules.** `{g2, g3} = {0, 0}`: `1/u^2`, `-2/u^3`, `1/u`, `u`. A double root
  (`g2^3 = 27 g3^2`, `e = -3 g3/(2 g2)`): `e + 3e Csch[Sqrt[3 e] u]^2` (and its
  derivative, `-e u + Sqrt[3 e] Coth[..]`, `E^(-e u^2/2) u Sinc[I Sqrt[3 e] u]`) for
  `e > 0`, the `Cot` / `Csc` / `Sinc` forms for `e < 0` (`WeierstrassP[u, {3, 1}] =
  1 + 3/2 Cot[Sqrt[3/2] u]^2`). `u = 0` is a pole (`Power::infy`, `ComplexInfinity`;
  `sigma(0) = 0`). Parity: `wp` even, the others odd. Exact half-periods of a
  double root (`{Pi/Sqrt[6], DirectedInfinity[I]}`, `{-I Pi/(2 Sqrt[3]), Infinity}`)
  and of `{0, 0}`; the exact invariants of the square (`{Gamma[1/4]^8/(256 Pi^2
  w1^4), 0}`) and hexagonal (`{0, Gamma[1/3]^18/(4096 Pi^6 w1^6)}`) lattices.
  `WeierstrassInvariants` of dependent half-periods emits
  `EllipticReducedHalfPeriods::nind`. The four functions thread over a List `u`.
- **Machine kernel** `sf_machine_weierstrass` (scalar, packed / `NDArray` `u`,
  `Compile[]` at scalar and rank 1; `CompileDiagnostics` reports `Compiled ->
  True`), one plan per `(g2, g3)` (an 8-slot memo): roots of `4t^3 - g2 t - g3`
  (Cardano + Newton), a basis `w1 = K(m)/Sqrt[e1 - e3]`, `w3 = I K(1-m)/Sqrt[e1 - e3]`
  (complex AGM), the normalisation above, and a **certificate** -- the invariants
  recomputed from the theta constants of the reduced `tau` must reproduce
  `(g2, g3)` to `1e-11`, else the next root ordering is tried. The basis, `tau`,
  `1/(2 w1)`, the nome and `eta = zeta(1/2)` are then refined in double-double, so
  `z = u/(2 w1)` is reduced modulo the lattice to `~1e-32 |z|`. With `|q| <= 0.066`
  every theta series has at most 6 terms; the rounding of `v = Pi z0` to a double
  is restored to first order. `wp` and `wp'` come from the theta quotients, `zeta`
  from `2 eta z + Pi theta1'/theta1 - 2 Pi I n`, `sigma` from `theta1` times an
  exponential whose exponent is formed in double-double. A value that cancels by
  more than `2^2` (near a zero of the function) is recomputed by arb from the plan's
  double-double `z` and `tau` (`flint_wp_unit_d`, ~10 µs); a plan the double route
  cannot certify goes to arb whole. Real invariants and real `u` give a `Real`.
- **Arbitrary precision** `flint_num_weierstrass`: the lattice in ball arithmetic
  (`acb_poly_find_roots`, `acb_elliptic_k`, the same integer normalisation, the
  certificate by `acb_elliptic_invariants`), then `acb_elliptic_p` / `p_jet` /
  `zeta` / `sigma` on `(1, tau)`, on the shared precision ladder.
- **Derivatives.** `D` in `u`: `wp'`, `6 wp^2 - g2/2`, `-wp`, `sigma zeta`; in `g2`
  and `g3` the closed forms through `wp`, `wp'`, `zeta`, `sigma` over `g2^3 - 27
  g3^2` (all eight checked against central differences at 40 digits).
- **Accuracy** against mpmath at 50-60 digits (400 random points: real and complex
  invariants of magnitude `10^-6 .. 10^6`, `u` up to 45 periods out): median 2.7 /
  2.6 / 0.7 / 1.1 ulp (`wp`, `wp'`, `zeta`, `sigma`), max 22 / 33 / 11 / 32. The
  reference system on the same points: median 37 / 61 / 4.8 / 17 ulp, and wrong
  (100% error) on about 1%.
- **Performance** (2.2 GHz Intel; wolframscript 14.0): scalar 0.9 µs with the
  plan memoised, 20 µs for new invariants (wolframscript 97 µs / 365 µs); 10^6
  packed real `u` 0.15 s, complex 1.0 s (wolframscript 98 s / 100 s);
  `WeierstrassHalfPeriods` 20 µs (wolframscript 255 µs); 50 digits 0.22 ms
  (wolframscript 0.34 ms; an mpmath theta-function composition 4.1 ms). SciPy and
  mpmath have no Weierstrass functions.

Attributes: `NumericFunction`, `Protected` (`WeierstrassHalfPeriods`,
`WeierstrassInvariants`: `Protected`).

```mathematica
In[1]:= WeierstrassP[0.3, {1., 2.}]
Out[1]= 11.1162

In[2]:= WeierstrassHalfPeriods[{2., 3.}]
Out[2]= {1.19722, 0.59861 + 1.17514 I}

In[3]:= WeierstrassP[u, {12, -8}]
Out[3]= 1 + 3 Csch[Sqrt[3] u]^2

In[4]:= D[WeierstrassZeta[u, {g2, g3}], u]
Out[4]= -WeierstrassP[u, {g2, g3}]
```

## MathieuC and MathieuCPrime

`MathieuC[a, q, z]` gives the even Mathieu function with characteristic value `a`
and parameter `q`, the even solution of `y'' + (a - 2 q Cos[2 z]) y = 0`;
`MathieuCPrime[a, q, z]` is its `z`-derivative. Both are defined for all complex
`a`, `q`, `z`.

- **Definition and normalization.** The Floquet form `Sum[c_k Cos[(nu + 2k) z]]`,
  where `nu` is the characteristic exponent (`Cos[Pi nu] = y_e(Pi)` for the even
  solution with `y_e(0) = 1`) and `(a - (nu+2k)^2) c_k = q (c_{k-1} + c_{k+1})`:
  - generic `nu`: `Sum[c_k^2] = 1` (the plain square, so the function is analytic
    in `a` and `q`), with the sign fixed by `Re(c_0) > 0`, `c_0` the `Cos[nu z]`
    coefficient for the representative of `nu` nearest `Sqrt[a]`;
  - integer `nu` (`a` a characteristic value, `MathieuCharacteristicA[r, q]` for
    an integer `r`): the DLMF normalization `Integrate[ce^2, {z, 0, 2 Pi}] = Pi`,
    i.e. scipy's `mathieu_cem`. A machine `a` within ~2e-6 in `nu` of a
    characteristic value is treated as one;
  - `q = 0`: `Cos[Sqrt[a] z]`, and `1/Sqrt[2]` for `a = 0` (`ce_0`).

  For real `a`, `q` in a stability band (real `nu`), and at every characteristic
  value, this coincides with the reference CAS. In an instability band, or for
  complex `a` or `q` (complex `nu`), the reference uses a non-analytic
  normalization (`Sum[|c_k|^2] = 1` with a phase convention), so there the two
  differ by a `z`-independent constant factor. Both are the even solution.
  Example: `MathieuC[I, -2, 2.5 + I]` is `40.7538 - 38.419 I` here and
  `20.292 - 27.0428 I` there. For real `a`, `q` the value is real or purely
  imaginary on the real and imaginary axes (`MathieuC[4, 1, 0.5] = 0.632114 I`).
- **Exact values.** `MathieuC[a, 0, z] = Cos[Sqrt[a] z]`,
  `MathieuC[0, 0, z] = 1/Sqrt[2]`, `MathieuCPrime[a, 0, z] = -Sqrt[a] Sin[Sqrt[a] z]`,
  `MathieuCPrime[a, q, 0] = 0`. Even in `z`, with the superficial sign rule (as
  `Cos`): `MathieuC[a, q, -x] = MathieuC[a, q, x]`,
  `MathieuCPrime[a, q, -x] = -MathieuCPrime[a, q, x]`. All-exact numeric arguments
  stay unevaluated (`MathieuC[2, 1, 3]`); `N` evaluates them, and `NHoldFirst`
  keeps an exact `a`.
- **Numerics.**
  - Machine: `MathieuC[2, 1, 3.2] = -0.899531`.
  - Arbitrary precision, tracking the input precision:
    `MathieuC[2, 1, 3.2`20] = -0.899531163441461283`,
    `MathieuC[2, 1, 3.2`20 + 0.5`20 I] = -0.882739739030372925 + 0.156794630450643144 I`,
    `N[MathieuC[2, 1, 37/10], 50]` to 50 digits.
  - Large `|z|`: `N[MathieuC[3, 1, 10^15], 20] = -0.78255601989076339546`. Machine
    `MathieuC[3, 1, 10.^15]` is meaningless, as it is for any machine-precision
    computation of `nu 10^15`.
  - Large `Im z`: where the Floquet sum cancels (`MathieuC[-1, 2, 5. I]`, terms
    ~`e^40`, value ~`1e-3`), the machine path declines to the ball kernel, which
    answers correctly. The reference returns `~1e45` there; the value
    `N y_e(5 i)` is confirmed by an independent mpmath ODE integration.
- **Algorithm.** One plan per `(a, q)`, then a short sum per `z`:
  - `nu`: Hill's determinant, rewritten without the poles at `a = (2k)^2` (the
    sinc factor divided out analytically) and with its polynomially convergent
    tail in closed form, refined by Newton on the continued-fraction consistency
    function (every root is an equivalent exponent);
  - the dominant coefficient index comes from inverse iteration on the
    truncated Floquet matrix;
  - coefficients come from two Miller sweeps glued there;
  - the sum is a centred Horner in `e^(2iz)` (complex), or a single real Horner
    pass (real stable).

  The arbitrary-precision kernel (`src/special_functions/mathieu_arb.c`, arb balls
  in the precision ladder) does the same, with three refinements:
  - precision-doubling Newton with a half-precision derivative;
  - a final first-order coefficient correction instead of a last full pass;
  - a per-index precision profile: tail coefficients carry only the bits they
    contribute.

  Clenshaw and the Miller sweeps run on midpoints with a rigorous a-priori
  error bound, because ball radii grow like `1.5^N` in those recurrences while
  the true error stays at the working precision.
- **Performance.**
  - 500 digits: `MathieuC[2, 1, 3.7`500]` in 1.5 ms (reference 16 ms).
  - 10^4 digits: 0.15 s wall, 0.24 s CPU (reference 0.27 s). Above 4096 bits the
    two Miller sweeps, the two Clenshaw sums and the coefficient corrections each
    run as two worker threads.
  - Machine-precision buffers: `MathieuC[2, 1, v]` over 10^6 points in 15 ms.
  - Characteristic order: `ce_2(., 5)` over 10^6 points in 17 ms, against
    scipy `mathieu_cem` 920 ms.
  - Per point: about 80 ns real and 200 ns complex; a new `(a, q)` plan costs
    about 28 µs.
- **Surfaces.** The machine helper `sf_machine_mathieu` is shared bit-for-bit by
  the scalar path, the packed / NDArray map (an n-ary kernel, parallel; a
  complex-valued result widens the buffer), and `Compile[]` at scalar and array
  shapes, real and complex (`CompileDiagnostics[{{v, _Real, 1}}, MathieuC[2., 1., v]]`
  reports `Compiled -> True`). A real-typed compiled call whose value is
  imaginary hands the call to the interpreter.
- **Derivative.** `D[MathieuC[a, q, x], x] = MathieuCPrime[a, q, x]` and
  `D[MathieuCPrime[a, q, x], x] = (-a + 2 q Cos[2 x]) MathieuC[a, q, x]`, so
  `Table[D[MathieuC[a, q, x], {x, n}], {n, 4}]` is closed in `MathieuC` /
  `MathieuCPrime`. `d/da` and `d/dq` stay inert
  (`Derivative[1, 0, 0][MathieuC][a, q, x]`).
- **Series** follows from the derivative:
  - `Series[MathieuC[3, 2, x], {x, 0, 8}] = MathieuC[3, 2, 0] (1 + x^2/2 - 5/8 x^4 - 47/720 x^6 + 473/4480 x^8) + O[x]^9`;
  - at a generic point:
    `Series[MathieuC[a, q, x], {x, x0, 2}] = MathieuC[a, q, x0] + MathieuCPrime[a, q, x0] (x - x0) - 1/2 (a - 2 q Cos[2 x0]) MathieuC[a, q, x0] (x - x0)^2 + O[x - x0]^3`.
- **Listable.** `MathieuC[a, q, {z1, z2, z3}] = {MathieuC[a, q, z1], ...}`.
- **Messages.** `MathieuC[]` emits `MathieuC::argrx` (through the message funnel:
  `Check` and `Quiet` see it).

- **At `a = b_n(q)`** (an odd characteristic value) no periodic even solution
  exists, and the Floquet-normalized even function is identically `0`; that limit
  is returned: `MathieuC[3.917024772998471, 1., 0.7] = 0.` (`b_2(1)`). The
  reference returns a nonzero value there that does not satisfy the Mathieu
  equation.

Not yet provided: `MathieuCharacteristicExponent` as a builtin (the exponent is
computed internally), and asymptotic Series at Infinity. The characteristic
values are `MathieuCharacteristicA` / `MathieuCharacteristicB` (below).

Attributes: `Listable`, `NHoldFirst`, `NumericFunction`, `Protected`.

## MathieuCharacteristicA and MathieuCharacteristicB

`MathieuCharacteristicA[r, q]` and `MathieuCharacteristicB[r, q]` are the
characteristic values `a_r(q)`, `b_r(q)` of Mathieu's equation
`y'' + (a - 2 q cos 2z) y = 0`. For an integer `r`, `a_r` is the value with an
even periodic solution (period `π` for even `r`, `2π` for odd) and `b_r`
(`r >= 1`) the value with an odd one; `a_{-r} = a_r`, `b_{-r} = b_r`. For a
non-integer `r`, `a_r = b_r` is the value whose Floquet exponent is `r`
(continuous in `q` from `r^2`, even in `r` and in `q`).

```
In[1]:= MathieuCharacteristicA[0, 25.]
Out[1]= -40.2568
In[2]:= N[MathieuCharacteristicA[0, 25], 30]
Out[2]= -40.2567795465667872758017193249
In[3]:= MathieuCharacteristicB[2, 1.]
Out[3]= 3.91702
In[4]:= MathieuCharacteristicA[1.5, 2.] == MathieuCharacteristicB[1.5, 2.]
Out[4]= True
```

- Exact arguments stay symbolic except `q = 0` (`r^2`). `NHoldFirst` keeps `N`
  from making an integer order a Real, so `N[MathieuCharacteristicA[0, 25], 30]`
  numericalizes `q` alone. `MathieuCharacteristicB[0, q]` emits
  `MathieuCharacteristicB::zord` (there is no `b_0`).
- Method: the eigenvalue of the symmetric tridiagonal Fourier-coefficient matrix
  of the class (`cos 2kz`, `cos (2k+1)z`, `sin (2k+1)z`, `sin (2k+2)z`, or
  `e^{i (nu+2k) z}` over all `k` for a non-integer order), truncated where the
  eigenvalue's truncation error falls below `2^-(bits+30)`. The wanted value is
  the eigenvalue of fixed rank (the rank of `r^2` on the diagonal at `q = 0`):
  for real `q` the eigenvalues of an unreduced symmetric tridiagonal matrix are
  simple, so the rank never changes.
- Machine: Sturm-sequence bisection in doubles (`sf_mathieu_char_r`), the one
  helper behind the scalar path, the packed / `NDArray` map
  (`MathieuCharacteristicA[2, v]` over a buffer of `q`) and `Compile[]`
  (`CompileDiagnostics[{{x, _Real}}, MathieuCharacteristicA[2, x]]` compiles), so
  the surfaces agree bit for bit.
- Arbitrary precision: the machine value refined by Newton on `det(T - x)` in
  MPFR (the logarithmic derivative is a sum over the `LDL^T` pivots) with the
  precision doubling, then certified by two Sturm counts at `x -+ 2^-(bits+2)`;
  bisection when the certificate fails.
- Complex `q` or `r`: unevaluated (the eigenvalue curves meet at branch points
  off the real axis, so the rank argument does not apply).

Attributes: `Listable`, `NHoldFirst`, `NumericFunction`, `Protected`.

## MathieuS and MathieuSPrime

`MathieuS[a, q, z]` gives the odd Mathieu function with characteristic value `a`
and parameter `q`, the odd solution of `y'' + (a - 2 q Cos[2 z]) y = 0`;
`MathieuSPrime[a, q, z]` is its `z`-derivative. Both are defined for all complex
`a`, `q`, `z`.

- **Definition and normalization.** `Sum[c_k Sin[(nu + 2k) z]]` with the
  coefficients and exponent of `MathieuC`, so `MathieuC + I MathieuS` is the
  Floquet solution `Sum[c_k E^(I (nu + 2k) z)]`
  (`F(z + Pi) = E^(I Pi nu) F(z)`):
  - generic `nu`: `Sum[c_k^2] = 1`, `Re(c_0) > 0`. The representative of `nu`
    is the one nearest `Sqrt[a]`, taking `Re(nu) >= 0` on a tie (`a < 0`). The
    sine sum changes sign with `nu -> -nu`; the tie rule keeps
    `MathieuS -> Sin[Sqrt[a] z]` as `q -> 0`;
  - `a = b_n(q)`: `se_n` with the DLMF normalization (scipy's `mathieu_sem`):
    `MathieuS[3.917024772998471, 1., 0.5] = 0.763792` (`b_2(1)`);
  - `a = a_n(q)` (an even characteristic value): no periodic odd solution
    exists, and the Floquet-normalized odd function is identically `0`, its limit
    from both sides. `MathieuS` returns `0` there. The reference returns a
    nonzero value at `MathieuS[MathieuCharacteristicA[2, 1], 1, 1/2]`, but its
    residual in the equation is about `-4.5`;
  - `q = 0`: `Sin[Sqrt[a] z]`.

  In stability bands (real `nu`) and at every `b_n`, this agrees with the
  reference to the last digit. For complex `nu` (instability bands, complex `a`
  or `q`), Mathilda's value is the analytic normalization: real or purely
  imaginary for real `a`, `q`, `z` (`MathieuS[-5, 3, 0.7] = 3.15987 I`).
  For complex `a` the reference's `MathieuS` does not solve the equation
  (residual about `3` at `MathieuS[I, -2, 3/2 + 3 I/10]`), so
  `MathieuS[I, -2, 2.5 + I]` is `38.4941 + 40.642 I` here and
  `26.8675 + 19.7052 I` there.
- **Exact values.** `MathieuS[a, 0, z] = Sin[Sqrt[a] z]`,
  `MathieuSPrime[a, 0, z] = Sqrt[a] Cos[Sqrt[a] z]`, `MathieuS[a, q, 0] = 0`.
  Odd in `z`: `MathieuS[a, q, -x] = -MathieuS[a, q, x]`, while `MathieuSPrime` is
  even.
- **Numerics.**
  - `MathieuS[2, 1, 3.2] = -0.768098`;
  - `MathieuS[2, 1, 3.2`30] = -0.7680982888758991508884760160`;
  - `N[MathieuS[2, 1, 37/10], 50]`;
  - `N[MathieuS[3, 3, 10^12], 20] = -1.1132442186984920324`.

  Machine `MathieuS[3, 3, 10.^12]` is limited by the rounding of `nu` times
  `10^12`, as in any machine computation.

  Small `|z|` is relatively accurate: the odd functions (`se`, `ce'`) use the
  Taylor series of the equation for `|z| <= 1/(4 Sqrt[1 + |a| + 2|q|])`, where
  the Floquet sum has only absolute accuracy.
- **Performance** (10^6 machine points, wall clock):
  - `MathieuS[2, 1, v]`: 19 ms (reference about 126 s);
  - `se_2(., 5)` (`MathieuS[b_2(5), 5., v]`): 20 ms, against scipy
    `mathieu_sem` 765 ms;
  - 500 digits: 1.2 ms (reference 5.7 ms);
  - 10^4 digits: 0.15 s wall, 0.24 s CPU (reference 0.27 s).
- **Surfaces.** The same machine helper (`sf_machine_mathieu_fn`, sharing the
  plan cache with `MathieuC`) serves the scalar path, the packed / NDArray n-ary
  kernel, and `Compile[]` at scalar and array shapes, bit-for-bit.
- **Derivative and Series.** `D[MathieuS[a, q, x], x] = MathieuSPrime[a, q, x]`
  and `D[MathieuSPrime[a, q, x], x] = (-a + 2 q Cos[2 x]) MathieuS[a, q, x]`.
  `Series[MathieuS[3, 2, x], {x, 0, 8}] = MathieuSPrime[3, 2, 0] (x + x^3/6 - 47/120 x^5 + 113/5040 x^7) + O[x]^9`.
- **Listable**; `MathieuS[]` emits `MathieuS::argrx`.

Attributes: `Listable`, `NHoldFirst`, `NumericFunction`, `Protected`.

## MeijerG

`MeijerG[{{a1, …, an}, {an+1, …, ap}}, {{b1, …, bm}, {bm+1, …, bq}}, z]` is the
Meijer G-function `G^{m,n}_{p,q}(z | a; b)`, the Mellin–Barnes integral

    (1/2πi) ∫ Γ(b1+s)…Γ(bm+s) Γ(1-a1-s)…Γ(1-an-s)
              / (Γ(an+1+s)…Γ(ap+s) Γ(1-bm+1-s)…Γ(1-bq-s)) z^-s ds.

`MeijerG[alist, blist, z, r]` is the generalized form (`s → r s`, times `r`); on
the principal branch it equals `MeijerG[alist, blist, z^(1/r)]`, which is how
every path evaluates it. `MeijerG[alist, blist, z, 1]` is the ordinary function.

- **Canonical form.** The parameters are sorted within each of the four groups
  (`MeijerG[{{1},{2}},{{2/3,1/5},{4}},z]` → `MeijerG[{{1},{2}},{{1/5,2/3},{4}},z]`),
  as in the reference.
- **Consistency.** `m = n = 0`, or a pole collision `a_k − b_j ∈ ℤ⁺` (`k ≤ n`,
  `j ≤ m`, exact or inexact), gives `MeijerG::hdiv` and stays unevaluated — unless
  a Gamma pair cancels: for `b_j − a_k = N ∈ ℤ≥0` with `j ≤ m`, `k > n`,
  `Γ(b_j+s)/Γ(a_k+s) = (−1)^N Γ(1−a_k−s)/Γ(1−b_j−s)`, so the two move groups
  (the reference's rewrite: `MeijerG[{{31,4},{2,3}},{{5,3/2},{1/10}},3]` →
  `−MeijerG[{{2,4,31},{3}},{{3/2},{1/10,5}},3]`). A wrong argument count gives
  `hdiv` then `argt`. All messages route through `Quiet`/`Check`.
- **Automatic reduction** (Slater's theorem, calibrated against the reference):
  `G^{2,0}_{0,2}` → `2 z^((b+c)/2) BesselK[b−c, 2√z]`; `G^{2,1}_{1,2}` →
  `Γ(1−a+b) Γ(1−a+c) z^b HypergeometricU[1−a+b, 1+b−c, z]`; any one-term
  Slater expansion (`m = 1`, `p ≤ q`, or the dual `n = 1`, `p ≥ q`); and a
  two-term one whose pFq is a 0F1, 1F1 or 2F1 with the exponents not an integer
  apart. Each pFq is rewritten to `Exp`, a power, `BesselJ`/`BesselI`,
  `1F1(u; 2u)` → `BesselI`, `EllipticK`/`EllipticE`, or left as a
  `HypergeometricPFQ` divided by the `Gamma`s of its regularized lower
  parameters. For `p = q` a numeric `z` uses the expansion that converges there
  — for `m + n ≤ p` too, where the two expansions are not continuations of each
  other but each is `G` in its own disk (`MeijerG[{{1},{1}},{{1/2},{0}},4]` =
  `Sqrt[Pi]`); a symbolic `z` needs `m + n > p`. On `|z| = 1` only `z = 1` reduces,
  and only when every term's pF(p−1) converges there (parameter excess
  `Re(Σa − Σb) − 1` above 0 at `x = 1`, above −1 at `x = −1`, by Abel's theorem):
  `MeijerG[{{1/2},{1}},{{0,1/2},{}},1]` = `Pi/2`. Everything else — e.g. two 2F2
  terms, or the log case — stays `MeijerG`. The dual named shapes (BesselK /
  HypergeometricU in `1/z`) are left to `FunctionExpand`, as in the reference.
- **Gamma-pair cancellation.** Equal parameters in opposite groups — an `a_k`
  (`k > n`) equal to a `b_j` (`j ≤ m`), or an `a_k` (`k ≤ n`) equal to a `b_j`
  (`j > m`) — cancel in the Mellin–Barnes integrand and are removed:
  `MeijerG[{{},{0}},{{0,0},{}},z]` → `E^-z`,
  `MeijerG[{{0,0},{0}},{{0},{}},z]` → `MeijerG[{{0,0},{}},{{},{}},z]`. A
  cancellation that would leave `m = n = 0` gives `MeijerG::hdiv`.
- **FunctionExpand.** The named shapes in `z` or `1/z`
  (`FunctionExpand[MeijerG[{{0,0},{}},{{},{}},z]]` = `2 BesselK[0, 2/Sqrt[z]]/z`),
  else the named Slater expansion, else — `p = q`, `m + n ≤ p`, symbolic `z` —
  `Piecewise[{{inner, z == 1 || Abs[z] < 1}, {outer, Abs[z] > 1}}, Indeterminate]`
  (the `z == 1` clause when the inner terms converge there):
  `FunctionExpand[MeijerG[{{1},{1}},{{1/2},{0}},z]]` =
  `Piecewise[{{2 ArcSin[Sqrt[z]]/Sqrt[Pi], z == 1 || Abs[z] < 1}, {Sqrt[Pi], Abs[z] > 1}}, Indeterminate]`.
  At `z = 1` with `p = q` the Slater sum is used only when every term
  converges there (the automatic reduction's rule); before, a term
  diverging where `G` is finite gave `ComplexInfinity`
  (`MeijerG[{{1/3,2/3},{}},{{0,1/3},{}},1]`, 5.5009).
- **G^{2,2}_{2,2} in 1 − z** (v0.6223, FunctionExpand only, when the Slater
  route declines):
  `G^{2,2}_{2,2}(z | a₁, a₂; b₁, b₂) = ΠΓ(1 − aᵢ + b_j)/Γ(2 − a₁ − a₂ + b₁ + b₂) z^b₁ 2F1(1 − a₁ + b₁, 1 − a₂ + b₁; 2 − a₁ − a₂ + b₁ + b₂; 1 − z)`.
  - Both sides have the cut z ≤ 0. The 2F1 is analytic at z = 1, so this
    covers the log case (b₁ − b₂ an integer) and the unit argument.
  - Example: `FunctionExpand[MeijerG[{{1/2, 1/2}, {}}, {{0, 0}, {}}, 1]]` =
    `π²`, which gives `∫₀^∞ K₀(x)² dx = π²/4`.
  - Checked numerically at z = 7/10, 13/10, 1/2 and 4. It declines when a
    prefactor Gamma has a pole.
- **pFq parameters over one denominator** (v0.6223). Before reducing a
  `HypergeometricPFQ`, `FunctionExpand` puts each parameter through
  `Together`, so parameters that agree only in that form cancel
  (`1/2 (2 + 2 v)` against `1 + v`).
  - This exposed a bug in the Formula algorithm's final tidy: it turned
    `2F3((3+ν)/4, (1+ν)/4; 1/2, (2+ν)/2, (1+ν)/2; w)` into `Indeterminate`,
    while its ν → −ν mirror was fine.
  - `hr_tidy` now keeps the untidied closed form when the tidy is not
    finite. The form above is now a sum of two `BesselI[ν/2, ±√w]` terms.
- **G at z = 1 by residue summation** (`Hrep\`unitG`). For `p = q` and
  rational parameters the Mellin–Barnes integrand reduces, by the recurrence
  and the reflection formula, to `R(s) Π π/sin(π(c + s))` with `R`
  rational, when every Gamma pairs. `G(1)` is then minus the sum of the
  residues at the right poles (Abel's limit of the Slater series; the
  powers of `Log z` drop out at 1). The sines are periodic up to sign, so
  past the finitely many singular points of `R` the residues are a rational
  function of `k` times `(±1)^k`. That sum is Hurwitz ζ, or the digamma
  difference of the alternating sum, at rational arguments, which reduce to
  Catalan, π and Log forms. `R` must decay, and a harmonic part must cancel
  (sign +1); otherwise the sum declines. This covers the log case the Slater
  route leaves:
  `FunctionExpand[MeijerG[{{0,0,1/2,1/2},{}},{{0,0,1/2},{-1/2}},1]]` =
  `Pi (4 Catalan + Pi Log[2])`, the convolution of `ArcCot[x]` and
  `1/(1 + x)`. It is tried after the certificate walk.
  A definite `Integrate` answers through the same rewrite (its final pass,
  Refining such a Piecewise under the call's assumptions):
  `∫₀^∞ x J₀(ax)/(x²+c²) dx` = `BesselK[0, a c]`, `∫₀^∞ J₀(ax) sin(bx)/x dx` =
  `Pi/2` for `b > a > 0`, `∫₀^∞ e^−x sin x/x dx` = `Pi/4`.
- **z = 0.** `p > q`: `ComplexInfinity` when `n = 1`, else unevaluated; `p ≤ q`:
  `0` when every `Re b_h > 0` (`h ≤ m`), the coefficient when one `b_h = 0`,
  `ComplexInfinity` when some `Re b_h < 0` (with a non-zero coefficient). This
  answers a few cases the reference leaves unevaluated or returns as `Indeterminate`
  (`MeijerG[{{},{}},{{1/2,1/3},{}},0]` is `0`).
- **Numerics.** Machine (real and complex) through `sf_machine_meijerg.c`: the
  Slater sum in `w = z^(1/r)` (`p < q`, or `p = q` with `|w| < 1`) or in `1/w`
  (`p > q`, or `|w| > 1`), with the regularized pFq of every term from the shared
  hypergeometric kernels and the `z`-independent coefficients
  `∏ π/sin(π(b_j−b_h)) ∏Γ(1+b_h−a_j)/∏Γ(a_j−b_h)` precomputed once per parameter
  set. It declines — and the same point is re-evaluated in arb and rounded — for
  more than 6 bits of cancellation between terms, or when a term kernel declines,
  so a machine answer is always within a few ulp (25 reference values agree
  to ≤ 6×10⁻¹⁵).
  The **log case** (two of `b1…bm` an integer apart, where the `π/sin` factors
  have cancelling poles) has its own machine kernel (v0.446): the residue series
  of the Mellin–Barnes integrand, `G = Σ_k w^(β+k) Σ_j c[k][j] (log w)^j` per
  integer-spaced class with base `β`, `j` below the class size. The `c[k][j]` do
  not depend on `z`, so the plan tabulates them once (256 terms): the Laurent
  expansion about `β` from arb's Γ / 1/Γ series, then the exact step
  `F(s+1) = F(s) R(s)` (a rational `R`) in double-double. A point is summed in
  double, re-summed in double-double (`log w` and `w` included) when the measured
  cancellation exceeds 3 bits, and handed to arb only past 45 bits or the table's
  end (`|w|` near 1 for `p = q`, very large `|w|` for `p < q`). The dual's
  parameters `1 − b` are carried exactly, as `G` can amplify their rounding
  ~10²×. Machine log-case values agree with 30-digit references to ≤ 4.5×10⁻¹⁶
  (36 points: triple poles, equal `b`'s, complex parameters, `n > 0`, `p = q`
  both sides, `p > q`, `r ≠ 1`). `MeijerG[{{},{2}},{{1/2,3/2},{}},z]` on a 10⁶
  buffer went from 119 s (per-point arb) to ~0.3 s, and from 456 µs to ~3 µs per
  scalar call; a log-case plan is memoised, so repeated scalar calls build its
  table once.
  Arbitrary precision (`N`, or any MPFR argument) uses arb throughout: regularized
  `acb_hypgeom_pfq` terms, the log case by perturbing `b_j → b_j + (j+1)ε` and
  comparing two `ε` until they agree to the target, `p = q = 2` at `z = 1` (where
  each 2F1 term diverges but `G` does not) by the Gauss regular parts. The output
  precision is the lowest input precision; a machine Real anywhere gives a
  machine answer.
- **Buffers and Compile.** A packed List / NDArray `z` is evaluated in one pass
  with the plan built once (`nd_parallel_for`), bit-for-bit equal to the scalar
  path; `Compile[]` lowers MeijerG at scalar and rank-1 `z`, real and complex
  (`OP_KERNN` / the strip-mined `OP_VKERNN_*`), and a real-typed compiled call
  whose value is complex (`z < 0`) falls back to the interpreter rather than
  dropping the imaginary part. `CompileDiagnostics[{{v,_Real,1}}, MeijerG[{{1/3},{}},{{1/2},{1}},v]]`
  reports `Compiled -> True`.
- **Derivative.** `D[MeijerG[{{a},{a'}},{{b},{b'}},z],z] =
  MeijerG[{{−1, a−1},{a'−1}},{{b−1},{0, b'−1}},z]` (the reference's form); for the
  generalized form `(1/(r z)) MeijerG[{{0,a},{a'}},{{b},{1,b'}},z,r]`.
- **Series.** At 0 (`p ≤ q`) and at `Infinity` (`p ≥ q`) from the truncated Slater
  expansion, including the log case (exact Laurent expansion in `ε` of every
  Gamma factor, giving `Log[z]`, `EulerGamma` and `PolyGamma` terms):
  `Series[MeijerG[{{},{}},{{a,a+1},{}},z],{z,0,1}] = z^a (1 + (−1 + 2 EulerGamma +
  Log[z]) z + O[z]²)` (via the `BesselK` it reduces to). At a regular point,
  Taylor via `D`.
- **Performance.** Scalar machine calls 9.5–16 µs (the log case 0.46 ms, through
  arb); 10^6 packed elements 0.29 s; 50 digits 0.24–1.5 ms; 1000 digits 27 ms.
  The fastest engine in every MR / C / AP cell of the special-functions shootout
  (mpmath 2–234x slower, the reference 11–670x).
- **Threads** over a List in its third argument; not `Listable`.

Gaps: the general-term `SeriesCoefficient` at `x = 1`; the exponential (Stokes)
asymptotics at `Infinity` for `p < q`; `p = q ≥ 3` exactly on `|w| = 1` (no
continuation; stays unevaluated); `Hypergeometric*Regularized` heads do not exist,
so reduced forms carry `HypergeometricPFQ[…]/Gamma[…]` instead.

Attributes: `NumericFunction`, `Protected`.

```mathematica
In[1]:= N[MeijerG[{{},{2}},{{1/2,3/2},{}},3],50]
Out[1]= 0.0415911598753969836703303907497787586265297151590396

In[2]:= MeijerG[{{},{}},{{v},{-v}},z]
Out[2]= BesselJ[2 v, 2 Sqrt[z]]

In[3]:= MeijerG[{{1/2},{}},{{0},{0}},x]
Out[3]= Sqrt[Pi] BesselI[0, 1/2 x] E^(-1/2 x)

In[4]:= MeijerG[{{3+I},{2.}},{{2/3,1/5},{4}},3]
Out[4]= -0.664578 - 3.02966 I

In[5]:= MeijerG[{{1/3},{}},{{1/2},{1}},{0.1, 0.3, 0.5},1/2]
Out[5]= {0.0511285, 0.126112, 0.133903}

In[6]:= D[MeijerG[{{a1,a2},{a3,a4}},{{v},{-v}},z],z]
Out[6]= MeijerG[{{-1, -1 + a1, -1 + a2}, {-1 + a3, -1 + a4}}, {{-1 + v}, {0, -1 - v}}, z]
```

## MeijerGReduce

`MeijerGReduce[expr, x]` rewrites `expr` in terms of inactive Meijer G
functions of `x`, for `x > 0`:

- a factor whose G argument is `c x` becomes
  `Inactive[MeijerG][{{a1..an},{an+1..ap}}, {{b1..bm},{bm+1..bq}}, c x]`;
- a factor whose argument is `c x^q`, `q` an integer `≠ 1`, uses the
  four-argument form `Inactive[MeijerG][alist, blist, w, 1/q]` with
  `w = c^(1/q) x` (so `w^q = c x^q` exactly; `MeijerG[a, b, w, r]` is
  `MeijerG[a, b, w^(1/r)]`): `Sin[x]` → `Sqrt[Pi] Inactive[MeijerG][{{},{}},{{1/2},{0}}, x/2, 1/2]`;
- a fractional `q` keeps the three-argument form at `c x^q`
  (`Sin[Sqrt[x]]` → `Sqrt[Pi] Inactive[MeijerG][{{},{}},{{1/2},{0}}, x/4]`).

Sums are mapped term by term; in a product the constants and powers of `x` stay
outside and every other factor is converted on its own, so a product of two
transcendental factors is a product of two G functions (the reference sometimes
fuses such a product into one G; both activate to `expr`). A positive integer
power of a convertible factor is that power of its G; `BesselJ[ν, t]^2` has its
own single-G entry. `E^(a + b x)` splits into `E^a` times the G of `E^(b x)`.
`UnitStep[1 - x]` is `G^{1,0}_{1,1}(x | 1; 0)`, and together with a
`(1 - x^q)^(a-1)` or `Log[x^q]` factor (`q > 0`) it is a single G
(`(1 - x)^(a-1) UnitStep[1 - x]` → `Gamma[a] Inactive[MeijerG][{{},{a}},{{0},{}}, x]`).

The table is the one `Integrate[..., Method -> "Roach"]` decomposes integrands
with (`Roach`GForm1` in `src/internal/RoachIntegrate.m`, PBM vol. 3 sec. 8.4), so
the two cannot drift apart: `Exp`, `Sin`/`Cos`, `Sinh`/`Cosh`, `Log[1 + c x^q]`,
`ArcTan`/`ArcCot`, `(b + c x^q)^(-a)`, `Erf`/`Erfc`, `ExpIntegralEi[-c x^q]`,
`ExpIntegralE`, `SinIntegral`/`CosIntegral`, `FresnelS`/`FresnelC`,
`Gamma[a, t]`, `BesselJ`/`BesselY`/`BesselI`/`BesselK`, `StruveH`,
`EllipticK`/`EllipticE`, `BesselJ[ν, t]^2`. A factor with no G form is kept as it
is; when no factor of `expr` has one (`MeijerGReduce[foo[x], x]`, `Log[x]` with
no indicator) the call stays unevaluated. A constant or a bare power of `x` is
returned unchanged. A non-symbol second argument stays unevaluated; a wrong
argument count gives `MeijerGReduce::argr`/`argrx`.

`Activate[MeijerGReduce[f, x]]` equals `f` on `x > 0`; `tests/test_meijergreduce.c`
checks this numerically for 50 inputs. The module is loaded lazily on the first
call.

Attributes: `Protected`.

```
In[1]:= MeijerGReduce[x^2 BesselJ[0, x] + Log[1 + a x], x]
Out[1]= Inactive[MeijerG][{{1, 1}, {}}, {{1}, {0}}, a x] + x^2 Inactive[MeijerG][{{}, {}}, {{0}, {0}}, 1/2 x, 1/2]

In[2]:= MeijerGReduce[E^(-x) BesselJ[0, x], x]
Out[2]= Inactive[MeijerG][{{}, {}}, {{0}, {}}, x] Inactive[MeijerG][{{}, {}}, {{0}, {0}}, 1/2 x, 1/2]

In[3]:= Activate[MeijerGReduce[ArcTan[x], x]] /. x -> 0.7
Out[3]= 0.610726
```

## FunctionExpand

`FunctionExpand[expr]` (and `FunctionExpand[expr, assum]`, the assumptions
accepted but not yet used) expands special functions in `expr` into simpler
ones, bottom-up (`src/special_functions/functionexpand.c`):

- **MeijerG** → its named shape (`BesselK`, `HypergeometricU`) in `z` or `1/z`,
  else its Slater expansion — a sum of `HypergeometricPFQ` terms in `z`
  (`p ≤ q`) or `1/z` (`p > q`), each in a classical form where it has one — else,
  for `p = q` with `m + n ≤ p` and a symbolic `z`, the Piecewise of the inner and
  outer expansions (see `MeijerG`). The log case is returned unchanged.
- **HypergeometricPFQ** → a closed form by Roach's Formula algorithm (v0.6203,
  see *Formula representations* below):
  `FunctionExpand[Hypergeometric2F1[-3/2, -1/2, 1/2, z]]` =
  `1/2 Sqrt[1 - z] (2 + z) + 3/2 Sqrt[z] ArcSin[Sqrt[z]]`. Unchanged when the
  algorithm declines. Automatic evaluation of `HypergeometricPFQ` is not
  affected.
- **PolyGamma[n, q]** at a positive rational `q` (v0.6131): integer `q`, every
  order, as `n! ((-1)^(n+1) (Zeta[n+1] - HarmonicNumber[q-1, n+1]))`
  (`2 (9/8 - Zeta[3])` for `PolyGamma[2, 3]`); half-integer `q`, every order;
  digamma at any rational by the Gauss digamma theorem. See `PolyGamma`.
- **PolyLog[n, w]**, integer `n ≥ 2`, `w = E^(I Pi r)` an exact root of unity
  (spelled with `I`, `−1`, `(−1)^r`, `E^(I Pi r)` and products of them) →
  `Re + I Im` with `x = r/2 mod 1 = p/q`:
  the Bernoulli part (Re for even `n`, Im for odd `n`) is always
  `(2π)^n B_n(x)/(2 n!)` up to sign; the Clausen part is a multiple of `Zeta[n]`
  for a cosine sum with `q ∈ {1,2,3,4,6}` (Ramanujan sums:
  `ζ(n) Σ_{d|q} μ(q/d) d^(1−n)/φ(q)`), `±Catalan` for `n = 2, q = 4`, and
  otherwise the Hurwitz decomposition `q^(−n) Σ_j f(2πjp/q) ζ(n, j/q)` in
  `PolyGamma` terms (as the reference writes it). The Clausen part is built
  from the representative `p ≤ q/2`, so a conjugate pair cancels structurally:
  `FunctionExpand[PolyLog[3, I] + PolyLog[3, −I]]` = `−3 Zeta[3]/16`,
  `FunctionExpand[PolyLog[3, E^(2 Pi I/3)]]` = `2 I Pi^3/81 − 4 Zeta[3]/9`.
- **PolyLog[2, 2]** → `Pi^2/4 - I Pi Log[2]` and **PolyLog[3, 2]** →
  `7/8 Zeta[3] + Pi^2/4 Log[2] - I Pi Log[2]^2/2`, by the inversion formula with
  `Li_n(1/2)` closed, on the branch PolyLog evaluates (the value approached from
  below the cut, so the imaginary part is negative). Automatic evaluation leaves
  both as they are.
- **PolyLog[2, z]**, z an exact number one anharmonic map (`1 - z`, `1/z`,
  `z/(z - 1)`, `1/(1 - z)`, `(z - 1)/z`) away from a root of unity, `1/2` or `2`,
  by the reflection, Landen and inversion formulas; each candidate is accepted
  only when it matches `PolyLog[2, z]` to 30 digits (a wrong branch is off by a
  multiple of `2 Pi I Log` or `Pi^2`). `FunctionExpand[PolyLog[2, (1 + I)/2]]` =
  `I Catalan + 5/96 Pi^2 - I/8 Pi Log[2] - 1/8 Log[2]^2`,
  `FunctionExpand[PolyLog[2, 1 + I]]` = `I Catalan + Pi^2/16 + I/4 Pi Log[2]`;
  `PolyLog[2, 2 I]` stays.
- In the final pass of a definite integral only: a conjugate pair
  `c (PolyLog[3, z] + PolyLog[3, Conjugate[z]])` with `Re z = 1/2` or
  `Re 1/z = 1/2` becomes `2 c Re Li_3(z)` by Landen's three-term identity (and
  inversion), checked to 30 digits: `Li_3(1 + I) + Li_3(1 - I) =
  Pi^2 Log[2]/16 + 35 Zeta[3]/32`. The log of an exact Gaussian rational is split
  into real and imaginary parts there too.
- **Log[w]** → `I Pi r`, **Log[1 + w]** → `Log[2 Cos[Pi r/2]] + I Pi r/2`
  (`|r| < 1`) for the same roots of unity.
- **PolyGamma[m, x]**, `x ∈ {1/2, 1/4, 3/4}`, via `ζ(n, 1/2) = (2^n−1) ζ(n)` and
  `ζ(n, 1/4 | 3/4) = ((4^n − 2^n) ζ(n) ± 4^n β(n))/2` where the Dirichlet beta is
  closed (`n = 2`: Catalan; odd `n`: `EulerE`): `FunctionExpand[PolyGamma[1, 1/4]]`
  = `8 Catalan + Pi^2`; `PolyGamma[3, 1/4]` stays.

Every rule is checked to 30 digits against the unexpanded numerics in
`tests/test_meijerg.c` and `tests/test_polylog.c`. Everything else is returned
unchanged; a packed / NDArray buffer passes straight through.

A definite `Integrate` runs the same rules as its final pass
(`fe_reduce_definite_result`), keeping a MeijerG rewrite only when no MeijerG,
pFq or unresolved Piecewise is left, and a root-of-unity rewrite only when the
(expanded) result is strictly smaller with no new `PolyGamma`:
`∫₀¹ ArcTan[x]/x dx` = `Catalan`, `∫₀^{π/2} x Log[2 Cos x] dx` = `−7 Zeta[3]/16`.
Last, a value containing `Log` is rewritten with `Log[1/n]` → `−Log[n]` and
expanded, kept only when strictly smaller, so endpoint values such as
`Log[1/2]` and the `Log[2]` of `PolyLog[3, 1/2]` merge into one sum of
`Zeta[3]`, `Pi^2 Log[2]`, `Log[2]^k` and rationals.

Attributes: `Protected`.

```mathematica
In[1]:= FunctionExpand[MeijerG[{{1/2,1},{3,4}},{{1/5},{}},z]]
Out[1]= 1/12 Gamma[1/5] Sqrt[Pi] HypergeometricPFQ[{1/5}, {1/2, 3, 4}, 1/z] - (256 Gamma[7/10] HypergeometricPFQ[{7/10}, {3/2, 7/2, 9/2}, 1/z])/(1575 Sqrt[Pi] Sqrt[z])

In[2]:= FunctionExpand[Sin[x] + MeijerG[{{},{}},{{0},{}},z]]
Out[2]= Sin[x] + E^(-z)

In[3]:= FunctionExpand[MeijerG[{{},{2}},{{1/2,3/2},{}},z]]
Out[3]= MeijerG[{{}, {2}}, {{1/2, 3/2}, {}}, z]
```

**Point tables (v0.6491).** `PolyLog[2, z]` at the golden-ratio points φ⁻², φ⁻¹,
−φ⁻¹ and −φ (Landen). `EllipticK` at the singular moduli (2 ∓ √3)/4,
(√2 − 1)² and 1 − (√2 − 1)² (Chowla–Selberg). Each argument is screened against
a double value, and the match is confirmed exactly with `RootReduce`.

## Formula representations (Roach's algorithm)

`FunctionExpand[HypergeometricPFQ[{a1, …, ap}, {b1, …, bq}, z]]` and the
diagnostic head `` HypergeometricPFQ`Representation[{a…}, {b…}, z] `` write a
pFq in closed form with Kelly Roach's Formula algorithm (*Hypergeometric
Function Representations*, ISSAC '96). The code is in
`src/special_functions/hrep_core.c` (the algebra) and `hrep_certs.c` (the
lookup certificates).

- **Module.** With θ = w d/dw, F(a; b; w) satisfies
  `[w Π(θ + a_j) − θ Π(θ + b_k − 1)] F = 0`, of order d = max(p, q + 1).
  Every contiguous F(a + m; b + n; w) lies in the Q(params)(w)-module
  spanned by the companion basis `[f0, θ f0, …, θ^(d−1) f0]` of an origin
  F(a0; b0; w) = f0. θ acts on a coefficient vector through the origin's ODE,
  so a **certificate** is just `(a0, b0, f0)`: Roach's derivative matrix M is
  never entered by hand.
- **Shifts.** `A_i: a_i → a_i + 1` is `1 + θ/a_i` and `B_i: b_i → b_i − 1`
  is `1 + θ/(b_i − 1)`. Each inverse is one polynomial division of the ODE by
  the shift's linear θ-factor, and the remainder constant is Roach's
  definedness condition: `(1 − a_i) Π(b_k − a_i) ≠ 0` for lowering `a_i`,
  `b_i Π(a_j − b_i) ≠ 0` for raising `b_i`.
- **Strategy.** Parameters split into residue classes: symbolic part, plus
  rational part mod 1. A certificate matches when its classes hold the same
  upper and lower counts, and parameters pair in sorted order. Down-movers
  are walked in ascending target order and up-movers in descending order
  (Roach's Y-below-Y₀ rule). Every step is checked as it is taken, and an
  undefined step abandons that candidate for the next. When a static
  certificate is in the class it beats a family certificate, so that
  half-integer Bessel answers come out in `Sinh` and `Cosh`. After that,
  the shortest walk wins.
- **Derivative formula.** When an upper parameter exceeds a lower one by
  n ∈ Z⁺, the pair is removed and `Π_{m<n} (θ + b + m)/(b + m)` is applied to
  the solved lower-order function (Roach '96 §17).
- **Lerch Φ formula** (Roach '96 §18). When every upper parameter pairs by
  an integer difference with a lower parameter, or with the 1 of k!, the
  coefficient R(k) is rational. Apart gives `R = P(k) + Σ q/(k + r)^j`, and
  `F = P(θ) 1/(1 − z) + Σ q Φ(z, j, r)`. At a rational r, Φ is
  range-reduced into (0, 1] and written in PolyLog over the n-th roots of
  z. For j = 1 that sum takes the real Log/ArcTan form. A real negative z
  with odd n uses the real root, so the answer stays real. `LerchPhi` is
  Σ zᵏ/((a+k)²)^{s/2}, so it is used only for a symbolic r.
- **Integral formula** (generated certificates). When a lower parameter
  exceeds an upper one by an integer, the engine:
  1. solves the lower-order F(c; d; t) in closed form;
  2. integrates it, `F(a, c; a+1, d; w) = a w^(−a) ∫₀^w t^(a−1) F(c; d; t) dt`,
     with `a` shifted into (0, 1];
  3. walks from that certificate.

  Each `Integrate` probe has a 3-second deadline. A generated certificate,
  and a Lerch result, is checked against the series at w = 1/7 and −1/5
  (parameters at sample rationals) before it is used.
- **PFD duplication formula** (Roach '96 §19, generated certificates).
  Write Δ(c, n) = c/n, (c+1)/n, …, (c+n−1)/n, so Π_i ((c+i)/n)_l =
  (c)_{nl}/n^{nl}, and count the 1 of l! as a lower parameter. When the upper
  parameters split into groups Δ(c_j, n) and singles, the lower ones into
  groups Δ(d_j, n) and singles, and every single upper pairs with a single
  lower by an integer difference, then with R(l) = Π (u)_l/(v)_l (rational)
  and x = n^(#lower groups − #upper groups) w^(1/n),
  `F = (1/n) Σ_{k<n} H(ω^k x)`, ω = e^(2πi/n),
  `H(x) = Σ_m R(m/n) (c)_m/(d)_m x^m`.
  R(m/n) is again a Pochhammer ratio in m, so H is a pFq whose extra
  parameters differ in pairs by one: Roach's partial fraction decomposition
  of R, left in product form for the derivative and Lerch formulas to take
  apart. H is solved recursively. A pair parameter of H that would be a
  nonpositive integer is range-reduced by shifting the series,
  `F = Σ_{l<s} t_l w^l + t_s w^s F(a+s, 1; b+s, 1+s; w)`, which keeps every
  group and pair. For n ≥ 3 an exponential `e^((a + b i) x)` is put in the
  real form `e^(a x)(cos b x + i sin b x)`, so conjugate roots combine.
  The formula is a certificate generator (Roach's "subroutine of Lookup"):
  the search finds every grouping, the origin is the target with each group
  moved to the nearest exact Δ(c, n) (a group holding the 1 of l! is pinned;
  an origin with a cancelling pair is skipped), f0 is checked against the
  series, and the walk reaches the target. It runs after the derivative
  formula and before integration; the auxiliary H uses only the cheap
  routes until integration is allowed. `0F2(; 1/3, 2/3; z)` =
  `(e^(3 z^(1/3)) + 2 e^(−3 z^(1/3)/2) cos(3√3 z^(1/3)/2))/3` and its lattice
  (Roach's gallery 1F3(1/2; −1/3, −1/2, −2/3; z), in his printed form) and
  `0F3(; 1/4, 1/2, 3/4; z)` = `(cos 4z^(1/4) + cosh 4z^(1/4))/2` come out
  this way.
- **Transformations.** When the direct routes fail, 2F1 also tries its
  Euler and two Pfaff images and 1F1 its Kummer image. The cheap routes run
  on every image before any integration is attempted.
- **Values at z = 1, −1, 1/2.** The classical summation theorems come first:
  - Gauss (convergent), Kummer, Gauss's second theorem, Bailey;
  - Dixon, Watson, Whipple;
  - 3F2(a, b, 1; c, 2; 1), with its digamma limit at b = 1.

  These are tried over parameter permutations, together with their
  derivative-formula reductions (`Σ β_m z^m (A)_m/(B)_m F(A+m; B+m)`).
  Otherwise the closed form is substituted, and if that is not finite (a
  removable singularity) the one-sided limit is taken. The engine never
  returns `Indeterminate`.
- **Presentation.**
  - The 0F1 family is written in the entire kernel `w^(−ν/2) I_ν(2√w)`.
    Its orders are reduced to two adjacent ones (0 and 1 for integer
    orders), and a negated argument gives BesselJ exactly.
  - A constant result has its Gamma at rationals shifted into (0, 1] and
    reflected past 1/2, then Simplify is applied, capped at 2 seconds.
    A constant written with complex numbers (a Lerch Φ value at a negative
    argument takes the principal, complex root) is first normalised: roots
    of unity become `Cos + I Sin`, every argument of an elementary function
    is put in `a + b I` form, the whole is `ComplexExpand`-ed, and a zero
    imaginary part (`PossibleZeroQ`) is dropped. For a constant every step
    is exact; the result must agree with the input to 30 digits, or the
    input stands. `FunctionExpand[Hypergeometric2F1[1, 1, 2/3, 8/9]]` =
    `9 + 3 Sqrt[3] Pi − 3 Log[3]`.
  - Bessel functions at quarter turns are named by the Kelvin functions.
    When every kernel is a BesselI / BesselJ of integer order, the answer is
    a Laurent polynomial in w^(1/4) times entire functions of it, so any
    fourth root may stand for w^(1/4); at a negated argument w = −u the root
    `(−1)^(1/4) u^(1/4)` is taken. Then `I_n(y) = i^(−n) J_n(i y)` and
    `J_n(e^(3πi/4) Y) = ber_n Y + i bei_n Y`, `J_n(e^(5πi/4) Y) = ber_n Y − i bei_n Y`
    (DLMF 10.61.1) with `J_n(−y) = (−1)^n J_n(y)`, chosen by the phase of the
    argument's numeric coefficient. Integer Kelvin orders are reduced to 0
    and 1 by DLMF 10.63.1.
- **Probe region.** The internal Integrate, Limit and check messages are
  suppressed and do not trip an enclosing `Check`.
- **Certificates** (v0.6203):
  - 0F0;
  - 1F0 (family);
  - 0F1 at b = 1/2 and 3/2, plus the Bessel-I family;
  - 1F1: (1/2; 3/2) Erfi, (1; 2), and the families (a; a+1), (1; b) and
    (a; 2a);
  - 2F1: (1,1;2) Log, (1,1;3/2), (1/2,1/2;3/2) ArcSin, (1/2,1;3/2) ArcTanh,
    (±1/2,1/2;1) K/E, and the families (a, a+½; ½), (a, a+½; 3/2),
    (a, −a; ½), (a, 1−a; 3/2), (1, a; 2) and (1, a; a+1) LerchPhi;
  - 1F2: (1/2; 3/2, 3/2) Shi, plus the families (a; a+½, 2a) BesselI² and
    (1; 2, b) (the integral of 0F1);
  - 0F3 Kelvin families: (1/2, (ν+1)/2, (ν+2)/2) =
    Γ(ν+1)(y/2)^(−ν)(I_ν(y) + J_ν(y))/2 with y = 4w^(1/4), which is
    `KelvinBer[0, 4 z^(1/4)]` at `−z` for ν = 0 (and
    `KelvinBei[0, 4 z^(1/4)]/(4 Sqrt[z])` for 0F3(; 1, 3/2, 3/2; −z), the same
    lattice); and (ν+1, (ν+1)/2, (ν+2)/2) = Γ(ν+1)²(x/2)^(−2ν) J_ν(x) I_ν(x)
    with x⁴ = −64w, which is `ber_0² + bei_0²` of 2√2 z^(1/4) at z for ν = 0;
  - 3F2: Clausen's (1, 1, 1; 3/2, 2) = ArcSin(√z)²/z.
  - Meijer G (the G walk, Roach '97): besides the Bessel, Struve, incomplete
    Gamma, Erfc, Si/Ci and Log/ArcTan/EllipticK entries,
    `G^{3,0}_{0,4}(w | ; 0, 0, 1/2 ; 1/2) = 4 ker_0(4 w^(1/4))` and
    `G^{3,0}_{0,4}(w | ; 0, 1/2, 1/2 ; 0) = −4 kei_0(4 w^(1/4))`, from the
    Mellin transforms of ker_0 and kei_0 and Gamma duplication.

  `` HypergeometricPFQ`Certificates[] `` lists the table, and
  `tests/test_hrep.c` checks every entry against its series.
- **Symbolic parameters are generic.** A definedness constant that is a
  nonzero rational function of the parameters counts as nonzero, as for any
  CAS identity in symbolic parameters.
- **Presentation.** The answer is collected over its special-function
  kernels. Each coefficient is factored over one denominator, and Simplify
  is applied when the coefficient carries a radical.

```mathematica
In[1]:= FunctionExpand[HypergeometricPFQ[{-3/2}, {-1/2, 1/2}, z]]
Out[1]= (1 + 2 z) Cosh[2 Sqrt[z]] + Sqrt[z] Sinh[2 Sqrt[z]] - 4 z^(3/2) SinhIntegral[2 Sqrt[z]]

In[2]:= FunctionExpand[Hypergeometric2F1[-3/2, -1/2, 2, z]]
Out[2]= (4 (1 + 14 z + z^2) EllipticE[z])/(15 Pi z) + (4 (-1 + z) (1 + 7 z) EllipticK[z])/(15 Pi z)

In[3]:= FunctionExpand[HypergeometricPFQ[{-1/2, 1/2, 1}, {3/2, 5/2}, z]]
Out[3]= (3 (-1 + z))/(16 z) + (3 (-1 + z^2) Log[1 - Sqrt[z]])/(32 z^(3/2)) - (3 (-1 + z^2) Log[1 + Sqrt[z]])/(32 z^(3/2)) - (3 PolyLog[2, -Sqrt[z]])/(8 Sqrt[z]) + (3 PolyLog[2, Sqrt[z]])/(8 Sqrt[z])

In[4]:= FunctionExpand[HypergeometricPFQ[{1, 3}, {2, 4, 5}, -z]]
Out[4]= (6 (-12 + z^2))/z^3 - (72 BesselJ[0, 2 Sqrt[z]])/z^3 + (144 BesselJ[1, 2 Sqrt[z]])/z^(7/2)

In[5]:= FunctionExpand[HypergeometricPFQ[{-11/4, -5/2, -5/2}, {-7/2, 1/4}, 1]]
Out[5]= (355 Gamma[1/4]^2)/(56 Sqrt[2] Sqrt[Pi])

In[6]:= FunctionExpand[HypergeometricPFQ[{}, {1/2, 1/2, 1}, -z]]
Out[6]= KelvinBer[0, 4 z^(1/4)]

In[7]:= FunctionExpand[HypergeometricPFQ[{1/2}, {-1/3, -1/2, -2/3}, z]]
Out[7]= E^(3 z^(1/3)) (1/3 - z^(1/3) + (3 z^(2/3))/2 - 3 z) + ((2/3 + z^(1/3) - (3 z^(2/3))/2 - 6 z) Cos[(3 Sqrt[3] z^(1/3))/2] + (Sqrt[3] (2 z^(1/3) + 3 z^(2/3)) Sin[(3 Sqrt[3] z^(1/3))/2])/2)/E^((3 z^(1/3))/2)

In[8]:= FunctionExpand[Hypergeometric2F1[1, 1, 2/3, 8/9]]
Out[8]= 9 + 3 Sqrt[3] Pi - 3 Log[3]
```

Roach 2008 (*Solving Integrals with the Quantum Computer Algebra System*),
p. 6: all eleven pFq values close, in Roach's printed form or an equal
closed form (2F1(1, 1; 2/3; 8/9) = 3^(3/2) π + 9 − 3 log 3 through the
complex-constant tidy). Of the twelve ISSAC '96 gallery entries, eleven close
through `FunctionExpand`, ten of them in Roach's printed form (the Kelvin
0F3(; 1, 5/2, 4; z), the 1F4(1/3; −2/3, −1/2, 1/2, 1; z) in I_0, J_0, I_1, J_1
of 4z^(1/4), and the duplication-formula 1F3(1/2; −1/3, −1/2, −2/3; z)
included). The one that does not is 2F3(−1/2, 1; 1/4, 1/2, 3/4; z), Roach's
erf·erfi formula: its lattice needs the integral of the duplication-formula
function, `∫ y² (e^(y²) erf y − e^(−y²) erfi y) dy`, which `Integrate`
declines (it is elementary in erf y · erfi y only after integration by
parts), and no contiguity walk crosses the lower 1/2.

## JacobiP

`JacobiP[n, a, b, x]` gives the Jacobi polynomial `P_n^(a,b)(x)`, orthogonal on
`[-1, 1]` with weight `(1-x)^a (1+x)^b`:

    P_n^(a,b)(x) = Σ_{l=0}^n (n+a+b+1)_l (a+l+1)_(n-l) / (l! (n-l)!) ((x-1)/2)^l
                 = Binomial[n+a, n] 2F1(-n, n+a+b+1; a+1; (1-x)/2)        (DLMF 18.5.7)

The sum form has no division by a parameter, so it is the polynomial for every
`a`, `b` (including a negative integer `a` in `[-n, -1]`).

- **Exact polynomial.** An integer `n ≥ 0` with exact rational `a`, `b` gives the
  expanded polynomial with exact rational coefficients (`O(n)` coefficient
  operations, then Horner in `(x-1)/2`); an exact rational `x` evaluates it
  exactly. A symbolic `a` or `b` expands the sum form (degree ≤ 40):
  `JacobiP[1, a, b, x] = (a - b)/2 + (a + b + 2) x/2`. `JacobiP[0, a, b, x] = 1`
  and `JacobiP[n, a, b, 1] = Binomial[n + a, n]`. `a = b = 0` is `LegendreP`.
- **Numerics.** An integer degree with an inexact argument runs the three-term
  recurrence (DLMF 18.9.2) -- in `double` (`sf_machine_jacobi_p[_c]`) at machine
  precision, the helper the packed / `NDArray` kernel and `Compile[]` share (exact
  rational parameters enter as their nearest doubles), else in MPFR, real or
  complex (complex `a`, `b` too) -- never the expanded polynomial, whose
  coefficients cancel catastrophically near `x = 1`. A non-integer degree, or the
  degenerate recurrence `a + b ∈ {-2, -3, …}`, evaluates
  `Binomial[n+a, n] Hypergeometric2F1[…]`.
  `JacobiP[300, 3, 1/2, 0.9998] = 100955.529716…` (1 µs),
  `JacobiP[2.5, 1, 2, 0.3] = -0.814842`.
- **Derivative.** `D[JacobiP[n, a, b, x], x] = (n + a + b + 1)/2 JacobiP[n - 1,
  a + 1, b + 1, x]` (DLMF 18.9.15); derivatives in `n`, `a`, `b` stay inert.
- **Surfaces.** Scalar, `List`, packed `List`, `NDArray` and `Compile[]` (scalar
  and rank-1) agree bit for bit (`REG_N` four-argument kernel).

Attributes: `Listable`, `NumericFunction`, `Protected`.

## RiemannR

`RiemannR[x]` is Riemann's prime-counting function
`R(x) = Sum_k MoebiusMu[k]/k LogIntegral[x^(1/k)]`.

- Attributes: `Listable`, `NumericFunction`, `Protected`.
- Exact arguments: `RiemannR[1] = 1`, `RiemannR[Infinity] = Infinity`; every
  other exact argument stays unevaluated (`RiemannR[2]`, `RiemannR[0]`), as in
  the reference. `N[RiemannR[2], 30]` evaluates.
- Method: always the Gram series `1 + Sum_k L^k/(k k! Zeta[k + 1])` in
  `L = Log x` (principal Log, so `x < 0` and complex `x` are covered: it is
  entire in `L`). The Moebius sum converges only conditionally.
- Machine real `x >= 1/e`: `sf_riemannr_r`, the series in double-double with
  `Log x` itself in double-double (`frexp` + `log1p` + split `ln 2`): a plain
  `log` costs ~300 ulps at `x = 10^300`, because `dR/R ~ dLog x`. Within an ulp.
  Machine complex `x` with `|Log x| <= 3/2`: `sf_riemannr_c`. Both are the unary
  kernel the packed / `NDArray` map and `Compile[]` use (scalar == packed ==
  `NDArray` == `Compile`, bit for bit); outside those domains the kernel declines
  and the scalar path computes in MPFR at 53 bits and rounds.
- `0 < x < 1` cancels: the terms grow to `~e^|L|/|L|^(3/2)` before summing to
  `R(x)`, which is `O(1)` or smaller (`R(e^-u) ~ 21.9/u^3`). The MPFR path works
  at `p + guard + log2(largest term)` bits and is adaptive in the guard: an
  answer below its error budget is recomputed with enough extra bits for `p`
  correct bits relative to `|R|`. At the largest zero of R (SIAM 100-Digit
  Challenge D4, `x* = 1.8286...e-14828`, `Log x* = -34142.128...`) the terms reach
  `10^14820` and the sum is `10^-13`: about 49000 bits, 1.7 s per 60-digit value.
- `1/Zeta[k + 1]` is never a per-term full-precision zeta: `k < k0` (where the
  series is short and cheap) uses `mpfr_zeta_ui` at the bits term `k` needs; for
  `k >= k0` the Dirichlet series `1/Zeta[s] = Sum MoebiusMu[n] n^-s` is run as one
  recurrence per squarefree `n <= 1024`, each at the precision its own largest
  term needs, multiplying by the SAME rounded `L` the main recurrence uses (a
  separately rounded `L/n` would make each stream the tail of a series in another
  variable, whose truncated head no longer cancels).
- Derivative: `D[RiemannR[x], x]` is the inert `Derivative[1][RiemannR][x]` (as
  in the reference), but it has a numeric value at an inexact point --
  `R'(x) = 1/(x L) Sum L^k/(k! Zeta[k+1])`, from the same pass -- through a
  `Derivative` SubValue to `Internal`RiemannRPrime`. `FindRoot` sees that the
  derivative numericalises and takes Newton steps on it instead of finite
  differences; a one-entry memo makes `R` and `R'` at the same expensive point
  one pass.
  `FindRoot[RiemannR[Exp[t]] == 0, {t, -34000}, WorkingPrecision -> 60]` gives
  `t* = -34142.1281846064649926300460679835405356170611068821899648398` in ~17 s
  (6 Newton steps).

## RiemannSiegelTheta and RiemannSiegelZ

`RiemannSiegelTheta[t]` is the Riemann–Siegel theta function
`θ(t) = Im LogGamma(1/4 + i t/2) - t Log(π)/2` (for complex `t`, the analytic
continuation `(LogGamma(1/4 + i t/2) - LogGamma(1/4 - i t/2))/(2 i) - t Log(π)/2`);
`RiemannSiegelZ[t]` is Hardy's function `Z(t) = Exp(i θ(t)) Zeta(1/2 + i t)`,
real for real `t`, with `|Z(t)| = |Zeta(1/2 + i t)|`, so its real zeros are the
zeros of `Zeta` on the critical line (the 10th at `49.7738324776723021819…`).
`θ` is odd and `Z` even.

- `RiemannSiegelTheta[0] = 0`, `RiemannSiegelZ[0] = Zeta[1/2]`; other exact
  arguments stay symbolic.
- Machine, arbitrary-precision and complex arguments evaluate through arb's
  `acb_dirichlet_hardy_theta` / `acb_dirichlet_hardy_z` behind the bridge's
  precision ladder (so `Z` next to a zero is a tight ball around 0, and a large
  `t` climbs precision to keep its digits). A real argument answers a real
  number. The machine helper (`sf_machine_riemann_siegel_z[_c]`, `_theta[_c]`,
  arb at a 53-bit target rounded to `double`) is shared by the scalar path, the
  packed / `NDArray` kernels and `Compile[]`, so the surfaces agree bit for bit.
  Without FLINT the defining expressions over `LogGamma` and `Zeta` are evaluated
  numerically (and the packed / `Compile[]` kernels decline to the interpreter).
- `FindRoot[z[t], {t, 4977/100}]` with `z[t_?NumericQ] := RiemannSiegelZ[t]`
  gives the 10th zeta zero at machine precision and to 50 digits.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```
In[1]:= RiemannSiegelZ[{1., 20.5}]
Out[1]= {-0.736305, 0.599329}

In[2]:= RiemannSiegelTheta[10.]
Out[2]= -3.06707
```

## ZetaZero

`ZetaZero[k]` is the `k`-th nontrivial zero of the Riemann zeta function above
the real axis, `1/2 + I t_k` (`t_k` the `k`-th positive zero of
`RiemannSiegelZ`); `ZetaZero[-k]` is its conjugate.

```
In[1]:= N[ZetaZero[1]]
Out[1]= 0.5 + 14.1347 I
In[2]:= N[Im[ZetaZero[1000]], 30]
Out[2]= 1419.42248094599568646598903808
In[3]:= N[ZetaZero[10^6]]
Out[3]= 0.5 + 600270. I
```

- An exact index stays symbolic; `N` gives the value (the head is `NHoldFirst`,
  and `N` reaches it through `numeric_register_held_value`). A number that is
  not a nonzero integer emits `ZetaZero::intnz` and stays unevaluated.
- arb's `acb_dirichlet_hardy_z_zero`: the zero is isolated rigorously (Gram /
  Rosser blocks, Turing's method), so the `k`-th zero is certified to be the
  `k`-th, then refined at the requested precision (machine: 53 bits, lowered).
  `ZetaZero[10^6]` in ~75 ms. Without FLINT it stays unevaluated.
- No packed / `NDArray` / `Compile` surface: the argument is an exact index.

Attributes: `Listable`, `NHoldFirst`, `NumericFunction`, `Protected`.

## AiryAiZero, AiryAiPrimeZero, AiryBiZero and AiryBiPrimeZero

`AiryAiZero[k]`, `AiryAiPrimeZero[k]`, `AiryBiZero[k]` and `AiryBiPrimeZero[k]`
are the `k`-th zeros (`k = 1, 2, ...`) of `AiryAi`, `AiryAiPrime`, `AiryBi` and
`AiryBiPrime`; all are negative reals.

```
In[1]:= N[AiryAiZero[1]]
Out[1]= -2.33811
In[2]:= N[AiryAiZero[500], 30]
Out[2]= -177.009251046479206825466153442
```

- An exact index stays symbolic and `N` gives the value (`NHoldFirst`, as for
  `ZetaZero`); a non-integer or non-positive index stays unevaluated.
- arb's `arb_hypgeom_airy_zero` (the asymptotic expansion of the zero, DLMF
  9.9.6-9.9.9, refined by Newton with a certified enclosure), at the requested
  precision. Without FLINT they stay unevaluated.
- No packed / `NDArray` / `Compile` surface: the argument is an exact index.

Attributes: `Listable`, `NHoldFirst`, `NumericFunction`, `Protected`.

## LegendreP

`LegendreP[n, x]` gives the Legendre polynomial / function `P_n(x)`;
`LegendreP[n, m, x]` gives the associated Legendre function `P_n^m(x)` (type 1);
`LegendreP[n, m, a, x]` gives the Legendre function of type `a` (`a ∈ {1, 2, 3}`,
default `1`).

- **Exact polynomial.** An exact integer order `n` produces the explicit
  degree-`|n'|` polynomial in `x` with exact rational coefficients, built from
  the three-term recurrence `k P_k = (2k-1) x P_{k-1} - (k-1) P_{k-2}`. Negative
  orders use `P_{-1-n} = P_n`. For example `LegendreP[10, x] = -63/256 +
  3465/256 x² - 15015/128 x⁴ + 45045/128 x⁶ - 109395/256 x⁸ + 46189/256 x¹⁰`.
  `LegendreP[n, 1] = 1` for any order, and `LegendreP[2, 2] = 11/2`.
- **Arbitrary-precision argument, integer degree.** An exact integer order with
  an inexact argument (beyond the machine helper below: an MPFR `x`, or a degree
  above 1000) runs the three-term recurrence in MPFR, real or complex, at the
  argument's precision plus guard bits -- never the expanded polynomial, whose
  coefficients reach `~4^n/√n` while `|P_n| ≤ 1` on `[-1, 1]`:
  `LegendreP[1000, 0.999997`50]` substituted into the polynomial gave
  `2.4×10^329` in 2 s; the recurrence gives `-0.02359077123132523175718700766…`
  (all 50 digits) in 0.7 ms.
- **Numerics.** A non-integer order with an inexact argument evaluates through
  the Gauss series `P_n(x) = 2F1(-n, n+1; 1; (1-x)/2)` at machine or
  arbitrary (MPFR) precision, real or complex (requires `|(1-x)/2| < 1`).
  Examples: `LegendreP[2.5, 2] = 9.58312`, `LegendreP[3/2 + I, 1.5 - I] =
  5.20466 + 0.299479 I`, and `N[LegendreP[3/2, 2], 50]` is accurate to 50
  digits with precision tracking the input. An exact non-integer order with an
  exact argument (e.g. `LegendreP[3/2, 2]`) stays symbolic and numericalizes
  only under `N`.
- **Associated functions (type 1).** For integer `n` and integer `m ≥ 0` the
  default type uses the Rodrigues derivative form `(-1)^m (1-x²)^(m/2)
  d^m/dx^m P_n(x)`; it is `0` when `m > |n'|`. For example `LegendreP[2, 1, x] =
  -3 x Sqrt[1 - x²]` and `LegendreP[2, 2, 2] = -9`.
- **Types 2 and 3.** These multiply the regularized Gauss polynomial
  `C(x) = 2F1Reg(-n, n+1, 1-m, (1-x)/2)` by `(1+x)^(m/2) (1-x)^(-m/2)` (type 2)
  or `(1+x)^(m/2) (-1+x)^(-m/2)` (type 3): `LegendreP[2, 1, 2, z]` and
  `LegendreP[2, 1, 3, z]` give the corresponding branch forms.
- **Listable.** `LegendreP[{1, 2, 3}, x] = {x, -1/2 + 3/2 x², -3/2 x + 5/2 x³}`.
- **Machine numerics (integer degree).** `LegendreP[n, x]` with an
  integer-valued `n` (exact, or a machine real like `3.`; negative `n` via
  `P_{-n-1} = P_n`) and a machine real or complex `x` is evaluated by the
  three-term recurrence in double / double-complex
  (`sf_machine_legendre_p_c`, `src/special_functions/sf_machine_plog.c`), the
  same helper the `NDArray` kernel and `Compile[]` call, so all five surfaces
  agree bit-for-bit; the exact polynomial is not built for a machine argument.
  (`LegendreP[n, 1]` with an exact 1 stays the exact `1`.)
- **Performance over machine arrays.** When the argument is a packed-List or
  `NDArray` of machine reals/complexes, the associated and type-2/3 polynomials
  are evaluated by Horner in machine `double` rather than as a sum of monomials:
  the sum form computed each `x^k` as its own full-buffer `pow()`, so a
  degree-`n` associated `LegendreP` did `O(n)` `pow()` passes and grew
  pathological with `n` (degree 30 over 10⁶ points did not finish). Horner is
  `O(n)` fused multiply-adds, uses the same coefficients (result unchanged to
  rounding), and runs the three-argument form faster than SciPy's `lpmv` at
  every degree tested (≈1.4–1.9×; degree 30 over 10⁵ points 10.6 ms vs 19.8 ms).
  `LegendreQ`'s associated form shares the same path. The symbolic output of
  every form is unchanged (the Horner route is taken only for a numeric buffer).

- **Derivative.** `D[LegendreP[n, x], x] = ((1+n) P_{n+1}(x) + (-1-n) x P_n(x)) /
  (x² - 1)` for symbolic `n`; an integer order differentiates the explicit
  polynomial.

Non-integer associated/type forms, negative `m`, symbolic `Series` /
`SeriesCoefficient`, and analytic continuation of the numeric series for
`|(1-x)/2| ≥ 1` are left symbolic.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= LegendreP[3, x]
Out[1]= -3/2 x + 5/2 x^3

In[2]:= LegendreP[10, 2, x]
Out[2]= (1 - x^2) (3465/128 - 45045/32 x^2 + 675675/64 x^4 - 765765/32 x^6 + 2078505/128 x^8)

In[3]:= LegendreP[5, 0.3 + 0.2 I]
Out[3]= 0.594236 - 0.0178925*I

In[4]:= LegendreP[2.5, 2]
Out[4]= 9.58312

In[5]:= D[LegendreP[n, x], x]
Out[5]= ((1 + n) LegendreP[1 + n, x] + (-1 - n) x LegendreP[n, x])/(-1 + x^2)
```

- **General degree and order** (v0.477): a numeric call with an inexact
  argument that the series above do not reach -- a non-integer `n` or `m`, or a
  point outside their disks -- evaluates through arb's
  `acb_hypgeom_legendre_p/q`. Types 1 and 2 are both the Ferrers form (arb type
  0); type 3 is arb type 1. The type label may arrive as `2.` (as `N` leaves it).
  Checked against the reference inside and outside (-1, 1) and at complex `x`.

## LegendreQ

`LegendreQ[n, x]` gives the Legendre function of the second kind `Q_n(x)`;
`LegendreQ[n, m, x]` gives the associated Legendre function `Q_n^m(x)` (type 1);
`LegendreQ[n, m, a, x]` gives the Legendre function of type `a` (`a ∈ {1, 2, 3}`,
default `1`).

- **Exact closed form.** Unlike `P_n`, `Q_n` carries a logarithm even for
  integer order. An exact integer order `n ≥ 0` produces
  `Q_n(x) = P_n(x) L(x) + v_n(x)` with `L(x) = (1/2)(Log[1+x] - Log[1-x])`, where
  the polynomial part `v_n` satisfies the same three-term recurrence as `P_n`
  but seeded `v_0 = 0`, `v_1 = -1`. For example
  `LegendreQ[0, x] = (1/2)(Log[1+x] - Log[1-x])`,
  `LegendreQ[1, x] = -1 + x L(x)`, and
  `LegendreQ[5, x] = -8/15 + 49/8 x² - 63/8 x⁴ + L(x) (15/8 x - 35/4 x³ +
  63/8 x⁵)`. Negative integer orders are singular (`Q_n` diverges) and stay
  symbolic.
- **Numerics.** A non-integer order with an inexact argument on the cut
  (`|x| < 1`) evaluates through the two origin Frobenius series
  `Q_v(x) = Q_v(0) 2F1(-v/2,(v+1)/2;1/2;x²) + Q_v'(0) x 2F1((1-v)/2,(v+2)/2;3/2;x²)`
  at machine or arbitrary (MPFR) precision, real or complex. Examples:
  `LegendreQ[1/3, 0.5] = -0.0399533`, `LegendreQ[1/2, 0.5] = -0.265596`,
  `LegendreQ[1/3 - I, 0.5 + 0.45 I] = -0.346939 + 2.52382 I`, and
  `N[LegendreQ[3/2, 1/2], 50]` is accurate to 50 digits with precision tracking
  the input. An exact non-integer order at an exact non-zero argument (e.g.
  `LegendreQ[3/2, 2]`) stays symbolic and numericalizes only under `N`; at an
  exact zero it gives the special value
  `Q_v(0) = -(√π/2) Sin[v π/2] Gamma[(v+1)/2] / Gamma[v/2+1]`.
- **Machine numerics (integer degree).** `LegendreQ[n, x]` with an integer
  `n ≥ 0` and a machine real or complex `x` is evaluated in double /
  double-complex (`sf_machine_legendre_q`, `sf_machine_legendre_q_c` in
  `src/special_functions/sf_machine_plog.c`, shared with the `NDArray` kernel
  and `Compile[]`) by Christoffel's form `Q_n = P_n Q_0 - W_{n-1}`, `P_n` and
  `W_{n-1}` both run up the three-term recurrence (`W_{-1} = 0`, `W_0 = 1`) —
  running the recurrence on `Q` itself seeds it with `Q_1 = x Q_0 - 1`, which
  cancels near `|x| ≈ 0.8`. `Q_0` is the closed form's `L(x)` evaluated without
  cancellation (`Re = log1p(4x/((1-x)² + y²))/4`, `Im = (arg(1+x) - arg(1-x))/2`,
  Kahan 1987) and with the closed form's branch: on the real axis beyond `±1`,
  `Log[1 - x] = log(x - 1) + iπ` (x > 1), so `LegendreQ[2, 1.5] =
  0.063567 - 4.51604 I`. `x = ±1` declines to the exact path (`Indeterminate`).
- **Associated functions (type 1).** For integer `n ≥ 0` and integer `m ≥ 0`,
  `(-1)^m (1-x²)^(m/2) d^m/dx^m Q_n(x)`. Because `Q_n` is not a polynomial, these
  carry `ArcTanh`/`Log` terms, e.g. `LegendreQ[3, 1, x]` and `LegendreQ[2, 2,
  0.5] = 4.06927`. `m = 0` reduces to `Q_n`.
- **Types 2 and 3.** These apply the branch-split prefactors `(1+x)^(m/2)
  (1-x)^(-m/2)` (type 2) or `(1+x)^(m/2) (-1+x)^(-m/2)` (type 3) to the shared
  core `(-1)^m (1-x)^m d^m/dx^m Q_n(x)` — the identical construction used for
  `LegendreP` types 2/3 — so they equal type 1 on the cut but carry the branch
  structure of type `a`.
- **Derivative and Series.** `D[LegendreQ[n, x], x] = ((-1-n) x Q_n + (1+n)
  Q_{n+1}) / (x²-1)`; higher integer-order derivatives follow by repetition.
  `Series[LegendreQ[v, x], {x, 0, k}]` gives the origin expansion.
- **Listable.** `LegendreQ[{0, 1, 2}, x]` threads element-wise.

Numeric analytic continuation for `|x| ≥ 1`, non-integer / complex associated
forms, negative-integer order, the asymptotic (`x → ∞`) series, and the
symbolic-order (`{x, l}` with symbolic `l`) derivative are left symbolic.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= LegendreQ[2, x]
Out[1]= -3/2 x + (1/2 Log[1 + x] - 1/2 Log[1 - x]) (-1/2 + 3/2 x^2)

In[2]:= LegendreQ[2, 0.5]
Out[2]= -0.818663

In[3]:= LegendreQ[3, 2.]
Out[3]= 0.00487112 - 26.7035*I

In[4]:= LegendreQ[1/2, 0]
Out[4]= -(Gamma[3/4] Sqrt[Pi])/(2 Sqrt[2] Gamma[5/4])
```

- **General degree and order** (v0.477): a numeric call with an inexact
  argument that the series above do not reach -- a non-integer `n` or `m`, or a
  point outside their disks -- evaluates through arb's
  `acb_hypgeom_legendre_p/q`. Types 1 and 2 are both the Ferrers form (arb type
  0); type 3 is arb type 1. The type label may arrive as `2.` (as `N` leaves it).
  Checked against the reference inside and outside (-1, 1) and at complex `x`.

## EllipticK, EllipticF, EllipticE and EllipticPi

The Legendre elliptic integrals. **The parameter argument is `m = k²`, not the
modulus `k`** — the reference convention. This is the one place a reader
of these signatures goes silently wrong, so each is written out as its integral:

- `EllipticK[m]` — ∫₀^(π/2) dt / √(1 − m sin²t), complete, first kind.
- `EllipticF[φ, m]` — ∫₀^φ dt / √(1 − m sin²t), incomplete, first kind.
- `EllipticE[m]` — ∫₀^(π/2) √(1 − m sin²t) dt, complete, second kind.
- `EllipticE[φ, m]` — ∫₀^φ √(1 − m sin²t) dt, incomplete, second kind.
- `EllipticPi[n, m]` — ∫₀^(π/2) dt / ((1 − n sin²t) √(1 − m sin²t)), complete,
  third kind.
- `EllipticPi[n, φ, m]` — the same with upper limit φ, incomplete, third kind.

`EllipticE` and `EllipticPi` are arity-overloaded exactly as in the reference: one
argument is complete for `E`, two arguments are complete for `Π`. There is no
`EllipticK[φ, m]` — that spelling is `EllipticF`.

**Attributes**: `Listable`, `NumericFunction`, `Protected`.

**Features**:
- Exact reductions: `EllipticK[0] = π/2`, `EllipticK[1] = ComplexInfinity`,
  `EllipticE[0] = π/2`, `EllipticE[1] = 1`; `EllipticF[0, m] = 0`,
  `EllipticF[φ, 0] = φ`, `EllipticF[π/2, m] = EllipticK[m]`;
  `EllipticE[π/2, m] = EllipticE[m]`;
  `EllipticPi[0, m] = EllipticK[m]`, `EllipticPi[0, φ, m] = EllipticF[φ, m]`,
  `EllipticPi[n, 0, m] = 0`, `EllipticPi[n, π/2, m] = EllipticPi[n, m]`.
- Exact non-special arguments stay symbolic (`EllipticF[1/3, 1/2]`); a numeric
  value follows from an inexact argument or from `N`.
- **Numeric evaluation** routes through FLINT/Arb's `acb_elliptic_*`, which is
  rigorous, arbitrary-precision, defined on the whole complex plane, and already
  uses the parameter convention and branch placement above. That buys the three
  things a hand-rolled kernel most easily gets wrong:
  - the **quasi-periodic extension** off the principal strip
    (`F(φ + kπ | m) = F(φ | m) + 2k K(m)`, and likewise for `E` and `Π`), so φ of
    any size is handled;
  - **complex φ**, which `EllipticF[ArcSin[z], m]` produces routinely as soon as
    `|z| > 1` — the normal case for an elliptic pencil written in `ArcSin` form;
  - the **value past the pole** for `EllipticPi` with `n > 1`, where the path
    crosses `sin²t = 1/n`. It is not a real principal value — `Π(3/2 | 1/2)` is
    `−0.456720313453 − 2.72069904635 i`, mpmath and Arb agreeing — so the machine
    kernel declines there and Arb answers.
- **Complete integrals at machine precision** (`src/special_functions/sf_machine_ellk.c`).
  One kernel per (function, real | complex), called by the scalar builtin, the
  packed/NDArray element-wise path and `Compile[]` alike, so all five surfaces
  (scalar, List, packed, NDArray, Compile) agree bit for bit.
  - *Real `m < 1`*: Cephes' form (Moshier's `ellpk`/`ellpe`; the structure of
    A&S 17.3.34/17.3.36) `K = P(t) − log(t) Q(t)`, `E = P(t) − t log(t) Q(t)`,
    `t = 1 − m`, with degree-10 `P`, `Q` **re-fitted** for Mathilda by
    `tools/sf_accuracy/fit_ellk.py` (linear relative-error minimax by Lawson
    IRLS in mpmath, points down to `t = 10⁻³⁰⁰`; fit error ≤ 1.2·10⁻¹⁸).
    `m < 0` by `K(m) = K(m/(m−1))/√(1−m)`, `E(m) = √(1−m) E(m/(m−1))`
    (DLMF 19.7.5), i.e. the same polynomials at `t' = 1/t`. Accuracy vs mpmath:
    max 1.7 ulp (`K`), 1.3 ulp (`E`) on `[0, 0.95]`; ≤ 2.4 ulp on `[−50, 1)`.
  - *Complex `m`*: the AGM, `K = π/(2 M(1, √(1−m)))`,
    `E = K (1 − Σ 2ⁿ⁻¹ cₙ²)` (DLMF 19.8.5–6), in raw `double` re/im with the
    "right choice" of square root at each step (Cox 1984), which gives the
    principal branch with the cut on `[1, ∞)`; on the cut the value is the limit
    from below (`EllipticK[2.] = 1.31103 − 1.31103 I`, as the reference and
    mpmath). Max 3.4 ulp (`K`), 3.9 ulp (`E`) on `|Re m|, |Im m| ≤ 0.9`; ≤ 4.7
    ulp on `|Re m|, |Im m| ≤ 20`.
  - A machine **complex** argument, and a machine **real** `m > 1`, now return
    a machine `Complex`; they used to go through Arb and come back as 53-bit
    `EXPR_MPFR` numbers. The pole `K(1)` and `|m| > 10³⁰⁰` decline to the
    exact/Arb path.
  - Cost (net, per call): `EllipticK[0.6]` 0.19 → 0.12 µs, `EllipticE[0.6]`
    0.24 → 0.12 µs; complex 6.7 → 0.37 µs (`K`) and 19.9 → 0.36 µs (`E`).
    10⁶-element buffer: `K` 12.5 → 4.0 ms, `E` 18.5 → 3.3 ms (SciPy 16.7 / 16.5).
- **Machine kernels (incomplete).** `EllipticE[φ, m]` and `EllipticF`
  carry `double` kernels built on Carlson's symmetric forms `R_F` and `R_D`, used
  by the packed/NDArray element-wise paths. They cover the real principal domain
  and **decline** outside it (`m > 1`, or `1 − m sin²φ < 0`), which abandons the
  buffer so the List path answers exactly through Arb — slower, never wrong.
  Measured against mpmath on identical machine inputs, the buffer now agrees with
  the scalar path to 3 ulp (`K`, `F`) and 6 ulp (`E`); `R_F` and `R_D` share one
  duplication loop, since every caller that wants `R_D` wants `R_F` at the same
  arguments and the two recurrences walk an identical sequence.
  `EllipticPi` carries one too, on Carlson's `R_J` (and `R_C`, which `R_J`
  needs): the **complete** `EllipticPi[n, m]` is a binary (`REG_B`) kernel at
  41 ns/element and 2 ulp, and the **incomplete** `EllipticPi[n, φ, m]` is an
  n-ary (`REG_N`) one at 1 ulp, carrying the quasi-period
  `Π(n; φ+kπ|m) = Π(n; φ|m) + 2k Π(n|m)`. The element-wise NDArray layer tops
  out at arity 2, so only the complete form rides that path; `Compile[]` lowers
  both, the three-argument shape through `OP_KERNN`. Like the others the kernel
  **declines** outside the real principal domain: `n > 1` is not a real
  principal value at all (`Π(3/2 | 1/2)` is `−0.45672 − 2.72070 i`), so the
  honest machine answer there is to hand back to Arb.
  Because `packed_aware` is a property of the symbol and not of one arity,
  registering the two-argument kernel also stopped the transparency gate
  materialising packed `List`s for the three-argument form — which is why the
  `ndarray_delist_and_reeval` fallback in all three argument positions had to
  land first, and must stay.
- **Exact values.** `EllipticK[-1] = Γ(1/4)²/(4√(2π))` and
  `EllipticK[1/2] = 8π^(3/2)/Γ(-1/4)²` (the lemniscatic singular values; `E` has
  no closed form at either and keeps none). `EllipticPi[n, 0] = π/(2√(1-n))`,
  `EllipticPi[n, φ, 0] = ArcTanh[√(n−1) tan φ]/√(n−1)`. At infinity:
  `EllipticK[∞] = 0`, `EllipticE[∞] = EllipticE[ComplexInfinity] =
  ComplexInfinity`, `EllipticF[φ, ∞] = 0`,
  `EllipticPi[∞, m] = EllipticPi[n, ∞] = 0`. All three incomplete forms are
  **odd in the amplitude** — `EllipticF[−φ, m] = −EllipticF[φ, m]`, likewise `E`
  and `Π` — because each integrand is even in `t`; the fold uses the same
  superficial-negativity test the trig heads use, so `−2x` folds and `−x−y` does
  not.

- **`Series`.** The complete integrals expand at `m = 0` from a dedicated
  kernel: `K = (π/2) Σ aₖ mᵏ` and `E = (π/2) Σ aₖ mᵏ/(1−2k)` with `a₀ = 1`,
  `aₖ = aₖ₋₁((2k−1)/(2k))²`. It has to be dedicated — the generic
  Taylor-via-`D` path evaluates `d/dm K` at `m = 0`, which is `(π/2 − π/2)/0`
  → `Indeterminate`, and the expansion is then abandoned, so
  `Series[EllipticK[m], {m, 0, 1}]` did not previously emit even the leading
  `π/2`. Composition through an inner series works; an expansion about a regular
  point still takes Taylor-via-`D`. The incomplete forms expand in `φ` by the
  generic path, as before.

- **`Interval`.** The complete integrals thread by certified monotonicity —
  `K` increasing and `E` decreasing, both only below `m = 1`, beyond which the
  value is complex. The amplitude slot of `EllipticF`, `EllipticE[φ,m]` and
  `EllipticPi` threads through the general derivative certifier. The `m` slot
  declines: its derivative reproduces the head, so the certifier can never
  bottom out there.

- **Machine numbers.** An inexact machine scalar is answered by the `double`
  kernel and comes back as a machine `Real`, not a 53-bit `EXPR_MPFR` — so
  `Precision[EllipticK[0.5]]` is `MachinePrecision`, `MachineNumberQ` is `True`,
  and an elliptic-produced list packs for its consumers. An `EXPR_MPFR`
  argument still routes to Arb, so `SetPrecision[m, 16]` is never silently
  answered in `double`.

- **The `m = 1` continuation.** `E(φ | 1) = ∫₀^φ |cos t| dt`, which is `Sin[φ]`
  only on the principal strip `|φ| ≤ π/2`; off it the value continues as
  `Sin[φ − kπ] + 2k`. The reduction therefore fires only where it is *provably*
  valid — decided by evaluating `Abs[φ] ≤ π/2`, so an exact rational, a `Real`,
  an `MPFR` and a symbol carrying assumptions all get an answer. A symbolic
  amplitude is undecidable and stays symbolic rather than becoming wrong:
  `EllipticE[2, 1]` is left alone, and `N[EllipticE[2, 1]]` is `1.0907025731743`
  (`= 2 − Sin[2]`), agreeing with its neighbour at `m = 1 − 10⁻²⁵`.

- **Inexactness and poles.** The numeric path runs *before* the remaining exact
  reductions, so an inexact argument gets an inexact answer — `EllipticK[0.]` is
  `1.5707963267948966`, not the exact `π/2`, and `EllipticPi[0., 1/2]` is
  `1.8540746773013719`, not the symbolic `EllipticK[1/2]`. The exact reductions
  remain underneath as the fallback, which is what still reduces
  `EllipticE[π/2, 0.5]` (Arb reads numbers, not a `Times`). Poles are checked
  *ahead* of the numeric path, in every spelling of the argument, because Arb
  returns a non-finite ball there and the bridge renders that as "unevaluated":
  `EllipticK[1]`, `EllipticK[1.]` and `EllipticK[SetPrecision[1, 30]]` all give
  `ComplexInfinity`, and so does the complete `EllipticPi[1, m]` — whose
  divergence is real (`Π(1−ε | 1/2)` grows as `1/√ε`). The *incomplete* third
  kind has no pole at `n = 1`: `EllipticPi[1, 1, 1/2]` is `1.73199154202`.

- **Derivatives.** The φ-derivatives are the integrands, which is what makes a
  numeric verification of an antiderivative built from these kernels close:
  `D[EllipticF[φ, m], φ] = 1/√(1 − m sin²φ)`,
  `D[EllipticE[φ, m], φ] = √(1 − m sin²φ)`,
  `D[EllipticPi[n, φ, m], φ] = 1/((1 − n sin²φ) √(1 − m sin²φ))`. Also
  `D[EllipticK[m], m] = (E(m) − (1−m) K(m)) / (2m(1−m))`,
  `D[EllipticE[m], m] = (E(m) − K(m)) / (2m)` and
  `D[EllipticE[φ, m], m] = (E(φ,m) − F(φ,m)) / (2m)`.
  `D[EllipticF[φ, m], m] = −E(φ,m)/(2(m−1)m) − F(φ,m)/(2m) +
  sin(2φ)/(4(m−1)√(1−m sin²φ))`, and for the **complete** third kind
  `D[EllipticPi[n, m], n] = (nE + (m−n)K + (n²−m)Π)/(2(m−n)(n−1)n)` and
  `D[EllipticPi[n, m], m] = (E/(m−1) + Π)/(2(n−m))`. Each of those three was an
  inert placeholder until it had been checked against a central difference at 30
  digits.

  The `n`- and `m`-derivatives of the **incomplete** `EllipticPi` remain
  deliberately inert, and not for want of trying: a least-squares fit over the
  natural candidate basis (E, F, Π over `(n−m)`, `(m−1)`, `(n−1)`, plus the
  `sin(2φ)/√(…)` boundary term) does not recover them — residual 1.5 relative,
  with no simple rational coefficients — so the closed form involves terms that
  basis does not span. An inert derivative is honest where a guessed one corrupts
  every caller silently; `BesselJ` treats its order the same way. One visible
  consequence: `Interval[]` threading for the incomplete `EllipticPi` in its `m`
  slot stays symbolic, since the certifier differentiates to get a sign.
- Wrong arity emits `EllipticK::argx` / `EllipticF::argrx` (fixed arity) or
  `EllipticE::argt` / `EllipticPi::argt` (the overloaded pair) and stays
  unevaluated.
- Without FLINT (`USE_FLINT` undefined) the numeric path falls back to the
  machine kernels (Cephes-form / AGM for the complete integrals, Carlson for
  the incomplete ones) and otherwise declines, leaving the
  call symbolic.

```mathematica
In[1]:= EllipticK[0]
Out[1]= 1/2 Pi

In[2]:= N[EllipticK[1/2], 25]
Out[2]= 1.8540746773013719184338503

In[3]:= D[EllipticF[phi, m], phi]
Out[3]= 1/Sqrt[1 - m Sin[phi]^2]

In[4]:= N[EllipticPi[3/2, ArcSin[Sqrt[2] Sqrt[1/(1 + 7/5)]], 1/2], 20]
Out[4]= 0.987739699728522080208 - 2.7206990463513267759*I

In[5]:= EllipticK[2.]
Out[5]= 1.31103 - 1.31103*I

In[6]:= EllipticE[0.3 + 0.4 I]
Out[6]= 1.46251 - 0.175161*I
```

## Surd

`Surd[x, n]` is the real-valued n-th root of a real `x` (v0.477; new): `Sign[x]
Abs[x]^(1/n)` for odd `n`, `x^(1/n)` for even `n` and `x ≥ 0` — unlike
`x^(1/n)`, which is the principal complex root for `x < 0`.

- Attributes: `Listable`, `NumericFunction`, `Protected`.
- Exact `x`: `Surd[8, 3] = 2`, `Surd[-8, 3] = -2`, `Surd[-2, 3] = -2^(1/3)`,
  `Surd[Sqrt[2], 3] = 2^(1/6)`, `Surd[8, -3] = 1/2`; a symbolic `x` stays.
- Machine `x`: `sf_surd_r` (sqrt / cbrt / pow with a Newton step), shared by the
  scalar path, the packed / `NDArray` map and `Compile`.
- `N` keeps the integer index exact, as it keeps an integer `Power` exponent:
  `N[Surd[x, 3], 30]` is `Surd[x, 3]` (it was `Surd[x, 3.`30]`, an invalid
  call that emitted `Surd::int` at every evaluation -- thousands of times
  inside a `WorkingPrecision` `NIntegrate`, whose integrand is numericalized
  first).
- Messages, as in the reference: `Surd::noneg` (even root of a negative,
  `Indeterminate`), `Surd::indet` (`n = 0`, `Indeterminate`), `Surd::preal` (a
  non-real `x`, unevaluated), `Surd::int` (a non-integer `n`, unevaluated),
  `::argr` / `::argrx` for the argument count.
- Deliberate difference: an integer-valued machine `n` (`Surd[8., 3.]`) is
  accepted as that integer; the reference rejects it with `Surd::int`. The
  packed kernel receives every operand as a double, so the scalar path accepts it
  too, to keep the surfaces in agreement.

## CubeRoot

`CubeRoot[x]` is the real cube root of a real `x`, the same as `Surd[x, 3]`.

- Attributes: `Listable`, `NumericFunction`, `Protected`.
- Exact `x`: `CubeRoot[8] = 2`, `CubeRoot[-8] = -2`, `CubeRoot[-27/8] = -3/2`,
  `CubeRoot[-2] = -2^(1/3)`, `CubeRoot[16] = 2 2^(1/3)`.
- Symbolic `x`: `Surd[x, 3]` (as in the reference, `CubeRoot[x^3] = Surd[x^3, 3]`).
- Machine `x`: `sf_cuberoot_r` (= `sf_surd_r(x, 3)`, `cbrt`), registered as
  CubeRoot's own unary kernel, so the scalar path, the packed / `NDArray` map and
  `Compile[]` (scalar and rank-1) agree bit for bit. Arbitrary precision: through
  `Surd`, at the argument's precision.
- A non-real number emits `CubeRoot::preal` and stays unevaluated
  (`CubeRoot[I]`); a packed buffer with such an element goes to the List path.

## HaarWavelet, DaubechiesWavelet, WaveletPhi, WaveletPsi, WaveletFilterCoefficients

`WaveletPhi[wave, x]` / `WaveletPsi[wave, x]`: the scaling and wavelet function
of `HaarWavelet[]` or `DaubechiesWavelet[n]` (`n` vanishing moments, `2n` taps;
`DaubechiesWavelet[] = DaubechiesWavelet[2]`). Options `MaxRecursion ->
Automatic` (= 8) and `WorkingPrecision -> MachinePrecision`. src/wavelets.c.

- Haar: the exact step functions -- `phi = 1` on `[0, 1)`; `psi = 1` on
  `(0, 1/2)`, `-1` on `(1/2, 1)`, else 0 (`psi(0) = psi(1/2) = 0`), returned as
  Reals; a symbolic `x` gives the reference's `Piecewise`.
- Daubechies: the cascade's exact values on the dyadic grid of spacing
  `2^-MaxRecursion` -- `phi` at the integers from the two-scale eigenproblem
  (normalised `Sum phi(k) = 1`), each finer level from `phi(x) = 2 Sum h_k
  phi(2x - k)`, and `psi(x) = 2 Sum g_k phi(2x - k)`, `g_k = (-1)^k h_(1-k)`
  (support `[1 - n, n]`) -- LINEARLY interpolated between grid points, which is
  what the reference returns: `WaveletPhi[DaubechiesWavelet[2], 1.3] =
  0.34067377391167836263936642769158893255...` (40 digits agree). A symbolic `x`
  gives `Piecewise[{{InterpolatingFunction[..][x], 0 <= x <= 2n - 1}}, 0]` (an
  order-1 `Interpolation` of the grid), so it can serve as an `NDSolve` initial
  condition.
- The filter: spectral factorisation of `P(y) = Sum_(k<n) Binomial[n-1+k, k] y^k`
  -- its roots by Durand-Kerner in MPFR complex, `z + 1/z = 2 - 4y`, the root
  with `|z| < 1` kept, `h` the reversed, normalised (`Sum h = 1`) coefficients of
  `(1 + z)^n Prod (z - z_j)`. This reproduces the reference's ordering
  (`{0.3415, 0.5915, 0.1585, -0.0915}` for `n = 2`; `n = 3` to 25 digits).
  `WaveletFilterCoefficients[wave, "PrimalLowpass" | "PrimalHighpass"]` returns
  it as `{{k, h_k}, ...}` (highpass indices `2 - 2n .. 1`).
- Precision: the smaller of `WorkingPrecision` and the precision of `x`; an exact
  `x` takes the `WorkingPrecision`, so the default is a machine number
  (`WaveletPhi[w, 1/3]` is a Real), as in the reference.
- The grid is memoised per `(n, MaxRecursion, precision)`; `NIntegrate` and
  `NDSolve` sample `phi` thousands of times.
- Deliberate difference: a List `x` is mapped element-wise (the reference
  returns a `Piecewise` over the whole list).
- Attributes: `Protected` (all five heads). No packed / `NDArray` / `Compile`
  surface: the first argument is a wavelet object, not a number.

## KroneckerDelta

`KroneckerDelta[n1, n2, ...]` is 1 when every `ni` is equal and 0 when any two
are unequal (v0.477; new). `KroneckerDelta[]` is 1 and `KroneckerDelta[n]` is
`KroneckerDelta[n, 0]`.

- Attributes: `NumericFunction`, `Orderless`, `Protected` (not `Listable`).
- Identical arguments give 1 (`KroneckerDelta[x, x]`); numbers and numeric
  quantities are compared by `Equal` (`KroneckerDelta[1/2, 0.5] = 1`,
  `KroneckerDelta[Pi, 3] = 0`); one provably unequal pair settles the whole call
  (`KroneckerDelta[x, 1, 2] = 0`). A symbolic pair that is not identical is
  never guessed (`KroneckerDelta[x, y]` stays).

## BesselJZero

`BesselJZero[n, k]` is the k-th positive zero of `BesselJ[n, x]`.

- Attributes: `Listable`, `NumericFunction`, `Protected`.
- An exact call stays unevaluated (the Hadamard-product recogniser in `Product`
  relies on it). With an inexact argument — and so under `N`, which leaves
  `k = 2.` — it evaluates (v0.477; previously never) for a real order `n > -1`
  and a positive integer `k`.
- Machine: `sf_bessel_j_zero_r`, McMahon's expansion when `k > n + 20`, otherwise
  a sign-change scan of `J_n` from the origin in steps of 1/4 (consecutive zeros
  are more than 2 apart), then Newton on `J_n / J_n'`. Shared by the scalar path,
  the packed / `NDArray` map and `Compile`. Within 3 ulp of the reference.
- Arbitrary precision: the machine zero refined by Newton in arb with the
  precision doubling (`N[BesselJZero[0, 1], 30]`). A step uses
  `J_n' = J_{n-1} - (n/x) J_n` (two Bessel evaluations, not three) and one step
  per doubled precision from 96 bits plus one at the target, since the machine
  seed already carries ~50 correct bits: `N[BesselJZero[0, 1000], 30]` went from
  ~1 ms to ~0.3 ms.

## BesselYZero

`BesselYZero[n, k]` is the k-th positive zero of `BesselY[n, x]`, for a real
order `n >= 0`.

- Attributes: `Listable`, `NumericFunction`, `Protected`. Exact calls stay
  unevaluated; an inexact argument (and so `N`) evaluates, as for `BesselJZero`.
- The same machinery as `BesselJZero` (one implementation, `bz_zero_r`):
  McMahon's expansion with `beta = (k + n/2 - 3/4) Pi` when `k > n + 20`,
  otherwise a sign-change scan from `max(n, 0.05)` (no zero of `Y_n` lies below
  `n`), then Newton. Machine helper `sf_bessel_y_zero_r`, shared by the scalar
  path, the packed / `NDArray` map and `Compile`; arbitrary precision by Newton in
  arb (`flint_num_bessel_y_zero`).

```
In[1]:= N[BesselYZero[0, 1]]
Out[1]= 0.893577
In[2]:= N[BesselYZero[1, 2], 30]
Out[2]= 5.42968104079413513277200519085
```

## HeavisideTheta and DiracDelta

Defined by rules in `src/internal/distributions.m`. `HeavisideTheta[x]` is 1 for
positive and 0 for negative numeric `x`, and `HeavisideTheta[0] = 0` (the causal,
left-continuous convention DSolve's impulse responses rely on; the reference
leaves it unevaluated). `DiracDelta[x]` is 0 for a nonzero real numeric `x`.

The multivariate forms are products (v0.477): `HeavisideTheta[x1, x2, ...]` is 0
when any `xi` is negative (or 0) and 1 when all are positive, and drops its
positive arguments otherwise (`HeavisideTheta[1, x] = HeavisideTheta[x]`);
`DiracDelta[x1, x2, ...]` is 0 as soon as one `xi` is nonzero. Mixed exact and
inexact arguments (`HeavisideTheta[1., Sqrt[2]]`) decide by sign like any other.

**DiracDelta normal form** (v0.6140), per argument, as the reference writes it:
- **Sign.** An argument with more terms of coefficient exactly -1 than of exactly 1
  (numbers aside) is negated: `DiracDelta[1 - x1 - x2]` is `DiracDelta[-1 + x1 + x2]`
  and `DiracDelta[-Sin[x]]` is `DiracDelta[Sin[x]]`. `DiracDelta[a - b]`,
  `DiracDelta[-a + b]` and `DiracDelta[1 - 2 x]` stay, and so does an inexact
  argument (`DiracDelta[1.5 - x]`).
- **Scale.** A numeric factor comes out as 1/|c|: `DiracDelta[2 x]` is
  `DiracDelta[x]/2`, `DiracDelta[2.5 x]` is `0.4 DiracDelta[x]`,
  `DiracDelta[4 x + 6 y]` is `DiracDelta[2 x + 3 y]/2` (the rational gcd of an exact
  sum, constants included), and `DiracDelta[2 x, y]` is `DiracDelta[x, y]/2`.
- `HeavisideTheta` is not normalised (`HeavisideTheta[-x]` stays). FeynCalc's
  `FCFeynmanParameterJoin` writes its delta function in this form.
