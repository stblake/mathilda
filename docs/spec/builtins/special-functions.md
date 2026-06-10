# Special Functions

Higher transcendental functions: the gamma function `Gamma`, the error function `Erf` and its complement `Erfc`, the digamma/polygamma family `PolyGamma`, the log-gamma function `LogGamma`, the Pochhammer symbol (rising factorial) `Pochhammer`, the Riemann/Hurwitz zeta function `Zeta` (with the inert Stieltjes constants `StieltjesGamma`), the Bernoulli numbers and polynomials `BernoulliB`, the Euler numbers and polynomials `EulerE`, the polylogarithm `PolyLog`, and the hypergeometric family `Hypergeometric0F1`, `Hypergeometric1F1`, `Hypergeometric2F1`, and the generalized `HypergeometricPFQ`.

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
  - Machine-precision real → libm `tgamma`; machine-precision complex →
    Lanczos approximation, e.g. `Gamma[2.3 + I] = 0.719141 + 0.540614 I`.
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
  - Machine-precision real → libm `erf`, e.g. `Erf[0.95] = 0.820891`,
    `Erf[1.5] = 0.966105`.
  - Arbitrary precision (MPFR) real → `mpfr_erf`, output precision tracking the
    input, e.g.
    `N[Erf[3/2], 50] = 0.96610514647531072706697626164594785868141047925764`
    and `Erf[0.95`100]`.
  - **Complex** (machine *and* arbitrary precision) → the cancellation-aware
    Maclaurin series erf(z) = (2/√π) e^(−z²) Σ 2ⁿ z^(2n+1)/(1·3···(2n+1))
    (DLMF 7.6.2), evaluated in MPFR with `|z|²/ln2` guard bits so even
    machine-precision complex results carry full accuracy, e.g.
    `Erf[1.5 - I] = 1.0784 + 0.0279637 I`,
    `N[Erf[1/2 + I], 30] = 1.20484755831421800270211268210 + 1.02440088160844588172486045441 I`.
    A `double complex` series is the fallback for `USE_MPFR=0` builds.
- The two-argument form `Erf[z0, z1]` reduces to `erf(z1) − erf(z0)` only when
  both reduce to something concrete, e.g. `Erf[1.5, 2] = 0.0292171`,
  `Erf[-Infinity, Infinity] = 2`; exact/symbolic pairs such as `Erf[2, 3]` and
  `Erf[a, b]` stay unevaluated.
- Derivative: `D[Erf[z], z] = (2/Sqrt[Pi]) E^(−z²)` (chain rule applies), so
  the origin Taylor series follows from the generic `D`-based fallback, e.g.
  `Series[Erf[x], {x, 0, 5}]` begins `2/Sqrt[Pi] x − 2/(3 Sqrt[Pi]) x^3 + …`.
- All other arguments (symbolic `Erf[x]`, exact `Erf[2]`) stay unevaluated.

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
  - Machine-precision real → libm `erfc`, e.g. `Erfc[0.95] = 0.179109`,
    `Erfc[1.5] = 0.0338949`.
  - Arbitrary precision (MPFR) real → `mpfr_erfc`, which is cancellation-free
    even for large positive z (where `1 − erf(z)` would lose all significance),
    output precision tracking the input, e.g.
    `N[Erfc[3/2], 50] = 0.033894853524689272933023738354052141318589520742363`.
  - **Complex** (machine *and* arbitrary precision) → `1 − erf(z)`, with erf(z)
    from the cancellation-aware DLMF 7.6.2 series evaluated in MPFR; the
    complement is formed at working precision (with `|z|²/ln2` guard bits) before
    rounding, so even machine-precision complex results carry full accuracy, e.g.
    `Erfc[1.5 - I] = -0.0783992 - 0.0279637 I`. A `double complex` series is the
    fallback for `USE_MPFR=0` builds.
- Derivative: `D[Erfc[z], z] = -(2/Sqrt[Pi]) E^(−z²)` (chain rule applies), so
  the origin Taylor series follows from the generic `D`-based fallback, e.g.
  `Series[Erfc[x], {x, 0, 3}]` begins `1 − 2/Sqrt[Pi] x + …`.
- All other arguments (symbolic `Erfc[x]`, exact `Erfc[2]`) stay unevaluated.

## InverseErf

- `InverseErf[s]` — the inverse error function: the z solving s = erf(z).
- `InverseErf[z0, s]` — the inverse of the generalized error function
  `Erf[z0, z]`, i.e. `InverseErf[s + Erf[z0]]`.

**Attributes**: `Listable`, `NumericFunction`, `Protected`.

`InverseErf` is odd in s and (per Mathematica) produces explicit numerical
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
  - Machine-precision real → a Winitzki seed polished by Newton on libm `erf`,
    e.g. `InverseErf[0.6] = 0.595116`, `InverseErf[1/{2.,3.,4.,5.}] =
    {0.476936, 0.30457, 0.225312, 0.179143}`.
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
  - Machine-precision real → a Winitzki seed polished by Newton on libm `erfc`,
    e.g. `InverseErfc[0.6] = 0.370807`, `InverseErfc[1/{2.,3.,4.,5.}] =
    {0.476936, 0.68407, 0.81342, 0.906194}`.
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

## PolyGamma

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

- `PolyGamma[z]` — the digamma function ψ(z) = Γ′(z)/Γ(z). Always rewrites to
  the two-argument form `PolyGamma[0, z]`.
- `PolyGamma[n, z]` — the n-th polygamma ψ⁽ⁿ⁾(z) = dⁿ/dzⁿ ψ(z).

Behaviour:

- **Special values.** Non-positive integer arguments are poles:
  `PolyGamma[0]`, `PolyGamma[n, -k]` → `ComplexInfinity`. At infinity,
  `PolyGamma[0, Infinity]` → `Infinity` and `PolyGamma[n, Infinity]` → `0` for
  `n ≥ 1`. `PolyGamma[ComplexInfinity]` and `PolyGamma[n, Indeterminate]`
  → `Indeterminate`.
- **Exact at positive integers.** `PolyGamma[0, m] = H_{m-1} − EulerGamma`
  (a rational minus Euler's constant), e.g. `PolyGamma[5] = 25/12 - EulerGamma`,
  `PolyGamma[100]` a large exact rational minus `EulerGamma`. For **odd** order
  `n ≥ 1` the value closes via ζ(n+1) into a rational plus a rational multiple of
  `Pi^(n+1)`: `PolyGamma[1, 1] = Pi^2/6`, `PolyGamma[3, 1] = Pi^4/15`,
  `PolyGamma[3, 5] = -22369/3456 + Pi^4/15`. For **even** order `n ≥ 1` the value
  involves ζ(odd) and stays symbolic, e.g. `PolyGamma[2, 1]`, `PolyGamma[4, 1]`.
- **Negative order.** `PolyGamma[-1, z] = LogGamma[z]` (the log-gamma function;
  see below). Orders ≤ −2 stay unevaluated.
- **Numeric.** Inexact real arguments evaluate at machine precision (`Real`) or
  arbitrary precision (`MPFR`, precision tracked from the input):
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

### LogGamma

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
- **Numerics.** Machine real via `lgamma`; arbitrary-precision real via MPFR
  `lgamma`; machine complex via a Lanczos log-series; arbitrary-precision
  complex via a Stirling series with argument reduction (the continuous
  branch — its imaginary part grows past π where `Log[Gamma]` would wrap).
  Negative real arguments return the complex value `log|Γ(z)| - Pi Ceiling[-z] I`.
- **Derivative.** `D[LogGamma[z], z] = PolyGamma[0, z]`; higher derivatives
  raise the `PolyGamma` order. Produced by `PolyGamma[-1, z]`.

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
  Mathematica expresses through `PolyGamma`, are not yet implemented.)

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
  `2F1(1,1;2) -> -Log[1-z]/z`.
- Numeric evaluation at machine, arbitrary (MPFR), and complex precision by
  direct series summation, with output precision tracking the input. The
  series is summed only where it converges — `p <= q` (entire) and `p == q+1`
  with `|z| < 1`; outside that regime (and `p > q+1`) the call stays
  unevaluated (analytic continuation beyond the unit disk is not yet
  implemented).
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

## Hypergeometric0F1


`Hypergeometric0F1[b, z]` is the confluent hypergeometric function ₀F₁,
equal to `HypergeometricPFQ[{}, {b}, z]`.

```mathematica
In[1]:= Hypergeometric0F1[1/2, z]
Out[1]= Cosh[2 Sqrt[z]]
```

## Hypergeometric1F1


`Hypergeometric1F1[a, b, z]` is Kummer's confluent hypergeometric function
₁F₁, equal to `HypergeometricPFQ[{a}, {b}, z]`.

```mathematica
In[1]:= Hypergeometric1F1[a, b, 0]
Out[1]= 1
```

## Hypergeometric2F1


`Hypergeometric2F1[a, b, c, z]` is the Gauss hypergeometric function ₂F₁,
equal to `HypergeometricPFQ[{a, b}, {c}, z]`.

```mathematica
In[1]:= Hypergeometric2F1[1, 1, 2, z]
Out[1]= -Log[1 - z]/z
```

## Zeta

- `Zeta[s]` — the Riemann zeta function ζ(s) = Σ_{k≥1} k⁻ˢ (Re s > 1; analytic
  continuation elsewhere).
- `Zeta[s, a]` — the Hurwitz (generalized) zeta function ζ(s, a) = Σ_{k≥0} (k+a)⁻ˢ.

Attributes: `Listable`, `NumericFunction`, `Protected`.

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
  `Zeta[s, 1] = Zeta[s]`.
- **Numerics.** Real one-argument zeta uses MPFR's `mpfr_zeta` (machine or
  arbitrary precision). All complex inputs, and all two-argument (Hurwitz)
  inputs, use one Euler–Maclaurin complex-MPFR kernel (Riemann is the `a = 1`
  case). Precision tracks the input, e.g.
  `N[Zeta[3], 50] = 1.2020569031595942853997381615114499907649862923405…`,
  `N[Zeta[5/4], 50] = 4.5951118258429433806853780396946256522810297806048…`,
  `N[Zeta[1/2 + I/2]] = -0.459303 - 0.961254 I`,
  `Zeta[-1.5 + I, 2.5 - I] = 0.0184868 + 1.67553 I`.
- **Derivatives.** `D[Zeta[s, a], a] = -s Zeta[1+s, a]` (the `s`-derivative is the
  generic `Derivative[1,0][Zeta][s, a]`); `D[Zeta[s], s] = Derivative[1][Zeta][s]`.
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
```

## StieltjesGamma

`StieltjesGamma[n]` is the n-th Stieltjes constant γₙ, the coefficients of the
Laurent expansion of `Zeta` about s = 1:
ζ(s) = 1/(s-1) + Σ_{n≥0} ((-1)ⁿ/n!) γₙ (s-1)ⁿ.

Attributes: `Listable`, `Protected`. It is inert: `StieltjesGamma[0]` reduces to
`EulerGamma`, and higher indices stay symbolic (they have no elementary closed
form). It appears in the `Series` expansions of `Zeta`.

```mathematica
In[1]:= StieltjesGamma[0]
Out[1]= EulerGamma
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
- `PolyLog[n, p, z]` — the Nielsen generalized polylogarithm S_{n,p}(z). Accepted
  for surface compatibility but left unevaluated (no closed-form engine).

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
  - The dilogarithm and trilogarithm at ½:
    `PolyLog[2, 1/2] = Pi²/12 - Log[2]²/2`,
    `PolyLog[3, 1/2] = Log[2]³/6 - Pi² Log[2]/12 + 7 Zeta[3]/8`.
- **Numerics.** With at least one inexact operand, machine- and
  arbitrary-precision (MPFR) real and complex inputs evaluate numerically:
  a direct power series for `|z| ≤ 1/2` (a fast pure-real path for real order and
  real −1 < z < 1), and the Jonquière / zeta expansion
  Liₛ(z) = Γ(1−s)(−ln z)^{s−1} + Σ_{k≥0} ζ(s−k) (ln z)ᵏ/k! (with the
  integer-order variant carrying the Hₙ₋₁ − ln(−ln z) term) for
  1/2 < |z| with |ln z| < 2π. The required ζ values for a complex order use the
  ζ functional equation in the left half-plane to stay well-conditioned.
  Examples: `PolyLog[2, 0.9] = 1.29971`, `PolyLog[0, 5.0] = -1.25`,
  `N[PolyLog[1, 1/3], 50] = 0.40546510810816438197801311546434913657199042346249`,
  `PolyLog[.2 + I, .5 - I] = -0.0898526 - 0.595865 I`. The branch cut runs along
  the real axis from 1 to ∞; values on it are taken continuous from below,
  matching Mathematica (e.g. Im Li₂(x) = −π Log[x] for real x > 1, so
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
```
