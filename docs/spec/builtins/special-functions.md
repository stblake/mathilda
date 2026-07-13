# Special Functions

Higher transcendental functions: the gamma function `Gamma`, the error function `Erf`, its complement `Erfc` and the imaginary error function `Erfi`, the digamma/polygamma family `PolyGamma`, the log-gamma function `LogGamma`, the harmonic numbers `HarmonicNumber`, the Pochhammer symbol (rising factorial) `Pochhammer`, the Riemann/Hurwitz zeta function `Zeta` (with the inert Stieltjes constants `StieltjesGamma`), the Hurwitz zeta function `HurwitzZeta`, the Bernoulli numbers and polynomials `BernoulliB`, the Euler numbers and polynomials `EulerE`, the polylogarithm `PolyLog`, the Lerch transcendent `LerchPhi`, the hypergeometric family `Hypergeometric0F1`, `Hypergeometric1F1`, `Hypergeometric2F1`, and the generalized `HypergeometricPFQ`, the Airy functions `AiryAi` and `AiryBi`, the Lambert W function `ProductLog`, and the Legendre polynomials and associated Legendre functions `LegendreP`.

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

## Beta

- `Beta[a, b]` — the Euler beta function B(a, b) = Γ(a)Γ(b)/Γ(a+b)
  = ∫₀¹ tᵃ⁻¹(1−t)ᵇ⁻¹ dt.
- `Beta[z, a, b]` — the incomplete beta function B_z(a, b)
  = ∫₀ᶻ tᵃ⁻¹(1−t)ᵇ⁻¹ dt (branch cut along z < 0).
- `Beta[z0, z1, a, b]` — the generalized incomplete beta
  B(z0, z1; a, b) = Beta[z1, a, b] − Beta[z0, a, b].

Attributes: `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.
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
- Numeric `Beta[a, b]`: machine and arbitrary-precision (MPFR) real and complex
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
  - Real (machine *and* arbitrary precision) → the all-positive Maclaurin series
    erfi(x) = (2/√π) Σ x^(2n+1)/(n!(2n+1)). For real x every term shares x's
    sign, so the partial sums climb monotonically to the result — no
    cancellation, no e^(x²) prefactor — evaluated in MPFR with output precision
    tracking the input, e.g. `Erfi[2.5] = 130.396`, `Erfi[0.5] = 0.614952`,
    `N[Erfi[1/2], 50] = 0.61495209469651098083968118562364139305134561789540`.
  - **Complex** (machine *and* arbitrary precision) → erfi(z) = -i erf(i z),
    reusing the cancellation-aware erf series (DLMF 7.6.2) in MPFR with
    `|z|²/ln2` guard bits, so even machine-precision complex results carry full
    accuracy, e.g. `Erfi[1.5 - I] = -0.70136 - 1.84683 I`,
    `N[Erfi[1/2 + I], 30] = 0.187973467223383313628263810077 + 0.950709728318957173804611826379 I`.
    A `double complex` series is the fallback for `USE_MPFR=0` builds.
- Derivative: `D[Erfi[z], z] = (2/Sqrt[Pi]) E^(z^2)` (positive exponent, vs
  `Erf`'s E^(−z²); chain rule applies), so the origin Taylor series follows from
  the generic `D`-based fallback, e.g. `Series[Erfi[x], {x, 0, 7}]` begins
  `2/Sqrt[Pi] x + 2/(3 Sqrt[Pi]) x^3 + …`.
- All other arguments (symbolic `Erfi[x]`, exact `Erfi[2]`) stay unevaluated.

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
- Derivative: `D[ExpIntegralEi[z], z] = E^z/z` (chain rule applies, e.g.
  `D[ExpIntegralEi[x^2], x] = (2 E^x^2)/x`); the origin Taylor series at a regular
  point follows from the generic `D`-based fallback.
- Wrong arity emits `ExpIntegralEi::argx` and stays unevaluated.

## LogIntegral

- `LogIntegral[z]` — the logarithmic integral li(z), the principal value of
  ∫₀^z dt/ln t.

**Attributes**: `Listable`, `NumericFunction`, `Protected`.

`LogIntegral` is computed from the identity **li(z) = Ei(Log z)**, reusing the
`ExpIntegralEi` numeric kernels; the principal `Log` supplies the ±i π jump, so
the branch cut runs along (−∞, +1).

**Features**:
- Exact special values: `LogIntegral[0] = 0`, `LogIntegral[1] = -Infinity`,
  `LogIntegral[Infinity] = Infinity`; `ComplexInfinity` and `Indeterminate` map
  to `Indeterminate`.
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

## Sinc

- `Sinc[z]` — the cardinal sine Sin[z]/z, with the removable singularity filled
  in as `Sinc[0] = 1`.

`Sinc` is entire and even. Real arguments use `mpfr_sin(z)/z`; complex arguments
use `sin(z)/z` via the shared `ncpx` toolkit. Provided as a first-class head so
`SinIntegral` differentiates to it.

- Exact special values: `Sinc[0] = 1`, `Sinc[±Infinity] = 0`;
  `ComplexInfinity` maps to `Indeterminate`.
- Numeric evaluation at machine or arbitrary (MPFR) precision, tracking input
  precision: `Sinc[2.] = 0.454649`,
  `N[Sinc[2], 45] = 0.454648713412840847698009932955872421351127485`. Complex:
  `Sinc[1. + I] = 0.966711 - 0.331747 I`.
- Symbolic arguments stay unevaluated (`Sinc[x]`, `Sinc[-x]`; WL does not
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
- **Exact at rational arguments (Gauss digamma theorem).** For order `n = 0` and
  a positive rational `p/q`, `PolyGamma[0, p/q]` closes into an elementary form
  `−EulerGamma − ln(2q) − (π/2) cot(πp/q) + 2 Σ cos(2πnp/q) ln(sin(nπ/q))`
  (integer part folded back via `ψ(x+1) = ψ(x) + 1/x`), e.g.
  `PolyGamma[0, 1/2] = -EulerGamma - Log[4]`,
  `PolyGamma[0, 3/4] - PolyGamma[0, 1/4] = Pi`,
  `PolyGamma[0, 2/3] - PolyGamma[0, 1/3] = Pi/Sqrt[3]`.
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
  `2F1(1,1;2) -> -Log[1-z]/z`. The central-binomial / arcsin family also closes:
  `2F1(1,1;1/2;z) -> 1/(1-z) + Sqrt[z] ArcSin[Sqrt[z]]/(1-z)^(3/2)`,
  `2F1(1,1;3/2;z) -> ArcSin[Sqrt[z]]/(Sqrt[z] Sqrt[1-z])`,
  `2F1(2,1;3/2;z) -> 1/(2(1-z)) + ArcSin[Sqrt[z]]/(2 Sqrt[z] (1-z)^(3/2))`.
  A table-backed very-well-poised reduction closes the classic
  central-binomial-cubed Ramanujan `1/Pi` series:
  `4F3({1/2,1/2,1/2,5/4}, {1/4,1,1}, -1) -> 2/Pi` (general `1/Pi` summation is an
  open problem, so this is a recognizer for the known class, not universal).
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
  `Zeta[s, 1] = Zeta[s]` and `Zeta[s, 1/2] = (2^s - 1) Zeta[s]`
  (e.g. `Zeta[2, 1/2] = Pi^2/2`, `Zeta[3, 1/2] = 7 Zeta[3]`).
- **Numerics.** Real one-argument zeta uses MPFR's `mpfr_zeta` (machine or
  arbitrary precision). All complex inputs, and all two-argument (Hurwitz)
  inputs, use one Euler–Maclaurin complex-MPFR kernel (Riemann is the `a = 1`
  case). Precision tracks the input, e.g.
  `N[Zeta[3], 50] = 1.2020569031595942853997381615114499907649862923405…`,
  `N[Zeta[5/4], 50] = 4.5951118258429433806853780396946256522810297806048…`,
  `N[Zeta[1/2 + I/2]] = -0.459303 - 0.961254 I`,
  `Zeta[-1.5 + I, 2.5 - I] = 0.0184868 + 1.67553 I`. The two-argument kernel uses
  the symmetric power `((a+k)^2)^(-s/2)` (= `|a+k|^-s` on the real axis), so for
  `Re(a) < 0` it differs from the principal-branch `HurwitzZeta`:
  `Zeta[3, -1.5] = 16.7107` while `HurwitzZeta[3, -1.5] = 0.118102`.
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
- **Numerics.** Real and complex, machine- and arbitrary-precision arguments
  evaluate through one Euler–Maclaurin complex-MPFR kernel; a numeric pole gives
  `ComplexInfinity`. Precision tracks the input, e.g.
  `HurwitzZeta[3, .2] = 125.739`, `HurwitzZeta[7., 5] = 1.84949×10⁻⁵`,
  `HurwitzZeta[2.3, 8 + I] = 0.0544701 - 0.00944852 I`,
  `N[HurwitzZeta[1/3, 8/7], 50] = -1.1389367444490991746548674334535727810961919460755…`.
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

## LerchPhi

- `LerchPhi[z, s, a]` — the Lerch transcendent Φ(z, s, a) = Σ_{k≥0} zᵏ/(k+a)ˢ
  (|z| < 1; analytic continuation elsewhere, branch cut z ∈ [1, ∞)). For Re a < 0
  the principal value uses the symmetric power ((k+a)²)^(−s/2) and excludes any
  term with k + a = 0.

`LerchPhi` is the common generalization of `Zeta`, `HurwitzZeta` and `PolyLog`:
`LerchPhi[1, s, a] = Zeta[s, a]` and `z LerchPhi[z, s, 1] = PolyLog[s, z]`.
Implemented in `src/special_functions/lerchphi.c`, registered via
`lerchphi_init()` in `core_init()`.

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
  - A non-positive integer `a = -m` drops the singular `k = m` term (Wolfram's
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
  For `|z| > 1` (off the branch cut `[1, ∞)`, with `s` not a positive integer,
  `a` not a non-positive integer, and `|Log z| < 2π`) the value comes from the
  Lerch/Erdélyi continuation
  `Φ(z,s,a) = z^(-a)[Γ(1-s)(-Log z)^(s-1) + Σ_{n≥0} HurwitzZeta(s-n,a)(Log z)^n/n!]`,
  summed with optimal truncation: `LerchPhi[5. + I, I, I + 2] =
  -0.581502 + 0.384767 I` (cross-checked against `z LerchPhi[z,s,1] =
  PolyLog[s,z]`). Outside that domain — integer `s`, or `z` on the cut
  `[1, ∞)` (e.g. `LerchPhi[2, 3, -1.5]`, which needs the logarithmic confluent
  form) — `|z| > 1` stays symbolic.
- **Series at `z = 0`.** `Series[LerchPhi[z, s, a], {z, 0, n}]` returns the
  defining power series, coefficient of `z^k` being `((k+a)^2)^(-s/2)`:
  `Series[LerchPhi[z, s, a], {z, 0, 5}] // Normal =
  (a^2)^(-s/2) + ((1+a)^2)^(-s/2) z + …`.
- **Derivatives.** `D[LerchPhi[z, s, a], z] = (LerchPhi[z, -1+s, a]
  - a LerchPhi[z, s, a])/z` and `D[LerchPhi[z, s, a], a] = -s LerchPhi[z, 1+s, a]`;
  the `s`-derivative is the generic `Derivative[0, 1, 0][LerchPhi][z, s, a]`.
- **Diagnostics.** Fewer than three arguments emit `LerchPhi::argrx`; a
  non-option beyond position 3 emits `LerchPhi::nonopt`; both leave the call
  unevaluated.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= LerchPhi[z, s, 1]
Out[1]= PolyLog[s, z]/z

In[2]:= LerchPhi[0.5, 3, 2.5]
Out[2]= 0.0794983
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
- **Numerics.** With an inexact argument, machine- and arbitrary-precision
  (MPFR) real and complex inputs evaluate numerically. The unified core routes
  on `r = |z|` and the requested precision `P`:
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

Attributes (both heads): `Listable`, `NumericFunction`, `Protected`,
`ReadProtected`.

```mathematica
In[1]:= AiryAi[0]
Out[1]= 1/(3^(2/3) Gamma[2/3])

In[2]:= AiryAi[1.8]
Out[2]= 0.0470362
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

Attributes (both heads): `Listable`, `NumericFunction`, `Protected`,
`ReadProtected`.

```mathematica
In[1]:= AiryBi[0]
Out[1]= 1/(3^(1/6) Gamma[2/3])

In[2]:= AiryBi[1.8]
Out[2]= 2.59587
```

## BesselJ

- `BesselJ[n, z]` — the Bessel function of the first kind Jₙ(z), the solution of
  z² y″ + z y′ + (z² − n²) y = 0 that is regular at the origin. For non-integer
  order it has a branch cut along the negative real z axis (from the (z/2)ⁿ
  factor).

Implemented in `src/special_functions/besselj.c`. The complex-MPFR arithmetic
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
- **Numerics.** When some argument is inexact, machine- and arbitrary-precision
  (MPFR) real and complex order and argument evaluate numerically. The core
  routes on `r = |z|` and the requested precision `P`:
  - integer order with a real argument takes MPFR's correctly-rounded
    `mpfr_jn` fast path (J₋ₙ = (−1)ⁿ Jₙ);
  - a **power series** for small/moderate `r` (DLMF 10.2.2,
    `t_k = t_{k-1}·(−(z/2)²)/(k(ν+k))`, one `Gamma(ν+1)` via `mpfr_gamma` for
    real order or the `Gamma` builtin for complex order, with ~`2r/ln2` guard
    bits to absorb the alternating-series cancellation);
  - an **asymptotic series** for large `r` away from the negative-real axis
    (DLMF 10.17.3, summed to the optimal smallest-term truncation).

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
  `src/internal/bessel.m`.
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

Attributes: `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

Not yet implemented: `Limit[BesselJ[n, x], x -> Infinity]` is left unevaluated
(the limit engine does not yet consume the oscillatory asymptotic leading term);
the symbolic-index general term `SeriesCoefficient[BesselJ[0, x], {x, 0, n}]`
(the Piecewise closed form) — concrete-index `SeriesCoefficient` works.

```mathematica
In[1]:= BesselJ[0, 5.2]
Out[1]= -0.11029

In[2]:= D[BesselJ[n, x], x]
Out[2]= 1/2 (BesselJ[-1 + n, x] - BesselJ[1 + n, x])
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
  even-order reflection `BesselK[-n, x] = BesselK[n, x]` live in
  `src/internal/bessel.m`.
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

Attributes: `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

On the branch cut (negative real axis) the principal-branch value is returned
consistently at all precisions; this convention can differ from other systems'
large-order machine-precision heuristics.

```mathematica
In[1]:= BesselK[0, 4.0]
Out[1]= 0.0111597

In[2]:= D[BesselK[n, x], x]
Out[2]= -1/2 (BesselK[-1 + n, x] + BesselK[1 + n, x])
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
  `src/internal/bessel.m`.
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
  Mathematica's `Normal` form, valid for symbolic order too.
- **Indefinite integrals** (recurrence forms, `src/internal/CRCMathTablesIntegrals.m`):
  `∫ xᵖ⁺¹ Iₚ(x) dx = xᵖ⁺¹ Iₚ₊₁(x)`, `∫ x⁻ᵖ⁺¹ Iₚ(x) dx = x⁻ᵖ⁺¹ Iₚ₋₁(x)`,
  `∫ I₁(x) dx = I₀(x)`.
- **Listable.** `BesselI[{0, 1, 2}, 1.] = {1.26607, 0.565159, 0.135748}`.

Attributes: `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

Not yet implemented: `Limit[BesselI[n, x], x -> Infinity]` (= ∞) and at
`I Infinity` (= 0) are left unevaluated — the Limit engine does not distribute
over the `Plus` of the two-exponential asymptotic at infinity (the `Series`
itself is exact).

```mathematica
In[1]:= BesselI[0, 2.0]
Out[1]= 2.27959

In[2]:= D[BesselI[n, x], x]
Out[2]= 1/2 (BesselI[-1 + n, x] + BesselI[1 + n, x])
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
  `src/internal/bessel.m`.
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

Attributes: `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

```mathematica
In[1]:= BesselY[0, 2.5]
Out[1]= 0.49807

In[2]:= D[BesselY[n, x], x]
Out[2]= 1/2 (BesselY[-1 + n, x] - BesselY[1 + n, x])
```

## Hyperfactorial

`Hyperfactorial[n]` gives the hyperfactorial `prod_{k=1}^{n} k^k`
(`H(0) = H(1) = 1`).  Exact via GMP for a non-negative integer `n`;
non-positive-integer, non-integer or symbolic orders are left
unevaluated.  Used by `Product` to recognise `Product[k^k, {k, 1, n}]`.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= Hyperfactorial[4]
Out[1]= 27648
```

## BarnesG

`BarnesG[z]`, the Barnes G-function, with `G(1) = G(2) = 1` and
`G(z+1) = Gamma[z] G(z)`.  For a positive integer `n`,
`G(n+1) = prod_{k=1}^{n-1} k!` (the superfactorial, exact via GMP), and
`G(m) = 0` for non-positive integer `m`.  Non-integer orders are left
unevaluated.  Used by `Product` to recognise
`Product[Gamma[i], {i, 1, n-1}]` → `BarnesG[n]`.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= BarnesG[5]
Out[1]= 12
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

Attributes: `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

```mathematica
In[1]:= ProductLog[1.0]
Out[1]= 0.567143

In[2]:= ProductLog[-1/E]
Out[2]= -1
```

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

Non-integer associated/type forms, negative `m`, symbolic `Series` /
`SeriesCoefficient`, `D[]` rules, and analytic continuation of the numeric
series for `|(1-x)/2| ≥ 1` are left symbolic.

Attributes: `Listable`, `NumericFunction`, `Protected`.

```mathematica
In[1]:= LegendreP[3, x]
Out[1]= -3/2 x + 5/2 x^3

In[2]:= LegendreP[10, 2, x]
Out[2]= (1 - x^2) (3465/128 - 45045/32 x^2 + 675675/64 x^4 - 765765/32 x^6 + 2078505/128 x^8)
```
