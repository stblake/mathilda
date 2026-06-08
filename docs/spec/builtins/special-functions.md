# Special Functions

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
  arguments; the `a`-derivative is the generic `Derivative[1,0][Gamma][a,z]`,
  and `D[Gamma[z], z]` stays `Derivative[1][Gamma][z]` since `PolyGamma`
  is not yet implemented).
- All other arguments (e.g. `Gamma[1/3]`, `Gamma[x]`, `Gamma[a, z]`,
  exact non-integer `Gamma[3/2, z]`, exact complex `Gamma[3/2, I]`) stay
  unevaluated.

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

## EulerGamma

`EulerGamma` is the Euler–Mascheroni constant γ, the limiting value of
`HarmonicNumber[n] - Log[n]` as `n → ∞`, with numerical value
≈ 0.5772156649015328606.

**Features**:
- Attributes `Constant`, `Protected`. `Attributes[EulerGamma] = {Constant,
  Protected}`; the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[EulerGamma]` is
  `True` and `D[EulerGamma, x] = 0`.
- `N[EulerGamma]` gives the machine value `0.577216`; `N[EulerGamma, prec]`
  gives any precision (MPFR `mpfr_const_euler`), e.g.
  `N[EulerGamma, 50] = 0.57721566490153286060651209008240243104215933593992`.
- Participates in exact numeric work, e.g.
  `Round[1/EulerGamma^100] = 734833795660954410469466`, and digit/continued-
  fraction extraction, e.g. `RealDigits[EulerGamma, 10, 50, -10^4]` returns
  decimal digits 10000–10049.

The constant value lives in the central numeric constant table (`src/numeric.c`);
the symbol's identity (attributes) is registered in `src/eulergamma.c`.

## Pi

`Pi` is π, with numerical value ≈ 3.14159.

**Features**:
- Attributes `Constant`, `Protected`. `Attributes[Pi] = {Constant, Protected}`;
  the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[Pi]` is `True` and
  `D[Pi, x] = 0`.
- `N[Pi]` gives the machine value `3.14159`; `N[Pi, prec]` gives any precision
  (MPFR `mpfr_const_pi`), e.g.
  `N[Pi, 50] = 3.1415926535897932384626433832795028841971693993751`.
- Participates in exact numeric work, e.g.
  `Round[Pi^100] = 51878483143196131920862615246303013562686760680406`.

## E

`E` is the exponential constant *e* (base of natural logarithms), with
numerical value ≈ 2.71828.

**Features**:
- Attributes `Constant`, `Protected`. `Attributes[E] = {Constant, Protected}`;
  the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[E]` is `True` and
  `D[E, x] = 0`.
- `N[E]` gives the machine value `2.71828`; `N[E, prec]` gives any precision
  (MPFR `mpfr_exp` of 1), e.g.
  `N[E, 50] = 2.71828182845904523536028747135266249775724709369996`.
- Participates in exact numeric work, e.g.
  `Round[E^100] = 26881171418161354484126255515800135873611119`.

## Catalan

`Catalan` is Catalan's constant *G*, the sum `Sum_{k>=0} (-1)^k (2 k + 1)^-2`,
with numerical value ≈ 0.915966.

**Features**:
- Attributes `Constant`, `Protected`. `Attributes[Catalan] = {Constant,
  Protected}`; the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[Catalan]` is `True` and
  `D[Catalan, x] = 0`.
- `N[Catalan]` gives the machine value `0.915966`; `N[Catalan, prec]` gives any
  precision (MPFR `mpfr_const_catalan`), e.g.
  `N[Catalan, 50] = 0.915965594177219015054603514932384110774149374281673`.

## GoldenRatio

`GoldenRatio` is the golden ratio φ = `(1 + Sqrt[5])/2`, the positive root of
`x^2 == x + 1`, with numerical value ≈ 1.61803.

**Features**:
- Attributes `Constant`, `Protected`. `Attributes[GoldenRatio] = {Constant,
  Protected}`; the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[GoldenRatio]` is `True`
  and `D[GoldenRatio, x] = 0`.
- `N[GoldenRatio]` gives the machine value `1.61803`; `N[GoldenRatio, prec]`
  gives any precision, e.g.
  `N[GoldenRatio, 50] = 1.61803398874989484820458683436563811772030917980576`.

## Degree

`Degree` gives the number of radians in one degree, with numerical value
`Pi/180` ≈ 0.0174533. Multiply by `Degree` to convert from degrees to radians,
so `30 Degree` represents 30°.

**Features**:
- Attributes `Constant`, `Protected`. `Attributes[Degree] = {Constant,
  Protected}`; the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[Degree]` is `True` and
  `D[Degree, x] = 0`.
- Used in arguments of trigonometric functions to express angles in degrees,
  e.g. `30 Degree` is `π/6`; the trig value evaluates numerically under `N`,
  e.g. `N[Cos[30 Degree]] = 0.866025` (the exact symbolic form
  `Cos[30 Degree]` is left unevaluated).
- `N[Degree]` gives the machine value `0.0174533`; `N[Degree, prec]` gives any
  precision, e.g.
  `N[Degree, 50] = 0.0174532925199432957692369076848861271344287188854173`.

The constant values for `Pi`, `E`, `Catalan`, `GoldenRatio`, and `Degree` all
live in the central numeric constant table (`src/numeric.c`); their identities
(attributes `Constant`, `Protected`) are stamped in `numeric_init`.

## GoldenAngle

`GoldenAngle` is the golden angle `(3 - Sqrt[5]) Pi` = `2 Pi / GoldenRatio^2`,
with numerical value ≈ 2.39996 radians ≈ 137.5°.

**Features**:
- Attributes `Constant`, `Protected`. `Attributes[GoldenAngle] = {Constant,
  Protected}`; the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[GoldenAngle]` is `True`
  and `D[GoldenAngle, x] = 0`.
- `N[GoldenAngle]` gives the machine value `2.39996`; `N[GoldenAngle, prec]`
  gives any precision (computed in MPFR from the closed form `(3 - Sqrt[5])
  Pi`), e.g.
  `N[GoldenAngle, 50] = 2.3999632297286533222315555066336138531249990110581`.

## Glaisher

`Glaisher` is the Glaisher–Kinkelin constant *A*, with numerical value
≈ 1.2824271291. It satisfies `Log[A] == 1/12 − Zeta'(−1)`, where `Zeta` is the
Riemann zeta function.

**Features**:
- Attributes `Constant`, `Protected`. `Attributes[Glaisher] = {Constant,
  Protected}`; the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[Glaisher]` is `True`
  and `D[Glaisher, x] = 0`.
- `N[Glaisher]` gives the machine value `1.28243`; `N[Glaisher, prec]` gives any
  precision, e.g.
  `N[Glaisher, 50] = 1.282427129100622636875342568869791727767688927325`.

  Arbitrary precision is computed from `ln A = (γ + ln(2π))/12 − ζ'(2)/(2π²)`,
  with `ζ'(2)` evaluated by Euler–Maclaurin summation of `−Σ ln(n)/n²` (the
  even Bernoulli numbers obtained from `Zeta[2k]`). Verified to 250 digits.

## Khinchin

`Khinchin` is Khinchin's constant *K* (also "Khintchine's constant"), with
numerical value ≈ 2.6854520011. For almost every real number, the partial
quotients `a_k` of its continued-fraction expansion have the same limiting
geometric mean `K = lim (a_1 a_2 ⋯ a_k)^(1/k)`, with closed form
`K = Product_{s>=1} (1 + 1/(s(s+2)))^(log_2 s)`.

**Features**:
- Attributes `Constant`, `Protected`. `Attributes[Khinchin] = {Constant,
  Protected}`; the symbol cannot be reassigned.
- Propagated as an exact, unevaluated symbol; `NumericQ[Khinchin]` is `True`
  and `D[Khinchin, x] = 0`.
- `N[Khinchin]` gives the machine value `2.68545`; `N[Khinchin, prec]` gives any
  precision, e.g.
  `N[Khinchin, 50] = 2.6854520010653064453097148354817956938203822939945`.

  Arbitrary precision uses the geometrically convergent zeta series
  `ln K · ln 2 = Σ_{n>=1} (ζ(2n) − 1)/n · Σ_{k=1}^{2n−1} (−1)^(k+1)/k`
  (the Bailey–Borwein–Crandall form). Verified to 250 digits.

The constant values for `GoldenAngle`, `Glaisher`, and `Khinchin` live in the
numeric constant table (`src/numeric.c`); their MPFR fillers compute
`GoldenAngle` from its closed form, and `Glaisher`/`Khinchin` from the series
above. Their `Constant`/`Protected` attributes are stamped in `numeric_init`.

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
