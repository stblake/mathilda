# Mathematical Constants

Exact symbolic constants with arbitrary-precision numeric values via `N`. Each is a protected symbol that stays unevaluated in exact computation and is filled to the requested precision on demand: the circle constant `Pi`, Euler's number `E`, the degree `Degree` (= Pi/180), the Euler–Mascheroni constant `EulerGamma`, `Catalan`'s constant, the `GoldenRatio` and `GoldenAngle`, and the `Glaisher` and `Khinchin` constants.

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
the symbol's identity (attributes) is registered in `src/special_functions/eulergamma.c`.

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

## I

`I` is the imaginary unit, Sqrt[-1].

**Features**:
- Attribute `Protected`. `Attributes[I] = {Protected}`; the symbol cannot be
  reassigned.
- Carries the OwnValue `Complex[0, 1]`, so `I` evaluates to the imaginary unit;
  `I^2 = -1` and complex numbers are written `a + b I`.
- `N[I] = 0. + 1. I`.

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
