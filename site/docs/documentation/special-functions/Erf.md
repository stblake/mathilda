# Erf

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Erf[z]`**

gives the error function erf(z) = (2/Sqrt\[Pi\]) Integral\_0^z e^(-t^2) dt.

**`Erf[z0, z1]`**

gives the generalized error function erf(z1) - erf(z0).

**`Erf[0] = 0, Erf[Infinity] = 1, Erf[-Infinity] = -1. An entire function,`**

<details>
<summary>Notes</summary>

odd in z. Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[Erf\[z\], z\] = (2/Sqrt\[Pi\]) E^(-z^2). Listable.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (7)

```mathematica
In[1]:= Erf[0]
Out[1]= 0

In[2]:= Erf[Infinity]
Out[2]= 1

In[3]:= Erf[-z]
Out[3]= -Erf[z]

In[4]:= N[Erf[1], 40]
Out[4]= 0.84270079294971486934122063508260925929605

In[5]:= N[Erf[1 + I], 20]
Out[5]= 1.31615128169794764488 + 0.190453469237834686284*I

In[6]:= Series[Erf[x], {x, 0, 7}]
Out[6]= 2/Sqrt[Pi] x + -2/3/Sqrt[Pi] x^3 + 1/5/Sqrt[Pi] x^5 + -1/21/Sqrt[Pi] x^7 + O[x]^8

In[7]:= D[Erf[x^2], x]
Out[7]= (4 x E^(-x^4))/Sqrt[Pi]
```

## Algorithm

Mathilda -- the error function.

```text
  Erf[z]       error function       erf(z)  = (2/sqrt(pi)) Int_0^z e^-t^2 dt
  Erf[z0, z1]  generalized error    erf(z1) - erf(z0)
```

erf is an entire function (no branch cuts) and odd in z. Evaluation is layered so each kind of argument takes the cheapest route:

```text
  exact special values   ->  0, +-1, DirectedInfinity[+-I], ...
  symbolic odd argument  ->  Erf[-x] = -Erf[x]
  machine real           ->  libm   erf
  arbitrary real         ->  MPFR   mpfr_erf
  complex (any precision) ->  the cancellation-aware Maclaurin series
                              (DLMF 7.6.2) evaluated in MPFR with guard
                              bits, so even machine-precision complex
                              results carry full accuracy. A double-complex
                              series is the fallback for USE_MPFR=0 builds.
  everything else        ->  stays symbolic (return NULL)

The series  erf(z) = (2/sqrt(pi)) e^-z^2 Sum_{n>=0} t_n,
  t_0 = z,  t_n = t_{n-1} (2 z^2)/(2n+1),
```

is convergent for every z. For real z all terms share a sign (no cancellation); for complex z the partial sums can reach magnitude ~e^|z|^2 before the e^-z^2 prefactor brings them back, so the MPFR path adds |z|^2/ln2 guard bits to absorb that cancellation exactly.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

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

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [D](../../calculus/D/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_erf.c`](https://github.com/stblake/mathilda/blob/main/tests/test_erf.c)
- Tests: [`tests/test_erfc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_erfc.c)

## Notes & additional examples

### Notes

`Erf[z]` is the error function `(2/Sqrt[Pi]) Integral_0^z e^(-t^2) dt`, an entire
odd function with the exact values `Erf[0] = 0` and `Erf[±Infinity] = ±1`. Real
and complex arguments evaluate numerically at machine or arbitrary (MPFR)
precision — the complex path uses a DLMF series so `N[Erf[1 + I], 20]` is correct
to the requested digits. The Maclaurin series and the chain-rule derivative
`D[Erf[z], z] = (2/Sqrt[Pi]) E^(-z^2)` are both built in, and `Erf` is `Listable`.
