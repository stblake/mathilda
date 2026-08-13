# Erfi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Erfi[z]`**

gives the imaginary error function erfi(z) = -I Erf\[I z\] = (2/Sqrt\[Pi\]) Integral\_0^z e^(t^2) dt.

**`Erfi[0] = 0, Erfi[Infinity] = Infinity, Erfi[I Infinity] = I. An entire`**

<details>
<summary>Notes</summary>

function, odd in z. Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[Erfi\[z\], z\] = (2/Sqrt\[Pi\]) E^(z^2). Listable.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

```mathematica
In[1]:= Erfi[0]
Out[1]= 0

In[2]:= Erfi[-x]
Out[2]= -Erfi[x]

In[3]:= N[Erfi[1], 30]
Out[3]= 1.650425758797542876025337729561

In[4]:= Series[Erfi[x], {x, 0, 7}]
Out[4]= 2/Sqrt[Pi] x + 2/3/Sqrt[Pi] x^3 + 1/5/Sqrt[Pi] x^5 + 1/21/Sqrt[Pi] x^7 + O[x]^8

In[5]:= D[Erfi[x], x]
Out[5]= (2 E^x^2)/Sqrt[Pi]
```

## Algorithm

Mathilda -- the imaginary error function.

```text
  Erfi[z]   imaginary error function   erfi(z) = erf(i z)/i
                                               = (2/sqrt(pi)) Int_0^z e^t^2 dt
```

erfi is an entire function (no branch cuts) and odd in z. Evaluation is layered so each kind of argument takes the cheapest route:

```text
  exact special values   ->  0, +-Infinity, +-I (the imaginary-axis limits
                             are FINITE: erfi(i y) -> i as y -> +Infinity)
  symbolic odd argument   ->  Erfi[-x] = -Erfi[x]
  machine / arbitrary real -> the all-positive Maclaurin series
                             erfi(x) = (2/sqrt(pi)) Sum x^(2n+1)/(n!(2n+1)),
                             evaluated in MPFR. For real x every term shares
                             x's sign, so the partial sums climb monotonically
                             to the result -- no cancellation, no e^x^2
                             prefactor, pure-real arithmetic.
  complex (any precision)  -> erfi(z) = -i erf(i z), reusing the
                             cancellation-aware erf series (DLMF 7.6.2) in
                             MPFR with |z|^2/ln2 guard bits, so even
                             machine-precision complex results carry full
                             accuracy. Double-precision series is the
                             USE_MPFR=0 fallback.
  everything else         ->  stays symbolic (return NULL)
```

There is no libm erfi and no mpfr_erfi, so both numeric kernels are hand-rolled here.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

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

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Erf](../../special-functions/Erf/), [D](../../calculus/D/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_cherry_ei.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_ei.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_erfi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_erfi.c)
- Tests: [`tests/test_intrischnorman.c`](https://github.com/stblake/mathilda/blob/main/tests/test_intrischnorman.c)

## Notes & additional examples

### Notes

`Erfi[z] = -I Erf[I z] = (2/Sqrt[Pi]) Integral_0^z e^(t^2) dt` is the imaginary
error function, an entire odd function with `Erfi[0] = 0`,
`Erfi[Infinity] = Infinity`, and `Erfi[I Infinity] = I`. Compared with `Erf`, the
sign of every Maclaurin coefficient is positive, reflecting the `+t^2` in the
integrand. Real and complex arguments evaluate numerically at machine or
arbitrary (MPFR) precision, the derivative is `(2/Sqrt[Pi]) E^(z^2)`, and `Erfi`
is `Listable`.
