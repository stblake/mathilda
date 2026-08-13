# Erfc

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Erfc[z]`**

gives the complementary error function erfc(z) = 1 - erf(z).

**`Erfc[0] = 1, Erfc[Infinity] = 0, Erfc[-Infinity] = 2. An entire`**

**`D[Erfc[z], z] = -(2/Sqrt[Pi]) E^(-z^2). Listable.`**

<details>
<summary>Notes</summary>

function. Real inputs evaluate via libm/MPFR erfc (cancellation-free); complex inputs via 1 - erf(z) at machine or arbitrary (MPFR) precision.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

```mathematica
In[1]:= Erfc[0]
Out[1]= 1

In[2]:= N[Erfc[2], 40]
Out[2]= 0.0046777349810472658379307436327470713891081

In[3]:= N[Erfc[1 + I], 25]
Out[3]= -0.31615128169794764488027107 - 0.19045346923783468628410886*I

In[4]:= Series[Erfc[x], {x, 0, 5}]
Out[4]= 1 + -2/Sqrt[Pi] x + 2/3/Sqrt[Pi] x^3 + -1/5/Sqrt[Pi] x^5 + O[x]^6

In[5]:= D[Erfc[Sqrt[x]], x]
Out[5]= -E^(-x)/(Sqrt[Pi] Sqrt[x])
```

## Algorithm

Mathilda -- the complementary error function.

```text
  Erfc[z]   complementary error function   erfc(z) = 1 - erf(z)
```

erfc is an entire function (no branch cuts). It is the complement of erf; unlike erf it has no symmetry that simplifies Erfc[-x] (erfc(-x) = 2 - erfc(x), which Mathilda leaves unexpanded). Evaluation is layered so each kind of argument takes the cheapest, most accurate route:

```text
  exact special values    ->  1, 0, 2, DirectedInfinity[-+I], ...
  machine real            ->  libm   erfc
  arbitrary real          ->  MPFR   mpfr_erfc   (cancellation-free even for
                              large positive z, where 1 - erf(z) would lose
                              all significance)
  complex (any precision) ->  1 - erf(z), with erf(z) from the
                              cancellation-aware Maclaurin series (DLMF
                              7.6.2) evaluated in MPFR with guard bits; the
                              complement is formed at working precision
                              before rounding, so even machine-precision
                              complex results carry full accuracy. A
                              double-complex series is the USE_MPFR=0
                              fallback.
  everything else         ->  stays symbolic (return NULL)

The erf series  erf(z) = (2/sqrt(pi)) e^-z^2 Sum_{n>=0} t_n,
  t_0 = z,  t_n = t_{n-1} (2 z^2)/(2n+1),
```

is convergent for every z. For complex z the partial sums can reach magnitude ~e^|z|^2 before the e^-z^2 prefactor brings them back, so the MPFR path adds |z|^2/ln2 guard bits to absorb that cancellation exactly.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

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

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Erf](../../special-functions/Erf/), [D](../../calculus/D/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_erfc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_erfc.c)
- Tests: [`tests/test_fullsimplify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fullsimplify.c)
- Tests: [`tests/test_gruntz.c`](https://github.com/stblake/mathilda/blob/main/tests/test_gruntz.c)
- Tests: [`tests/test_gruntz_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_gruntz_stress.c)

## Notes & additional examples

### Notes

`Erfc[z] = 1 - Erf[z]` is the complementary error function, with `Erfc[0] = 1`,
`Erfc[Infinity] = 0`, and `Erfc[-Infinity] = 2`. Real inputs evaluate through the
cancellation-free `erfc` of libm / MPFR — important in the right tail, where the
`N[Erfc[2], 40]` example keeps full precision instead of losing it to a `1 - Erf`
subtraction. Complex inputs route through `1 - Erf[z]` at machine or arbitrary
precision. The derivative is `D[Erfc[z], z] = -(2/Sqrt[Pi]) E^(-z^2)`, and `Erfc`
is `Listable`.
