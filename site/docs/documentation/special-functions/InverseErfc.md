# InverseErfc

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`InverseErfc[s]`**

gives the inverse complementary error function: the z solving s = Erfc\[z\].

**`InverseErfc[0] = Infinity, InverseErfc[1] = 0, InverseErfc[2] = -Infinity.`**

<details>
<summary>Notes</summary>

Numerical values are given only for real s in \[0, 2\], at machine or arbitrary (MPFR) precision; D\[InverseErfc\[z\], z\] = -(Sqrt\[Pi\]/2) E^(InverseErfc\[z\]^2). Listable.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

```mathematica
In[1]:= InverseErfc[1]
Out[1]= 0

In[2]:= InverseErfc[2]
Out[2]= -Infinity

In[3]:= N[InverseErfc[1/1000000], 40]
Out[3]= 3.4589107372795000221509276359575695199155

In[4]:= D[InverseErfc[z], z]
Out[4]= -1/2 Sqrt[Pi] E^InverseErfc[z]^2

In[5]:= InverseErf[1 - 3/10] == InverseErfc[3/10]
Out[5]= True
```

## Algorithm

Mathilda -- the inverse complementary error function.

```text
  InverseErfc[s]   inverse of erfc: the z solving s = erfc(z)
```

Since erfc maps the real line onto (0, 2) (decreasing from 2 at -Infinity to 0 at +Infinity), explicit numerical values are produced only for *real* s in [0, 2]; out-of-domain (s < 0 or s > 2) and complex inputs stay symbolic. Evaluation mirrors InverseErf -- no complex machinery:

```text
  exact special values   ->  InverseErfc[0]=Infinity, InverseErfc[1]=0,
                             InverseErfc[2]=-Infinity
  machine real (0<s<2)   ->  Winitzki seed + Newton polish on libm erfc
  arbitrary real (0<s<2) ->  Newton with precision doubling on mpfr_erfc
  everything else        ->  stays symbolic (return NULL)
```

Mathematically InverseErfc[s] = InverseErf[1 - s], but we do NOT route through InverseErf: for small s (large z) forming 1 - s and inverting erf near 1 loses all significance to cancellation. Instead Newton iterates on erfc directly -- both libm and MPFR ship cancellation-free erfc -- with

```text
  f(z) = erfc(z) - s,  f'(z) = -(2/sqrt(pi)) e^{-z^2}, i.e.
  z <- z + (erfc(z) - s) (sqrt(pi)/2) e^{z^2}.
```

Newton converges quadratically; the MPFR path doubles the working precision each step so the final full-precision erfc dominates the cost.

The derivative D[InverseErfc[z], z] = -(sqrt(pi)/2) e^{InverseErfc[z]^2} lives in src/calculus/deriv.c; Series follows from it via the generic Taylor-via-D fallback. (erfc is not odd, so unlike InverseErf there is no auto-applied symmetry rewrite.)

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

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

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [InverseErf](../../special-functions/InverseErf/), [D](../../calculus/D/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_interval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interval.c)
- Tests: [`tests/test_inverfc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_inverfc.c)

## Notes & additional examples

### Notes

`InverseErfc[s]` returns the `z` solving `Erfc[z] == s`, with
`InverseErfc[0] = Infinity`, `InverseErfc[1] = 0`, `InverseErfc[2] = -Infinity`.
Numerical values are produced only for real `s` in `[0, 2]`, at machine or
arbitrary (MPFR) precision. It is the natural function for accurate evaluation of
extreme normal-distribution quantiles, where `1 - Erf` would lose all precision
to cancellation.
