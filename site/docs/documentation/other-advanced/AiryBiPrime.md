# AiryBiPrime

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
AiryBiPrime[z]
    gives the derivative Bi'(z) of the Airy function AiryBi.
AiryBiPrime[0] = 3^(1/6)/Gamma[1/3], AiryBiPrime[+Infinity] = Infinity. Real
and complex inputs evaluate numerically at machine or arbitrary (MPFR)
precision; D[AiryBiPrime[z], z] = z AiryBi[z]. Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= AiryBiPrime[0]
Out[1]= 3^(1/6)/Gamma[1/3]
```

```mathematica
In[1]:= N[AiryBiPrime[0], 40]
Out[1]= 0.44828835735382635791482371039882839086616
```

```mathematica
In[1]:= D[AiryBiPrime[z], z]
Out[1]= z AiryBi[z]
```

```mathematica
In[1]:= AiryBiPrime[1.0 + 1.0 I]
Out[1]= 0.0756628 + 0.783701*I
```

### Notes

`AiryBiPrime[z]` is the derivative `Bi'(z)`. Its exact origin value is
`3^(1/6)/Gamma[1/3]`, and a further derivative satisfies the Airy equation in
the form `D[AiryBiPrime[z], z] == z AiryBi[z]`. Complex arguments evaluate to
machine precision and, under `N[..., n]`, to arbitrary MPFR precision;
`AiryBiPrime` is Listable.
