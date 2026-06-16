# AiryAiPrime

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
AiryAiPrime[z]
    gives the derivative Ai'(z) of the Airy function AiryAi.
AiryAiPrime[0] = -1/(3^(1/3) Gamma[1/3]), AiryAiPrime[+Infinity] = 0. Real
and complex inputs evaluate numerically at machine or arbitrary (MPFR)
precision; D[AiryAiPrime[z], z] = z AiryAi[z]. Listable.
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
In[1]:= AiryAiPrime[0]
Out[1]= -1/(3^(1/3) Gamma[1/3])
```

```mathematica
In[1]:= N[AiryAiPrime[0], 40]
Out[1]= -0.25881940379280679840518356018920396347907
```

```mathematica
In[1]:= D[AiryAiPrime[z], z]
Out[1]= z AiryAi[z]
```

```mathematica
In[1]:= AiryAiPrime[1.0 + 1.0 I]
Out[1]= -0.130628 + 0.163068*I
```

### Notes

`AiryAiPrime[z]` is the derivative `Ai'(z)`. Its exact origin value is
`-1/(3^(1/3) Gamma[1/3])`, and differentiating once more recovers the Airy
equation in the form `D[AiryAiPrime[z], z] == z AiryAi[z]`. Complex arguments
evaluate to machine precision (and to arbitrary precision under `N[..., n]`);
`AiryAiPrime` is Listable.
