# AiryAi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
AiryAi[z]
    gives the Airy function Ai(z), the solution of y'' = z y that decays
    as z -> +Infinity.
AiryAi[0] = 1/(3^(2/3) Gamma[2/3]), AiryAi[+-Infinity] = 0. An entire
function of z. Real and complex inputs evaluate numerically at machine
or arbitrary (MPFR) precision; D[AiryAi[z], z] = AiryAiPrime[z]. Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= AiryAi[0]
Out[1]= 1/(3^(2/3) Gamma[2/3])

In[2]:= AiryAi[1.8]
Out[2]= 0.0470362
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
