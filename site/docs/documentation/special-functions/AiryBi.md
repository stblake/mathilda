# AiryBi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
AiryBi[z]
    gives the Airy function Bi(z), the solution of y'' = z y that grows
    exponentially as z -> +Infinity.
AiryBi[0] = 1/(3^(1/6) Gamma[2/3]), AiryBi[Infinity] = Infinity,
AiryBi[-Infinity] = 0. An entire function of z. Real and complex inputs
evaluate numerically at machine or arbitrary (MPFR) precision;
D[AiryBi[z], z] = AiryBiPrime[z]. Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= AiryBi[0]
Out[1]= 1/(3^(1/6) Gamma[2/3])

In[2]:= AiryBi[1.8]
Out[2]= 2.59587
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
