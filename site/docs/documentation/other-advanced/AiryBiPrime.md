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
