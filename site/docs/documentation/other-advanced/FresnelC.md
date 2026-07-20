# FresnelC

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FresnelC[z]
    gives the Fresnel integral C(z) = Integral_0^z Cos[Pi t^2/2] dt.
An entire, odd function with no branch cuts. FresnelC[0] = 0,
FresnelC[+-Infinity] = +-1/2, FresnelC[+-I Infinity] = +-I/2.
Real and complex inputs evaluate numerically at machine or arbitrary (MPFR)
precision; D[FresnelC[z], z] = Cos[Pi z^2/2]. Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
