# SinhIntegral

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SinhIntegral[z]
    gives the hyperbolic sine integral Shi(z) = Integral_0^z Sinh[t]/t dt.
An entire, odd function with no branch cuts. SinhIntegral[0] = 0,
SinhIntegral[+-Infinity] = +-Infinity, SinhIntegral[+-I Infinity] = +-I Pi/2.
Real and complex inputs evaluate numerically at machine or arbitrary (MPFR)
precision; D[SinhIntegral[z], z] = Sinh[z]/z. Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
