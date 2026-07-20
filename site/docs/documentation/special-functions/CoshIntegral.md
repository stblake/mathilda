# CoshIntegral

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
CoshIntegral[z]
    gives the hyperbolic cosine integral
    Chi(z) = EulerGamma + Log[z] + Integral_0^z (Cosh[t] - 1)/t dt.
Has a logarithmic singularity at 0 and a branch cut on (-Infinity, 0].
CoshIntegral[0] = -Infinity, CoshIntegral[Infinity] = Infinity,
CoshIntegral[+-I Infinity] = +-I Pi/2.
Real and complex inputs evaluate numerically at machine or arbitrary (MPFR)
precision; D[CoshIntegral[z], z] = Cosh[z]/z. Listable.
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
