# Beta

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Beta[a, b]
    is the Euler beta function B(a, b) = Gamma(a) Gamma(b) / Gamma(a+b).
Beta[z, a, b]
    is the incomplete beta function Integral_0^z t^(a-1) (1-t)^(b-1) dt.
Beta[z0, z1, a, b]
    is the generalized incomplete beta Beta[z1, a, b] - Beta[z0, a, b].
Exact for rational arguments (a positive integer gives a rational via
Pochhammer); non-positive integer poles give ComplexInfinity. Machine
and arbitrary-precision (MPFR) real and complex inputs evaluate
numerically. The incomplete form reduces through Hypergeometric2F1.
Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
