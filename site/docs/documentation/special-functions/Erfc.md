# Erfc

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Erfc[z]
    gives the complementary error function erfc(z) = 1 - erf(z).
Erfc[0] = 1, Erfc[Infinity] = 0, Erfc[-Infinity] = 2. An entire
function. Real inputs evaluate via libm/MPFR erfc (cancellation-free);
complex inputs via 1 - erf(z) at machine or arbitrary (MPFR) precision.
D[Erfc[z], z] = -(2/Sqrt[Pi]) E^(-z^2). Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Exact special values: `Erfc[0] = 1`, `Erfc[Infinity] = 0`,

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
