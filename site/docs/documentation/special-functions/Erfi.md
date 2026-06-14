# Erfi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Erfi[z]
    gives the imaginary error function erfi(z) = -I Erf[I z]
    = (2/Sqrt[Pi]) Integral_0^z e^(t^2) dt.
Erfi[0] = 0, Erfi[Infinity] = Infinity, Erfi[I Infinity] = I. An entire
function, odd in z. Real and complex inputs evaluate numerically at
machine or arbitrary (MPFR) precision; D[Erfi[z], z] = (2/Sqrt[Pi]) E^(z^2).
Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Exact special values: `Erfi[0] = 0`, `Erfi[Infinity] = Infinity`,

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
