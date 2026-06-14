# InverseErfc

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
InverseErfc[s]
    gives the inverse complementary error function: the z solving s = Erfc[z].
InverseErfc[0] = Infinity, InverseErfc[1] = 0, InverseErfc[2] = -Infinity.
Numerical values are given only for real s in [0, 2], at machine or
arbitrary (MPFR) precision; D[InverseErfc[z], z] =
-(Sqrt[Pi]/2) E^(InverseErfc[z]^2). Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Exact special values: `InverseErfc[0] = Infinity`, `InverseErfc[1] = 0`,

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
