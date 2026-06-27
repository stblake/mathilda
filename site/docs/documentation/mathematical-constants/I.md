# I

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
I
    is the imaginary unit Sqrt[-1].
I represents the imaginary unit; I^2 evaluates to -1 and complex numbers
are written a + b I. It has attribute Protected, and N[I] is 0. + 1. I.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Attribute `Protected`. `Attributes[I] = {Protected}`; the symbol cannot be

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)
