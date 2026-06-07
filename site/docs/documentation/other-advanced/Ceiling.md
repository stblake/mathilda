# Ceiling

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Ceiling[x]
    gives the smallest integer greater than or equal to x.
Ceiling[x, a]
    gives the smallest multiple of a greater than or equal to x.
Ceiling is Listable. Exact inputs return exact integers; Real / MPFR
inputs are rounded toward +Infinity at the input precision.
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
