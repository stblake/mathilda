# DataType

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
DataType[a]
    gives the element data type of the NDArray a as a string: one of
    "float64", "float32", "complex64", or "complex32". Set the type
    when constructing with NDArray[list, DataType -> "float32"]; the
    four types map onto BLAS's s/d/c/z precisions. Returns unevaluated
    for a non-NDArray argument.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
