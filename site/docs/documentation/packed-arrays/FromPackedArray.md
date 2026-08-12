# FromPackedArray

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FromPackedArray[expr] is FromNDArray[expr]: it returns expr with any dense buffer storage undone, so a packed List becomes an ordinary List of separate elements and an NDArray[...] becomes the nested List of its entries. Provided under Mathematica's name for the same operation; see FromNDArray.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/pack.c`](https://github.com/stblake/mathilda/blob/main/src/pack.c)
- Specification: [`docs/spec/builtins/packed-arrays.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/packed-arrays.md)
