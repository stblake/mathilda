# $AutoArrayPacking

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
$AutoArrayPacking
    controls whether Mathilda stores large lists of machine numbers as
    dense buffers (packed arrays). True by default; set it to False to
    build every list one element at a time.

A packed list is an ordinary List -- same Head, printed form, elements,
ordering and pattern matches -- and only NDArrayQ tells the two apart.
So this changes storage and speed, not answers. Does not affect
ToNDArray or ToPackedArray, which are explicit requests, nor the
explicit NDArray[...] head.

Reads back False in a session started with the environment variable
MATHILDA_NO_PACK set. Only True or False is accepted.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
