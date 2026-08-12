# FromNDArray

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FromNDArray[expr] returns expr with any dense buffer storage undone: a packed List becomes an ordinary List of separate elements, and an NDArray[...] becomes the nested List of its entries. Anything else is returned unchanged. Inverse of ToNDArray.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= NDArrayQ[FromNDArray[ToNDArray[{1., 2., 3.}]]]
Out[1]= False

In[2]:= FromNDArray[NDArray[{1., 2.}]]
Out[2]= {1.0, 2.0}
```

## Implementation notes

**Attributes:** `Protected`.

## See also

[Normal](../../power-series/Normal/)

## References

- Source: [`src/pack.c`](https://github.com/stblake/mathilda/blob/main/src/pack.c)
- Specification: [`docs/spec/builtins/packed-arrays.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/packed-arrays.md)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
