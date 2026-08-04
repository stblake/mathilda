# ToPackedArray

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ToPackedArray[list] is ToNDArray[list]: it returns list stored as a dense machine-precision buffer. The result is still a List -- same Head, same printed form, same elements -- but NDArrayQ gives True for it. Provided under Mathematica's name for the same operation; see ToNDArray.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ToPackedArray[{1., 2., 3.}] === ToNDArray[{1., 2., 3.}]
Out[1]= True
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/pack.c`](https://github.com/stblake/mathilda/blob/main/src/pack.c)
- Specification: [`docs/spec/builtins/packed-arrays.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/packed-arrays.md)
