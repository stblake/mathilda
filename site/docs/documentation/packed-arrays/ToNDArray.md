# ToNDArray

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ToNDArray[list] returns list stored as a dense machine-precision buffer. The result is still a List -- same Head, same printed form, same elements -- but NDArrayQ gives True for it. ToNDArray[list, DataType -> "float64"] forces the element type. Returns list unchanged when it is not rectangular, is empty, or holds anything other than uniformly Integer or uniformly Real machine values. Unlike automatic packing it ignores the size threshold.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= DataType[ToNDArray[{1, 2, 3}, DataType -> "float64"]]
Out[1]= "float64"

In[2]:= NDArrayQ[ToNDArray[{1., 2.5}, DataType -> "int64"]]
Out[2]= False
```

```mathematica
In[1]:= DataType[ToNDArray[{True, False, True}]]
Out[1]= "bool"

In[2]:= Positive[ToNDArray[{-1, 0, 2}]]      (* a numeric buffer -> a bool one *)
Out[2]= {False, False, True}

In[3]:= Sin[ToNDArray[{True, False}]]        (* not numeric: delists to symbolic *)
Out[3]= {Sin[True], Sin[False]}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/pack.c`](https://github.com/stblake/mathilda/blob/main/src/pack.c)
- Specification: [`docs/spec/builtins/packed-arrays.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/packed-arrays.md)
