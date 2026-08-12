# ToNDArray

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ToNDArray[list] returns list stored as a dense machine-precision buffer. The result is still a List -- same Head, same printed form, same elements -- but NDArrayQ gives True for it. ToNDArray[list, DataType -> "float64"] forces the element type. Returns list unchanged when it is not rectangular, is empty, or holds anything other than uniformly Integer or uniformly Real machine values. Unlike automatic packing it ignores the size threshold.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Scope (3)

```mathematica
In[1]:= DataType[ToNDArray[{True, False, True}]]
Out[1]= "bool"

In[2]:= Positive[ToNDArray[{-1, 0, 2}]]      (* a numeric buffer -> a bool one *)
Out[2]= {False, False, True}

In[3]:= Sin[ToNDArray[{True, False}]]        (* not numeric: delists to symbolic *)
Out[3]= {Sin[True], Sin[False]}
```

### Options (2)

```mathematica
In[4]:= DataType[ToNDArray[{1, 2, 3}, DataType -> "float64"]]
Out[4]= "float64"

In[5]:= NDArrayQ[ToNDArray[{1., 2.5}, DataType -> "int64"]]
Out[5]= False
```

## Implementation notes

**Attributes:** `Protected`.

## See also

[Rational](../../arithmetic/Rational/), [Complex](../../arithmetic/Complex/), [DataType](../../other-advanced/DataType/), [List](../../other-advanced/List/)

## References

- Source: [`src/pack.c`](https://github.com/stblake/mathilda/blob/main/src/pack.c)
- Specification: [`docs/spec/builtins/packed-arrays.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/packed-arrays.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
- Tests: [`tests/test_ndarray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
