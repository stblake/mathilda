# ToNDArray

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ToNDArray[list] returns list stored as a dense machine-precision buffer. The result is still a List -- same Head, same printed form, same elements -- but NDArrayQ gives True for it. ToNDArray[list, DataType -> "float64"] forces the element type. A mix of Integer and Real machine values is widened to a Real buffer (the integers become doubles); an all-Integer list stays Integer. Returns list unchanged when it is not rectangular, is empty, or holds any non-machine value. Unlike automatic packing it ignores the size threshold.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= a = ToPackedArray[{1, 2, 3.}]
Out[1]= {1.0, 2.0, 3.0}

In[2]:= {PackedArrayQ[a], DataType[a]}
Out[2]= {True, "float64"}
```

### Scope (3)

```mathematica
In[3]:= DataType[ToNDArray[{True, False, True}]]
Out[3]= "bool"
```

A numeric buffer -> a bool one

```mathematica
In[4]:= Positive[ToNDArray[{-1, 0, 2}]]
Out[4]= {False, False, True}
```

Not numeric: delists to symbolic

```mathematica
In[5]:= Sin[ToNDArray[{True, False}]]
Out[5]= {Sin[True], Sin[False]}
```

### Options (2)

```mathematica
In[6]:= DataType[ToNDArray[{1, 2, 3}, DataType -> "float64"]]
Out[6]= "float64"

In[7]:= NDArrayQ[ToNDArray[{1., 2.5}, DataType -> "int64"]]
Out[7]= False
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Rational](../../arithmetic/Rational/), [Complex](../../arithmetic/Complex/), [Real](../../other-advanced/Real/), [ToPackedArray](../../packed-arrays/ToPackedArray/), [DataType](../../other-advanced/DataType/), [List](../../other-advanced/List/)

- Source: [`src/pack.c`](https://github.com/stblake/mathilda/blob/main/src/pack.c)
- Specification: [`docs/spec/builtins/packed-arrays.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/packed-arrays.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
- Tests: [`tests/test_ndarray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
