# NDArray

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NDArray[nested_list]
    Packs a rectangular, machine-precision (Integer/Real) nested
    list into a dense N-dimensional array (numpy ndarray style).
    Visibly distinct from List: Head, ListQ, and printing never
    treat an NDArray as a List. Dimensions gives its shape,
    ArrayDepth its rank, Length its leading-axis length. Builtins
    that recognize NDArray (Dot, Plus, Times) use a fast C-level
    path; results that would need a non-machine-precision entry
    auto-degrade to an ordinary nested List.
NDArray[nested_list, DataType -> "float32"]
    Packs at the given element type: "float64" (default), "float32",
    "complex64", or "complex32". DataType[a] gives an array's type.
    A ragged (non-rectangular) list is rejected with an
    NDArray::ragged warning; an empty or non-machine-precision
    list stays unevaluated.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= NDArray[{{1, 2}, {3, 4}}]
Out[1]= NDArray[{{1.0, 2.0}, {3.0, 4.0}}]

In[2]:= Dimensions[NDArray[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]]
Out[2]= {2, 2, 2}

In[3]:= Depth[NDArray[{{1, 2}, {3, 4}}]]
Out[3]= 3

In[4]:= Dot[NDArray[{{1, 2}, {3, 4}}], NDArray[{{5, 6}, {7, 8}}]]
Out[4]= NDArray[{{19.0, 22.0}, {43.0, 50.0}}]

In[5]:= NDArray[{{1, 2}, {3, 4}}] + NDArray[{{5, 6}, {7, 8}}]
Out[5]= NDArray[{{6.0, 8.0}, {10.0, 12.0}}]

In[6]:= NDArrayQ[NDArray[{1, 2, 3}]]
Out[6]= True

In[7]:= NDArray[{{1, x}, {3, 4}}]
Out[7]= NDArray[{{1, x}, {3, 4}}]
```

## Implementation notes

- `Protected`.
- A ragged (non-rectangular) `list` — unequal sublist shapes, or a mix of

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/ndarray.c`](https://github.com/stblake/mathilda/blob/main/src/ndarray.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
