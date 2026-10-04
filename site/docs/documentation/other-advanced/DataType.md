# DataType

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DataType[a]`**

gives the element data type of the NDArray a as a string: one of "float64", "float32", "complex64", or "complex32". Set the type when constructing with NDArray\[list, DataType -\> "float32"\]; the four types map onto BLAS's s/d/c/z precisions. Returns unevaluated for a non-NDArray argument.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

The default packed element type

```mathematica
In[1]:= DataType[NDArray[{1.0, 2.0, 3.0}]]
Out[1]= "float64"
```

DataType is also the option keyword that set it

```mathematica
In[2]:= DataType[NDArray[{1, 2, 3}, DataType -> "float32"]]
Out[2]= "float32"
```

A plain list has no dtype, so this stays symbolic

```mathematica
In[3]:= DataType[{1, 2, 3}]
Out[3]= DataType[{1, 2, 3}]
```

DataType is the one NDArray option, defaulting to float64

```mathematica
In[4]:= Options[NDArray]
Out[4]= {DataType -> "float64"}
```

## Implementation notes

**Algorithm.** `DataType` plays two roles, both wired in `ndarray.c`.

As a **function**, `builtin_datatype` reads the element type of a packed NDArray.
With one argument it checks `is_ndarray(arg)`; on a hit it returns
`ndt_to_string(arg->data.ndarray.dtype)` — one of the strings `"float64"`,
`"float32"`, `"complex64"`, `"complex32"`, or `"bool"`. On anything that is not an
NDArray (an ordinary `List`, a scalar) it returns `NULL`, so `DataType[...]` stays
symbolic rather than guessing a type.

As an **option keyword**, `DataType` is the name on the left of the
`DataType -> "float32"` rule that `NDArray[list, DataType -> "..."]` reads when
packing, and it is the single entry in `Options[NDArray]` (default
`DataType -> "float64"`). It is registered `Protected`; its docstring lives
centrally in `info.c`.

**Data structures.** The dtype is a small enum tag stored in the NDArray's
`ndarray` struct (it selects the buffer's element width and the BLAS s/d/c/z
precision). The function emits a freshly allocated `EXPR_STRING`.

**Complexity / limits.** `O(1)` — a tag read, no buffer traversal. Only the packed
`NDArray[...]` surface carries a dtype; a plain list has none, which is why
`DataType` declines on one.

**Attributes:** `Protected`.

## References

- Source: [`src/ndarray.c`](https://github.com/stblake/mathilda/blob/main/src/ndarray.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_fourier.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fourier.c)
- Tests: [`tests/test_ndarray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)

## Notes & additional examples

### Notes

`DataType` is both a reader and an option keyword. As a function, `DataType[a]`
gives the element type of a packed `NDArray` `a` as a string — `"float64"`,
`"float32"`, `"complex64"`, `"complex32"` or `"bool"`; on a non-array it declines
and stays symbolic. As an option, `DataType -> "float32"` is how `NDArray[list,
...]` chooses the packed element type, and it is the sole entry of
`Options[NDArray]`. The five type names map onto BLAS's s/d/c/z precisions.
