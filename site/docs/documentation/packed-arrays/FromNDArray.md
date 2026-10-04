# FromNDArray

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FromNDArray[expr] returns expr with any dense buffer storage undone: a packed List becomes an ordinary List of separate elements, and an NDArray[...] becomes the nested List of its entries. Anything else is returned unchanged. Inverse of ToNDArray.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= NDArrayQ[FromNDArray[ToNDArray[{1., 2., 3.}]]]
Out[1]= False

In[2]:= FromNDArray[NDArray[{1., 2.}]]
Out[2]= {1.0, 2.0}
```

### Applications (3)

Undoes the buffer: an ordinary List of elements

```mathematica
In[3]:= FromNDArray[ToNDArray[{1., 2., 3.}]]
Out[3]= {1.0, 2.0, 3.0}
```

No longer packed

```mathematica
In[4]:= NDArrayQ[FromNDArray[ToNDArray[{1., 2., 3.}]]]
Out[4]= False
```

Also unpacks a visible NDArray into a nested List

```mathematica
In[5]:= FromNDArray[NDArray[{1., 2.}]]
Out[5]= {1.0, 2.0}
```

## Implementation notes

**Algorithm.** `builtin_fromndarray` undoes buffer storage. It takes one argument: anything
that is not an `EXPR_NDARRAY` (as reported by `is_ndarray`, which is true for *both* the
packed-`List` surface and a visible `NDArray[...]`) is returned unchanged by `expr_copy`;
otherwise the buffer is expanded into a nested `List` of separate element nodes by
`ndarray_to_nested_list`. So a packed `List` becomes an ordinary `List` (`NDArrayQ` then
`False`), and an `NDArray[...]` becomes the nested `List` of its entries. It is the inverse of
`ToNDArray` and does the same job as `Normal` on these two forms.

**Data structures.** Reads the dense row-major buffer and rebuilds one boxed `Expr` leaf per
element (`Integer` from an `int64` buffer, `Real` from `float64`, `True`/`False` from `bool`),
reassembling the nested `List` from the stored `int64` dims.

**Complexity / limits.** `O(n)` — it allocates an expression node per element, which is
exactly the per-element cost that packing exists to avoid, so this is the point at which a
large array stops being cheap to hold. Pure representation change: the values are identical to
the buffer's.

**Attributes:** `Protected`.

## References

**See also:** [Normal](../../data-structures/Normal/)

- Source: [`src/pack.c`](https://github.com/stblake/mathilda/blob/main/src/pack.c)
- Specification: [`docs/spec/builtins/packed-arrays.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/packed-arrays.md)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

`FromNDArray[expr]` undoes buffer storage. A packed `List` becomes an ordinary `List` of
separate element nodes (so `NDArrayQ` is then `False`), and an `NDArray[...]` becomes the
nested `List` of its entries. Anything else is returned unchanged. It is the inverse of
`ToNDArray`, and does the same job as `Normal` on these two forms.

The value is identical to the buffer's — only the representation changes — but it allocates one
expression node per element, which is exactly the per-element cost packing exists to avoid, so
this is the point at which a large array stops being cheap to hold.
