# FromPackedArray

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FromPackedArray[expr] is FromNDArray[expr]: it returns expr with any dense buffer storage undone, so a packed List becomes an ordinary List of separate elements and an NDArray[...] becomes the nested List of its entries. Provided under Mathematica's name for the same operation; see FromNDArray.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= FromPackedArray[ToPackedArray[{1., 2., 3., 4.}]]
Out[1]= {1.0, 2.0, 3.0, 4.0}

In[2]:= NDArrayQ[FromPackedArray[ToPackedArray[{1., 2., 3., 4.}]]]
Out[2]= False

In[3]:= FromPackedArray[NDArray[{1., 2., 3.}]]
Out[3]= {1.0, 2.0, 3.0}
```

### Applications (3)

Mathematica's name; undoes the buffer

```mathematica
In[4]:= FromPackedArray[ToPackedArray[{1., 2., 3., 4.}]]
Out[4]= {1.0, 2.0, 3.0, 4.0}
```

The result is an ordinary List

```mathematica
In[5]:= NDArrayQ[FromPackedArray[ToPackedArray[{1., 2., 3., 4.}]]]
Out[5]= False
```

Unpacks a visible NDArray too

```mathematica
In[6]:= FromPackedArray[NDArray[{1., 2., 3.}]]
Out[6]= {1.0, 2.0, 3.0}
```

## Implementation notes

**Algorithm.** `FromPackedArray` is Mathematica's name (`Developer`​`FromPackedArray`) for
`FromNDArray` and is **the same C builtin** (`builtin_fromndarray`) registered under a second
name — an aliased builtin, not a DownValue that rewrites to the other, so it costs no extra
evaluation pass and cannot be shadowed. It undoes both forms of buffer storage: a packed
`List` becomes an ordinary `List` of separate element nodes (`NDArrayQ` then `False`), and a
visible `NDArray[...]` becomes the nested `List` of its entries (`ndarray_to_nested_list`).
Anything that is not an `EXPR_NDARRAY` is returned unchanged. It is the inverse of
`ToPackedArray`.

**Data structures.** Reads the dense row-major buffer and allocates one boxed `Expr` leaf per
element (`Integer`/`Real`/`True`/`False` by dtype), rebuilding the nested `List` from the
stored `int64` dims.

**Complexity / limits.** `O(n)`, allocating an expression node per element — the per-element
cost packing avoids, so this is where a large array stops being cheap to hold. The value is
unchanged; only the representation differs.

**Attributes:** `Protected`.

## References

**See also:** [FromNDArray](../../packed-arrays/FromNDArray/), [ToPackedArray](../../packed-arrays/ToPackedArray/), [List](../../other-advanced/List/)

- Source: [`src/pack.c`](https://github.com/stblake/mathilda/blob/main/src/pack.c)
- Specification: [`docs/spec/builtins/packed-arrays.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/packed-arrays.md)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

`FromPackedArray` is the Wolfram Language's name (`Developer`​`FromPackedArray`) for
`FromNDArray`, and is the **same builtin registered twice** — not a rule that rewrites to the
other — so it costs no extra evaluation pass and cannot be shadowed. It undoes both forms of
buffer storage: a packed `List` becomes an ordinary `List` of separate elements (so `NDArrayQ`
is then `False`), and a visible `NDArray[...]` becomes the nested `List` of its entries.
Anything else is returned unchanged. It is the inverse of `ToPackedArray`.

The value is unchanged; only the representation differs. Like `FromNDArray`, it allocates one
expression node per element — the per-element cost that packing avoids.
