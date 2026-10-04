# ToNDArray

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ToNDArray[list] returns list stored as a dense machine-precision buffer. The result is still a List -- same Head, same printed form, same elements -- but NDArrayQ gives True for it. ToNDArray[list, DataType -> "float64"] forces the element type. A mix of Integer and Real machine values is widened to a Real buffer (the integers become doubles); an all-Integer list stays Integer. Returns list unchanged when it is not rectangular, is empty, or holds any non-machine value. Unlike automatic packing it ignores the size threshold.`**

## Examples (12)

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

### Applications (5)

A real list stored as a dense buffer, still printed as a List

```mathematica
In[8]:= ToNDArray[{1., 2., 3., 4.}]
Out[8]= {1.0, 2.0, 3.0, 4.0}
```

Only NDArrayQ reveals that it is packed

```mathematica
In[9]:= NDArrayQ[ToNDArray[{1., 2., 3., 4.}]]
Out[9]= True
```

An all-integer list packs to an int64 buffer

```mathematica
In[10]:= DataType[ToNDArray[{1, 2, 3, 4}]]
Out[10]= "int64"
```

A mixed integer/real list widens to a float buffer on request

```mathematica
In[11]:= ToNDArray[{1, 2, 3.}]
Out[11]= {1.0, 2.0, 3.0}
```

An int64 request that would round reals is refused; the list is returned unchanged

```mathematica
In[12]:= ToNDArray[{1., 2.5}, DataType -> "int64"]
Out[12]= {1.0, 2.5}
```

## Implementation notes

**Algorithm.** `builtin_tondarray` first strips an optional trailing `DataType -> "..."` rule
with `pack_take_dtype` (rightmost wins); after that there must be exactly one positional
argument. An argument that is already a packed `List` with no conflicting `DataType` is
returned by `expr_copy` (a no-op); a visible `NDArray[...]` is turned back into a nested
`List` first (`ndarray_to_nested_list`), so `ToNDArray` is also a way to say "same values, as
a `List`". The list is then packed by `pack_force_coerce` → `pack_build` with `min_elems = 0`
(no size threshold, unlike automatic packing) and `coerce = true`:

1. `pack_sniff` walks the nested list once, checking rectangularity and folding every leaf's
   class (`PK_INT`/`PK_REAL`/`PK_BOOL`) into one; with `coerce = true` a mix of machine
   `Integer` and `Real` folds to `PK_REAL` (the integers will widen to `double`), whereas the
   automatic path (`coerce = false`) declines that mix to keep `1 === 1.` observable.
2. The dtype is `int64` (all-integer), `float64` (real, or coerced mixed), or `bool`
   (all-`True`/`False`). An explicit `DataType` may *widen* an exact list to a float buffer
   but is refused when it would round a `Real` into an `int64` slot (`want == NDT_INT64 &&
   cls != PK_INT`), request a complex buffer, or mismatch bool-vs-number — in which case the
   original list is returned unchanged.
3. `pack_flatten` writes the leaves row-major into a freshly `malloc`'d buffer.

A list that is ragged, empty, or holds any `Rational`/`BigInt`/arbitrary-precision/`Complex`/
symbolic element simply comes back unchanged (never an error).

**Data structures.** The result is an `EXPR_NDARRAY` with `present_as = NDA_HEAD_LIST` — the
*packed-list* surface, so `Head` is still `List` and only `NDArrayQ`/`PackedArrayQ` reveal the
buffer. Storage is a dense row-major buffer of `ndt_elem_size(dt) * n` bytes carrying rank and
`int64` dims (up to `NDARRAY_MAX_RANK`).

**Complexity / limits.** `O(n)` in two passes (sniff then flatten) and one allocation. Ignores
the automatic-packing size threshold. Complex dtypes are rejected (no faithful round trip
yet); the packing contract is representation-only, so every downstream operation either
answers exactly as the ordinary list does or falls back to it.

**Attributes:** `Protected`.

## References

**See also:** [Rational](../../arithmetic/Rational/), [Complex](../../arithmetic/Complex/), [Real](../../other-advanced/Real/), [ToPackedArray](../../packed-arrays/ToPackedArray/), [DataType](../../other-advanced/DataType/), [List](../../other-advanced/List/)

- Source: [`src/pack.c`](https://github.com/stblake/mathilda/blob/main/src/pack.c)
- Specification: [`docs/spec/builtins/packed-arrays.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/packed-arrays.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
- Tests: [`tests/test_ndarray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)

## Notes & additional examples

### Notes

`ToNDArray[list]` returns `list` stored as a dense machine-precision buffer. The result is
still a `List` — same `Head`, same printed form, same elements — so the packing is invisible
except to `NDArrayQ` / `PackedArrayQ`. Unlike automatic packing it ignores the size threshold,
so even a tiny list packs.

The element type is inferred: all-`Integer` → `"int64"`, all-`Real` → `"float64"`,
all-`True`/`False` → `"bool"`. A list that **mixes** machine `Integer` and `Real` widens to a
`"float64"` buffer (the integers become doubles) — the one point where the explicit request
diverges from automatic packing, which leaves `{1, 2, 3.}` an ordinary list because turning
`1` into `1.` would be an observable head change. `ToNDArray[list, DataType -> "..."]` forces
the type, and may widen an exact list to a float buffer, but never rounds a `Real` into an
`int64` slot (that request returns the list unchanged).

It returns `list` unchanged — not an error — whenever it cannot pack: a ragged (non-
rectangular) list, an empty list, or any `Rational`/`BigInt`/arbitrary-precision/`Complex`/
symbolic element. On an already-packed list it is a no-op; on a visible `NDArray[...]` it
restates the same values as a `List`.
