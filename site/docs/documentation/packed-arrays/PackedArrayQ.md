# PackedArrayQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PackedArrayQ[expr]`**

Gives True if expr is a packed array -- a List stored as a dense machine-precision buffer -- else False. Provided under Mathematica's name; unlike NDArrayQ, a visible NDArray\[...\] object gives False.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= PackedArrayQ[ToPackedArray[{1., 2., 3.}]]
Out[1]= True

In[2]:= PackedArrayQ[NDArray[{1., 2., 3.}]]
Out[2]= False

In[3]:= NDArrayQ[NDArray[{1., 2., 3.}]]
Out[3]= True

In[4]:= PackedArrayQ[{a, b, c}]
Out[4]= False
```

### Applications (3)

A packed List

```mathematica
In[5]:= PackedArrayQ[ToPackedArray[{1., 2., 3.}]]
Out[5]= True
```

A visible NDArray is not a packed List

```mathematica
In[6]:= PackedArrayQ[NDArray[{1., 2., 3.}]]
Out[6]= False
```

An ordinary unpacked List

```mathematica
In[7]:= PackedArrayQ[{a, b, c}]
Out[7]= False
```

## Implementation notes

**Algorithm.** `builtin_packedarrayq` is a one-argument predicate that returns the symbol
`True` when its argument is a packed array — a `List` the system stores as a dense
machine-precision buffer — and `False` otherwise. The test is exactly `is_packed_list`, so it
is deliberately **narrower** than `NDArrayQ`: `NDArrayQ` is `True` for either surface of an
`EXPR_NDARRAY` (the packed-`List` surface *or* a visible `NDArray[...]` object), whereas
`PackedArrayQ` is `True` only for the packed-`List` surface. A visible `NDArray[...]` is a
distinct atom (`AtomQ` true, `ListQ` false) rather than a `List` that happens to be packed, so
it gives `False` — matching Mathematica (`Developer`​`PackedArrayQ`), which has no visible
`NDArray` head at all. The two predicates differ only on that one input.

For the answer to be correct the head must itself be **packed-aware**: `PackedArrayQ` is on
both the `AWARE` and `INT64_OK` lists in `src/pack.c`. Otherwise the evaluator's transparency
gate would materialise a packed argument into a plain `List` before the builtin runs, and it
could only ever observe an already-unpacked list and answer `False` for every packed input.

**Data structures.** None of its own — it inspects the argument's tag (`EXPR_NDARRAY` with the
packed-`List` presentation) and allocates only the returned `True`/`False` symbol.

**Complexity / limits.** `O(1)`. `Protected`.

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [NDArrayQ](../../other-advanced/NDArrayQ/), [AtomQ](../../expression-information/AtomQ/), [ListQ](../../expression-information/ListQ/), [NDArray](../../linear-algebra/NDArray/)

- Source: [`src/ndarray.c`](https://github.com/stblake/mathilda/blob/main/src/ndarray.c)
- Specification: [`docs/spec/builtins/packed-arrays.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/packed-arrays.md)
- Tests: [`tests/test_logexp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_logexp.c)

## Notes & additional examples

### Notes

`PackedArrayQ[expr]` gives `True` when `expr` is a packed array — a `List` stored as a dense
machine-precision buffer — and `False` otherwise. It is the Wolfram Language's name
(`Developer`​`PackedArrayQ`) for that test.

It is deliberately **narrower** than `NDArrayQ`. `NDArrayQ` is `True` for either surface of the
internal array object — a packed `List` *or* a visible `NDArray[...]`. `PackedArrayQ` is `True`
only for the packed-`List` surface: a visible `NDArray[...]` is a distinct atom (`AtomQ` is
`True`, `ListQ` is `False`), not a `List` that happens to be packed, so it gives `False` —
matching the Wolfram Language, which has no visible `NDArray` head at all. The two predicates
differ only on that one input.

Like `NDArrayQ`, `PackedArrayQ` is packed-aware, so the evaluator hands it a packed argument
intact (for `int64` buffers as well as `float64`) instead of unpacking it first — without
that it could only ever observe an already-unpacked `List` and would answer `False` for every
packed input. `Protected`.
