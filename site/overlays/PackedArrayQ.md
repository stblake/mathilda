### Worked examples

```mathematica
In[1]:= PackedArrayQ[ToPackedArray[{1., 2., 3.}]]  (* a packed List *)
```

```mathematica
In[1]:= PackedArrayQ[NDArray[{1., 2., 3.}]]  (* a visible NDArray is not a packed List *)
```

```mathematica
In[1]:= PackedArrayQ[{a, b, c}]  (* an ordinary unpacked List *)
```

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
