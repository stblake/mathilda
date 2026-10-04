### Worked examples

```mathematica
In[1]:= NDArrayQ[NDArray[{1, 2, 3}]]  (* a visible NDArray object *)
```

```mathematica
In[1]:= NDArrayQ[Range[10000]]  (* a big machine-number list packs transparently *)
```

```mathematica
In[1]:= NDArrayQ[{1, 2, 3}]  (* a small plain list was never packed *)
```

```mathematica
In[1]:= NDArrayQ[5]  (* a scalar is not an array *)
```

### Notes

`NDArrayQ[expr]` gives `True` when `expr` is an NDArray object, and `False` otherwise. It is
an `O(1)` representation test that recognises **both** packed-array surfaces: the visible
`NDArray[...]` head and a transparently packed `List` (the dense machine-number buffer
Mathilda builds behind an ordinary list head).

That dual recognition is what makes it the definitive way to see through the packing gate: a
packed list looks exactly like a plain `List` in every other respect, so only `NDArrayQ`
confirms a value is on the buffer fast path — hence `Range[10000]` is `True` while the small
`{1, 2, 3}` is `False`. Its sibling `PackedArrayQ` tests only the packed-`List` surface, so
it would answer `False` on a visible `NDArray[...]`.
