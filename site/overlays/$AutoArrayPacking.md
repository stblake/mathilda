### Worked examples

```mathematica
In[1]:= $AutoArrayPacking  (* True by default -- automatic packing is on *)
```

```mathematica
In[1]:= NDArrayQ[Range[10000]]  (* a big machine-number list packs, so this is True *)
```

```mathematica
In[1]:= $AutoArrayPacking = False; NDArrayQ[Range[10000]]  (* with packing off, the same list is a plain List *)
In[2]:= $AutoArrayPacking = True  (* restore the default *)
```

### Notes

`$AutoArrayPacking` is a boolean system variable: `True` (the default) lets Mathilda store
large lists of machine numbers as dense packed buffers; `False` builds every list one
element at a time. It is a read/write OwnValue tied to the real packing flag, and only
`True` or `False` is accepted — anything else draws a `$AutoArrayPacking::flagset` message
and is rolled back.

It changes storage and speed, not answers: a packed list is an ordinary `List` — same head,
printed form, elements, ordering and pattern matches — and only `NDArrayQ` distinguishes
the two, which is why toggling it flips the `NDArrayQ` results above. It does not affect the
explicit `ToNDArray`/`ToPackedArray` or the visible `NDArray[...]` head, and it reads back
`False` in a session started with `MATHILDA_NO_PACK` set.
