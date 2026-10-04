### Worked examples

```mathematica
In[1]:= ToNDArray[{1., 2., 3., 4.}]  (* a real list stored as a dense buffer, still printed as a List *)
```

```mathematica
In[1]:= NDArrayQ[ToNDArray[{1., 2., 3., 4.}]]  (* only NDArrayQ reveals that it is packed *)
```

```mathematica
In[1]:= DataType[ToNDArray[{1, 2, 3, 4}]]  (* an all-integer list packs to an int64 buffer *)
```

```mathematica
In[1]:= ToNDArray[{1, 2, 3.}]  (* a mixed integer/real list widens to a float buffer on request *)
```

```mathematica
In[1]:= ToNDArray[{1., 2.5}, DataType -> "int64"]  (* an int64 request that would round reals is refused; the list is returned unchanged *)
```

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
