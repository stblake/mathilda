### Worked examples

```mathematica
In[1]:= FromPackedArray[ToPackedArray[{1., 2., 3., 4.}]]  (* Mathematica's name; undoes the buffer *)
```

```mathematica
In[1]:= NDArrayQ[FromPackedArray[ToPackedArray[{1., 2., 3., 4.}]]]  (* the result is an ordinary List *)
```

```mathematica
In[1]:= FromPackedArray[NDArray[{1., 2., 3.}]]  (* unpacks a visible NDArray too *)
```

### Notes

`FromPackedArray` is the Wolfram Language's name (`Developer`​`FromPackedArray`) for
`FromNDArray`, and is the **same builtin registered twice** — not a rule that rewrites to the
other — so it costs no extra evaluation pass and cannot be shadowed. It undoes both forms of
buffer storage: a packed `List` becomes an ordinary `List` of separate elements (so `NDArrayQ`
is then `False`), and a visible `NDArray[...]` becomes the nested `List` of its entries.
Anything else is returned unchanged. It is the inverse of `ToPackedArray`.

The value is unchanged; only the representation differs. Like `FromNDArray`, it allocates one
expression node per element — the per-element cost that packing avoids.
