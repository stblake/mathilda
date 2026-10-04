### Worked examples

```mathematica
In[1]:= PadRight[{a, b, c}, 5]  (* pad with zeros on the right to length 5 *)
```

```mathematica
In[2]:= PadRight[{{a, b}, {c}}, {3, 5}]  (* a dimension list builds a full rectangular array *)
```

```mathematica
In[3]:= PadRight[{{a, b}, {c}}]  (* Automatic pads a ragged array to rectangular *)
```

```mathematica
In[4]:= PadRight[{a, b, c}, 2]  (* a shorter length keeps the leading elements *)
```

### Notes

`PadRight` is the exact mirror of `PadLeft`, padding on the right instead of the
left. `PadRight[list, n]` pads with `0`; `PadRight[list, n, x]` repeats `x`, and a
list padding is tiled cyclically. A length shorter than the list keeps its
*first* `n` elements; a negative length pads on the left. `PadRight[list, {n1,
n2, ...}]` builds a full nested array, and `PadRight[list]` (or `Automatic`) pads
a ragged array to the smallest enclosing rectangle. The head of `list` need not
be `List`. A rank-1 packed/`NDArray` buffer takes a native fast path; padding an
exact `0` into a float buffer yields a mixed list that cannot repack. There is no
`Compile[]` lowering.
