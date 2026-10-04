### Worked examples

```mathematica
In[1]:= PadLeft[{a, b, c}, 5]  (* pad with zeros on the left to length 5 *)
```

```mathematica
In[2]:= PadLeft[{a, b, c}, 10, {x, y, z}]  (* a list padding is tiled cyclically *)
```

```mathematica
In[3]:= PadLeft[{1, 2, 3}, 2]  (* a shorter length keeps the trailing elements *)
```

```mathematica
In[4]:= PadLeft[{{a, b}, {c}}, {3, 4}]  (* a dimension list builds a full nested array *)
```

### Notes

`PadLeft[list, n]` returns a length-`n` list padded on the left with `0`;
`PadLeft[list, n, x]` repeats a given element, and `PadLeft[list, n, {x1, ...}]`
tiles a list of pad elements cyclically so the padding phase is continuous. A
length shorter than the list keeps its *last* `n` elements (left padding, so the
front is dropped); a negative length pads on the opposite side. The dimension-list
form `PadLeft[list, {n1, n2, ...}]` builds a full nested array with length `ni` at
level `i`, and the bare `PadLeft[list]` pads a ragged array to rectangular. The
head of `list` need not be `List`. A rank-1 packed/`NDArray` buffer takes a native
fast path; there is no `Compile[]` lowering.
