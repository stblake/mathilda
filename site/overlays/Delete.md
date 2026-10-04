### Worked examples

```mathematica
In[1]:= Delete[{a, b, c, d}, 2]
```

```mathematica
In[1]:= Delete[{a, b, c, d}, -1]  (* negative indices count from the end *)
```

```mathematica
In[1]:= Delete[{a, b, c, d}, {{1}, {3}}]  (* a list of positions deletes several at once *)
```

```mathematica
In[1]:= Delete[f[a, b, c], 2]  (* works on any head *)
```

### Notes

`Delete[expr, n]` removes the element at position `n`, rebuilding the enclosing
expression without it; negative indices count from the end, and a position path
`{i, j, ...}` reaches a nested element. A list of positions
`{{p1}, {p2}, ...}` deletes each of them in one call. `Delete` works on any head,
not just `List`. Unlike `Drop`, which removes a count or a strided range, `Delete`
targets positions explicitly; out-of-range positions leave the structure
unchanged.
