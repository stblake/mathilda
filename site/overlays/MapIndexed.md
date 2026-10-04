### Worked examples

```mathematica
In[1]:= MapIndexed[f, {a, b, c}]  (* the index arrives as a LIST, {1}, {2}, ... *)
```

```mathematica
In[1]:= MapIndexed[#2 &, {a, b, c}]  (* keep only the position of each element *)
```

```mathematica
In[1]:= MapIndexed[{#1, First[#2]} &, {x, y, z}]  (* pair each element with its integer index *)
```

```mathematica
In[1]:= MapIndexed[f, {{a, b}, {c, d}}, {2}]  (* at level 2 the position is a two-index list *)
```

### Notes

`MapIndexed[f, expr]` is `Map` that also hands `f` the position of each element as
a second argument: `{f[e1, {1}], f[e2, {2}], ...}`. The position is always a list
of indices, so at deeper levels it has one entry per level (`{2}` makes the
positions `{i, j}`). Inside a pure function `#1` is the element and `#2` its
position — `First[#2]` recovers the plain integer index.

A level spec selects which parts are wrapped (default `{1}`), and `Heads -> True`
also indexes heads. Over an association the values are mapped and positioned by
`Key[k]`.
