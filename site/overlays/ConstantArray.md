### Worked examples

```mathematica
In[1]:= ConstantArray[Pi, 4]  (* any expression may be the repeated element *)
```

```mathematica
In[1]:= ConstantArray[1, {2, 2}]  (* a nested rectangular array *)
```

```mathematica
In[1]:= ConstantArray[x, {2, 3}]  (* a symbolic element is copied verbatim *)
```

### Notes

`ConstantArray[c, n]` is a flat list of `n` copies of `c`; `ConstantArray[c, {n1,
..., nk}]` is an `n1 x ... x nk` nested array. The element `c` is copied
verbatim, so it may be a symbol, a number, or a compound structure such as a
matrix. A dimension of `0` yields an empty list at that level; dimensions must be
non-negative machine integers.

When `c` is a machine number over a rectangular shape, the whole result is a
packed array written straight into a buffer, so building a large constant array
costs a memory fill rather than one boxed element per cell.
