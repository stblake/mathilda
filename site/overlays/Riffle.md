### Worked examples

```mathematica
In[1]:= Riffle[{1, 2, 3, 4}, 0]  (* one separator in every gap *)
```

```mathematica
In[1]:= Riffle[Range[5], x]  (* interleave a symbol between the elements *)
```

```mathematica
In[1]:= Riffle[{a, b, c, d, e}, {1, 2, 3}]  (* a separator list cycles through the gaps *)
```

```mathematica
In[1]:= StringJoin[Riffle[{"a", "b", "c"}, ", "]]  (* the comma-join idiom *)
```

### Notes

`Riffle[list, x]` places `x` in each of the gaps between consecutive elements,
and `Riffle[list, {x1, ..., xk}]` cycles through the separators left to right. A
list of `n` elements has `n - 1` gaps, so separators go only *between* elements —
never before the first or after the last — and the output has `2n - 1` slots.

A single element (or an empty list) has no gaps, so the list is returned
unchanged. The head of the first argument is preserved, so
`Riffle[f[a, b], x]` gives `f[a, x, b]`. Interleaving a separator and then
`StringJoin`-ing is the usual way to build a delimited string.
