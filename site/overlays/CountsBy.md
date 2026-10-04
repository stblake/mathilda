### Worked examples

```mathematica
In[1]:= CountsBy[{1, 2, 3, 4, 5}, EvenQ]
```

```mathematica
In[1]:= CountsBy[{-2, -1, 0, 1, 2}, Sign]
```

```mathematica
In[1]:= CountsBy[Range[10], Mod[#, 3] &]
```

### Notes

`CountsBy[list, f]` returns `<|f[x] -> count|>` — a histogram keyed by the value
of `f` applied to each element, in first-appearance order of the `f`-values. It
is `Counts` after bucketing by `f`: `CountsBy[list, Sign]` tallies how many
elements are negative, zero, and positive. Over an association the values are
tallied by `f`.
