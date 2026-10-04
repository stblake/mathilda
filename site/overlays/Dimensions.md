### Worked examples

```mathematica
In[1]:= Dimensions[{{1, 2}, {3, 4}}]  (* a full rectangular matrix *)
```

```mathematica
In[2]:= Dimensions[{{a, b, c}, {d, e}, {f}}]  (* ragged below level 1, so only the outer length is counted *)
```

```mathematica
In[3]:= Dimensions[{{{{a, b}}}}]  (* four nested levels *)
```

```mathematica
In[4]:= Dimensions[{{{{a, b}}}}, 2]  (* capped at the first two levels *)
```

```mathematica
In[5]:= Dimensions[x]  (* an atom has no parts *)
```

### Notes

`Dimensions[expr]` reports the shape only down to the level at which `expr` stops
being a full array — a level counts only when every sub-piece there shares the
same head and length, so the second example stops at `{3}` rather than inventing
a second dimension for the ragged rows. The scan uses a fixed per-level stack
(`DIMENSIONS_MAX_DEPTH = 64`) and takes its reference head from the top-level
expression, so a nested `List` of `List`s is measured as a tensor. An atom
returns the empty `List` `{}`, and `Dimensions[expr, n]` truncates the result to
the first `n` levels.
