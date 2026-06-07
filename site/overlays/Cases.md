---
status: Stable
references:
---
### Worked examples

```mathematica
In[1]:= Cases[{1, a, 2, b, 3}, _Integer]
Out[1]= {1, 2, 3}
```

```mathematica
In[1]:= Cases[{1, 2, 3, 4}, x_ /; x > 2]
Out[1]= {3, 4}
```

```mathematica
In[1]:= Cases[{f[1], g[2], f[3]}, f[x_] -> x]
Out[1]= {1, 3}
```

```mathematica
In[1]:= Cases[{{1, 2}, {3, 4}}, _Integer, 2]
Out[1]= {1, 2, 3, 4}
```

### Notes

`Cases[list, pattern]` returns the elements that match `pattern`. With a
`pattern -> rhs` rule it instead returns the transformed values (`f[x_] -> x`
extracts the argument of every `f`). Conditional patterns (`x_ /; x > 2`) filter
on a predicate. A trailing level specification (here `2`) descends into nested
lists and collects matches from those deeper levels — without it `Cases` only
inspects level 1.
