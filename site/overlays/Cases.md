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

```mathematica
In[1]:= Cases[{1, 2, 3, 4, 5}, _?PrimeQ]
Out[1]= {2, 3, 5}
```

```mathematica
In[1]:= Cases[{f[1, 2], f[3], g[4], f[5, 6]}, f[a_, b_] -> a + b]
Out[1]= {3, 11}
```

```mathematica
In[1]:= Cases[{x^2, y^3, z, w^4}, p_^n_ -> {p, n}]
Out[1]= {{x, 2}, {y, 3}, {w, 4}}
```

### Notes

`Cases[list, pattern]` returns the elements that match `pattern`. With a
`pattern -> rhs` rule it instead returns the transformed values (`f[x_] -> x`
extracts the argument of every `f`). Conditional patterns (`x_ /; x > 2`) filter
on a predicate. A trailing level specification (here `2`) descends into nested
lists and collects matches from those deeper levels — without it `Cases` only
inspects level 1.

`PatternTest` (`_?PrimeQ`) filters by an arbitrary predicate — here keeping only
the primes. Rules can do real structural work: `f[a_, b_] -> a + b` matches only
the two-argument `f` heads and returns the sums of their arguments, skipping
`f[3]` (wrong arity) and `g[4]` (wrong head). The pattern `p_^n_ -> {p, n}`
deconstructs every power into a `{base, exponent}` pair while ignoring the
non-power element `z`, illustrating how `Cases` doubles as a structural query
and extraction tool.
