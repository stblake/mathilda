### Worked examples

```mathematica
In[1]:= ReplaceAt[{a, b, c, d}, x_ -> X, 2]
Out[1]= {a, X, c, d}
```

```mathematica
In[1]:= ReplaceAt[{a, b, c, d}, x_ -> X, -1]
Out[1]= {a, b, c, X}
```

```mathematica
In[1]:= ReplaceAt[{{a, b}, {c, d}}, x_ -> X, {2, 1}]
Out[1]= {{a, b}, {X, d}}
```

```mathematica
In[1]:= ReplaceAt[1 + x + x^2 + x^3, e_ :> D[e, x], {2}]
Out[1]= 2 + x^2 + x^3
```

```mathematica
In[1]:= ReplaceAt[{1, 2, 3, 4, 5}, n_ :> n^2, {2 ;; 4}]
Out[1]= {1, 4, 9, 16, 5}
```

### Notes

`ReplaceAt` applies rules only at explicitly named positions, leaving the rest
of the expression untouched. Positions may be a single index, a part list
`{i, j, ...}`, a list of positions, `All`, or a `Span`. Negative indices count
from the end and `0` targets the head. Because the rule sees the *part* as a
whole expression, position-targeted rewrites can do real work: applying
`e_ :> D[e, x]` at part `{2}` of `1 + x + x^2 + x^3` differentiates just that
one summand (here the `x` term, whose derivative `1` merges into the leading
constant), while `n_ :> n^2` over the span `2 ;; 4` squares a contiguous slice.
