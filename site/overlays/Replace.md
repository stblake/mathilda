### Worked examples

```mathematica
In[1]:= Replace[x^2, x^2 -> done]
Out[1]= done
```

```mathematica
In[1]:= Replace[x, {x -> 1, _ -> 0}]
Out[1]= 1

In[2]:= Replace[w, {x -> 1, _ -> 0}]
Out[2]= 0
```

```mathematica
In[1]:= Replace[a + b + c + d, x_ + y_ -> {x, y}]
Out[1]= {a, b + c + d}
```

```mathematica
In[1]:= Replace[{{1, 2}, {3, 4}}, x_Integer :> x^2, {2}]
Out[1]= {{1, 4}, {9, 16}}
```

```mathematica
In[1]:= Replace[{1, {2, {3, {4}}}}, x_Integer :> x^2, {-1}]
Out[1]= {1, {4, {9, {16}}}}
```

### Notes

`Replace[expr, rules]` matches `expr` only at the **top level** (unlike
`ReplaceAll`, which descends into every subexpression). The first matching
rule in the list wins; a fall-through `_ -> 0` rule acts as a default. A
`levelspec` third argument restricts where matching happens: `{2}` targets
parts two levels deep, and the negative spec `{-1}` targets the leaves, so
`x_Integer :> x^2` squares every atom in a nested list. Because matching
respects flat/orderless heads, `x_ + y_` against a four-term sum binds `x` to
the first term and `y` to the rest.
