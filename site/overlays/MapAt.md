### Worked examples

```mathematica
In[1]:= MapAt[f, {a, b, c, d}, 2]
Out[1]= {a, f[b], c, d}
```

```mathematica
In[1]:= MapAt[f, {a, b, c, d}, -1]
Out[1]= {a, b, c, f[d]}
```

```mathematica
In[1]:= MapAt[f, {{a, b}, {c, d}}, {2, 1}]
Out[1]= {{a, b}, {f[c], d}}
```

```mathematica
In[1]:= MapAt[f, {a, b, c, d}, {{1}, {3}}]
Out[1]= {f[a], b, f[c], d}
```

```mathematica
In[1]:= MapAt[Framed, {1, 2, 3, 4, 5}, {{1}, {-1}}]
Out[1]= {Framed[1], 2, 3, 4, Framed[5]}
```

### Notes

`MapAt[f, expr, n]` applies `f` to the element at position `n`; negative `n` counts from the end and `0` targets the head. A position list like `{2, 1}` selects a part deep in a nested expression, while a list of positions `{{1}, {3}}` applies `f` at several places at once (last two examples). Positions may also use `All` or `Span`. Repeated positions apply `f` more than once at that part.
