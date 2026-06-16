### Worked examples

```mathematica
In[1]:= Take[{a, b, c, d, e}, 3]
Out[1]= {a, b, c}
```

```mathematica
In[1]:= Take[{a, b, c, d, e}, -2]
Out[1]= {d, e}
```

```mathematica
In[1]:= Take[Range[10], {2, 8, 2}]
Out[1]= {2, 4, 6, 8}
```

```mathematica
In[1]:= Take[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, 2, 2]
Out[1]= {{1, 2}, {4, 5}}
```

```mathematica
In[1]:= Take[Table[Fibonacci[n], {n, 1, 15}], {3, 15, 3}]
Out[1]= {2, 8, 34, 144, 610}
```

### Notes

`Take[list, n]` takes the first `n` elements, `Take[list, -n]` the last `n`, and
`Take[list, {m, n}]` (optionally `{m, n, s}` with a step) an inclusive index
range. Indices are 1-based and negative indices count from the end; `UpTo[n]`,
`All`, and `None` are also accepted. Multiple specifications act level by level,
so `Take[mat, 2, 2]` extracts the top-left 2x2 sub-block of a matrix. `Take`
operates on any expression, not only `List`; out-of-range requests are left
unevaluated.
