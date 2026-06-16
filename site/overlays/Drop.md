### Worked examples

```mathematica
In[1]:= Drop[{a, b, c, d, e}, 2]
Out[1]= {c, d, e}
```

```mathematica
In[1]:= Drop[{a, b, c, d, e}, -2]
Out[1]= {a, b, c}
```

```mathematica
In[1]:= Drop[{a, b, c, d, e}, {2, 4}]
Out[1]= {a, e}
```

```mathematica
In[1]:= Drop[{a, b, c, d, e, f, g}, {2, 7, 2}]
Out[1]= {a, c, e, g}
```

```mathematica
In[1]:= Drop[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, {2}, {2}]
Out[1]= {{1, 3}, {7, 9}}
```

### Notes

`Drop` is the complement of `Take`. A plain count drops from the front
(`Drop[list, n]`) or, with a negative count, from the back. The `{m, n}` form
drops a contiguous block, `{m, n, s}` drops a strided slice, and `{m}` drops a
single element. Multiple level specifications drop along successive list
dimensions, so the `{2}, {2}` example deletes the second row and the second
column of a matrix in one call. Indices are 1-based and negative indices count
from the end; out-of-range requests leave the expression unevaluated.
