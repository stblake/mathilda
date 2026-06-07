### Worked examples

```mathematica
In[1]:= MapAll[f, {a, {b, c}}]
Out[1]= f[{f[a], f[{f[b], f[c]}]}]

In[2]:= f //@ {a, b}
Out[2]= f[{f[a], f[b]}]

In[3]:= MapAll[g, 1 + x]
Out[3]= g[g[1] + g[x]]

In[4]:= Map[f, {a, {b, c}}]
Out[4]= {f[a], f[{b, c}]}
```

### Notes

`MapAll` applies `f` to every subexpression including atomic leaves, equivalent to `Map[f, expr, {0, Infinity}]`; its operator form is `f //@ expr`. Unlike `Map`, which only touches the first level (compare In[1] vs In[4]), `MapAll` reaches all levels and wraps the whole expression too.
