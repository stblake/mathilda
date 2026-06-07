---
status: Stable
references:
---
### Worked examples

```mathematica
In[1]:= Range[5]
Out[1]= {1, 2, 3, 4, 5}
```

```mathematica
In[1]:= Range[2, 10, 2]
Out[1]= {2, 4, 6, 8, 10}
```

```mathematica
In[1]:= Range[0, 1, 1/4]
Out[1]= {0, 1/4, 1/2, 3/4, 1}
```

### Notes

`Range[n]` produces the integers `1` through `n`. `Range[n, m]` runs from `n` to
`m` in unit steps, and `Range[n, m, d]` uses step `d`. The step may be an exact
rational, in which case the result stays exact rather than being converted to
floating point. `Range` is the most direct way to build the index list that
`Map`, `Select`, or `Fold` then consume.
