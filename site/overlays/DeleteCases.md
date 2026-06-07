### Worked examples

```mathematica
In[1]:= DeleteCases[{1,2,3,2},2]
Out[1]= {1, 3}

In[2]:= DeleteCases[{1,a,2,b,3},_Integer]
Out[2]= {a, b}

In[3]:= DeleteCases[{1,2,3,4,5},_?OddQ]
Out[3]= {2, 4}
```

### Notes

Removes every element matching the pattern; combine with `_Integer` or `_?OddQ` to filter by head or predicate. The default level is `{1}` (top-level elements).
