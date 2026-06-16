### Worked examples

```mathematica
In[1]:= 1 < 2 < 3
Out[1]= True

In[2]:= 1 < 5 < 3
Out[2]= False

In[3]:= 2 <= 2 < 5
Out[3]= True
```

```mathematica
In[1]:= 1 < x < 3 < 5
Out[1]= 1 < x < 3
```

```mathematica
In[1]:= 2 < 3 < x < 1
Out[1]= 3 < x < 1
```

### Notes

A chained comparison such as `a < b <= c` parses to the canonical `Inequality` form and holds only when every adjacent pair holds. Mixed operators (e.g. `<=` and `<`) may be combined in a single chain. When some pairs are undecidable (because a value is symbolic), the decidable-and-true pairs are dropped and only the residual chain is returned — `1 < x < 3 < 5` collapses the true `3 < 5` and keeps `1 < x < 3`.
