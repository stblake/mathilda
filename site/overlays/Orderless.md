### Worked examples

```mathematica
In[1]:= SetAttributes[g, Orderless]
Out[1]= Null

In[2]:= g[3, 1, 2]
Out[2]= g[1, 2, 3]

In[3]:= g[c, a, b]
Out[3]= g[a, b, c]
```

### Notes

`Orderless` marks a head as commutative, so its arguments are automatically sorted into canonical order. `Plus` and `Times` carry this attribute by default; canonical ordering is also used when matching patterns against `Orderless` heads.
