### Worked examples

```mathematica
In[1]:= SetAttributes[f, Flat]
Out[1]= Null

In[2]:= f[f[a, b], c]
Out[2]= f[a, b, c]

In[3]:= f[a, f[b, f[c, d]]]
Out[3]= f[a, b, c, d]
```

### Notes

`Flat` marks a head as associative, so nested calls with the same head are flattened into a single argument sequence. `Plus` and `Times` carry this attribute by default; the flattening is also accounted for during pattern matching.
