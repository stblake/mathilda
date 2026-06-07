### Worked examples

```mathematica
In[1]:= ArcCos[0]
Out[1]= 1/2 Pi

In[2]:= ArcCos[1/2]
Out[2]= 1/3 Pi

In[3]:= ArcCos[1]
Out[3]= 0

In[4]:= N[ArcCos[0.5]]
Out[4]= 1.0472
```

### Notes

`ArcCos[z]` gives the principal inverse cosine, in `[0, Pi]` for real `z` in `[-1, 1]`. `ArcCos` is Listable.
