### Worked examples

```mathematica
In[1]:= Ceiling[7/2]
Out[1]= 4

In[2]:= Ceiling[-2.7]
Out[2]= -2

In[3]:= Ceiling[17, 5]
Out[3]= 20
```

### Notes

`Ceiling[x]` rounds toward `+Infinity`; the two-argument `Ceiling[x, a]` gives the smallest multiple of `a` at least `x`. Exact inputs return exact integers.
