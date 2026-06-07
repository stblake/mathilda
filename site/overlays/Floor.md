### Worked examples

```mathematica
In[1]:= Floor[7/2]
Out[1]= 3

In[2]:= Floor[-2.3]
Out[2]= -3

In[3]:= Floor[17, 5]
Out[3]= 15
```

### Notes

`Floor[x]` rounds toward `-Infinity`; the two-argument `Floor[x, a]` gives the greatest multiple of `a` not exceeding `x`. Exact inputs return exact integers.
