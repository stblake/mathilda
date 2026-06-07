### Worked examples

```mathematica
In[1]:= ArcSin[1/2]
Out[1]= 1/6 Pi

In[2]:= ArcSin[1]
Out[2]= 1/2 Pi

In[3]:= N[ArcSin[0.5]]
Out[3]= 0.523599

In[4]:= ArcSin[Sin[x]]
Out[4]= ArcSin[Sin[x]]
```

### Notes

`ArcSin[z]` gives the principal inverse sine, in `[-Pi/2, Pi/2]` for real `z` in `[-1, 1]`. The inverse-of-forward composition `ArcSin[Sin[x]]` is deliberately *not* folded to `x`, since that holds only on the principal branch. `ArcSin` is Listable.
