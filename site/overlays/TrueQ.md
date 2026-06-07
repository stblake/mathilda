### Worked examples

```mathematica
In[1]:= TrueQ[True]
Out[1]= True

In[2]:= TrueQ[x]
Out[2]= False

In[3]:= TrueQ[1 < 2]
Out[3]= True
```

### Notes

`TrueQ` yields `True` only when its argument evaluates to the symbol `True`; anything else (including unevaluated symbolic expressions) yields `False`, making it useful as a safe guard in conditions.
