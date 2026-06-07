### Worked examples

```mathematica
In[1]:= FactorTermsList[6 x^2 + 4 x]
Out[1]= {2, 2 x + 3 x^2}

In[2]:= FactorTermsList[2 x^2 + 4 x + 2]
Out[2]= {2, 1 + 2 x + x^2}
```

### Notes

`FactorTermsList[poly]` returns `{overall numerical factor, polynomial with that factor removed}`. The numerical content (here 2) is pulled out, leaving the primitive part; it does not factor the remaining polynomial further.
