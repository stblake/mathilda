### Worked examples

```mathematica
In[1]:= f[n_] := n^2
Out[1]= Null

In[2]:= f[4]
Out[2]= 16

In[3]:= f[a + b]
Out[3]= (a + b)^2
```

### Notes

`:=` (`SetDelayed`) holds the right-hand side and re-evaluates it each time the rule fires, with the pattern bindings substituted in. The definition itself returns `Null`.
