### Worked examples

```mathematica
In[1]:= 2 == 2
Out[1]= True

In[2]:= 1 == 1.
Out[2]= True

In[3]:= a == b
Out[3]= a == b
```

### Notes

Unlike `SameQ`, `Equal` (`==`) tests mathematical equality, so `1 == 1.` is `True`. When equality cannot be decided, the call stays unevaluated as a symbolic equation.
