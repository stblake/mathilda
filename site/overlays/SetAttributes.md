### Worked examples

```mathematica
In[1]:= SetAttributes[g, Orderless]
Out[1]= Null

In[2]:= Attributes[g]
Out[2]= {Orderless}

In[3]:= g[3, 1, 2]
Out[3]= g[1, 2, 3]
```

### Notes

`SetAttributes[s, attr]` adds an attribute to `s`. Here `Orderless` makes the evaluator sort `g`'s arguments into canonical order on every call.
