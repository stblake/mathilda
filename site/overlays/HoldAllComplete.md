### Worked examples

```mathematica
In[1]:= SetAttributes[h, HoldAllComplete]
Out[1]= Null

In[2]:= h[1+1]
Out[2]= h[1 + 1]

In[3]:= Attributes[h]
Out[3]= {HoldAllComplete}
```

### Notes

`HoldAllComplete` is the strongest hold attribute: it not only suppresses argument evaluation but also blocks `Sequence` flattening, `Unevaluated` stripping, and `Evaluate`, so the arguments are passed through untouched.
