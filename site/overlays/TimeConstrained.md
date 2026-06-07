### Worked examples

```mathematica
In[1]:= TimeConstrained[1 + 1, 5]
Out[1]= 2

In[2]:= TimeConstrained[2^10, 5]
Out[2]= 1024
```

### Notes

`TimeConstrained[expr, t]` returns the result of `expr` if it finishes within `t` seconds; otherwise evaluation is aborted (returning `$Aborted` by default).
