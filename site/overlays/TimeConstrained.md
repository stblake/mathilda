### Worked examples

```mathematica
In[1]:= TimeConstrained[1 + 1, 5]
Out[1]= 2

In[2]:= TimeConstrained[2^10, 5]
Out[2]= 1024
```

```mathematica
In[1]:= TimeConstrained[Integrate[x^2 Exp[x], x], 10]
Out[1]= 2 E^x + x^2 E^x - 2 x E^x
```

```mathematica
In[1]:= TimeConstrained[Solve[x^2 - 3 x + 2 == 0, x], 10]
Out[1]= {{x -> 1}, {x -> 2}}
```

### Notes

`TimeConstrained[expr, t]` returns the result of `expr` if it finishes within `t` seconds; otherwise evaluation is aborted (returning `$Aborted` by default).
