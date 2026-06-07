### Worked examples

```mathematica
In[1]:= Default[h] = 1
Out[1]= 1

In[2]:= h[x_, y_.] := {x, y}
Out[2]= Null

In[3]:= h[5]
Out[3]= {5, 1}

In[4]:= h[5, 9]
Out[4]= {5, 9}
```

### Notes

`Default[f] = v` sets the value used for a missing optional argument matched by `_.` (`Optional[Blank[]]`). When the optional argument is supplied it is used directly.
