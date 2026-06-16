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

```mathematica
In[1]:= Default[q] = 0
Out[1]= 0

In[2]:= q[a_, b_.] := {a, b}
Out[2]= Null

In[3]:= q[7]
Out[3]= {7, 0}
```

### Notes

`Default[f] = v` sets the value used for a missing optional argument matched by `_.` (`Optional[Blank[]]`). When the optional argument is supplied it is used directly, otherwise the stored default is substituted before the rule fires. This is the mechanism that lets `Plus`- and `Times`-style patterns absorb absent terms; for an ad hoc head you register the default yourself, as the two examples (additive default `0`, multiplicative-style default `1`) show. To attach a per-position default use the `Default[f, i]` form.
