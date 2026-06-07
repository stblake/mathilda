### Worked examples

```mathematica
In[1]:= f[x_, y_ : 1] := x + y

In[2]:= f[3]
Out[2]= 4

In[3]:= f[3, 10]
Out[3]= 13

In[4]:= g[x_, y_ : 0, z_ : 0] := {x, y, z}

In[5]:= g[1, 2]
Out[5]= {1, 2, 0}
```

### Notes

`Optional[p, def]` (surface syntax `p : def`) lets a pattern argument be omitted, supplying `def` in its place — the standard way to give function definitions default-valued parameters (Out[2], Out[5]).
