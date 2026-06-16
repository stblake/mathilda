### Worked examples

```mathematica
In[1]:= f = Interpolation[{1, 4, 9, 16}]
Out[1]= InterpolatingFunction[{{1, 4}}, <>]

In[2]:= f[2]
Out[2]= 4

In[3]:= f[2.5]
Out[3]= 6.25
```

```mathematica
In[1]:= g = Interpolation[Table[{x, Sin[x]}, {x, 0., 6., 0.5}]]
Out[1]= InterpolatingFunction[{{0.0, 6.0}}, <>]

In[2]:= g[1.5]
Out[2]= 0.997495
```

```mathematica
In[1]:= p = Interpolation[{{0, 0}, {1, 1}, {2, 4}, {3, 9}}, InterpolationOrder -> 2]
Out[1]= InterpolatingFunction[{{0, 3}}, <>]

In[2]:= p[1.5]
Out[2]= 2.25
```

### Notes

`Interpolation[data]` builds an `InterpolatingFunction` (default piecewise
cubic) over the given samples; here `{1, 4, 9, 16}` are the values at
`x = 1, 2, 3, 4`, so evaluating recovers `x^2`. The returned object prints with
only its domain shown and is applied like an ordinary function, `f[x]`.
