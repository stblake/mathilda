### Worked examples

```mathematica
In[1]:= f = Interpolation[{1, 4, 9, 16}]
Out[1]= InterpolatingFunction[{{1, 4}}, <>]

In[2]:= f[2.5]
Out[2]= 6.25
```

```mathematica
In[1]:= g = Interpolation[Table[{x, Sin[x]}, {x, 0., 6., 0.5}]]
Out[1]= InterpolatingFunction[{{0.0, 6.0}}, <>]

In[2]:= g[1.5]
Out[2]= 0.997495

In[3]:= Sin[1.5]
Out[3]= 0.997495
```

```mathematica
In[1]:= d = Interpolation[{1, 4, 9, 16, 25}]; dd = d'; dd[2.5]
Out[1]= 5.0
```

### Notes

`InterpolatingFunction[domain, table]` is the approximate-function object
returned by `Interpolation`; in standard output only the domain
(`{{xmin, xmax}, ...}`) is shown and the data table is elided as `<>`. Apply it
to a coordinate, `f[x]`, to get the tensor-product piecewise-polynomial value;
arguments outside the domain are extrapolated with a warning.
