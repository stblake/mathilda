### Worked examples

```mathematica
In[1]:= f = Interpolation[{1, 4, 9, 16}]
Out[1]= InterpolatingFunction[{{1, 4}}, <>]

In[2]:= f[2.5]
Out[2]= 6.25
```

### Notes

`InterpolatingFunction[domain, table]` is the approximate-function object
returned by `Interpolation`; in standard output only the domain
(`{{xmin, xmax}, ...}`) is shown and the data table is elided as `<>`. Apply it
to a coordinate, `f[x]`, to get the tensor-product piecewise-polynomial value;
arguments outside the domain are extrapolated with a warning.
