### Worked examples

```mathematica
In[1]:= FindMaximum[Sin[x], {x, 1}]
Out[1]= {1.0, {x -> 1.5708}}
```

```mathematica
In[1]:= FindMaximum[x (10 - x), {x, 0}]
Out[1]= {25.0, {x -> 5.0}}
```

```mathematica
In[1]:= FindMaximum[Sin[x] Sin[2 y], {{x, 1}, {y, 1}}]
Out[1]= {1.0, {x -> 1.5708, y -> 0.785398}}
```

### Notes

`FindMaximum[f, {x, x0}]` returns `{fmax, {x -> xmax, ...}}`. Internally it
maximises by minimising `-f`, so the same Brent (1D) and BFGS quasi-Newton
(n-D) machinery as `FindMinimum` applies. The first example recovers the
peak of `Sin` at `x = π/2`; the multivariate case locates a saddle-free
maximum of the product `Sin[x] Sin[2 y]` at `(π/2, π/4)`.
