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

```mathematica
In[1]:= FindMaximum[10 - (x - 3)^2 - (y + 1)^2, {{x, 0}, {y, 0}}, Method -> "LBFGSB"]
Out[1]= {10.0, {x -> 3.0, y -> -1.0}}
```

### Notes

`FindMaximum[f, {x, x0}]` returns `{fmax, {x -> xmax, ...}}`. Internally it
maximises by minimising `-f`, so the same Brent (1D) and BFGS quasi-Newton
(n-D) machinery as `FindMinimum` applies. The first example recovers the
peak of `Sin` at `x = π/2`; the multivariate case locates a saddle-free
maximum of the product `Sin[x] Sin[2 y]` at `(π/2, π/4)`.

Every `FindMinimum` method is available, including `Method -> "LBFGSB"`
(limited-memory BFGS with bound constraints; aliases `"LBFGS"`,
`"LimitedMemoryBFGS"`), shown in the last example locating the peak of a
concave paraboloid at `(3, -1)`.
