### Worked examples

```mathematica
In[1]:= FindMinimum[x^2 - 4 x + 7, {x, 0}]
Out[1]= {3.0, {x -> 2.0}}
```

```mathematica
In[1]:= FindMinimum[Cos[x] + x/5, {x, 0, 10}]
Out[1]= {-0.391749, {x -> 2.94023}}
```

```mathematica
In[1]:= FindMinimum[(1 - x)^2 + 100 (y - x^2)^2, {{x, -1}, {y, 1}}]
Out[1]= {3.46541e-23, {x -> 1.0, y -> 1.0}}
```

```mathematica
In[1]:= FindMinimum[Gamma[x], {x, 1.5}]
Out[1]= {0.885603, {x -> 1.46163}}
```

### Notes

`FindMinimum[f, {x, x0}]` performs a local search from the start `x0`,
returning `{fmin, {x -> xmin, ...}}`. The third example is the notorious
Rosenbrock banana valley: BFGS quasi-Newton drives the iterate into the
curved trough and locates the global minimum `(1, 1)` to machine precision.
The last example finds the minimum of the Gamma function on the positive
axis (a root of the digamma function) at `x ≈ 1.4616`.
