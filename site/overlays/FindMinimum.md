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

```mathematica
In[1]:= FindMinimum[(1 - x)^2 + 100 (y - x^2)^2, {{x, -1.2}, {y, 1}}, Method -> "LBFGSB"]
Out[1]= {4.71192e-22, {x -> 1.0, y -> 1.0}}
```

```mathematica
In[1]:= FindMinimum[(x - 2)^2 + (y - 3)^2, {{x, 0, 0, 1}, {y, 0, 0, 1}}, Method -> "LBFGSB"]
Out[1]= {5.0, {x -> 1.0, y -> 1.0}}
```

### Notes

`FindMinimum[f, {x, x0}]` performs a local search from the start `x0`,
returning `{fmin, {x -> xmin, ...}}`. The third example is the notorious
Rosenbrock banana valley: BFGS quasi-Newton drives the iterate into the
curved trough and locates the global minimum `(1, 1)` to machine precision.
The Gamma example finds the minimum of the Gamma function on the positive
axis (a root of the digamma function) at `x ≈ 1.4616`.

`Method -> "LBFGSB"` selects **limited-memory BFGS with bound constraints**
(aliases `"LBFGS"`, `"LimitedMemoryBFGS"`; a Mathilda extension). It keeps only
the last 10 correction pairs — `O(m·n)` per step rather than the full-memory
QuasiNewton's `O(n²)` — so it scales to large variable counts, and it handles
box bounds by an active-set projection. The two `LBFGSB` examples above show an
unconstrained solve (the Rosenbrock valley, reached to machine precision) and a
box-constrained solve whose unconstrained minimum `(2, 3)` lies outside the unit
box, so the constrained optimum is pinned to the corner `(1, 1)` with value `5`.
