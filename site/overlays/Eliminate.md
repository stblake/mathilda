### Worked examples

```mathematica
In[1]:= Eliminate[{x + y == 2, x - y == 0}, y]
Out[1]= x == 1

In[2]:= Eliminate[{a == b + c, d == a - c}, c]
Out[2]= d == b
```

### Notes

`Eliminate[eqns, vars]` removes `vars` from a system of polynomial equations
over the rationals (via a lexicographic Gröbner basis with an elimination
block), returning the relations among the remaining variables. It yields
`True` if the elimination ideal is empty and `False` if the system is
inconsistent; non-polynomial systems return unevaluated with `Eliminate::nlin`.
