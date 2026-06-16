### Worked examples

```mathematica
In[1]:= Eliminate[{x + y == 2, x - y == 0}, y]
Out[1]= x == 1

In[2]:= Eliminate[{a == b + c, d == a - c}, c]
Out[2]= d == b
```

```mathematica
In[1]:= Eliminate[{a == x + y, b == x y}, {x, y}]
Out[1]= True

In[2]:= Eliminate[{p == x + 1/x, q == x^2 + 1/x^2}, x]
Out[2]= 2 + q == p^2
```

```mathematica
In[1]:= Eliminate[{x == a Cos[t], y == a Sin[t]}, t]
Out[1]= x^2 + y^2 == a^2
```

```mathematica
In[1]:= Eliminate[{u == Exp[x], v == Exp[2 x]}, x]
Out[1]= v == u^2
```

### Notes

`Eliminate[eqns, vars]` removes `vars` from a system of polynomial equations
over the rationals (via a lexicographic Gröbner basis with an elimination
block), returning the relations among the remaining variables. It yields
`True` if the elimination ideal is empty and `False` if the system is
inconsistent; non-polynomial systems return unevaluated with `Eliminate::nlin`.
