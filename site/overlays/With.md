### Worked examples

```mathematica
In[1]:= With[{a = 2, b = 3}, a^2 + b^2]
Out[1]= 13
```

The bound values are substituted *before* the body evaluates, so injecting an exact algebraic constant lets the symbolic engine verify a defining identity — here the golden ratio satisfying `phi^2 - phi - 1 == 0`:

```mathematica
In[1]:= With[{phi = (1 + Sqrt[5])/2}, Simplify[phi^2 - phi - 1]]
Out[1]= 0
```

`With` localizes the value cleanly for numeric work too, fixing a parameter and driving a high-precision computation:

```mathematica
In[1]:= With[{n = 20}, N[Sum[1/k^2, {k, 1, n}], 30]]
Out[1]= 1.596163243913023316640878872058
```

Bindings flow through complex arithmetic transparently:

```mathematica
In[1]:= With[{x = 1 + I}, x^2]
Out[1]= 2*I
```

### Notes

`With[{x = x0, ...}, expr]` replaces each `x` by its value `x0` throughout `expr` and then evaluates the result. Unlike `Module`, the substitution is literal and immediate, making `With` ideal for inlining constants and parameters.
