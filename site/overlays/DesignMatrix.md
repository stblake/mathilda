### Worked examples

```mathematica
In[1]:= DesignMatrix[{{0, 1}, {1, 0}, {3, 2}, {5, 4}}, {1, x}, x]
Out[1]= {{1, 0}, {1, 1}, {1, 3}, {1, 5}}
```

```mathematica
In[1]:= DesignMatrix[{{1, 1}, {2, 8}, {3, 27}}, {1, x, x^2, x^3}, x]
Out[1]= {{1, 1, 1, 1}, {1, 2, 4, 8}, {1, 3, 9, 27}}
```

```mathematica
In[1]:= DesignMatrix[{{1, 1, 5}, {2, 4, 6}, {3, 9, 2}}, {1, x, y, x*y}, {x, y}]
Out[1]= {{1, 1, 1, 1}, {1, 2, 4, 8}, {1, 3, 9, 27}}
```

```mathematica
In[1]:= DesignMatrix[{{1, 2}, {2, 5}}, {1, Sin[x]}, x, WorkingPrecision -> 40]
Out[1]= {{1.0, 0.84147098480789650665250232163029899962254}, {1.0, 0.90929742682568169539601986591174484270222}}
```

### Notes

`DesignMatrix[data, {f1, ..., fn}, vars]` builds the matrix whose rows are the
basis functions `f_i` evaluated at each data coordinate — exactly the matrix
`Fit` assembles internally before solving the normal equations. The data shapes
match `Fit`: with a single variable, each row is `{x_k, y_k}` (or just `{x_k}`)
and only the leading coordinate(s) are substituted, so the response column is
ignored when forming the design entries. A polynomial basis `{1, x, x^2, x^3}`
therefore produces a Vandermonde matrix. For several predictors the basis may mix
the variables freely (`{1, x, y, x*y}`). Entries are kept exact unless
`WorkingPrecision -> MachinePrecision` or a digit count is supplied, in which
case each entry is converted to an approximate number — useful when the basis
contains transcendental functions such as `Sin`.
