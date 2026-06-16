### Worked examples

```mathematica
In[1]:= PolynomialQ[x^3 + 2 x + 1, x]
Out[1]= True
```

```mathematica
In[1]:= PolynomialQ[Sin[x] + x, x]
Out[1]= False
```

```mathematica
In[1]:= PolynomialQ[x^2 y + x y^2 + 1, {x, y}]
Out[1]= True
```

```mathematica
In[1]:= PolynomialQ[x^2 + y/x, x]
Out[1]= False
```

### Notes

`PolynomialQ[expr, var]` tests whether `expr` is a polynomial in `var` — i.e.
expands to a sum of products of non-negative integer powers of the variables
with variable-free coefficients. `Sin[x] + x` fails because of the
transcendental term, and `x^2 + y/x` fails because `y/x = y x^(-1)` carries a
negative power of `x`. The multivariate form `PolynomialQ[expr, {x, y}]`
requires polynomiality in all listed variables simultaneously. Note that the
test is syntactic up to expansion and does not cancel rational expressions: a
fraction such as `(x^2 - 1)/(x - 1)` is reported `False` even though it equals
the polynomial `x + 1`.
