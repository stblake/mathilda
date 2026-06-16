### Worked examples

```mathematica
In[1]:= FactorTerms[6 x^2 + 4 x]
Out[1]= 2 (2 x + 3 x^2)
```

```mathematica
In[1]:= FactorTerms[2 x^2 + 4 x + 2]
Out[1]= 2 (1 + 2 x + x^2)
```

```mathematica
In[1]:= FactorTerms[3 x^2 y + 6 x y^2, x]
Out[1]= 3 y (x^2 + 2 x y)
```

### Notes

`FactorTerms[poly]` pulls out the overall numerical content (the integer GCD of the coefficients) without touching the polynomial structure — unlike `Factor`, it never splits `1 + 2 x + x^2` into `(1 + x)^2`. With a second argument `FactorTerms[poly, x]` extracts the content with respect to `x`, i.e. the factor that does not depend on `x`: here `3 y` is the part free of `x`, leaving `x^2 + 2 x y`.
