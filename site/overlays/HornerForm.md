### Worked examples

```mathematica
In[1]:= HornerForm[1 + x + x^2 + x^3]
Out[1]= 1 + x (1 + x (1 + x))

In[2]:= HornerForm[a x^3 + b x^2 + c x + d, x]
Out[2]= d + x (c + x (b + a x))
```

Sparse and mixed-sign polynomials nest just as cleanly, skipping absent
degrees:

```mathematica
In[1]:= HornerForm[3 x^4 - 2 x^3 + x - 7]
Out[1]= -7 + x (1 + x^2 (-2 + 3 x))
```

For multivariate input, choose the recursion variable; the other variables
ride along inside the coefficients:

```mathematica
In[1]:= HornerForm[1 + 2 x + 3 x^2 y + 4 x y^2, x]
Out[1]= 1 + x (2 + 3 x y + 4 y^2)
```

A rational function can be nested numerator and denominator independently:

```mathematica
In[1]:= HornerForm[(x^2 + 1)/(x^3 - x + 2), x, x]
Out[1]= (1 + x^2)/(2 + x (-1 + x^2))
```

### Notes

The nested (Horner) form evaluates a degree-n polynomial in n multiplications and n additions instead of the naive 2n. Pass an explicit variable as the second argument when the coefficients are themselves symbolic.
