### Worked examples

```mathematica
In[1]:= Discriminant[x^2 - 5 x + 6, x]
Out[1]= 1
```

```mathematica
In[1]:= Discriminant[a x^2 + b x + c, x]
Out[1]= b^2 - 4 a c
```

```mathematica
In[1]:= Discriminant[x^3 + p x + q, x]
Out[1]= -4 p^3 - 27 q^2
```

```mathematica
In[1]:= Discriminant[x^4 + 1, x]
Out[1]= 256
```

### Notes

`Discriminant[poly, var]` is, up to sign and leading-coefficient scaling,
`Resultant[poly, D[poly, var], var] / lc[poly, var]`, and it vanishes exactly when
`poly` has a repeated root in `var`. The first example has distinct roots `2` and
`3`, so its discriminant is a nonzero constant. The second recovers the familiar
`b^2 - 4 a c` from the quadratic formula. The third gives the classical
depressed-cubic discriminant `-4 p^3 - 27 q^2`, whose sign distinguishes three real
roots from one real and two complex roots. The fourth shows `x^4 + 1`, which has
four distinct complex roots and discriminant `256`.
