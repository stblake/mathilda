### Worked examples

```mathematica
In[1]:= PolynomialMod[7 x^2 + 5 x + 3, 3]
Out[1]= 2 x + x^2
```

```mathematica
In[1]:= PolynomialMod[x^4, x^2 + 1]
Out[1]= 1
```

```mathematica
In[1]:= PolynomialMod[x^5 + x^3 + x, x^2 - 2]
Out[1]= 7 x
```

### Notes

`PolynomialMod` has two modes. With an integer modulus each coefficient is
reduced to its canonical residue in `{0, ..., m-1}`, so `7 x^2 + 5 x + 3`
modulo `3` becomes `x^2 + 2 x` (the constant `3` vanishes). With a polynomial
modulus the input is reduced as a polynomial: working in `Q[x]/(x^2+1)` the
relation `x^2 ≡ -1` gives `x^4 ≡ 1`, and in `Q[x]/(x^2-2)` the relation
`x^2 ≡ 2` collapses `x^5 + x^3 + x = 4x + 2x + x = 7 x`. Unlike
`PolynomialRemainder`, the leading coefficient of the polynomial modulus is
not normalised to one.
