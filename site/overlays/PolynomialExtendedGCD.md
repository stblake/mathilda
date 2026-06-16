### Worked examples

```mathematica
In[1]:= PolynomialExtendedGCD[x^2 - 1, x^3 - 1, x]
Out[1]= {-1 + x, {-x, 1}}
```

```mathematica
In[1]:= PolynomialExtendedGCD[x^4 + x^3 + x^2 + x + 1, x^2 + 1, x]
Out[1]= {1, {1, -x - x^2}}
```

```mathematica
In[1]:= PolynomialExtendedGCD[x^7 - 1, x^5 - 1, x]
Out[1]= {-1 + x, {-x - x^3, 1 + x^3 + x^5}}
```

### Notes

`PolynomialExtendedGCD[a, b, x]` returns `{g, {s, t}}` where `g` is the
polynomial GCD and `s`, `t` are the Bezout cofactors satisfying
`s a + t b == g`. The first example certifies `(-x)(x^2-1) + 1*(x^3-1) = x - 1`,
the GCD. When the inputs are coprime the GCD is `1` and the cofactors give a
constructive proof of coprimality — exactly the data needed to invert one
polynomial modulo another (e.g. building inverses in `Q[x]/(b)`). For the two
cyclotomic-flavoured inputs `x^7-1` and `x^5-1`, whose only common root is the
trivial seventh-and-fifth root `x = 1`, the GCD is the linear factor `x - 1`.
