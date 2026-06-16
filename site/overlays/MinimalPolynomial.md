### Worked examples

```mathematica
In[1]:= MinimalPolynomial[Sqrt[2], x]
Out[1]= -2 + x^2
```

```mathematica
In[1]:= MinimalPolynomial[(1 + Sqrt[5])/2, x]
Out[1]= -1 - x + x^2
```

```mathematica
In[1]:= MinimalPolynomial[Sqrt[2] + Sqrt[3], x]
Out[1]= 1 - 10 x^2 + x^4
```

```mathematica
In[1]:= MinimalPolynomial[Cos[2 Pi/5], x]
Out[1]= -1 + 2 x + 4 x^2
```

```mathematica
In[1]:= MinimalPolynomial[Sqrt[2 + Sqrt[2]], x]
Out[1]= 2 - 4 x^2 + x^4
```

### Notes

`MinimalPolynomial[s, x]` returns the lowest-degree integer polynomial in `x`,
with positive leading coefficient and content 1, that has the algebraic number
`s` as a root. The first two examples recover the defining polynomials of `Sqrt[2]`
and the golden ratio (`x^2 - x - 1`). The real power shows up with compound
algebraic numbers: `Sqrt[2] + Sqrt[3]` is degree 4 over the rationals, and
`MinimalPolynomial` finds its quartic `x^4 - 10 x^2 + 1` by eliminating the
radicals with resultants — not by numerical root-finding. It also handles
algebraic constants beyond plain radicals: `Cos[2 Pi/5]` (a root of a cyclotomic
relation) yields `4 x^2 + 2 x - 1`, and the nested radical `Sqrt[2 + Sqrt[2]]`
yields the quartic `x^4 - 4 x^2 + 2`. The input may be built from rationals,
radicals, the imaginary unit, roots of unity, and `Root[]` objects.
