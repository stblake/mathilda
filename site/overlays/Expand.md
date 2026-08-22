---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), Ch. 3 (normal forms and the distributive expansion of polynomials)."
---
### Worked examples

```mathematica
In[1]:= Expand[(x + 1)^3]   (* a binomial power: the coefficients are the binomial coefficients *)
Out[1]= 1 + 3 x + 3 x^2 + x^3
```

```mathematica
In[1]:= Expand[(a + b)(c + d)]   (* a product of sums distributes into every pairwise term *)
Out[1]= a c + b c + a d + b d
```

```mathematica
In[1]:= Expand[(1 + x + y)^3]   (* several variables: the multinomial, not just the binomial, case *)
Out[1]= 1 + 3 x + 3 x^2 + x^3 + 3 y + 6 x y + 3 x^2 y + 3 y^2 + 3 x y^2 + y^3
```

```mathematica
In[1]:= Expand[(x + 2)^2 (x - 1)]   (* numeric coefficients are folded, so terms collect and cancel *)
Out[1]= -4 + 3 x^2 + x^3
```

```mathematica
In[1]:= Expand[(1 + x)^10]   (* the work is polynomial in the exponent, not exponential *)
Out[1]= 1 + 10 x + 45 x^2 + 120 x^3 + 210 x^4 + 252 x^5 + 210 x^6 + 120 x^7 + 45 x^8 + 10 x^9 + x^10
```

```mathematica
In[1]:= Expand[(x + 1)^2/(y + 1)]   (* only the numerator expands: Expand leaves denominators alone *)
Out[1]= (1 + 2 x + x^2)/(1 + y)
```

### Notes

`Expand` applies the distributive law to products and integer powers,
producing a flat sum of monomials in canonical order (ascending total
degree in the leading variable). Like terms are combined automatically, so
`(x + 2)^2 (x - 1)` collapses the `x^1` coefficient to zero and it drops out
of the result. `Expand` only multiplies out — it does not factor or cancel —
and a second argument `Expand[expr, patt]` leaves alone any parts free of the
pattern `patt`.
