---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), Ch. 3 (normal forms and the distributive expansion of polynomials)."
---
### Worked examples

```mathematica
In[1]:= Expand[(x + 1)^3]
Out[1]= 1 + 3 x + 3 x^2 + x^3
```

```mathematica
In[1]:= Expand[(x + 1)^4]
Out[1]= 1 + 4 x + 6 x^2 + 4 x^3 + x^4
```

```mathematica
In[1]:= Expand[(a + b)(c + d)]
Out[1]= a c + b c + a d + b d
```

```mathematica
In[1]:= Expand[(x + 2)^2 (x - 1)]
Out[1]= -4 + 3 x^2 + x^3
```

### Notes

`Expand` applies the distributive law to products and integer powers,
producing a flat sum of monomials in canonical order (ascending total
degree in the leading variable). Like terms are combined automatically, so
`(x + 2)^2 (x - 1)` collapses the `x^1` coefficient to zero and it drops out
of the result. `Expand` only multiplies out — it does not factor or cancel —
and a second argument `Expand[expr, patt]` leaves alone any parts free of the
pattern `patt`.
