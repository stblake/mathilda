### Worked examples

```mathematica
In[1]:= CoefficientRules[x^2 + 2 x y + y^2, {x, y}]  (* exponent vector -> coefficient *)
```

```mathematica
In[2]:= CoefficientRules[3 x^2 + 1, x]  (* a single variable *)
```

```mathematica
In[3]:= CoefficientRules[1 + x y + x^3, {x, y}, "DegreeLexicographic"]  (* order by total degree *)
```

### Notes

`CoefficientRules[poly, {x1, ..., xk}]` gives a sparse `{expvec -> coeff, ...}`
view of `poly`: each rule's left side is the length-`k` integer exponent vector
of a monomial and its right side the coefficient. The variables default to
`Variables[poly]`, and an optional third argument sets the monomial order — the
six named orders (`"Lexicographic"` default, `"DegreeLexicographic"`,
`"DegreeReverseLexicographic"`, and `"Negative"` variants) or an explicit weight
matrix — with `Modulus -> m` reducing coefficients modulo `m`. It works whether
or not `poly` is expanded, and `FromCoefficientRules` is the exact inverse.
