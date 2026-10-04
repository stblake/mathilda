### Worked examples

```mathematica
In[1]:= PolynomialSqrt[x^2 + 2 x + 1]  (* (1 + x)^2 *)
```

```mathematica
In[1]:= PolynomialSqrt[(x^2 - 1)^2]  (* returned in factored form *)
```

```mathematica
In[1]:= PolynomialSqrt[x^2 + 1]  (* not a perfect square *)
```

### Notes

`PolynomialSqrt[p]` returns a polynomial `s` with `s^2 == p` when `p` is a perfect
square — every non-constant irreducible factor must have even multiplicity, and the
numeric content is carried through `Sqrt` — and `$Failed` otherwise. `PolynomialSqrt[p,
x]` treats `p` as a polynomial in `x`, so any factor free of `x` counts as constant
content.

The result is given in **factored** form: `PolynomialSqrt[(x^2 - 1)^2]` is
`(-1 + x)(1 + x)`, not the expanded `x^2 - 1`. Every success carries an exact
certificate — `Expand[s^2 - p]` must be zero before `s` is returned — so a near-miss is
reported as `$Failed` rather than an approximate root. `PolynomialSqrt` is `Protected`.
