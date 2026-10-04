### Worked examples

```mathematica
In[1]:= MonomialList[x^2 + 2 x y + y^2, {x, y}]  (* the individual monomials, with coefficients *)
```

```mathematica
In[2]:= MonomialList[(x + y)^2, {x, y}]  (* expands first, so an unexpanded form works *)
```

```mathematica
In[3]:= MonomialList[1 + x y + x^3, {x, y}, "DegreeLexicographic"]  (* highest total degree first *)
```

### Notes

`MonomialList[poly, {x1, ..., xk}]` gives the list of monomials of `poly`, each
with its coefficient, so their sum is `poly`. The variables default to
`Variables[poly]` (or `All`), and an optional order argument sorts the monomials:
the default `"Lexicographic"`, the degree orders
(`"DegreeLexicographic"`/`"DegreeReverseLexicographic"`), their `"Negative"`
variants, or an explicit weight matrix; `Modulus -> m` reduces coefficients. The
polynomial is expanded first, so an unexpanded form such as `(x + y)^2` is
accepted. `CoefficientRules` is the same decomposition rendered as `expvec ->
coeff` rules.
