---
status: Stable
references:
  - "Cox, Little & O'Shea, \"Ideals, Varieties, and Algorithms\" (Springer) — Buchberger's algorithm, monomial orders, and elimination."
  - "von zur Gathen & Gerhard, \"Modern Computer Algebra\" (3rd ed.), Ch. 21 (Gröbner bases)."
---
### Worked examples

```mathematica
In[1]:= GroebnerBasis[{x^2 + y^2 - 1, x - y}, {x, y}]
Out[1]= {-1 + 2 y^2, x - y}
```

```mathematica
In[1]:= GroebnerBasis[{x y - 1, x - y}, {x, y}]
Out[1]= {-1 + y^2, x - y}
```

```mathematica
In[1]:= GroebnerBasis[{x^2 - y, x^3 - z}, {x, y, z}]
Out[1]= {y^3 - z^2, -y^2 + x z, x y - z, x^2 - y}
```

```mathematica
In[1]:= GroebnerBasis[{x + y + z, x y + y z + z x, x y z - 1}, {x, y, z}]
Out[1]= {-1 + z^3, y^2 + y z + z^2, x + y + z}
```

### Notes

`GroebnerBasis[polys, vars]` computes a Gröbner basis of the ideal under the
default Lexicographic monomial order via Buchberger's algorithm. The lex order
triangularizes the system: in the circle-meets-line example the basis isolates
`-1 + 2 y^2` (a univariate consequence in the last variable) plus the linear
relation. The symmetric-functions system reduces to `z^3 = 1`, exposing the
roots of unity, while the implicitization example returns the full reduced
basis describing the twisted cubic. Set
`MonomialOrder -> DegreeReverseLexicographic` for systems with three or more
variables, where lex order can be exponentially slow; a third argument selects
elimination variables, and Equal equations are accepted in place of bare
polynomials.
