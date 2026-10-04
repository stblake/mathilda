### Worked examples

```mathematica
In[1]:= Product[k, {k, 1, n}]  (* a symbolic finite product is a factorial *)
```

```mathematica
In[1]:= Product[k^2, {k, 1, 5}]  (* a finite numeric range is multiplied out *)
```

```mathematica
In[1]:= Product[(k + 1)/k, {k, 1, n}]  (* a telescoping rational product *)
```

```mathematica
In[1]:= Product[2^k, {k, 1, n}]  (* a polynomial-exponential (geometric) product *)
```

```mathematica
In[1]:= Product[1 - 1/k^2, {k, 2, Infinity}]  (* a convergent infinite product *)
```

### Notes

`Product` is the multiplicative analogue of `Sum`. It is `HoldAll`, so the index
is localised and the iterator bounds are not evaluated against an outer binding.
A finite numeric range (or an explicit list of values) is multiplied out
directly, with an empty product giving `1`; a symbolic, indefinite, or convergent
infinite product is handed to a closed-form method cascade.

The cascade tries, cheapest-first, the telescoping (Gamma-free rational),
rational (Pochhammer / Gamma), geometric (`base^k`) and q-product families, plus
several infinite-product specialists. The method can be pinned with
`Method -> "Telescoping" | "Rational" | "Geometric" | "QProduct"`, and
convergence testing for infinite products can be disabled with
`VerifyConvergence -> False`. Multiple iterators `Product[f, s1, s2]` form nested
products, so an inner bound may depend on an outer index. When no method applies
the `Product[...]` is returned unevaluated.
