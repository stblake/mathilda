### Worked examples

```mathematica
In[1]:= CoefficientArrays[{x + 2 y - 3, 3 x - y + 1}, {x, y}]  (* {constant vector, coefficient matrix} *)
```

```mathematica
In[1]:= CoefficientArrays[1 + 2 x + 3 x^2, x]  (* grouped by degree: {const, {linear}, {{quadratic}}} *)
```

### Notes

`CoefficientArrays[polys, vars]` returns the coefficient arrays of a polynomial system,
grouped by total degree: element `d + 1` holds every degree-`d` coefficient, with a
degree-`d` coefficient placed at the index tuple naming its variables in non-decreasing
order. For a linear system this means `{b, m} = CoefficientArrays[eqs, vars]` gives the
constant vector `b` and coefficient matrix `m` directly — the first example returns
`{{-3, 1}, {{1, 2}, {3, -1}}}`.

A `List` of polynomials is a system, carrying the equation index as the leading axis; a
single polynomial (second example) carries no leading axis, so its degree-`d` slice is a
rank-`d` array. The arrays are dense `List`s, not `SparseArray` objects, and `Modulus -> p`
reduces the coefficients. A high degree in many variables makes the dense result explode and
is declined with a `CoefficientArrays::dense` message.
