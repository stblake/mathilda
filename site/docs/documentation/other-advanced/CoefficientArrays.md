# CoefficientArrays

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`CoefficientArrays[polys, vars] gives the coefficient arrays of the polynomial system polys in the variables vars, grouped by total degree: element d + 1 holds every degree-d coefficient, so {b, m} = CoefficientArrays[eqs, vars] is the constant vector and coefficient matrix of a linear system. A degree-d coefficient sits at the index tuple naming its variables in non-decreasing order, every other permutation being 0; a list of polynomials carries the equation index as the leading axis. CoefficientArrays[poly, vars] treats a non-List first argument as a single polynomial. Option Modulus -> p reduces the coefficients. The arrays are dense Lists rather than SparseArrays.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

{constant vector, coefficient matrix}

```mathematica
In[1]:= CoefficientArrays[{x + 2 y - 3, 3 x - y + 1}, {x, y}]
Out[1]= {{-3, 1}, {{1, 2}, {3, -1}}}
```

Grouped by degree: {const, {linear}, {{quadratic}}}

```mathematica
In[2]:= CoefficientArrays[1 + 2 x + 3 x^2, x]
Out[2]= {1, {2}, {{3}}}
```

## Implementation notes

**Algorithm.** `builtin_coefficientarrays` (`src/poly/monomials.c`) decomposes each
polynomial into `(exponent-vector, coefficient)` terms (the shared monomial front end also
used by `MonomialList`/`CoefficientRules`), reduces coefficients modulo an optional
`Modulus -> m`, finds the maximum total degree `D`, and lays the terms into `D + 1` nested
**dense** arrays grouped by total degree: result element `d + 1` holds every degree-`d`
coefficient. A degree-`d` coefficient sits at the rank-`d` index tuple naming its variables
in non-decreasing order, every other permutation of that tuple being `0`. So `{b, m} =
CoefficientArrays[eqs, vars]` is the constant vector `b` and the coefficient matrix `m` of a
linear system. A `List` first argument is treated as a **system** — the equation index
becomes the leading axis, raising each array's rank by one and forcing one shared column
order across equations; a non-`List` argument is a single polynomial with no leading axis.

**Data structures.** The arrays are built as nested dense `List`s of `Expr` (zero-filled by
`ca_make_zero`, then overwritten by `ca_set` at a computed multi-index), not `SparseArray`
objects — `Normal` is the identity on them. Any visible `NDArray`/packed-`List` argument is
unpacked first, since this head is not on the `AWARE` list.

**Complexity / limits.** The dense layout costs `O(nvars^d)` entries at degree `d`, so a
high degree in many variables blows up; the routine guards this with a `CoefficientArrays
::dense` message (routed through `mth_message`) and declines when the dense result would
exceed its entry cap. `FromCoefficientRules` is the structural inverse of the
`CoefficientRules` form rather than of these arrays.

**Attributes:** `Protected`.

## References

- Source: [`src/poly/monomials.c`](https://github.com/stblake/mathilda/blob/main/src/poly/monomials.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

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
