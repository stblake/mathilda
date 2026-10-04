---
source: src/poly/monomials.c
---
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
