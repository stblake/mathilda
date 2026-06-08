# SubresultantPolynomials

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SubresultantPolynomials[poly1, poly2, var]
    gives the list of subresultant polynomials {S_0, ..., S_m} of
    poly1 and poly2 with respect to var, where m = Exponent[poly2,
    var].  The list has length m + 1, its first element is
    Resultant[poly1, poly2, var], and the coefficient of var^j in S_j
    is the j-th principal subresultant coefficient.  Requires
    Exponent[poly1, var] >= Exponent[poly2, var] and exact
    coefficients.  Computed by a subresultant polynomial-remainder
    sequence, with a determinant-polynomial fallback for algebraic
    coefficients.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= SubresultantPolynomials[(x - 1)^2 (x - 2) (x - 3), (x - 1) (x - 4)^2, x]
Out[1]= {0, -36 + 36 x, 38 - 49 x + 11 x^2, -16 + 24 x - 9 x^2 + x^3}

In[2]:= SubresultantPolynomials[a x^3 + b x^2 + c x + d, 3 a x^2 + b x + c, x]
Out[2]= {4 a^2 c^3 + 2 a b^3 d - 18 a^2 b c d + 27 a^3 d^2, -2 a b c + 9 a^2 d - 2 a b^2 x + 6 a^2 c x, c + b x + 3 a x^2}

In[3]:= SubresultantPolynomials[2 x^7 + 3 x^3 + 5 x - 1, 7 x^6 + 8 x - 9, x]
Out[3]= {-183782157189, -761749829 + 3208696817 x, -3143546 + 11222638 x + 3838135 x^2, -21609 + 163611 x - 49392 x^2 + 64827 x^3, 0, -49 + 371 x - 112 x^2 + 147 x^3, -9 + 8 x + 7 x^6}

In[4]:= First[%] - Resultant[2 x^7 + 3 x^3 + 5 x - 1, 7 x^6 + 8 x - 9, x]
Out[4]= 0
```

## Implementation notes

**Algorithm.** `builtin_subresultantpolynomials` returns the list of subresultant *polynomials* `{S_0, ..., S_m}` of `poly1`, `poly2` w.r.t. `var`, where `m = Exponent[poly2, var]` (requires `deg p1 ≥ deg p2`, exact coefficients). `S_0` is `Resultant[p1, p2, var]`, and the coefficient of `var^j` in `S_j` is the j-th principal subresultant coefficient (matching `Subresultants[...][[j+1]]`).

By the fundamental theorem of subresultants each `S_j` is either zero or a scalar multiple of a single member of the subresultant polynomial remainder sequence. The implementation reuses the same Bronstein gamma/beta/delta PRS as `Subresultants` and `Resultant`, then classifies each output index: a *regular* index (`j == deg(R_p)` at a strict-drop step) gives `S_j = (psc_j / lc(R_p)) · R_p` — the chain member rescaled so its leading coefficient is `psc_j = clc_p^{δ_p}` (leaving `S_j = R_p` when `δ_p == 1`); a *defective* index (across a degree gap, `δ_p > 1`) is computed directly from the determinant-polynomial definition (small Sylvester minor since these sit high in the chain); all other indices are zero. For algebraic-number coefficients the whole list is built from the determinant-polynomial definition, mirroring `Resultant`/`Subresultants`.

**Data structures.** Descending coefficient arrays `Expr**` (`desc_coeffs`); the determinant-polynomial path builds truncated shifted polynomials via `trunc_shift_poly`. Output is a `List` of polynomial `Expr*`.

**Complexity / limits.** Same coefficient-growth control as the subresultant PRS; defective and algebraic-coefficient indices fall back to small-minor determinant evaluation.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- W. S. Brown and J. F. Traub, "On Euclid's Algorithm and the Theory of Subresultants", J. ACM 18(4), 1971.
- M. Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (Springer, 2005).
- Source: [`src/poly/subresultantpoly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/subresultantpoly.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
