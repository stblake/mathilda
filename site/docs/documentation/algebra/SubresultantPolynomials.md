# SubresultantPolynomials

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SubresultantPolynomials[poly1, poly2, var]`**

gives the list of subresultant polynomials {S\_0, ..., S\_m} of poly1 and poly2 with respect to var, where m = Exponent\[poly2, var\].  The list has length m + 1, its first element is Resultant\[poly1, poly2, var\], and the coefficient of var^j in S\_j is the j-th principal subresultant coefficient.  Requires Exponent\[poly1, var\] \>= Exponent\[poly2, var\] and exact coefficients.  Computed by a subresultant polynomial-remainder sequence, with a determinant-polynomial fallback for algebraic coefficients.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= SubresultantPolynomials[(x - 1)^2 (x - 2) (x - 3), (x - 1) (x - 4)^2, x]
Out[1]= {0, -36 + 36 x, 38 - 49 x + 11 x^2, -16 + 24 x - 9 x^2 + x^3}

In[2]:= First[%] - Resultant[2 x^7 + 3 x^3 + 5 x - 1, 7 x^6 + 8 x - 9, x]
Out[2]= 183782157188
```

### Applications (3)

```mathematica
In[3]:= SubresultantPolynomials[x^2 - 1, x^2 - 4, x]
Out[3]= {9, -3, 1}

In[4]:= SubresultantPolynomials[x^3 + x + 1, x^2 + 1, x]
Out[4]= {1, 1, 1 + x^2}

In[5]:= Resultant[x^2 - 1, x^2 - 4, x]
Out[5]= 9
```

## Algorithm

==================================================================== subresultantpoly.c -- Subresultant polynomials.

```text
  SubresultantPolynomials[poly1, poly2, var]
    generates the list of subresultant polynomials {S_0, ..., S_m} of
    poly1 and poly2 with respect to var, where m = Exponent[poly2, var].
    The list has length m + 1, its first element is
    Resultant[poly1, poly2, var], and the coefficient of var^j in S_j is
    the j-th principal subresultant coefficient (i.e. Subresultants[
    poly1, poly2, var][[j+1]]).  Requires Exponent[poly1, var] >=
    Exponent[poly2, var] and exact coefficients.
```

Each subresultant polynomial S_j is, by the fundamental theorem of subresultants, either zero or a scalar multiple of a single member of

```text
the subresultant polynomial remainder sequence (PRS).  We therefore
```

reuse the Bronstein gamma/beta/delta PRS (the same recurrence as Subresultants in subresultants.c and Resultant in poly.c) and classify every output index:

```text
  * Regular index  (j == deg(R_p) for a strict-drop chain step p):
      S_j = (psc_j / lc(R_p)) * R_p, the chain member rescaled so its
      leading coefficient is the j-th principal subresultant coefficient
      psc_j = clc_p^{delta_p}.  (delta_p == 1 leaves S_j = R_p.)
  * Defective index (j == deg(R_{p-1}) - 1 across a degree gap,
      delta_p > 1): a lower-degree polynomial computed directly from the
      determinant-polynomial definition.  Such indices sit high in the
      chain, so the associated Sylvester minor is small.
  * Otherwise: S_j = 0 (gap interior, or below the last chain degree).
```

For algebraic-number coefficients (where the pseudo-remainder chain bloats) we skip the PRS and build the whole list from the determinant- polynomial definition, mirroring Resultant / Subresultants.

Memory convention matches the rest of Mathilda: every helper returning Expr* hands fresh ownership to the caller; builtin_subresultantpolynomials leaves its input `res` alive for the evaluator to free. ====================================================================

## Implementation notes

**Algorithm.** `builtin_subresultantpolynomials` returns the list of subresultant *polynomials* `{S_0, ..., S_m}` of `poly1`, `poly2` w.r.t. `var`, where `m = Exponent[poly2, var]` (requires `deg p1 ≥ deg p2`, exact coefficients). `S_0` is `Resultant[p1, p2, var]`, and the coefficient of `var^j` in `S_j` is the j-th principal subresultant coefficient (matching `Subresultants[...][[j+1]]`).

By the fundamental theorem of subresultants each `S_j` is either zero or a scalar multiple of a single member of the subresultant polynomial remainder sequence. The implementation reuses the same Bronstein gamma/beta/delta PRS as `Subresultants` and `Resultant`, then classifies each output index: a *regular* index (`j == deg(R_p)` at a strict-drop step) gives `S_j = (psc_j / lc(R_p)) · R_p` — the chain member rescaled so its leading coefficient is `psc_j = clc_p^{δ_p}` (leaving `S_j = R_p` when `δ_p == 1`); a *defective* index (across a degree gap, `δ_p > 1`) is computed directly from the determinant-polynomial definition (small Sylvester minor since these sit high in the chain); all other indices are zero. For algebraic-number coefficients the whole list is built from the determinant-polynomial definition, mirroring `Resultant`/`Subresultants`.

**Data structures.** Descending coefficient arrays `Expr**` (`desc_coeffs`); the determinant-polynomial path builds truncated shifted polynomials via `trunc_shift_poly`. Output is a `List` of polynomial `Expr*`.

**Complexity / limits.** Same coefficient-growth control as the subresultant PRS; defective and algebraic-coefficient indices fall back to small-minor determinant evaluation.

**Attributes:** `Protected`.

## References

**See also:** [Subresultants](../../algebra/Subresultants/)

- W. S. Brown and J. F. Traub, "On Euclid's Algorithm and the Theory of Subresultants", J. ACM 18(4), 1971.
- M. Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (Springer, 2005).
- Source: [`src/poly/subresultantpoly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/subresultantpoly.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
- Tests: [`tests/test_subresultantpolynomials.c`](https://github.com/stblake/mathilda/blob/main/tests/test_subresultantpolynomials.c)

## Notes & additional examples

### Notes

`SubresultantPolynomials[p1, p2, x]` returns the full chain of subresultant
*polynomials* `{S_0, ..., S_m}` with `m = Exponent[p2, x]`, so the list has
length `m + 1`. Its first element is `Resultant[p1, p2, x]`, and the coefficient
of `x^j` in `S_j` is the `j`-th principal subresultant coefficient (the value
returned by `Subresultants`). It requires `Exponent[p1, x] >= Exponent[p2, x]`
and exact coefficients, and is computed by a subresultant polynomial-remainder
sequence with a determinant-polynomial fallback for algebraic coefficients.
