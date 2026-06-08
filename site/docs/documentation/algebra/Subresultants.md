# Subresultants

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Subresultants[poly1, poly2, var]
    gives the list of principal subresultant coefficients of poly1 and
    poly2 with respect to var.  The list has length
    Min[Exponent[poly1, var], Exponent[poly2, var]] + 1, its first
    element is Resultant[poly1, poly2, var], and the first k entries
    vanish exactly when the polynomials share k roots (with
    multiplicity).  Computed by a subresultant polynomial-remainder
    sequence, or a Sylvester-minor determinant for algebraic
    coefficients.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Subresultants[2x^7 + 3x^3 - 7x + 1, 3x^5 - 17x + 21, x]
Out[1]= {273612691817, 68946901, 1299537, 16641, 0, 9}

In[2]:= Subresultants[(x - 1)(x - 2)(x - 3), (x - 1)(x - 4)(x - 5), x]
Out[2]= {0, 12, -4, 1}

In[3]:= Subresultants[(x - 1)^5 (x - 2)(x - 3), (x - 1)^4 (x - 4)(x - 5), x]
Out[3]= {0, 0, 0, 0, 144, 18, 1}

In[4]:= Subresultants[a x^3 + b x^2 + c x + d, x^3 - 5 b x - 7 a, x]
Out[4]= {-343 a^6 + 245 a^4 b^2 - 49 a^2 b^3 + 147 a^3 b c - 175 a^3 b^2 c + 35 a b^3 c - 70 a^2 b c^2 - 7 a c^3 - 147 a^4 d - 35 a^2 b^2 d + 125 a^2 b^3 d - 25 b^4 d + 21 a b c d + 50 a b^2 c d + 5 b c^2 d - 21 a^2 d^2 - 10 b^2 d^2 - d^3, -7 a^2 b + 25 a^2 b^2 - 5 b^3 + 10 a b c + c^2 - b d, -b, 1}

In[5]:= Length[Subresultants[x^50 + a, x^20 + b, x]]
Out[5]= 21
```

## Implementation notes

**Algorithm.** `builtin_subresultants` returns the list of principal subresultant coefficients (PSCs) of `poly1`, `poly2` w.r.t. `var`. The list has length `min(deg p1, deg p2) + 1`; element 0 equals `Resultant[p1, p2, var]`, and the first k PSCs vanish exactly when the polynomials share k roots (with multiplicity). Both arguments are checked polynomial in `var` (`internal_polynomialq`), then expanded.

The efficient path is `subresultants_prs`: the Bronstein subresultant polynomial remainder sequence with the gamma/beta/delta recurrence (`gamma_i = (-lc_{i-1})^{δ_{i-1}} γ_{i-1}^{1-δ_{i-1}}`, `beta_i = -lc_{i-1} γ_i^{δ_i}`, `r_{i+1} = prem(r_{i-1}, r_i)/β_i`) — the same recurrence as `Resultant` in `src/poly/poly.c`. It orients so `deg A ≥ deg B`, runs the chain via `pseudo_rem_standard` and `poly_divide_by_scalar`, retains each member's degree and leading coefficient, and extracts PSCs by the fundamental theorem of subresultants using the cumulative recurrence `s_p = lc(R_p)^{δ_p} / s_{p-1}^{δ_p-1}`, placing each at its degree index. Equal input degrees set `psc_L = 1` (empty top minor); a swap-correction sign `(-1)^{(n-j)(m-j)}` is applied when the inputs were oriented. Intermediate gamma/beta and PSC arithmetic goes through `internal_cancel` + `expr_expand`.

For algebraic-number coefficients (Sqrt, cube roots — `subres_has_algebraic`), where the pseudo-remainder chain bloats, it falls back to `subresultants_determinant`: `PSC_j = Det(M_j)`, where `M_j` is the Sylvester matrix restricted to the first `m-j` `poly1`-shift rows, `n-j` `poly2`-shift rows, and `n+m-2j` columns, evaluated via the `Det` builtin and `expr_expand`. The same determinant path is the fallback when the PRS path returns NULL.

**Data structures.** Descending coefficient arrays `Expr**` (`desc_coeffs`, leading coeff at index 0), built via the bulk `get_all_coeffs_expanded` extractor when possible. The PRS keeps parallel growable arrays `R[]` (chain members), `deg[]`, `lc[]`. Matrices are `List`-of-`List` `Expr*`.

**Complexity / limits.** The subresultant PRS keeps coefficient growth polynomial (the whole point of the Bronstein recurrence vs. naive PRS). The determinant fallback is `O(min(n,m))` minors, each a Sylvester-minor determinant. Identically-zero input is left unevaluated.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- W. S. Brown and J. F. Traub, "On Euclid's Algorithm and the Theory of Subresultants", J. ACM 18(4), 1971.
- M. Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (Springer, 2005).
- Source: [`src/poly/subresultants.c`](https://github.com/stblake/mathilda/blob/main/src/poly/subresultants.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
