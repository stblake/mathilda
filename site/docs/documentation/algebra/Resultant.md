# Resultant

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Resultant[p, q, var]
    gives the resultant of p and q as polynomials in var: the unique
    integer / polynomial scalar that vanishes iff p and q share a
    root in var.  Computed via a Sylvester-matrix determinant or, in
    the exact path, a subresultant pseudo-remainder sequence.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Resultant[x^2 - 2x + 7, x^3 - x + 5, x]
Out[1]= 265

In[2]:= Resultant[x^3 - 5x^2 - 7x + 14, x^3 - 8x^2 + 9x + 58, x]
Out[2]= 0
```

## Implementation notes

**Algorithm.** `builtin_resultant` first verifies both arguments are polynomials in the
given variable (via `internal_polynomialq`), then delegates to `resultant_internal`. That
routine factors the problem multiplicatively: a `Times` or `Power` argument is split using
`Resultant(fg, h) = Resultant(f, h) Resultant(g, h)` and `Resultant(f^k, h) =
Resultant(f, h)^k`, recursing on each factor. After expanding both inputs it handles the
degenerate degree-0 cases (`Resultant(c, q) = c^deg(q)`, etc.) directly.

For the general case it tries Bronstein's **subresultant polynomial remainder sequence**
(`resultant_subresultant`) first. The chain works entirely in the coefficient ring `D[x]`,
using `pseudo_rem_standard` (the full `lc(B)^(deg A − deg B + 1)` pseudo-remainder Bronstein's
identities require) plus a single *scalar* exact division by `β_i ∈ D` per step. This makes
only `O(min(n, m))` pseudo-remainder calls. It deliberately bails out (returns NULL) when the
inputs carry algebraic-number coefficients — `subres_has_algebraic` detects `Power[base,
Rational[a,b]]` (b > 1) subterms such as `Sqrt[N]`, which cause the PRS to bloat geometrically
— and on a size-budget overflow.

When the subresultant path declines, `resultant_internal` falls back to the **Sylvester
matrix** definition: it bulk-extracts all coefficients (`get_all_coeffs_expanded`), builds the
`(n+m)×(n+m)` matrix as a `List` of `List` rows (m shifted copies of P's coefficients atop n
shifted copies of Q's), and computes the determinant by constructing a `Det[...]` call and
running it through `evaluate`, then expands the result. The Sylvester path needs `O(n^3)` exact
polynomial-coefficient divisions in `Det`'s Bareiss elimination, collapsing to Laplace
expansion over rings (like `Q(α)`) where those divisions can't be certified — which is exactly
why the subresultant PRS is preferred.

**Data structures.** Polynomials are ordinary `Expr` trees in expanded form; coefficients are
extracted into `Expr**` arrays. The Sylvester matrix is materialised as nested `List`
expressions and handed to the linalg `Det` builtin.

`Discriminant` reuses `resultant_internal` as `(-1)^(n(n-1)/2)/a_n · Resultant(p, p', x)`.

- `Protected`, `Listable`.
- Computes the resultant of polynomials `poly1` and `poly2` with respect to the variable `var`.
- The resultant is independent of common roots and vanishes exactly when the polynomials have roots in common.
- Default algorithm is Bronstein's subresultant PRS (Symbolic Integration I, p.24): a linear chain of pseudo-remainders with scalar exact divisions in the coefficient ring, avoiding the (n+m)x(n+m) Sylvester matrix construction and its O(n^3) Bareiss reduction.  For Z/Q coefficients this is materially faster than the matrix path, and on inputs with symbolic coefficients it sidesteps the O(n!) Laplace expansion that the matrix path falls back to when Bareiss exact-division certification fails.
- Inputs containing algebraic-number coefficients (e.g. `Sqrt[N]`, cube roots — any `Power[X, Rational[a,b]]` with `b > 1`) are routed to the Sylvester+Det path instead, because the subresultant chain bloats geometrically when `Power[base, k/m]` forms can't be combined with their `Times[base^q, Sqrt[base]]` equivalents by `Plus` alone.
- A size-budget guard inside the subresultant path falls back to Sylvester+Det for any pathological input where chain elements exceed ~30x the input leaf-count.
- Automatically preserves multiplicativity (e.g., $Res(A \cdot B, Q) = Res(A, Q) Res(B, Q)$ and $Res(A^k, Q) = Res(A, Q)^k$).
- **FLINT acceleration** (when built with FLINT): plain rational inputs are computed via `fmpq_mpoly_resultant`, which avoids the subresultant PRS coefficient growth on higher-degree / multivariate inputs. FLINT's convention matches the classical output exactly; inputs with algebraic-number or otherwise non-rational coefficients fall through to the paths above. The same kernel is exposed directly as `` FLINT`Resultant `` (see the FLINT` context section).

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- von zur Gathen & Gerhard, "Modern Computer Algebra" (3rd ed.), Ch. 6 (resultants and the Sylvester matrix).
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), Ch. 7 (subresultant PRS).
- M. Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (Springer, 2005) — SubResultant, p. 24.
- W. S. Brown and J. F. Traub, "On Euclid's Algorithm and the Theory of Subresultants", JACM 18(4), 1971.
- J. J. Sylvester, dialytic elimination / the Sylvester matrix.
- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Resultant[x^2 - 1, x^2 - 4, x]
Out[1]= 9
```

```mathematica
In[1]:= Resultant[x^2 - 2, x^2 - 3, x]
Out[1]= 1
```

```mathematica
In[1]:= Resultant[x^2 + a, x + b, x]
Out[1]= a + b^2
```

```mathematica
In[1]:= Resultant[x^2 - y, x^2 + y, x]
Out[1]= 4 y^2
```

```mathematica
In[1]:= Resultant[x^2 + a x + b, 2 x + a, x]
Out[1]= -a^2 + 4 b
```

```mathematica
In[1]:= Resultant[x^3 + p x + q, 3 x^2 + p, x]
Out[1]= 4 p^3 + 27 q^2
```

```mathematica
In[1]:= Factor[Resultant[x^2 + y^2 - 1, x + y - 1, x]]
Out[1]= 2 y (-1 + y)
```

### Notes

`Resultant[p, q, x]` returns a scalar in the remaining variables that
vanishes exactly when `p` and `q` share a common root in `x`. For
`x^2 - 1` and `x^2 - 4` (roots `±1` and `±2`) the value is the nonzero `9`,
confirming no shared root, whereas eliminating `x` from `x^2 - y` and
`x^2 + y` gives `4 y^2`, which is zero only at the shared-root locus `y = 0`.
The computation uses either a Sylvester-matrix determinant or, on the exact
path, a subresultant pseudo-remainder sequence; coefficients may themselves
be polynomials in the other variables, as in `a + b^2`.
