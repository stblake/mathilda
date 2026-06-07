---
references:
  - "M. Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (Springer, 2005) — SubResultant, p. 24."
  - "W. S. Brown and J. F. Traub, \"On Euclid's Algorithm and the Theory of Subresultants\", JACM 18(4), 1971."
  - "J. J. Sylvester, dialytic elimination / the Sylvester matrix."
source: src/poly/poly.c
---
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
