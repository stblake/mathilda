---
references:
  - "K. O. Geddes, S. R. Czapor, G. Labahn, *Algorithms for Computer Algebra* (Kluwer, 1992) — polynomial remainder sequences / subresultants."
  - "W. S. Brown, J. F. Traub, \"On Euclid's Algorithm and the Theory of Subresultants\", JACM 18 (1971)."
source: src/poly/poly.c
---
**Algorithm.** `SubresultantPolynomialRemainders[a, b, x]` returns the
polynomial remainder chain `{a, b, R_2, R_3, ...}` in `K(coeffs)[x]`.
`builtin_subresultantpolynomialremainders` (src/poly/poly.c) expands both
inputs, swaps them if needed so `deg a ≥ deg b` (the documented chain
orientation), seeds the chain with `{a, b}`, and then repeatedly appends
`R = pseudo_rem(prev, cur, x)` — the **pseudo-remainder** — until the remainder
is zero or constant (degree ≤ 0). The chain array grows by doubling.

Note (per the source comment) this is a **pseudo-remainder PRS, not the
Lazard-scaled subresultant chain**: it does not divide out the subresultant
content scaling factors. This is deliberate — its sole consumer, the
Lazard–Rioboo–Trager log part of rational integration (src/calculus/intrat.c),
only uses each chain element's degree in `x` and its primitive part in the
auxiliary variable `t`, and both are invariant under content scaling, so the
simpler pseudo-remainder chain is a correct, cheaper substrate.

**Data structures.** A dynamically grown `Expr**` array of polynomial trees,
returned wrapped in a `List`. Each step relies on `pseudo_rem`,
`get_degree_poly`, `is_zero_poly`, and `expr_expand`.

**Complexity / limits.** Length of the chain is `O(deg b)` reductions; because
content is not removed, intermediate coefficients can swell relative to a true
subresultant PRS, but element degrees and `t`-primitive parts (all the consumer
needs) are exact. Requires the variable to be a symbol.
