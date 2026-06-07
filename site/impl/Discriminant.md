---
source: src/poly/poly.c
references:
  - "B. L. van der Waerden, *Algebra*, vol. 1 (Springer)."
---
**Algorithm.** `builtin_discriminant` (in `src/poly/poly.c`) computes `Disc_x(p)` from the standard resultant identity. After confirming `p` is a polynomial in `x` (`PolynomialQ`) and expanding it, it reads the degree `n` (returning 0 for degree 0/1). It then forms the derivative `p'` with `poly_derivative`, computes `R = Res_x(p, p')` via `resultant_internal`, and applies the closed form

  Disc = (-1)^{n(n-1)/2} · Res(p, p') / a_n,

where `a_n` is the leading coefficient (obtained with `get_coeff`). The quotient is built as `Times[sign·R, a_n^{-1}]` and the final result is `expr_expand`ed.

**Data structures.** Plain `Expr*` polynomials throughout; the heavy lifting is the subresultant/Euclidean `resultant_internal` (see `Resultant`). Degree and coefficient queries use `get_degree_poly`/`get_coeff`.

**Complexity / limits.** Dominated by the resultant computation, which is polynomial in `n` and the coefficient sizes but can be expensive for large multivariate inputs.
