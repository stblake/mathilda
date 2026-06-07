---
references:
  - "D. Kozen and S. Landau, \"Polynomial Decomposition Algorithms\", J. Symbolic Computation 7 (1989)."
source: src/poly/poly.c
---
**Algorithm.** `builtin_decompose` (`src/poly/poly.c`) finds polynomials with f(x) = g(h(x)), deg(h) ≥ 2. It first verifies the input is a polynomial in x via `internal_polynomialq` (returns `NULL`/unevaluated otherwise) and then calls the recursive driver `decompose_recursive`. That driver expands f, takes its degree n, and applies two reductions in order: (1) a common-degree shortcut — if the gcd `d` of the degrees of all nonzero monomials exceeds 1, every term is divisible by x^d, so it substitutes y = x^d and recurses; (2) trial composition (the Kozen–Landau scheme) — for each proper divisor s of n it builds the candidate inner polynomial of degree s by matching the high-order coefficients of `a_n*(x^s)^r` against f (solving for each coefficient with `exact_poly_div` over the non-x variables), then divides f by successive powers of the candidate inner h with `poly_div_rem`, accepting the decomposition only if every remainder is x-free and the final quotient is zero. The outer factor list `Lg` from the recursion is concatenated with the inner h to form the returned `List`. If no decomposition is found it returns `{f}`.

**Data structures.** `Expr*` polynomial trees throughout; coefficients are extracted with `get_coeff`/`get_all_coeffs_expanded` and variables collected via `collect_variables`. Exact divisions use the multivariate `exact_poly_div` and `poly_div_rem` helpers in the same file.

**Complexity / limits.** Univariate decomposition over the polynomial's coefficient domain; multivariate coefficients are tolerated through the variable-list machinery. Cost is dominated by repeated `expr_expand` and `poly_div_rem` over the divisors of n.
