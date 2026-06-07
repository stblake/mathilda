---
source: src/poly/facpoly_squarefree.inc
references:
  - "D. Y. Y. Yun, \"On square-free decomposition algorithms\", SYMSAC 1976."
  - "K. O. Geddes, S. R. Czapor, G. Labahn, *Algorithms for Computer Algebra* (Kluwer, 1992)."
---
**Algorithm.** `builtin_factorsquarefree` (in `src/poly/facpoly_squarefree.inc`, compiled into `facpoly.c`) computes the square-free decomposition `p = ∏ p_i^i` with each `p_i` square-free and pairwise coprime, using **Yun's algorithm**. After a `Rationalize`/`Numericalize` round-trip for inexact inputs, the dispatcher `factor_square_free_dispatcher` expands and routes to `factor_square_free_poly`.

`factor_square_free_poly` works recursively on the variable list. It first extracts the polynomial content (`poly_content`) and recurses on it; on the primitive part `pp` it runs Yun's recurrence: set `A = pp`, `B = gcd(A, A')`, `C = A/B`, `D = A'/B − C'`, then iterate `P_i = gcd(C, D)` (the degree-`i` square-free factor), `C ← C/P_i`, `D ← D/P_i − (C/P_i)'`, accumulating `P_i^i`. Derivatives are computed by the local `poly_deriv` and GCDs by `poly_gcd_internal`; exact divisions by `exact_poly_div`. An F4 "cheap squarefree" pre-check (`sqfree_cheap_check`, a univariate-substitution probe) short-circuits the expensive `gcd(pp, pp')` when `pp` is provably square-free. A final `exact_poly_div` recovers any leading/missing scalar factor so the product round-trips to the input.

**Data structures.** `Expr*` polynomials manipulated through the generic evaluator (`eval_and_free`) and the poly helpers (`get_degree_poly`, `poly_gcd_internal`, `exact_poly_div`); factors collected in a growable `Expr**` buffer with their multiplicities encoded as `Power[P_i, i]`.

**Complexity / limits.** Dominated by the GCD chain; the cheap pre-check avoids one full multivariate GCD on the common square-free case. Square-free decomposition is a prerequisite stage for full `Factor`.
