---
source: src/poly/facpoly_factorterms.inc
---
**Algorithm.** `builtin_factorterms` (in `src/poly/facpoly_factorterms.inc`, compiled into `facpoly.c`) pulls out the content of a polynomial without factoring the polynomial part. It threads over `List`/equation/inequality/logic heads (`ft_is_threading_head`), then calls the shared engine `ft_compute_list` and multiplies the resulting factor list back into a single `Times`.

`ft_compute_list` first `Together`-normalises and splits into numerator/denominator. It collects and sorts the numerator's variables, then: (1) extracts the **numerical content** via `ft_content_wrt_set` (content with respect to *all* variables over an empty ground ring, i.e. the integer GCD of coefficients) and divides it out with `ft_divide_out` (exact polynomial division, falling back to symbolic `Times[poly, content^{-1}]`); (2) for a 2-arg call `FactorTerms[poly, {x_1,…,x_k}]`, peels the content with respect to progressively smaller variable subsets, where `ft_content_wrt_set` recursively computes the multivariate `poly_gcd_internal` of the coefficients of each monomial in the chosen variables over the shrinking ground ring; (3) appends the final residue (re-multiplied by `1/den` to round-trip rational inputs).

**Data structures.** `Expr*` polynomials throughout; `poly_gcd_internal` (multivariate polynomial GCD) is the workhorse for symbolic content; variable lists are `Expr**` sorted with `compare_expr_ptrs`.
