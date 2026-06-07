---
source: src/poly/facpoly_factorterms.inc
---
**Algorithm.** `builtin_factortermslist` (in `src/poly/facpoly_factorterms.inc`) is the list-returning sibling of `FactorTerms`. It calls the same engine `ft_compute_list` but returns the factor `List` directly instead of multiplying it out: `{c_0, c_1, …, c_k, residue}`, where `c_0` is the numerical content, `c_1..c_k` are the contents extracted with respect to progressively smaller subsets of the requested variables, and `residue` is the remaining polynomial part. Unlike `FactorTerms` it does not auto-thread (its result shape is already a list). See `FactorTerms` for the content-extraction algorithm.

**Data structures.** `Expr*` factors gathered into a single `List`; content computation uses the recursive multivariate-GCD helper `ft_content_wrt_set`.
