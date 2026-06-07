---
source: src/poly/poly.c
---
**Algorithm.** `builtin_collect` (in `src/poly/poly.c`) groups an expression by powers of one or more keys, delegating to the recursive worker `collect_internal`. For each key it expands the expression with respect to that key (`expr_expand_patt`), except when the key is itself a `Plus` — there expansion would distribute the subterm and is skipped so e.g. `Collect[a(c+s)+b(c+s), c+s]` stays grouped. Each summand is decomposed into base–power form; terms are bucketed by the exponent at which the key appears (single-base keys group by the rational/symbolic exponent ratio, multi-factor monomial keys by the integer multiplicity from `get_k`). The collected coefficients are summed and, if a third head argument `h` is given, wrapped with `h`. It threads over `List`, equations and inequalities (skipping operator slots in `Inequality`), and recurses across multiple keys.

**Data structures.** Base–power lists (as in `Coefficient`) for term decomposition; results are rebuilt with `internal_times`/`Plus` and re-evaluated through `eval_and_free`.
