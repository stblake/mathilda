---
source: src/funcprog.c
---
`builtin_select` filters the arguments of a compound expression by a predicate.
It iterates the args of `list` (any head, not only `List`), and for each element
builds `crit[elem]` and runs `evaluate()`; the element is kept only when the
result is exactly the symbol `True`. The optional third argument caps the number
of kept elements (`n_max`), stopping the scan early once reached. The surviving
elements are reassembled under the original head via `expr_new_function`. Returns
`NULL` (unevaluated) when the first argument is an atom or when the count
argument is non-integer. Each predicate test allocates a copied call and frees it
plus its evaluated result, so memory is bounded per element.
