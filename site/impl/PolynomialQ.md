---
source: src/poly/poly.c
---
**Algorithm.** `builtin_polynomialq` is a structural predicate. It normalises the second
argument into a variable set (a single symbol or a `List` of symbols) and calls
`is_polynomial`, which recurses over the expression tree: an expression is a polynomial in the
given variables iff every node is one of the variables, a sub-expression free of all the
variables (a degree-0 constant), a `Plus`/`Times` whose arguments are all polynomials, or a
`Power[base, k]` with a non-negative integer exponent `k` and polynomial base. Anything else
containing a variable in a non-polynomial position (e.g. `Sin[x]`, `1/x`, `x^(1/2)`) returns
`False`. Returns the symbol `True` or `False`.
