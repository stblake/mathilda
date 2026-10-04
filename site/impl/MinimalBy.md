---
source: src/sort.c
---
**Algorithm.** `builtin_minimal_by` is `maximal_minimal_by(res, 1)`, the mode-1
twin of `MaximalBy`. It evaluates the key `f[e]` for every element (for an
association, `f` of each *value*) into a parallel `keys[]` array, finds the least
key in one pass by `expr_compare`, then collects in a second pass every element
whose key equals that minimum, in original order — so all ties are returned. The
result keeps the input's head (a list stays a list; an association returns the
matching entries). `MinimalBy[f]` is the operator form.

The ranking is Mathilda's canonical order (`expr_compare`), so the keys may be any
expressions, not only numbers.

**Data structures.** An owned `keys[]` array of the *n* evaluated key `Expr`s and
an owned `out[]` array of the tied elements, both freed once the result is built.
The input argument array is borrowed.

**Complexity / limits.** `O(n)` key evaluations and `O(n)` comparisons; it selects
the minimum rather than sorting. The empty collection returns itself.
