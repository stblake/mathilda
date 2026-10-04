---
source: src/sort.c
---
**Algorithm.** `builtin_maximal_by` is `maximal_minimal_by(res, 0)` (mode 1 is
`MinimalBy`). It evaluates the key `f[e]` for every element `e` — for an
association, `f` of each *value* — into a parallel `keys[]` array, finds the best
key in one pass by `expr_compare` (greatest for `MaximalBy`, least for
`MinimalBy`), then makes a second pass collecting **every** element whose key
equals the best, in original order. So all ties are returned, and the result keeps
the collection's head (a `List` stays a list, an association returns the matching
entries). `MaximalBy[f]` with one argument is the operator form `MaximalBy[#1,
f] &`.

The comparison is Mathilda's canonical order (`expr_compare`), the same total
order `Sort` uses, so keys need not be numeric — any expressions are ranked.

**Data structures.** An owned `keys[]` array of the *n* evaluated key `Expr`s, and
an owned `out[]` array of the tied elements; both freed after the result node is
built. The input argument array is borrowed.

**Complexity / limits.** `O(n)` key evaluations and `O(n)` comparisons — it finds
an extreme, so it does **not** sort. The empty collection returns itself.
