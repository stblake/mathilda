---
source: src/funcprog.c
---
**Algorithm.** `builtin_tuples` (in `src/funcprog.c`) enumerates the Cartesian product. `Tuples[{l1, ..., lk}]` takes one tuple from each list; `Tuples[list, n]` is the n-ary product of one list with itself; `Tuples[list, {n1, ..., nd}]` produces n1·…·nd-element tuples reshaped into a d-dimensional array. All three normalize to an array of source-list pointers and call the recursive `tuples_rec`.

`tuples_rec` is a straightforward odometer: it iterates the current list's elements in order, recursing one position deeper for each, and emits a completed tuple when all positions are filled. This yields **lexicographic order with the last list varying fastest** (the standard Cartesian-product order). The tuple head is taken from the input expression's head (so `Tuples` over non-`List` heads is structure-preserving). For the `{n1,...,nd}` form, each flat tuple is then folded into nested sublists by `reshape_rec`. Results grow in a geometrically resized `Expr**` buffer; an empty source list short-circuits to no output.
