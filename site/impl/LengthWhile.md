---
source: src/funcprog.c
---
**Algorithm.** `builtin_lengthwhile` returns the length of the leading run of
elements for which `crit` gives `True`, stopping at the first failure. It tries a
compiled-predicate fast path first: `pred_run_length` compiles `crit` and scans a
machine buffer directly, so for a numeric predicate nothing is boxed and the answer
is just the count. Failing that, a visible `NDArray` is materialised and repacked
(`ndstruct_delist_repack`), and otherwise `leading_run_length` walks the elements one
by one (for an association it tests the *values*).

**Data structures.** A packed numeric buffer on the compiled path; plain `Expr`
element pointers on the interpreted path. Only a count is returned, so nothing is
copied.

**Complexity / limits.** O(k) where `k` is the run length — the scan stops at the
first element that fails `crit`, so a predicate that fails immediately is O(1).
`TakeWhile` returns the run itself rather than its length.
