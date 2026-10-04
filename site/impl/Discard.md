---
source: src/assoc_ops.c
---
**Algorithm.** `builtin_discard` is the complement of `Select`: it keeps the
elements of a non-atomic expression for which the criterion does *not* give
`True`. `Discard[expr, crit]` tests every element; `Discard[expr, crit, n]`
discards at most the first `n` matches and keeps the rest regardless. The result
reuses the original head, so it works on any non-atomic expression, and on an
association it tests the values while keeping the keys.

**Data structures.** A single pass builds a `keep` array of element copies; the
criterion is applied with `eval_call1` (evaluate `crit[elem]`) and the verdict
read with `is_true_sym`. A visible `NDArray` is de-listed first via
`ops_delist_visible`.

**Complexity / limits.** `O(n)` evaluations of the criterion. A non-atomic first
argument is required (`Discard::normal`), and the optional limit must be a
non-negative integer or `Infinity` (`Discard::innf`); both diagnostics route
through the message funnel so `Quiet[]`/`Check[]` behave. `crit` is arbitrary, so
there is no packed fast path.
