---
source: src/funcprog.c
---
**Algorithm.** `builtin_select_first` returns the first element of a collection
for which a predicate holds, the one-element analogue of `Select`. A
compiled-predicate fast path (`pred_find_first`) scans a packed real buffer and
returns the first passing element directly, or `PRED_NO_MATCH` to fall to the
default/`Missing` branch — cheap when a match comes early, a full scan otherwise.
A visible `NDArray` is materialised and re-dispatched (`ndstruct_delist_repack`);
an association is handled over its values (`assoc_apply_over_values`).

The fallback loops the elements, building and evaluating `pred[e]`, and returns a
copy of the first element whose test is `True`. A `Throw` from the predicate
propagates. If nothing matches, `SelectFirst[list, pred, default]` returns
`default` and the two-argument form returns `Missing["NotFound"]`.

**Data structures.** The argument array is borrowed; each `pred[e]` is a freshly
built, evaluated, then freed `Expr`, and the match is a copy of the stored
element. No result list is accumulated — the scan stops at the first hit.

**Complexity / limits.** `O(1)` when an early element matches, `O(n)` predicate
evaluations when the match is late or absent.
