---
source: src/funcprog.c
---
**Algorithm.** `builtin_any_true` is `all_any_none_true(res, 1)`, sharing one
quantifier with `AllTrue` (mode 0) and `NoneTrue` (mode 2). A bool `NDArray` with
an identity-shaped predicate (`TrueQ`/`Identity`) is answered by scanning the raw
byte buffer for any set byte — `np.any`. Otherwise the compiled-predicate path
`pred_quantify` runs, then a visible `NDArray` is materialised and re-dispatched
(`ndstruct_delist_repack`) and an association is quantified over its values.

The fallback walks the elements building and evaluating `test[e]`; `AnyTrue`
short-circuits to `True` on the first `test[e]` that is `True`. A result that is
neither `True` nor `False` makes the whole call stay unevaluated (`return NULL`),
and a `Throw` from the test propagates. If nothing matches, the result is `False`
— so the empty list gives `False`.

**Data structures.** The argument array is borrowed; each `test[e]` is a built,
evaluated, then freed `Expr`. The bool-buffer path reads only the `unsigned char`
NDArray buffer, closing the sign-predicate → quantifier pipeline without
materialising a `True`/`False` symbol per element.

**Complexity / limits.** `O(n)` tests with early exit (`O(1)` when an early
element matches); the bool scan is `O(n)` bytes, no allocation. An indeterminate
predicate value leaves the call unevaluated.
