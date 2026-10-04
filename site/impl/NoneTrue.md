---
source: src/funcprog.c
---
**Algorithm.** `builtin_none_true` is `all_any_none_true(res, 2)`, the `mode = 2`
case of the quantifier it shares with `AllTrue` (mode 0) and `AnyTrue` (mode 1).
`NoneTrue` is the logical complement of `AnyTrue`: a bool `NDArray` with an
identity-shaped predicate (`TrueQ`/`Identity`) scans the raw byte buffer and
returns `False` if any byte is set, `True` otherwise. The compiled-predicate path
`pred_quantify`, visible-`NDArray` materialisation, and association-over-values
handling are all shared with the siblings.

The fallback walks the elements building and evaluating `test[e]`;
`NoneTrue` short-circuits to `False` on the first `test[e]` that is `True`. A
result that is neither `True` nor `False` leaves the whole call unevaluated
(`return NULL`), and a `Throw` from the test propagates. With no match the result
is `True`, so the empty list gives `True`.

**Data structures.** The argument array is borrowed; each `test[e]` is built,
evaluated, read for `True`/`False`, and freed. The bool-buffer path reads only the
`unsigned char` NDArray buffer.

**Complexity / limits.** `O(n)` tests with early exit; the bool scan is `O(n)`
bytes. An indeterminate predicate value leaves the call unevaluated, matching
Wolfram.
