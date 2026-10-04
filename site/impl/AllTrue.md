---
source: src/funcprog.c
---
**Algorithm.** `builtin_all_true` is `all_any_none_true(res, 0)`, the `mode = 0`
case of a shared quantifier (mode 1 is `AnyTrue`, mode 2 is `NoneTrue`). Three
paths, cheapest first. A bool `NDArray` with an identity-shaped predicate (`TrueQ`
or `Identity`, which hand a boolean element back unchanged) is a single
early-exit scan of the raw byte buffer — `np.all` — returning `False` on the first
unset byte. Otherwise the compiled-predicate fast path `pred_quantify` is tried,
then a visible `NDArray` is materialised and re-dispatched (`ndstruct_delist_repack`)
and an association is quantified over its values (`assoc_apply_over_values`).

The fallback walks the elements, building and evaluating `test[e]` for each. It
short-circuits — `AllTrue` returns `False` on the first `test[e]` that is `False`
— and if any result is neither `True` nor `False` the whole call is left
unevaluated (`return NULL`), matching Wolfram rather than guessing. An in-flight
`Throw` from the test propagates out immediately. With no short-circuit, the empty
and all-pass cases give `True`.

**Data structures.** The collection's argument array is borrowed; each `test[e]`
is a freshly built and evaluated `Expr` that is read for `True`/`False` and freed.
The bool-buffer path touches only the `unsigned char` NDArray buffer, so a 10⁶-element
`True`/`False` array never boxes a single symbol.

**Complexity / limits.** `O(n)` element tests with early exit; the bool scan is
`O(n)` over bytes with no allocation. An indeterminate predicate value leaves the
call unevaluated, so partially-symbolic input flows through unchanged.
