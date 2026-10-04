---
source: src/ndarray.c
---
**Algorithm.** `builtin_ndarrayq` (`src/ndarray.c`) is a one-argument predicate: it returns
`True` when `is_ndarray(arg)` holds and `False` otherwise. `is_ndarray` recognises *both*
packed-array surfaces Mathilda uses — a visible `NDArray[...]` object and a transparently
packed `List` (the dense machine-number buffer that `$AutoArrayPacking` builds behind an
ordinary list head) — so `NDArrayQ` is the one test that sees through the transparency gate.

**Data structures.** None beyond the input `Expr`. The function reads the argument's tag and
packed-buffer metadata through `is_ndarray` and constructs a fresh `True`/`False` symbol.
It is registered `Protected` in `ndarray_init`; a non-unary call returns `NULL`. Its sibling
`PackedArrayQ` tests only the packed-`List` surface (`is_packed_list`), so a visible
`NDArray[...]` makes `PackedArrayQ` `False` while `NDArrayQ` is `True`.

**Complexity / limits.** `O(1)` — it inspects the representation, not the contents, so it
never materialises a packed buffer. Because a packed list is otherwise indistinguishable
from a plain `List` (same head, printed form, elements, ordering, matches), `NDArrayQ` is
the intended way to confirm that a value is actually on the buffer fast path: a *small*
plain list that was never packed gives `False`, while `Range[10000]` under the default
`$AutoArrayPacking` gives `True`.
