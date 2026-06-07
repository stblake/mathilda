---
source: src/purefunc.c
---
`Slot[n]` (`#`, `#n`) is an inert marker — `builtin_slot` always returns `NULL` so the node stays unevaluated on its own. Substitution happens only when a pure `Function` is applied: `substitute_slots` walks the function body and replaces each `Slot[n]` with a copy of the `n`-th argument (when `1 ≤ n ≤ arg_count`), stopping recursion at any nested `Function` so inner pure functions are not captured by the outer slots. `ATTR_PROTECTED`.
