---
source: src/core.c
---
`builtin_own_values` (1-arg, symbol only) calls `symtab_get_own_values(sym)` to retrieve the symbol's `Rule*` list — the immediate `x = value` assignments — and hands it to the shared `rules_to_list` helper, which deep-copies each `pattern`/`replacement` pair into a `Rule[lhs, rhs]` node and gathers them into a `List`. The list is returned unevaluated so the stored values are reported verbatim rather than re-evaluated.
