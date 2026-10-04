---
source: src/core.c
---
**Algorithm.** `builtin_valueq` (`src/core.c`) inspects the symbol table rather
than evaluating its argument — `ValueQ` carries `HoldAll`, so the argument reaches
the builtin unevaluated. A bare `EXPR_SYMBOL` is valued iff its `SymbolDef`
carries a non-empty `own_values` chain (an immediate or delayed assignment). A
compound `f[...]` is valued iff the `SymbolDef` of its head `f` carries any
`down_values` — regardless of whether the particular arguments would match a rule.
Everything else (numbers, strings, undefined symbols) returns `False`.

**Data structures.** Two `symtab_lookup` probes at most, reading the `own_values`
/ `down_values` linked lists hanging off the `SymbolDef`; no rule matching and no
tree walk.

**Complexity / limits.** `O(1)` (a hash lookup plus a null test on a list head).
Because the test is existence-of-a-rule, not applicability, `ValueQ[f[a, b]]` is
`True` as soon as `f` has any `DownValue`. Wrong arity routes a `ValueQ::argx`
diagnostic through `builtin_arg_error`. Attributes `HoldAll`, `Protected`.
