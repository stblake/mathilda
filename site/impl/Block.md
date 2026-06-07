---
source: src/modular.c
---
**Algorithm.** `builtin_block` (in `src/modular.c`) implements **dynamic scoping** by save/restore of global symbol state — it does *not* rename anything (contrast `Module`). `Block` carries `HoldAll | Protected` (set in `src/attr.c`). For each local `x` (or `x = init`) the handler looks up the symbol's `SymbolDef`, saves its current `own_values` chain and `attributes` into a `SavedVar` record, then clears `own_values` to `NULL` so the name is unbound inside the block. If an initializer is given it is evaluated (in the surrounding scope) and installed as a fresh `OwnValue` on the *same* global symbol.

The body is then `evaluate`d directly — any reference to `x` anywhere in the call tree (even inside functions defined elsewhere) sees the block-local value, which is the defining property of dynamic scope. `Return[v]`/`Return[v, Block]` is trapped via `eval_classify_return`. On exit the handler frees the OwnValues created during the block and restores each symbol's saved `own_values` and `attributes`, so the global binding is left exactly as it was — even if the body threw or returned early.

**Data structures.** A `SavedVar[]` array (name, saved `Rule*` own-value list, saved attribute bitmask) parallel to the variable list.
