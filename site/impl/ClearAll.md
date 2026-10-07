---
source: src/core.c
---
**Algorithm.** `builtin_clear_all` (`src/core.c`) is one of the four
symbol-management heads driven by `core_apply_symbol_action`, which walks the
argument list and applies a per-symbol action to each spec — a bare symbol, a
string naming a symbol, or a flat `List` of such specs (so `ClearAll[{a, b}]`
works), with the name read out by `core_symbol_name_of`. The action here is
`core_clear_all_one`.

For each name, `core_clear_all_one` first skips any `Protected`
symbol (which is what shields every builtin), then does the full erase that
distinguishes `ClearAll` from `Clear`: `symtab_clear_symbol` drops the
OwnValues/DownValues, the attribute word is zeroed (bumping the rule epoch so
the evaluation cache is invalidated), and the docstring (usage message) is
freed. The C builtin function pointer, if any, is left intact — the `Protected`
guard already keeps `ClearAll` away from builtins. The head returns `Null`.

**Attributes & limits.** `ClearAll` carries `HoldAll | Protected`, so the
symbols arrive unevaluated rather than being replaced by their current values.
Non-symbol/non-string specs are silently ignored.
