---
source: src/repl_hooks.c
---
A REPL session hook, not a builtin. `repl_hooks_init` (`src/repl_hooks.c`) merely touches the symbol so `?$Pre` works; no default OwnValue is installed, so out of the box the hook is a no-op. Each REPL iteration `repl.c` calls `repl_apply_pre(parsed)`; if `hook_is_set("$Pre")` (i.e. `symtab_get_own_values("$Pre")` is non-empty) it builds `$Pre[expr]` and runs it through the standard `evaluate()` via `hook_call_eval`, applied after parsing but before the main evaluation. Because the wrapped expression is evaluated before `$Pre` sees it unless `$Pre` is assigned a `HoldAll` head, its effect is usually indistinguishable from `$Post`.
