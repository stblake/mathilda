---
source: src/repl_hooks.c
---
A REPL session hook, not a builtin. Registered (docstring only) in `repl_hooks_init` (`src/repl_hooks.c`). `repl.c` calls `repl_apply_pre_print(out)` just before printing; if an OwnValue is set, `hook_call_eval` builds and evaluates `$PrePrint[expr]`. Crucially this is display-only: `Out[n]` is assigned the unmodified post-`$Post` result above, and only the rendered form reflects the `$PrePrint` value. Unset = identity.
