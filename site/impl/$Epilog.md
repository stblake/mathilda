---
source: src/repl_hooks.c
---
A REPL session hook, not a builtin. Registered (docstring only) in `repl_hooks_init` (`src/repl_hooks.c`). Unlike the per-line hooks, `$Epilog` is a bare symbol evaluated once at session teardown: `repl.c` calls `repl_apply_epilog()` on `Quit[]` or EOF, which — if an OwnValue is set — evaluates the symbol `$Epilog` (not a call) for its side effects and discards the result. Unset = no-op.
