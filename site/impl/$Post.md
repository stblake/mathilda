---
source: src/repl_hooks.c
---
A REPL session hook, not a builtin. Registered (docstring only, no default OwnValue) in `repl_hooks_init` (`src/repl_hooks.c`). In each REPL cycle `repl.c` calls `repl_apply_post(result)` on the evaluator's output; when an OwnValue is set the helper builds `$Post[expr]` and runs it through `evaluate()` via `hook_call_eval`. This happens after evaluation and before `Out[n]` is stored, so `$Post` can transform the visible result. Unset = identity.
