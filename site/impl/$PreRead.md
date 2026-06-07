---
source: src/repl_hooks.c
---
A REPL session hook, not a builtin, and the only one operating on raw text rather than expressions. Registered (docstring only) in `repl_hooks_init` (`src/repl_hooks.c`). Before parsing, `repl.c` passes the input line through `repl_apply_pre_read`; if `$PreRead` has an OwnValue, the helper wraps the line as a `String`, evaluates `$PreRead[str]`, and expects a `String` back. A non-string return triggers a `$PreRead::strret` diagnostic and the original input is used. When unset (or NULL input) the original string is returned via a local `hooks_strdup`.
