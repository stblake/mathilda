---
source: src/core.c
---
`builtin_symbol` (`src/core.c`) converts a string to a symbol. It validates the name with `symbol_name_is_valid` (each backtick-delimited context segment must start with a letter or `$` and continue with alphanumerics/`$`), runs it through `context_resolve_name` to apply the current context, and returns an `EXPR_SYMBOL`. Invalid names emit a `Symbol::symname` message and return `NULL`.
