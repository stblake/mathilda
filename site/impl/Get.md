---
source: src/readwrite.c
---
**Algorithm.** `builtin_get` reads a Mathilda source file and evaluates it expression by expression, returning the last value. It opens the file (`Get::noopen` + `$Failed` on failure), slurps the entire contents into a `malloc`'d buffer, then walks the buffer with the parser's `parse_next_expression(&ptr)` — the same Pratt parser used by the REPL — `evaluate`ing each parsed expression and keeping the last non-`NULL` result (defaulting to `Null` for an empty file). Parsing stops when `parse_next_expression` returns `NULL` at end-of-input. This is the mechanism `init.m` uses to load the internal `.m` bootstrap files. `ATTR_PROTECTED`.
