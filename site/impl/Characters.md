---
source: src/picostrings.c
---
`builtin_characters` requires a single `EXPR_STRING` argument (else returns `NULL`). It allocates an `Expr*` array of `strlen(str)` entries and, for each byte, builds a length-1 `EXPR_STRING` via `expr_new_string` from a two-char `buf`, wrapping the lot in a `List`. It is byte-oriented (one element per `char`), not Unicode-codepoint aware. Carries `ATTR_LISTABLE | ATTR_PROTECTED`.
