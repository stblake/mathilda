---
source: src/picostrings.c
---
`builtin_stringlength` checks for a single `EXPR_STRING` argument and returns `expr_new_integer((int64_t)strlen(arg->data.string))` — a byte count, not a codepoint count. Non-string arguments leave the call unevaluated (`NULL`). `ATTR_LISTABLE | ATTR_PROTECTED`, so it threads element-wise over a list of strings automatically.
