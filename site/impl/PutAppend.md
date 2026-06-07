---
source: src/readwrite.c
---
`builtin_putappend` shares `put_common` with `Put`, differing only in fopen mode: it opens the (last-argument) filename in `"a"` mode, appending each preceding expression — serialized via `expr_to_string` (the standard printer) plus a trailing `\n` — to the end of the file, creating it if absent and never truncating an existing one. `PutAppend["file"]` writes nothing. Open failure prints `PutAppend::noopen` and returns `$Failed`; success returns `Null`. `ATTR_PROTECTED`. The infix `expr >>> "file"` lowers to `PutAppend[expr, "file"]`.
