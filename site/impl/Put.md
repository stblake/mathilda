---
source: src/readwrite.c
---
`builtin_put` (and `builtin_putappend`) share `put_common`, which treats the **last** argument as the filename and writes every preceding expression to it, one per line. `Put` opens the file in `"w"` mode (truncating); each `expr_i` is serialized with `expr_to_string` (the standard printer) followed by `\n`. The single-argument form `Put["file"]` runs the loop zero times, leaving an empty file. Open failure prints `Put::noopen` and returns `$Failed`; success returns the symbol `Null`. The output is intended to be re-readable with `Get`. `ATTR_PROTECTED`. The infix `expr >> "file"` lowers to `Put[expr, "file"]`.
