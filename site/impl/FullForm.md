---
source: src/print.c
---
`FullForm` is an unevaluated display wrapper: the builtin `builtin_fullform` (`src/print.c`) returns `NULL`, leaving `FullForm[expr]` intact. Rendering is done by the printer — `print_standard` detects the `FullForm` head and calls `expr_print_fullform`, which writes the raw tree as `head[arg, ...]` with no infix/operator sugar. `ToString[expr, FullForm]` reuses the same path via `expr_to_string_fullform`.
