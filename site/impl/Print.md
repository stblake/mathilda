---
source: src/print.c
---
`builtin_print` (`src/print.c`) calls `print_standard` on each argument in turn (no separators), emits a trailing newline, and returns the symbol `Null`.
