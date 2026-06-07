---
source: src/print.c
---
`InputForm` is an unevaluated display wrapper: `builtin_inputform` (`src/print.c`) returns `NULL`, leaving `InputForm[expr]` intact. The printer's `print_standard` detects the `InputForm` head and renders the argument in a re-parseable form (a printer flag toggles InputForm-specific formatting); `ToString[expr, InputForm]` routes through the same standard printer.
