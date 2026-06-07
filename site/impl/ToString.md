---
source: src/core.c
---
`builtin_tostring` (`src/core.c`) renders an expression to a string. The optional second argument selects the form: `FullForm` uses `expr_to_string_fullform`; `TeXForm` wraps in `TeXForm[...]` and prints; `InputForm`/`StandardForm`/`OutputForm` (and the default) use the standard printer `expr_to_string`. All formatting is shared with the `src/print.c` printer.
