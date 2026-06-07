---
source: src/core.c
---
`builtin_toexpression` (`src/core.c`) feeds a string argument to `parse_expression` (the Pratt parser, `src/parse.c`) and returns the parsed tree for the evaluator to reduce. An optional second argument (`InputForm`/`FullForm`/`StandardForm`) is accepted but ignored, since the parser is form-agnostic; an optional third argument is a head `h` wrapped around the result (commonly `Hold`). A parse failure returns `$Failed`; non-string input returns `NULL`. The symbol is `Listable`.
