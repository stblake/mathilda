---
source: src/core.c
---
`builtin_information` (`src/core.c`) looks up the symbol's docstring with `symtab_get_docstring` and returns it as a string. If none exists it returns a string `No information available for symbol "..."` using `context_display_name` for the shortened name. (The interactive `?name` syntax routes to the same docstring store.)
