---
source: src/core.c
---
`builtin_clear` (`src/core.c`) iterates its arguments and, for each that is a symbol, calls `symtab_clear_symbol(name)` to remove that symbol's OwnValues and DownValues (its rules/assignments) while leaving the symbol itself, its attributes, and any builtin binding intact. Non-symbol arguments are ignored. Returns `Null`. It carries `ATTR_HOLDALL` so the symbols are not evaluated to their current values before being cleared.
