---
source: src/context.c
---
`builtin_context` (`src/context.c`) has three forms. `Context[]` returns `context_current()` (the `g_current` string, default `"Global`"`). `Context[sym]` / `Context["name"]` (the symbol form is held — `ATTR_HOLDFIRST`) report the symbol's context: if the name already carries a backtick prefix it is peeled off and returned via `context_prefix_len`; otherwise the symbol is looked up with `symtab_lookup` and reported as `"System`"` when it is a builtin (`def->builtin_func != NULL`), else `"Global`"`. An unknown string-form name emits `Context::notfound` and leaves the call unevaluated; an unassigned bare symbol defaults to `"Global`"`.
