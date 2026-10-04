---
source: src/modular.c
---
**Algorithm.** `builtin_unique` (`src/modular.c`) generates fresh, never-before-used
symbols. `Unique[]` uses the prefix `"$"`; `Unique[x]` or `Unique["x"]` uses the
symbol's name or the string as prefix (`unique_prefix_of`); `Unique[{a, b, ...}]`
returns a list of fresh symbols, one per element, all sharing one numeric suffix
(the Wolfram behaviour). The suffix is drawn from the module counter
`module_number`, the same monotone source `Module` uses for its `x$n` temporaries.

`unique_make_batch` guarantees freshness by *scanning*: it advances
`module_number` until `prefix<n>` names no existing symbol for every prefix in
the batch, then takes that value and post-increments. Each generated symbol is
created with `expr_new_symbol` and tagged `ATTR_TEMPORARY`, and the user-visible
`$ModuleNumber` OwnValue is re-synced from the counter (`unique_sync_module_number`),
so `Unique` and `Module` share one numbering. An unusable prefix (not a symbol or
string) makes the batch fail and the head returns `NULL`.

**Attributes & limits.** `Unique` has no special attributes — its single argument
is evaluated normally. It accepts zero or one argument; the list form yields one
suffix-sharing symbol per element.
