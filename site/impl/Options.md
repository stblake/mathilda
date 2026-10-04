---
source: src/options_builtin.c
---
**Algorithm.** `builtin_options` (`src/options_builtin.c`) returns an object's
option list, always as a freshly built `List` (never `NULL` for a recognised
shape). `options_of_object` branches on the first argument: a **symbol** yields a
copy of its registered defaults via `symtab_get_options` (or `{}` if it has
none); a **compound expression** yields just the option rules explicitly present
among its own arguments (each `Rule`/`RuleDelayed` whose left side is a
symbol/string, detected by `is_option_rule`); anything else yields `{}`.

The two-argument forms `Options[obj, name]` and `Options[obj, {names}]` build the
full list and then select the matching whole rules with `lookup_rule`, whose name
comparison is context-insensitive (`strip_context` drops any `` ` ``-qualified
prefix). Selected rules are returned in the requested order; a name with no
setting simply contributes nothing.

**Data structures & limits.** A symbol's defaults are stored on
`SymbolDef.default_options` as `List[Rule[name, val], ...]` — the
`DefaultValues`-equivalent — reached through `symtab_get_options`/`_set_options`.
Every element handed back is `expr_copy`'d, so the stored list is never aliased
or mutated. `Options` is `Protected`.
