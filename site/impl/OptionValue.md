---
source: src/options_builtin.c
---
**Algorithm.** `builtin_optionvalue` (`src/options_builtin.c`) resolves a single
option setting. It accepts `OptionValue[f, name]`, `OptionValue[f, opts, name]`,
and a four-argument `... , Hold` form; the bare `OptionValue[name]` has no
context at top level and is left unevaluated (it is only meaningful inside a
fired rule, where `optionvalue_inject_context` rewrites it to the three-argument
form carrying the enclosing head and its explicit options). The name must be a
symbol or string.

Resolution (`ov_resolve`) checks **explicit options first**, then defaults
derived from `f`: `lookup_value` scans the explicit `opts` (a single rule or a
list), and `defaults_lookup` interprets `f` as `Automatic` (no defaults), a
symbol (→ its `Options`), a single rule, or a list of such specs, first match
winning. All name comparisons strip any context prefix (`strip_context`) so a
symbol and its qualified string name compare equal. An unresolved option returns
`NULL` (the call stays unevaluated); the `Hold` form wraps the resolved value in
`Hold[...]`.

**Data structures & limits.** Values are read from the same
`List[Rule[name, val]]` option representation used by `Options`, and every
returned value is `expr_copy`'d. `OptionValue` is `Protected`.
