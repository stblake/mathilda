---
source: src/options_builtin.c
---
**Algorithm.** `builtin_setoptions` (`src/options_builtin.c`) redefines
individual default options of a symbol. The first argument must be a symbol. It
takes a working
copy of the symbol's current option rules (from `symtab_get_options`) as a flat
vector and, for each trailing `name -> value` rule, finds the existing option of
that name (context-insensitive match) and **replaces it in place**, preserving
the option's original position.

`SetOptions` cannot *add* an option: a name not already in the symbol's defaults
raises `SetOptions::optnf` and the call aborts leaving the stored options
unchanged; a malformed (non-rule) argument returns `NULL`. On success the new
list is stored with `symtab_set_options`, the evaluation cache is dropped
(`eval_clock_bump`, since option changes can alter results), and a copy of the
updated option list is returned.

**Data structures & limits.** Options live on `SymbolDef.default_options` as
`List[Rule[...]]`; the working vector holds `expr_copy`'d rules so nothing
aliased is mutated before the atomic store. `SetOptions` is `Protected`.
