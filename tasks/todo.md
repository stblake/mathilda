# Task: RootReduce should thread over Rule (Solve output)

## Bug
`RootReduce` applied to a `Solve` result `{{u -> val, ...}}` emits
`RootReduce::argx: called with 0 arguments` and leaves `RootReduce[u -> val]`
unevaluated. `Reduce` output (`u == Root[...] && ...`) reduces fine because
`RootReduce` threads over `Equal`/`And`.

## Root cause
`src/rootreduce.c:is_option_rule` classifies *any* `Rule` with a symbol LHS as a
trailing option. So `u -> val` (u a symbol) is eaten as an option → `npos == 0`
→ `argx`. It should only treat `Method -> ...` (a registered option name) as an
option.

## Plan
- [x] Reproduce the bug (`RootReduce[u -> 3]` → argx)
- [x] Fix `is_option_rule`: only Rules whose LHS is a *registered option name*
      of RootReduce (i.e. `Method`) are options. Consult `symtab_get_options`.
- [x] Thread over `Rule` (not `RuleDelayed`): map RootReduce over the parts,
      rebuild, evaluate once (extracted shared `thread_parts` helper).
- [x] Update file header comment + docstring.
- [x] Update `docs/spec/builtins/algebra.md` (threads over rules too).
- [x] Changelog note in `docs/spec/changelog/2026-08-24.md`.
- [x] Add tests to `tests/test_rootreduce.c`; build & run.
- [x] Verify Method option still works; verify no argx; verify Solve example.

## Review

Root cause: `is_option_rule` classified *any* symbol-LHS `Rule` as a trailing
option, so a Solve entry `u -> value` was eaten as an option → `npos == 0` →
`RootReduce::argx`. Fix (all in `src/rootreduce.c`):

1. `is_option_rule` now recognises only a `Rule`/`RuleDelayed` whose LHS names a
   *registered* option of RootReduce — `Method`, read from `Options[RootReduce]`
   via `symtab_get_options` (new helper `is_rootreduce_option_name`). Any other
   symbol LHS is a positional argument.
2. Extracted the generic "map RootReduce over the parts, rebuild, evaluate once"
   tail of `thread_relational` into a shared `thread_parts` helper, and reused it
   to thread over an immediate `Rule`. `RuleDelayed` is *not* routed through
   `thread_parts` (it holds its RHS → would strand a `RootReduce`); its RHS is
   still reduced cleanly by the existing `rr_thread_coeffs` coefficient path.

Result: `Solve[...] // RootReduce` now returns the *same* `Root` objects as
`Reduce[...] // RootReduce`. `Method`, `argx`, `mtd` diagnostics unchanged.
All `rootreduce_tests` pass (incl. new `test_rootreduce_rule`).
