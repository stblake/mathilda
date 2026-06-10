# Task: Implement ClearAll, Remove, Protect, Unprotect

## Plan
- [ ] sym_names.{h,c}: declare + intern SYM_ClearAll, SYM_Remove, SYM_Protect, SYM_Unprotect
- [ ] symtab.{h,c}: add `symtab_remove_symbol(name)` (full table-entry deletion)
- [ ] core.c: implement builtin_clear_all, builtin_remove, builtin_protect, builtin_unprotect
  - HoldAll; respect Locked/Protected guards per WL semantics
  - Protect/Unprotect return list of changed symbol-name strings (WL-faithful)
  - ClearAll/Remove return Null
  - Support symbol, string, and List[...] arg forms
- [ ] core.c core_init(): register builtins + attributes
  - ClearAll: HoldAll | Protected
  - Remove:   HoldAll | Locked | Protected
  - Protect:  HoldAll | Protected
  - Unprotect:HoldAll | Protected
- [ ] core.h: prototypes
- [ ] info.c: terse docstrings for all four
- [ ] tests/: extensive unit tests + CMake wiring
- [ ] docs/spec/builtins + changelog update
- [ ] Build clean (c99 -Wall -Wextra), run tests, valgrind no new leaks

## Notes
- Protection enforced via get_attributes() (static floor | def->attributes); Unprotect
  manipulates def->attributes, fully effective for user symbols (consistent w/ ClearAttributes).
- Remove guarded against Protected/Locked so builtins can't be deleted.

## Review

Done. All four builtins implemented and verified.

- **sym_names**: added SYM_ClearAll, SYM_Remove, SYM_Protect, SYM_Unprotect (decl + intern).
- **symtab**: new `symtab_remove_symbol()` — unlinks + frees the SymbolDef (rules,
  docstring), leaves the interned name (owned by sym_intern), bumps eval clock.
- **core.c**: builtin_clear_all / _remove / _protect / _unprotect + shared helpers.
  - ClearAll/Remove return Null; Protect/Unprotect return List of changed names (strings) — WL-faithful.
  - Guards via get_attributes(): ClearAll/Remove skip Locked|Protected; Protect/Unprotect skip Locked.
  - Accept symbol, string, or flat List of specs.
- **Attributes**: ClearAll/Protect/Unprotect = HoldAll|Protected; Remove also Locked.
- **info.c**: terse docstrings for all four.
- **Tests**: tests/test_clearall_remove_protect.c, 26 cases, all pass (Debug build, rc=0).
  Only stderr = 2 expected `wrsym` lines from the protection-block tests.
- **valgrind**: 12,864B/402 blocks "definitely lost" == macOS baseline (chop_tests: 12,800/400);
  no Mathilda/src frames, none of my functions appear in any leak stack. No new leaks.
- **Regression**: symtab_tests pass. core_tests aborts at test_gcd_lcm — PRE-EXISTING test bug
  (expects standard-form string but passes is_fullform=1; REPL gives the right value). Unrelated.
- Full `make` clean (-std=c99 -Wall -Wextra), binary links.

Known limitation (consistent w/ existing ClearAttributes): Unprotect on a builtin listed in
attr.c's static `builtin_attrs` table clears the def bit but get_attributes() re-ORs the static
floor, so enforcement still treats it as Protected. Fully effective for user symbols (the WL
Unprotect[f];...;Protect[f] workflow). Not changed — out of scope, would touch get_attributes globally.
