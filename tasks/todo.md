# Task: Implement FirstPosition

Delegating to Position (shares Position's codebase). Attributes: HoldRest, Protected.

## Plan checklist

- [x] `src/patterns.h` — declare `builtin_first_position`
- [x] `src/patterns.c` — implement `builtin_first_position` + register in `patterns_init`
- [x] `src/sym_names.h` — `extern SYM_FirstPosition`
- [x] `src/sym_names.c` — define + intern `SYM_FirstPosition`
- [x] `src/info.c` — docstring for FirstPosition
- [x] `src/options_builtin.c` — `{ "FirstPosition", "True" }` Heads default
- [x] `tests/test_firstposition.c` — new test file (23 cases)
- [x] `tests/CMakeLists.txt` — register `firstposition_tests`
- [x] `docs/spec/builtins/pattern-matching.md` + `data-structures.md` — FirstPosition entries (Position's canonical doc is in pattern-matching.md, not structural-manipulation.md)
- [x] `docs/spec/changelog/2026-09-14.md` — changelog note
- [x] Build clean (`make -j`, `make check-c99`)
- [x] Run `firstposition_tests` (23/23) + `patterns_tests` (no regression)
- [x] REPL spot-check all 15 spec examples — all match
- [x] leak check (`leaks --atExit`): 0 leaks
- [x] Rebuild code-review graph

## Review

Implemented `FirstPosition[expr, pattern, default, levelspec]` (attributes
`HoldRest, Protected`) by **delegating to `Position`** with the first-match cap
(`Position[expr, pattern, levelspec, 1, Heads->opt]`), taking `[[1]]`, and falling
back to the held `default` (evaluated only when returned) or `Missing["NotFound"]`.
Associations without a levelspec use the 2-arg `Position` form so its value→`Key[...]`
remapping fires. This maximally shares Position's code (levelspec parsing, `Heads`,
traversal/ordering, association handling) and is Wolfram-faithful because Position's
ordering is already asserted by `test_position`.

Verified: all 15 documented examples reproduce the spec exactly; 23-case unit suite
passes; no `Position` regression; 0 memory leaks; clean build + C99 gate.

Design note: `FirstPosition` is a structural/pattern-search head (not numeric), so —
like `Position` — it carries no ND-kernel / packed-aware / `Compile[]` obligation and
reuses `patterns_delist_visible` for visible `NDArray` inputs. Known limitation
(shared with `Position`): a levelspec over an association is not `Key`-remapped.
