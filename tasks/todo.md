# Task: Move Solve family into `src/solve/` and split `solveint.c`

Plan: `/Users/user/.claude/plans/the-source-code-for-temporal-haven.md`
Pure structural refactor — no behaviour change.

## Phase 1 — Baseline ✅
- [x] Build `Mathilda` on current tree (clean; already up to date)
- [x] Build + run 9 solve-family suites — ALL PASS (baseline_tests.txt)
- [x] Capture REPL corpus (18 clusters, 1.4s) → baseline_corpus.out

## Phase 2 — Move only (no split) ✅
- [x] `git mv` all 10 `solve*.{c,h}` into `src/solve/` (git renames)
- [x] makefile: added `$(wildcard $(SRC_DIR)/solve/*.c)` + `-I./src/solve`
- [x] tests/CMakeLists.txt: repointed 10 paths; added `include_directories(../src/solve)`
- [x] `make clean && make -j` clean (1m14s, no diagnostics); tests rebuilt
- [x] Corpus byte-identical; all 9 suites pass identical to baseline

## Phase 3 — Split solveint.c ✅
- [x] Created `solveint_internal.h` (5 macros, 3 structs, 11 static-inline helpers, 49 protos)
- [x] Created `solve_common.c` (32 shared funcs) + 17 `solveint_*.c` method files
- [x] `solveint.c` reduced 5667→461 lines (dispatcher + special-forms + fermat + init)
- [x] Added 17 new `.c` files to tests/CMakeLists.txt COMMON_SRC
- [x] Clean build; **pure-move verified** (all 136 bodies byte-identical); corpus + tests identical

## Phase 4 — Docs/comments ✅
- [x] Updated `src/linalg/hnf.c` + `hnf.h` + `test_solve_integers.c` comments
- [x] SPEC.md §2 layout tree: added `src/solve/`
- [x] Changelog note in `docs/spec/changelog/2026-08-17.md`
- [x] Path refs in SOLVE_INTEGERS.md, spec/builtins/*, benchmark README/REPORT/run.py

## Phase 5 — Verify ✅
- [x] `make check-c99` clean; full clean rebuild from scratch clean (1m16s)
- [x] All 9 solve suites green, identical to baseline
- [x] Behaviour corpus byte-identical (baseline vs move vs split vs final)
- [x] valgrind: def-lost 13,632→13,480 (pre vs post split, identical within noise) — no new leaks
- [x] Rebuild code-review graph

## Phase 6 — (user request) fix context.c:39 leak ✅
- [x] Root cause: `context_shutdown()` declared + written but **never called** (dead code);
      context state (g_current / $ContextPath / package-frame save-slots via ctx_strdup)
      lingers to exit. Frame accounting itself is correct (instrumented: 2 push / 2 pop, balanced).
- [x] Fix: call `context_shutdown()` on all post-init exit paths in `repl.c main()` (+`#include "context.h"`).
      context.c itself untouched.
- [x] Verified: 9 context blocks now freed, 0 context.c blocks lost; $Context/$ContextPath
      + BeginPackage/Begin/End/EndPackage work; corpus identical.

## Review
**Outcome:** Pure structural refactor delivered. `src/solve/` now holds the whole Solve
family; `solveint.c` split 5667→461 lines across `solve_common.c` + 16 `solveint_<method>.c`
+ `solveint_internal.h`. Proven a pure move: all 136 function bodies byte-identical. Build
(make + CMake) rewired; docs updated. Plus a bonus fix: wired up the never-called
`context_shutdown()` surfaced by the valgrind pass. Not yet committed (awaiting user).
