# Task: Implement AlgebraicNumber and ToNumberField

Faithful WL recreation. FLINT `qqbar` backend (declines cleanly without FLINT).
Plan: /Users/user/.claude/plans/let-s-implement-tonumberfield-tonumberfi-sorted-rain.md

## Phase 0 — Skeleton & registration
- [ ] SYM_AlgebraicNumber, SYM_ToNumberField (3-site in sym_names.{h,c})
- [ ] src/poly/algebraicnumber.{c,h} — builtin stub returning NULL + _init
- [ ] src/poly/tonumberfield.{c,h} — builtin stub returning NULL + _init
- [ ] Wire _init() in core.c (near minpoly/rootreduce block)
- [ ] Attributes: AlgebraicNumber NHOLDALL|PROTECTED, ToNumberField PROTECTED
- [ ] Docstrings in src/info.c (no examples)
- [ ] tests/CMakeLists.txt: add both .c to COMMON_SRC + test target
- [ ] Build clean, parses/prints, no-op

## Phase 1 — Recognition in flint_qqbar.c
- [ ] to_qqbar + is_constant_algebraic + collect_atoms handle AlgebraicNumber
- [ ] Verify RootReduce[AN]→Root; value-preservation

## Phase 2 — AlgebraicNumber canonicalisation
- [ ] lc/φ/M reduction, rational branch, collapse, over-length fold, empty→0, malformed→NULL
- [ ] expr_eq fixpoint guard
- [ ] Verify examples 1,4–12 + idempotence

## Phase 3 — numericalize branch (numeric.c)
- [ ] Build Σ ci θ^i, recurse. Verify N[...,50], Round, Less

## Phase 4 — is_numeric_quantity (core.c) + Re/Im (complex.c)
- [ ] Verify Re[real]→self, Im→0, complex Equal/Unequal

## Phase 5 — Field arithmetic (plus/times/power)
- [ ] algnum combine pre-pass. Verify examples 2,13–16 + mismatched-gen sum
- [ ] Regression: plus/times/power/rootreduce suites

## Phase 6 — ToNumberField + docs
- [ ] All arg forms, decline a∉Q(θ), ToNumberField[2,1/2]→2
- [ ] docs/spec + changelog + full-suite run + valgrind

## Review — COMPLETE

All phases done and verified. Every WL spec example reproduced exactly.

**Files changed**
- `src/sym_names.{h,c}` — SYM_AlgebraicNumber, SYM_ToNumberField (3-site).
- `src/poly/flint_qqbar.{c,h}` — to_qqbar/is_constant_algebraic/collect_atoms
  recognise AlgebraicNumber, Re/Im/Abs/Conjugate and E^(I Pi r); Expr-level entry
  points (canonicalise, to-number-field, common-field, field arithmetic).
- `src/poly/algebraicnumber.{c,h}`, `src/poly/tonumberfield.{c,h}` — new builtins.
- `src/numeric.c` — N branch for AlgebraicNumber. `src/core.c` — is_numeric_quantity
  clause + init wiring. `src/complex.c` — Re/Im/Abs real-algebraic hook.
- `src/plus.c`, `src/times.c`, `src/power.c` — field-arithmetic combine passes
  (placed after the int64 fast path — no hot-path cost).
- `src/info.c` — docstrings. Docs: `docs/spec/builtins/algebra.md` + changelog.
- Tests: `tests/test_algebraicnumber.c` (+ CMake target, COMMON_SRC entries).

**Verification**
- `tests/test_algebraicnumber.c`: all pass (spec examples pinned + RootReduce
  value-preservation cross-checks all == 0).
- Regressions green: core, eval, evaluate, numeric, numeric_stress/domain/complex,
  complexexpand, power*, root_numeric, core_algebra, simplify, expand, comparisons,
  boolean, rootreduce, numberfield.
- `make check-c99`: clean. `make check-packed-aware`: my heads not flagged
  (pre-existing FirstPosition failure is unrelated, from commit 34bbbf37).
- Valgrind: no new lost blocks vs macOS baseline (definitely/indirectly/possibly
  lost byte-identical to the `1+1` baseline run).

**Known limitations (documented, not bugs)**
- `RootReduce[Re[complex algebraic]]` yields a valid but differently-normalised
  representation than WL (same value — cross-checked == 0). WL documents that
  AlgebraicNumber representations are non-unique.
- `Automatic` and `All` both use the qqbar minimal primitive element.
