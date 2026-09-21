# Task: Implement NumberFieldIntegralBasis (+ AlgebraicIntegerQ)

Plan: /Users/user/.claude/plans/continuing-on-from-our-magical-alpaca.md

## Phase A — Engine (src/poly/flint_qqbar.c + .h)
- [ ] Add `#include "numberfield.h"` + `#include "numberfield_internal.h"`
- [ ] `flint_qqbar_integral_basis(const Expr* a)` engine fn
- [ ] `flint_qqbar_algebraic_integer_q(const Expr* x)` engine fn
- [ ] FLINT-off stubs for both
- [ ] Declare both in flint_qqbar.h

## Phase B — WL surface files
- [ ] src/poly/numberfieldintegralbasis.{c,h}
- [ ] src/poly/algebraicintegerq.{c,h}

## Phase C — Registration & wiring
- [ ] sym_names.h/.c: SYM_NumberFieldIntegralBasis, SYM_AlgebraicIntegerQ
- [ ] core.c: forward-declare + call both inits
- [ ] info.c: docstrings for both
- [ ] version.h: bump + tag

## Phase D — Docs
- [ ] docs/spec/builtins/algebra.md: two sections
- [ ] docs/spec/changelog/2026-09-21.md: entry
- [ ] tools/nd_fastpath_sweep.py: SKIP_EXPLOSIVE += NumberFieldIntegralBasis

## Phase E — Tests
- [ ] tests/test_numberfieldintegralbasis.c
- [ ] tests/CMakeLists.txt: COMMON_SRC + add_executable/add_test
- [ ] optional leak script

## Verification
- [ ] make -j; make check-c99
- [ ] build + run tests
- [ ] REPL smoke test of all examples
- [ ] valgrind vs baseline
- [ ] audits green
- [ ] graph refresh + memory note

## Review

Status: DONE (pending user commit).

- Engine `flint_qqbar_integral_basis` + `flint_qqbar_algebraic_integer_q` in
  `src/poly/flint_qqbar.c` (reuse `to_qqbar`/`algint_generator`/`poly_to_algnum`
  + `nf_field_create`/`nf_ok_basis`). Round-2 HNF canonicalised to lower-triangular
  (1-first, ascending degree) via rot180 of the column-reversed `fmpz_mat_hnf`.
- Thin WL wrappers `numberfieldintegralbasis.c` (Listable|Protected) and
  `algebraicintegerq.c` (Protected). Registered in core.c; SYM_ in sym_names;
  docstrings in info.c; version 0.157.
- Docs: algebra.md two sections; changelog 2026-09-21.md; SKIP_EXPLOSIVE.
- Tests: test_numberfieldintegralbasis.c — ALL PASS. Sibling suites
  (numberfield/algebraicnumber/algebraicnumberpolynomial) still pass.
- Verification: build clean; check-c99 clean; check-packed-aware OK; valgrind
  definitely-lost byte-identical to baseline (no new leaks).
  check-compile-coverage RED is pre-existing (Image*/PackedArrayQ), my heads not
  implicated.
- Every user example reproduced exactly (via N[]); AlgebraicIntegerQ of every
  basis element and integer combination is True.

TODO for user: `git add` + commit (v0.157) + `git tag v0.157`.
