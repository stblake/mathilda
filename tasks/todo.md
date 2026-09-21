# Task: Implement AlgebraicNumberDenominator

Plan: /Users/user/.claude/plans/continuing-on-from-our-crispy-stallman.md

Definition: AlgebraicNumberDenominator[a] = smallest positive integer n with n*a
an algebraic integer. NOTE: not qqbar_denominator (that is the min-poly leading
coeff, which over-counts, e.g. AlgebraicNumber[Sqrt[2],{1/5,1}] -> 5 not 25).
Algorithm: per-prime valuation over the primitive integer minimal polynomial.

## Phase A — Engine (src/poly/flint_qqbar.c + .h)
- [ ] Add `#include <flint/fmpz_factor.h>`
- [ ] `flint_qqbar_algebraic_number_denominator(const Expr* x, Expr** out)` engine fn
- [ ] FLINT-off stub (return -1)
- [ ] Declare in flint_qqbar.h

## Phase B — WL surface files
- [ ] src/poly/algebraicnumberdenominator.{c,h}  (::nalg message; Listable, Protected)

## Phase C — Registration & wiring
- [ ] sym_names.h/.c: SYM_AlgebraicNumberDenominator
- [ ] core.c: forward-declare + call init
- [ ] info.c: docstring
- [ ] version.h: bump 0.157 -> 0.158

## Phase D — Docs
- [ ] docs/spec/builtins/algebra.md: section
- [ ] docs/spec/changelog/2026-09-21.md: entry

## Phase E — Tests
- [ ] tests/test_algebraicnumberdenominator.c
- [ ] tests/CMakeLists.txt: COMMON_SRC + add_executable/add_test

## Phase F — Verify
- [x] make -j; make check-c99
- [x] run unit test
- [x] REPL smoke test all prompt examples
- [x] valgrind

## Review
All phases complete (v0.158). AlgebraicNumberDenominator[a] = smallest n with n*a
an algebraic integer, via per-prime valuation over the minimal polynomial (NOT
qqbar_denominator, which over-counts: AlgebraicNumber[Sqrt[2],{1/5,1}] -> 5 not 25).

Verification results:
- Build: clean (gcc-16, -std=c99 -Wall -Wextra), binary relinked.
- make check-c99: exit 0.
- Unit test (tests/test_algebraicnumberdenominator.c, ~40 cases): all passed.
- REPL smoke test: every example from the prompt matches exactly, incl. the
  AlgebraicNumber divergence case (5) and the ToNumberField round-trip
  (alpha = AlgebraicNumber[Sqrt[3], {-1, 1}], AlgebraicIntegerQ[alpha] -> True).
- valgrind: byte-identical to a no-builtin baseline (13,440 def-lost in both =
  known macOS baseline noise); zero leak stacks reference the new code.
- make check-packed-aware: OK. Sibling algebraic tests: no regression.

Files: src/poly/flint_qqbar.{c,h} (engine + stub), src/poly/algebraicnumberdenominator.{c,h}
(surface), src/sym_names.{c,h}, src/core.c, src/info.c, src/version.h (0.157->0.158),
docs/spec/builtins/algebra.md, docs/spec/changelog/2026-09-21.md,
tools/nd_fastpath_sweep.py (SKIP_EXPLOSIVE), tests/test_algebraicnumberdenominator.c,
tests/CMakeLists.txt.

Not done (awaiting user): git commit / tag v0.158 / push (release policy).
