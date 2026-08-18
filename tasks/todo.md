# Task: Diophantine solving in `Solve[..., Integers]`

Plan: `/Users/user/.claude/plans/we-do-not-have-tidy-sparrow.md`

## Phase 1 — core bounded engine (headline deliverable) — DONE

- [x] `src/solveint.h` — module contract
- [x] `src/solveint.c` — Stage A / B (bound propagation incl. lower-bound sign
      deduction) / C (recursive elimination + exact leaf) + verify + emit
- [x] Reuse `expr_to_mpoly`/`mpoly_*`, GMP integer roots, `internal_together`
- [x] Wire into `builtin_solve` as an `Integers` pre-pass
- [x] Integer-restriction filter in `solve_finish`
- [x] `solveint_init()` from `core_init()`
- [x] **Meet-in-the-middle** for separable additive equations (Euler 5s -> 0.8s)
- [x] `tests/test_solve_integers.c` + registered in CMake (8 tests pass)
- [x] Docs + changelog updated
- [x] `make check-c99` clean; solve suites pass; valgrind = baseline noise only

Solved: orig, 1, 3, 8, 10, 12, 14 (7 of the 16 incl. the motivating case).
Deferred cases (2, 4, 5, 6, 7, 9, 11, 13, 15) correctly return unevaluated.

## Phase 2 — divisor-factoring / reciprocal special forms — DONE
- [x] Bilinear divisor solver `(a u + c)(a v + b) = bc - ad` (Pythagorean),
      reached by linear elimination trying each keep-pair.
- [x] Unit-fraction recursion `Σ 1/x_i == R` with ordering (Egyptian).
- [x] Wired as `si_try_special_forms` before the unbounded decline.
- [x] Tests added; solve suites pass; check-c99 clean; valgrind = baseline.
  Pythagorean (z>0) -> 3 triangles; Egyptian 4/2027 -> 73 decompositions.

## Phase 2b (future) — MITM for n=3 additive (ex 11), linear-Diophantine
##   parametric/lattice output (ex 15, x+y==10)
## Phase 3 (future) — Pell via continued fractions (ex 4)
## Phase 4 (future) — Thue/Mordell/hyperelliptic/Catalan (ex 5,6,7,13)
