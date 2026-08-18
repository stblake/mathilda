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

## Phase 3 — Pell via continued fractions — DONE
- [x] si_pell_detect (x^2 - D y^2 == +/-1), si_pell_cf (CF of sqrt(D)),
      orbit generation up to bound, 4 sign variants + verify.
- [x] Negative Pell handled (base iff period odd); unsolvable -> {}.
- [x] Tests added; solve/poly suites pass; check-c99 clean; valgrind = baseline.

## Perf — special forms tried first + numeric verification — DONE
- [x] Pythagorean z>0: 7.45s -> ~0.5ms (was routed to leaf search, not divisor).
- [x] Fixed latent mpoly_to_expr array leak.

## Phase 2b — linear Diophantine parametric + bounded gcd test — DONE
- [x] Parametric family for unconstrained multivar linear (gcd staircase).
- [x] Bounded: gcd unsolvability -> {} (ex 16).

## Phase 2c — LLL lattice enumeration for solvable bounded linear — DONE
- [x] LatticeReduce-reduced basis; coefficient box = value box projected
      through the pseudoinverse (B B^T)^{-1} B; exact per-candidate check.
- [x] Few-solution large-coeff boxes enumerated (2000-pt case < 1s);
      dense/intractable boxes declined. Verified complete (gap-free 2-var
      progression), no dups, leak-clean.

## Phase 2d — odd-power-sum divisor method — DONE
- [x] Divisor method for separable sums of odd powers: s=x+y | m, power sum
      p_e(s,p) via Newton recurrence -> integer roots p -> (x,y). Outer vars
      enumerated => O(N*factoring), not O(N^2).
- [x] Sum of three cubes ==42 over |.|<10^5 in ~7s; ==3 finds all knowns;
      two-cube taxicab; general over odd e (5th, 7th). m-size guard declines
      when factoring would be impractical (5th powers over |.|<10^5).
- [x] Verified complete; tests added; suites pass; check-c99; valgrind baseline.

## Phase 4 (future) — Thue/Mordell/hyperelliptic/Catalan/general-N Pell;
##   Booker-style O(N) sum-of-three-cubes; MITM sieving refinements
