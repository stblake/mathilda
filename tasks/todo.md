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

## Phase 4 — exponential Diophantine + bounded elliptic coverage — DONE (partial)
- [x] Exponential Diophantine (variable exponents): bounded-box enumeration +
      Catalan/Mihailescu (x^a-y^b==+/-1 -> unique 3^2-2^3). si_solve_exponential
      runs before the MPoly stage. Handles the Catalan example.
- [x] Confirmed + tested bounded Mordell/hyperelliptic (y^m==f(x)) resolve via
      the leaf search: y^2=x^3-10000 -> {25,75}, y^2=x^3-2 -> {3,5}, etc.
- Still deferred (genuinely research-grade, unbounded): Mordell-Weil/elliptic
  integral points with no box (ex 5/7 unbounded), indefinite Thue equations
  (ex 6), general-N Pell, Booker O(N) sum-of-three-cubes.

## Benchmark vs sympy + correctness fix (2026-08-18)

**Benchmark** `benchmarks/87-diophantine-integers/` — head-to-head against
sympy `diophantine` across 19 well-known families. `cases.py` is the single
source of truth; `run.py` generates `diophantine.m`, runs it, runs sympy + a
same-box Python search, writes `REPORT.md` + `results.json`. Result: **Mathilda
19/19 exact** (all counts verified); sympy answers 3 directly/comparably, a few
more only parametrically or after manual elimination, and `NotImplementedError`
on 11 (all cubic/exponential forms + every system).

**Fix found while benchmarking** — `Solve[x^2+y^2==25,{x,y},Integers]` returned
`{}` (silent wrong answer). Root cause: the bounder only bounded provably-
non-negative variables, so an unconstrained sum of even powers was declined and
Solve fell through to the generic path → `{}`. Added `derive_even_only_bounds`
(sign-symmetric even-only vars get `[-B,B]`) + a one-line `urest==0` tightening
(fixes `x^2+y^2==0`). Verified vs sympy (12/28/72 exact). New test
`test_sum_of_two_squares`. check-c99 clean, all solve suites pass, valgrind =
baseline.

## Diophantine campaign — correctness + 10 new methods (2026-08-18) — DONE

Plan: `~/.claude/plans/cosmic-riding-hoare.md`. Driven by a ~18-case stress test.

- **P0 correctness** (`src/solve.c`): `y^2==x^3-2 → {}` was a silent wrong answer
  (unbounded nonlinear curve solved parametrically, then Integers-filtered to
  `{}`). `is_single_multivar_equation` + a `solve_finish` guard now leave such a
  lone multivar equation **unevaluated** after solveint declines. Also `y==x^2`,
  `y^2==x^3+1`, ... Regression bank added (curves, verified negatives, the
  already-passing `y^3==x^5-x+1` / `x^2+y^3==z^7` / nine-cubes==239).
- **Tranche A** (bounded engine): multi-leaf staged elimination (`si_solve_multileaf`
  → Euler brick), ordering-reduced estimate (`si_longest_chain`) + int64 fast leaf
  (`si_leaf_roots_i64` + coef cache → `2Σsq=(Σ)^2`, Markov-Hurwitz), non-poly
  power-leaf (`si_solve_bounded_powerleaf` → Brocard `n!+1==m^2`).
- **Tranche B/C/D** (closed form): conic `Y^2=AX^2+BX+C` (`si_solve_conic` →
  `n^2+n+41`), unbounded Pell parametric (`si_solve_pell_parametric`), homogeneous
  linear ray (`si_solve_linear_system_ray` + Bareiss det), fixed-base exponential
  (`3^m-2^n=1`), PTE→{} via Newton (`si_solve_power_sum_equal`), general unbounded
  imaginary Mordell (`si_solve_mordell`: k<0 squarefree, k≡2,3 mod4, 3∤h via
  reduced-form class number → `(3,±5)`, `(17,±70)` for k=-13, `{}` for k=-5;
  verified vs brute force for all 49 engaged k in [-150,-2]).
- **Honest declines** (unevaluated, never wrong): Elkies `x^4+y^4+z^4==w^4` @1e7,
  Cassels `3x^3+4y^3+5z^3`, Pyth-area (Fermat), Ramanujan-Nagell `x^2+7==2^n`
  (needs the descent in Z[(1+√-7)/2]).
- 11 new tests in `test_solve_integers.c`; all solve suites + corpus 97/97 pass;
  `make check-c99` clean; valgrind = startup baseline. Docs + changelog updated.

## Tier 1 primitives campaign (2026-08-18) — matching Mathematica's Diophantine dispatcher

Driven by a review of Mathematica's ~25 engines vs Mathilda's current coverage.
The bounded search + closed-form layer is strong (19/19 vs sympy); the gaps are
all in the *unbounded / algebraic-number-theory* half.

**Held-out measurement first (the benchmark over-fit).** `benchmarks/87` is the
single source of truth AND was co-designed with each method (Phase 2d added the
sum-of-three-cubes cases with the divisor method, Phase 3 the Pell cases with the
CF solver, ...), so 19/19 measures "method works on the example it was written
for", not coverage. A held-out corpus (15 equations from standard references,
none in `cases.py`) run cold found: the bounded engine is genuinely complete on
every boxed Thue/Mordell/Pell-N/Legendre case, but **two silent wrong answers**
in unbounded *linear systems*, plus honest declines on factorable-leading-form
BQF.

- [x] **P0(a) correctness** — underdetermined linear system over Integers →
      silent wrong `{}`. Fixed: `solvelinsys` declines the underdetermined
      integer case; `src/solve.c` `solution_set_is_parametric` guard converts any
      unproven parametric `{}` to unevaluated. Test added; 8 suites pass;
      check-c99 clean. (Changelog: "underdetermined linear SYSTEM ... {}".)
- [x] **P0(b) HNF** — DONE. `HermiteDecomposition` builtin (`src/linalg/hnf.c`,
      row HNF + unimodular transform, `u.m==r`) + integer linear-system solver
      (`si_solve_linear_system_hnf`): particular solution by forward substitution
      with divisibility test (failure ⇒ `{}` proof), kernel lattice as `C[k]`.
      `{x+2y+3z==10, x-y+z==2}` → `{{x->18+5C[1], y->8+2C[1], z->-8-3C[1]}}`;
      `2x+2y==3 && x-y==0` → `{}` (precise, supersedes P0a's unevaluated). Tests
      added (74 solve-int cases + HermiteDecomposition invariants); check-c99
      clean; valgrind delta = 0. Docs + changelog updated.
- [x] **P1 Runge / factorable binary quadratic** — DONE.
      `si_solve_factorable_conic` (`src/solveint.c`): a 2-var degree-2 equation
      with a cross term and perfect-square discriminant `δ=B²-4AC>0` completes to
      a difference of squares `(2kU)²-V²=W` and divisor-enumerates `W` —
      exhaustive, so `{}` is a proof. `x²+xy-2y²==4` → 6 points;
      `(x-y)(x+2y)==15` → `{}` (mod-3 obstruction, NOT a decline — the earlier
      "should have solutions" read was wrong); `2x²+3xy-2y²==7` → 2 points
      (non-unit leading coeffs the old conic couldn't do). Generalises
      `si_solve_conic`. Test `test_factorable_conic`; check-c99 clean; valgrind
      delta 0. Docs + changelog updated.
- [x] **P2 general-N Pell (unbounded)** — DONE. `si_solve_genpell_parametric`
      (`src/solveint.c`): unbounded `x²-Dy²==N && x>0 && y>0`, any `N≠+1`
      (incl. `N=-1`) → one `ConditionalExpression` family per solution class.
      Fundamental unit via CF of √D; Nagell bound `y≤u√(|N|/(2(t±1)))` bounds the
      fundamental search; each found (±x,±y) advanced into the positive orthant
      then reduced by ε⁻¹ to the minimal class rep (dedup); orbit
      `(a+b√D)(t+u√D)^C[1]`, C[1]≥0. Exhaustive → `{}` is a proof (incl. negative
      Pell over even-CF-period D). **Validated against brute force over ~30
      (D,N) pairs** (multi-class, negative Pell, D=61 large unit, unsolvable).
      13-assertion `test_generalized_pell`; check-c99 clean; valgrind delta 0.
      Docs + changelog updated.
      - Bug found & fixed during validation: `(1,1)` and `(5,3)=(1,1)·ε` for
        D=3,N=-2 were emitted as two overlapping families — the reduction only
        *advanced* into the orthant, never *reduced* a positive solution to
        minimal. Added the ε⁻¹ reduction loop.
- [x] **Systemic**: held-out validation gate — DONE.
      `benchmarks/87-diophantine-integers/heldout.py` (~20 equations from
      standard references, none in `cases.py`) + `validate.py`
      (`make check-diophantine-heldout`): runs Mathilda cold, cross-checks every
      answer against an independent Python brute-force oracle over the same box,
      FAILS (nonzero exit) on any silent wrong answer (a `{}`/finite/parametric
      result the oracle contradicts). Needs only the binary (no sympy); writes
      `HELDOUT_REPORT.md`. Status: OK 18 / DECLINE 2 / WRONG 0; detector verified
      to fail loud via a negative control. Guards all four Tier-1 features
      against regression.
