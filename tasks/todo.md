# Task: Solve every solvable Wikipedia Diophantine equation in `Solve[..., Integers]`

Plan: `/Users/user/.claude/plans/let-s-ensure-that-we-compiled-bubble.md`

Closes 6 gaps found by empirically testing Mathilda 0.074 against the Wikipedia
"Diophantine equation" page. All prerequisites (JacobiSymbol, LLL, number-field
layer, Baker/Thue machinery) already exist.

## Work items (approved order: cheapest correctness win first)

### 3 — L: Erdős–Straus / symmetric-Egyptian (small edit) ✅ DONE
- [x] `si_solve_reciprocal`: symmetric case → canonical order + emit permutations
- [x] Verified: `4/5` → 12 perms; ordered form unchanged; `1=1/x+1/y+1/z` → 10 tuples

### 4 — M: Frye one-sided-box robustness (small edit) ✅ DONE
- [x] Derive magnitude bound `|w|<=(Σ hi^4)^{1/4}`; positive search + emit BOTH
      signs filtered by si_verify (sound; also enables unconstrained w)
- [x] Verified: witness window unchanged; negative-only box → -422481 (0.003s);
      tight positive window keeps only +; engagement on `w<500000` confirmed

### 1 — Ternary quadratic / Legendre solver (flagship, closes D/F/G) ✅ CODE DONE
- [x] New `src/solve/solveint_ternary.c` — `si_solve_ternary_quadratic`; hooked after genpell
- [x] Legendre proof (per-prime Legendre) → `{{0,0,0}}` when it fails
- [x] Witness + chord/tangent + tangent-family; 16-branch sign/swap orbit (deduped)
- [x] Single-representation gate (≤1 prime ≡1 mod4); multi-rep declines (k=65)
- [x] Verified end-to-end vs brute: Pythagorean 113/113, 2z² 97/97 COMPLETE+SOUND
- [ ] Generalise `benchmarks/87/validate.py` to multi-parameter family enumeration (test phase)

### 2 — Ramanujan–Nagell unbounded solver (closes K) ✅ CODE DONE
- [x] New `src/solve/solveint_rn.c` — Route A (imag-quadratic → Lucas → BHV n≤40)
- [x] Promoted `SIExpTerm`/`si_exp_collect` to shared (header)
- [x] Hooked in `solveint.c` after `si_solve_exponential`
- [x] Gate (b=2, D≡7 mod 8, squarefree, h=1) + Lucas cross-check
- [x] Verified: `2^n-7==x^2` → {3,4,5,7,15}; alt spelling; declines b=3, D=15

## Verification
- [x] C tests in `tests/test_solve_integers.c` (4 new, property-based) — ALL PASS
- [x] Held-out cases in `benchmarks/87/heldout.py` (7 new) + validate.py multi-param fix;
      `make check-diophantine-heldout` GREEN (No silent wrong answers)
- [x] `make check-c99` (exit 0)
- [x] valgrind delta 0 (identical leak totals vs baseline; no new-file frames)
- [x] Docs: `solutions-of-equations.md`, changelog `2026-08-17.md`, SOLVE_INTEGERS.md "Done"
- [x] End-to-end: re-run all Wikipedia gap equations → all return Target column

## Review

All 6 gaps closed; all Wikipedia "solvable" Diophantine equations now solve.

**Delivered:**
1. `si_solve_ternary_quadratic` (new `solveint_ternary.c`) — Legendre proof;
   trivial-only `{{0,0,0}}` when unsolvable; complete sign/swap-orbit + tangent
   family for solvable single-representation k. Closes Pythagorean (D), `3z²` (F),
   `2z²` (G).
2. `si_solve_ramanujan_nagell` (new `solveint_rn.c`) — Route A (imag-quadratic →
   Lucas → BHV n≤32). Closes `2ⁿ−7==x²` (K).
3. `si_solve_reciprocal` symmetric-permutation path — closes Erdős–Straus without
   ordering (L).
4. `si_solve_biquadrate_frye` one-sided box + both-sign emission — closes Frye
   engagement (M).

**Verification:** 4 new C tests + 7 held-out cases (all green, 0 wrong answers),
`make check-c99` clean, valgrind delta 0, full `solve_integers_tests` /
`thue_tests` / `solve_tests` pass.

**Key correctness findings (non-obvious):**
- Ternary single-witness parametrization is INCOMPLETE for k with ≥2 primes ≡1
  (mod 4) (multiple sum-of-two-squares reps → multiple classes, e.g. 65=5·13);
  gate to single-representation, decline otherwise (safe).
- Chord parametrization couples x/y signs and gives 2·P0 (not P0) at the tangent;
  need the full sign/swap orbit + a separate tangent-line family `C·P0` for
  completeness (validated 16→deduped branches vs brute force).
- Frye: forcing `w>0` on a one-sided box would silently drop `w<0` solutions
  (wrong); instead bound `|w|` from the summand box, search positive, emit both
  signs filtered by `si_verify`.
