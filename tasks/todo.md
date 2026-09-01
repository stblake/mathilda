# Task: Complete Kovacic Case 1 (apparent singularities) — closed form over series

Goal: `DSolve` returns the elementary closed form for Legendre/Chebyshev/Gegenbauer/Jacobi
(currently a truncated series, and Chebyshev hangs). Fully algorithmic via a complete classical
Kovacic Case 1 for arbitrary rational `r` (apparent singularities, higher-order + complex poles).
All engine work in `src/calculus/dsolve_kovacic.c`.

## Steps

- [x] 1. Read `dsolve_common.h` helper API (dsolve_analyze_roots enumerates poles incl. complex).
- [x] 2. Factor out `solve_monic_P(theta, R2, d, x, counter)` (R2 precomputed per θ).
- [x] 3. Pole locals via Limit[(x-c)^2 r] (order ≤2 + complex); [√r]_∞=0 for δ≥2.
- [x] 4. Add `kovacic_case1_general`: enumerate poles+signs, dedup degenerate, search d,
        collect y1 candidates, pick cleanest by LeafCount, y-level reduction of order.
- [x] 5. Wire into `dsolve_kovacic_try` after pure-Riccati Case 1, for all rational r (degd≥1).
- [x] 5b. Fixes: nu==0 guard in solve_ansatz (no more Solve::ivar); drop numeric_verify (exact
         by construction); y-level assembly (avoid √·ArcTanh Simplify blow-up); skip final
         realify (avoid Gegenbauer 7s). All class members now <0.5s.
- [x] 6. Case-2 hang: root-caused to the σ-solve (ds_solve) on complex-pole (degree-≥2
        factor) systems; guarded by "all irreducible factors linear" gate → declines to series.
        Also made Case 1c decline fast via numeric degree test (complex α).
- [x] 7. REPL-verified: Legendre n=1..5, Chebyshev n=2/3, Gegenbauer, complex poles — fast+correct.
- [x] 8. Regression: single-pole, poly-r, Hermite series, non-Liouvillian declines — all intact.
- [x] 9. Tests: 6 new cases in `tests/test_dsolve.c` (Legendre 1/2, Chebyshev, complex poles,
        auto-closed-form, no-hang). All 3 DSolve suites (incl. m5 Case-2) pass.
- [x] 10. Gates: `make check-c99` clean; valgrind — no new leak (only pre-existing
         solve_ansatz→Together baseline ≈13.6KB); broad 9-ODE hang sweep clean.
- [ ] 11. Docs: calculus.md, changelog 2026-08-31, DSOLVE_PLAN.md §1c, docstring; memory note.

## Review

**Outcome:** `DSolve[(1-x^2)y''-2x y'+2y==0, y, x]` now returns the elementary closed
form `C[1] x + C[2](1 - x ArcTanh[x])` (Legendre P₁/Q₁) instead of a truncated series —
matching Maple/Mathematica. Fully algorithmic via the completed Kovacic Case 1, not a
Legendre pattern.

**What changed** (all in `src/calculus/dsolve_kovacic.c`):
- `kovacic_case1_general` — classical Case 1 apparent-singularity construction for any
  rational `r`: poles (incl. complex) via `dsolve_analyze_roots`, local exponents
  `α_c=(1±√(1+4b_c))/2`, monic `P` of the classical degree `d=α_∞−Σα_c`, `z1=P Exp[∫θ]`.
- `solve_monic_P`, `numeric_nonneg_int` helpers; `solve_ansatz` nu==0 short-circuit.
- Case-2 guard: skip when a denominator factor has degree ≥ 2 (avoids the σ-solve hang).
- Wired as Case 1c after pure Case 1, before Case 2.

**Key engineering lessons** (saved to memory
`project_dsolve_kovacic_case1_apparent_singularity`): y-level reduction of order (not
z-level — avoids a 50 s √·ArcTanh Simplify blow-up); numeric degree test (complex α
symbolic Simplify hangs); skip final realify (Gegenbauer 17 s → 0.1 s); cleanest-y1 by
LeafCount; Case-2 complex-pole guard (never TimeConstrained in-engine).

**Results:** Legendre n=1..5, Chebyshev n=2/3, Gegenbauer, complex poles — all closed
form, all < 0.5 s, residuals ~1e-14. Hermite/Laguerre correctly stay series (2nd sol
non-elementary). `(x³+1)y''+xy'+y` no longer hangs (1.1 s → series).

**Verification:** all 3 DSolve suites pass (6 new tests added); `make check-c99` clean;
valgrind — no new leak (only the pre-existing ~13.6 KB `solve_ansatz→Together`
baseline); broad 9-ODE hang sweep clean.
