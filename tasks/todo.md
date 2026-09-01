# DSolve / Kovacic stress-test hardening (Tiers A+B+C) — DONE

Plan: `/Users/user/.claude/plans/the-following-set-of-moonlit-wirth.md`

## Tier A — Robustness & correctness
- [x] A1. Early complexity gate in `dsolve_kovacic_try` (denom degree / pole count) — In[12] hang FIXED (fast series)
- [x] A2. Silence `Integrate::nonelem` on speculative Kovacic integrals — In[5] message GONE
- [x] A3. Frobenius expansion at a shifted ordinary point — In[7]→O[x-2]^7, session-2 In[1]→O[x-1]^7 (was bare-uneval); warnings muted
- [x] A4. Absorb cosmetic constant factor from `second_solution` into C[2] — In[2] now clean Euler form

## Tier B — Kovacic Case-1 completion (full classical algorithm)
- [x] B1. `[√r]_c` for even-order poles ≥ 4 (`kovacic_pole_data`, lifted `mult>2` cap)
- [x] B2. Non-constant `[√r]_∞` for r growing at ∞ (`kovacic_inf_data`, x→1/y; b uses `([√r])²`)
- [x] B3. (s∞, mask) enumeration + degree bound + `solve_monic_P`; **two independent families first**, then reduction of order
- [x] B4. Hang guard: decline reduction-of-order integral in Erf/Erfi/… (`has_hang_special`)
  - NOTE: In[1] (Kovacic's example) finds y1 but y2 is genuinely non-elementary → graceful series
  - NOTE: session-2 In[1] has irrational pole exponents (√2) → not Case-1 → graceful series
  - PROVEN on `y''−(x²+3/(4x²))y==0` (growth at ∞, old δ≥2 rejected it)

## Tier C — Bessel/Airy change-of-variable recognizer (In[11])
- [x] C1. Recognize `P==0, Q == A x^m` → `√x Z_{1/(m+2)}(...)`; J/Y for A>0, I/K for A<0. In[11] FIXED
  - NOTE: pipe protocol needs INTEGER ids; `"id":"11"` (string) is silently dropped (test-harness gotcha, cost a detour)
  - NOTE: permissive verify accepted a WRONG (sign-flipped) Bessel form — always numerically spot-check special-fn recognizers

## Verification
- [x] Regression tests: `test_dsolve.c` (Bessel-reducible, high-degree no-hang), `test_dsolve_m5_stress.c` (growth closed-forms, growth no-hang) — all pass
- [x] Manual NDJSON checks on all 14 corpus cases — closed form where elementary, graceful series otherwise, no hangs, clean stderr
- [x] `make check-c99` passes (rc=0); valgrind is Linux-CI (macOS baseline noise); new paths run 1.1s clean
- [x] Docs: `docs/spec/builtins/calculus.md` (Kovacic / SpecialFunctionForm / PowerSeries rows) + `docs/spec/changelog/2026-08-31.md`
- [x] Existing suites (`dsolve_tests`, `dsolve_m5_stress_tests`) still pass; stress back to 7s

## Review
- **Root causes were shared, not per-case.** Buckets: hang (missing degree gate),
  bare-unevaluated (Frobenius only expanded at x=0), fell-to-series (Case-1
  incompleteness OR genuinely non-Liouvillian), cosmetic (unabsorbed constant),
  leaked warning (speculative integral).
- **Reality check on the corpus:** In[1] and session-2 In[1] — the two "crown
  jewel" targets — turned out NOT to have elementary general solutions (In[1]'s 2nd
  solution is non-elementary; session-2 In[1] has irrational pole exponents, not
  Case-1). Tier B is still a correct, verified generalization (proven on a
  constructed growth-at-∞ case), but its corpus impact is the guard behaviour
  (fast, clean declines) rather than new closed forms for those two.
- **Two soundness lessons captured to memory:** (1) the permissive verify accepts a
  wrong special-function candidate whose residual zero_test can't disprove — always
  numerically spot-check a recognizer; (2) the NDJSON pipe silently drops string
  `id`s (test-harness gotcha).
- **Elegance:** kept the pinned-method contracts intact (shifted-point fallback is
  cascade-only; growth reduction-of-order declines special functions instead of
  hanging), so no existing behaviour regressed.
