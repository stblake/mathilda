# Task: Simplify-under-Assumptions assumed-trig stress campaign

Baseline v0.331: **22/50**. Target **49/50** (only #20 out of scope).
Plan: `/Users/user/.claude/plans/similar-to-previous-campaigns-eager-glacier.md`

## Mechanisms (one commit each, bump + tag)

- [ ] **M1** (v0.332) Sign oracle: `assume_sign_nonneg/nonpos_deep` → Reduce bridge + secant decomp.
      Wire into `try_simp_abs` + `try_simp_sqrt_of_square`. → 11,12,17,19,41,42,43,44,45 (Abs)
- [ ] **M2** (v0.333) Generalized `Sqrt[c·f²]→Sqrt[c]·Abs[f]` (Times radicand). → 13,46
- [ ] **M3** (v0.334) Symbolic-integer π-periodicity `Sin[θ+kπ]→(-1)ᵏSin[θ]`. → 34,35,38,39
- [ ] **M4** (v0.335) Sqrt-local radicand prep (half-angle, 1±Cos[2x], Sec²−1). → 14,15,16,18,47
- [ ] **M5** (v0.336) Inverse-trig-of-trig branch reduction (containment + reflection). → 48,49,50
- [ ] **M6** (v0.337) Conjugate Schwarz reflection `Conjugate[f[z]]→f[Conjugate[z]]`. → 3
- [ ] **M7** (v0.338) Reduce-verified inverse multi-angle. → 28,29,30

## Verification per commit
- [ ] stress.m 50-case harness (monotonic, no regressions)
- [ ] existing simplify suites (test_simplify/invtrig/logexp/fullsimplify/trigrat/simp/trigexp_zero)
- [ ] held-out generality cases (not in the 50)
- [ ] new regression tests added
- [ ] make check-messages / check-c99 clean; build clean
- [ ] docs: changelog 2026-10-05.md + simplification.md

## Review
(to be filled at end)
