# Task: Simplify-under-Assumptions assumed-trig stress campaign

Baseline v0.331: **22/50**. Target **49/50** (only #20 out of scope).
Plan: `/Users/user/.claude/plans/similar-to-previous-campaigns-eager-glacier.md`

## Mechanisms (one commit each, bump + tag)

- [x] **M1** (v0.332) Sign oracle: `assume_sign_nonneg/nonpos_deep` → Reduce bridge + secant decomp.
      DONE: 22→33/50 (cleared 11,12,17,18,19,28,41-45), sound, no regressions, leak-clean, tagged.
- [x] **M2** (v0.333) Generalized `Sqrt[c·f²]→Sqrt[c]·Abs[f]` (Times radicand). DONE 33→35/50.
- [x] **M3** (v0.334) Symbolic-integer π-periodicity `Sin[θ+kπ]→(-1)ᵏSin[θ]`. DONE 35→39/50.
- [x] **M4** (v0.335) Sqrt-local radicand prep (half-angle, 1±Cos[2x], Sec²−1). DONE 39→43/50.
- [x] **M5** (v0.336) Inverse-trig-of-trig branch reduction (containment + reflection). DONE 43→46/50.
- [x] **M6** (v0.337) Conjugate Schwarz reflection `Conjugate[f[z]]→f[Conjugate[z]]`. DONE 46→47/50.
- [x] **M7** (v0.338) Reduce-verified inverse multi-angle. DONE 47→49/50 (cleared 28,29,30).

## Verification per commit
- [ ] stress.m 50-case harness (monotonic, no regressions)
- [ ] existing simplify suites (test_simplify/invtrig/logexp/fullsimplify/trigrat/simp/trigexp_zero)
- [ ] held-out generality cases (not in the 50)
- [ ] new regression tests added
- [ ] make check-messages / check-c99 clean; build clean
- [ ] docs: changelog 2026-10-05.md + simplification.md

## Review

**Outcome: 22/50 → 49/50** (v0.331 → v0.338, 7 commits, each tagged). Only #20
(`Log[Sec+Tan] = ArcTanh[Sin]`, a one-off gudermannian identity) remains — it was
explicitly scoped out (no general algorithm).

Per-mechanism: M1 sign oracle 22→33 · M2 Sqrt extraction →35 · M3 π-periodicity
→39 · M4 Sqrt-local radicand prep →43 · M5 inverse-of-trig branch reduction →46 ·
M6 Conjugate reflection →47 · M7 Reduce-verified multi-angle →49.

All changes are general, sound mechanisms (not case pattern-matches): every one
has paired soundness tests confirming it DECLINES where the identity is not
universally valid (sign-changing/unbounded regions, wrong branch, missing
reality). No `PossibleZeroQ`; all `0`s come from exact rewrites or Reduce proofs.

Verification: 21 simplify/trig/refine/reduce suites green (hard_fail=0); the only
soft-FAILs (1 simplify, 2 simp, 1 logexp, all NDEBUG-elided) are **pre-existing**
and outside this campaign's domain (confirmed `a^p^q` at v0.331 via worktree;
`Coth` TrigToExp and a `(x^2)^(3/2)` printer cosmetic likewise). `make
check-messages` / `check-c99` clean. valgrind: zero leaks trace to any new
function. New regression + soundness tests in test_logexp_simplify.c (M1-M4) and
test_invtrig_simplify.c (M5-M7).

**Pre-existing NDEBUG-masked issues surfaced AND FIXED (v0.339):**
1. `Simplify[(a^p)^q, a>0]` → `a^(p q)` was UNSOUND (needs `p ∈ Reals`). FIXED:
   replaced the two context-free string rules with a structural rule gated on
   `prov_nn(base) && prov_re(inner_exp)`. Genuine code bug.
2. `TrigToExp/ExpToTrig[Coth[x]]` — stale test golden values from an old sign-flip
   bug since fixed in code; corrected the expectations (code was already right).
3. `(x^2)^(3/2)` printer — test expected ambiguous `x^2^(3/2)` (reparses as
   `x^(2√2)`); corrected the expectation (printer was already right).

All three asserts were `-DNDEBUG`-elided, so they'd been failing silently; the
14-suite sweep now shows 0 soft-FAILs (was 4).
