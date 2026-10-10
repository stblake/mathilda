# DSolve corpus wave M66 — §2.2.39 (Problems 3801–3900)

Plan: `/Users/user/.claude/plans/let-s-continue-our-implementation-keen-hollerith.md`

## Phase A — land & audit the corpus
- [x] Fetch upstream §2.1.39 `Ch2.S1.SS39.htm` (browser-UA curl) → scratchpad
- [x] Convert → `DSolve_test_status/DE_examples_2239.m` (label 2.2.39)
- [x] `make check-corpus-indvar` green; eyeballed — clean
- [x] 100 records; 8 scalar (0 IVP) + 92 systems (53 2×2, 32 3×3, 7 4×4)

## Phase B — measure baseline
- [x] Build `dsolve_corpus_tests`
- [x] Baseline run → 95/100, 5 UNEVAL (3804,3805,3832,3886,3891), 0 FAIL/crash
- [x] Non-PASS diagnosed from inside cascade (pinned methods, profiler)

## Phase C — diagnose & fix (FAIL > crash > UNEVAL)
- [x] Classified: 3805 = latency (Log·trig ∫ falls to PMS 9.4s); other 4 hard residues
- [x] Fix: new `Integrate` LogByParts stage (`integrate_logbyparts.c`) — 3805 9.4→0.35s
- [x] §2.2.39 → 96/100; integrate/risch suite (35) green, 0 regressions
- [x] A/B kill-switch: 4 failing dsolve unit/stress tests are PRE-EXISTING on HEAD

## Phase D — lock in & document
- [x] Add `dsolve_corpus_2_2_39_tests` gate at baseline 4
- [x] STATUS.md section + wave-history M66; DSOLVE_PLAN.md M66; README rows (37/38/39)
- [x] Changelog `docs/spec/changelog/2026-10-05.md`; spec `calculus.md` LogByParts note
- [x] Bump `src/version.h` → 0.343
- [x] Targeted regression on VoP/trig/Log sections (19/25/27/32/34/35/36/38/39)
- [x] Regenerate `reports/2.2.39.{tsv,md}` post-fix → 96/100, 4 non-PASS
- [x] valgrind: new file in 0/31581 leak records, no invalid access; check-c99 ✓ / check-messages ✓
- [ ] Commit + tag v0.343

## Review

**Result:** §2.2.39 (Problems 3801–3900) landed at **95 → 96/100, 0 FAIL, 0 crash**,
gate `dsolve_corpus_2_2_39_tests` at baseline **4**.

**One general fix:** new `Integrate` stage `LogByParts`
(`src/calculus/integrate_logbyparts.c`) — recognises `c Log[g] K` with a
trig/hyperbolic cofactor and does one integration by parts, closing the `Log·trig`
family to `Log·trig + Si/Ci` instead of the ~4.7–12 s `ParallelMixedSpecial`
search. `3805` (`y''+4y==Log[x]`): 9.4 s → 0.35 s. Placed after the cheap stages
(cascade `if(!result)` short-circuit ⇒ sees only what they declined, perturbs
nothing already fast) and before the heavy tail; recursive sub-integrals raise
`g_integrate_no_special` so a non-closing residual declines fast; exact `Simplify`
diff-back is the sole acceptance test.

**Regression-free (A/B kill-switch proven):** integrate/risch suite 35/35 green;
every affected corpus section identical with/without the stage; valgrind — new
file in 0 leak records, no invalid access.

**Residue (4, honest hard classes):** 3804 (variable-coeff Airy-inhomogeneous),
3832/3891 (variable-coeff 2×2/3×3 systems), 3886 (const-coeff 3×3, irreducible
cubic-`Root` spectrum churn).

**Pre-existing, NOT this wave (flagged for a drift-reconciliation pass):** gates
`dsolve_corpus_2_2_19_tests` (8 vs stale baseline 5) and `_2_2_35_tests` (3 vs 2)
have drifted across M40–M65; unit/stress `dsolve_tests`, `dsolve_m34/m62/m63_stress`
are red on HEAD. All A/B-confirmed independent of M66; baselines left untouched
(raising without root-cause could mask a real regression).
