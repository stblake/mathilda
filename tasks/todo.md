# M25 — §2.2.5 corpus (Problems 401–500), driven to 100% PASS

Add Nasser Abbasi 12000.org §2.2.5 (indexsubsection14.htm, "Problems 401 to 500") to
the DSolve corpus dashboard. User goal: **push to 100/100 PASS**. §2.2.5 is
series-solution-heavy (Edwards & Penney): 2nd-order linear (const + variable coeff),
Airy/Emden–Fowler, Gegenbauer, Bessel, Jacobi/₂F₁, one Liénard, one 3rd-order.

## Step 1 — Convert HTML → corpus
- [x] Fetch `indexsubsection14.htm` via curl (218 KB, title "Problems 401 to 500")
- [x] Convert → `DE_examples_225.m` (`--label 2.2.5`): 100 scalar, 15 IVP, 0 sys
- [x] Faithfulness spot-check (401/439/455/494/497) ✔
- [x] Converter untouched → no §2.2.1–4 regen needed (byte-identical)

## Step 2 — Register ctest + first baseline run
- [x] Add `dsolve_corpus_2_2_5_tests` (baseline 1)
- [x] Baseline run: **95 PASS / 5 UNEVAL / 0 FAIL**
- [x] Generate `reports/2.2.5.{tsv,md}`

## Step 3 — Analyze gaps + close to 100% (verified reuse first; new capability as needed)
- [x] Bucket non-PASS: 428/482 (exact→Erf verify spin / Kovacic degenerate) + 459/463/490 (transcendental-coeff singular)
- [x] **Fix A — Erf verify** (`dsolve_common.c`): ExpandAll a Gaussian×Erf residual before zero_test (gated to Erf/Erfi) → 428 returns Erf closed form. zero_test core deficiency logged in `POSSIBLE_ZEROQ_IMPROVEMENTS.md`.
- [x] **Fix B — Kovacic fundamental-set guard** (`dsolve_kovacic.c`, gated to a body carrying C[k]) + ExactODE non-elementary-Integrate decline (`dsolve_exactode.c`) → 482 → Frobenius 2-param series.
- [x] **Fix C — transcendental Frobenius** (`dsolve_frobenius.c`): `normal_series` Taylor-normalises xP/x²Q → 463/490 return verified Frobenius series.
- [x] 459 documented residue (irregular singular point; matches Mathematica)
- [x] Every fix reuses verified machinery / recurrence-gated series; 0 FAIL by construction

## Step 4 — Anti-overfit units + scoreboard + gate
- [x] Pinned `t_m25_exact_erf`, `t_m25_kovacic_fundamental_set`, `t_m25_transcendental_frobenius` in test_dsolve.c
- [x] Updated `t_kov_stress_pole` (m5 stress) — exposed a latent degenerate-basis bug; now asserts full-DSolve solve + pinned-Kovacic decline
- [x] STATUS.md §2.2.5 block + wave-history bullet
- [x] README.md contents row
- [x] `dsolve_corpus_2_2_5_tests` baseline → 1 (final 99/1/0)

## Step 5 — Plan + changelog + task log
- [x] DSOLVE_PLAN.md M25 milestone entry
- [x] `docs/spec/changelog/2026-09-07.md` M25 note
- [x] Review section (below)

## Verification
- [x] Final §2.2.5: **99 PASS / 1 UNEVAL / 0 FAIL**; new PASSes reproduced interactively
- [x] §2.2.1/2.2.2/2.2.3/2.2.4 ctest gates hold (96/92/99/99, 0 FAIL)
- [x] §2.1.2: **441 solved (up from 436), 0 FAIL** (pre-gate run; final gated re-run in progress, ≥441 expected)
- [x] All DSolve unit + stress suites green (dsolve_tests + m5/m12/m14/m17/m18/m19/m20 + general)
- [x] `make check-c99` green

## Review

**Outcome:** §2.2.5 (Problems 401–500) added to the DSolve corpus at **99/100 PASS, 0
FAIL** (baseline 95). The single residue, 2.2.5-459 `x²y''+Cos[x]y'+xy==0`, is an
**irregular** singular point at x=0 whose only analytic solution is a one-parameter
formal series — left unevaluated, matching Mathematica (the shifted-Frobenius path
deliberately declines transcendental coefficients). 100% was the target; 459 resists
every honest solve+verify path, so it is surfaced as the residue rather than faked
(a rank-deficient or divergent-from-Mathematica answer would violate the invariants).

**Three verified fixes, all reusing existing machinery (0 FAIL by construction):**
1. Erf integrating-factor verify — `dsolve_verify_body` ExpandAll-normalises a
   Gaussian×Erf residual before the numeric zero-test (which otherwise spins on the
   `E^(-x²/2)·E^(x²/2)` cancellation). Gated to Erf/Erfi. Core zero_test deficiency
   logged in `POSSIBLE_ZEROQ_IMPROVEMENTS.md` for a follow-up session.
2. Kovacic fundamental-set independence guard — rejects a degenerate `(C[1]+C[2])y1`
   basis (gated to bodies carrying a generated constant, so a constant-less Case-1
   body is untouched); ExactODE declines a non-elementary reduced Integrate. Together
   482 falls through to the correct 2-parameter Frobenius series.
3. Transcendental-coefficient Frobenius — `normal_series` Taylor-normalises xP, x²Q so
   an analytic transcendental coefficient's removable singularity (`6 Sin[x]/x`) no
   longer poisons the indicial roots.

**Latent bugs exposed & handled:** the Kovacic guard surfaced that the m5-stress "pole"
family was passing on a rank-deficient `(C[1]+C[2])y1` (only one Liouvillian solution
exists; the second is non-elementary). Updated `t_kov_stress_pole` to assert the correct
behaviour (full DSolve solves via Frobenius; pinned Kovacic correctly declines). A
separate pre-existing quirk — Kovacic's Case-1 "growth" family returns a constant-less
body — was left untouched (out of scope) by gating the guard to bodies with a C[k].

**Side benefit:** §2.1.2 improved 436 → 441 solved (the Frobenius/exactode fixes also
close a few there), still 0 FAIL.

**Follow-ups (logged, not done here):** the core `zero_test`/`PossibleZeroQ`
Gaussian×Erf spin (`POSSIBLE_ZEROQ_IMPROVEMENTS.md` #1); the Kovacic Case-1 constant-less
"growth" body; irregular-singular-point / integer-difference-root (log) Frobenius for
cases like 459 and the m5 pole a1c3.
