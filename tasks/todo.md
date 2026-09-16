# DSolve M51 — §2.2.32 corpus (Problems 3101–3200)

Milestone: M51 = §2.2.32. Version bump 0.148 → 0.149.
Definition of done: 0 FAIL, 0 crash; UNEVAL OK only for no-closed-form residue.
Source: https://12000.org/my_notes/solving_ODE/current_version/Ch2.S2.SS32.htm
Naming: DE_examples_2232.m, --label 2.2.32, ctest dsolve_corpus_2_2_32_tests.

Section character: ALL constant-coefficient LINEAR ODEs (2nd/3rd/high-order
nonhomogeneous) with polynomial / exponential / sinusoid forcing + resonance,
plus VoP-only forcing (Sec, Tan, Csc, Log). Core Mathilda strength
(UndeterminedCoefficients + LinearConstantCoefficients VoP). Expect high pass rate.

## Prereq — M50 finalization (DONE)
- [x] §2.1.2 standalone regression check: 647 ≤ 655 (improved to 557 PASS)
- [x] dsolve_corpus_2_2_31_tests gate green (14.75s); main binary rebuilt

## Phase 1 — Fetch & convert
- [x] Fetch Ch2.S2.SS32.htm — title "Problems 3101 to 3200" confirmed
- [x] Convert → DE_examples_2232.m (100 scalar [8 IVP], 0 systems)
- [x] Spot-audit: all 100 records parse. Note 3168 `y'+P(x)y==Q(x)` — converter
      renders P(x)/Q(x) as P*x/Q*x (juxtaposition), a solvable const-linear ODE (not
      the intended arbitrary-function form); watch how it measures. Symbolic params
      in 3129/3130/3155/3184 (n,k,a).

## Phase 2 — Baseline measurement
- [x] Baseline (clean corpus): 98 PASS / 2 UNEVAL / 0 FAIL / 0 crash
- [x] reports/2.2.32.md + .tsv
- [x] MANDATORY BAR met: 0 FAIL / 0 crash

## Phase 3 — Root-cause investigation (fix REVERTED)
- [x] Ranked gap: 3155 + 3165 (both VoP-forcing 2nd-order linear)
- [x] Probed: 3165 CORRECT but cold DSolve >8s (Tan^2 VoP). 3155 was a WRONG answer
      (Sec[ax]/a^2 + Re/Im mush) masked as UNEVAL by leaked-C rule — residual -74.7.
- [x] ROOT CAUSE: dsolve_homog_basis realifies complex roots to Exp[Re x](Cos,Sin)[Im x]
      but only ComplexExpands Re/Im when NUMERIC. Symbolic root ±√(-a²) → unevaluated
      Re/Im → garbage basis that doesn't back-substitute → wrong VoP particular.
- [x] Tried SOLVER fix: symbolic complex root → complex-EXPONENTIAL basis Exp[r x].
      Correct (3155 fixed, y''+9y==0 unchanged) BUT the √(-a²) exponential form slows
      the 2nd-order cascade → 13 symbolic-coeff §2.1.2 cases time out at 8s cold
      (557→544 PASS, non-PASS 647→660 > 655 baseline). NET REGRESSION.
- [x] REVERTED the fix (documented M51 note in dsolve_common.c). Deferred: needs a
      narrower fix that doesn't add cascade latency.
- [x] Residue 4 (reverted/original code, stable): 3155 (masked-wrong, deferred),
      3165 (cold >8s), 3161/3164 (correct VoP at the 8s boundary, timing-flaky)

## Phase 4 — Record & land (no solver change)
- [x] tests/CMakeLists.txt — dsolve_corpus_2_2_32_tests gate (baseline 4, robust)
- [x] STATUS.md §2.2.32 block + M51 wave-history line (rewritten: no solver change)
- [x] README.md row (rewritten)
- [x] DSOLVE_PLAN.md M51 entry + const-coeff catalog "known bug (deferred)" note
- [x] changelog M51 entry (investigated + reverted)
- [x] src/version.h 0.148 → 0.149
- [x] memory: symbolic-complex-root wrong-answer + deferred-fix lesson
- [x] make -j + make check-c99 clean (reverted binary)
- [ ] §2.1.2 reverted-binary re-confirm ≤655 (running; revert is byte-identical to main)

## Review
- §2.2.32 = M51. All const-coeff linear ODEs w/ forcing. 96/100, 0 FAIL, 0 crash.
  NO solver change lands.
- Key finding: a masked WRONG answer (3155) traced to the symbolic-complex-root
  homogeneous basis. The obvious exp-basis fix is CORRECT but regresses 13 §2.1.2 cases
  by slowing the hot 2nd-order cascade past the 8s cold budget → reverted, deferred.
- Lessons (saved to memory): (1) the verifier's leaked-C[k]→UNFIT rule MASKS wrong
  general solutions as UNEVAL — numerically check residual SIGN to distinguish
  correct-slow from wrong-masked. (2) A homogeneous-basis change is on the hot 2nd-order
  path — always measure cold §2.1.2, not just the target case.
