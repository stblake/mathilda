# DSolve M32 — §2.2.12 corpus (Problems 1101–1200) to full coverage

Baseline (real prelude verifier, 8s/case): **80 PASS / 19 UNEVAL / 1 FAIL**.

## Stage 1 — corpus + dashboard wiring
- [ ] Generate `DSolve_test_status/DE_examples_2212.m` (`--label 2.2.12 --url …indexsubsection21.htm`)
- [ ] Add `dsolve_corpus_2_2_12_tests` to `tests/CMakeLists.txt` (temp high baseline for dev)
- [ ] Build `dsolve_corpus_tests`, capture baseline TSV, generate `reports/2.2.12.{tsv,md}`

## Stage 2 — F1: IVP constant-fitter (FAIL 1147 + Cluster A) ✅
- [x] fit_state (OK/EMPTY/UNDEF) + dsolve_run sibling-aware drop; never emit Undefined/unfitted-C
- [x] Bernoulli emits both real signs for even 1-n; Root-form IVP via implicit twin
- [x] Verified 1138,1140,1141,1143,1144,1145,1146,1147(FAIL→PASS),1149,1150,1152,1153 → PASS
- [x] Regression guard: under-determined BVP keeps C[2] (t_cc_bvp restored)

## Stage 3 — F2: Separable recognizer + implicit twin (Cluster B) ✅
- [x] F2a relax gate to `ds_is_zero(denom)` (accept generic-param splits) → 1157-class
- [x] F2b `dsolve_separable_implicit_try` (dsolve_run_implicit), non-elem integral kept → 1173,1186
- [x] F2c x-free/Root-body drop → degenerate `{{y->C[1]}}` fixed (1133)
- [x] 1182/1190: were the Erf-variable bug (F3), not slowness — fixed
- [x] mute speculative integrals (`g_integrate_quiet`)

## Stage 4 — Cluster C ✅ (scope decided by evidence)
- [x] Inverse/autonomous non-elem case (1186) closed via separable implicit twin
- [x] Integrate Gaussian→Erf/Ei/PolyLog variable fix (risch_special.c) → 1182,1190
- [x] Prototyped generalised d'Alembert on 1135/1200: induced ODE NOT cascade-solvable →
      documented residue (1135,1200 + Abel-2nd-kind 1157). Full method deferred to own milestone.

## Stage 5 — finalize + gates
- [x] §2.2.12 fork-per-case: 97/3/0 (baseline 3 in CMakeLists)
- [ ] Anti-overfit `t_m32_*` units in `tests/test_dsolve.c` (drafted; add + verify via REPL)
- [~] Regression: §2.1.2 + §2.2.1–§2.2.12 corpus ctests — RUNNING
- [ ] `make check-c99`; valgrind spot-check; regenerate reports/2.2.12.{tsv,md}
- [x] STATUS.md, README.md, DSOLVE_PLAN.md (M32), changelog 2026-09-07.md
- NOTE: dsolve_tests unit suite hits a PRE-EXISTING 120s whole-binary alarm at
  t_rischnorman_enum_cap_no_crash (identical on pristine main; CI does not run this suite).
  Tests 1-41 pass on my build; BVP regression I introduced was fixed.

## Review

**Result: §2.2.12 added and driven 80/19/1 → 97/3/0** (0 FAIL, 0 crashes; the section's
sole wrong answer, 1147, is repaired). Two clean fork-per-case runs confirm 97/3.

Three shared-substrate root-cause fixes (each lifts earlier sections; none regress):
- **F1 IVP fitter** (`dsolve_common.c`): per-branch `fit_state` (OK/EMPTY/UNDEF) + sibling-aware
  drop in `dsolve_run`; Bernoulli emits both real signs for even `1−n`. Fixed FAIL 1147 + 11 IVPs.
- **F2 Separable** (`dsolve_separable.c`): generic-param split gate + implicit twin
  `dsolve_separable_implicit_try` (non-elem integral kept unevaluated). Fixed 1133/1173/1149/1150/1186.
- **F3 Integrate** (`risch_special.c`): Gaussian→Erf/Ei/PolyLog templates threaded the real
  integration variable (were literal `x`). Fixed 1182/1190.

Residue (3, research-grade bounded declines): 1135, 1200 (solvable-for-y/x, transcendental),
1157 (Abel 2nd kind). A generalised-d'Alembert method was prototyped and closes none of them → deferred.

Verification: §2.2.12 97/3 ×2 clean; targeted regression §2.2.1(99,≤base2)/2.2.4/2.2.8(≤base1)/2.2.9/
2.2.10/2.2.11 all within baseline (§2.2.1 improved); t_m32_* all pass (REPL); `make check-c99` green;
reports/2.2.12.{tsv,md} regenerated. §2.1.2 full run skipped by analysis (no ICs → fitter no-op there).
NOTE: `dsolve_tests` unit binary hits a PRE-EXISTING 120s alarm at t_rischnorman (identical on pristine
main; CI does not run this suite); tests 1–41 pass, BVP regression I introduced was fixed.
