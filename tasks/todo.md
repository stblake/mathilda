# DSolve corpus wave M46 — §2.2.27 (Problems 2601–2700)

## Phase A — Corpus generation & registration
- [x] Generate `DSolve_test_status/DE_examples_2227.m` via converter (100 records)
- [x] Converter-invariance gate: converter untouched → prior sections invariant
- [x] Register `dsolve_corpus_2_2_27_tests` ctest in `tests/CMakeLists.txt`

## Phase B — Baseline measurement
- [x] Build `dsolve_corpus_tests`, run §2.2.27 → baseline 97/100 (0 FAIL, 1 CRASH, 2 UNEVAL)
- [x] Bucket report `reports/2.2.27.md`

## Phase C — Root-cause fixes (one general fix; no overfit)
- [x] CRASH 2690 root-caused: prelude `dsFreeParams` collected relational operators
      (LessEqual/Less) from a CHAINED-inequality Piecewise condition
      (Inequality[0,LessEqual,t,Less,2]) and substituted numbers → corrupt expr →
      SIGSEGV. Fix: exclude relational/logical/piecewise operators. 2690 CRASH→PASS.
- [x] Anti-overfit / regression: §2.2.15/2.2.16/2.2.25 (the only piecewise sections)
      hold at baseline; non-piecewise sections provably unchanged. 1430 flake noted.
- [x] §2.2.27 → 98/100, baseline set to 2 (2621/2641 sympy=False residue)

## Phase D — Documentation, gates, version
- [x] STATUS.md §2.2.27 block + M46 wave-history bullet
- [x] README.md section row
- [x] DSOLVE_PLAN.md M46 milestone
- [x] docs/spec/changelog/2026-09-14.md note (top, newest-first)
- [x] src/version.h 0.143 → 0.144
- [x] Verify: gate passes (baseline 2, 52.8s); check-c99 exit 0; REPL spot-checks

## Review

**Result: §2.2.27 = 98/100, 0 FAIL, 0 crash.** (baseline 97/100 had 1 CRASH + 2 UNEVAL.)

**One general root-cause fix** (`dsolve_corpus_prelude.m`, `dsFreeParams`): a chained
inequality in a Piecewise/step condition (`0<=t<2` ≡ `Inequality[0, LessEqual, t, Less, 2]`)
puts comparison operators in argument position; the verifier collected `LessEqual`/`Less` as
free parameters and substituted numbers → corrupted Piecewise → SIGSEGV on differentiation
(2690). Fix = exclude relational/logical/piecewise operator symbols (never ODE parameters).
Monotone-safe; 2690 CRASH→PASS. **No C code changed** (only version.h macro).

**Regression:** the only sections with chained-inequality Piecewise conditions
(§2.2.15/2.2.16/2.2.25) all hold at baseline; every other section is provably unaffected
(no such operators in their residuals). §2.2.15-1430 is a pre-existing flaky near-8s system
solve (unrelated).

**Residue 2** (both sympy=False, bounded declines): 2621 (exact → non-elementary integrating
factor), 2641 (transcendental coefficient). Not chased — genuine no-closed-form cases.

**Notable:** the CAS itself is robust on legitimate chained-inequality Piecewise D/Integrate;
Emden–Fowler 13/13 (Airy/Bessel), systems 4/4, ortho-poly 7/7 all solved out of the box.

**Files:** DE_examples_2227.m + reports/2.2.27.{tsv,md} (new); dsolve_corpus_prelude.m,
tests/CMakeLists.txt (new gate, baseline 2), STATUS.md, README.md, DSOLVE_PLAN.md,
changelog, version.h (edited).
