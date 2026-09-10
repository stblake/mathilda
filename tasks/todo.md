# M34 — DSolve §2.2.14 corpus (Problems 1301–1400) → full coverage

Baseline **90/100** → **99/100, 0 FAIL** (residue 1360 = forced Duffing, no CAS solves it).

## Step 0 — Land the corpus
- [x] `DSolve_test_status/DE_examples_2214.m` (100 records)
- [x] `dsolve_corpus_2_2_14_tests` in `tests/CMakeLists.txt` (baseline → 1)
- [x] STATUS.md `## Section 2.2.14` block + wave-history (M33 + M34 backfilled)
- [x] `reports/2.2.14.{md,tsv}`

## Fixes (all verified against corpus + m34 stress)
- [x] **A/B — robust VoP + numeric-zero verify short-circuit** (`dsolve_common.c`) → 1337/1341/1350/1354
- [x] **C — bounded-Kovacic complex-pole gate** (`dsolve_kovacic.c`, alarm-free `time()`) → 1392/1393; forced-g inert VoP → 1350
- [x] **D — exact plain-symbol first integral + SeriesData IVP Normal-fit + 2nd-order FIT_UNDECIDED fall-through** (`dsolve_exactode.c`, `dsolve_common.c`) → 1384/1381
- [x] **E — IC-point Frobenius series** (`dsolve_frobenius.c`) → 1385

## Anti-overfit + hygiene
- [x] `test_dsolve_m34_stress.c` (6 forward-generator families A–F + pinned) — PASSES
- [x] `t_m34_corpus_cases` in `tests/test_dsolve.c` — PASSES (runs before the pre-existing slow test)
- [x] version 0.131 → 0.132; DSOLVE_PLAN.md M34 bullet; changelog `2026-09-07.md`
- [x] `make check-c99` — exit 0
- [x] DSolve stress suites (m5/m12/m14/m17/m18/m19/m20) — all PASS
- [ ] Full corpus regression sweep (`ctest -R dsolve_corpus_`) — RUNNING
- [ ] valgrind spot-check VoP path (1337/1350)

## Notes
- `dsolve_tests` (206) SIGALRMs on this machine — **PRE-EXISTING** (clean HEAD/M33 also SIGALRMs,
  exit 142, real 120.7s): the `alarm(120)` watchdog trips on a borderline suite dominated by the
  pre-existing >30s `t_rischnorman_enum_cap_no_crash` (a first-order Abel `DSolve` — code my changes
  do not touch). Both versions die at the same test. Not an M34 regression.

## Review
- M34 closed 9 genuine gaps (90→99). Residue 1360 is the honest ceiling (no CAS has a closed form).
