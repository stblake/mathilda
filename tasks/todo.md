# M28 — DSolve §2.2.8 corpus (Problems 701–800)

## Plan
Replicate the M27 corpus-milestone process for §2.2.8 (Edwards & Penney, Table 2.33,
100 scalar first-order ODEs). Generate corpus → baseline → root-cause-fix tractable
residues → gate at honest residue → docs. 0 FAIL invariant.

## Tasks
- [x] 1. Generate `DE_examples_228.m` (100 records, 701–800, 20 IVP, 0 systems; spot-checked faithful).
- [x] 2. Register gate; run baseline → **98/100 PASS, 0 FAIL, 0 crash**; 2 UNEVAL (757, 783). reports/2.2.8.{tsv,md} written.
- [x] 3. Triage: 757 = tractable cascade-hang (Bernoulli); 783 = honest residue (needs general y=_G(x,y') method; all 3 CAS solve it but via a method we lack).
- [x] 4. Root-cause fix: `dsolve_bernoulli.c` `bern_Y_nonalgebraic` fast-decline on transcendental-in-y → 757 solves via Linearizable. Re-run → **99/100, 0 FAIL**.
- [x] 5. Final gate baseline set to 1 in CMakeLists (with rationale comment).
- [x] 6. Added `t_m28_bernoulli_hang_trig_substitution` (757 solves+verifies, a-family, genuine-Bernoulli 752 guard). dsolve_tests: All passed.
- [x] 7. Docs: STATUS.md §2.2.8 block + M28 bullet; README row; DSOLVE_PLAN.md M28; changelog 2026-09-07.md.
- [~] 8. Verify: check-c99 ✅ exit 0; dsolve_tests ✅; §2.2.x + §2.1.2 corpus regression + stress suites — IN PROGRESS.

## Review
### Outcome
- §2.2.8 (Problems 701–800): **99/100 PASS, 0 FAIL, 0 crashes**. Sole residue 783 (quartic-in-y, `y=_G(x,y')`; future `SolvableForY` method).
- One root-cause engine fix (Bernoulli fast-decline on transcendental-in-y), reusing verified machinery + the existing reconstruction/verify gate → no wrong answers possible.
- 783 documented as honest residue (matches the M24–M27 "leave the research-grade residue" pattern).
### Files
- New: `DSolve_test_status/DE_examples_228.m`, `reports/2.2.8.{tsv,md}`.
- Edited: `src/calculus/dsolve_bernoulli.c`, `tests/CMakeLists.txt`, `tests/test_dsolve.c`, `DSolve_test_status/STATUS.md`, `DSolve_test_status/README.md`, `DSOLVE_PLAN.md`, `docs/spec/changelog/2026-09-07.md`.
