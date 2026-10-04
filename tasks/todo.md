# DSolve campaign — land §2.2.37 (Problems 3601–3700)

Full plan: `~/.claude/plans/let-s-continue-our-implementation-magical-peacock.md`
Scope: land + one general-fix wave, 0-FAIL invariant, general fixes only.

## Phase A — acquire & convert
- [x] curl `Ch2.S1.SS37.htm` (browser UA) → scratch; 317 KB LaTeXML "3601 to 3700", ~100 records
- [x] convert → `DE_examples_2237.m`: 100 records (21 IVP, 0 systems)
- [x] spot-check: no leftover LaTeX, labels 3601–3700, 20 linear/10 sep/8 Bernoulli/…

## Phase B — audit & gate
- [x] `make check-corpus-indvar` green (4904 records, 38 files)
- [x] added `dsolve_corpus_2_2_37_tests` to `tests/CMakeLists.txt`

## Phase C — baseline measurement
- [x] baseline: **96 PASS, 4 UNEVAL, 0 FAIL, 0 crash, 0 timeout**
- [x] 4 gaps diagnosed: 3666(n=Pi), 3668(n=√3) Bernoulli; 3662 verify-artifact; 3650 Root-IVP

## Phase D — one general-fix wave
- [x] diagnosed from inside cascade: `DSolve`Bernoulli` declines symbolic/irrational exponent
- [x] general fix: per-term exponent + `Y^n→W` abstraction (`src/calculus/dsolve_bernoulli.c`)
- [x] re-measure §2.2.37: **98 PASS, 2 UNEVAL, 0 FAIL** (+2: 3666, 3668)
- [x] Bernoulli unit regression (n=2,3,1/2,-1,x-coeff) all PASS
- [~] full corpus regression (`ctest -R dsolve_corpus`) — RUNNING (on §2.1.2 now)
- [x] added held-out unit test `t_m64_bernoulli_irrational_exponent` (test_dsolve.c)

## Phase E — land
- [x] ratchet §2.2.37 baseline 100000→2; reports/2.2.37.{md,tsv} regenerated
- [x] STATUS.md section block + M64 wave-history; M64 milestone in DSOLVE_PLAN.md
- [x] changelog (docs/spec/changelog/2026-09-28.md); calculus.md Bernoulli row updated
- [x] bump src/version.h 0.266→0.267
- [ ] rebuild dsolve_tests + run (confirm new unit test passes); valgrind spot-check
- [ ] commit + tag v0.267 (when user asks)

## Review
Baseline 96→98/100 on one general fix (Bernoulli irrational exponent), 0 FAIL.
Residue 2 (3650 Root-IVP, 3662 verify-artifact) documented. Pending: full
regression result + unit-test rebuild + valgrind, then commit/tag.
