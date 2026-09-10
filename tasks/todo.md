# DSolve M33 — §2.2.13 corpus (Problems 1201–1300) to full coverage

Section: 100 records, 32 IVPs, 0 systems. 52 second-order (incl. 6 Euler–Cauchy
"Emden"), 48 first-order. Hard residue: Abel-2nd-type (M13-deferred), a few
`[F(x)*G(y),0]` symmetry forms, one dAlembert.

## Stage 1 — corpus + dashboard wiring
- [ ] Generate `DSolve_test_status/DE_examples_2213.m` (`--label 2.2.13`, from indexsubsection22.htm)
- [ ] Add `dsolve_corpus_2_2_13_tests` to `tests/CMakeLists.txt` (temp high baseline for dev)
- [ ] Build `dsolve_corpus_tests`, capture baseline TSV, generate `reports/2.2.13.{tsv,md}`

## Stage 2 — root-cause fixes (measured-biggest-first, 0 FAIL invariant)
- [ ] Fix every FAIL (wrong answer) — trace to shared substrate
- [ ] Fix tractable UNEVAL via existing methods' substrate
- [ ] Pinned units `t_m33_*` per fix in `tests/test_dsolve.c`
- [ ] Defer research-grade Abel residue as bounded declines (document in STATUS)

## Stage 3 — lower gate + no regression
- [ ] Lower ctest baseline (argv[3]) to final non-PASS count
- [ ] All existing dsolve corpus gates green; `make check-c99` green

## Stage 4 — docs
- [ ] DSOLVE_PLAN.md M33 block; STATUS.md §2.2.13 + wave line; README row
- [ ] docs/spec/changelog/2026-09-07.md `## DSolve M33 …`; version.h 0.130→0.131

## Review

**Result: §2.2.13 baseline 91/100 → 99/100, 0 FAIL, 0 regression.** Sole residue: 1203
(Abel-2nd-kind class B, M13-deferred). Version 0.130 → 0.131.

Five shared-substrate root-cause fixes (each verified; each lifts earlier sections too):
1. `dsolve_exact.c` — potential built from the cleaner of `∫M dx`/`∫N dy` (Path 1/2, no
   hot-path Simplify → 1201); syntactic + Together denominator-clearing candidates for
   inexact rational forms (1233 negative-exp, 1216/1238); `mu(y)` tried when `mu(x)` yields
   nothing (1214).
2. `dsolve_common.c` — implicit verify `Together`s the residual before the zero-test
   (kills a Csc/Cot false-negative that rejected 1214).
3. `dsolve_separable.c` — `sep_find_split` numeric pre-filter fast-declines non-separable
   transcendental RHS (~15 s → ~0 s; unblocks 1201/1233 under the 8 s corpus timeout).
4. `dsolve_common.c` — `FIT_UNDECIDED`: a first-order IVP whose fit bubbles unevaluated
   declines so the cascade reaches Exact (closes exact/homogeneous overlap 1205/1231);
   gated `nfun==1 && max_order==1` (BVP free-constant + 2nd-order series IVPs untouched).
5. Regression caught+gated mid-wave: the ungated FIT_UNDECIDED declined correct 2nd-order
   IVPs in §2.2.5/§2.2.6/§2.2.11 — the order gate fixed it.

Verification: `dsolve_corpus_2_2_13_tests` ctest passes at baseline 1; §2.2.1/§2.2.2/…/§2.2.12
+ §2.1.2 gates held at baseline (§2.2.2, §2.2.7 incidentally improved); `t_m33_*` units all
True; `make check-c99` green; clean rebuild, no debug/backup artifacts.
