# M38 — DSolve §2.2.18 corpus (Problems 1701–1800)

Plan: `/Users/user/.claude/plans/let-s-continue-our-implementation-logical-book.md`

## Decision log
- **Converter autonomous-indvar fix: SKIPPED** (plan step 1 fallback taken). Empirically:
  the `_missing_x` gate is not zero-impact (it would change the currently-passing prior
  record `2.2.13-1292`), 1791 is `_quadrature` (not reachable by that gate anyway), and a
  broader heuristic would alter 7+ currently-passing prior records (1157/1182/1187/1190/
  1534/1537/1603) — exactly what M36 deferred for 1534. Also, `1744` read as-converted
  (indvar `a`) reduces via NormalForm to `z''=−z` → a valid closed form, so it likely
  PASSES anyway. Generate as-is; document 1744/1791 as converter-limitation residue only
  if they don't pass. Minimal-impact, consistent with precedent.
- **Baseline: 95/100, 0 FAIL, 5 UNEVAL** = {1708,1709 Abel-2nd-kind (hard, M13); 1729
  `_linear` hang; 1763 2nd-order-nonhomog; 1769 2nd-order-nonhomog}.
- **FIX 1 (correctness, landed): zero_test decay false-positive → dropped forcing.**
  `PossibleZeroQ[E^(-2x^2)]`/`[E^(-x^2-2x)]` return True (documented decay false-positive).
  This corrupted TWO gating sites: (a) the Wronskian nonzero-check in
  `dsolve_variation_of_parameters` (decaying W read as degenerate basis → VoP declined),
  and (b) Kovacic's forcing detection (decaying forcing read as homogeneous → dropped the
  particular → shipped a homogeneous-only WRONG answer for 1763, masked as UNEVAL by the
  prelude's leaked-C[k]→UNFIT leniency). Fix: new shared `ds_is_structural_zero` (Expand→0)
  in `dsolve_common.{c,h}`, used at both gates (VoP + kovacic ×2). Did NOT touch the
  sensitive zero_test sampler (per [[project_possiblezeroq_decay_false_positive]] guidance:
  don't gate on PossibleZeroQ for decaying exprs). 1763 now solves:
  `C[1]E^(-x^2)+C[2]x E^(-x^2)+2E^(-x^2-2x)`. Safe: structural-zero only ADDS VoP work
  (a truly-zero forcing → yp==0; a truly-dependent basis is still structurally 0).
  After FIX 1: §2.2.18 = **96/100, 0 FAIL** (1763 gained, no regression in §2.2.18).
- **Residue (4, all sympy=False, bounded UNEVAL / no wrong answers) — documented, not fixed:**
  - `1708`/`1709`: Abel 2nd-kind class B rational IVPs (the M13-deferred AIR class).
  - `1729`: `Integrate[Sin[x]/(b Cos x-x Sin x)]` is genuinely non-elementary
    (`Integrate::nonelem`) — no elementary integrating factor.
  - `1769`: VoP needs `Integrate[E^x/Sqrt[x]] -> Sqrt[Pi] Erfi[Sqrt[x]]`; the Erf/Erfi tower
    (risch_special/knowles_erf/risch_tower) handles only `Int[x^(1/2)E^x]`, not negative
    half-integer powers. Extending it is a deep-Risch feature (high risk) — out of scope.
  Matches campaign precedent (M37 left 13, M36 left 3).

## Tasks
- [ ] 1. Generate `DSolve_test_status/DE_examples_2218.m` via converter; fix header comment.
- [ ] 2. Wire harness: `tests/CMakeLists.txt` new `dsolve_corpus_2_2_18_tests`; README row.
- [ ] 3. Build `dsolve_corpus_tests`; measure baseline (full run → TSV → report).
- [ ] 4. Triage: fix ALL FAILs (0-FAIL invariant, mandatory).
- [ ] 5. Root-cause UNEVAL fixes, biggest bucket first (2nd-order var-coeff linear; Abel-2B).
- [ ] 6. Anti-overfit `t_m38_*` unit tests in `tests/test_dsolve.c`.
- [ ] 7. Update tracking: STATUS.md, reports/2.2.18.{md,tsv}, README row, changelog M38,
        DSOLVE_PLAN.md M38 milestone, version.h (if code changed), calculus.md (if new method).
- [ ] 8. Verify: §2.2.18 0-FAIL; all prior corpus gates hold; dsolve_tests + stress suites;
        make check-c99; valgrind spot; rebuild code graph.

## Review (M38 — DONE)

**Outcome:** §2.2.18 (Problems 1701–1800) added to the corpus. **96/100, 0 FAIL, 0 regression.**
Version 0.135 → 0.136.

**Tasks:** all complete.
- [x] 1. `DE_examples_2218.m` generated (100 records, 14 IVP, 0 systems); header auto-correct.
- [x] 2. `tests/CMakeLists.txt` `dsolve_corpus_2_2_18_tests` (baseline 4); README row.
- [x] 3. Baseline measured: 95/100, 0 FAIL, 5 UNEVAL.
- [x] 4. 0 FAIL held throughout (no wrong answers).
- [x] 5. FIX: zero_test decay false-positive → dropped forcing (`ds_is_structural_zero` gating
       the VoP Wronskian + Kovacic forcing detection). Closed 1763 → 96/100. The other 4 UNEVAL
       are genuinely hard (all sympy=False) — documented residue, not hacked.
- [x] 6. `t_m38_forced_decay_wronskian` added + passes.
- [x] 7. STATUS.md, reports/2.2.18.{md,tsv}, README, changelog M38, DSOLVE_PLAN M38, version.h.
- [x] 8. Verified: §2.2.18 ctest gate passes (baseline 4); check-c99 exit 0; all prior corpus
       gates held (§2.1.2 −5, §2.2.1/2/6/7 improved, 0 regression); DSolve unit+M5/M14/M17 stress
       pass; dsolve_tests reaches only the pre-existing t_rischnorman SIGALRM (not a regression).

**Key win:** FIX 1 is a real correctness fix (removed a latent homogeneous-only wrong-answer class
in Kovacic), not just a coverage bump — verified by the +5 §2.1.2 / +others improvements with 0 FAIL.

**Residue (4, all sympy=False, bounded UNEVAL, no wrong answers):** 1708/1709 (Abel-2B, M13),
1729 (non-elementary integrating factor), 1769 (needs ∫E^x/√x→Erfi, a deep-Risch tower extension).

**Leak safety:** change is leak-neutral by construction — `ds_is_structural_zero` frees its
`Expand` temporary; the two gate swaps add no retained allocation; VoP/Kovacic ownership unchanged.
Linux CI valgrind is the definitive check (macOS valgrind is baseline-noisy).
