# DSolve M19 — confluent Whittaker / ₁F₁ recognizer (+ corpus re-baseline)

Plan: `/Users/user/.claude/plans/twinkly-chasing-hammock.md`
Method: confluent equation → normal form `z''=r z` with one finite double pole +
nonzero const at ∞ → Whittaker (κ,μ) → emit `e^{cz} z^μ ₁F₁[…]` (verifiable),
new branch in `specialform_reduced_basis` (`src/calculus/dsolve_specialform.c`).

## Step 0 — corpus re-baseline (stale report predates M17/M18)
- [x] full §2.1.2 re-run: TRUE baseline **424/1000 scalar (42.4%), 0 FAIL**, 576 UNEVAL
- [~] regenerate reports/2.1.2.{md,tsv}; record in STATUS.md (pending M19 re-run finish)

## Step 1 — Whittaker/₁F₁ recognizer ✅ DONE
- [x] specialform_whittaker_basis(): single (x−x0)^2 pole via squarefree trick; num/den deg 2
- [x] b0=g|x0, b1=g'|x0, b2=lead; c=2√(−b2), μ=√(1/4−b0), κ=b1/c (Qc-convention)
- [x] emit W1/W2 = Exp[-z/2] z^(1/2±μ) Hypergeometric1F1[1/2±μ−κ, 1±2μ, z] (→auto PFQ)
- [x] gate: decline IntegerQ[2μ]; keep symbolic μ; sf_num_ok in both P==0 and pre-pass
- [x] pre-pass Whittaker fallback @ wider ≤200 leaf budget (cheap sf_ct)

## Step 2 — SCOPE DECISION: P==0-only (correctness-first)
- [x] pre-pass (P≠0) Whittaker gained 13 (97/101/104/…) BUT regressed 5 (94/470/472/806/811):
      recovery factor Exp[-∫P/2] shares finite-pole base w/ Whittaker z^(1/2±μ) → composed
      candidate stacks same-base radical powers → verify/zero_test/1F1-numeric $IterationLimit
      → starves Frobenius series fallback (those cases' baseline PASS) → REGRESSION.
- [x] leaf-count / μ-complexity CANNOT separate gains from regressions (data-dependent).
- [x] DECISION: run Whittaker on y'-free (P==0) surface ONLY. 0 regressions; base self-verifies
      against reduced eq. P≠0 confluent family (~13 cases) = documented future work.
- [x] verified: 94/470/472/806/811 back to SERIES-PASS; 102/568 solve; 97/101/104 → baseline ABORT
      (UNEVAL, not regression — never solved at baseline)

## Step 3 — tests + gates + docs ✅ DONE
- [x] test_dsolve_m19_stress.c (nf/shifted/energy/prepass grids) + CMake — pass
- [x] t_m19_* units in test_dsolve.c (confluent/prepass/integer-2μ decline) — pass
- [x] all DSolve ctest suites (dsolve_tests + m5/m12/m14/m17/m18/m19) + check-c99 green
- [x] valgrind: whittaker path leak-flat (decline A==B, no growth under ×8)
- [x] DSOLVE_PLAN.md M19 + §1c row; calculus.md; changelog; docstring
- [~] STATUS.md re-baseline + M19 line; lower argv[3] gate baseline (pending re-run)

## Review

Landed **M19**: the confluent Whittaker/₁F₁ recogniser (`specialform_whittaker_basis`
in `dsolve_specialform.c`) + a **§2.1.2 re-baseline** (the scoreboard was stale, pre-
M17/M18). A y'-free `y''+Q y==0` whose `Q` has a single finite double pole + rank-1
irregular point at ∞ → verifiable `Exp[-z/2] z^(1/2±μ) ₁F₁[1/2±μ-κ, 1±2μ, z]`, self-
verified against the reduced equation (inert `WhittakerM/W` never emitted).

**Measured: baseline 424/1000 (42.4%, the true post-M17/M18 number), M19 427 (+3: 102,
568, 611), 0 FAIL, 0 regressions.** Gate 612→576. All DSolve ctest suites (dsolve_tests +
m5/m12/m14/m17/m18/m19 stress) + check-c99 green; Whittaker path valgrind leak-flat.

Key decision — **P==0-only scope (correctness over yield).** The P≠0 pre-pass Whittaker
(with the normal-form recovery factor) gained ~13 (97/101/104/…) but **regressed ~5**
(94/470/472/806/811): the recovery factor `Exp[-∫P/2]` shares the finite-pole base with
the Whittaker `z^(1/2±μ)`, so the composed candidate stacks same-base symbolic-radical
powers whose verify/`zero_test`/`HypergeometricPFQ`-numeric hits `$IterationLimit` and
starves the Frobenius series fallback those cases relied on. Leaf-count and μ-complexity
cannot separate gains from regressions (data-dependent evaluator fragility). Restricting
to the P==0 surface (no recovery, no stacking) keeps 0 regressions; the P≠0 confluent
family (~13 cases) is documented future work needing an evaluator-robustness fix.
