# DSolve corpus status

Scoreboard for `DSolve` against the 12000.org *Solving ODEs* corpora. **Update
after every method-wave** (see `README.md` and `../DSOLVE_PLAN.md` M15+). The
`ctest` gate for each section is the *baseline* (argv[3] in
`tests/CMakeLists.txt`) — the checked-in non-PASS high-water mark; every wave
must lower it.

Verdicts: PASS = verified closed form · UNEVAL = declined/timeout/`{}` ·
FAIL = wrong branch (numeric back-substitution) · SKIP = system.

---

## Section 2.1.2 — "Problems not solved, but were solved by Maple and Mathematica"

Corpus: `DE_examples_2.m` — 1204 records (1000 scalar + 204 systems).
`ctest -R dsolve_corpus_2_1_2_tests` · gate baseline **576**.

| Date | Scalar solved | Solve % | Gap (non-PASS) | Notes |
|------|--------------:|--------:|---------------:|-------|
| 2026-09-06 (M15 baseline) | 385 / 1000 | 38.5% | 615 | 0 FAIL, 6 crashes. Infrastructure landed. |
| 2026-09-06 (crash fixes)  | 388 / 1000 | 38.8% | 612 | 0 FAIL, crashes 6→2 (both non-reproducing). |
| 2026-09-06 (M16)          | 396 / 1000 | 39.6% | 604 | +8 (cv_num_ok symbolic-param verify + Pöschl-Teller 2F1 recognizer). |
| 2026-09-07 (M18 Stage 1)  | ~403 / 1000 | ~40.3% | ~597 | +7 projected on reducible-μ targets, **0 FAIL** (μ(x,y) Cheb-Terrab & Roche 1999 + `TrigToExp[Coth]` fix). Full re-run was pending — see next row. |
| 2026-09-07 (**M19 re-baseline**) | **424 / 1000** | **42.4%** | **576** | First full re-run since M16 — M17 + M18 measured together. The 396→424 jump is M17/M18 (the reports were stale). Honest post-M17/M18 baseline. 0 FAIL. |
| 2026-09-07 (**M19**)      | **427 / 1000** | **42.7%** | **573** | **+3** (2.1.2-102, -568, -611), **0 FAIL, 0 regressions**. Confluent Whittaker/₁F₁ recogniser on the y'-free (P==0) surface. Gate baseline 612→**576** (573 measured non-PASS + 3 margin for intermittent fork-harness crashes). |
| 2026-09-07 (**M20**)      | **432 / 1000** | **43.2%** | **568** | **+6** (2.1.2-402, -371, -372, -376, -403, -424), **0 FAIL, 0 real regressions**. `PolynomialShiftSubstitution` (`dsolve_polyshift.c`): the radical `[F(x),G(x)]`-symmetry sub-cluster of the 1st-order symmetry gap, `u=φ(x)+c y` → separable → implicit first integral. (The one P→U, 2.1.2-879, is a **load-flaky timeout** — a 2nd-order Frobenius case that PASSes in 5.9 s in isolation, under the 8 s fork limit, and is untouched by polyshift; effective +6 → 433 on a clean run.) Gate baseline 576→**572** (568 measured non-PASS + margin for the flaky fork-timeout cluster 879/208/872/983). |

### Gap by bucket (baseline, ranked)

| Bucket | tot | PASS | UNEVAL | pass% | target method |
|---|--:|--:|--:|--:|---|
| 2nd_linear            | 414 | 203 | 211 | 49.0 | M16 landed Legendre-symbolic + Pöschl-Teller; residue = Heun / parabolic-cylinder / Gegenbauer-Möbius / power→Bessel |
| 3rd_high_linear       | 142 |  32 | 110 | 22.5 | OperatorFactor Beke / 2nd-order right factors |
| 2nd_reducible_mu      | 102 |  39 |  63 | 38.2 | **M18 Stage 1 (μ(x,y)) +7**; residue needs μ(x,y')/μ(y,y') Stages 2/3 (Lemma-3 Cases C–F) |
| 1st_Abel              |  68 |   4 |  64 |  5.9 | Abel Invariant Rational (AIR, revive M13) |
| 1st_with_symmetry     |  65 |  27 |  38 | 41.5 | targeted Lie symmetry ansätze |
| 3rd_high_reducible    |  34 |   7 |  27 | 20.6 | higher-order missing-x/y + μ reduction |
| 1st_solvable_for_yx   |  42 |  17 |  25 | 40.5 | general solve-for-y/x (extend Lagrange) |
| NONE                  |  21 |   0 |  21 |  0.0 | (unclassified — inspect) |
| rational_misc         |  12 |   2 |  10 | 16.7 | (mixed) |
| Liouville_2nd         |   9 |   2 |   7 | 22.2 | Liouville transformation gaps |
| 1st_Riccati           |   7 |   0 |   7 |  0.0 | Riccati linearization gaps |
| 1st_other             |   6 |   0 |   6 |  0.0 | (inspect) |
| Emden_Fowler          |  10 |   5 |   5 | 50.0 | Emden–Fowler / elliptic quadrature |
| 2nd_other, 1st_Bernoulli | 4 | 1 | 3 | — | (inspect) |
| **ortho_poly**        |  54 |  54 |   0 | 100 | ✅ already solved (Kovacic) |
| **elliptic**          |   7 |   7 |   0 | 100 | ✅ already solved |
| **Bessel_special**    |   2 |   2 |   0 | 100 | ✅ |
| **Lienard**           |   1 |   1 |   0 | 100 | ✅ |

Full per-case results: `reports/2.1.2.tsv`; bucketed report: `reports/2.1.2.md`.

### Crashes — all 6 FIXED (2026-09-06)

Two root-cause fixes; all 6 baseline crashes now run clean (298/876/879 → **SOLVED**,
208/269/587 → clean decline). Both fixes protect every caller, not just DSolve.

| Case | Was | Now | Root cause / fix |
|---|---|---|---|
| 2.1.2-298 | SIGSEGV | ✅ SOLVED | see fix (A) |
| 2.1.2-876 | SIGSEGV | ✅ SOLVED | see fix (A) |
| 2.1.2-879 | SIGSEGV | ✅ SOLVED | see fix (A) |
| 2.1.2-587 | SIGILL  | decline | see fix (A) (transcendental coeff) |
| 2.1.2-208 | flaky   | decline | (A)-adjacent; stable now |
| 2.1.2-269 | SIGABRT | decline | see fix (B) |

- **(A) `dsolve_linear_normalize` transcendental-coefficient gate**
  (`src/calculus/dsolve_common.c`, new `ds_is_rational_in`): the polynomial/rational
  normaliser was handing coefficients like `a(λe^{λx}−a e^{2λx})` or `(1+e^{t²/2})²`
  to `PolynomialGCD`/`Cancel`, whose FLINT rational-canonicaliser emitted a
  non-finite content and then recursed on `GCD(−∞, 1)` → stack overflow. Gated to
  rational-in-x coefficients only (its documented domain); transcendental-coeff
  equations aren't const-coeff/Euler anyway, so they now fall through to
  Kovacic/Frobenius (which SOLVE 298/876/879).
- **(B) `expr_to_mpolyq` zero-base-negative-power guard**
  (`src/poly/flint_bridge.c`): `base^(−k)` with `base` converting to the zero mpoly
  called `fmpz_mpoly_q_inv(0)`, which hard-ABORTs FLINT (SIGABRT). Now declines so
  the caller falls back.

Post-fix re-run: 388/1000 solved, **612 non-PASS**, and 2 remaining "crashes"
(2.1.2-208, 2.1.2-983) that **do not reproduce in isolation** (3 clean trials each:
208 declines, 983 is a 4th-order Euler that times out) — intermittent, fork-layout
-dependent, effectively declines. Left counted in the 612 baseline; chase if they
become reproducible.

---

## Section 2.1.3

Corpus: `DE_examples_3.m` — pending fetch/convert. Gate:
`dsolve_corpus_2_1_3_tests` (to be added).

| Date | Scalar solved | Solve % | Gap | Notes |
|------|--------------:|--------:|----:|-------|
| — | — | — | — | awaiting section3.html |

---

## Wave history

- **M15 (2026-09-06)** — infrastructure: converter (`tools/latex_ode_to_mathilda.py`),
  corpora, self-verifying fork-per-case harness, this dashboard. 2.1.2 baseline 615.
- **M15 crash hardening (2026-09-06)** — two root-cause fixes (normalize
  transcendental gate; `expr_to_mpolyq` zero-inverse guard); 298/876/879 now solve.
  2.1.2 baseline 615 → **612** (388 solved). All DSolve suites + check-c99 green.
- **M16 (2026-09-06)** — 2nd-order-linear change-of-variable + recognizer wave.
  (1) `cv_num_ok` now instantiates free parameters, so M14's transform verifies
  symbolic-degree solutions (symbolic-k Legendre-Cot now solves). (2) New
  Pöschl-Teller / trig-potential recognizer in `SpecialFunctionForm` →
  numerically-verified `Hypergeometric2F1`. 388→**396 solved**, gap 612→604,
  0 wrong answers; all DSolve suites + check-c99 green. Anti-overfit families
  `t_m16_legendre_symbolic`, `t_m16_poschl_teller`.
- **M17 (2026-09-07)** — 2nd-order-linear affine→Gauss ₂F₁ recognizer + Liouville
  normal-form pre-pass in `SpecialFunctionForm`; the largest bucket (2nd-order linear).
- **M18 Stage 1 (2026-09-07)** — reducible-μ integrating factor μ(x,y) (Cheb-Terrab &
  Roche 1999, `dsolve_ifactor.c`): Φ degree-≤2 poly in y', Case A closed-form μ / Case B
  linear-ν-ODE μ; before `SecondOrderSymmetry`, linearity-gated. +7 reducible-μ solves
  (0 FAIL). Also: symbolic-parameter numeric-verify gate; `TrigToExp[Coth]` sign-bug fix.
  Anti-overfit `test_dsolve_m18_stress.c`; units `t_m18_*`.
- **M19 (2026-09-07)** — confluent Whittaker/₁F₁ recogniser + **§2.1.2 re-baseline**. The
  scoreboard was stale (pre-M17/M18); the first full re-run puts the honest baseline at
  **424/1000 (42.4%), 0 FAIL** (the 396→424 jump is M17+M18). New
  `specialform_whittaker_basis()` in `dsolve_specialform.c`: a y'-free `y''+Q y==0` with a
  single finite double pole + rank-1 irregular point at ∞ → verifiable
  `Exp[-z/2] z^(1/2±μ) ₁F₁[1/2±μ-κ, 1±2μ, z]` (self-verified vs the reduced equation;
  inert `WhittakerM/W` never emitted). **424→427 (+3: 102, 568, 611), 0 FAIL, 0
  regressions.** Gate 612→**576**. Anti-overfit `test_dsolve_m19_stress.c`; units `t_m19_*`.
  *Scope note:* Whittaker runs on the **P==0 surface only** — the P≠0 pre-pass (recovery
  factor) gained ~13 but regressed ~5 (94/470/472/806/811): the recovery factor shares the
  finite-pole base with `z^(1/2±μ)`, stacking same-base radical powers whose verify hits
  `$IterationLimit` and starves the Frobenius fallback. Data-dependent, not gate-able by
  size → restricted to P==0 for 0 regressions.
- **Next** — the **P≠0 confluent family** (~13 cases: 97/101/104 and kin) is the biggest
  Whittaker residue, pending an evaluator-robustness fix for the same-base symbolic-radical
  verify (`zero_test` / `HypergeometricPFQ`-numeric `$IterationLimit`). Then parabolic-
  cylinder/Hermite (`POLY_r ndeg=2`) and trig/hyperbolic-potential-with-`y'`. Also still open:
  M18 Stages 2/3 (μ(x,y'), μ(y,y') — BLOCKED on a non-elementary first-order solver),
  3rd/high-order operator factoring (needs non-`TimeConstrained` bounding), Abel Invariant
  Rational (deferred M13).
