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
| 2026-09-08 (**M21** side-effect) | **436 / 1000** | **43.6%** | **564** | **+4, 0 FAIL, 0 regression.** Not a §2.1.2-targeted wave — the M21 §2.2.1 shared fixes (scalar-Solve IC fit + `NthAlgebraic` denominator-clearing) also close 4 §2.1.2 first-order cases, and the `dsFreeParams` verifier fix (numeric back-substitution was vacuous) surfaced **no** new FAIL. Gate baseline kept at **572** (margin for the flaky fork cluster; not lowered since §2.1.2 was not the focus). |

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

## Section 2.2.1 — "Problems 1 to 100" (Table 2.19, sorted by problem number)

Corpus: `DE_examples_221.m` — 100 records, **all scalar, 63 IVPs** (4 symbolic ICs
`y(a)=b`, 3 swapped-variable `x=x(y)`). Elementary first-order + simple 2nd-order
(quadrature / linear / separable / homogeneous / Riccati). **Zero overlap** with §2.1.2
(that corpus is the hard "SymPy-failed" residue with no ICs). This is the first corpus
section carrying **initial conditions**: the equation slot is the DSolve-native list
`{ode, ic...}` and the harness verifies the ODE residual **and every IC** (M21).
`ctest -R dsolve_corpus_2_2_1_tests` · gate baseline **4**.

| Date | Scalar solved | Solve % | Gap (non-PASS) | Notes |
|------|--------------:|--------:|---------------:|-------|
| 2026-09-08 (M21 baseline) | 86 / 100 | 86.0% | 14 | 0 FAIL. First IC-verifying run (prelude `dsFreeParams` `Heads->True` bug fixed → numeric verify is real). |
| 2026-09-08 (**M21**)      | **96 / 100** | **96.0%** | **4** | **+10, 0 FAIL, 0 regression.** Three solver fixes (below). Gate baseline **4**. |

**M21 fixes** (all in the scalar first-order cascade / fit substrate):
1. **Swapped-variable `A/y'==B`** (#98/99/100) — `NthAlgebraic` now clears a
   top-derivative-bearing denominator (`Numerator[Together[·]]`) and recurses on the
   cleared ODE, so the `x=x(y)` spelling solves. (Worked around an `Exponent` bug — it
   returns 0 when any funcapp is present — via `FreeQ`.)
2. **Transcendental-inverse IC fit** (#29/30/33/34/40/61) — `dsolve_fit_constants` uses
   Solve's **scalar** form for a single-condition/single-constant fit; only the scalar
   form applies inverse-function inversion (a constant inside `Sqrt`/`Log`/`^(3/2)`/Airy
   ratio now fits, where the list form bubbled back and leaked `C[1]`).
3. **`ConditionalExpression` principal-branch collapse** (#60) — a multivalued `Tan`
   inverse fits as a `ConditionalExpression` over an integer family; take the principal
   branch → `y'=1+y², y(0)=0` gives `Tan[x]`.

**Residue (4, all bounded UNEVAL, 0 wrong answers):**

| Case | ODE | Why |
|---|---|---|
| 2.2.1-35 | `y'=Log[1+y²], y(0)=0` | non-elementary `∫dy/Log[1+y²]`; DSolve doesn't spot the equilibrium `y≡0` (bounded timeout, matches Maple/Mma) |
| 2.2.1-47 | `y'=4(xy)^(1/3)` | homogeneous class-G; general solution is a degree-12 `Root` object the verify can't confirm |
| 2.2.1-48 | `y'=2x Sec[y]` | separable, but the `ArcSin[Cos[2](…)]` inversion is slow (>8 s prelude limit) and carries constant artifacts |
| 2.2.1-67 | `y'=6 e^{2x−y}, y(0)=0` | `Solve[E^y==Q, y]` reuses `C[1]` as the Log branch-index while `C[1]` is already the integration constant → collision (a `Solve` generated-constant bug; separate follow-up) |

Full per-case results: `reports/2.2.1.tsv`; bucketed report: `reports/2.2.1.md`.

**Verifier fix (affects all sections):** the corpus prelude's `dsFreeParams` used
`Cases[…, Heads->True]`, collecting operator heads (`Plus`, `Times`, `Tan`, `Sec`) as
"parameters" and substituting numbers for them — so every residual became
non-numericizable and scored a vacuous `UNK` (trusted). Removing `Heads->True` makes the
numeric back-substitution real for the first time. §2.1.2 re-verified: **0 new FAIL**.

---

## Section 2.2.2 — "Problems 101 to 200" (Table 2.19, sorted by problem number)

Corpus: `DE_examples_222.m` — 100 records, **all scalar, 9 IVPs**. Continuation of
§2.2.1 (same elementary Table 2.19): linear / separable / homogeneous (classes A/C/G) /
Bernoulli / exact / Riccati / d'Alembert, plus a handful of 2nd-order missing-x/missing-y.
**Zero overlap** with §2.1.2. `ctest -R dsolve_corpus_2_2_2_tests` · gate baseline **8**.

| Date | Scalar solved | Solve % | Gap (non-PASS) | Notes |
|------|--------------:|--------:|---------------:|-------|
| 2026-09-08 (baseline) | 82 / 100 | 82.0% | 18 | 0 FAIL. Converter-fix generation (below). |
| 2026-09-08 (**M22**)  | **92 / 100** | **92.0%** | **8** | **+10, 0 FAIL, 0 regression.** Homogeneous-correctness + Exact-transcendental + FOS-implicit waves. Gate baseline **8**. |

**Converter fix (`tools/latex_ode_to_mathilda.py`, benefits every section):**
`is_condition_row` matched `<main>·(…)` multiplication (`y²(y'x+y)`, `x(5−x)`) as a
`y(P)` initial condition and dropped 8 ODE rows (109/119/129/166/175–178). Anchored to
the LHS (must be *solely* `y(…)`/`y'(…)`). Regenerating §2.2.1 is byte-for-byte identical
→ no regression; +6 of the 8 immediately PASS.

**M22 solver waves** (all reuse the verified implicit first-integral substrate
`dsolve_run_implicit`; 0 FAIL by construction):
1. **Homogeneous correctness** (`dsolve_homogeneous.c`) — retired latent WRONG answers
   (117 `√(x²+y²)`, 112 `E^(y/x)`) and a hang (107). Reduced RHS now via `F(1,v)` not
   `F(x,v·x)` (radicals collapse under `x→1`, no spurious `x`); `homog_exp_log_invert`
   gated to the log-sum case; explicit bodies numerically verified (drop the spurious
   inverse-branch → implicit fallback); `$rad` placeholder leaks rejected.
2. **Exact transcendental** (`dsolve_exact.c`, `dsolve.c`) — implicit potential
   `F(x,y)==C[1]` for exact ODEs with transcendental `M,N` (140/141/142/182/195),
   wired right after explicit Exact so it preempts the downstream hang. A Linearizable
   Bernoulli-shape recursion gate (skip a reduced eqn with a transcendental of `u`,
   e.g. `E^(u+E^u)`) removes the pre-Exact hang on 141.
3. **FirstOrderSubstitution implicit** (`dsolve_fos.c`, `dsolve.c`) — inert-integral
   relation for `y'==f[a x+b y+c]` with arbitrary `f` (159), as Mathematica returns.
4. **Bernoulli mixed-radical gate** (`dsolve_bernoulli.c`) — decline `(x y)^p` (not the
   pure `B(x)y^n` form; its detector spun), letting Homogeneous own 107.

**Residue (8, all bounded UNEVAL, 0 wrong answers):**

| Case | ODE | Why |
|---|---|---|
| 2.2.2-133 | `2x Sin y Cos y y'=4x²+Sin²y` | `y=G(x,y')` trig; slow/undecidable |
| 2.2.2-160 | `y'+p x y=q x yⁿ` | Bernoulli with a **symbolic** exponent `n` (genuinely hard) |
| 2.2.2-165 | `y'=Sin[x−y]` | correct explicit form but slow (>8 s prelude limit); implicit would be cleaner (future) |
| 2.2.2-170 | `r y''=(1+y'²)^{3/2}` | elastica/catenary; reduced 1st-order ODE the cascade cannot close |
| 2.2.2-175/176 | `x'=3x(5−x), x(0)=8/2` | logistic IVP; solves interactively (~1 s) but flaky under the forked cold-cache 8 s limit |
| 2.2.2-177/178 | `x'=4x(7−x)`, `x'=7x(x−13)` | **pre-existing** general-solution hang on the autonomous quadratic (separate follow-up) |

Full per-case results: `reports/2.2.2.tsv`; bucketed report: `reports/2.2.2.md`.

---

## Section 2.2.3 — "Problems 201 to 300" (Table 2.19, sorted by problem number)

Corpus: `DE_examples_223.m` — 100 records, **all scalar, 35 IVPs**. Continuation of
§2.2.2 (same elementary Table 2.19), but skewed toward higher-order **constant-
coefficient linear** (2nd/3rd/4th order, homogeneous + forced, real/repeated/complex
roots), **Euler–Cauchy / Emden–Fowler**, and a handful of elementary first-order
(exact / separable / linear / homogeneous class A/C/G / Bernoulli / d'Alembert).
**Zero overlap** with §2.1.2. `ctest -R dsolve_corpus_2_2_3_tests` · gate baseline **1**.

| Date | Scalar solved | Solve % | Gap (non-PASS) | Notes |
|------|--------------:|--------:|---------------:|-------|
| 2026-09-08 (baseline) | 98 / 100 | 98.0% | 2 | 0 FAIL. Converter unchanged (byte-identical §2.2.1/§2.2.2 regen). Two gaps: 204 (radical exact hang), 232 (Emden–Fowler). |
| 2026-09-08 (**M23**)  | **99 / 100** | **99.0%** | **1** | **+1, 0 FAIL, 0 regression.** Exact radical-potential → implicit first integral. Gate baseline **1**. |

**M23 solver fix** (reuses the verified implicit first-integral substrate; 0 FAIL by
construction):
- **Exact radical potential → implicit** (`dsolve_exact.c`) — 204's potential
  `F = 6 x^(3/2) y^(4/3) − 10 x^(6/5) y^(3/2)` carries fractional powers of `y`, so the
  explicit `Solve[F==C[1], y]` did not terminate (it hung the whole solve). The explicit
  Exact entry is now gated to a **rational-in-`y`** potential (`ds_is_rational_in`), so a
  radical/transcendental potential falls through to the existing implicit entry
  `dsolve_exact_implicit_try`, which returns `F(x, y[x]) == C[1]` verbatim (as Maple and
  Mathematica do), verified by the implicit-function rule. The gate never demotes an
  invertible case: a rational-in-`y` potential (`t_exact_xayb`'s Laurent-in-`y` form,
  `2xy+1+x²y'==0`) still solves explicitly.

**Residue (1, bounded UNEVAL, 0 wrong answers):**

| Case | ODE | Why |
|---|---|---|
| 2.2.3-232 | `y y'' == 6 x^4` | Emden–Fowler `_with_linear_symmetries`; the (correct) scaling-symmetry reduction `r=y/x³, s=ln x` lands on the autonomous `r r''+5r r'+6r²==6`, whose first-order reduction `r p p'==6−5rp−6r²` is an **Abel equation of the 2nd kind** — non-elementary for the cascade and squarely in the deferred-M13 (Abel Invariant Rational) territory. Declines cleanly (no wrong answer). |

Full per-case results: `reports/2.2.3.tsv`; bucketed report: `reports/2.2.3.md`.

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
- **M20 (2026-09-08)** — first-order symmetry gap Stage 1: `PolynomialShiftSubstitution`
  (`dsolve_polyshift.c`), the radical `[F(x),G(x)]`-symmetry sub-cluster. §2.1.2
  432/1000 (+6). See DSOLVE_PLAN.md M20 and the §2.1.2 block above.
- **M21 (2026-09-08)** — §2.2.1 corpus (Problems 1–100), first IVP-carrying section.
  IC verification in the harness; 3 solver fixes. **96/100 (+10, 0 FAIL).** See §2.2.1 block.
- **M22 (2026-09-08)** — §2.2.2 corpus (Problems 101–200). Homogeneous-correctness /
  exact-transcendental / FOS-implicit waves. **92/100 (+10, 0 FAIL).** See §2.2.2 block.
- **M23 (2026-09-08)** — §2.2.3 corpus (Problems 201–300). Exact radical-potential →
  implicit first integral (`dsolve_exact.c`, `ds_is_rational_in` gate). **99/100 (+1,
  0 FAIL, 0 regression);** sole residue 232 (Emden–Fowler → Abel 2nd kind). Anti-overfit
  unit `t_m23_exact_radical`. See §2.2.3 block.
- **Next** — the **P≠0 confluent family** (~13 cases: 97/101/104 and kin) is the biggest
  Whittaker residue, pending an evaluator-robustness fix for the same-base symbolic-radical
  verify (`zero_test` / `HypergeometricPFQ`-numeric `$IterationLimit`). Then parabolic-
  cylinder/Hermite (`POLY_r ndeg=2`) and trig/hyperbolic-potential-with-`y'`. Also still open:
  M18 Stages 2/3 (μ(x,y'), μ(y,y') — BLOCKED on a non-elementary first-order solver),
  3rd/high-order operator factoring (needs non-`TimeConstrained` bounding), Abel Invariant
  Rational (deferred M13).
