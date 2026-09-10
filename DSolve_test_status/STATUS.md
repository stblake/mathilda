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
`ctest -R dsolve_corpus_2_1_2_tests` · gate baseline **655** (M27: systems now verified).

| Date | Solved | Solve % | Gap (non-PASS) | Notes |
|------|-------:|--------:|---------------:|-------|
| 2026-09-06 (M15 baseline) | 385 / 1000 | 38.5% | 615 | 0 FAIL, 6 crashes. Infrastructure landed. |
| 2026-09-06 (crash fixes)  | 388 / 1000 | 38.8% | 612 | 0 FAIL, crashes 6→2 (both non-reproducing). |
| 2026-09-06 (M16)          | 396 / 1000 | 39.6% | 604 | +8 (cv_num_ok symbolic-param verify + Pöschl-Teller 2F1 recognizer). |
| 2026-09-07 (M18 Stage 1)  | ~403 / 1000 | ~40.3% | ~597 | +7 projected on reducible-μ targets, **0 FAIL** (μ(x,y) Cheb-Terrab & Roche 1999 + `TrigToExp[Coth]` fix). Full re-run was pending — see next row. |
| 2026-09-07 (**M19 re-baseline**) | **424 / 1000** | **42.4%** | **576** | First full re-run since M16 — M17 + M18 measured together. The 396→424 jump is M17/M18 (the reports were stale). Honest post-M17/M18 baseline. 0 FAIL. |
| 2026-09-07 (**M19**)      | **427 / 1000** | **42.7%** | **573** | **+3** (2.1.2-102, -568, -611), **0 FAIL, 0 regressions**. Confluent Whittaker/₁F₁ recogniser on the y'-free (P==0) surface. Gate baseline 612→**576** (573 measured non-PASS + 3 margin for intermittent fork-harness crashes). |
| 2026-09-07 (**M20**)      | **432 / 1000** | **43.2%** | **568** | **+6** (2.1.2-402, -371, -372, -376, -403, -424), **0 FAIL, 0 real regressions**. `PolynomialShiftSubstitution` (`dsolve_polyshift.c`): the radical `[F(x),G(x)]`-symmetry sub-cluster of the 1st-order symmetry gap, `u=φ(x)+c y` → separable → implicit first integral. (The one P→U, 2.1.2-879, is a **load-flaky timeout** — a 2nd-order Frobenius case that PASSes in 5.9 s in isolation, under the 8 s fork limit, and is untouched by polyshift; effective +6 → 433 on a clean run.) Gate baseline 576→**572** (568 measured non-PASS + margin for the flaky fork-timeout cluster 879/208/872/983). |
| 2026-09-08 (**M21** side-effect) | **436 / 1000** | **43.6%** | **564** | **+4, 0 FAIL, 0 regression.** Not a §2.1.2-targeted wave — the M21 §2.2.1 shared fixes (scalar-Solve IC fit + `NthAlgebraic` denominator-clearing) also close 4 §2.1.2 first-order cases, and the `dsFreeParams` verifier fix (numeric back-substitution was vacuous) surfaced **no** new FAIL. Gate baseline kept at **572** (margin for the flaky fork cluster; not lowered since §2.1.2 was not the focus). |
| 2026-09-09 (**M27**) | **446 / 1000 sc + 107 / 204 sys** | **44.6% sc / 52.5% sys** | **651** | **Systems now VERIFIED by back-substitution (was skipped): 553/1204 total, 0 FAIL.** Scalar solved **436 → 446 (+10)** — a side-effect of the M27 Solve periodicity-index fix (fresh mint index + integer-family collapse) closing inverse-function first-order cases. Systems scored for the first time: 107/204. Non-PASS 651 = 554 scalar + 97 systems (incl. 2 flaky fork-timeout crashes). Gate baseline **572 → 655** (651 + 4 margin). |

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
| 2026-09-09 (M31 re-baseline) | **98 / 100** | **98.0%** | **2** | +2, 0 FAIL. Side-effect of M31's shared prelude precision fix (large-cancellation residuals now verify). Gate baseline **4 → 2**. |

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
| 2026-09-09 (M31 re-baseline) | **93 / 100** | **93.0%** | **7** | +1, 0 FAIL. Side-effect of M31's shared prelude precision fix. Gate baseline **8 → 7**. |

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

## Section 2.2.4 — "Problems 301 to 400" (Table 2.19, sorted by problem number)

Corpus: `DE_examples_224.m` — 100 records, **all scalar, 24 IVPs**. Continuation of
§2.2.3 (same elementary Table 2.19), skewed toward higher-order **constant-coefficient
linear** (2nd/3rd/**5th** order, homogeneous + forced, real/repeated/complex roots),
**missing-x / missing-y** reductions, nonlinear **`_with_linear_symmetries`**, and a
handful of Euler / Emden–Fowler / exact / quadrature. **Zero overlap** with §2.1.2.
`ctest -R dsolve_corpus_2_2_4_tests` · gate baseline **1**.

| Date | Scalar solved | Solve % | Gap (non-PASS) | Notes |
|------|--------------:|--------:|---------------:|-------|
| 2026-09-08 (baseline) | 93 / 100 | 93.0% | 7 | Two **converter** fixes (below). 6 UNEVAL + **1 FAIL** (387 — a harness verify bug, not a wrong answer). |
| 2026-09-08 (**M24**)  | **99 / 100** | **99.0%** | **1** | **+6, FAIL→0, 0 regression.** One harness fix + two solver fixes + two converter fixes. Gate baseline **1**. |
| 2026-09-09 (M31 re-baseline) | **100 / 100** | **100.0%** | **0** | +1, 0 FAIL. Side-effect of M31's shared prelude precision fix. Gate baseline **1 → 0**. |

**Converter fixes (`tools/latex_ode_to_mathilda.py`; §2.2.1/2/3 regenerate byte-for-byte
identical → no regression):**
1. **Imaginary unit `i`** (309/310/311, `y''+2 i y'+3 y=0`, `y''=(-2+2 i√3)y`) — a
   constant-coefficient **complex** ODE with no explicit independent variable made
   `detect_symbols` pick the imaginary unit `i` as the indep var *and* keep it as a plain
   symbol. `i`/`I` are now excluded from indep-var candidates (like `e`) and a standalone
   `i` maps to Mathilda's `I`.
2. **`y^{(n)}` derivative notation** (336/340/343, 5th-order) — `y^{(5)}` was converted to
   `y[x]^((5))` (a **power** of y) instead of the 5th derivative; the mains substitution now
   recognises the parenthesized-order superscript and emits `y'''''[x]` (`Derivative[5]`).

**M24 fixes:**
1. **Harness verify variable-capture** (`dsolve_corpus_prelude.m`, benefits **every
   section**) — retired the sole FAIL (387, `m x''+k x==F0 Cos[om t]`). `dsResidVerdict`
   swept the residual over a loop variable `k`, which **collided with the ODE parameter
   `k`** (the spring constant): the loop forced `k` to 0…5 instead of its generic sample
   value, so a *correct* fitted solution produced a bogus nonzero residual → a FALSE
   "BAD". The sweep index / value holder are now `$`-prefixed (`$dsSweep`/`$dsVal`) — names
   the converter can never emit. Strictly more correct: it can only fix false FAILs.
2. **Trig-power/product forcing** (`dsolve_undetcoeff.c`, 326/362/363/365) —
   `UndeterminedCoefficients` now `TrigReduce`-linearises the forcing (`Sin[x]^2 →
   (1−Cos2x)/2`, `Cos[x]^3 → (3Cosx+Cos3x)/4`, `Sin[3x]Sin[x] → (Cos2x−Cos4x)/2`,
   `x Cos[x]^3 → (3x Cosx + x Cos3x)/4`) into first-harmonic sinusoids — each a UC
   function — so a trig power/product forcing solves tidily instead of declining to the
   (hanging) variation-of-parameters fallback. TrigReduce preserves value and never turns a
   UC function into a non-UC one, so it can only help.
3. **Numeric complex roots concretized** (`dsolve_common.c` `dsolve_homog_basis`, 312) —
   `y'''==y`'s complex cube roots were emitted as `Re[-(-1)^(1/3)]`/`Im[-(-1)^(1/3)]` (Re/Im
   do not auto-evaluate on a radical power), which blocked the IVP constant-fit. The basis
   builder now `ComplexExpand`s the real/imag parts of a **numeric** complex root (gated by
   `NumericQ`, so a symbolic-parameter root — where ComplexExpand could introduce Abs/Sign —
   is untouched), so `y'''==y, y(0)=1, y'(0)=0, y''(0)=0` fits to a concrete solution.

**Residue (1, bounded UNEVAL, 0 wrong answers):**

| Case | ODE | Why |
|---|---|---|
| 2.2.4-381 | `(x²−1)y''−2x y'+2y = x²−1` | variable-coefficient Legendre-type (`_with_linear_symmetries`, **SymPy-failed**); the homogeneous solves via Kovacic (`y1=x`) but the nonhomogeneous particular needs variation-of-parameters on a Kovacic basis with a rational forcing — a genuine new capability, not a reuse tweak. Declines cleanly (Maple solves it; SymPy does not). |

Full per-case results: `reports/2.2.4.tsv`; bucketed report: `reports/2.2.4.md`.

---

## Section 2.2.5 — "Problems 401 to 500" (Table 2.19, sorted by problem number)

Corpus: `DE_examples_225.m` — 100 records, **all scalar, 15 IVPs**. Continuation of
§2.2.4, but a **series-solution-heavy** chunk (Edwards & Penney, Ch. 8): 2nd-order
linear (constant- and variable-coefficient), Airy/Emden–Fowler, Gegenbauer (already
Kovacic), Bessel, Jacobi/₂F₁, one Liénard, one 3rd-order. **Zero overlap** with §2.1.2.
`ctest -R dsolve_corpus_2_2_5_tests` · gate baseline **1**.

| Date | Scalar solved | Solve % | Gap (non-PASS) | Notes |
|------|--------------:|--------:|---------------:|-------|
| 2026-09-08 (baseline) | 95 / 100 | 95.0% | 5 | 0 FAIL. Converter unchanged (byte-identical §2.2.1–4 regen). 5 UNEVAL: 428/482 (exact→Erf verify spin / Kovacic degenerate basis) + 459/463/490 (transcendental-coeff singular). |
| 2026-09-08 (**M25**)  | **99 / 100** | **99.0%** | **1** | **+4, 0 FAIL, 0 regression.** Three verified fixes (below). Gate baseline **1**. |

**M25 fixes** (all reuse verified machinery / return series the recurrence gates; 0 FAIL
by construction):
1. **Erf integrating-factor verify (`dsolve_common.c`, 428).** `y''+x y'+y==0` is exact
   → first-order linear `y'+x y==C[2]` → Erf closed form, but `dsolve_verify_body`'s
   `zero_test` spun on the Gaussian×Erf residual (`E^(-x²/2)·E^(x²/2)` products that never
   combine defeat the numeric precision ladder — `POSSIBLE_ZEROQ_IMPROVEMENTS.md` #1). The
   verify now `ExpandAll`-normalises an Erf/Erfi residual first (gated to Erf/Erfi so every
   other verify path is unchanged), so the exact Erf form returns.
2. **Kovacic fundamental-set guard (`dsolve_kovacic.c` + `dsolve_exactode.c`, 482).**
   `2x y''+(1-2x²)y'-4x y==0` has one Liouvillian solution `√x E^(x²/2)` (second is
   non-elementary); Kovacic's coincident-exponent path collapsed to the rank-deficient
   `(C[1]+C[2])√x E^(x²/2)` (it verifies but is not a general solution). A final
   independence guard rejects a degenerate basis → cascade falls through to Frobenius →
   correct two-parameter series. (ExactODE also now declines when its reduced sub-solve
   leaves a non-elementary `Integrate`, so 482 reaches Kovacic/Frobenius rather than hanging.)
3. **Transcendental-coefficient Frobenius (`dsolve_frobenius.c`, 463/490).** At a regular
   singular point, forming `xP = x·P` leaves a removable singularity when P,Q carry an
   analytic transcendental (`6 Sin[x]/x` is 6 at 0 but substitutes to `6 Sin[0]/0 =
   Indeterminate`), so the indicial roots came out garbage and Frobenius declined.
   `frobenius_regsing` now Taylor-normalises `xP`, `x²Q` first (a no-op for genuine
   polynomials), so analytic transcendental coefficients yield a verified Frobenius series.

**Residue (1, bounded UNEVAL, 0 wrong answers):**

| Case | ODE | Why |
|---|---|---|
| 2.2.5-459 | `x² y''+Cos[x] y'+x y==0` | **irregular** singular point at x=0 (`x·P = Cos[x]/x` is not analytic); the only analytic solution is a one-parameter formal power series (the second has an essential singularity). A transcendental-coefficient equation Mathilda leaves unevaluated, **matching Mathematica** (the shifted-Frobenius path deliberately declines transcendental coefficients rather than expand about an arbitrary point). Declines cleanly (no wrong answer). |

Full per-case results: `reports/2.2.5.tsv`; bucketed report: `reports/2.2.5.md`.

---

## Section 2.2.6 — "Problems 501 to 600" (Table 2.29, Edwards & Penney 6th ed.)

Corpus: `DE_examples_226.m` — 100 records, **74 scalar (47 IVP) + 26 systems**.
A forced-linear chunk: constant-coefficient
2nd/high-order IVPs with **general forcing f(t)** and **DiracDelta impulses**, plus
variable-coefficient series/Bessel/Emden–Fowler/Liénard and the special Riccati
`y'=x²+y²`. `ctest -R dsolve_corpus_2_2_6_tests` · gate baseline **5** (M27
re-baseline: systems are now verified, not skipped).

| Date | Solved | Solve % | Gap (non-PASS) | Notes |
|------|-------:|--------:|---------------:|-------|
| 2026-09-09 (baseline) | 57 / 74 sc | 77.0% | 17 | 0 FAIL. 17 UNEVAL: the forcing family 561–575 (general f + DiracDelta) all declined, plus 524 (slow Bessel) and 555 (singular-reduction hang, via timeout). Systems skipped. |
| 2026-09-09 (**M26**)  | **72 / 74 sc** | **97.3%** | **2** | **+15, 0 FAIL, 0 regression.** The whole forcing family now solves (below). |
| 2026-09-09 (**M27**)  | **95 / 100** | **95.0%** | **5** | Systems now VERIFIED by back-substitution (was skipped): **72/74 scalar + 23/26 systems**, 0 FAIL. Gate baseline **2 → 5** (the 3 new system residues exceed the 8 s budget). |

**M26 fixes** (Green's-function variation of parameters, no Laplace transform):
1. **Definite-integral variation of parameters (`dsolve_common.c`).** When the
   indefinite VoP integral does not close (arbitrary or impulse forcing), emit the
   causal Duhamel convolution `x_p = Integrate[Σ_i basis_i(t)·cof_i(s)·g(s)/(a_n W(s)),
   {s,0,t}]` over a fresh dummy. The kernel vanishes on the diagonal, so a zero-IC IVP
   fits its constants to 0. `TrigReduce`+`Expand` normalise the kernel so a resonant
   cos/sin forcing closes to the clean `t Sin[w t]` form and an exponential kernel
   integrates termwise. Covers 561–563, 572–575 (arbitrary f).
2. **DiracDelta sifting under a definite integral (`integrate_dirac.c`, new).**
   `∫ DiracDelta[αx+β] h(x) dx over [lo,hi] → h(x0)/|α|·B`, with a HeavisideTheta step
   for a symbolic upper limit (causal convention, full step at the base point) and a
   1 / 0 / ½ factor for numeric limits; a mixed integrand (1+δ, t+δ, δ+cos) is split by
   linearity. Covers 564–571 (impulse forcing) and fixes standalone `Integrate[δ·f]`.
3. **HeavisideTheta / DiracDelta value + derivative rules (`distributions.m`,
   `deriv.c`).** `H(0)=0`, `H(x>0)=1`, `H(x<0)=0`, `δ(x≠0)=0`, and `H' = DiracDelta`,
   so an IC fit at the base point resolves a shifted step/impulse and the causal
   response satisfies its pre-impulse conditions. `dsolve_constcoeff.c` gates the
   "homogeneous" branch off a DiracDelta forcing (it samples numerically to 0).
4. **Undefined-function guards (`dsolve_common.c`, `integrate.c`).** Variation of
   parameters skips the indefinite attempt for an arbitrary/undefined forcing (it never
   closes and can hang), and `integrate_definite` skips the improper/parametric methods
   (residue/Ramanujan/differentiation-under-the-integral) on an undefined-function
   integrand — they cannot close a `K(t,s) f(s)` convolution and spent seconds churning
   (the exp-kernel case dropped 6 s → 0.3 s). `dsolve_verify_body` keeps (never rejects)
   a distributional residual (definite Integrate / DiracDelta / HeavisideTheta) rather
   than driving `zero_test` into a spin.

**Residue (2, bounded UNEVAL, 0 wrong answers):**

| Case | ODE | Why |
|---|---|---|
| 2.2.6-524 | `y''+x⁴ y==0` | Emden–Fowler; the correct closed form `√x BesselJ[1/6, x³/3] + √x BesselY[1/6, x³/3]` is returned but takes longer than the 8 s per-case `TimeConstrained` DSolve budget (a performance residue, not a gap). |
| 2.2.6-555 | `t x''+(t-2)x'+x==0`, x(0)=0 | Exact → the regular-singular first-order reduction `t x'+(t-3)x==C[2]` whose integrating-factor quadrature `∫E^t/t⁴` is non-elementary (`ExpIntegralEi`); a series residue. Declines via timeout (no wrong answer). |

Full per-case results: `reports/2.2.6.tsv`; bucketed report: `reports/2.2.6.md`.

---

## Section 2.2.7 — "Problems 601 to 700" (Table 2.31)

Corpus: `DE_examples_227.m` — 100 records, **50 scalar (23 IVP) + 50 systems** — the
first section that is **half systems**, and the wave (**M27**) that taught the corpus
harness to **verify systems** by back-substituting the multi-function solution (they
were skipped through M26). Scalar half is elementary first-order (25 separable, 14
quadrature, 8 linear, 2 homogeneous class-G, 1 Riccati); systems are 25 2-D / 18 3-D /
7 4-D, 47 constant-coefficient + 3 variable-coefficient.
`ctest -R dsolve_corpus_2_2_7_tests` · gate baseline **7**.

| Date | Solved | Solve % | Gap (non-PASS) | Notes |
|------|-------:|--------:|---------------:|-------|
| 2026-09-09 (baseline) | 90 / 100 | 90.0% | 10 | Before the two engine fixes: 0 FAIL. Systems verified for the first time (48/50 scalar + 42/50 systems). |
| 2026-09-09 (**M27**)  | **93 / 100** | **93.0%** | **7** | **+3, 0 FAIL, 0 regression.** 48/50 scalar + 45/50 systems. Two engine fixes below took 90→93. |
| 2026-09-09 (M31 re-baseline) | **94 / 100** | **94.0%** | **6** | +1, 0 FAIL (48/50 scalar + 46/50 systems). Side-effect of M31's shared prelude precision fix. Gate baseline **7 → 6**. |

**M27 fixes** (system verification + two root-cause engine bugs the corpus surfaced):
1. **Corpus harness verifies systems (`dsolve_corpus_prelude.m`).** `dsExplicitQ`
   generalised to a List function slot, and the system-skip removed: a system branch
   `{x->Function[…], y->Function[…], …}` is back-substituted per equation exactly like a
   scalar (the numeric machinery already substituted a whole rule-list). §2.1.2 (204
   systems) and §2.2.6 (26 systems) were re-baselined; §2.2.1–§2.2.5 (pure scalar) are
   unchanged.
2. **Separated-exponent Simplify hang (`dsolve_linsys.c`).** A constant-coefficient
   system as ordinary as `x'=-50x+20y, y'=100x-60y` (eigenvalues -10, -100) HUNG:
   `dsolve_linsys_tidy` Simplify-ed a body that is a sum of exponentials with
   widely-separated real decay rates, and Simplify's zero-test spins numericising
   `E^(-10 t)` against `E^(-100 t)` (catastrophic dynamic range). tidy now Expands any
   exponential body (it already did for complex/large ones); +2 systems (636, 650).
3. **Solve periodicity-index collision (`solveinv.c` + `dsolve_common.c`).** `y'=2x Sec[y]`
   returned a WRONG answer (masked as UNEVAL by the leaked→UNFIT rule): `DSolve\`Separable`
   feeds Solve an equation already carrying the integration constant `C[1]`, and Solve
   reused `C[1]` as the `2πk` inverse-trig periodicity index; once its
   `Element[C[1],Integers]` constraint was stripped the shared `C[1]` corrupted the
   solution at non-integer values. `solveinv` now **seeds its mint counter past every
   `C[k]` already in the equation** (fresh index), and `dsolve_extract_solutions` collapses
   the integer family (`Element[C[k],Integers]` only — not a range condition, which would
   zero the integration constant and broke `y'=3x²(1+y²), y(0)=1`) to its principal branch.
   +1 (684).

**Residue (7, bounded UNEVAL, 0 wrong answers):**

| Case | ODE / system | Why |
|---|---|---|
| 2.2.7-603/606/607 | constant-coeff 2-D/3-D systems, forced | Irrational/complex eigenvalues + polynomial/exponential forcing: the correct closed form exceeds the 8 s per-case DSolve budget (a performance residue). |
| 2.2.7-604/608 | variable-coefficient non-triangular systems | Genuinely coupled, variable-coefficient (`t x`, `E^t y` entries) — the honest engine gap (`LinearSystemVarCoeff` covers only the scalar-factor `A(x)=f(x)B` class). |
| 2.2.7-675 | `y'=Log[1+y²], y(0)=0` | Non-elementary separable (plus the missed equilibrium `y≡0`) — matches Mathematica. |
| 2.2.7-683 | `y'=4(x y)^(1/3)` | Homogeneous class-G; the implicit inversion path exceeds the 8 s budget (a performance residue). |

Full per-case results: `reports/2.2.7.tsv`; bucketed report: `reports/2.2.7.md`.

---

## Section 2.2.8 — "Problems 701 to 800" (Table 2.33, Edwards & Penney)

Corpus: `DE_examples_228.m` — 100 records, **100 scalar (20 IVP), 0 systems** — a return
to elementary first-order after the half-systems §2.2.7. Same territory as §2.2.1–§2.2.4:
23 linear, 16 separable, ~25 homogeneous (class A/G/C), 11 exact (incl. two fractional-
power potentials), 19 Bernoulli (several fractional-exponent), 3 quadrature, plus a few
`y'=F(ax+by+c)` substitution / Riccati forms.
`ctest -R dsolve_corpus_2_2_8_tests` · gate baseline **1**.

| Date | Solved | Solve % | Gap (non-PASS) | Notes |
|------|-------:|--------:|---------------:|-------|
| 2026-09-09 (baseline) | 98 / 100 | 98.0% | 2 | 0 FAIL, 0 crashes. Both residues in the `y=_G(x,y')` class. |
| 2026-09-09 (**M28**)  | **99 / 100** | **99.0%** | **1** | **+1, 0 FAIL, 0 regression.** One Bernoulli cascade-hang fix below took 98→99. |

**M28 fix** (one root-cause engine bug the corpus surfaced):
1. **Bernoulli hangs the cascade on a transcendental-in-`y` RHS (`dsolve_bernoulli.c`).**
   `2.2.8-757` `2 x Sin[y]Cos[y] y' == 4 x² + Sin[y]²` reduces (`u = Sin[y]²`) to the linear
   `u' − u/x == 4 x`, which `Linearizable` solves in one step. But solved for `y'` the RHS is
   a rational function of `Sin[y]`/`Cos[y]` — transcendental in `y`, not the Bernoulli form
   `A(x) y + B(x) y^n` — and the Bernoulli exponent detector (`Y F_Y`, `Cancel`, `ds_free_of`)
   spun on it for seconds, timing out the whole cascade (`$Aborted`) before `Linearizable` was
   reached. `dsolve_bernoulli.c` now declines immediately when `y` appears inside a non-`Power`
   function head (`Sin[y]`, `Exp[y]`, …) or in a power exponent — mirroring the existing
   `bern_mixed_radical` early-decline guard; a genuine Bernoulli (algebraic in `y`) is
   untouched. +1 (757).

**Residue (1, bounded UNEVAL, 0 wrong answers):**

| Case | ODE | Why |
|---|---|---|
| 2.2.8-783 | `y'=1+x²+y²+x²y⁴` | Quartic in `y` (beyond Abel). Maple/Mma/SymPy solve it only via the general "solve-for-`y` then differentiate" (Maple's `y=_G(x,y')`) method Mathilda lacks — declines cleanly. A future `SolvableForY` method would close it. |

Full per-case results: `reports/2.2.8.tsv`; bucketed report: `reports/2.2.8.md`.

---

## Section 2.2.9 — "Problems 801 to 900" (Edwards & Penney)

Corpus: `DE_examples_229.m` — 100 records, **100 scalar (35 IVP), 0 systems**. Dominated by
second-order linear: 38 constant-coefficient homogeneous, 35 constant-coefficient
nonhomogeneous (undetermined coefficients / variation of parameters), 11 Euler/Emden–Fowler,
plus 7 with `x(t)` as the dependent variable (`x''` notation, 862–868), 3 complex-coefficient
(`i` in the equation: 857/858/859), and 6 first-order (separable / homogeneous / Bernoulli /
Abel). Mathilda's strongest DSolve territory.
`ctest -R dsolve_corpus_2_2_9_tests` · gate baseline **0**.

| Date | Solved | Solve % | Gap (non-PASS) | Notes |
|------|-------:|--------:|---------------:|-------|
| 2026-09-09 (baseline) | 100 / 100 | 100.0% | 0 | 0 FAIL, 0 crashes. Fully covered out of the box by the existing const-coeff / Euler / VoP / first-order specialists. |
| 2026-09-09 (**M29**)  | **100 / 100** | **100.0%** | **0** | **0 FAIL, 0 regression.** No ODE-solver fix needed; the wave added piecewise rounding-function derivatives (below), upgrading 898's verification. |

**M29 feature** (a general engine addition the corpus made visible, not an ODE-solver fix):
1. **Piecewise derivatives of the integer-rounding functions (`src/calculus/deriv.c`).**
   `D` of `Floor`/`Ceiling`/`Round`/`IntegerPart`/`FractionalPart` now returns the
   Mathematica `Piecewise[{{v, cond}}, Indeterminate]` form (0 off the jump set, 1 for
   `FractionalPart`), composing with the chain rule. Problem `2.2.9-898`
   (`y''+9y == 2 Sec[3x]`) is solved by variation of parameters and its solution carries a
   `Floor` branch-tracking term; previously `D[Floor[u],x]` was the inert
   `Derivative[1][Floor][u]`, so the ODE residual never numericized and the harness passed
   898 only under the "non-numericizable ⇒ trust DSolve" path. The residual now reduces to a
   genuine numeric ~0 — 898 is verified, not merely trusted.

**Residue: none** (0 UNEVAL, 0 FAIL).

Full per-case results: `reports/2.2.9.tsv`; bucketed report: `reports/2.2.9.md`.

---

## Section 2.2.10 — "Problems 901 to 1000" (Edwards & Penney)

Corpus: `DE_examples_2210.m` — 100 records, **56 scalar (15 IVP) + 44 systems**. The most
systems-heavy section: constant-coefficient linear (2nd/3rd/4th order, homogeneous + forced,
real / repeated / complex roots), Euler–Cauchy (2nd & 3rd order, incl. a Bessel and two
Gegenbauer/Legendre-type equations), and **44 first-order linear systems** (2×2 / 3×3 / 4×4
constant matrices — 42 homogeneous, 2 with polynomial/exponential forcing).
`ctest -R dsolve_corpus_2_2_10_tests` · gate baseline **0**.

| Date | Solved | Solve % | Gap (non-PASS) | Notes |
|------|-------:|--------:|---------------:|-------|
| 2026-09-09 (baseline) | 98 / 100 | 98.0% | 2 | 0 FAIL, 0 crashes, 0 timeouts. Two UNEVAL: 924 (forced 2×2 system, irrational spectrum) and 907 (inhomogeneous Legendre-type var-coeff). |
| 2026-09-09 (**M30**)  | **100 / 100** | **100.0%** | **0** | **+2, 0 FAIL, 0 regression.** Two root-cause fixes below (both in the linear-ODE machinery). 56/56 scalar + 44/44 systems. |

**M30 fixes** (two root-cause bugs the section surfaced, both about *forced* linear ODEs):
1. **Forced constant-coefficient systems with an irrational spectrum (`src/calculus/dsolve_linsys.c`).**
   `2.2.10-924` (`x'=2x+4y+3eᵗ, y'=5x-y-t²`, eigenvalues `(1±√89)/2`) ran **>90 s**: the
   variation-of-parameters integral `∫e^{−λt}tᵐ dt` made `Integrate` rationalise the `1/λᵏ`
   coefficient of an irrational `λ` into a *hundreds-of-digit* integer and spin. The assembler
   now abstracts a **real-irrational** eigenvalue to a fresh symbol before the integral (and
   `Simplify`s the algebraic coefficients — safe once the exponents are symbolic, so the
   widely-separated-decay-rate `Simplify` hang cannot apply), then substitutes the eigenvalue
   back — **0.6 s**. **Rational** eigenvalues stay concrete so genuine resonance is still
   handled; **complex** eigenvalues stay concrete too (their `e^{at}Cos/Sin[bt]` real form is
   the existing `ComplexExpand` path — this is what keeps the complex-spectrum 4×4 `2.2.10-927`
   fast).
2. **Kovacic now closes an *inhomogeneous* variable-coefficient ODE (`src/calculus/dsolve_kovacic.c`).**
   `2.2.10-907` (`(x²−1)y″−2xy′+2y == x²−1`, a Legendre-type equation) declined: the Kovacic
   solver handled only the homogeneous equation. It now accepts a forcing
   (`dsolve_second_order_PQ_forced`), **de-obfuscates** the fundamental set it recovers
   (`Sqrt[−1+x²]·E^(−½Log[1+x]+3⁄2Log[−1+x])` is really `(x−1)²`; convert `Exp[c Log u]→uᶜ`,
   split radicands, `PowerExpand`), and adds the particular solution by **variation of
   parameters** over the cleaned basis, re-verifying the full solution numerically.

**Residue: none** (0 UNEVAL, 0 FAIL).

Full per-case results: `reports/2.2.10.tsv`; bucketed report: `reports/2.2.10.md`.

---

## Section 2.2.11 — "Problems 1001 to 1100" (Edwards & Penney)

Corpus: `DE_examples_2211.m` — 100 records, **59 scalar (13 IVP) + 41 systems**. A large leading
block of **41 first-order constant-coefficient linear systems** (1001–1041, 2×2 … 6×6 constant
matrices — defective / repeated / complex spectra), then scalar first-order quadrature /
separable / linear ("class A") and second-order linear — constant-coefficient (`missing_x`),
exact, Euler, Gegenbauer, Emden–Fowler, Liénard, Airy. The scalar half is solved entirely out of
the box (all 59/59); both gaps were systems.
`ctest -R dsolve_corpus_2_2_11_tests` · gate baseline **0**.

| Date | Solved | Solve % | Gap (non-PASS) | Notes |
|------|-------:|--------:|---------------:|-------|
| 2026-09-09 (baseline) | 98 / 100 | 98.0% | 2 | 0 FAIL, 0 crashes, 0 timeouts. Two UNEVAL, both systems: 1014 (coupled 3×3 DAG) and 1001 (4×4, spectrum {16,32,48,64}). |
| 2026-09-09 (**M31**)  | **100 / 100** | **100.0%** | **0** | **+2, 0 FAIL, 0 regression.** Two root-cause fixes below. 59/59 scalar + 41/41 systems. |

**M31 fixes** (the two systems the section surfaced):
1. **`Integrate` linearity over a distributed product (`src/calculus/integrate.c`, `try_linearity`).**
   `2.2.11-1014` (`{x1'=2x1, x2'=−7x1+9x2+7x3, x3'=2x3}`) is a DAG solved by `TriangularSystem`,
   which peels the sources `x1,x3` and asks the scalar engine to integrate the integrating-factor
   integrand `Integrate[e^{−9x}(7 C[k]e^{2x} − 7 C[j]e^{2x}), x]` — a product of an exponential
   with a **sum** of exponentials (`Times[c, Plus[…]]`, exponents uncombined). `Integrate` took
   its exponential-substitution path and returned, in **55 s**, a **branch-wrong** antiderivative
   (a spurious `(−1)^{1/9}` factor; it survived `dsolve_run_system`'s verify because the residual
   is zero-test-*undecidable*). The **root fix** is in `Integrate` itself: `try_linearity` now
   distributes a product over a sum factor — `Integrate` is linear, so `c(g+h) → cg+ch`,
   integrated term-by-term and committed only if every term closes elementary (otherwise the
   whole-integrand cascade still runs, so a sum that is elementary only as a whole is not lost).
   The exponentials then collapse (`e^{−9x}e^{2x}→e^{−7x}`, the clean path) — fast and correct.
   This also repairs the **direct** user-reported bug:
   `Integrate[e^{−9x}(a e^{2x} − b e^{2x}), x]` was `−(a−b)/7 · (e^{2x})^{−7/2}` in ~9 s (ugly, and
   genuinely branch-wrong for symbolic/funcapp coefficients), now `−(a−b)/7 · e^{−7x}` in ~4 ms.
   Guarded by `test_linearity_distributes_product` (`tests/test_integrate_dispatch.c`).
2. **Corpus verifier made cancellation-robust (`dsolve_corpus_prelude.m`, `dsResidVerdict`).**
   `2.2.11-1001` (4×4, eigenvalues `{16,32,48,64}`) solves **correctly** (`Simplify[residual]≡0`),
   but back-substitutes to a difference of `e^{64x}`-scale terms that, at the prelude's 20-digit
   numeric sweep over `x≈1.1…3`, looks large (catastrophic cancellation: `|resid|@20 ≈ 8.9·10⁴³`
   at `x=3`, but `@120 ≈ 2·10⁻⁵⁶`) → a false "BAD" → `UNFIT` → UNEVAL. The shared verifier now
   **re-evaluates a not-small residual sample at 200-digit precision**: a genuine nonzero stays
   nonzero, a cancellation artifact collapses to ~0. The change is **monotone** — it can only
   turn a spurious "not small" into "small", so it can lower a section's non-PASS count, never
   raise it, and never introduces a FAIL.

**Residue: none** (0 UNEVAL, 0 FAIL).

Full per-case results: `reports/2.2.11.tsv`; bucketed report: `reports/2.2.11.md`.

---

## Section 2.2.12 — "Problems 1101 to 1200" (Edwards & Penney)

Corpus: `DE_examples_2212.m` — 100 records, **100 scalar (38 IVP), 0 systems**. Elementary
first-order territory (36 separable / 18 linear / 15 linear "class A" / 13 quadrature /
Bernoulli / exact / homogeneous), a run of "Abel 2nd type / class A" that are really
homogeneous-degree-0 rational (solved as homogeneous), plus two solvable-for-y/x forms.
`ctest -R dsolve_corpus_2_2_12_tests` · gate baseline **3**.

| Date | Solved | Solve % | Gap (non-PASS) | Notes |
|------|-------:|--------:|---------------:|-------|
| 2026-09-09 (baseline) | 80 / 100 | 80.0% | 20 | **1 FAIL** (1147, a wrong answer), 19 UNEVAL. Not solved out of the box — the section exposes real defects. |
| 2026-09-09 (**M32**)  | **97 / 100** | **97.0%** | **3** | **+17, FAIL→0, 0 crashes, 0 regression.** Three root-cause fixes below. |

**M32 fixes** (all shared-substrate, so they lift earlier sections too — none regress):
1. **IVP constant-fitting (`src/calculus/dsolve_common.c`).** The fitter now DROPS a branch that
   an initial condition cannot be met on: it substituted `Undefined` (Solve had no consistent
   constant — `2.2.12-1147`, the section's sole **FAIL/wrong answer**), or it is a scalar
   unsatisfiable inverse branch (the wrong `±` / `Root` index) — but ONLY when a sibling branch
   actually fits, so a lone basis singularity is still kept and a legitimately UNDER-determined BVP
   (`y''+y==0, y[0]==0, y[π]==0 → C[2] Sin[x]`) keeps its free constant. Bernoulli now emits BOTH
   real signs for an even `1−n` root (`y'==(1−2x)/y, y[1]==−2` needs `−√`). This fixed the FAIL plus
   11 UNEVAL separable IVPs (`±`-branch selection, sign flip, and the Root-form cubic separables
   1149/1150, whose constant is fitted on the implicit first integral `G(x0,y0)` with no inversion).
2. **Separable recognizer + implicit twin (`src/calculus/dsolve_separable.c`).** Accept
   generic-parameter splits (`(a y+b)/(c y+d)`), and add `dsolve_separable_implicit_try` (mirroring
   Exact/ExactImplicit) that returns the first integral `∫dy/g == ∫f dx + C[1]` — keeping a
   non-elementary integral **unevaluated** — when the relation does not invert for y. Solves
   `Cot[t]y/(1+y)`, `Cos²x Cos²2y`, and the autonomous non-elementary `−2 ArcTan[y]/(1+y²)` (1186).
3. **Integrate Gaussian→Erf / Ei / PolyLog recognizer variable (`src/calculus/risch_special.c`).**
   The completing-the-square templates emitted the antiderivative in a **literal `x` for every
   integration variable** (`Integrate[E^(a^2/2), a]` was `… Erf[… x …]`), so any Bernoulli/linear
   solve over a non-`x` independent variable produced a stray-variable — hence non-verifying —
   answer. Now the real variable is threaded through the template. Fixed 1182/1190.

**Residue (3, research-grade — NOT wrong answers, all bounded declines):**
- `1135` `y'==(x−e^{−x})/(x+e^y)` (`y=_G(x,y')`) — solvable-for-y then differentiate; the induced
  `p`-ODE is transcendental and not cascade-solvable.
- `1200` `e^x sin y + 3y − (3x − e^x sin y)y' == 0` (`x=_G(y,y')`) — cannot even be solved for `x`
  (`e^x` + linear `x`).
- `1157` `y'==(a y+b)/(c y+d)` with `a` the independent variable — an Abel equation of the second
  kind (the M13-deferred `[_Abel]` class).

These three are the documented "solvable-for-y/x + Abel" gap (`DSOLVE_PLAN.md` M13 / M28 residue);
a general generalised-d'Alembert method was prototyped and confirmed to close none of them, so it is
deferred to its own milestone rather than half-built here.

Full per-case results: `reports/2.2.12.tsv`; bucketed report: `reports/2.2.12.md`.

---

## Section 2.2.13 — "Problems 1201 to 1300" (Edwards & Penney)

Corpus: `DE_examples_2213.m` — 100 records, **100 scalar (32 IVP), 0 systems**. A MIX of
first-order (exact / linear / separable / homogeneous / Abel / symmetry) and 2nd-order linear
(46 reducible-μ, 6 Euler–Cauchy "Emden–Fowler", plus const-coeff). Unlike the pure first-order
§2.2.8/§2.2.12, half the section is 2nd-order.
`ctest -R dsolve_corpus_2_2_13_tests` · gate baseline **1**.

| Date | Solved | Solve % | Gap (non-PASS) | Notes |
|------|-------:|--------:|---------------:|-------|
| 2026-09-09 (baseline) | 91 / 100 | 91.0% | 9 | 0 FAIL, 9 UNEVAL: exact-transcendental / rational-clear declines, a `mu=Sin y` hang, an exact/homogeneous IVP shadow, and 1 Abel-2nd-kind. |
| 2026-09-09 (**M33**)  | **99 / 100** | **99.0%** | **1** | **+8, FAIL→0, 0 regression.** Five root-cause fixes below. |

**M33 fixes** (all shared-substrate, so they lift earlier sections too — none regress):
1. **Exact potential Path-1/Path-2 + syntactic-denominator clearing + robust `mu(y)`
   (`src/calculus/dsolve_exact.c`).** Build the potential from whichever coefficient integrates
   cleanly (`∫M dx` or `∫N dy`, no hot-path `Simplify`) — the E^(x y) family (1201) uses Path 2.
   Clear an inexact rational form by its common denominator D, trying the SYNTACTIC denominators
   (handles a negative exponential `E^(-x)` that `Together` mis-factors — 1233) then the
   `Together` denominator (summed `1/x`,`1/y` — 1216/1238); condition-free only. Try `mu(y)`
   whenever `mu(x)` yields no factor (`mu=Sin y` — 1214).
2. **Implicit verify `Together`s the residual (`src/calculus/dsolve_common.c`).** A `mu=Sin y`
   equation's telescoping implicit residual has uncancelled `Csc`/`Cot` poles that the numeric
   zero-test FALSE-NEGATIVES, wrongly rejecting a correct branch (1214); combining over a common
   denominator first (value-preserving) settles it.
3. **Separable fast numeric pre-filter (`src/calculus/dsolve_separable.c`).** `sep_find_split`
   numerically samples `F - g·h` and skips the symbolic zero-test on a clearly-nonzero point — a
   non-separable transcendental RHS (1201/1233) fast-declines in ~0 s not ~15 s, the pre-Exact
   cost that pushed those past the 8 s per-case timeout. Never rejects a genuine split.
4. **First-order IVP undecided-fit fall-through (`src/calculus/dsolve_common.c`).** A scalar
   first-order IVP whose fit bubbles back unevaluated (constant unfitted) now DECLINES so the
   cascade continues — closing the exact/homogeneous overlap 1205/1231 (Homogeneous's
   transcendental log-form shadows Exact's fitting polynomial integral). Gated to
   `nfun==1 && max_order==1` so an under-determined BVP keeps its free constant and higher-order
   series IVPs are untouched.

**Residue (1, research-grade — NOT a wrong answer, a bounded decline):**
- `1203` `x Log x + x y + (y Log x + x y) y' == 0` — an Abel equation of the second kind
  (class B), the M13-deferred `[_Abel]` class (Nasser's own solver and SymPy also fail).

Full per-case results: `reports/2.2.13.tsv`; bucketed report: `reports/2.2.13.md`.

---

## Section 2.2.14 — "Problems 1301 to 1400" (Boyce & DiPrima)

Corpus: `DE_examples_2214.m` — 100 records, **99 scalar (25 IVP), 1 system**. The
first Boyce & DiPrima section: 2nd-order-linear dominated — 38 `_with_linear_symmetries`,
19 `_missing_x` (const-coeff), 13 nonhomogeneous, 12 Emden–Fowler (mostly Euler–Cauchy),
11 exact.
`ctest -R dsolve_corpus_2_2_14_tests` · gate baseline **1**.

| Date | Solved | Solve % | Gap (non-PASS) | Notes |
|------|-------:|--------:|---------------:|-------|
| 2026-09-10 (baseline) | 90 / 100 | 90.0% | 10 | 0 FAIL, 10 UNEVAL: VoP verify-hangs, exact/Kovacic cascade-hangs, a singular-IVP fit, a transcendental-coeff series, a forced Bessel, plus the forced-Duffing residue. |
| 2026-09-10 (**M34**)  | **99 / 100** | **99.0%** | **1** | **+9, FAIL→0, 0 regression.** Four root-cause fixes below. |

**M34 fixes** (all shared-substrate — lift earlier sections, none regress):
1. **Numeric-zero verify short-circuit + robust variation of parameters
   (`src/calculus/dsolve_common.c`).** (a) `dsolve_verify_body` keeps a branch whose
   residual is NUMERICALLY zero at a spread of clean real points, before the
   symbolic `zero_test_decide` whose precision ladder climbs for >8 s on a residual
   that IS zero but carries Log/ArcTan branch cuts (the VoP answer of `y''+y==Tan[x]`
   / `2 Sec[x/2]`, 1337/1341); the reject path (a decidably-nonzero residual) is
   unchanged. (b) `dsolve_variation_of_parameters` reserves the symbolic-limit
   definite convolution for DiracDelta forcing and keeps a per-term INDEFINITE
   Wronskian integral (inert when non-elementary) for every other forcing, matching
   Mathematica's integral form and never entering the parametric DiffUnderInt
   escalation that blows up on `Tan`/`Sec`/arbitrary `g` (1350/1354).
2. **Bounded Kovacic + complex-pole gate (`src/calculus/dsolve_kovacic.c`).** A
   per-call wall-clock budget + bounded coefficient `Solve`, and Case-1c declines a
   NON-REAL pole (the complex-conjugate pole pair of `(x^3+1)y''+4x y'+y==0`, whose
   `ds_simplify(theta)` on the complex radicals spins for many seconds) — a genuinely
   Heun equation with no Liouvillian solution, which then falls to the Frobenius
   ordinary-point series (1392/1393). Its forcing closure also accepts an
   arbitrary-`g` inert-Integrate VoP particular (skips the un-numericizable
   `numeric_verify`), closing the forced Bessel operator 1350.
3. **Bounded ExactODE sub-solve + SeriesData IVP fit (`src/calculus/dsolve_exactode.c`,
   `dsolve_common.c`).** The exact reduction's recursive first-order sub-solve is
   `TimeConstrained` (its integrating-factor quadrature `Integrate[E^(-Cos[x]),x]` is
   non-elementary AND slow to give up, 1384), so it declines to Frobenius; and
   `dsolve_fit_constants` now takes `Normal[body]` of a SeriesData body so a series
   IVP fits `a[0]=C[1], a[1]=C[2]` at the IC point. The FIT_UNDECIDED fall-through is
   extended to a 2nd-order IVP, so a special-function general solution singular at the
   IC point (Bessel at x=0, 1381) declines to the origin-centred series.
4. **Frobenius series about the IC point (`src/calculus/dsolve_frobenius.c`).**
   `dsolve_frobenius_shifted_try` prefers the IVP's IC point as the expansion center
   when it is an ordinary point, and lifts the rational-only gate on that path, so a
   TRANSCENDENTAL-coefficient equation with no closed form Taylor-expands about x0
   (`x^2 y''+(x+1)y'+3 Log[x] y==0, y[1]==2, y'[1]==0`, 1385).

**Residue (1, NOT a wrong answer — a bounded decline):**
- `1360` `u''+u'+u^3/5==Cos[t]`, `u[0]==2, u'[0]==0` — a forced Duffing oscillator,
  genuinely nonlinear, with no closed form in **Maple, Mathematica, or SymPy** (all ✗
  in the source table).

Full per-case results: `reports/2.2.14.tsv`; bucketed report: `reports/2.2.14.md`.

---

## Section 2.2.15 — "Problems 1401 to 1500" (Boyce & DiPrima)

Corpus: `DE_examples_2215.m` — 100 records, **39 scalar (18 IVP), 61 systems** (53 2×2 +
8 3×3 constant-coefficient linear). The scalar half is high-order constant-coefficient linear
(1462–1489) plus **step / piecewise / Heaviside-forced 2nd-order IVPs** (1492–1500, the Boyce
& DiPrima Laplace-transform chapter).
`ctest -R dsolve_corpus_2_2_15_tests` · gate baseline **2**.

| Date | Solved | Solve % | Gap (non-PASS) | Notes |
|------|-------:|--------:|---------------:|-------|
| 2026-09-10 (baseline) | 98 / 100 | 98.0% | 2 | 0 FAIL. **But 1492–1500 were FALSE passes**: DSolve returned an inert `Integrate[UnitStep[…]·Cos,t]` the prelude could not numericize (UNK → trusted). Genuine gap was 11. |
| 2026-09-10 (**M35**)  | **98 / 100** | **98.0%** | **2** | **0 FAIL, 0 regression.** 1492–1500 now GENUINELY solve (verified `Piecewise` closed forms) via the new `DSolve\`PiecewiseForcing`. Gate baseline **2**. |

**M35 additions.**
- **`DSolve\`PiecewiseForcing`** (`src/calculus/dsolve_piecewise.c`, new method + builtin).
  Linear ODE IVPs with `Piecewise`/`UnitStep`/Heaviside forcing, by **interval continuation**:
  split at the forcing's breakpoints, solve each interval (recursing the scalar cascade), match
  `y,…,y^(n-1)` across each breakpoint, assemble a verified `Piecewise`. No `LaplaceTransform`
  needed. Runs before `UndeterminedCoefficients`/`LinearConstantCoefficients`, gated to
  step/piecewise-forced linear IVPs; bounded (re-entry guard + `TimeConstrained` sub-solves +
  wall-clock deadline + decline memo) with its own per-interval + per-IC numeric verify. E.g.
  `y''+4y==Piecewise[{{1,0<=t<Pi}},0], y(0)=1,y'(0)=0` → `Piecewise[{{1/4+3/4Cos[2t],t<Pi},{Cos[2t],t>=Pi}},0]`.
- **`dsolve_verify_body` distributional-keep** (`src/calculus/dsolve_common.c`): a residual
  carrying `UnitStep`/`Piecewise` is accepted on construction (joins DiracDelta/HeavisideTheta/
  definite-Integrate) — the numeric probe cannot run on it and `zero_test` spuriously rejected
  the correct step-forced answer of 1497 (`y''+4y==Sin[t]−UnitStep[t−2π]Sin[t]`).
- **Converter** (`tools/latex_ode_to_mathilda.py`): `\left\{…cases…\right.` → `Piecewise[…]`
  (`otherwise`→default, `±∞` bounds dropped), `\le`/`\ge`/`\infty` mapping, `Heaviside`→`UnitStep`.

**Residue (2, NOT wrong answers — bounded declines, both also ✗ in SymPy):**
- `1463` `t(t−1)y''''+E^t y''+4t²y==0` — 4th-order, transcendental coefficient.
- `1469` `t y'''+2y''−y'+t y==0` — 3rd-order, variable coefficient (the `3rd_high_linear`
  operator-factoring bucket, future work).

Full per-case results: `reports/2.2.15.tsv`; bucketed report: `reports/2.2.15.md`.

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
- **M24 (2026-09-08)** — §2.2.4 corpus (Problems 301–400). Two converter fixes (imaginary
  unit `i`; `y^{(n)}` derivative notation — §2.2.1/2/3 byte-identical), one harness fix
  (`dsResidVerdict` sweep-variable/ODE-parameter `k` collision → false FAIL, benefits every
  section), two solver fixes (`TrigReduce` trig-power/product forcing in
  `UndeterminedCoefficients`; `ComplexExpand` of numeric complex roots so an IVP fits).
  **99/100 (+6, FAIL→0, 0 regression);** sole residue 381 (variable-coeff Legendre-type,
  SymPy-failed). Anti-overfit units `t_m24_trig_power_forcing`, `t_m24_complex_cuberoot_ivp`.
  See §2.2.4 block.
- **M25 (2026-09-08)** — §2.2.5 corpus (Problems 401–500), a series-solution-heavy chunk. Three
  verified fixes: (1) Erf integrating-factor verify — `dsolve_verify_body` `ExpandAll`-normalises a
  Gaussian×Erf residual before `zero_test`, gated to Erf/Erfi, so the exact `y''+x y'+y==0` returns
  its Erf closed form instead of hanging (`POSSIBLE_ZEROQ_IMPROVEMENTS.md` #1 logs the core
  deficiency); (2) Kovacic fundamental-set independence guard — a coincident-exponent basis that
  collapses to `(C[1]+C[2]) y1` is rejected so Frobenius returns the correct two-parameter series
  (plus ExactODE declines a non-elementary reduced `Integrate`); (3) transcendental-coefficient
  Frobenius — `frobenius_regsing` Taylor-normalises `xP`, `x²Q` so an analytic transcendental
  coefficient (`6 Sin[x]/x`, removable at 0) yields a verified series. **99/100 (+4, 0 FAIL,
  0 regression);** sole residue 459 (irregular singular point, matches Mathematica). Anti-overfit
  units `t_m25_exact_erf`, `t_m25_kovacic_fundamental_set`, `t_m25_transcendental_frobenius`.
  See §2.2.5 block.
- **M26 (2026-09-09)** — §2.2.6 corpus (Problems 501–600), a forced-linear chunk.
  Green's-function (definite-integral) variation of parameters, DiracDelta sifting under
  `Integrate` (`integrate_dirac.c`), HeavisideTheta/DiracDelta value+derivative rules
  (`distributions.m`). **72/74 scalar (+15, 0 FAIL).** See §2.2.6 block.
- **M27 (2026-09-09)** — §2.2.7 corpus (Problems 601–700), the first **half-systems**
  section, and the wave that made the corpus harness **verify systems** by
  back-substitution (`dsolve_corpus_prelude.m`: `dsExplicitQ` generalised to a List
  function slot, system-skip removed). §2.1.2 (204 systems) and §2.2.6 (26 systems)
  re-baselined; §2.2.1–§2.2.5 (pure scalar) unchanged. Two root-cause engine fixes the
  corpus surfaced: (2) `dsolve_linsys_tidy` Expands an exponential body instead of
  Simplifying it — Simplify's zero-test hangs on a sum of exponentials with
  widely-separated real decay rates (`x'=-50x+20y, y'=100x-60y`, eigenvalues -10,-100);
  (3) `solveinv` seeds its parameter-mint counter past every `C[k]` already in the
  equation and `dsolve_extract_solutions` collapses the integer periodicity family to its
  principal branch, fixing a WRONG `y'=2x Sec[y]` answer from Solve reusing `C[1]` as the
  `2πk` index. **§2.2.7 93/100 (48/50 scalar + 45/50 systems), 0 FAIL, 0 regression.**
  Anti-overfit units `t_m27_system_verify`, `t_m27_separable_inverse_constant`,
  `t_m27_ivp_family_intact`. See §2.2.7 block.
- **M28 (2026-09-09)** — §2.2.8 corpus (Problems 701–800, Table 2.33, Edwards & Penney),
  a return to elementary first-order (100 scalar, 20 IVP, 0 systems). One root-cause
  engine fix: `dsolve_bernoulli.c` now declines fast when `y` appears non-algebraically
  (inside a transcendental function or a power exponent) — the Bernoulli exponent detector
  used to spin for seconds on the trig-in-`y` RHS of `2 x Sin[y]Cos[y] y'==4x²+Sin[y]²`
  (757), timing out the whole cascade (`$Aborted`) before `Linearizable` could solve it.
  **§2.2.8 99/100, 0 FAIL, 0 regression.** Sole residue 783 (`y'=1+x²+y²+x²y⁴`, a
  quartic-in-`y` needing the general `y=_G(x,y')` solve-for-`y`-and-differentiate method).
  Anti-overfit unit `t_m28_bernoulli_hang_trig_substitution`. See §2.2.8 block.
- **M29 (2026-09-09)** — §2.2.9 corpus (Problems 801–900, Edwards & Penney): 100 scalar
  (35 IVP), 0 systems, dominated by second-order linear constant-coefficient (homogeneous +
  nonhomogeneous), Euler/Emden–Fowler, plus complex-coefficient and `x(t)`-dependent-variable
  forms. **Fully solved out of the box — §2.2.9 100/100, 0 FAIL, 0 regression, baseline 0.**
  No ODE-solver fix was required; instead the wave added a general engine feature the corpus
  made visible: **piecewise derivatives of `Floor`/`Ceiling`/`Round`/`IntegerPart`/
  `FractionalPart`** (`src/calculus/deriv.c`), returning the Mathematica
  `Piecewise[{{v, cond}}, Indeterminate]` forms and composing with the chain rule. This
  upgrades 898 (`y''+9y==2 Sec[3x]`, a variation-of-parameters solution carrying a `Floor`
  branch-tracking term) from a non-numericizable "trust DSolve" pass to a genuine numeric
  residual ~0. Anti-overfit units `t_m29_sec_floor_verifies` (test_dsolve.c) and
  `test_rounding_deriv` (test_deriv.c). See §2.2.9 block.
- **M30 (2026-09-09)** — §2.2.10 corpus (Problems 901–1000, Edwards & Penney): 56 scalar
  (15 IVP) + **44 systems** — the most systems-heavy section (constant-coefficient linear of
  every order, Euler–Cauchy incl. Bessel/Legendre-type, and 44 first-order linear systems up
  to 4×4). **§2.2.10 100/100, 0 FAIL, 0 regression, baseline 0.** Two root-cause fixes, both
  about *forced* linear ODEs: (1) a forced constant-coefficient system with a real-irrational
  spectrum (924, `(1±√89)/2`) ran >90 s because `Integrate` rationalised the `1/λᵏ`
  variation-of-parameters coefficient into a hundreds-of-digit integer — `dsolve_linsys.c` now
  abstracts a real-irrational eigenvalue to a symbol before the integral and substitutes it
  back (0.6 s; complex/rational spectra stay concrete, keeping the complex-4×4 927 fast);
  (2) Kovacic gained an inhomogeneous closure — it accepts a forcing, de-obfuscates its
  fundamental set (`Exp[c Log u]→uᶜ`, split radicands, `PowerExpand`), and adds a
  variation-of-parameters particular (`dsolve_kovacic.c` + `dsolve_second_order_PQ_forced`),
  solving the Legendre-type 907. Anti-overfit units `t_m30_linsys_irrational_forcing`,
  `t_m30_kovacic_inhomogeneous`. See §2.2.10 block.
- **M31 (2026-09-09)** — §2.2.11 corpus (Problems 1001–1100, Edwards & Penney): 59 scalar
  (13 IVP) + **41 first-order constant-coefficient linear systems** (2×2 … 6×6, defective /
  complex spectra). **§2.2.11 100/100, 0 FAIL, 0 regression, baseline 0** — the scalar half
  (separable / quadrature / first-order-linear / 2nd-order const-coeff + Euler + Gegenbauer +
  Emden–Fowler + Liénard + Airy) solved entirely out of the box; both gaps were systems. Two
  root-cause fixes: (1) `2.2.11-1014` (a DAG solved by `TriangularSystem`) asked `Integrate` for
  `Integrate[e^{−9x}(7 C[k]e^{2x}−7 C[j]e^{2x}), x]` — a product of an exponential with a SUM of
  exponentials — which took the exponential-substitution path → 55 s and a branch-wrong
  `(−1)^{1/9}` antiderivative (kept because its residual is zero-test-undecidable). The ROOT fix
  is in `Integrate`: `src/calculus/integrate.c::try_linearity` now distributes a product over a
  sum factor (`Integrate` is linear: `c(g+h)→cg+ch`, committed only if every term closes), so the
  exponentials collapse to `e^{−7x}` — and it repairs the user-reported direct-`Integrate` bug too
  (`Integrate[e^{−9x}(a e^{2x}−b e^{2x}), x]`, was slow + branch-wrong); (2) `2.2.11-1001` (4×4,
  eigenvalues {16,32,48,64}) solves
  **correctly** but back-substitutes to a difference of `e^{64x}`-scale terms the prelude's
  20-digit sweep misread as nonzero (catastrophic cancellation) — `dsolve_corpus_prelude.m` now
  re-checks a not-small residual sample at 200-digit precision (monotone: never introduces a
  FAIL). The shared prelude fix also lifted four earlier sections whose correct-but-
  cancellation-heavy answers now verify; their ctest baselines were tightened to match (see the
  re-baseline note below). Anti-overfit units `t_m31_triangular_exp_forcing`,
  `t_m31_linsys_large_eigenvalue` (`test_dsolve.c`) + `test_linearity_distributes_product`
  (`test_integrate_dispatch.c`, the direct-`Integrate` regression guard). See §2.2.11 block.
- **M32 (2026-09-09)** — §2.2.12 corpus (Problems 1101–1200, Edwards & Penney): 100 scalar
  first-order ODEs (38 IVP), 0 systems. **§2.2.12 80/19/1 → 97/3/0** (FAIL→0, 0 crashes, 0
  regression) — the first §2.2.x section NOT solved out of the box; it exposed a wrong answer
  and real defects. Three shared-substrate root-cause fixes (they lift earlier sections, none
  regress): (1) **IVP constant-fitting** (`dsolve_common.c`) drops a branch an IC cannot be met
  on — one fitted to `Undefined` (the FAIL 1147) or an unsatisfiable scalar `±`/`Root` inverse
  branch, but only when a sibling fits (a lone singularity and a genuinely under-determined BVP
  keep their constant), and Bernoulli emits both real signs for an even `1−n` root; (2)
  **Separable** (`dsolve_separable.c`) accepts generic-parameter splits and gains an implicit
  twin `dsolve_separable_implicit_try` (∫dy/g==∫f dx+C, non-elementary integral kept unevaluated)
  for the non-invertible / Root-form / autonomous-non-elementary separables; (3) the **Integrate
  Gaussian→Erf/Ei/PolyLog recognizer** (`risch_special.c`) emitted a literal `x` for every
  integration variable — now threads the real variable, fixing every non-`x`-variable
  Bernoulli/linear solve. Residue 3 (research-grade, bounded declines): 1135/1200 (solvable-for-y/x,
  transcendental) and 1157 (Abel 2nd kind). Anti-overfit units `t_m32_*` (`test_dsolve.c`). See
  §2.2.12 block.
- **Next** — the **P≠0 confluent family** (~13 cases: 97/101/104 and kin) is the biggest
  Whittaker residue, pending an evaluator-robustness fix for the same-base symbolic-radical
  verify (`zero_test` / `HypergeometricPFQ`-numeric `$IterationLimit`). Then parabolic-
  cylinder/Hermite (`POLY_r ndeg=2`) and trig/hyperbolic-potential-with-`y'`. Also still open:
  M18 Stages 2/3 (μ(x,y'), μ(y,y') — BLOCKED on a non-elementary first-order solver),
  3rd/high-order operator factoring (needs non-`TimeConstrained` bounding), Abel Invariant
  Rational (deferred M13).
