# Bronstein-faithfulness audit of the integrate code — consolidated findings

**Date:** 2026-07-13
**Reference:** Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (2004).
**Scope:** a line-by-line review of the existing transcendental Risch / rational
integration code against the book, per the directive "all existing integrate
code should be thoroughly reviewed against Bronstein's book for inconsistencies;
when in doubt, re-implement using the book."

Three independent audits (full detail in `RISCH_AUDIT_A1.md`, `A2.md`, `A3.md`).

## Headline

**No CRITICAL findings and no soundness defects anywhere.** Every closed form
the integrator emits is gated by an exact certificate (SolveAlways polynomial
identity, exact tower-derivation identity, or SPDE polynomial identity) before it
ships, so a wrong antiderivative cannot be produced. The only gaps are
**completeness** gaps — inputs that are declined (returned unevaluated) where a
more complete implementation would integrate them or *prove* non-integrability.
These align with the existing roadmap (P1/P3).

## A1 — Chapter 6 Risch DE stack (`integrate_risch_transcendental.c` L828–1272)

Verdict: **faithful, no CRITICAL/MAJOR.** `rde_spde`, the inlined
`RdeNormalDenominator` / `RdeBoundDegreeBase`, `rde_weak_normalizer`,
`rde_polyrischde_nocancel1/2`, and the `rde_base` driver were each transcribed
against the book and match — the base-field (C(x), D=d/dx) specialization
(SplitFactor = identity, no special part) is provably sound and complete for
that field.

- **F5 — RESOLVED (2026-07-13):** audited `flint_rde_base_solve_fg`
  (`src/poly/flint_bridge.c` L791–880) line-by-line against Bronstein §6.1–6.5.
  It is a faithful `fmpq_poly` port of the same stack the Expr path runs, and
  **every "0 = no solution" it returns is authoritative** (backed by a theorem),
  so it cannot silently drop an integrable in-scope input:
  * `RdeNormalDenominator` (base `dₙ=1`): `h = gcd(eₙ, eₙ')`, guard `eₙ ∤ h²`
    ⇒ decline. Authoritative by **Cor. 6.1.1(ii)** (a solution ⇒ `eₙ | dₙh²`).
    Outputs `a=h`, `b=hf−Dh`, `c=(h²/eₙ)·num_g` match eq. (6.2) exactly.
  * `RdeBoundDegreeBase` (inlined, L834–853) computes the book's
    `n = max(0, dc−max(db, da−1))` and, in the cancellation branch `db=da−1`,
    raises to `max(0, m, dc−db)` with `m=−lc(b)/lc(a)∈ℤ` — verified equal to the
    p.199 box (the base term already equals `dc−db` when `max(db,da−1)=db`). This
    is a proven **upper** bound (Cor. 6.3.1), so the ladder is never given too
    small an `n`.
  * `fq_spde` mirrors the **Thm 6.4.1** box (`n<0 ∧ c≠0` ⇒ no-solution;
    `g∤c` ⇒ no-solution; `ExtendedEuclidean` via `xgcd`+rem; reconstruct
    `q=αh+β`) — every decline is authoritative.
  * `fq_polyrischde_nocancel1/integrate` are the §6.5 top-down peels (base:
    `D=d/dx` always satisfies the non-cancellation hypothesis); they decline
    only when a required monomial degree falls outside `[0,n]`, i.e. genuinely
    no bounded solution. Each `nr=1` solution is exact **by construction**
    (residual `c` driven to 0), so no diff-back is needed.
  * Scope: returns `-1` (defer to Expr) when `f` is rational (`FD≠1`) or `g` is
    not univariate over ℚ — outside its declared authority. Correct.
  * Empirical: solves constructed RDEs (poly & rational-`g` repeated-pole,
    self-verifying `D[y]+fy−g≡0`) and declines exactly on Bronstein's own
    non-elementary **Ex. 6.1.1** (`eᵗ/t`) and **Ex. 6.3.2** (`e^{−t²}`).
  * Corroboration: the Expr fallback (`rde_base` L1176–1249) is a line-for-line
    mirror of the same stack — two independently-coded paths agree on the
    algorithm. No code change required.
- F1/F3/F4: behaviour-equivalent simplifications (documented). F2: an unreachable
  defensive guard. None actionable.

## A2 — Chapter 5 reductions & residue (`integrate_risch_transcendental.c`, `intrat.c`)

Verdict: **no CRITICAL.** Primitive-polynomial (§5.8), hyperexponential-
polynomial (§5.9, Laurent range exact), and the LRT log part are faithful. Two
**MAJOR completeness** gaps:

- **A2-1 — Hermite is an ansatz, not literal `HermiteReduce`** (`rt_hermite_try`
  L1866–1897, `rt_field_ratint` L3481–3534; §5.3). The code fits an undetermined-
  coefficient numerator whose x-coefficients are *polynomials* of heuristic degree
  `Nx = max(deg_x)+2`, rather than the book's squarefree + per-multiplicity
  `ExtendedEuclidean` loop with arbitrary *rational* k=C(x) coefficients.
  Integrands whose Hermite numerator is rational (or higher x-degree) are declined.
  → **RESOLVED (2026-07-13):** the literal quadratic `HermiteReduce` (§5.3, Thm
  5.3.1) is implemented in `src/calculus/risch_hermite.{c,h}` as
  `Risch\`HermiteReduce[f,t,deriv] -> {g,h,r}` with the exact invariant
  `f = D[g]+h+r`, handling arbitrary rational k=C(x) coefficients. It is **wired
  into `rt_field_ratint`, and the ~150-line undetermined-coefficient ansatz is
  deleted** — the live proper-part integrator is now Hermite (§5.3) + residue
  criterion (§5.6), with a tower diff-back self-verify. The full transcendental
  integrator suite still passes; rational-coefficient repeated-pole cases the
  ansatz declined now close. The two pre-existing `intrat.c` (Ch.2 rational)
  test failures are unrelated (confirmed identical at HEAD via A/B).
  **Extended to the exponential top (2026-07-13):** `rt_field_hyperexp_coupled`
  now delegates to `rt_field_hyperexp_hermite`, the literal hyperexponential
  pipeline (HermiteReduce §5.3 → residue logs §5.6 → coupling reconciliation
  `P = h + r − D_tower[L]` → IntegrateHyperexponentialPolynomial §5.9, each Laurent
  coefficient by a Risch DE). Its unified `SolveAlways` ansatz is likewise
  **deleted**; repeated exponential poles carrying rational lower-field
  coefficients now close (e.g. `∫ D[1/((1+x)(1+E^x)²)]`). Both proper-part paths
  (log/primitive and exponential) are now ansatz-free.
- **A2-2 — residue criterion discards the elementary `r_s` logs**
  (`intrat_log_part_core`, `intrat.c` L1759–1775; Thm 5.6.1). The gate requires
  the *entire* residue resultant to be constant; when the resultant mixes constant
  and non-constant roots (`r = r_s·r_n`), Bronstein still returns the elementary
  logs from `r_s`, but the code returns NULL and discards them (Bronstein's own
  Example 5.6.2 cannot be reproduced even partially). → Fix: split by κ_D, return
  the `r_s` logs, report `r_n` unevaluated.
  → **RESOLVED (2026-07-14):** `intrat_log_part_core` now performs the κ_D split
  (per-irreducible `intrat_freeq_all` partition of the squarefree resultant
  buckets), reconstructs the `r_n` remainder root-free via the resultant–gcd
  association `d_n = gcd(d, Res_z(r_n, a−z·D(d)))`, and returns
  `Integrate`PartialLogPart[logs, remainder]`. `Integrate`RischTranscendental`
  surfaces it as `logs + Integrate[remainder, x]`; a new FTC rule
  `D[Integrate[f,x],x]→f` closes the diff-back. Non-regressing (only the
  previously-always-declined mixed branch changes).
  → **Extended to the hyperexponential coupled path (2026-07-14, P3 tail):**
  `rt_field_hyperexp_hermite` now accepts a partial log part too — it splits the
  simple part `h = h_s + h_n`, reconciles the constant-residue `h_s` through the
  §5.9 Laurent step (`P = h_s + r − D[L_s]`), integrates it, and reports `h_n`
  unintegrated (so a depth-2 `∫ E^x/(1+E^(E^x)) + 1/(x+E^(E^x))` returns
  `E^x − Log[1+E^(E^x)] + Integrate[1/(x+E^(E^x)),x]` instead of declining).
  Gated by the exact tower-variable diff-back, so the κ_D reconstruction (not
  exact over a transcendental monomial) can only decline, never ship a wrong
  remainder. Tests in `tests/test_risch_residue_split.c`.
- A2-3: `b=0` is only consumed as "decline" — a *semi*-decision procedure (never a
  positive "provably non-elementary" verdict). OK-by-design for a CAS; the decision
  half is roadmap P3.

## A3 — Rational integration + Bronstein foundation modules

Verdict: **no CRITICAL/MAJOR — all faithful ports.**

- `intrat.c` Hermite (§2.2) and LRT log part (§2.5): exact port of Mack's linear
  box; the pseudo-remainder-vs-subresultant substitution is sound.
- `risch_field.c` normal/special (Def 3.4.2) and `PolynomialReduce` (§5.4): correct;
  field gcd over k=C(other vars) is correct by Gauss's lemma.
- `risch_canonical.c` `SplitFactor` / `SplitSquarefreeFactor` /
  `CanonicalRepresentation` (§3.5): match their boxes; `f_s`, `f_n` proper.
- `risch_structure.c` `RationalSpan` / `LogReducible` / `ExpReducible` (Cor 9.3.1):
  the Q-span decision is sound (no false "reducible").
- `risch_hypertangent.c` (Thm 5.10.2): `c = coeff(r,t)/(2α)` and the Dc≠0
  certificate match the box.

## Disposition

- **No fixes are required for correctness.** The code is sound.
- A2-1 (literal Hermite over the tower) is **RESOLVED** — both the log/primitive
  and exponential proper-part paths are now ansatz-free. The remaining MAJOR
  completeness gap A2-2 (residue `r_s` logs) and the decision half (A2-3) are
  folded into roadmap **P3**.
- **F5** (audit `flint_rde_base_solve_fg`) is **RESOLVED (2026-07-13)** —
  faithful `fmpq_poly` port, every decline authoritative; no defect, no code
  change. The Track-1 audit is now fully closed.
