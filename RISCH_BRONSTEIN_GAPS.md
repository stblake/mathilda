# Risch Transcendental Integrator — Gaps vs. Bronstein's Scheme

**Date:** 2026-07-14 *(refreshed from the 2026-07-12 original to reflect the P0/P1/P2/P5
landings on `main`)*
**Subject file:** `src/calculus/integrate_risch_transcendental.c` (5181 LoC) + `intrischnorman.c`
**Reference:** Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (2004)
**Forward target:** G. W. Cherry, *Integration in Finite Terms with Special Functions: The
Logarithmic Integral* (SIAM J. Comput. 1986) and *An Analysis of the Rational Exponential
Integral* (SIAM J. Comput. 1989) — for a future special-function generalization.
**Companion docs:** `RISCH_STATUS.md` (implementation-level case map, kept current),
`RISCH_AUDIT_FINDINGS.md` (Bronstein-faithfulness audit — no soundness defects),
`CHERRY_DESIGN.md` (the P5 architecture + landed substrate refactors).

> **Reading key.** Section/page references like *§5.10, p.163* are Bronstein book pages.
> Code references like `rt_hypertangent…:L4990` are function names / line numbers in
> `integrate_risch_transcendental.c` unless another file is named. Line numbers are as of
> this refresh (5181 LoC) and drift with edits — treat them as anchors, not addresses.

---

## 0. Executive summary

Mathilda's transcendental integrator is a genuine **recursive tower integrator** with an
explicit differential-tower data structure (`RtTower`), a real base-field Risch-DE solver
(`rde_base` with weak normalization, normal-denominator reduction, SPDE, principled cap-free
degree bounds), and Rothstein–Trager/LRT log parts. On the cases it targets it is *correct by
construction*. The 2026-07-12 original rated it a "heuristic parallel-Risch graft on a tower
skeleton." Since then the three structural divergences it identified have each been **materially
closed** — the Bronstein canonical/structure foundation now exists as tested modules on `main`,
the canonical split is wired into the recursive integrator, and tangents now integrate to **real**
output. What remains is genuine but narrower: the *live independence oracle* stays heuristic (by an
empirical finding, not an oversight), the *decision half* (proving non-integrability) is still
partial, and tangents integrate to real output via the direct §5.10 hypertangent case + real
reconstruction rather than a first-class `RT_TAN` tower monomial.

Status of the three original divergences:

1. **Canonical (normal/special) representation — NOW PRESENT AND WIRED.** `SplitFactor` /
   `SplitSquarefreeFactor` / `CanonicalRepresentation` (§3.5) are implemented in
   `src/calculus/risch_canonical.{c,h}` over the field arithmetic in `risch_field.{c,h}`, exposed
   as `Risch\`` builtins, and — as of the "canonical spine" landing — `rt_field_integrate`'s
   logarithmic case dispatches on the `f = f_p + f_s + f_n` split (`rt_canonical_split:L3558`)
   instead of the ad-hoc `PolynomialQuotient/Remainder` gate. The exponential-top path performs the
   split implicitly through its Hermite pipeline. The special-polynomial notion `S^irr` (`t` for
   exp, `t²+1` for tan) is realized inside `CanonicalRepresentation`'s `f_s`. **G-A1 is largely
   closed**; the residual is uniformity (per-case gates still exist alongside the spine) rather than
   absence.

2. **Risch Structure Theorem decision — NOW PRESENT (as a reusable oracle), deliberately NOT
   swapped into the live path.** `src/calculus/risch_structure.{c,h}` implements the ℚ-linear
   membership decision (`RationalSpan`) and the structure-theorem front-ends
   `Risch\`LogReducible` / `Risch\`ExpReducible` (Cor. 9.3.1) **and the real case**
   `Risch\`TanReducible` / `Risch\`ArcTanReducible` (Cor. 9.3.2, eqs. 9.14/9.15, with the disjoint
   `E`/`L`/`T`/`A` index-set partition), plus the standalone
   `Risch\`LogarithmicDerivativeOfRadical` (§5.12). A Phase-C experiment to replace
   `rt_class_primitive`/`rt_expand_logs` in the live integrator with this oracle **regressed** a
   Bronstein example and was reverted: the live tower already tolerates dependent generators via the
   tower-variable + diff-back mechanism, so G-A2's *live* impact was overstated. **G-A2/G-A3 are
   closed as decision facilities; the live independence path remains heuristic by choice.**

3. **Tangents — NOW REAL, via the direct §5.10 hypertangent case + real reconstruction of the
   exp route.** Two things landed. (a) A genuine **direct hypertangent family** integrator
   (`rt_hypertangent_case:L4990`, dispatched in the transcendental case) handles `Tan`/`Cot`
   (special `t²+1`, the §5.10 driver) and `Tanh`/`Coth` (special `t²−1`, the hyperbolic driver)
   over the real hypertangent monomial with **no Weierstrass / complex-exponential rewrite**,
   emitting real `ArcTan`/`Log[Cos]` — including irreducible-quadratic normal poles via the genuine
   §5.10 residue criterion (`rt_hypertan_family:L4811`, degree-≤2 factor gate). (b) For the general
   rational-trig integrands still routed through the exp tower, `rt_realify:L4580` (with `cx_reim`)
   reconstructs the **real** part of the I-laden complex-log output behind an exact numeric
   diff-back gate. Net: `∫tan x → −Log[Cos x]`, `Sec`/`Csc`/`Sec²`/`Sec³`/`1/(2+Cos)` etc. now come
   out real. **G-B4/G-B5 are substantially closed.** What remains (§2.3) is that this is *not* a
   first-class `RT_TAN` tower monomial: `RtKind` is still `{RT_LOG, RT_EXP}`, the RDE tangent
   branches and the deep Chapter-8 real recursion are absent, so tangents nested inside deeper
   towers, and `1/(a+b Sin[x])`, still fall outside.

The cross-cutting fourth gap is now **substantially closed** (2026-07-14, P3): Mathilda
**integrates and proves non-integrability** on the single-tower field class. `Risch\`ElementaryIntegralQ[f,x]`
is a sound Boolean decision predicate — `True`/`False`/unevaluated behind exact certificates only
(non-constant residue Thm 5.6.1(ii); a Ch.6 Risch DE with no rational solution; the §5.8 `Dc≠0`
certificate) — and `Integrate\`RischTranscendental` emits `Integrate::nonelem` on a proven-non-elem
unevaluated result. The mixed-resultant `r_s` partial-log-part landing (A2-2) now covers **both** the
primitive/log and the hyperexponential coupled paths (the coupled path splits `h = h_s + h_n` and
reconciles `h_s` through §5.9). What remains: `IntegrateNonLinearNoSpecial` (§5.11) is out of scope by
construction, and the constant-residue `Log` argument of an exp-monomial partial is unsimplified (a
`Simplify`/`LogToReal` form gap, not a correctness issue).

The remainder of this document maps every Bronstein subalgorithm to its Mathilda counterpart,
rates the (updated) gap, and closes with a Cherry-oriented roadmap whose P0/P1 prerequisites are
now satisfied.

---

## 1. Focus area A — Canonical representation & the structure theorems

### 1.1 What Bronstein requires (§3.4–3.5, §4.4, Ch. 9)

**Monomial (Def. 3.4.1).** `t` transcendental over `k` with `Dt ∈ k[t]`. Classified by
`δ(t) = deg_t(Dt)`: *primitive* (`Dt ∈ k`), *hyperexponential* (`Dt/t ∈ k`), *hypertangent*
(`Dt/(t²+1) ∈ k`), else *nonlinear*.

**Normal / special split (Def. 3.4.2, Thm 3.5.1).**
- `p` **normal**: `gcd(p, Dp) = 1` (its roots move under `D`).
- `p` **special**: `p | Dp` (its roots are fixed by `D` — e.g. `t` for an exponential,
  `t²+1` for a hypertangent).
- `SplitFactor(p) → (p_n, p_s)` computes this by gcd only.

**Canonical representation (§3.5).** Every `f ∈ k(t)` splits **uniquely** as
```
f = f_p + f_s + f_n     (polynomial part) + (special part b/d_s) + (normal part c/d_n)
```
via `CanonicalRepresentation` (`PolyDivide` + `SplitFactor` + `ExtendedEuclidean`). The whole
integration algorithm is defined on this decomposition: Hermite reduces `f_n`, the residue
criterion tests `f_n`'s remaining simple part, polynomial reduction handles `f_p`, and `f_s`
carries the special (exp/tan) poles.

**Structure Theorems (§9.3).** Given a tower `C(x)(t₁,…,tₙ)` with index sets
`E` (exponentials), `L` (logarithms), plus `T`/`A` (tangents/arctangents) in the real case
(Thm 9.3.2), the theorems reduce three questions to a **ℚ-linear system** (Cor. 9.3.1/9.3.2):
- (i) is `Da/a` the derivative of a tower element? (eq. 9.8/9.12)
- (ii) is `Db` the logarithmic derivative of a `K`-radical? (eq. 9.9/9.13)
- (iii/iv) the tangent analogues (eq. 9.14/9.15).
The Rothstein–Caviness theorem (§9.4) extends (i)/(ii) to log-explicit Liouvillian towers
(primitives allowed), which is exactly what recursion on primitives produces.

### 1.2 What Mathilda has (updated 2026-07-14)

| Bronstein construct | Mathilda status | Where |
|---|---|---|
| Monomial classification by `δ(t)`, `Dt/t`, `Dt/(t²+1)` | **Present in the foundation module** (`risch_field` classifies normal/special per monomial kind incl. hypertangent). In the *live tower* `RtKind` still distinguishes `RT_LOG`/`RT_EXP` only — tangents are handled by a **separate direct case**, not as a tower kind | `risch_field.c`; `RtKind:L3163`, `rt_hypertangent_case:L4990` |
| `SplitFactor` / `SplitSquarefreeFactor` (normal/special) | **Present.** `Risch\`SplitFactor` etc.; used by `CanonicalRepresentation` | `risch_canonical.c` |
| Special polynomials `S^irr` (e.g. `t`, `t²+1`) | **Present** inside `CanonicalRepresentation`'s special part `f_s` (`t^k` for exp, `t²+1` for tan) | `risch_canonical.c`, `rt_canonical_split:L3558` |
| `CanonicalRepresentation` (`f_p+f_s+f_n`) | **Present and WIRED** into the log/primitive top of the recursion; exp top splits implicitly via Hermite | `risch_canonical.c`; `rt_canonical_split:L3558`, `rt_field_integrate:L3766` |
| Structure Theorem ℚ-linear decision (Cor. 9.3.1/2) | **Present as a reusable oracle** (complex **and real/tangent** cases); **not swapped into the live independence path** (Phase-C finding: net-negative) | `risch_structure.c`; live path still `rt_class_primitive:L145`, `rt_expand_logs:L2542` |
| `IsLogarithmicDerivativeOfRadical` (§5.12) | **Present** as a standalone decision (`{n,u}` or `False`); the full radical-containing-log recursion is a pinned decline | `risch_structure.c` (`Risch\`LogarithmicDerivativeOfRadical`) |
| Triangular tower soundness (`Dcoef_i ∈ K_{i-1}`) | **Present** (a genuine, if narrower, structure check) | `rt_tower_build:L3207` |
| Rothstein–Trager resultant (§4.4, §5.6) | **Present** for the log part (single-kernel + tower), via `Integrate\`TranscendentalLogPart`; residue criterion now returns the elementary `r_s` logs of a **mixed** resultant | `rt_frac_lrt:L1681`, `rt_field_lrt_logpart:L3449`, `Integrate\`PartialLogPart` |

### 1.3 The gaps, concretely (updated)

- **G-A1 — Canonical splitting: LARGELY CLOSED.** `CanonicalRepresentation` exists and
  `rt_field_integrate`'s logarithmic case now dispatches on `f = f_p + f_s + f_n`
  (`rt_canonical_split`), behavior-preservingly (primitive monomials have `f_s ≡ 0`). The
  exponential top splits implicitly via its Hermite pipeline. **Residual:** the split is not yet the
  *sole* control-flow spine — several per-case denominator gates still run alongside it, so the
  "genuine mixed denominator" uniformity win is partial. This is now a *refactoring/uniformity*
  item, not an *absence*.

- **G-A2 — Structure-theorem independence decision: CLOSED as a facility; live impact reassessed.**
  The ℚ-linear system of Cor. 9.3.1/9.3.2 is implemented (`RationalSpan`, `LogReducible`,
  `ExpReducible`, `TanReducible`, `ArcTanReducible`). The original doc predicted large *live*
  completeness wins from swapping it in for `rt_class_primitive`/`rt_expand_logs`. The Phase-C
  experiment did exactly that and **regressed** the Bronstein `Log[x/Log[x]]` family (a composite
  log becomes a coupled two-log denominator the field integrator declines), so it was reverted. The
  live tower already tolerates ℚ-dependent generators through the tower-variable + diff-back
  mechanism. **Residual:** the oracle is available for callers that need a *proof* of dependence
  (e.g. Cherry), but the live path stays heuristic — an intentional trade, not a gap.

- **G-A3 — "New logarithm" / "logarithmic derivative of a radical": CLOSED for the integration
  form.** `Risch\`LogarithmicDerivativeOfRadical[f,x,mons]` returns `{n,u}` with `Du/u == n·f` (or
  `False`), exact/complete for the `f = Db` form that arises from integration. **Residual:** the
  full §5.12 recursion for radicals whose radicand contains log monomials (Bronstein's `x⁵ log x`
  case) is a documented, sound decline.

- **G-A4 — `RT_MAXK = 5` depth cap** (`:L3161`). **Unchanged.** Still a hard cap on tower depth in
  the live engine; a symptom of the ansatz approach (each level widens the `SolveAlways` system).
  Bronstein's scheme is depth-agnostic. Low practical impact but still a cap.

---

## 2. Focus area B — Extensions into tangents (the hypertangent case)

This was the largest original divergence. It is now **substantially closed for real output**,
though not via a first-class `RT_TAN` tower monomial.

### 2.1 Bronstein's direct real-tangent machinery

- **Hypertangent monomial (Def. 5.10.1, p.164).** `t` with `Dt/(t²+1) = a ∈ k`;
  `t = tan(∫a)`. Provided `√-1·Dt/(t²+1)` is not a logarithmic derivative of a `k(√-1)`-radical,
  `t` is a monomial over `k`, `Const` is preserved, and **the only special irreducible is
  `t²+1`** (Thm 5.10.1). Crucially, `√-1 ∉ k` is *maintained*, so tangents integrate **without**
  rewriting to complex exponentials.
- **Polynomial part (`IntegrateHypertangentPolynomial`, p.167).** `PolynomialReduce` to
  `deg ≤ 1`; the linear coefficient over `2a` yields the `c·log(t²+1)` (i.e. real
  `log(sec²)`/`arctan`) term. `Dc ≠ 0 ⇒ not elementary`.
- **Reduced part (`IntegrateHypertangentReduced`, p.169).** Peels poles at `t²+1` by solving a
  **coupled 2×2 system** `CoupledDESystem(0, 2mα, a, b)` per multiplicity `m`.
- **Coupled system (§8.4, p.264).** Writing `q* = y₁ + y₂√-1 ∈ k(√-1)`, the equation
  `Dq* + (b₀ - (…)√-1)q* = …` splits into a **real 2×2 system over `k`** (eq. 8.10):
  ```
  D(y₁,y₂)ᵀ + [[b₀, nη-b₂],[b₂-nη, b₀]]·(y₁,y₂)ᵀ = (z₁,z₂)ᵀ
  ```
  solved by `CoupledDECancelTan`/`CoupledDESystem` recursively over `k` — never leaving the
  real field. Real `arctan` output falls out. (`Example 8.4.1` walks a full case.)
- **RDE support (§6.2, §6.6).** `RdeSpecialDenomTan` (special part at `t²+1`) and
  `PolyRischDECancelTan` (cancellation) are the hypertangent branches of the Risch-DE pipeline.

### 2.2 What Mathilda now does

**Two paths, both producing real output:**

1. **Direct hypertangent family (`rt_hypertangent_case:L4990`, dispatched in the transcendental
   case at `:L5073`).** For `Tan[u]`/`Cot[u]` it builds the real hypertangent monomial (special
   `t²+1`) and runs `IntegrateHypertangentPolynomial` + the §5.10 reduced/residue driver
   (`Risch\`IntegrateHypertangent`, in `risch_hypertangent.c` / `risch_coupled.c`); for
   `Tanh[u]`/`Coth[u]` it uses the hyperbolic analogue (special `t²−1`). The normal-pole gate
   (`rt_hypertan_family:L4811`) admits any denominator whose **irreducible factors over ℚ have
   degree ≤ 2** (`rt_max_irr_degree` via `Factor`), and the driver realizes the quadratic-algebraic
   Rothstein–Trager residues as real `ArcTan` — **no Weierstrass, no complex exponentials.**
2. **Exp-tower route with real reconstruction (`rt_trig_frontend` → `rt_realify:L4580`).** General
   rational functions of the circular kernels still route through `TrigToExp` → the exponential
   machinery (now including `rt_exp_ratreduce_case`, closing the Laurent forms with `E^(−Ix)` in the
   denominator — `Sec`, `Sec²`, `Csc`, `1/(2+Cos)`, `1/(5+4Cos)`). The correct-but-I-laden log part
   is then reduced to its **real** part by `cx_reim` + `rt_realify` behind an exact numeric
   diff-back gate (a slip can only decline, never emit a wrong form). A two-argument
   `ArcTan[u,v] = atan2` derivative rule was added to `deriv.c` so the real forms verify.

Real results now include `∫tan x = −Log[Cos x]`, `Sec → (Log[2+2Sin]−Log[2−2Sin])/2`,
`1/(2+Cos) → real ArcTan`, `Cos/(1+Cos) → x − Sin/(1+Cos)`, `Tan[x]/(3+Tan[x]²)`,
`1/(3+Tan[x]²) → x/2 − ArcTan[Tan x/√3]/(2√3)`, and `Sec³`/`Csc³`/`Sec⁴`.

### 2.3 The gaps (updated)

- **G-B1 — `RT_TAN` tower monomial kind: STILL ABSENT.** `RtKind` remains `{RT_LOG, RT_EXP}`
  (`:L3163`). Tangents are handled by the *separate* `rt_hypertangent_case`, not as a first-class
  tower monomial the recursion can nest. (The structure module `risch_structure.c` does carry an
  `ArcTan` monomial kind for the ℚ-decision, but that is not the `RtTower` kind.) Consequence: a
  tangent monomial *inside* a deeper tower (`k` itself carrying further tangents) is not modeled.
- **G-B2 — Coupled differential system: PRESENT (base field).** `Risch\`CoupledDESystem` /
  `Risch\`IntegrateHypertangentReduced` (`risch_coupled.c`) implement the §8.1/§8.4 reduction to a
  single Risch DE over `C(i)(x)`, reproducing Bronstein Examples 5.10.3, 8.4.1, 5.10.2. **Residual:**
  the full Chapter-8 real recursion `CoupledDECancelTan` for coupled systems over a field `k` that
  itself carries further tangent monomials (deep towers) is not implemented.
- **G-B3 — `RdeSpecialDenomTan` / `PolyRischDECancelTan`: STILL ABSENT.** The RDE pipeline still
  lacks the tangent special-denominator and tangent cancellation branches (§6.2 Lemma 6.2.4 / §6.6).
  Not needed by the current base-field driver but required for the general tower-nested case.
- **G-B4 — Complex, non-real output: CLOSED for the covered families.** `∫tan x` now yields
  `−Log[Cos x]` (direct case) and the exp-route families are realified. **Residual:** a handful of
  forms (`1/(a+b Sin[x])`) still decline, and some real reconstructions rest on a *numeric* rather
  than symbolic diff-back because `Simplify` cannot reduce the emitted multiple-angle form — a
  `Simplify` deficiency tracked in `SIMPLIFY_DEFICIENCIES.md`, not an integrator gap.
- **G-B5 — `arctan` in the answer basis: CLOSED.** Real `ArcTan` now falls out of both the direct
  §5.10 residue criterion (irreducible-quadratic normal poles) and the realified exp route.

**Net:** tangents now integrate to **real** closed forms via a genuine §5.10 hypertangent driver
(base field) plus real reconstruction of the exp route — the original "highest-value chunk." The
remaining tangent work is the *general nested-tower* case: an `RT_TAN` tower kind, the RDE tangent
branches, and the deep `CoupledDECancelTan` recursion.

---

## 3. Case-by-case algorithm map (Bronstein Ch. 5–8 → Mathilda)

Ratings: **✔ present** · **◑ partial / ad-hoc** · **⌁ heuristic-substitute** · **✗ absent**.

| Bronstein algorithm (§, p.) | Purpose | Mathilda | Rating |
|---|---|---|---|
| `CanonicalRepresentation` (§3.5, p.100) | `f_p+f_s+f_n` split | **`risch_canonical_representation`; wired into `rt_field_integrate` log top (`rt_canonical_split`)** | ✔ |
| `SplitFactor` / `SplitSquarefreeFactor` (§3.5) | normal/special | **`risch_canonical.c`** | ✔ |
| `HermiteReduce` (§5.3, p.139) | kill repeated normal poles | `risch_hermite_reduce` (literal quadratic algorithm), wired into `rt_field_ratint` (ansatz removed) | ✔ |
| `PolynomialReduce` (§5.4, p.141) | reduce `f_p` for nonlinear `t` | **`risch_field_polynomial_reduce` (`Risch\`PolynomialReduce`)** + implicit in exp/log poly cases | ✔ |
| Liouville's Theorem (§5.5) | shape of the answer | assumed implicitly | ✔ (implicit) |
| `ResidueReduce` / residue criterion (§5.6, p.151) | log part **+ decide non-integrability** | log part ✔; **mixed-resultant `r_s` logs now returned** (`Integrate\`PartialLogPart`, Thm 5.6.1); Boolean decision half ✗ | ◑ |
| Integration of reduced funcs (§5.7) | reduce to `k[t]`/`k(t)` | canonical spine + per-case gates | ◑ |
| `IntegratePrimitive(Polynomial)` (§5.8, p.158) | primitive top | `rt_log_poly_case`, `rt_int_primitive_poly` | ✔ |
| `IntegrateHyperexponential(Polynomial)` (§5.9, p.161) | exp top | `rt_exp_poly_case`, `rt_int_hyperexp_poly` | ✔ |
| **`IntegrateHypertangent(Reduced/Polynomial)` (§5.10)** | **tangent top** | **`Risch\`IntegrateHypertangent{Polynomial,Reduced}` (`risch_hypertangent.c`/`risch_coupled.c`); live via `rt_hypertangent_case:L4990`** | ✔ (base field) |
| `IntegrateNonLinearNoSpecial` (§5.11, p.172) | nonlinear, `S^irr=∅` | — | ✗ |
| In-field integration (§5.12, p.175): `Du=f`, `Du/u=f`, radical | recognize (log-)derivatives | **`Risch\`LogarithmicDerivativeOfRadical`** + SolveAlways gates | ◑→✔ |
| **RDE** `WeakNormalizer` (§6.1, p.183) | strip spurious residues | `rde_weak_normalizer` | ✔ |
| `RdeNormalDenominator` (§6.1, p.185) | normal part of denom | present | ✔ |
| `RdeSpecialDenomExp` (§6.2, p.190) | special part, exp | in exp RDE path | ◑ |
| **`RdeSpecialDenomTan` (§6.2, p.192)** | special part, tan | — | ✗ |
| `RdeBoundDegree{Prim,Exp,NonLinear}` (§6.3) | degree bound by type | **`rt_rde_var_bound:L795` — principled, cap-free, monomial-type-aware, resonance-widened; unit-tested** | ✔ |
| `SPDE` (§6.4, p.203) | Rothstein SPDE | `rde_spde` | ✔ |
| `PolyRischDENoCancel{1,2}` (§6.5) | non-cancellation solve | `rde_polyrischde_nocancel1`, `…_integrate` | ✔ |
| `PolyRischDECancelPrim` (§6.6, p.212) | cancellation, primitive | partial (integrate path) | ◑ |
| **`PolyRischDECancelTan` (§6.6, p.215)** | cancellation, tangent | — | ✗ |
| `ParamRde*` (§7.1) | `Dy+fy=Σcᵢgᵢ` | approximated by multi-unknown `SolveAlways` | ⌁ |
| `LimitedIntegrate` (§7.2, p.246) | `f=Dv+Σcᵢwᵢ` | `rt_limited_integrate` / `rt_limited_field_integrate` (narrow) | ◑ |
| `ParametricLogarithmicDerivative` (§7.3, p.253) | new-log / radical decision | partial via `Risch\`LogarithmicDerivativeOfRadical` | ◑ |
| **`CoupledDESystem` (§8, p.259/265)** | 2×2 real system | **`Risch\`CoupledDESystem` (base field, `risch_coupled.c`)** | ✔ (base) |
| **`CoupledDECancelTan` (§8)** | recursion over tan-carrying `k` | — | ✗ |

Where Mathilda scores ✔ it is genuinely solid (the base-field Risch-DE stack is real and even has a
FLINT fast path; the foundation modules are Bronstein-audited). The remaining ✗/⌁ rows now cluster
narrowly around **(a) the deep tangent recursion** (`RT_TAN` tower kind, RDE tangent branches,
`CoupledDECancelTan`), **(b) the residue-criterion Boolean decision half**, and **(c) the
first-class parametric layer** (`ParamRde`).

---

## 4. The methodological gap: SolveAlways ansatz vs. deterministic reduction

Bronstein's algorithm is a **chain of exact reductions**: Hermite → polynomial reduction →
residue criterion → per-type polynomial integration, each step either shrinking the problem
deterministically or *proving* non-integrability. Degree bounds and the SPDE make every solve a
finite, certain computation.

Mathilda's tower cases (`rt_log_tower_case`, `rt_exp_tower_case`, `rt_recursive_tower_case`,
`rt_hyperexp_case`) still **build a bounded-degree ansatz `Q` and call `SolveAlways`**, then
(for the flat-tower cases) **verify by differentiating back** (`rt_verify_antideriv`). This remains
a Norman–Moses-flavored parallel-Risch shape *wearing a tower's clothes* (distinct from the actual
`Integrate\`RischNorman` engine in `intrischnorman.c`). Two things have improved since the original
writing:

- **Degree bounds are no longer heuristic.** The original called `rt_rde_var_bound` a heuristic. It
  is now Bronstein's exact leading-degree bound (`RdeBoundDegree`), monomial-type-aware and
  **cap-free** (the old cap-at-5 / cap-at-10 / `nmono ≤ 60` ceilings are gone), including the
  resonance-widening sub-case, and is **unit-tested directly** (`test_rde_var_bound`). So
  "completeness bounded by a too-small heuristic degree" is largely retired — a bound can now only
  ever *decline*, never truncate a genuinely solvable coefficient.
- **The canonical spine is threaded** into the log/primitive top, so that path is closer to the
  deterministic reduction chain than the pure ansatz.

Implications that persist:
- **Cannot fully prove non-integrability.** A `NULL` from `SolveAlways` still means "no solution in
  my ansatz box," not "no elementary antiderivative." The mixed-resultant partial-log-part landing
  (A2-2) returns the elementary logs and defers the rest; the `Dc≠0` tangent non-elementarity
  certificate *is* present in `risch_hypertangent.c`. But the general residue-criterion Boolean
  (§5.6 Thm 5.6.1(ii)) and `IntegrateNonLinearNoSpecial` are still absent.
- **The diff-back check still does real work** in the flat-tower and realification paths — a tell
  that those are search-with-verification, not decision procedures. (The recursive-tower case
  attempts an exact tower-variable identity `D_tower[Q] == F` before falling back to diff-back.)

This is not "wrong" — it is a legitimate, sound engineering choice. The remaining deterministic-core
work is the residue Boolean decision half (P3) and, for tangents, the RDE tangent branches (P2
residual).

---

## 5. Forward target — generalizing to Cherry's special functions

The user intends to extend this integrator toward **Cherry's algorithms** for integration in
terms of the **logarithmic integral `Li`/dilogarithm** (1986) and the **error function `Erf` /
exponential integral `Ei`** (1989). The prerequisites this section flagged (P0 canonical rep, P1
structure decision) are **now on `main`**; the P5 design and substrate refactors have landed
(`CHERRY_DESIGN.md`).

### 5.1 What Cherry adds (the extended Liouville theorem)

Both papers generalize Liouville's theorem to allow special-function terms in the antiderivative:

- **Log integral (1986, Thm 2.2 / 5.4).** For a Liouvillian `F` over an algebraically closed
  `C`, if `γ` has an integral in an **li-elementary** extension then
  ```
  γ = w₀' + Σ cᵢ (wᵢ'/wᵢ) + Σ dᵢ (uᵢ'/uᵢ),   with vᵢ' = uᵢ'/uᵢ   (⇒ ∫ = … + Σ dᵢ Li(uᵢ))
  ```
- **Rational exponential (1989, Thm 2.2 / 2.3).** For `∫ (f/g) eᵃ dx`,
  ```
  g = y' + f'y + Σ cᵢ(uᵢ'/uᵢ) + Σ dᵢ ūᵢ'   (⇒ … + Σ cᵢ e^{-a}Ei(ūᵢ) + Σ dᵢ e^{-a}Erf(ūᵢ))
  ```
  with the sharp structural bounds: **error-function arguments are quadratic over `C` and there
  are at most two of them** (Thm 3.2) — because `f+β=g²` has ≤2 solutions `β`.

Notably, **Cherry deliberately avoids the general Risch DE** `S'−fS=g` and instead uses
**undetermined coefficients + linear systems** (1989, p.893). So the dependency surface is
*lighter* than full Bronstein — but it is very much the **canonical/structure foundation** of
§1, not the tower-ansatz machinery. Both papers reduce to the single ansatz `∫γ = v + Σ kᵢ SF(aᵢ)`;
`CHERRY_DESIGN.md` unifies them as `rt_tower_solve` with extra basis terms whose arguments come from
two new candidate generators.

### 5.2 Cherry's dependency list vs. Mathilda (updated)

| Cherry subroutine | Needed for | Mathilda status | Gap |
|---|---|---|---|
| Squarefree + partial-fraction decomposition over the tower | both | present (`Apart`, poly subsystem, FLINT) | — |
| Hermite reduction | both (desingularize denominators) | **literal `risch_hermite_reduce`** | ✔ |
| Resultant / Rothstein–Trager (log residues, factor args) | Li args, Ei args | present (`TranscendentalLogPart`) + **`Integrate\`RothsteinTragerResultant` (R3)** | ✔ |
| Canonical `f = f_p + f_s + f_n` split by top monomial (Lemma 5.1) | both | **present + wired** (`risch_canonical_representation`) | ✔ (was G-A1) |
| Structure decision `vᵢ = rᵢaₙ + …` (Roca79 / Cor. 9.3) | both, recursion | **present** (`LogReducible`/`ExpReducible`/`Tan…`/`LogarithmicDerivativeOfRadical`) | ✔ (was G-A2/G-A3) |
| **Σ-decomposition** (Thm 4.4) — Li argument generator | Li | **absent** — the Li candidate-generator | **build** |
| **Completing-square decision** `f+β=g²` (Thms 3.1–3.3) — Erf argument generator | Erf | **`PolynomialSqrt` (R5) landed**; the ≤2/quadratic `β²−Uβ+V` search itself still to build | **build** |
| **Risch Main Thm part (b)** → reduce `cᵢ,dᵢ` to a linear system over constants | both | `rt_tower_solve` is the reduce-to-linear-system solver and accepts special-function basis terms in `Q` (contract documented, R1) | ◑→✔ |
| Special-function output menu / registry | both | **`RtSpecialForm` registry (R2)** — the four `rt_try_*` recognizers are registry entries; `rt_special_case` loops it | ✔ (seam ready) |

### 5.3 What Mathilda has toward Cherry — and the remaining delta

Mathilda still ships the **ad-hoc pattern recognizers** that emit Cherry's target functions —
`rt_try_erf:L327`, `rt_try_ei:L385`, `rt_try_li:L432`, `rt_try_dilog:L526` — but they are now
**registry entries** in `RT_SPECIAL_FORMS[]:L572`, the seam Cherry's decision procedures grow into
(R2). They remain the **narrow, enumerated special cases** of Cherry's general decision procedures:
- `rt_try_erf`'s "`f'/f` is degree-1" test is a one-off instance of the completing-square `f+β=g²`
  search (whose `PolynomialSqrt` primitive, R5, now exists).
- `rt_try_li`'s `c·w^{p-1}w'/Log[w]` template is a single Σ-decomposition shape.

Generalizing per Cherry means **replacing these templates with the two candidate-generators
(Σ-decomposition; completing-square `β²−Uβ+V`) feeding the `rt_tower_solve` linear-system coefficient
solve**, on top of the (now-present) canonical representation and structure oracle.

### 5.4 Consequence for sequencing (updated)

Cherry's prerequisites — **(a)** canonical normal/special representation, **(b)** the
structure-theorem ℚ-decision / logarithmic-derivative decision, and **(c)** a general
"reduce to a linear system over constants" solver — are **now satisfied**: (a) and (b) are on
`main`, and (c) is `rt_tower_solve` with its documented special-function-basis extension point. The
five substrate refactors (`CHERRY_DESIGN.md §3`: `PolynomialSqrt`, `RothsteinTragerResultant`, the
`RtSpecialForm` registry, the `rt_tower_solve` contract, the recursion hook) landed
behavior-preservingly. **What remains for P5 is the two argument-generators + the extended-Liouville
ansatz assembly**, not foundation work. The tangent work (§2) is orthogonal to Cherry and largely
done for the base field.

---

## 6. Prioritized roadmap (status refreshed 2026-07-14)

Ordered by (leverage × how many downstream gaps it unlocks). Each item is scoped to a coherent
increment; none requires touching `src/external/`.

**P0 — Canonical representation & special polynomials. ✅ DONE (modules + wiring).**
`SplitFactor`/`CanonicalRepresentation` (§3.5) are implemented in `risch_canonical.{c,h}` over
`risch_field.{c,h}`, exposed as `Risch\`` builtins, and the `f = f_p + f_s + f_n` split is **wired
into** `rt_field_integrate`'s logarithmic case (`rt_canonical_split`), replacing the ad-hoc
`PolynomialQuotient/Remainder` gate behavior-preservingly; the exp top splits implicitly via Hermite.
The special-irreducible set `S^irr` (`t` for exp, `t²+1` for tan) lives in `f_s`. Foundation coverage:
278 assertions across `test_risch_{field,canonical,structure,hypertangent}.c` plus
`test_risch_canonical_wiring.c`. Closes **G-A1**. *Residual (refactor, not blocker):* make the
canonical split the sole control-flow spine, retiring the remaining per-case denominator gates.

**P1 — Structure-theorem decision + logarithmic-derivative decision. ✅ DONE (A+B), live swap
declined by finding (C).**
`risch_structure.{c,h}` implements `RationalSpan` and the structure-theorem front-ends for the
**complex** case (`LogReducible`/`ExpReducible`, Cor. 9.3.1) **and the real case**
(`TanReducible`/`ArcTanReducible`, Cor. 9.3.2, eqs. 9.14/9.15, with the disjoint `E`/`L`/`T`/`A`
index-set partition and an `ArcTan` monomial kind), plus the standalone
`LogarithmicDerivativeOfRadical` (§5.12, exact for the `f = Db` integration form). Builtins + tests
(`test_risch_structure{,_real}.c`, `test_risch_logderiv_radical.c`), leak-clean. **Phase C** (swap
the oracle into the live independence path) was attempted and **reverted** — it regressed
`Log[x/Log[x]]`; the live tower already tolerates dependent generators via tower-variable +
diff-back, so G-A2's live impact was overstated. P1 ships as the **reusable `Risch\`` decision
builtins**. Closes **G-A2/G-A3** as facilities.

**P2 — Hypertangent extension (the tangent focus area). ✅ DONE for the base field + real output;
deep-tower recursion remains.**
Landed: (a) `Risch\`PolynomialReduce` + `IntegrateHypertangentPolynomial` (§5.4/§5.10 polynomial
part, `risch_hypertangent.c`); (b) `Risch\`CoupledDESystem` + `IntegrateHypertangentReduced` (§8.1/
§8.4 base-field coupled solver over `C(i)(x)`, `risch_coupled.c`, reproducing Bronstein Ex. 5.10.3 /
8.4.1 / 5.10.2); (d) the **live** direct hypertangent driver `rt_hypertangent_case` routing
`Tan`/`Cot` (special `t²+1`) and `Tanh`/`Coth` (special `t²−1`) **instead of** Weierstrass, with an
irreducible-quadratic normal-pole gate → real `ArcTan`; plus **real reconstruction** (`rt_realify` +
`cx_reim`) of the exp-tower route's I-laden output and a two-arg `ArcTan[u,v]` derivative. Real now:
`∫tan x = −Log[Cos x]`, `Sec`/`Csc`/`Sec²`/`Sec³`/`Sec⁴`, `1/(2+Cos)`, `1/(5+4Cos)`, `Cos/(1+Cos)`,
`1/(3+Tan²)`. Closes **G-B2 (base), G-B4, G-B5.** *Remaining (G-B1/G-B3 + deep recursion):* a
first-class `RT_TAN` tower monomial kind, the RDE tangent branches `RdeSpecialDenomTan` /
`PolyRischDECancelTan` (§6, c), and the full Chapter-8 `CoupledDECancelTan` recursion for tangents
nested in deeper towers; and the residual declines (`1/(a+b Sin[x])`).

**P3 — Residue-criterion decision half + tight degree bounds. ✅ DONE (Boolean verdict shipped).**
*Done (degree bounds):* Bronstein's exact cap-free `RdeBoundDegree` (`rt_rde_var_bound`,
monomial-type-aware, resonance-widened, unit-tested) across the RDE and all flat-tower/Hermite ansatz
sites — the "heuristic `Nx/Ntop`" concern is retired. *Done (partial log part, A2-2, Thm 5.6.1):*
mixed Rothstein–Trager resultants return the elementary `r_s` logs (`Integrate\`PartialLogPart`, κ_D
split) plus an unevaluated `Integrate[remainder,x]` and an FTC rule `D[Integrate[f,x],x]→f`.
*Done (2026-07-14, the Boolean decision half):* `Risch\`ElementaryIntegralQ[f,x]` — a sound decision
predicate returning `True` (elementary antiderivative exhibited), `False` (**proven** non-elementary:
non-constant residue Thm 5.6.1(ii) via a decide-mode `Integrate\`TranscendentalLogPart` marker; a
Ch.6 Risch DE with no rational solution; or the §5.8 `Dc≠0` certificate), or unevaluated outside the
tower field scope. `Integrate\`RischTranscendental` emits `Integrate::nonelem` on a proven-non-elem
unevaluated result. **The RDE integrators are now genuine decision procedures**: `rt_field_rde`'s
base case routes *every* `C(x)` RDE (rational RHS included — the former `rt_is_poly(p,x)` guard
diverted `E^x/x`, `E^(x²)` to a bounded ansatz) to the complete `rde_base`, so every `NULL` is an
authoritative "no rational solution", not a scope decline. Decides `E^x/x`/`E^(x²)`/`1/Log[x]`/
`E^(E^x)/(1+E^(E^x))` → `False`; elementary siblings → `True`. *Remaining:* `IntegrateNonLinearNoSpecial`
(§5.11) — out of scope by construction (the live tower builds only log/exp monomials, no nonlinear
no-special monomial arises). *The hyperexponential coupled path now also returns the partial log part*
for a mixed resultant (`rt_field_hyperexp_hermite` splits `h = h_s + h_n`, reconciles `h_s` through
§5.9, reports `h_n` unintegrated — depth-≥2 towers; gated by the tower diff-back), so the last A2-2
residual is closed. The residual form quality (the constant-residue `Log` argument is unsimplified for
exp monomials — a `Simplify`/`LogToReal` gap, not a correctness issue) is tracked separately.

**P4 — Parametric layer. ◑ NARROW.**
`LimitedIntegrate` exists in a narrow form (`rt_limited_integrate` / `rt_limited_field_integrate`);
`ParametricLogarithmicDerivative` is partially covered by `Risch\`LogarithmicDerivativeOfRadical`.
`ParamRde` (§7.1) is still approximated by multi-unknown `SolveAlways` rather than a first-class
facility. Directly reused by the primitive case and by Cherry's coefficient solve.

**P5 — Cherry special-function generalization (the forward target). 📐 DESIGNED + SUBSTRATE LANDED;
engines to build.**
`CHERRY_DESIGN.md` sets the architecture (`∫γ = v + Σ kᵢ SF(aᵢ)` as `rt_tower_solve` + argument
generators behind the `RtSpecialForm` registry, on the P0/P1 foundation). Five behavior-preserving
substrate refactors landed: **R5** `PolynomialSqrt` (completing-square primitive), **R3**
`Integrate\`RothsteinTragerResultant` (Ei/log argument generator), **R2** the `RtSpecialForm`
registry (the four `rt_try_*` recognizers are entries), **R1** the documented `rt_tower_solve`
extension-point contract, **R4** the lower-field recursion hook (deferred to the Cherry module). *To
build:* (a) Σ-decomposition (Li args, 1986 Thm 4.4); (b) the completing-square `β²−Uβ+V` search
(Erf args, 1989 Thms 3.1–3.3, ≤2/quadratic bounds); (c) the extended-Liouville ansatz assembly that
subsumes and retires `rt_try_{erf,ei,li,dilog}` as recognizers, replacing them with decision
procedures.

**Dependency graph:** `P0 → P1 → {P3, P4}` ✅ satisfied; `P4 → P5`, `P1 → P5` (P1 ✅, P4 narrow);
`P2` parallel — base field ✅, deep-tower recursion outstanding.

---

## 7. Quick reference — where the pieces live now

| Capability | Home | Status |
|---|---|---|
| `SplitFactor` / `CanonicalRepresentation` / `S^irr` | `risch_canonical.c`; wired via `rt_canonical_split:L3558` | ✅ |
| Structure-theorem ℚ-linear decision (complex + real) | `risch_structure.c` (`Risch\`{Log,Exp,Tan,ArcTan}Reducible`) | ✅ (facility) |
| `IsLogarithmicDerivativeOfRadical` | `risch_structure.c` (`Risch\`LogarithmicDerivativeOfRadical`) | ✅ (integration form) |
| Hypertangent kind + detection (tower) | `RtKind:L3163` — still `{RT_LOG,RT_EXP}`; direct case `rt_hypertangent_case:L4990` | ◑ (no `RT_TAN` tower kind) |
| `CoupledDESystem` | `risch_coupled.c` (`Risch\`CoupledDESystem`, base field) | ✅ (base) |
| `RdeSpecialDenomTan` / `PolyRischDECancelTan` / `CoupledDECancelTan` | RDE block — not yet | ✗ |
| `IntegrateHypertangent{Polynomial,Reduced}` | `risch_hypertangent.c` / `risch_coupled.c`; live via `rt_hypertangent_case` | ✅ |
| Real reconstruction of the exp-trig route | `rt_realify:L4580` + `cx_reim` | ✅ |
| Residue Boolean decision half | `Risch\`ElementaryIntegralQ` (§5.6 residue + Ch.6 RDE + §5.8 Dc≠0 certificates); `Integrate::nonelem` message; RDE base case now fully authoritative (`rde_base`) | ✅ |
| Cherry Σ-decomposition / completing-square / linear-coeff solve | new `src/calculus/cherry_*.c` (to build); substrate: `PolynomialSqrt`, `RothsteinTragerResultant`, `RtSpecialForm` registry | 📐 |

---

*Sources: implementation inventory of `integrate_risch_transcendental.c` (5181 LoC) /
`intrischnorman.c` and the `risch_{field,canonical,structure,hermite,hypertangent,coupled}` modules;
the git history 2026-07-12…14 (P0 foundation → canonical spine → P1 structure theorem → P2 real
tangents → P5 substrate); `RISCH_STATUS.md`, `RISCH_AUDIT_FINDINGS.md`, `CHERRY_DESIGN.md`;
Bronstein 2nd ed. §3.4–3.5, §4.4, Ch. 5–9; Cherry 1986 (Log Integral) & 1989 (Rational Exponential
Integral).*
