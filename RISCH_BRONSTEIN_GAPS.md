# Risch Transcendental Integrator — Gaps vs. Bronstein's Scheme

**Date:** 2026-07-12
**Subject file:** `src/calculus/integrate_risch_transcendental.c` (4612 LoC) + `intrischnorman.c`
**Reference:** Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (2004)
**Forward target:** G. W. Cherry, *Integration in Finite Terms with Special Functions: The
Logarithmic Integral* (SIAM J. Comput. 1986) and *An Analysis of the Rational Exponential
Integral* (SIAM J. Comput. 1989) — for a future special-function generalization.

> **Reading key.** Section/page references like *§5.10, p.163* are Bronstein book pages.
> Code references like `rt_hypertangent…:L4386` are function names / line numbers in
> `integrate_risch_transcendental.c` unless another file is named.

---

## 0. Executive summary

Mathilda's transcendental integrator is a genuine **recursive tower integrator** with an
explicit differential-tower data structure (`RtTower`), a real base-field Risch-DE solver
(`rde_base` with weak normalization, normal-denominator reduction, SPDE, degree bounds), and
Rothstein–Trager/LRT log parts. On the cases it targets it is *correct by construction*. But
measured against Bronstein's canonical scheme it diverges in three structural ways that bound
its completeness and, more importantly, its ability to *decide* integrability:

1. **No canonical (normal/special) representation.** Bronstein's entire Chapter 5 operates on
   the `f = f_p + f_s + f_n` splitting (§3.5) computed by `SplitFactor`. Mathilda has **no
   normal/special split and no notion of the special polynomials `S`**. Every case instead
   guesses a bounded ansatz and calls `SolveAlways`, gated by structural certificates and (for
   the flat-tower cases) a diff-back check. This works but is a *heuristic parallel-Risch
   graft on a tower skeleton*, not the deterministic reduction chain.

2. **No Risch Structure Theorem decision.** Mathilda decides monomial independence with an
   ad-hoc **rational-commensurability test on exponents** (`rt_class_primitive`) plus
   **`Log` expansion** (`rt_expand_logs`). Bronstein's Structure Theorems (§9.3, Cor. 9.3.1/2)
   reduce independence and the "is this a new logarithm / logarithmic-derivative-of-a-radical"
   questions to a **ℚ-linear system** over the constants. Mathilda builds no such system, so it
   misses non-multiplicative ℚ-linear relations and cannot run the residue criterion's
   decision half.

3. **Tangents go through complex exponentials, not the hypertangent extension.** This is the
   single largest focus-area gap. Bronstein integrates `tan`/`arctan` **directly over a real
   hypertangent monomial** (`Dt/(t²+1) ∈ k`, §5.10) via a **coupled 2×2 differential system**
   (§8.4), producing real `arctan` output. Mathilda's `rt_trig_frontend` does
   `TrigToExp → integrate → ExpToTrig`, which is exactly the complex-exponential route Bronstein
   built the hypertangent case to avoid — and the file itself documents the resulting I-laden,
   un-real-simplified output as a known limitation (`:L4394`).

A fourth, cross-cutting gap: Mathilda **integrates but rarely proves non-integrability**. It
returns `NULL` ("I can't") where Bronstein's `ResidueReduce` and the polynomial-case theorems
return a Boolean that *proves* "no elementary integral exists." The decision-procedure property
— the thing that makes it Risch rather than a heuristic — is only partially present.

The remainder of this document maps every Bronstein subalgorithm to its Mathilda counterpart,
rates the gap, and closes with a Cherry-oriented roadmap.

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

### 1.2 What Mathilda has

| Bronstein construct | Mathilda status | Where |
|---|---|---|
| Monomial classification by `δ(t)`, `Dt/t`, `Dt/(t²+1)` | **Partial.** Distinguishes `RT_LOG` vs `RT_EXP` only; **no hypertangent kind, no nonlinear kind** | `RtKind`, `RtTower` `:L3099` |
| `SplitFactor` (normal/special) | **Absent.** No normal/special concept anywhere | — |
| Special polynomials `S^irr` (e.g. `t`, `t²+1`) | **Absent** (`t`-at-0 multiplicity is counted ad-hoc in the exp cases) | `rt_var_mult_at_zero`, hyperexp cases |
| `CanonicalRepresentation` (`f_p+f_s+f_n`) | **Absent.** Cases split num/den and gate "den free of `t`" / "den = `c·t^M`" by hand | `rt_field_integrate:L3799` |
| Structure Theorem ℚ-linear decision (Cor. 9.3.1/2) | **Absent.** Replaced by rational-exponent commensurability + `Log` expansion | `rt_class_primitive:L142`, `rt_expand_logs:L2463` |
| Triangular tower soundness (`Dcoef_i ∈ K_{i-1}`) | **Present** (a genuine, if narrower, structure check) | `rt_tower_build:L3295` |
| Rothstein–Trager resultant (§4.4, §5.6) | **Present** for the log part (single-kernel + tower), via `Integrate`TranscendentalLogPart` | `rt_frac_lrt:L1656`, `rt_field_lrt_logpart:L3378` |

### 1.3 The gaps, concretely

- **G-A1 — No canonical splitting.** Because there is no `f = f_p + f_s + f_n`, each case
  re-derives an ad-hoc denominator classification and a bounded ansatz. Consequences: (a) the
  same reduction logic is re-implemented per case with per-case gates; (b) genuine mixed
  denominators (a normal factor *and* a special factor together, with a repeated normal pole)
  fall between the cases; (c) the code cannot state "this is the reduced part" cleanly, which is
  the precondition for the residue criterion's decision property.

- **G-A2 — No structure-theorem independence decision.** `rt_class_primitive` only detects
  **multiplicative** ℚ-relations among exponential exponents (`e^{x/2}, e^{x/3} → e^{x/6}`).
  It does **not** detect additive ℚ-relations that the structure theorem catches, e.g. that
  `log(x)` and `log(2x)` differ by a constant, that `log(x²+x)` = `log x + log(x+1)` couples to
  existing logs, or that a proposed new `exp`/`log` is dependent via eq. 9.8/9.9. `rt_expand_logs`
  patches the log side syntactically (product/power laws) but is not the ℚ-linear-system decision
  and will miss anything not exposed by structural expansion.

- **G-A3 — No decision of "new logarithm" / "logarithmic derivative of a radical."** This is
  Bronstein's §5.12 in-field integration and Cor. 9.3.1(ii). It is the gate that both the
  hyperexponential special-denominator step and Cherry's methods need. Mathilda approximates it
  inside `SolveAlways` ansätze (folding would-be new logs back), but has **no standalone
  `IsLogarithmicDerivativeOfRadical` decision**.

- **G-A4 — `RT_MAXK = 5` depth cap** (`:L3095`). A hard cap on tower depth. Bronstein's scheme
  is depth-agnostic; the cap is a symptom of the ansatz approach (each level widens the
  `SolveAlways` system).

---

## 2. Focus area B — Extensions into tangents (the hypertangent case)

This is where Mathilda most diverges from Bronstein's design intent.

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

### 2.2 What Mathilda does instead

`rt_trig_frontend` (`:L4386`):
```
TrigToExp(f)  →  {exp-poly, frac, hyperexp, expsum} cases  →  ExpToTrig
```
i.e. **Weierstrass/complex-exponential rewriting**. `Tan[u]` becomes a rational function of
`E^{iu}` (or `Tan[u/2]` kept as an exp-derived generator). The exponential cases then integrate
it.

### 2.3 The gaps

- **G-B1 — No hypertangent monomial kind.** `RtKind` has only `RT_LOG`/`RT_EXP`. There is no
  `Dt/(t²+1) ∈ k` classification, so `tan` is never treated as a first-class real monomial.
- **G-B2 — No coupled differential system.** Nothing implements `CoupledDESystem` /
  `CoupledDECancelTan` (§8, eq. 8.1/8.10). This is the core missing primitive; without it the
  hypertangent reduced case cannot be done Bronstein-style at all.
- **G-B3 — No `RdeSpecialDenomTan` / `PolyRischDECancelTan`.** The Risch-DE pipeline
  (`rde_base:L1136`) implements the primitive/hyperexponential/base machinery but **not** the
  tangent special-denominator or tangent cancellation branches (§6.2 Lemma 6.2.4 / §6.6).
- **G-B4 — Complex, non-real output (documented).** Because of the `TrigToExp` route, results
  are I-laden and not simplified back to real closed forms — e.g. `∫tan x` yields
  `I x - Log[1 + E^{2Ix}]` instead of `-Log[Cos x]` (`:L4394`). Bronstein's whole §5.10 exists to
  avoid exactly this.
- **G-B5 — No `arctan` in the answer basis.** Bronstein emits real `arctan` from the `t²+1`
  special factor and from the structure-theorem tangent decision (eq. 9.14/9.15). Mathilda's log
  part is `Log`-only; `arctan` appears only accidentally via `TranscendentalLogPart`'s real
  ArcTan reconstruction on the exponential rewrite, and the tangent structure decision is absent.

**Net:** Mathilda can *sometimes* integrate tangents (via the exp rewrite) but does not implement
the hypertangent extension, the coupled system, or the real-output guarantee that define
Bronstein's treatment of tangents. This is the highest-value chunk to build for parity.

---

## 3. Case-by-case algorithm map (Bronstein Ch. 5–8 → Mathilda)

Ratings: **✔ present** · **◑ partial / ad-hoc** · **�’ heuristic-substitute** · **✗ absent**.

| Bronstein algorithm (§, p.) | Purpose | Mathilda | Rating |
|---|---|---|---|
| `CanonicalRepresentation` (§3.5, p.100) | `f_p+f_s+f_n` split | — | ✗ |
| `SplitFactor` / `SplitSquarefreeFactor` (§3.5) | normal/special | — | ✗ |
| `HermiteReduce` (§5.3, p.139) | kill repeated normal poles | **`risch_hermite_reduce` (literal quadratic algorithm), wired into `rt_field_ratint` (ansatz removed)** | ✔ |
| `PolynomialReduce` (§5.4, p.141) | reduce `f_p` for nonlinear `t` | only implicit in exp/log poly cases | ◑ |
| Liouville's Theorem (§5.5) | shape of the answer | assumed implicitly | ✔ (implicit) |
| `ResidueReduce` / residue criterion (§5.6, p.151) | log part **+ decide non-integrability** | log part ✔ (`rt_frac_lrt`, `rt_field_lrt_logpart`); **decision half ✗** | ◑ |
| Integration of reduced funcs (§5.7) | reduce to `k[t]`/`k(t)` | per-case gates | ◑ |
| `IntegratePrimitive(Polynomial)` (§5.8, p.158) | primitive top | `rt_log_poly_case:L622`, `rt_int_primitive_poly:L3867` | ✔ |
| `IntegrateHyperexponential(Polynomial)` (§5.9, p.161) | exp top | `rt_exp_poly_case:L1364`, `rt_int_hyperexp_poly:L3966` | ✔ |
| **`IntegrateHypertangent(Reduced/Polynomial)` (§5.10)** | **tangent top** | `TrigToExp` rewrite only | ✗ (see §2) |
| `IntegrateNonLinearNoSpecial` (§5.11, p.172) | nonlinear, `S^irr=∅` | — | ✗ |
| In-field integration (§5.12, p.175): `Du=f`, `Du/u=f`, radical | recognize (log-)derivatives | folded into SolveAlways gates | ◑ |
| **RDE** `WeakNormalizer` (§6.1, p.183) | strip spurious residues | `rde_weak_normalizer:L1050` | ✔ |
| `RdeNormalDenominator` (§6.1, p.185) | normal part of denom | `:L1173` | ✔ |
| `RdeSpecialDenomExp` (§6.2, p.190) | special part, exp | in exp RDE path | ◑ |
| **`RdeSpecialDenomTan` (§6.2, p.192)** | special part, tan | — | ✗ |
| `RdeBoundDegree{Prim,Exp,NonLinear}` (§6.3) | degree bound by type | `rde_base` bound + `rt_rde_var_bound:L4092` (heuristic) | ◑ |
| `SPDE` (§6.4, p.203) | Rothstein SPDE | `rde_spde:L909` | ✔ |
| `PolyRischDENoCancel{1,2}` (§6.5) | non-cancellation solve | `rde_polyrischde_nocancel1:L985`, `…_integrate:L1022` | ✔ |
| `PolyRischDECancelPrim` (§6.6, p.212) | cancellation, primitive | partial (integrate path) | ◑ |
| **`PolyRischDECancelTan` (§6.6, p.215)** | cancellation, tangent | — | ✗ |
| `ParamRde*` (§7.1) | `Dy+fy=Σcᵢgᵢ` | approximated by multi-unknown `SolveAlways` | �’ |
| `LimitedIntegrate` (§7.2, p.246) | `f=Dv+Σcᵢwᵢ` | `rt_limited_integrate` / `rt_limited_field_integrate:L3932` (narrow) | ◑ |
| `ParametricLogarithmicDerivative` (§7.3, p.253) | new-log / radical decision | — (folded into ansatz) | ✗ |
| **`CoupledDESystem` / `CoupledDECancelTan` (§8, p.259/265)** | 2×2 real system | — | ✗ |

Where Mathilda scores ✔ it is genuinely solid (the base-field Risch-DE stack is real and even
has a FLINT fast path). The ✗/�’ rows cluster around **(a) the canonical/structure foundation**,
**(b) the tangent/coupled-system column**, and **(c) the parametric decision problems** — and
those three clusters are causally linked: (c) needs (a), and the tangent column needs the coupled
system which needs (b).

---

## 4. The methodological gap: SolveAlways ansatz vs. deterministic reduction

Bronstein's algorithm is a **chain of exact reductions**: Hermite → polynomial reduction →
residue criterion → per-type polynomial integration, each step either shrinking the problem
deterministically or *proving* non-integrability. Degree bounds and the SPDE make every solve a
finite, certain computation.

Mathilda's tower cases (`rt_log_tower_case`, `rt_exp_tower_case`, `rt_recursive_tower_case`,
`rt_hyperexp_case`) instead **build a bounded-degree ansatz `Q` and call `SolveAlways`**, then
(for the flat-tower cases) **verify by differentiating back** (`rt_verify_antideriv`). This is
essentially the **Norman–Moses parallel-Risch heuristic wearing a tower's clothes** — even though
the file correctly notes it is *distinct* from the actual `Integrate`RischNorman` engine in
`intrischnorman.c`.

Implications:
- **Completeness is bounded by the ansatz.** If the degree bound (`rt_rde_var_bound`, the
  `Nx`/`Ntop` heuristics) is too small, a genuinely elementary integral is missed and returned as
  `NULL`. Bronstein's per-type degree bounds (§6.3) are provably tight; Mathilda's are heuristic.
- **Cannot prove non-integrability.** A `NULL` from `SolveAlways` means "no solution in my
  ansatz box," not "no elementary antiderivative." Bronstein's residue criterion + degree bounds
  yield a real proof. Mathilda's `rt_try_*` special-function recognizers and the tangent
  polynomial `Dc≠0` test (which Bronstein *does* use to prove non-elementarity) are absent.
- **The diff-back check is doing real work.** Because the ansatz can over-count freedom (spurious
  `SolveAlways` solutions), the flat-tower cases *must* verify. That is a tell that these cases
  are search-with-verification, not decision procedures. (The recursive-tower case does better: it
  attempts an exact tower-variable identity `D_tower[Q] == F` at `:L4336` before falling back to
  diff-back.)

This is not "wrong" — it is a legitimate, pragmatic engineering choice, and the structural
certificates keep it sound. But it is the reason the three focus-area gaps exist: a deterministic
Bronstein core would replace the ansatz with reductions and gain both completeness and the
decision property.

---

## 5. Forward target — generalizing to Cherry's special functions

The user intends to extend this integrator toward **Cherry's algorithms** for integration in
terms of the **logarithmic integral `Li`/dilogarithm** (1986) and the **error function `Erf` /
exponential integral `Ei`** (1989). This section maps what Cherry needs onto what exists.

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
§1, not the tower-ansatz machinery.

### 5.2 Cherry's dependency list vs. Mathilda

| Cherry subroutine | Needed for | Mathilda status | Gap |
|---|---|---|---|
| Squarefree + partial-fraction decomposition over the tower | both | present (`Apart`, poly subsystem, FLINT) | — |
| Hermite reduction | both (desingularize denominators) | `rt_hermite_try` / classical | ◑ usable |
| Resultant / Rothstein–Trager (log residues, factor args) | Li args, Ei args | present (`TranscendentalLogPart`) | ✔ |
| **Σ-decomposition** (Thm 4.4) — split a rational fn into irreducible pieces that become `Li` arguments | Li | **absent** — this is the Li candidate-generator | **build** |
| **Completing-square decision** `f+β=g²` (Thms 3.1–3.3) — `Erf` argument generator | Erf | **absent** (only pattern-matched in `rt_try_erf:L323`) | **build** |
| **Risch Main Thm part (b)** → reduce coefficients `cᵢ,dᵢ` to a **linear system over constants** | both | **absent as a general facility** (only per-case `SolveAlways`) | **build** |
| Logarithmic-derivative / structure decision (`u'/u` vs `u'`; which monomial) | both, recursion | weak (`rt_expand_logs`, commensurability) | **build** (= G-A2/G-A3) |
| Factored, normalized transcendental tower | both | `RtTower` exists but lacks canonical/special split | **extend** (= G-A1) |

### 5.3 What Mathilda already has toward Cherry — and why it isn't enough

Mathilda ships **ad-hoc pattern recognizers** that emit exactly Cherry's target functions:
`rt_try_erf` (`:L323`), `rt_try_ei` (`:L381`), `rt_try_li` (`:L428`), `rt_try_dilog` (`:L522`).
These are hand-written templates (e.g. `K·E^{ax²+bx+c}` with `f'/f` linear → `Erf`). They are the
**narrow, enumerated special cases** of Cherry's **general decision procedures**:
- `rt_try_erf`'s "`f'/f` is degree-1" test is a one-off instance of Cherry's completing-square
  `f+β=g²` search.
- `rt_try_li`'s `c·w^{p-1}w'/Log[w]` template is a single Σ-decomposition shape.

Generalizing per Cherry means **replacing these templates with the two candidate-generators
(Σ-decomposition; completing-square) feeding a general linear-system coefficient solve**, on top
of a proper canonical representation.

### 5.4 Consequence for sequencing

Cherry sits on the **same foundation** whose absence drives §1's gaps. Concretely, the
prerequisites for a faithful Cherry layer are: **(a) canonical normal/special representation
(G-A1)**, **(b) the structure-theorem ℚ-linear decision / logarithmic-derivative decision
(G-A2/G-A3)**, and **(c) a general "reduce to a linear system over constants" solver** (Cherry's
use of Risch Main Thm part (b), which is also Bronstein's `LimitedIntegrate` / `ParamRde` flavor).
The tangent work (§2) is *orthogonal* to Cherry and can proceed in parallel. Building Cherry on
today's `SolveAlways` ansatz without (a)–(c) would reproduce the current template approach at
larger scale, not Cherry's decision procedure.

---

## 6. Prioritized roadmap

Ordered by (leverage × how many downstream gaps it unlocks). Each item is scoped to a coherent
increment; none requires touching `src/external/`.

**P0 — Canonical representation & special polynomials (unblocks almost everything).**
Implement `SplitFactor`/`CanonicalRepresentation` (§3.5) on the existing `Expr`/poly layer and
thread `f = f_p + f_s + f_n` through `rt_field_integrate`. Add the special-irreducible set
`S^irr` per monomial kind (`t` for exp, `t²+1` for tan). Closes **G-A1**; is the substrate for
P1, P3, P4.
> **Status (2026-07-12): the standalone modules are done; wiring remains.**
> `src/calculus/risch_field.{c,h}` — monomial derivation; field gcd/divexact/**num-den/
> divmod/xgcd/diophantine** in `k[t]` over `k=C(x,…)`; normal/special classification.
> `src/calculus/risch_canonical.{c,h}` — `SplitFactor`, `SplitSquarefreeFactor`,
> **`CanonicalRepresentation`** (`f_p+f_s+f_n`, §3.5 p.103). Exposed as the `Risch\``
> builtins; `tests/test_risch_canonical.c` covers Bronstein Examples 3.5.1/3.5.2, the unique
> exp/hypertangent canonical decompositions, field division by an `x`-coefficient, and
> per-monomial classification. Clean `-Wall -Wextra`, leak-clean.
>
> **Test coverage (2026-07-13):** the field primitives `FieldGCD`/`DivExact`/`NumDen`/
> `ExtendedEuclidean`/`Diophantine` are now exposed as `Risch\`` builtins and covered
> directly by a new `tests/test_risch_field.c` (property round-trips: Bézout `u a+v b=g`,
> `a=q b+r`, `b dn+c ds=r`; degree bounds; field-unit behaviour; robustness). Total Bronstein
> foundation coverage: **278 assertions across 34 functions** in `test_risch_{field,canonical,
> structure,hypertangent}.c` (up from 88/14). **Remaining:** thread the
> `f_p+f_s+f_n` split (and the special-polynomial set `S^irr` per monomial kind) through the
> recursive integrator `rt_field_integrate`, replacing its ad-hoc per-case denominator gates.

**P1 — Structure-theorem decision + logarithmic-derivative decision.**
Build the ℚ-linear system of Cor. 9.3.1/9.3.2 (and §5.12 `IsLogarithmicDerivativeOfRadical`).
Replace `rt_class_primitive`+`rt_expand_logs` as the independence oracle where it matters. Closes
**G-A2/G-A3/G-A4-depth**; prerequisite for both a real decision procedure and Cherry.
> **Status (2026-07-12): core done (complex case).** `src/calculus/risch_structure.{c,h}` —
> `risch_rational_span` (the ℚ-linear membership decision via `Numerator[Together[…]]` +
> `SolveAlways`) and the structure-theorem front-ends `Risch\`LogReducible` / `Risch\`ExpReducible`
> (Cor. 9.3.1(i)/(ii)), returning the rational coefficient vector or `False`. Builtins +
> `tests/test_risch_structure.c` (log(x²)=2log x, exp(2x)=exp(x)², exp(x+log x)=x·eˣ, and
> genuine-new cases), leak-clean. `warn_svars` now honours the mute counter. **Remaining:** the
> **real/tangent case** (Cor. 9.3.2, eqs. 9.14/9.15 — the `T`/`A` index sets), the standalone
> `IsLogarithmicDerivativeOfRadical` wrapper (`n·f = Du/u`), and replacing `rt_class_primitive`/
> `rt_expand_logs` in the live integrator with this oracle.

**P2 — Hypertangent extension (the tangent focus area).**
(a) Add the `RT_TAN` monomial kind (`Dt/(t²+1) ∈ k`) and its detection. (b) Implement
`CoupledDESystem`/`CoupledDECancelTan` (§8.4). (c) Implement `RdeSpecialDenomTan` +
`PolyRischDECancelTan` (§6). (d) `IntegrateHypertangent{Polynomial,Reduced}` (§5.10). Route `tan`
through this **instead of** `TrigToExp`, yielding real `arctan`/`log(cos)` output. Closes
**G-B1..G-B5**. Independent of P0/P1 in principle, but cleaner on top of P0.
> **Status (2026-07-12): polynomial part done.** `risch_field_polynomial_reduce` (§5.4, the
> nonlinear-monomial reduction) and `src/calculus/risch_hypertangent.{c,h}`
> `IntegrateHypertangentPolynomial` (§5.10 d, polynomial part) land with builtins
> `Risch\`PolynomialReduce` / `Risch\`IntegrateHypertangentPolynomial` and
> `tests/test_risch_hypertangent.c` (Example 5.10.1, real `∫tan x`, the `Dc≠0` non-elementarity
> certificate, `a=2` scaling), leak-clean — the first real tangent integration with no complex
> exponentials.
>
> **Status (2026-07-13): reduced part + coupled system done.**
> `src/calculus/risch_coupled.{c,h}` — `CoupledDESystem` (b) via the §8.1 reduction to a single
> Risch DE over `C(i)(x)` (`Risch\`CoupledDESystem`, built on the newly-exposed base-field
> `Risch\`RischDE` = `rde_base`), and `IntegrateHypertangentReduced` (§5.10, p.169) peeling `t²+1`
> poles per multiplicity via Eq. (5.20) (`Risch\`IntegrateHypertangentReduced`). Reproduces
> Bronstein Examples 5.10.3 `(c,d)`, 8.4.1 `(s1,s2)`, and 5.10.2 (`∫ sin x / x` → non-elementary).
> Extensive stress-tested `tests/test_risch_coupled.c`, leak-clean. This is the base-field (single
> `tan`-over-`C(x)`) coupled solver — correct and real-valued for the common case. **Remaining:**
> the full Chapter-8 real recursion `CoupledDECancelTan` (for coupled systems that arise over a
> field `k` itself carrying further tangent monomials — deep towers), the RDE tangent branches
> `RdeSpecialDenomTan` / `PolyRischDECancelTan` (c), the `RT_TAN` detection + `IntegrateHypertangent`
> driver (Hermite + Residue + reduced + polynomial) (a/d), and routing `tan` through this in the
> live integrator instead of `TrigToExp`.

**P3 — Residue-criterion decision half + tight degree bounds.**
Return Booleans that *prove* non-elementarity (§5.6 Thm 5.6.1(ii), the `Dc≠0` tangent test,
`IntegrateNonLinearNoSpecial`). Replace heuristic `Nx/Ntop` bounds with `RdeBoundDegree{Prim,Exp,
Tan,NonLinear}` (§6.3). Converts the tower cases from search-with-verify toward decision
procedures; retires the diff-back fallback where a real certificate now exists.

**P4 — Parametric layer.**
`ParamRde` (§7.1), `LimitedIntegrate` (§7.2), `ParametricLogarithmicDerivative` (§7.3) as
first-class facilities. Directly reused by the primitive case and by Cherry's coefficient solve.

**P5 — Cherry special-function generalization (the forward target).**
On P0/P1/P4: (a) Σ-decomposition (Li args, 1986 Thm 4.4); (b) completing-square decision
(Erf args, 1989 Thms 3.1–3.3, with the ≤2 / quadratic bounds); (c) general linear-system
coefficient solve; (d) the extended-Liouville ansatz. Subsumes and retires `rt_try_{erf,ei,li,
dilog}` as recognizers, replacing them with decision procedures.

**Dependency graph:** `P0 → P1 → {P3, P4}`; `P4 → P5`; `P1 → P5`; `P2` parallel (best after P0).

---

## 7. Quick reference — where the missing pieces would live

| New capability | Natural home |
|---|---|
| `SplitFactor` / `CanonicalRepresentation` / `S^irr` | new `rt_canonical_*` block near `RtTower` (`:L3099`) |
| Structure-theorem ℚ-linear decision | new `rt_structure_decide_*`; replaces callers of `rt_class_primitive:L142` |
| `IsLogarithmicDerivativeOfRadical` | new; used by RDE special-denom + Cherry |
| Hypertangent kind + detection | extend `RtKind`/`rt_tower_build:L3141` |
| `CoupledDESystem` | new file `src/calculus/coupled_de.c` (parallels `rde_base`) |
| `RdeSpecialDenomTan` / `PolyRischDECancelTan` | extend the RDE block (`:L829–1272`) |
| `IntegrateHypertangent*` | new case alongside `rt_exp_poly_case`/`rt_log_poly_case` |
| Σ-decomposition / completing-square / linear-coeff solve | new `src/calculus/cherry_*.c`; retire `rt_try_{erf,ei,li,dilog}` |

---

*Sources: implementation inventory of `integrate_risch_transcendental.c`/`intrischnorman.c`;
Bronstein 2nd ed. §3.4–3.5, §4.4, Ch. 5–9; Cherry 1986 (Log Integral) & 1989 (Rational
Exponential Integral).*
