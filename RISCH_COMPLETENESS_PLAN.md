# Risch Transcendental — Implementation Plan to a Complete Bronstein Decision Procedure

**Date:** 2026-07-14
**Reference:** Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (2004).
**Subject:** `src/calculus/integrate_risch_transcendental.c` (~5450 LoC) + the `risch_*` modules.
**Companion docs:** `RISCH_STATUS.md` (case map), `RISCH_BRONSTEIN_GAPS.md` (gap survey),
`RISCH_AUDIT_FINDINGS.md` (soundness audit — no defects).

> This is the *build spec* that turns Mathilda's sound, wide-coverage integrator into a
> genuine **decision procedure** on the transcendental (Liouvillian, non-nonlinear) tower:
> for every elementary input it either returns an elementary antiderivative or **proves**
> none exists. It is written to be executed in small, independently-landable increments,
> each green + committed before the next.

---

## 0. Definition of "done"

A complete Bronstein transcendental decision procedure over the tower `C(x)(t₁,…,tₙ)` of
logarithmic / hyperexponential / hypertangent monomials satisfies, at **every** level `k(t)`:

- `HermiteReduce → PolynomialReduce → ResidueReduce` run literally (they already do at the
  recursion level — see §1);
- the level-integrator (`IntegratePrimitive` / `IntegrateHyperexponential` /
  `IntegrateHypertangent`) either returns `g` with `f − Dg ∈ k`, **or proves non-elementary**;
- its subproblem — **RischDE over `k`** (exp), **LimitedIntegrate over `k`** (primitive),
  **CoupledDESystem over `k`** (tangent) — is solved **literally and recursively over the
  tower field**, not by a bounded ansatz.

The one honest scope boundary Bronstein himself draws (§5.2 p.136): a **nonelementary
primitive monomial** or a **general nonlinear monomial with `Sⁱʳʳ ≠ ∅`** has no general
algorithm. Neither arises from a log/exp/tan tower, so both are **out of scope by
construction** for elementary integrands (Gap 5, document-only).

---

## 1. Current architecture (verified against source, 2026-07-14)

An independent source audit (not the self-authored docs) established:

**Literal / decision-procedure (correct-by-construction, exact `Together`-zero verify):**
- `rt_recursive_tower_case → rt_field_integrate → rt_field_ratint_hermite` /
  `rt_field_hyperexp_hermite`: literal `HermiteReduce` (`risch_hermite_reduce`) +
  residue-criterion LRT (`Integrate\`TranscendentalLogPart`) + hyperexp-poly Laurent.
  **The `RISCH_BRONSTEIN_GAPS.md §4 "ansatz wearing tower clothes" narrative is stale for
  this path.**
- Base-field `C(x)` RischDE `rde_base`: faithful Ch.6 port (`WeakNormalizer` →
  `RdeNormalDenominator` → `RdeBoundDegreeBase` → `SPDE` → `PolyRischDENoCancel1` / the
  `b=0` antidiff case), FLINT fast path, audited (A1/F5). Degree bounds `rt_rde_var_bound`
  exact + cap-free.
- Structure theorems (Cor 9.3.1 complex + 9.3.2 real) as a tested oracle; P3 decision
  (`Risch\`ElementaryIntegralQ` + `Integrate::nonelem`).

**Bounded ansatz + `SolveAlways` (declines rather than decides) — the completeness gaps:**
- **`rt_field_rde` general branch (L4122–4280):** the per-Laurent-coefficient RischDE over a
  tower field `Kₗ = C(x, t₀…t_{L-1})` is a `q = h/pd` polynomial-numerator ansatz solved by
  `SolveAlways`. **This is Gap 1 — the keystone.**
- `rt_log_tower_case` / `rt_exp_tower_case`: pure ansatz+diff-back flat cases, redundant with
  the recursive path, with local caps (`nl≤4`) and the hard `RT_MAXK=5` depth cap (Gap 4).
- `LimitedIntegrate` (§7.2) approximated by `SolveAlways` (Gap 2).
- No `RT_TAN` tower kind; tangents base-field-only (Gap 3).

---

## 2. Gap 1 — Recursive Risch DE over the tower (THE KEYSTONE)

### 2.1 Target

Replace the `rt_field_rde` general branch (L4122–4280) with a literal recursive Chapter-6
solver. `rt_int_hyperexp_poly` (L4063) already calls `rt_field_rde(pⱼ, ip, T, L, x)` per
Laurent power; that call site is unchanged — only the callee's non-base branch is rebuilt.

```
solve  D_tower[y] + f·y = g   for  y ∈ K_L = C(x, t₀…t_{L-1})    (f, g ∈ K_L)
```

### 2.2 The derivation abstraction (the crux)

The base stack hardcodes `D = d/dx` (`rde_dx`, `rt_degree(·,x)`, poly ops in `x`). The
recursive stack must operate on `k[τ]` where `τ = t_m` is the current top monomial and
`k = K_m = C(x, t₀…t_{m-1})` is treated as a **coefficient field** (lower tower vars are
transcendental constants for the `τ`-polynomial algebra). Introduce a small derivation
context so the same algorithm boxes run at every level:

```c
typedef struct {
    RtTower* T;      /* the tower                                           */
    long     m;      /* current top monomial index (τ = T->t[m]); m<0 = base C(x) */
    Expr*    x;      /* base variable                                        */
    /* derived: */
    Expr*    tau;    /* T->t[m]  (the polynomial variable)                   */
    RtKind   kind;   /* T->kind[m] : RT_LOG / RT_EXP / (RT_TAN, Gap 3)       */
    Expr*    Dtau;   /* D_tower[τ] = Dcoef_m·(1 | τ | τ²+1) per kind         */
    Expr*    Dlog;   /* T->Dcoef[m] : log u'/u ; exp w' ; tan a             */
} RdeCtx;
```

Primitive operations, each dispatching on `ctx->m < 0` (base → existing `rde_*`) vs.
recursive (poly-in-`τ` with lower-field coefficients):

| Op | Base (`m<0`) | Recursive (`k[τ]`) |
|---|---|---|
| `rde_deg(e, ctx)`   | `rt_degree(e, x)` | `rt_degree(e, τ)` |
| `rde_coeff(e,i,ctx)`| `rt_coeff(e, x, i)` | `rt_coeff(e, τ, i)` |
| `rde_D(e, ctx)`     | `d/dx` | `rt_tower_deriv(e, T, x)` restricted to `k[τ]` |
| `rde_Dk(c, ctx)`    | — | derivation of a **coefficient** `c ∈ k`: recurse `rde_D` at `m-1` |
| gcd / quot / rem / xgcd in `τ` | `PolynomialGCD[·,x]` etc. | over `k[τ]` — see §2.6 risk |

`Dtau` per kind: RT_LOG → `Dlog` (∈ k, `δ(τ)=0`, primitive); RT_EXP → `Dlog·τ`
(`Dτ/τ ∈ k`, hyperexponential); RT_TAN → `Dlog·(τ²+1)` (`Dτ/(τ²+1) ∈ k`, hypertangent,
Gap 3).

### 2.3 Function-by-function spec (Bronstein box → Mathilda function)

All new functions are file-local statics in `integrate_risch_transcendental.c` (they reuse
`rde_add/sub/mul/quot/rem/gcd`, `rt_tower_deriv`, `rt_rde_var_bound`, `rt_dec_nonelem`),
unless §4 moves the RDE stack to a `risch_rde.{c,h}` module.

- **`rde_tower(f, g, RdeCtx* ctx) -> Expr*` (NULL = no solution).** Top driver. `ctx->m<0` →
  `rde_base(f,g,x)`. Else: WeakNormalizer → RdeNormalDenominator → RdeSpecialDenom* →
  RdeBoundDegree* → SPDE → PolyRischDE(NoCancel|Cancel). Every NULL is authoritative →
  `rt_dec_nonelem()`.
- **`rde_weak_normalizer_field(f, ctx) -> Expr*` (§6.1 p.183).** Generalize
  `rde_weak_normalizer` to `k(τ)`: `SplitFactor` (via `Risch\`SplitFactor`) for the normal
  part, `Resultant_τ(a − z Dτ·(dτ), d₁)` for positive-integer residues. Common case → `1`.
- **`rde_normal_denominator_field(f, g, ctx) -> {a,b,c,h}` or "no solution" (§6.1 p.185).**
  `dₙ, eₙ = SplitFactor`-normal parts of `denom(f), denom(g)`; `h = gcd(eₙ, eₙ')/gcd(…)`;
  guard `eₙ | dₙh²` else no-solution (Cor 6.1.1(ii), authoritative). Mirrors the base logic
  already inlined in `rde_base` L1199–1231, lifted to `k(τ)`.
- **`rde_special_denominator(a, b, c, ctx) -> {ā,b̄,c̄,h}` (§6.2).** Dispatch:
  - RT_LOG (primitive p.188): identity, `k⟨τ⟩ = k[τ]`, `h=1`.
  - RT_EXP (`RdeSpecialDenomExp`, p.190): special `τᵏ`; compute the `τ`-adic valuations and
    clear the special part; may consult `ParametricLogarithmicDerivative` (Gap 2 dependency,
    heuristic ok initially).
  - RT_TAN (`RdeSpecialDenomTan`, p.192): special `τ²+1` (Gap 3).
- **`rde_bound_degree(a, b, c, ctx) -> long` (§6.3).** Reuse `rt_rde_var_bound(deg_τ(c/…),
  deg_τ(b/a), deriv_lowers=(kind==RT_LOG), m_res)`; the Prim cancellation branch needs a
  limited-integration test (`Dt = mη + Dz/z`) — heuristic/`m_res` ok initially, exact via
  Gap 2.
- **`rde_spde_field(a, b, c, ctx, n, RdeSpde* out) -> int` (§6.4 p.203).** Generalize
  `rde_spde`: same Rothstein degree-reducing recursion, but degree/xgcd/derivation via the
  `RdeCtx` ops (§2.2). `ExtendedEuclidean` over `k[τ]` (§2.6). Terminates when `deg_τ(a)=0`.
- **`rde_polyrischde_nocancel(b, c, ctx, n) -> Expr*` (§6.5 pp.208–210).** Top-down
  leading-coefficient matching; the "c ← c − Dp − bp" update uses `rde_D` (= `D_tower`).
  Three sub-cases by `deg_τ(b)` vs `δ(τ)`:
  - deg(b) large (p.208): pure iterative solve in `k[τ]`.
  - deg(b) small (p.209) / `δ(τ)≥2` degenerate (p.210): reduces the leftover to an RDE
    `Dy + b₀y = c₀` over `k` → **recurse `rde_tower(b₀, c₀, ctx@m-1)`**.
- **`rde_polyrischde_cancel(b, c, ctx, n) -> Expr*` (§6.6).** Cancellation cases:
  - `PolyRischDECancelPrim` (b∈k, p.212): recursively calls RischDE over `k` on the leading
    coefficient each pass → recurse `rde_tower(…, ctx@m-1)`.
  - `PolyRischDECancelExp` (p.213): analogous; the coefficient RischDE is over `k`.
  - `PolyRischDECancelTan` (p.215): Gap 3 (needs CoupledDESystem).

### 2.4 Increment breakdown (each lands green + committed)

- **1a — Scaffold + primitive (log) non-cancellation.** `RdeCtx`, `rde_tower`, the base
  passthrough, `rde_normal_denominator_field` + `rde_spde_field` + `rde_polyrischde_nocancel`
  for **RT_LOG** top only, with the coefficient recursion bottoming at `rde_base`. Wire
  `rt_field_rde`'s general branch to call `rde_tower` for the log case; keep the ansatz as a
  fallback initially (so nothing regresses), gated behind the new path.
- **1b — Exponential non-cancellation** (`RdeSpecialDenomExp` + RT_EXP `NoCancel`). This is
  the common tower case (nested exp Laurent coefficients).
- **1c — Cancellation cases** `PolyRischDECancelPrim` / `…Exp` (§6.6) — recurse into
  `rde_tower` over `k`. Closes the resonance/coupled-coefficient integrands the ansatz
  declined.
- **1d — WeakNormalizer + full RdeNormalDenominator** over `k(τ)` (spurious-residue strip);
  needed for correctness on integrands with positive-integer residues in the tower.
- **1e — Retire the `rt_field_rde` ansatz fallback.** Once 1a–1d cover the ansatz's domain
  (verified by the existing suite + new tests), delete L4122–4280's `SolveAlways` block; the
  general branch becomes `return rde_tower(i·Dcoef_L, p, ctx@L-1)`. Every NULL now
  authoritative.
- **1f — Decision wiring.** Confirm each `rde_tower` NULL flows to `rt_dec_nonelem()` so
  `ElementaryIntegralQ` decides the newly-covered class; add non-elementary siblings to the
  strict-unevaluated guard.

### 2.5 Tests (per increment)

New file `tests/test_risch_rde_tower.c` (register in `tests/CMakeLists.txt` COMMON_SRC +
target). Each increment adds:
- **Solve cases:** construct `D[Y] = F` for a known elementary `Y` at depth ≥2 whose
  coefficient RDE lands **off** the old ansatz box (rational lower-field coefficient,
  high τ-degree), assert `Simplify[D[∫F]−F]===0` (exact) — e.g.
  `∫ D[E^(E^x)/(1+Log[x]), x]`, and the deg-≥6 nested-exp Laurent analogues currently only
  reachable via the ansatz.
- **Decision cases:** a non-elementary sibling asserts `ElementaryIntegralQ→False` (RDE
  no-solution now authoritative) and `Integrate` leaves it unevaluated with `Integrate::nonelem`.
- **Regression:** the whole `integrate_risch_transcendental_tests` suite stays green; A/B any
  pre-existing failures at HEAD before attributing.
- **Invariant:** correct-by-construction — the recursive path emits only behind an exact
  `Together`-zero identity in the tower vars (same gate as `rt_field_ratint_hermite`); a
  mis-reduction can only decline. Valgrind: no new leak blocks over the `Sin[1.0]` baseline.

### 2.6 Risk register (Gap 1)

- **R1 — Polynomial algebra over `k[τ]` with a transcendental coefficient field.** `gcd`,
  `quot`, `rem`, `PolynomialExtendedGCD` must treat `τ` as the variable and `{x, t₀…t_{m-1}}`
  as field constants (rational functions). Mathilda's `PolynomialGCD[·, τ]` /
  `PolynomialExtendedGCD[·,·,τ]` over `Q(x, t…)` is the first choice; where it is too slow or
  incorrect over the multivariate coefficient field, route through the **FLINT tower engine**
  (`flint_bridge`, `project_flint_extension_engine`). **Mitigation:** 1a starts with the log
  case where coefficients are often polynomial; add a `Together`-normalize + exact-identity
  gate so any algebra slip declines, never ships wrong.
- **R2 — `RdeSpecialDenomExp` valuation logic** (τ-adic order at `τ=0`) is the subtlest box;
  build it against Bronstein Example 6.2.1/6.2.2 as unit tests before wiring.
- **R3 — Termination/perf** of the recursion at depth `RT_MAXK`; the base FLINT fast path
  keeps the leaves cheap. Keep the ansatz fallback until 1e proves parity.

---

## 3. Gaps 2–5 (after the keystone)

### Gap 2 — `LimitedIntegrate` first-class (§7.2)
`IntegratePrimitivePolynomial` (§5.8) needs a genuine limited-integration decision on the
leading coefficient (`f = Dv + Σ cᵢwᵢ`). Implement `LimitedIntegrateReduce` (p.248) + the
§7.1 parametric polynomial solve (`ParamRde` reduce → `LinearConstraints` → `ConstantSystem`).
Replaces the `SolveAlways` approximation in `rt_limited_field_integrate` (L3993). Also
sharpens `RdeBoundDegreePrim`'s cancellation test (Gap 1 R-mitigation) and
`ParametricLogarithmicDerivative` (§7.3). New: `src/calculus/risch_param.{c,h}`.
Tests: Bronstein Ex 7.2.x; `∫ Log[x]/(x+1)`-family primitives currently declined.

### Gap 3 — Tangent tower monomial (`RT_TAN`)
- Extend `RtKind` → `{RT_LOG, RT_EXP, RT_TAN}`; teach `rt_tower_build` to collect `Tan`/`Cot`
  (special `τ²+1`) and `Tanh`/`Coth` (special `τ²−1`) as tower monomials with
  `Dcoef = a` (`Dτ/(τ²+1)=a`).
- `rde_special_denominator` RT_TAN branch (`RdeSpecialDenomTan`, §6.2 p.192);
  `rde_polyrischde_cancel` RT_TAN branch (`PolyRischDECancelTan`, §6.6 p.215) → wire the
  existing `Risch\`CoupledDECancelTan` (`risch_coupled.c`) into the RDE stack.
- `rt_field_integrate` RT_TAN dispatch → `IntegrateHypertangent{Reduced,Polynomial}`
  (already exist in `risch_hypertangent.c`/`risch_coupled.c`) at every tower level.
- Closes tangents nested in deeper towers and `1/(a+b Sin[x])`.
Tests: extend `test_risch_hypertangent.c`; nested `∫ Tan[Log[x]]/x`-family.

### Gap 4 — Retire flat-ansatz cases + remove `RT_MAXK`
Once the recursive literal path (Gaps 1+3) subsumes them, delete `rt_log_tower_case` /
`rt_exp_tower_case` (or demote to a cheap fast-path with a comment), remove the `nl≤4`
local caps and the `#define RT_MAXK 5` depth cap (replace fixed arrays with dynamic
allocation in `RtTower`). The literal recursion is depth-agnostic. Diff the full suite.

### Gap 5 — §5.11 scope boundary (document-only)
`IntegrateNonLinearNoSpecial` and the nonelementary-primitive case have no general algorithm
(Bronstein §5.2 p.136) and never arise from a log/exp/tan tower. Add an explicit note to
`RISCH_STATUS.md` and a `Integrate` comment; no code. This is where the decision procedure's
completeness genuinely ends, and it is complete for **all elementary (Liouvillian
non-nonlinear) integrands**.

---

## 4. File layout decision

**Option A (minimal blast radius, chosen for 1a–1c):** grow the recursive RDE stack as
file-local statics in `integrate_risch_transcendental.c`, beside `rde_base`. Reuses every
helper + the decision hooks; no header churn.

**Option B (chosen at/around 1e, once the stack is large):** extract the whole RDE stack
(`rde_*`, `rde_tower`, `RdeCtx`, `RdeSpde`) into `src/calculus/risch_rde.{c,h}`, exposing
`rde_tower` + `rde_base`. Requires threading the decision context (pass `RtDecision*` instead
of the file-local `g_rt_decision`) and adding to `tests/CMakeLists.txt` COMMON_SRC. Do this
as a **behavior-preserving refactor commit** once the code justifies the module.

---

## 5. Sequencing & dependency graph

```
Gap 1 (keystone)  1a → 1b → 1c → 1d → 1e → 1f     ← start here
        │
        ├─ Gap 2 (LimitedIntegrate)  sharpens 1's RdeBoundDegreePrim + IntegratePrimitive
        ├─ Gap 3 (RT_TAN)            parallel; reuses 1's rde_tower + risch_coupled
        └─ Gap 4 (retire flat/caps)  after 1 (+3) subsume the ansatz cases
Gap 5  document-only, any time.
```

**Global invariants (non-negotiable, from `RISCH_STATUS.md §1):**
1. Recursive Risch, never parallel-Risch fallback.
2. Correct-by-construction: emit only behind an exact certificate (tower-var `Together`-zero
   identity); the heuristic diff-back is a backstop, not the gate, on the literal path.
3. Every NULL is either an authoritative "no elementary integral" (→ `rt_dec_nonelem`) or a
   dispatch-to-sibling decline — never a silent bounded-ansatz giveup once 1e lands.
4. `make -std=c99 -Wall -Wextra` clean; valgrind no new blocks over baseline; docs + weekly
   changelog updated per `CLAUDE.md`.

---

*Execution note: proceed continuously through 1a→1f, committing each green increment; do not
re-check scope between increments (the user has approved the plan). Re-plan only if an
increment reveals a structural obstacle (per `CLAUDE.md` "if something goes sideways, STOP and
re-plan").*
