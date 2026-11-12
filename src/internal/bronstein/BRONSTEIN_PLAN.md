# Plan: Bronstein's algorithm for mixed algebraic–transcendental integrands (simple radicals), plus Trager and Bronstein in the Automatic cascade

> **Execution note.** The first implementation step copies this file verbatim to
> `src/internal/bronstein/BRONSTEIN_PLAN.md`. That copy is the living plan. Tick its phase
> checkboxes, append a "Progress log" entry per session, and record every design revision
> in a dated section, as `docs/notes/TRAGER_PLAN.md` does. A fresh session starts by
> reading that file, then `tasks/lessons.md`, `SPEC.md` and `docs/notes/TRAGER_PLAN.md`.

## Context

Mathilda has three relevant pieces today:
- A complete decision procedure for **algebraic** integrands over one simple radical. This
  is Trager 1984, `src/internal/trager/`, run as `Method -> "Trager"` (v0.6081+).
- A **transcendental** Risch decision procedure, in C (`src/calculus/integrate_risch_*.c`,
  `risch_*.c`).
- The Risch–Norman heuristic **ParallelMixedTower** (`src/internal/mixed/ParallelMixed.m`)
  for mixed towers. It is incomplete: a failed ansatz gives no verdict.

Nothing *decides* a mixed integrand such as `(x+1)/((x Log[x]+x^2) Sqrt[x+Log[x]])` or
`Log[1+E^x]^(1/3)/(1+Log[1+E^x])`. Bronstein's 1990 paper, *Integration of Elementary
Functions* (JSC 9, 117–173), extends Trager's rational approach to `K = k(θ)(y)` with θ a
log or exp monomial. Restricted to n-th-root extensions it needs no Puiseux expansions
except explicit ones at the special places (§"Global approach", Prop 2.15, `S_log2`/`S_exp2`).
Trager's thesis, Appendix A (pp.73–80), gives the same reductions in elementary form for a
simple radical `y^n = P(θ)`.

**Goal.**
- Implement Bronstein's algorithm in Mathilda's own language under `src/internal/bronstein/`.
  Algebraic extensions are limited to one simple radical. It is exposed as
  `Integrate[f, x, Method -> "Bronstein"]` / `Integrate`Bronstein`.
- Put both Bronstein **and Trager** into the Automatic cascade (user decision, 2026-11-11).

**Scope (assumed; the scope question went unanswered, so this follows the recommendation).**
1. **Shape A, "radical on top"** (Phases 1–8): `K = k(θ)(y)`, `k = Q(x)(t_1..t_r)` a
   transcendental Log/Exp tower, θ = `t_{r+1}` the top monomial, `y^n = P ∈ k[θ]` with P
   actually depending on θ. This is Bronstein Part 2 together with Trager Appendix A. It
   needs only the transcendental Risch DE over k, which already exists in C.
2. **Shape B, "radical at the bottom"** (Phases 9–11): `K = Q(x)(y)(θ_1..θ_r)`,
   `y^n = g(x)`. This is transcendental Risch over an algebraic base. It needs a new Risch
   DE over `Q(x,y)` (Bronstein Part 3, pp.160–170) with explicit binomial Puiseux series.
3. **Out of scope** (clean decline):
   - nested radicals (Bronstein's `y⁴ − 2θy² + θ² − θ` example);
   - two non-Kummer-splittable radicals;
   - a radical in the middle of the tower with monomials above it that depend on it;
   - Tan/ArcTan kernels (Phase 12 may add a TrigToExp front end);
   - `η'` a Q-combination of logarithmic derivatives in the exponential special-place
     residue (Bronstein §"Linear combinations of logarithms", pp.155–157), deferred.

## Ground truth from exploration (read before coding)

**Papers.**
- Bronstein 1990: the PDF page is the journal page − 116 (p.117 is PDF p.1).
- Trager 1984: the printed page is the PDF page − 1.

Sections used:

| Topic | Reference |
|---|---|
| Normal/special places | [B] Def 1.3–1.4, Lemma 1.5–1.8 pp.121–124. Log case: every squarefree P is normal, and the only special place is θ=∞. Exp case: θ | P is special, so θ=0 and θ=∞ are special |
| Integral basis | [B] p.128. Radical basis `w_i = y^i/d_i` ([T] p.30); diagonal, hence normal at every place, ∞ included (Prop 1.14 iv) |
| Canonical representation, `f⁻` | [B] p.130 |
| Hermite at normal places | [B] pp.132–134, system (S), Lemma 2.1, Prop 2.2. Simple radical: [T] App A §2.1–2.2 pp.76 |
| Polynomial part / special places | [T] App A §2.3 (log, eqs 3–5; `d'≠0` subcase p.78) and §2.4 (exp, eq 6 p.79); [B] "Reduction at special places" pp.150–152 (`S_exp2`), Prop 2.18 |
| Residues at normal places | [B] Def 2.3, Prop 2.4–2.6: `R(z) = res_θ(prim_z(res_y(G − zHD', F)), D)` |
| Resultant criterion (non-elementary certificate) | [B] Prop 2.9: R must be `α·P(z)` with P of constant coefficients. Examples p.138–139 |
| Residues at special places | Log: Def 2.10 + Prop 2.11 (coefficient of θ⁻¹ in the Puiseux expansion at ∞, solved by `(E_0)`). Exp: Def 2.21, Lemma 2.19–2.20, Prop 2.22 |
| Divisors and log part | Log: Prop 2.12–2.13 pp.141–142. Exp: Prop 2.23–2.24 |
| Algebraic part, global approach | [B] Prop 2.15, `S_log2` p.146; `S_exp2` p.152. For an n-th root `M_ij = 0` (i≠j), so the coordinates decouple into Risch DEs |
| Torsion | [B] p.129–130 (good reduction, Prop 1.16–1.17); [T] Ch 6 pp.64–70 |
| Risch DE over an algebraic field (Shape B) | [B] Part 3 pp.160–170: Prop 3.3–3.5 bounds, `(R_2)`, Prop 3.7, Lemma 3.8, Puiseux-based solve |
| Component split | [T] App A §1 p.74: f is integrable iff each `g_i y^i` is; the rational part of `∫g_i y^i` is `v_i y^i` |

**Worked examples (use them as tests, exact answers).**
- [B] p.134: `∫(x+1)/((x Log[x]+x^2) Sqrt[x+Log[x]]) = −2/Sqrt[x+Log[x]]`. Hermite only.
- [B] p.139: `Log[1+E^x]^(1/3)/(1+Log[1+E^x])` is **not elementary** (resultant criterion,
  tower `x, γ=e^x, θ=log(1+γ)`, `y³=θ`).
- [B] p.147 (**corrected 2026-11-12**; the surface form first written here was a mis-transcription that
  does not match the answer): use the paper's θ-form
  `(((3x+1)Log[x] − x³ + x²)Sqrt[x+Log x] − (2x² − x − 1)Log[x] − 2x³ + x² + x)/(x Log[x]² − (x³ − 2x²)Log[x] − x⁴ + x³)`
  `= 2 Sqrt[x+Log x] + 2 Log[x + Sqrt[x+Log x]]`. Log case: residue at ∞, divisor `p₂ − p_∞`.
- [B] p.158: `(3(x+E^x)^(1/3) + (2x²+3x)E^x + 5x²)/(x (x+E^x)^(1/3)) = 3x(x+E^x)^(2/3) + 3 Log[x]`.
  Exp case: special-place reduction via `S_exp2`.
- [B] p.138: the nested example `Log[x]+Sqrt[Log x + Sqrt[Log x]]` must **decline** (nested
  radical).

**Mathilda pieces to reuse.**

`src/internal/trager/` (`BeginPackage["Trager`"]`, private functions):
- Field layer: `fieldQ/fieldNF/fieldGF`, `kc`, `kDivRem`, `kExtGCD`, `kGCD`, `kRowEchelon`
  (`Field.m`).
- `curveBasis` (`Basis.m`); `af*` arithmetic (`Arithmetic.m`, where `afD` is hard-wired to
  `D[·,x]`).
- `hermiteReduce`/`hermiteSolve` (`Hermite.m`, `D[·,x]` hard-wired).
- `componentCurve`/`reducedCurveQ` (`Reduce.m`); `canonicalForm`/`radicalPowers` (`Validate.m`).
- `residueIdeal`/`residueData` (`Ideal.m`, which takes `Ep` as **data**, so it can be reused
  with `Ep = Risch`Derivation[E]`).
- `hnf`/`normalizeAtInfinity` (`Normalize.m`); `principalGenerator`/`principalQ`.
- `torsionOrderSpec`/`torsionScreen`/`goodPrimes` (`TorsionModP.m`).
- `shiftCurve`/`unshift` (`Shift.m`); `realForm`/`presentable` (`LogToArcTan.m`, `Pipeline.m`).
- Failure protocol `tragerFail`/`catchTrager`/`throwNonElementary`/`throwDecline`
  (`Common.m`).
- Algebraic modules (Field/Normalize/Ideal/PrincipalGen) are field-generic through the field
  record `F` and a polynomial-variable argument. Only `fieldQ`/`fieldNF`/`fieldGF` exist,
  and each dispatches to FLINT kernels (`nativeOpts`). **There is no rational-function
  coefficient field Q(x, t…)**, which Shape A needs: the curve is over `k = Q(x,…)` in the
  variable θ.

Risch builtins callable from `.m`. All are Protected and documented, and all take surface
expressions with Log/Exp kernels:
- `Risch`RischDE[f, g, x]` solves `y' + f y = g` over C(x) or the inferred tower
  (`integrate_risch_transcendental.c:745`). **It stays unevaluated BOTH on "no solution" and
  on "out of scope"**, so it cannot by itself back a non-elementary certificate.
- `Risch`SPDE`, `Risch`PolyRischDENoCancel` (base field only).
- `Risch`Derivation[p, {var -> Dvar,…}]`, `Risch`SplitFactor`, `Risch`SplitSquarefreeFactor`,
  `Risch`CanonicalRepresentation`, `Risch`PolyDivide`, `Risch`FieldGCD`,
  `Risch`ExtendedEuclidean`, `Risch`Diophantine`, `Risch`NormalQ`, `Risch`SpecialQ`
  (`risch_canonical.c:387+`).
- `Risch`HermiteReduce` (`risch_hermite.c`).
- `Risch`LogReducible`, `Risch`ExpReducible`, `Risch`RationalSpan`,
  `Risch`LogarithmicDerivativeOfRadical` (structure theorems, `risch_structure.c:461+`).
- `Integrate`RischTranscendental[f, x]`, `Risch`ElementaryIntegralQ[f, x]` (True/False/unevaluated).

**Integrate glue** (`src/calculus/integrate.c`):
- Trager lazy load and builtin at :1457–1530 (`trager_lazy_load`, `integrate_trager_common`,
  `try_trager`). `try_trager` declines on a List result, i.e. the
  `{"not elementary"|"failed", reason}` lists.
- Method enum :1749–1801 (`METHOD_TRAGER`…); `Integrate::method` message :3208; dispatch
  :3507; registration :3617.
- Cascade at :3330–3420, with stage profiling enum `IntegrateProfStage` :2680 (`IP_*`,
  `ip_stage_name`). The order is: … `IP_RISCHTRANS` → `IP_CRCTABLE` → `IP_POLYLOG` →
  `IP_HOLONOMIC` → `IP_BESSELHYPER` → `IP_LOGBYPARTS` → `IP_PMT` (gated
  `!g_integrate_no_special && !has_pseudoelliptic_radical`) → `IP_GOURSAT` → `IP_PMS` (last;
  the only special-function stage) → `IP_ROACH`.
- Non-elementary announcement: `integrate_announce_nonelementary` (used at :1297 for the PMT
  certificate, `pmt_is_nonelementary_certificate` :1170).
- PMT cascade budgets pattern: `pmt_in_cascade`, and `$ParallelMixedTimeBudget` Block'd by
  the C side (see the ParallelMixed.m header :100–110).

**Tests.**
- `tests/trager_test_utils.h`: `tr_init()` loads init.m and the package and puts
  `Trager`Private`` on `$ContextPath`; `tr_eq`, `TR_ANTIDERIV`. Suites
  `tests/test_trager_*.c` are registered explicitly at `tests/CMakeLists.txt:6480+`.
- FriCAS harness: `tools/trager_fricas_stress.py` (nightly ctest `trager_fricas_stress`);
  FriCAS at `/opt/local/bin/fricas`.

**Memories to honour.**
- `.m` package context path: `$ContextPath={"System`"}` + Begin private, else `a_` captures.
- A new Integrate method lives in `.m`, with C only for lazy-load wiring.
- No WL templates in C.
- Cascade ordering (`feedback_integrate_cascade_ordering`): a new cascade stage costs what it
  declines.
- A budget must not stop the cascade.
- Decline-with-fallback is a silent wrong answer.
- PossibleZeroQ is the verifier: no numeric fallback, and a PZQ false negative gets fixed in
  PZQ.
- Version race: check tags right before a bump.
- No reference-system product name in source or commits.

## Design

### Key simplification (both shapes)

[T] App A §1 splits `f = Σ a_i y^i` into components. Each `a_i y^i` is integrated on its own
reduced curve `Y^m = G` with `Y = y^i/s`; reuse the logic of `Trager`Private`componentCurve`,
lifted to `k[θ]`. On that curve the integrand sits in **one coordinate**, `c·Y` with
`c ∈ k(θ)`, `w_1 = Y`, `d_1 = 1`. Then:
- The algebraic part is `b·Y`, and `D(bY) = (Db + s b) Y` with
  `s = (1/m) DG/G` (D the tower derivation).
- So the whole algebraic-part question, `S_log2`/`S_exp2` restricted to the one live
  coordinate, is the **Risch DE `Db + s b = c̃` over k(θ)**.
- Its solution is unique, because the homogeneous solution `1/Y ∉ k(θ)`.
- This replaces Bronstein's Puiseux-based Prop 2.14/2.17 algebraic-part systems entirely.

Hermite/special-place reductions are still run first, because they produce certificates
early and shrink the residue computation.

### Per-component pipeline, Shape A (Phases 2–7)

Input: `c Y`, `Y^m = G(θ)`, G ∈ k[θ] with every θ-irreducible multiplicity < m, plus the
tower record `T`.

1. **Hermite at normal finite places** ([B] (S), [T] App A §2.1–2.2).
   - Split the denominator with `Risch`SplitFactor` into normal × special.
   - For each normal squarefree factor `V^j`, j ≥ 2, solve
     `A ≡ B·U·(V s − (j−1) DV) (mod V)` over k[θ].
   - Invertibility mod V holds for the reasons in the header of trager `Hermite.m`, with
     `'` replaced by D.
2. **Special places.**
   - *Exp*, θ=0 ([B] `S_exp2` p.152; [T] eq. before (6), p.79): to drop `θ^k`, solve the
     Risch DE in **k**, `b₀' + (s₀ − kη') b₀ = c₀` (`s₀` the θ⁰ part).
   - *Exp*, polynomial part / θ=∞ ([T] eq. (6)): top-down degree reduction, one Risch DE in
     k per degree.
   - *Log*, polynomial part ([T] eqs. (3)–(5)): the leading-coefficient equation
     `c_{j+m} = b_{j+1}·w + b_j'`, with `b_{j+1}` constant, is **limited integration in k**
     against the single log θ.
   - *Log*, `d' ≠ 0` subcase (the radicand content `d ∈ k` is non-constant, [T] p.78): the
     coefficient equation `c_{k+m} = b_k' − (d'/(n d)) b_k` is a Risch DE in k.
   - A no-solution answer is a **certificate** (Prop 2.18 / [T]) **only if the Risch DE
     reported a certified no-solution**. Otherwise decline.
3. **Residues at normal places** (Prop 2.6):
   `R(z) = res_θ(N^m G − z^m E'^m, E)` up to the `prim_z` content, with `E' = DE`. This
   generalises trager `componentResiduePolynomial`.
   - Resultant criterion (Prop 2.9): when the monic primitive part of R in z does not have
     constant coefficients (free of x, t_i, θ), the result is a NonElementary certificate.
4. **Residues at special places.**
   - *Log*, θ=∞ (Def 2.10): take `a` = the coefficient of θ⁻¹ in the Puiseux expansion of
     `cY` at each place over ∞. For `Y^m = G` this is explicit:
     `Y = lc(G)^{1/m} θ^{deg G/m} (G/(lc θ^deg))^{1/m}`, a binomial series in `θ^{-1}`;
     there are `gcd(m, deg G)` places with ramification `r = m/gcd`.
   - Then solve `a = v' + c η'/η` (limited integration in k) and set `τ = −r c`.
     No solution ⇒ certificate (Prop 2.11).
   - *Exp*, θ=0 and θ=∞ (Def 2.21): values `v_p(f)` from the explicit Laurent expansions.
     `τ` comes from `v_p(f) − v_{p̄}(f) = cη' + Σ c_i d_i'/d_i`.
   - Decline when η' is a Q-combination of log derivatives (`Risch`LogReducible` /
     `Risch`RationalSpan` decides that).
5. **Q-basis and divisors** (Prop 2.13/2.24).
   - Take a Q-basis `r_1..r_l` of the τ's. Set `D_i = Σ_p q_ip p` over places of `K/k`.
   - Finite places: trager `residueIdeal` over the coefficient field `k(α)` (the new `Fn`
     field type, below).
   - Places over θ=∞ (and θ=0 in the exp case): move them to finite places with the purely
     algebraic shift `θ = a + 1/z` (trager `Shift.m`; divisors do not see the derivation),
     then use the same ideal code.
6. **Torsion by specialization**, a new step. The divisor lives on a curve over the function
   field `k = Q(x,t…)`, where the mod-p Weil bound does not apply directly ([B] p.130 notes
   the problem).
   - Specialise every transcendental of k (x, t_1..t_r, treated as independent
     indeterminates) to random small rationals. The specialization must be good: G's θ-factor
     degrees and multiplicities, genus, distinctness of the residue places, and nonvanishing
     of denominators and leading coefficients are all preserved.
   - On the specialised curve over `Q(α)`, run Trager's `torsionOrderSpec`/`torsionScreen`.
   - Soundness: in characteristic 0 specialization of a good model is injective on torsion,
     so the generic order is either `N_spec` or ∞.
   - Non-torsion at a good specialization ⇒ **NonElementary certificate**.
   - Otherwise test principality of `N_spec·D_i` **generically** over `k(α)[θ]` (trager
     `principalGenerator` with field `Fn`). Not principal ⇒ certificate; principal ⇒ the
     generator `w_i`.
7. **Log part and algebraic part.**
   - Log part: `L = (1/h) Σ r_i log w_i`. Verify the coordinate-1 identity exactly, as trager
     `logSumResidual` does.
   - Residual `c̃ Y = cY − (coordinate 1 of L')`: solve `Db + s b = c̃` over k(θ) by the
     certified Risch DE (tower k(θ) passed as surface kernels).
   - Coordinate-0 residual: it must be `Σ d_j v_j'/v_j` with `v_j ∈ C̄k` (Prop 2.13).
     Integrate it with `Integrate`RischTranscendental` and accept it only if the result is
     logs of k-elements; otherwise decline.
8. **Reassemble.**
   - Substitute θ → surface kernel and Y → the user's radical (trager `componentExpression`,
     `reduceRadicalPowers`).
   - Final `PossibleZeroQ[D[F,x] − f]`. On failure throw an Internal error and decline.
   - Real form via trager `presentable`/`realForm`.

### Front end (`Tower.m`)

- Collect Log/Exp kernels innermost-first and decide reducibility with
  `Risch`LogReducible`/`Risch`ExpReducible`; rewrite reducible kernels.
- Build the derivation rules `{t_i -> Dt_i}`.
- Detect radicals with trager `radicalPowers`/`detectRadical` (the B11 unification of
  `Power[_, 1/2]` and `x^(1/4)`). For a compositum, apply `KummerSplit` first, as trager
  `Pipeline.m:228` does.
- Classify: no radical → hand to `Integrate`RischTranscendental` (return its result
  verified); no kernel → Trager; radicand depends on the top monomial → **Shape A**; radicand
  in Q(x) with kernels above → **Shape B** (Phase 9+; decline until then); otherwise decline.
- Flatten a radicand linear (or Möbius) in x, as `(a x + b)^(1/n)`: substitute
  `x = (u^n − b)/a`, giving a purely transcendental tower in u. Hand to
  RischTranscendental and back-substitute. This is a cheap, large-coverage pre-pass and a
  sound identity.
- θ is reordered to be the top monomial on which the radicand depends. If the radicand
  depends on several monomials, θ is the highest one, and the lower ones belong to k.

### Coefficient field `Fn` (an edit to `src/internal/trager/Field.m`)

- Add `fieldFn[vars, ext]`: rational functions in `vars` (x, t_i) over Q or `Q(α)`.
  - `kc` = Together (with `Extension` for `Q(α)`), `kZeroQ`, `kDivRem` =
    `PolynomialQuotientRemainder` in the polynomial variable with Together'd coefficients,
    `kExtGCD` = `PolynomialExtendedGCD`, linear algebra via RowReduce.
  - No FLINT path: `nativeQ` is False.
- Everything else in Normalize/Ideal/PrincipalGen stays untouched.
- Spike first (Phase 1): run `hnf`, `normalizeAtInfinity`, `principalGenerator` on a genus-0
  curve over Q(x) in θ, and confirm `Fn` is enough. If a module hard-codes Q, generalise that
  module in place (backward compatible; the Trager suites gate it). **Do not fork the trager
  modules.**

### Package layout: `src/internal/bronstein/`

Every file gets a header block (purpose, [B]/[T] section and page references) and a
`::usage` docstring on every function. Function names carry a `br` prefix, so they never
collide with Trager's.

| File | Contents |
|---|---|
| `Bronstein.m` | `LoadModule["trager/Trager.m"]`; `BeginPackage["Bronstein`", {"Mathilda`"}]`; public `BronsteinIntegrate::usage`; `Begin["`Private`"]`; append `"Trager`Private`"` to `$ContextPath` (verify in Phase 1 that short trager names resolve; otherwise use fully qualified names); `LoadModule` siblings in order |
| `Common.m` | Failure protocol. Reuse Trager's `tragerFail` tags `NonElementary`/`Decline`/`Internal`; `brOption[...]` zero-arg config (budgets). No top-level variables (they land in Global`) |
| `Tower.m` | Surface → tower record `<|"x","kernels","rules","theta","kind"("Log"|"Exp"),"eta","lower"|>`; shape classification; linear-radicand flattening |
| `Poly.m` | Polynomial ops in θ over k: content/primitive part in k, squarefree and normal/special split via `Risch`SplitFactor`/`SplitSquarefreeFactor`, the derivation `brD` = `Risch`Derivation` |
| `Curve.m` | Component split over k[θ] (lifting `componentCurve`, including the content `d ∈ k` and the gcd-of-exponents reduction); `brCurve` record `{m, G, d, s}`; regularity / absolute irreducibility ([T] Ch 3 §3, with the `d'≠0` "no constant s with s·d an n-th power" check, [T] p.78) |
| `RDE.m` | `brRischDE[f, g, T]` → `{"Solution", y}` / `{"NoSolution"}` / `{"Unknown"}`, wrapping the certified C entry point below |
| `LimitedIntegrate.m` | `brLimitedIntegrate[a, w, T]`: find a constant λ and `v ∈ k` with `a = v' + λ w` (specialised to one log). It uses `Integrate`RischTranscendental` plus `Risch`LogReducible` normalisation of log arguments; certified "no" only from certified sub-results |
| `Hermite.m` | Normal-place Hermite in coordinate 1 (generalised `hermiteSolve` with D) |
| `Special.m` | Log polynomial-part reduction, exp θ=0 and polynomial reduction, the `d'≠0` subcase |
| `Puiseux.m` | Explicit binomial expansions of `cY` at the places over θ=∞ (and θ=0 for exp); ramification indices; place enumeration as `ω`-orbits |
| `Residues.m` | `R(z)`, the constant-coefficient certificate, finite residues, special-place τ |
| `Divisors.m` | Q-basis of τ, divisor ideals over `fieldFn[...]` (finite places via `residueIdeal`, special places via shift) |
| `Specialize.m` | Good-specialization selection and checks; torsion via Trager's TorsionModP on the specialised curve; generic principality of `N·D` |
| `LogPart.m` | Assembly of `(1/h) Σ r_i log w_i`, exact coordinate-1 verification |
| `AlgebraicPart.m` | `Db + s b = c̃` over k(θ); coordinate-0 residual |
| `Pipeline.m` | `BronsteinIntegrate[f, x]` → antiderivative / `{"not elementary", reason}` / `{"failed", reason}`; component loop; reassembly; final PZQ check; `presentable` |
| `RDEAlgebraic.m` (Phase 10) | Risch DE over `Q(x)(y)`, [B] Part 3 |
| `BaseAlgebraic.m` (Phase 11) | Shape B driver |

`makefile` install already copies `src/internal/**/*.m`; check that the subdirectory is
picked up.

### C changes (glue only)

1. `Risch`RischDECertified[f, g, x]` in `integrate_risch_transcendental.c` returns `y`,
   `"NoSolution"` (only on a path that proves no solution: SPDE `$Failed` / bound
   exhaustion), or `$Failed` (out of scope or unsupported).
   - Audit `rde_solve_tower`/`rde_base` to separate those paths (Phase 1). This is a small
     tri-state wrapper, not new algorithm code.
   - Docstring, Protected, a test in `tests/test_risch_rde_tower.c`.
2. `integrate.c`:
   - Add `METHOD_BRONSTEIN`, `"Bronstein"` in `method_from_string` and in the
     `Integrate::method` string, plus dispatch.
   - Add `bronstein_lazy_load()` (`mathilda_load_module("bronstein/Bronstein.m")`) and
     builtin `Integrate`Bronstein` (Protected, docstring, `tc_async_defer_push`,
     `mth_msg_suppress_push`, cloned from `integrate_trager_common`).
   - Add `try_bronstein(f, x, cascade)`.
   - Integrate docstring method list.
3. **Cascade** (user decision): two new stages, `IP_TRAGER` ("15b Trager") and
   `IP_BRONSTEIN` ("15c Bronstein"), between `IP_LOGBYPARTS` and `IP_PMT`, skipped under
   `g_integrate_no_special`.
   - Gates, cheap and structural:
     - Trager: exactly one radical class after Kummer, no Log/Exp/trig kernels, not
       flattenable by the linear-radical stages already run.
     - Bronstein: ≥ 1 non-flattenable radical and ≥ 1 Log/Exp kernel of x, no trig/foreign
       kernels.
   - Budgets: `$TragerCascadeBudget` / `$BronsteinCascadeBudget`, public zero-arg-safe
     symbols Block'd by C in the cascade, `Infinity` for an explicit Method. Initial values
     are set by measurement (Phase 8). A budget expiry is a decline that lets the cascade
     continue.
   - **A certified non-elementary verdict** from either stage calls
     `integrate_announce_nonelementary`, as PMT does at :1297, and sets a cascade flag that
     skips the remaining **elementary-only** stages (PMT, Goursat). PMS (special functions)
     and Roach still run, since an elliptic/Erf answer is still wanted.
   - Explicit `Method -> "Trager"/"Bronstein"` keeps today's behaviour: returns
     unevaluated + `Integrate::nonelem` on a certificate.

### Testing policy (user requirement: extensive unit tests for EVERY routine, .m and .c)

- **Every `.m` function** written or modified in this campaign has its own unit tests in
  the suite for its module. That covers every `br*` function in `src/internal/bronstein/`
  and every Trager function generalised or added (e.g. the `fieldFn` arms of `kc`,
  `kDivRem`, `kExtGCD`, `kRowEchelon`). "Own tests" means:
  - ≥ 3 direct calls on hand-checked inputs, covering the typical case, the edge cases
    (degree 0, constant content, m = 2 vs m ≥ 3, log vs exp, empty or zero input) and the
    decline/throw paths, each asserted on its exact failure tag;
  - plus a **property test** where an invariant exists (e.g. `D[alg] + rem == f`, a
    Bezout identity, `R(z)` roots equal the residues at sampled places, a divisor of a
    generator equal to the requested divisor, specialization commuting with arithmetic),
    run over a seeded random family (≥ 20 instances).
- **Every C routine** added or changed has a direct test:
  - `Risch`RischDECertified`: each of the three outcomes, plus agreement with
    `Risch`RischDE` on a corpus.
  - `bronstein_lazy_load`, `integrate_bronstein_common`, `try_bronstein`: decline rules,
    message routing, `Quiet`/`Check`.
  - The cascade gates are exposed for testing as a Protected internal builtin
    `Integrate`CascadeGateQ["Trager"|"Bronstein", f, x]`, so each gate is unit-tested
    directly on ≥ 30 accept/reject integrands, not only through the cascade.
  - The budget Block and the nonelem skip flag.
- **Coverage gate.** `tools/check_bronstein_test_coverage.py`, run by `make
  check-bronstein-coverage` and also as a push-tier ctest:
  - it parses every `name::usage` in `src/internal/bronstein/*.m` (and a listed set of
    trager names touched by this work) and every new C builtin name;
  - it fails when a name is not called in some `tests/test_bronstein_*.c` /
    `tests/test_trager_*.c` / `tests/test_risch_*.c`;
  - assert-empty, with an `EXEMPT` list (one-line reason each) only for trivial accessors.
- A routine is not done until its tests are in, green and covered by the gate. Each phase's
  checkbox below includes its routines' tests.

### Tests (an extensive suite per module, all on the push tier)

- `tests/bronstein_test_utils.h`, cloned from `trager_test_utils.h` (`br_init` loads the
  package, puts `Bronstein`Private`` on `$ContextPath`; hard asserts; `BR_ANTIDERIV` via
  PossibleZeroQ).
- Suites, each registered explicitly in `tests/CMakeLists.txt` next to the trager block:

| Suite | Coverage |
|---|---|
| `test_bronstein_tower.c` | kernel collection, reducibility rewrites (`Log[x^2]`, `E^(2x)` vs `E^x`), shape classification of 30+ integrands, linear-radicand flattening identity |
| `test_bronstein_field_fn.c` | `fieldFn` ops; HNF / normalize-at-∞ / principal generator over Q(x) in θ (property: same module, Bezout transform certified) |
| `test_bronstein_curve.c` | component split over k[θ], content and `d'≠0` detection, basis `s`, irreducibility |
| `test_bronstein_rde.c` | tri-state wrapper: solutions, certified no-solutions, `Unknown` on out-of-scope |
| `test_bronstein_hermite_special.c` | [B] p.134 exactly; `D[alg] + remainder == f`; exp θ=0 and polynomial reductions ([B] p.158 system `(S_exp2)`, solution `(0,0,3x)`); log polynomial part; Prop 2.18 certificate |
| `test_bronstein_residues.c` | `R(z)` of [B] p.139 (`−z³ − (1+e^x)³/e^{3x}` → certificate), p.147 (`P(x) z³ (z−2)`), p.158 (`P(x) z³`); Puiseux at ∞ (log) with τ_∞ = −2 for p.147; exp special τ |
| `test_bronstein_divisors_torsion.c` | divisor `p₂ − p_∞` of p.147 principal with generator `x + y`; specialization goodness checks and rejections; a non-torsion example certified; candidate-order verification |
| `test_bronstein_integrate.c` | End to end through `Integrate[…, Method -> "Bronstein"]`: the four paper examples; families `E^x Sqrt[1+E^x]`, `E^x/Sqrt[1+E^(2x)]`, `Sqrt[1+E^x]`, `1/(x Sqrt[Log x])`, `Log[x]/(x Sqrt[1+Log[x]])`, `(x+E^x)^(1/3)` variants. Certificates: `Sqrt[Log[x]]`, `1/Sqrt[Log[x]]`, `Sqrt[E^x+x]`-type. Declines: nested radical, two radicals, Tan kernel. Silent unevaluated plus `Integrate::nonelem` capture; `Quiet`/`Check` behaviour; definite pass-through |
| `test_integrate_cascade_trager_bronstein.c` | Automatic cascade now closes Trager-only cases (e.g. `1/(x Sqrt[1-x^3])`, `(x^4+1)^(-1/4)`) and Bronstein cases. The nonelem certificate skips PMT/Goursat but not PMS. Budget expiry continues the cascade |

**FriCAS stress** (nightly):
- `tools/bronstein_fricas_stress.py`, modelled on `tools/trager_fricas_stress.py`.
- Corpus `tests/data/bronstein_stress_corpus.m`: ~250 integrands. Hand families over
  `θ ∈ {Log x, Log(1+x), E^x, E^(x²)}` with radicands linear/quadratic in θ, n = 2..4, plus
  a seeded random generator `R(x,θ)·y^k`.
- Classes AGREE-ELEM, AGREE-NONELEM, MATHILDA-ONLY, FRICAS-ONLY, DISAGREE (hard fail),
  WRONG (hard fail), TIMEOUT.
- A `BASELINE` ratchet on FRICAS-ONLY; ctest `bronstein_fricas_stress` (nightly by name).

### Non-trivial integral generator and benchmark (user requirement)

**Goal.** Machine-generate integrals that Bronstein's (or Trager's) algorithm computes and
that are **non-trivial**. An integral is non-trivial if at least one of these holds; record
which:
- **(N1)** No obvious substitution or by-parts closes it.
- **(N2)** The antiderivative needs algebraic numbers absent from the integrand: its
  constant field is a proper extension of the integrand's.
- **(N3)** FriCAS cannot compute it (it fails, returns unevaluated, or times out).

**Construction: by differentiation, so every item ships with a verified reference
antiderivative.** The generator is `tools/integral_generator/generate.m`, run by
`./Mathilda -file` (no wolframscript), seeded and deterministic.

1. Pick a shape: Trager (`y^n = g(x)`); Bronstein Shape A (`y^n = P(θ)`,
   `θ ∈ {Log[a(x)], E^(b(x))}`, optionally over a lower monomial); Shape B later. Pick
   n ∈ {2,3,4,6} and small random polynomials with the radicand m-th-power-free. Keep the
   genus mix deliberate, with genus 0 and genus ≥ 1 both present.
2. Build `F = v + Σ_α Σ_j c_j(α) Log[u_j(α)]`:
   - **Algebraic part:** `v = b(θ)·y^i` with random `b ∈ k(θ)`.
   - **Log part (N2 by design):** take an irreducible `q(z) ∈ Q[z]` of degree 2–4 and,
     for its roots α, logands `u(α) ∈ Q(α)(x,θ)[y]`. Then
     `F_log = RootSum[q, α ↦ α Log[u(α)]]`; its derivative is a trace, rational over Q, so
     the integrand carries no α while the answer needs Q(α). Variants: conjugate pairs
     giving ArcTan; logands `p + q·y` from polynomial **Pell / norm equations**
     `p² − G q² = c·(monomial)` over `k[θ]`, which give genus ≥ 1 torsion divisors (the
     pseudo-elliptic class), searched with the Trager torsion machinery; logands whose
     divisor has torsion order 2–6.
3. `f = D[F, x]`, then **disguise it**: Together, reduce the powers of y modulo `y^n`,
   expand the numerator and denominator, and rationalise the denominator into the
   `(A + B y)/(C + D y)` surface form, so the structure of F is not visible. Reject the item
   if f's LeafCount exceeds a cap, so the integrands stay human-sized.
4. **Filter.**
   - **(N1):** run every cheap cascade stage up to and including `LogByParts` (everything
     before the Trager/Bronstein stages, through an `Integrate`StageOnly[stage, f, x]`
     internal hook or the explicit `Method ->` of each). If any closes it, it is trivial
     for N1.
     - Also try the substitution heuristics directly: `u = θ`, `u = y`, `u = radicand`
       (DerivativeDivides), and one by-parts split with `u = Log`-factor.
   - **(N2):** compare the constant fields: the algebraic numbers in the minimal
     RootSum/Log form of F against those in f (`ToNumberField` / free `Root`/radical
     constants).
   - **(N3):** FriCAS pipe (`fricas -nosman`, the `tools/trager_fricas_stress.py`
     pattern) with a timeout; record "fails / unevaluated / timeout / closes".
   - Keep an item if (N1 ∨ N2 ∨ N3), **and** `Method -> "Bronstein"`/`"Trager"` closes
     it, **and** the result passes `PossibleZeroQ[D[F_result, x] − f]`.
   - Items Mathilda does not yet close are kept in a separate `pending` list. They are
     targets, not failures.
5. **Optional certified class:** also emit items Bronstein/Trager **prove non-elementary**
   (perturb a generated f by a holomorphic differential or a non-constant residue), tagged
   NONELEM. They are cross-checked against FriCAS; a disagreement is a hard failure.

**Outputs.**
- Frozen corpus `tests/data/generated_integrals.m`: records
  `{id, shape, genus, n, f, Fref, tags -> {N1,N2,N3}, fricas -> status}`, regenerated only
  deliberately (`make gen-integrals SEED=…`).
- A generator unit suite `test_integral_generator.m`/`.c`: determinism under seed; every
  record's `D[Fref] == f`; disguise is identity-preserving; each filter verified on
  hand-made trivial and non-trivial items.

**Benchmark.** `benchmarks/256-generated-algebraic-integrals/`:
- `generated_integrals.m`: Mathilda, Automatic and the explicit methods, timed, each
  answer PZQ-checked.
- `generated_integrals.input`: FriCAS, picked up automatically by `run_all.py` (it already
  supports `.input` halves).
- `generated_integrals.py`: the sympy baseline.
- `README.md`: the generation method, the N1/N2/N3 definitions, and per-tag counts.
- Rows are labelled by corpus id and tag. `make bench-gap` reports Mathilda vs FriCAS vs
  sympy: closes, time, and N3 items closed only by Mathilda.
- A nightly ctest `generated_integrals_corpus` asserts zero WRONG over the frozen corpus,
  plus a ratcheting `BASELINE` of closed counts per tag.

## Phases (continuous execution; tests gate each phase; commit and bump per substantive phase)

- [x] **0. Bookkeeping.** Copy this plan to `src/internal/bronstein/BRONSTEIN_PLAN.md`. Add
  a `tasks/todo.md` section. Read `tasks/lessons.md` and the memories named above.
- [x] **1. Spikes and substrate.**
  - (a) `fieldFn` in trager `Field.m` plus the HNF / principality spike over Q(x).
  - (b) `$ContextPath` bridging from `Bronstein`Private`` to `Trager`Private``.
  - (c) Audit `Risch`RischDE` paths and add `Risch`RischDECertified`.
  - (d) Confirm `Risch`SplitFactor`/`Derivation` on tower expressions with θ as a symbol.
  - Suites `field_fn`, `rde`. Commit, v+0.0001.
- [x] **2. Tower and Curve.** `Tower.m`, `Poly.m`, `Curve.m`, linear-radicand flattening.
  Suites `tower`, `curve`.
- [x] **3. Hermite and special places.** `Hermite.m`, `Special.m`, `LimitedIntegrate.m`.
  Suite `hermite_special`. [B] p.134 and p.158 close end to end through a temporary driver.
- [x] **4. Residues.** `Puiseux.m`, `Residues.m` (normal, log-∞, exp-0/∞). Suite `residues`.
  The p.139 certificate.
- [x] **5. Divisors and torsion.** `Divisors.m`, `Specialize.m`. Suite `divisors_torsion`.
- [x] **6. Log and algebraic parts.** `LogPart.m`, `AlgebraicPart.m`. The p.147 example
  closes.
- [x] **7. Pipeline and C method glue.** `Pipeline.m`, `Bronstein.m`; `Method -> "Bronstein"`,
  `Integrate`Bronstein`; suite `integrate`. Docs: `docs/spec/builtins/calculus.md` (Method
  list plus an `Integrate`Bronstein` subsection: algorithm, scope, decline rules, examples),
  the changelog `docs/spec/changelog/2026-11-09.md` (or the current week's), the Integrate
  docstring. Commit, bump, tag.
- [x] **8. Cascade.**
  - Trager and Bronstein stages, gates, budgets, the nonelem skip flag; suite `cascade`.
  - **Measure**: stage profile (`IP_*` timings) over `ctest -R integrate` and the
    integration corpora, before and after. The decline cost per non-matching integrand must
    stay small; set the budgets from the data.
  - `ctest -LE nightly` stays green, apart from expectation updates where a previously
    unevaluated integral now closes (each verified by PZQ).
  - Commit, bump, tag.
- [ ] **9. FriCAS stress** harness and corpus. Triage every WRONG/DISAGREE to root cause; set
  `BASELINE`.
- [ ] **9b. Integral generator and benchmark.**
  - Build `tools/integral_generator/generate.m`, the filters (N1/N2/N3), the FriCAS
    driver, the generator unit suites, the frozen corpus (~300 items: Trager, Shape A log,
    Shape A exp; genus 0 and ≥ 1), `benchmarks/256-generated-algebraic-integrals/`
    (`.m`/`.input`/`.py`/README), the nightly ctest with `BASELINE`, and the `make
    gen-integrals` target.
  - Triage every pending or WRONG item to its root cause.
  - Extend the generator to Shape B after Phase 11.
- [ ] **10. Shape B substrate.** `RDEAlgebraic.m`: Risch DE `z' + f z = Σ c_i g_i` over
  `Q(x)(y)`, `y^n = g(x)`, following [B] Part 3:
  - order bounds at normal places (Prop 3.3), infinite places (Prop 3.5, the log analogue
    for `d/dx`), the integral equation `(R₂)`, Prop 3.7;
  - the linear system from explicit binomial Puiseux expansions at x=∞ and at branch places;
  - Lemma 3.8 constant reduction.
  - Uses Trager's `af*` arithmetic and integral basis. Suite `test_bronstein_rde_algebraic.c`.
- [ ] **11. Shape B driver.** `BaseAlgebraic.m`: transcendental Risch over `k₀ = Q(x,y)`
  (Hermite in θ over k₀, the Rothstein–Trager resultant with k₀-norms, polynomial/special
  parts via `RDEAlgebraic`, base integration via Trager `TragerIntegrateCurve`). Examples:
  `Log[x] Sqrt[x^2+1]`, `E^x/Sqrt[x^2+1]` (certificate), `x E^x/Sqrt[1+x^2]` variants.
  Suites, stress corpus extension, docs, bump.
- [ ] **12. (Optional)** a TrigToExp front end with real form; the
  linear-combination-of-logs exp-residue case ([B] pp.155–157).
- [ ] **13. Close-out.** Review section in `tasks/todo.md`; lessons; valgrind; book index if
  the calculus chapter is touched (`cd book && make pdf`).

## Verification

- `make -j` clean; `make check-c99`, `make check-messages` (the new C glue routes through
  `mth_message`).
- `cd tests/build && cmake .. && make -j && ctest -R "bronstein|trager|risch" -LE nightly`
  all green, then the full `ctest -LE nightly -j4`.
- REPL spot checks:
  - `Integrate[(x+1)/((x Log[x]+x^2) Sqrt[x+Log[x]]), x]` → `−2/Sqrt[x+Log[x]]` (Automatic
    and Method → "Bronstein").
  - The p.147 and p.158 examples → the paper's answers (PZQ-verified).
  - `Integrate[Log[1+E^x]^(1/3)/(1+Log[1+E^x]), x]` → unevaluated plus `Integrate::nonelem`.
  - `Integrate[1/(x Sqrt[1-x^3]), x]` (Automatic) now closes via the Trager stage.
- `ctest -R "bronstein_fricas_stress|trager_fricas_stress|generated_integrals_corpus"`:
  zero WRONG/DISAGREE; report the agreement table and timings.
- `make check-bronstein-coverage`: every routine has tests.
- `python3 benchmarks/run_all.py --only 256`: the report shows per-tag closes for
  Mathilda/FriCAS/sympy.
- Stage-profile table before and after Phase 8: total cascade time on the existing
  integration suites must not regress beyond the agreed budget.
- `valgrind --leak-check=full ./Mathilda -file` on a short script of the examples: no growth
  beyond the startup baseline.

## Progress log

(append one dated entry per session: what was done, test status, the next step, open questions)

### 2026-11-12 — Phase 0 + 1 (session 1)

Done:
- **(a) `fieldFn`** in `trager/Field.m`. Spikes over Q(x) pass: the genus-0
  generator is a k-multiple of u1/u2, genus-1 `P − Q` is non-principal, and
  Q(√2)(x) works. Trager's generic modules needed only:
  - `nativeQ[F]` at every FLINT call site;
  - `afPolynomialQ` via `fieldScalarQ`;
  - `idealInverse`/`idealReduce` arms for `Fn` (`commonDenFactorFn`).
- **(b) Context bridge.** `Bronstein.m` appends `"Trager`Private`"` to
  `$ContextPath` inside `Bronstein`Private``, and short Trager names resolve
  there. Only `br*` names are created in `Bronstein`Private``. Beware: a test
  script that calls a Trager name without that context creates a `Global`` symbol.
- **(c) `Risch`RischDECertified` + `Risch`RischDECertifiedReason`.** This is a
  certification scope (`rde_cert_begin`/`end`/`mark` in
  `integrate_risch_rde.h`) with a mark on every non-authoritative NULL. Findings:
  1. A **real `Risch`RischDE` bug**: the polynomial loops did not canonicalise
     `c`, giving a false "no solution" (fixed with `rde_canon`).
  2. `rde_tower` has **no weak normalisation**, so a certification-only check
     was added.
  3. The tower degree bounds are exact only outside the cancellation cases;
     those are marked.

  `RDE.m` (`brRischDE`) wraps it all as {"Solution",y} / {"NoSolution"} /
  {"Unknown",why}.
- **(d)** `Risch`Derivation`/`SplitFactor`/`SplitSquarefreeFactor`/`CanonicalRepresentation`/`NormalQ`/`SpecialQ`
  work on θ as a plain symbol, including a 2-level tower `{x->1, t->t, th->t/(1+t)}`.
- Suites: `bronstein_field_fn_tests` (100 assertions) and `bronstein_rde_tests` (49).
  `tests/bronstein_test_utils.h` has `BR_EQ/BR_TRUE/BR_THROWS/BR_ANTIDERIV`;
  `BR_KEEP_GOING=1` reports all mismatches in one run.

Notes for the next session:
- Printed strings carry quotes in `BR_EQ` expectations (`"\"Fn\""`); builtin
  docstrings are checked with `symtab_get_docstring`, not `::usage`.
- Certified "NoSolution" is still conservative. Many tower equations come back
  `$Failed` ("RdeBoundDegreePrim: cancellation case", "tower: f not weakly
  normalized"). If Phase 3/6 needs more certificates:
  - weak-normalise `f` in `brRischDE` (in `.m`) before the call;
  - implement the limited-integration bound (`LimitedIntegrate.m`).
- Next: Phase 2 (`Tower.m`, `Poly.m`, `Curve.m`).

### 2026-11-12 — Phase 2 (session 1)

Done:
- **`Tower.m`** (`brTower`): kernels → symbols `brT<i>` innermost-first.
  - Exp relations come from `Risch`ExpReducible`. Integer relations are
    rewritten (an identity); a fractional one re-bases the exp monomial
    (`brRebaseExp`, `E^b = (E^(b/q))^q`).
  - A reducible log declines, because it is not an identity on every branch.
  - The radical is found by `brDetectRadical`. The shape is one of Rational /
    Transcendental / Algebraic / A / B; for A, θ is moved to the top, and a
    monomial above θ that depends on it declines.
  - `brFlattenLinear` handles a Möbius radicand in x.
- **`Poly.m`**: content in k, `brSplit`/`brSplitSquarefree`/`brNormalQ`/`brSpecialQ`
  (θ-free counts as special: `Risch`SpecialQ[1]` is False),
  `brSpecialPlaces`.
- **`Curve.m`**: `brCanonicalForm` over k(θ) (Trager's power-basis inverse),
  `brComponentCurve` (Trager's `componentCurve` with θ as the curve variable,
  plus the two constant-field-extension declines), `brCurve`
  (`s = (1/m) DG/G`, content, Trager basis), `brReducedCurveQ`.
- Suites: `bronstein_tower_tests` (107 assertions), `bronstein_curve_tests` (27).
- Pitfall hit again: `Return` inside `With` inside `Module` only leaves the
  `With`. Use `Which`/`If` values instead.
- Deferred: the Kummer split for a compositum declines for now; plan it in
  Phase 12.

### 2026-11-12 — Phase 3 (session 1) and a design revision

**Design revision: no `Special.m`.** Re-reading [B] pp.140–146 changed two
things.
- The special-place residue of the log case (Def 2.10) needs only the θ⁻¹
  coefficient of the Puiseux expansion of the *normal* f. Nothing has to be
  reduced at ∞ first.
- For an n-th root the whole algebraic part is the coordinate system `S_log2`
  (Prop 2.15), and its solution is a POLYNOMIAL in θ (Prop 2.14).

So the [T] App A special-place reductions are subsumed by `AlgebraicPart.m`:
- **(1)** `Risch`RischDECertified` over the whole tower. It decides the exp case
  (p.158 certified).
- **(2)** When (1) is undecided and θ is a log, `brPolyRDELog`. The C log case is
  always "RdeBoundDegreePrim: cancellation case", because `s = DG/(mG)` has
  `deg b = deg a − 1`. `brPolyRDELog` is the top-down coefficient solve, with:
  - N = ⌊max(deg c + 1, 1 − deg G/m)⌋;
  - one Risch DE in k per coefficient (β = lcG'/(m lcG));
  - when lcG = const·u^m, a homogeneous solution h = 1/u, so each coefficient is
    a parametric integration in k with a free constant, and the lower
    equations become a linear system for the constants.

**`LimitedIntegrate.m`** (`brIntegrateInK`, `brSplitLogs`, `brLogBasisReduce`,
`brParamIntegrate`, `brLimitedIntegrate`):
- each piece is integrated by `Integrate`RischTranscendental`;
- new logs are reduced to a k-independent basis by `Risch`LogReducible`;
- the constants then solve a linear system.

A "no" is certified only from `ElementaryIntegralQ` False, LogReducible False or
an inconsistent system. ArcTan/RootSum answers give Unknown.

`Tower.m` gained `brFromSurface` (C answers back to symbols), `brLowerTower` (the
record of k) and a `"monos"` key.

Results:
- p.134 closes by Hermite alone;
- p.158's algebraic part is `3x`, and a temporary driver closes both end to end
  (PZQ-verified);
- p.147's special-place quantity gives λ = 1, v = x − x².

Suite `bronstein_hermite_special_tests` (76 assertions).

**Harness trap.** A `br_run` string that both creates a tower and mentions
`brT1` is parsed before `Bronstein`Private`brT1` exists, so the parser makes
`Global`brT1`, which shadows it for the rest of the run. Create the towers in
their own `br_run` (the suite asserts `Context[brT1]`).

Next: Phase 4.
- `Residues.m`: R(z) and the Prop 2.9 criterion.
- Puiseux at θ = ∞ for the log case. For a component Y (i ≥ 1), m ∤ deg G gives
  τ = 0; m | deg G with a constant lc(G) gives limited integration. A
  non-constant non-m-th-power lc declines for now.
- Then read [B] pp.148–158 for the exp special places (Def 2.21) before writing
  them.

### 2026-11-12 — Phase 4 (session 1)

- `Puiseux.m`: Laurent coefficients at θ = 0 and ∞ (`brSeriesCoefficients`,
  power-series division, with no `Series` call), binomial series (Miller
  recurrence), and place data at 0 and ∞ (count, ramification, lc/G0).
- **`Special.m` is back, for the EXPONENTIAL case only.** Def 2.21 needs a
  *reduced* f.
  - In one coordinate each step is a Risch DE in k for the extreme Laurent
    coefficient: `c1' + (s(0) − μη')c1 = a` at 0, and
    `c1' + (s(∞) + dη')c1 = a` at ∞.
  - It is valid even at ramified places, with fractional exponents: the
    homogeneous solution is not in k̄, and Lemma 2.16 applies.
  - "NoSolution" certifies (Prop 2.18).
  - p.158 reduces to 0, so its algebraic part `3x` is found here.
- `Residues.m`:
  - **R(z)** uses Trager's branch-free `Bp` and the derivative of the full
    normal denominator. Prop 2.9 gives the certificate, and p.139 is certified.
  - **Log, θ = ∞:** τ = 0 when m ∤ deg G. Otherwise limited integration of the
    θ⁻¹ coefficient, giving τ_ζ = −ζκ^{1/m}λ.
  - **Exp, special places:** values at unramified places, and τ_p = εr(c_p − c_p̄)
    with p̄ the ζ = 1 place above 0. This agrees with the true divisor up to
    div(θ), which is principal. Coordinate 0 absorbs the resulting θ factor of
    the generator.
  - The algebraic constants γ are supported only as κ^{1/m}·u, u ∈ k; otherwise
    the stage declines.
- **p.147 ground truth corrected.** The surface integrand in the plan was
  mis-transcribed; its residue at the y = ±x places was non-constant, which
  looked like a "false certificate". The paper's θ-form is right. The
  y-component has residues ±1, i.e. `log((x+y)/(x−y))`, and the rest of
  `2 log(x+y)` is coordinate 0.
- Pitfall: `FirstPosition[list, pat]` has default level spec `{0, Infinity}` and
  can match the whole list. Use an explicit scan.
- Suite `bronstein_residues_tests` (82 assertions).

Next: Phase 5.
- `Divisors.m`: finite places via `residueIdeal` over `fieldFn`.
- Special places by a shift (log ∞: θ = 1/u; exp 0/∞). Initially support only
  divisors supported at finite normal places, and decline when a special
  residue is nonzero; record it.
- `Specialize.m`: torsion by specialisation.

### 2026-11-12 — Phases 5–7 (session 1), v0.6518

- **`Divisors.m`**:
  - residue field = Trager's `residueField` and `qBasis` (over Q), lifted to
    `fieldFn[{x, lower}, {t, mu, theta}]`;
  - divisor ideals = Trager's `divisorIdeal` with the tower-derivative `Ep`;
  - principality = Trager's `principalGenerator`.
- **Special-place divisors (genus 0).** Shift frame θ = a + 1/u, Ŷ = Y u^k:
  - unramified special places are point ideals (u − u0, Ŷ − w);
  - ramified ones use Trager's `placeIdealAt`;
  - finite residue places are rewritten in the frame (`residueIdeal` with
    `(N̂, Ê)`);
  - the field also contains the roots of `z^m − κ`.

  Genus ≥ 1 with nonzero special residues declines; it would need the special
  places carried through the specialised torsion computation.
- **`Specialize.m`**: torsion order by a good specialisation. The goodness checks
  are deg G, the squarefree profile, Bp degree/squarefree/coprime, denominators,
  and the residue-place degrees. Trager's `torsionOrderSpec` gets a context with
  `N̄ ≡ N·Bp'/(DE) mod Bp`, so its d/dθ frames carry the same places. Then
  generic principality of N·δ.
  - Order 3 over Q(x) is found (`y² = θ³ + x²`).
  - The lifted Nagell–Lutz point is certified.
- **`LogPart.m`**: L = Σ(β/N) log u, projected to weight 1 as L − (1/m)Σ log Norm(u);
  m ≥ 3 with nonzero weight-≥2 coordinates declines.
- **`Pipeline.m`**: `BronsteinIntegrate`, with these hand-offs:
  - no radical → transcendental Risch;
  - no kernel → Trager;
  - Shape B with a Möbius radicand → flattened;
  - otherwise Shape B declines.

  Components are assembled and the final check is `antiderivativeQ`.
- **C**: `METHOD_BRONSTEIN`, `Integrate`Bronstein`, `try_bronstein`.
- **Exp special values with γ = G0^{1/m} ∉ C·k** (`brRadicalValueCheck`). By the
  Galois weight of γ, the η-coefficient is 0, and the condition becomes "∫γw is
  pure logs" (Def 2.21) or "∫γ a1 is purely algebraic" (Def 2.10), decided by
  Trager over Q(x). For a non-empty lower tower it declines.
- **FALSE CERTIFICATE found and fixed.** An unqualified `TragerIntegrateCurve`
  stayed unevaluated because `Trager`` was not on the package path; see
  `tasks/lessons.md` and memory `feedback-certificate-from-positive-evidence`.
  Two other certificate sites (`brExpValueIntegral`, `brResultantCriterion`) now
  validate shapes too.
- Suites: `bronstein_divisors_torsion_tests` (51), `bronstein_integrate_tests` (68,
  including 16 random derivatives with 0 wrong verdicts).

Next: Phase 8 (cascade).
- The patch is prepared in this session: gates `Integrate`CascadeGateQ`, budgets
  `$TragerCascadeBudget`/`$BronsteinCascadeBudget` (default 10 s), and a
  `cascade_nonelem` local flag that skips PMT/Goursat. The suite is
  `tests/test_integrate_cascade_trager_bronstein.c`.
- Then measure the stage profile and run the full ctest.
- Presentation polish: `Log[-E^x]` → x; `Log[x^6]/2` → 3 Log[x].

### 2026-11-12 — Phase 8 (session 1), v0.6519

- The cascade stages "15b Trager" and "15c Bronstein" are in place (see the
  `integrate.c` comment at `try_cascade_decider`):
  - `Integrate`CascadeGateQ` exposes the gates;
  - budgets `$TragerCascadeBudget`/`$BronsteinCascadeBudget` (10 s);
  - a local `cascade_nonelem` flag skips PMT/Goursat;
  - the cascade's held `nonelem` is dropped when PMS answers (existing design).
- `ctest -LE nightly` 647/647.
- **Open:** the before/after stage profile was killed by system memory pressure.
  The `N::meprec` flood came from a parametric DDI integrand in another stage.
  Rerun it on a reduced corpus (no symbolic parameters, one process per
  integrand), then set the budgets from the data.

## Stress campaign and open-gap plan (2026-11-12, after v0.6520)

### What the stress testing found

Tooling: `tools/bronstein_stress/`.
- `gen.m` is a seeded generator of KNOWN-elementary `f = D[F]` across six families
  (`alg`, `logxy`, `asinh`, `genus1`, `algres`, `ratb`). `FAMILY=random` gives
  random `R(x,θ)·Y^k` of unknown verdict.
- `run.m` runs a batch; `stress.py` is the memory-safe driver, one process per
  batch, with timeouts.
- `fricas_check.py` asks FriCAS 1.3.6 about every certificate and decline.

| Run | Cases | Result |
|---|---|---|
| Method→Bronstein, seeds 1–3 | 45 | 41 OK, 3 DECLINE, 1 TIMEOUT |
| Method→Bronstein, seeds 10–39 | 600 | 571 OK, **10 WRONG**, 12 TIMEOUT, 7 DECLINE |
| the 29 non-OK of that run, after v0.6520 fixes | 29 | 15 OK, 2 DECLINE (genus-1 special residue), 12 TIMEOUT |
| random family + FriCAS, seeds 100–107 | 160 | 160 certificates; FriCAS agrees on all 160 (unevaluated or "implementation incomplete") |
| plain Integrate (cascade), seeds 200–207, 60 s | 120 | 109 OK, 0 WRONG, **11 TIMEOUT** (cases Method→Bronstein closes in 0.02–4 s) |

Bugs fixed (v0.6520):
- a wrong answer from an unverified presentation (`Log[2 + 2x]` → `Log[1 + x]`
  inside a radical);
- missing geometry for exponential special places that have value 0;
- ArcTan in a special-place value integral.

Earlier, in v0.6518, a false certificate came from an unresolved package name.

**FriCAS observations.** FriCAS 1.3.6 answers "implementation incomplete (constant
residues)" on many mixed integrands, e.g. `Sqrt[Log[x]]` and `1/Sqrt[Log[x]]`. So it
can neither close nor certify them. These are N3 material for the benchmark.

### Open gaps, by priority

**P0 — soundness hardening** (before anything else):
1. **Verified-form invariant in Trager too.** `TragerIntegrate` also returns
   `presentable[...]`, whose fast path is unverified. Trager's integrands have no
   log inside a radical, so the hazard is latent, but it must be closed. Re-verify
   or fall back, as `brPresentable` does. Add a regression test.
2. **Certificate-site audit.** There are 12 `throwNonElementary` sites. Each must be
   reached only from a validated positive answer (a shape check). Each site gets a
   unit test that feeds a malformed or unevaluated sub-answer and asserts
   Decline/Internal, never NonElementary.
3. **Generator coverage:**
   - constant factors in every kernel argument (`Log[c x + c]`, `E^(c x + d)`);
   - two-level towers (`Log[1 + E^x]`, `E^(Log[x]^2)`-type);
   - a non-constant `lc(G)` and content `d ∈ k` (`x θ + 1`, `x^2 θ^2 + 1`);
   - exponential cube roots with special residues;
   - a **perturbed** family (an elementary f plus a non-elementary piece that keeps
     the residues constant). This exercises the deeper certificates: algebraic
     part, special reduction, torsion. Each such certificate is checked against
     FriCAS where FriCAS decides.
4. **Nightly ctests**:
   - `bronstein_stress`: fixed seeds, assert zero WRONG/FALSE_CERT/CRASH, with a
     ratcheting baseline of closes;
   - `bronstein_fricas_stress`: zero DISAGREE.
5. **Coverage gate** `tools/check_bronstein_test_coverage.py` with
   `make check-bronstein-coverage`. This is required by the testing policy and not
   built yet: every `::usage` name and new builtin must be called from a suite.

**P1 — user-visible quality:**
6. **Performance of m ≥ 3 answers.** This covers all 12 stress timeouts
   (`E^(x²)` cube and fourth roots). The components take about 3 s; the time goes
   to assembling, verifying and presenting huge ζ-laden logs. Plan:
   - (a) shrink generators: remove k-content and normalise coordinates;
   - (b) emit conjugate log families as one `RootSum` (Trager's per-orbit form)
     instead of m explicit ζ-logs;
   - (c) verify per component (the identity in k(θ)(Y) is already exact),
     leaving a lighter final PossibleZeroQ;
   - (d) a presentation budget.

   Acceptance: no stress timeout at 30 s.
7. **Presentation.** All re-verified:
   - `Log[-E^x]` → `x`;
   - `Log[x^6]/2` → `3 Log[x]`;
   - for m = 2, the pair `c Log[u] − (c/2) Log[N(u)]` → `(c/2) Log[u/σu]`;
   - drop constant logs;
   - real ArcTan forms.
8a. **Cascade grind before the decision procedures.** All 11 cascade timeouts are
    integrands Bronstein closes in under 5 s. The stage profile names
    **DerivativeDivides** (stage 12): 89.7 s of a 90 s run, one aborted call,
    before stages 15b/15c are reached. Fix options, cheapest first:
    - (a) gate DerivativeDivides' Eliminate/Solve search off when
      `Integrate\`CascadeGateQ["Bronstein"|"Trager", f, x]` holds. This is the
      pattern already used to keep PMT off pseudo-elliptic shapes;
    - (b) move the gated Trager/Bronstein stages ahead of DerivativeDivides for
      gate-matching integrands only.

    Either must be measured: a reordering can expose latent grind elsewhere, and
    the closed-form shape of integrands that DerivativeDivides closes today must be
    kept. Acceptance: the cascade stress run (seeds 200–207) with no timeout at
    60 s, and `ctest -LE nightly` green.
8. **Cascade budgets.** Rerun the stage profile on a reduced corpus (no symbolic
   parameters, one process per integrand). Set `$TragerCascadeBudget` /
   `$BronsteinCascadeBudget` from the data.

**P2 — coverage inside Shape A** (each a decline today, never a wrong answer):
9. **Special-place residues on genus ≥ 1.** Two stress declines. Carry the special
   places through the specialised torsion computation: Trager's
   `torsionOrderWith[ctx, mk]` accepts an ideal *builder*. Pass one that also
   builds the point and ramified-place ideals on the shifted frame, mod p. Then
   test generic principality on the frame (code exists for genus 0).
10. **Log-like atoms with algebraic coefficients** (ArcTan, RootSum, complex logs)
    in `brLimitedIntegrate`/`brParamIntegrate`, i.e. the log-case Def 2.10 and the
    coefficient integrations. Generalise `brLogBasisReduce` to logs over `k(t)`.
11. **Undecided Risch DE sub-cases in the exponential case** (`RdeSpecialDenomExp`
    / `RdeBoundDegreeExp` parametric-log-derivative marks). Write `brPolyRDEExp`,
    a top-down Laurent-polynomial solve analogous to `brPolyRDELog`. Its leading
    equations never cancel: the homogeneous solution needs θ^(j+g/m) ∈ k̄.
12. **m ≥ 3 log parts with weight-≥2 coordinates.** Not seen in the stress runs so
    far. The general projection `(1/m) Σ ζ^{-k} σ^k L` needs ζ adjoined to the
    residue field.
13. **Reducible logarithms** (`Log[x^2]` beside `Log[x]`). Today they decline.
    Rewrite derivative-equivalently with an explicit constant
    (`Log[x^2] = 2 θ + c`, with c recovered numerically and checked exactly) and let
    the final verification decide.
14. **Compositum of radicals.** Use the Kummer split, as Trager's pipeline does.
15. **Radical γ at special places over a non-empty lower tower.** This needs
    integration of an algebraic function over k, which is Shape B. It follows #16.

**P3 — scope extensions** (the original plan's later phases):
16. **Shape B** (Phases 10–11): the radical over Q(x) below the tower. This also
    unlocks #15.
17. **Integral generator and benchmark 256** (Phase 9b). **Done 2026-11-12, as the
    integration shootout** (`benchmarks/256-integration-shootout/`,
    `tools/integration_shootout/`): 1151 cases in five classes, every applicable
    Mathilda method beside FriCAS, Maxima, sympy and Wolfram, each answer checked by
    its own system. Nightly ctest `integration_shootout_stress` with a ratchet
    baseline. Shape B items join it after Phase 11.
18. **Phase 12**: a TrigToExp front end with real forms; the [B] pp.155–157 case
    where η' is a linear combination of logarithmic derivatives.

**Order of work:** 1–5 → 8a → 8 → 6–7 → 9 → 10 → 11 → 13–14 → 17 → 16 (with 15) → 12 → 18.

### 2026-11-12 — integration shootout (v0.6526)

First full run: Mathilda `Automatic` closed 870 of 1151 (Wolfram 908, FriCAS 666,
Maxima 577, sympy 443), with 5 wrong answers and no false certificates. Fixed:
- the 5 wrong answers: the cubic Legendre form in ParallelMixedSpecial had a
  piecewise sign (`CubicData` now orders the roots for a constant sign);
- the one AUTO-GAP (`lccontent-2-3`): Bronstein's optional presentation overran the
  cascade budget and lost the verified answer (`brPresentationTime`);
- a garbage answer (`$Failed + tragerFail[..]`): a `$RecursionLimit` overflow in
  ParallelMixedSpecial left the evaluator degraded for the rest of the cascade
  (contained at each stage now).

Final run: 878 closed, 0 WRONG / FALSE_CERT / DISAGREE / AUTO-GAP. Targets it
leaves (see `REPORT.md`):
- **A, algebraic, 57% vs Wolfram 83%:** torsion-heavy (0/18), genus hyper/super,
  Chebyshev and binomial families. Wolfram's extra closes there are mostly
  hypergeometric, which no Mathilda method produces for algebraic integrands.
- **B, radical at the bottom, 48% vs 72%:** Bronstein declines this shape, so it
  needs Shape B (#16).
- **M random/perturbed:** 0/33 `mrandom` for every system (most are presumably
  non-elementary, unproven).
- **Why ParallelMixedSpecial overflows `$RecursionLimit`** on
  `9/(x (-3 + 2x + Log[1+x^2])^(3/2) + x^3 (..)^(3/2))`: contained, not yet
  explained.

### 2026-11-12 — cascade grind fixed (v0.6521)

Plan item 8a is done. The fixes:
- DerivativeDivides has a 0.2 s budget per candidate in the cascade (user
  direction), with its sub-integral bounded.
- Trager and Bronstein moved ahead of Linearity, DerivativeDivides and
  RischTranscendental: stages "08b"/"08c".
- `g_integrate_depth`, `g_integrate_no_special` and `dd_depth` are registered
  unwind counters. A nested-timeout depth leak was silencing `Integrate::nonelem`
  for the whole session.

Cascade stress (seeds 200–207): 120/120 OK.

Next: the exact per-component verification (the Method→Bronstein m ≥ 3 timeouts
are 98 s spent in the final PossibleZeroQ), then the P0 items.

### 2026-11-12 — exact verification, P0 items 1–2 (v0.6522)

- `brComponentIdentityQ`: exact per-component verification on Y^m = G, replacing
  the global PZQ (97 s → 3.4 s).
- A core `With`/`Module`/`Block` initializer-Throw bug was found through the
  Bronstein certificates and fixed in `src/modular.c`.
- P0 #1 (Trager presentation re-verified) and P0 #2 (certificate-site audit:
  11/12 already read explicit answers, and the constant-system "no" is now
  shape-checked) are done.
- Stress (600): 592 OK, 2 DECLINE (genus-1 special residues, gap #9), 6 SLOWCHECK
  (5–6k-leaf m ≥ 3 answers, gap #6/#7). Cascade 120/120.


### 2026-11-12 — generator content removal (v0.6523, gap #6/#7 part 1)

- `brNormalizeGenerator` (LogPart.m) divides each log generator u by the
  k(θ)-content of its coordinates (lcm of the denominators over the gcd of the
  numerators). The weight-1 term Log[u] − (1/m) Log[Norm u] does not change when u
  is scaled by g ∈ k(θ), because Log[g] − (1/m) Log[g^m] is a constant. So this only
  shrinks the answer, and the exact per-component verification still runs on the
  normalised form.
- Worst m = 4 case: log arguments of 342/375 leaves become 60/75; the answer goes
  from ~6k to 1.8k leaves.
- Stress (600): 598 OK, 2 DECLINE (gap #9), **0 SLOWCHECK**.
- Still open under #6/#7: RootSum per-orbit presentation; `Log[-E^x]` → x.

### 2026-11-12 — P0 #3–#5: generator coverage, nightly ctests, coverage gate (v0.6524)

- **#5 done**: `tools/check_bronstein_test_coverage.py` / `make check-bronstein-coverage`
  / push-tier ctest `bronstein_test_coverage`. It covers 140 routines, with 1 exempt
  (an inert symbol). It found `brRadicalValueCheck` returning from inside `With`
  (a rational γ declined instead of passing); fixed and tested.
- **#4 done**: nightly ctests `bronstein_stress` (600 cases, `--max-open 8`),
  `bronstein_cascade_stress` (Automatic, `--max-open 0`) and `bronstein_fricas_stress`
  (random plus perturbed, zero DISAGREE).
- **#3 done**: four new families (constk, tower2, lccontent, expcube), plus a
  perturbed near-miss half in the random family. On their first run they found:
  - **a false certificate**: the `Risch`RischDECertified` canonical-arithmetic bug,
    with a non-constant θ-leading coefficient (fixed in C);
  - `Integrate`RischTranscendental` declining a single-kernel proper part plus a
    new k-log (fixed: single-kernel field recursion);
  - norm terms left as Log[const·p] and never merged (fixed: `brSplitLogs`,
    `brCollectLogTerms`; expcube went from 61 s and an unverifiable PZQ to 1.3 s).
- Stress (600, 10 families): 592 OK, 0 WRONG/FALSE_CERT, 8 DECLINE:
  - 2 genus-1 special residues (#9);
  - 3 `TorsionBeyondWeilBound` (lccontent, new);
  - 2 special values over a radical extension of k (#15);
  - 1 non-constant content over a radical constant field.
- Still open under #6/#7: presentation of conjugate pairs (an ArcTan pair summing
  to a constant survives the real form). Next: #8 cascade budgets, then the
  lccontent declines.
