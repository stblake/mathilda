# Review: Transcendental Risch Integration & a Roadmap for Cherry's Special-Function Extension

*Date: 2026-07-11 · Scope: `src/calculus/integrate_risch_transcendental.c` (4092 LoC),
`src/calculus/intrischnorman.c` (pmint, ~4300 LoC), the `Integrate` dispatch cascade,
and the special-function surface in `src/special_functions/`.*

---

## 1. Executive summary

Mathilda has **two** transcendental integrators, layered in the `Integrate` Automatic
cascade in this order (`src/calculus/integrate.c:783-812`):

| Stage | Symbol | Engine | Nature |
|------|--------|--------|--------|
| 2  | `Integrate`RischNorman` | `intrischnorman.c` | Parallel-Risch / Norman **heuristic** (Bronstein's *pmint*) |
| 2b | `Integrate`RischTranscendental` | `integrate_risch_transcendental.c` | Recursive Risch **decision procedure** (Bronstein/Roach lineage) |

`RischTranscendental` is the mathematically principled one: every branch fires behind an
exact structural certificate (`SolveAlways`, exact resultants), and every emitted closed
form is additionally **diff-back verified** by `rt_verify_antideriv`
(`integrate_risch_transcendental.c:2075-2089`, `Simplify[D[result,x] - f] === 0`). It is
correct by construction and cannot ship a spurious answer.

**The current special-function capability is real but shallow.** It consists of four
*ad-hoc pattern recognizers* — not Cherry's decision procedure — collected in
`rt_special_case` (`:510-517`) and tried **only at the top level, only on the whole
integrand**:

- `rt_try_erf` (`:279-322`) — `K·E^(ax²+bx+c)` → `Erf`/`Erfi`
- `rt_try_ei`  (`:337-382`) — `(M·E^(ax+b))/(cx+d)` → `ExpIntegralEi`
- `rt_try_li`  (`:384-420`) — `c·w^(p-1)·w'/Log[w]` → `LogIntegral`
- `rt_try_dilog` (`:478-508`) — `K·Log[1+px]/x` → `-K·PolyLog[2,-px]`

Each is diff-back verified against a genuine derivative rule that already exists in
`deriv.c:339-391` (`Erf`, `Erfi`, `Erfc`, `ExpIntegralEi`, `LogIntegral`) and
`deriv.c:1172` (`PolyLog`). So they are *correct where they fire* — but they fire on a
vanishingly small set of shapes.

**Incorporating Cherry means replacing these four pattern-matchers with a genuine
decision procedure** that (a) decides *whether* `Erf`/`Li`/`Ei` terms are needed, (b)
computes their arguments, and (c) does so recursively inside the tower, reusing the
Laurent-ansatz + RDE machinery that already exists for the elementary case.

---

## 2. Architecture as built

### 2.1 `RischTranscendental` — the recursive decision procedure

Top-level flow (`rt_integrate`, `:4042-4050`): rational base → transcendental → special.
The transcendental dispatcher (`rt_transcendental_case`, `:4000-4022`) tries ten cases in
order: log-poly, exp-poly, fractional (Rothstein–Trager), Hermite, hyperexponential,
exp-sum, log-tower, exp-tower, recursive-tower, trig-frontend.

**Bronstein algorithms present** (verified by structural map):

| Algorithm | Location | Status |
|-----------|----------|--------|
| Hermite reduction (log & exp kernels) | `rt_hermite_try :1479-1678` | ✅ |
| Rothstein–Trager log part (rational residues) | `rt_frac_try :1166-1357` | ✅ |
| Lazard–Rioboo–Trager (algebraic residues, pure resultant) | `rt_frac_lrt :1359-1476`; tower: `rt_field_lrt_logpart :3007-3068` | ✅ |
| Risch differential equation (RDE) + degree bound + resonance | `rt_rde_var_bound :726-743`, `rt_solve_rde_rational :793-903`, `rt_solve_rde :904-977`, `rt_resonance_int :754-765` | ✅ |
| Coupled hyperexponential (unified Laurent + log ansatz via `SolveAlways`) | `rt_hyperexp_case :1719-1948`; tower: `rt_field_hyperexp_coupled :3243-3425` | ✅ |
| Tower construction + commensurate-exponent reduction | `rt_tower_build :2778-2956` (`RT_MAXK=5`) | ✅ |
| Field-level recursion to `C(x)` base | `rt_field_integrate :3431-3494` → `Integrate`BronsteinRational` | ✅ |
| Trig via `TrigToExp`/`ExpToTrig` | `rt_trig_frontend :3967-3992` | ✅ (real-form gap, below) |

This is a serious, largely faithful recursive Risch engine. The degree bounds are
**derived, not capped** (`rt_rde_var_bound` implements Bronstein's `RdeBoundDegree` with
resonance widening) — consistent with the project's "no arbitrary caps in a decision
procedure" rule.

**Structurally deferred** (author's own comment, `:2726-2728`):
> "The genuine algorithm's proper-rational part (tower Hermite + Rothstein–Trager over
> K_{n-1}) and its general field Risch-DE are deferred to a later increment; a nonzero
> proper part or a non-base RDE declines cleanly here."

So at tower depth ≥ 2 the engine handles the polynomial/Laurent part with
coefficient-recursion but **declines a proper-rational part over an inner field**. This is
the single largest elementary-completeness gap; it also constrains where Cherry terms can
appear (see §5).

**`ParametricLogarithmicDerivative` / `LogarithmicDerivativeOfARadical`**: not present as
named routines. Their role (deciding when an exponential's integrating factor closes) is
partly subsumed by the `SolveAlways`-certified ansätze, but the general primitive is
missing — relevant to full RDE completeness and to Cherry's exp-case argument search.

### 2.2 `RischNorman` (pmint) — the heuristic

A near-line-for-line C port of Bronstein's *pmint*: `convert_to_tan` → collect
indeterminates → vector field → `splitFactor`/deflation → monomial ansatz →
undetermined-coefficient linear solve → log-candidate sum, with a `K=I` retry
(`intrischnorman.c` pipeline at `:3549-3954`).

Two findings matter for this review:

1. **Its `getSpecial` table is gutted.** In full pmint, `getSpecial` seeds the ansatz with
   special-function Darboux candidates (Erf, Ei, li, dilog). The C port's
   `get_special_all` (`:3213-3249`) contains **only `Tan`, `Tanh`, and `LambertW`**. So the
   heuristic path contributes *nothing* toward Erf/Ei/Li — all special-function capability
   today lives in `RischTranscendental`'s four recognizers.

2. **It is heuristic, with hard caps** (`PMINT_MAX_MONOMIALS=5000`, `MAX_INDETS=32`,
   `BUDGET_SEC=4.0`, solve-row/leaf caps) and — despite a header comment — **no active
   post-hoc differentiation check** on the final result. This is acceptable for a heuristic
   that runs *before* the verified engine, but it means pmint must never be the thing we
   extend for *correctness-critical* special-function work. Cherry belongs in
   `RischTranscendental`.

### 2.3 Special-function surface today

| Function | Builtin | Derivative rule | Cherry-relevant |
|----------|:------:|:---------------:|:---------------:|
| `Erf`, `Erfc`, `Erfi` | ✅ | ✅ `deriv.c:339-363` | core (error-fn paper) |
| `ExpIntegralEi` | ✅ | ✅ `deriv.c:384-389` | core (log-integral paper) |
| `LogIntegral` | ✅ | ✅ `deriv.c:390-392` | core (log-integral paper) |
| `PolyLog` (Li₂ = dilog) | ✅ | ✅ `deriv.c:1172` | extension (Baddoura) |
| `Gamma`, `LogGamma`, `PolyGamma` | ✅ | ✅ | arithmetic only |
| `SinIntegral` (Si), `CosIntegral` (Ci) | ❌ | ❌ | needed (imaginary-Ei) |
| `SinhIntegral` (Shi), `CoshIntegral` (Chi) | ❌ | ❌ | needed (real-Ei) |
| `FresnelS`, `FresnelC` | ❌ | ❌ | needed (imaginary-erf) |
| `ExpIntegralE` (Eₙ) | ❌ | ❌ | optional |

The derivative infrastructure for the *core* Cherry functions already exists — this is the
good news. The trigonometric/hyperbolic integral siblings (`Si/Ci/Shi/Chi/Fresnel`), which
Cherry's algorithm produces naturally as the imaginary-argument images of `Ei`/`Erf`, are
absent.

---

## 3. Correctness assessment

**Strong.** The `RischTranscendental` design is sound and defensively verified:

- Every special-function emission passes through `rt_verify_antideriv` before returning,
  and the required derivative rules genuinely exist, so the gate has teeth (I confirmed
  `deriv.c:339-392`). A recognizer that mis-fires produces a non-zero `Simplify[D-f]` and
  is rejected.
- The elementary branches are certified by `SolveAlways`/exact resultants; the "residue
  must be constant" class of Risch bug (see `[[project_solvealways_nonconstant_residue]]`)
  is guarded structurally.
- Degree bounds are derived (`rt_rde_var_bound`), not magic constants — no completeness is
  silently thrown away.

**Caveat to keep honest:** pmint (`RischNorman`) runs *first* in the cascade, is a heuristic
with size caps, and has **no** post-hoc verification. Because it precedes the verified
engine, a wrong pmint answer would be returned to the user without a diff-back check. This
is a pre-existing risk, orthogonal to Cherry, but worth a follow-up: gate pmint's output
with the same `Simplify[D-f]===0` check `RischTranscendental` uses. Cheap insurance.

---

## 4. Cherry's theory (what we are actually adding)

Two papers, both titled *"Integration in Finite Terms with Special Functions"*:

- **The Error Function** — G. W. Cherry, *J. Symbolic Computation* **1** (1985) 283–302.
- **The Logarithmic Integral** — G. W. Cherry, *SIAM J. Comput.* **15** (1986) 1–21.
  (Foundation: Cherry's 1983 Delaware PhD thesis under B. F. Caviness.)

Liouville's theorem says an elementary `f` has an elementary integral **iff**
`∫f = v₀ + Σ cᵢ log(vᵢ)` with `v₀, vᵢ` in the field and `cᵢ` constant. Cherry extends the
**allowed form of the answer** by one family of new terms:

- **Error-function extension.** `∫f = v₀ + Σ cᵢ log(vᵢ) + Σ dⱼ·erf(uⱼ)`,
  where `erf(u)' = (2/√π)·u'·e^(−u²)`. An `erf(u)` term is admissible exactly when
  `e^(−u²)` is a *monomial of the tower* — i.e. `−u²` equals (a rational multiple of) one
  of the exponential arguments. So error-function terms are exactly the "the exponent is a
  perfect quadratic" case, generalized from a single Gaussian to any Laurent-in-`θ`
  coefficient.

- **Logarithmic-integral extension.** `∫f = v₀ + Σ cᵢ log(vᵢ) + Σ dⱼ·li(uⱼ)`,
  where `li(u)' = u'/log(u)` and `Ei(u) = li(e^u)`. A `li`/`Ei` term is admissible when the
  integrand has a `1/log`-type (log tower) or `e^{monomial}/monomial` (exp tower) component
  whose residue closes.

**The algorithmic content — and why it fits Mathilda cleanly:** in both cases Cherry
reduces "does a special-function term exist, and what is its argument?" to the *same two
primitives the elementary Risch case already uses*:

1. a **Laurent/partial-fraction decomposition in the top monomial `θ`** (already:
   `rt_hyperexp_case`, `rt_field_hyperexp_coupled`), and
2. a **Risch differential equation** per coefficient to absorb the residual elementary part
   (already: `rt_solve_rde*`), plus a small algebraic side-condition ("is `−kL` a perfect
   square in the field?" for erf; "does the residue close?" for li/Ei).

In other words, Cherry is not a new engine — it is **an augmentation of the existing
hyperexponential/RDE ansatz with extra candidate basis terms**, decided rather than
pattern-matched. This is the crucial insight for the roadmap: we already own the hard
machinery.

Beyond Cherry proper: **Baddoura** (*Integration in finite terms with elementary functions
and dilogarithms*, J. Symbolic Comput. 1994/2011) extends this to `Li₂`, which is the
principled home for the current `rt_try_dilog` shape.

---

## 5. Gap analysis — current recognizers vs. Cherry's decision procedure

| Capability | Now | Cherry |
|-----------|-----|--------|
| Gaussian `e^(quadratic)` | ✅ single term only | ✅ any `p(x)·e^(quadratic)`, `p` polynomial |
| Gaussian **× rational** (`e^(−x²)/(x−a)`, `x·e^(−x²)·log x`) | ❌ | ✅ via Laurent + RDE, mixes `erf` with elementary terms |
| `e^(linear)/(linear)` → `Ei` | ✅ | ✅ |
| `e^(poly)/(poly)`, higher-degree denominator → `Ei` | ❌ | ✅ |
| `1/log(x)`, `w^(p-1)w'/log w` → `li` | ✅ narrow | ✅ general `R(x)/log(w)` |
| `Si/Ci/Shi/Chi` (imaginary-argument `Ei`) | ❌ (functions don't exist) | ✅ |
| `Fresnel S/C` (imaginary-argument `erf`) | ❌ (functions don't exist) | ✅ |
| **Sums / combinations** in one integrand | ❌ (single pattern) | ✅ (linear ansatz) |
| Special functions at an **inner tower level** | ❌ (`rt_special_case` is top-level only, `:4047`) | ✅ (part of the field recursion) |

The two structural gaps — (i) special functions only recognized on the *outermost* whole
integrand, and (ii) no proper-rational part over an inner field (`:2726-2728`) — are
**linked**. Cherry terms genuinely want to be produced *inside* `rt_field_integrate`, at
the same place the deferred tower Hermite/Rothstein–Trager belongs. Landing the deferred
elementary work and landing Cherry are naturally the same refactor of the field-level
integrator.

---

## 6. Development roadmap

Ordered by dependency and value. Each phase is independently shippable and preserves the
correct-by-construction + diff-back-verified contract.

### Phase 0 — Verification & scaffolding (prerequisite, low risk)
- Add the missing special-function **builtins** and, critically, their **derivative rules**
  so the diff-back gate can verify any new emission:
  `SinIntegral` (`sin(x)/x`), `CosIntegral` (`cos(x)/x`), `SinhIntegral`, `CoshIntegral`,
  `FresnelS` (`sin(πx²/2)`), `FresnelC` (`cos(πx²/2)`). Model on `erf.c` / `expintegralei.c`;
  wire derivative rules into `deriv.c` beside the existing block (`:339-392`); register in
  the symbol table with docstrings + attributes; add to `tests/CMakeLists.txt` COMMON_SRC
  (see `[[project_tests_common_src_list]]`).
- Add the connection identities as `.m` rules (`SinIntegral[x] = (Ei(ix)−Ei(−ix))/(2i)`
  etc.) so simplification can move between the `Ei`/`Erf` core and the trig siblings.
- **Insurance (independent):** gate `RischNorman`'s output with `Simplify[D−f]===0` before
  returning (§3 caveat).

### Phase 1 — Cherry error-function decision procedure
Generalize `rt_try_erf` (`:279-322`) from a single-Gaussian matcher into an **erf-augmented
hyperexponential ansatz**:
- In the exponential/`rt_hyperexp_case` path, after the Laurent-in-`θ` decomposition, add
  `erf` candidates for every Laurent exponent `kL` for which `−kL` is a perfect square in
  the field (argument `u = √(−kL)`), plus the elementary ansatz.
- Solve the coupled system with the existing `SolveAlways` extraction. Residuals fall back
  to `rt_solve_rde*`. Emit `v₀ + Σ cᵢ log vᵢ + Σ dⱼ erf(uⱼ)`; keep the `rt_verify_antideriv`
  gate as the safety net.
- Coverage unlocked: `∫ x e^(−x²) log x`, `∫ e^(−x²)/(x−a)` (as `erf` + elementary),
  `∫ e^(x²) dx` families with polynomial prefactors, and — via `K=I` — `FresnelS/C`.

### Phase 2 — Cherry logarithmic-integral decision procedure
Generalize `rt_try_ei` / `rt_try_li` (`:337-420`) into a **residue-style `li`/`Ei`
decision** on the top monomial:
- Log tower (`θ = log L`): detect `R(x)/log(w)` structure via the partial-fraction/residue
  path already used by `rt_frac_try`; the `li` argument is the closing residue.
- Exp tower (`θ = e^L`): `e^L·R/monomial` → `Ei`; reuse the Laurent machinery.
- Emit imaginary-argument images as `Si/Ci/Shi/Chi` using the Phase-0 identities so
  `∫ sin(x)/x`, `∫ cos(x)/x`, `∫ e^x/x^k` close in the natural named form.

### Phase 3 — Land the deferred field-level work + inner-tower special functions
- Implement the deferred **tower Hermite + Rothstein–Trager over K_{n-1}** and **general
  field Risch-DE** (`:2726-2728`). Add `ParametricLogarithmicDerivative` /
  `LogarithmicDerivativeOfARadical` as the missing RDE primitives.
- Fold the Phase 1/2 Cherry candidates into `rt_field_integrate` (`:3431-3494`) so special
  functions can be produced at inner tower levels, not just the outermost expression.
- This is the largest phase and subsumes the current top-level `rt_special_case` into the
  recursion.

### Phase 4 — Dilogarithm (Baddoura) and polish
- Promote `rt_try_dilog` to Baddoura's `Li₂` decision procedure.
- Close the trig **real-form simplification gap** (`:3975-3982`) so `∫ tan x` etc. come back
  real; this is a `Simplify` half-angle/log-of-product task, not a Risch task.

---

## 7. Concrete first steps (recommended)

1. **Phase 0 builtins** — highest value-to-risk. `Si/Ci/Shi/Chi/Fresnel` with derivative
   rules unblock every subsequent phase and are individually testable. Start with
   `SinIntegral`/`CosIntegral` (derivatives `sin x/x`, `cos x/x`).
2. **Phase 1 erf ansatz** — the biggest capability jump for the least new theory, because
   the Laurent + `SolveAlways` + RDE substrate already exists; we are adding basis terms,
   not an engine.
3. Defer Phase 3 until 1–2 are proven; it is the deepest refactor and is gated on the same
   field-level code as the already-deferred elementary work.

Throughout: keep the `rt_verify_antideriv` diff-back gate on every emission (it is the
reason this subsystem is trustworthy), add corpus cases to
`tests/test_integrate_risch_transcendental.c` (`test_special_functions`, `:467-498`) with
`assert_rm_num` numeric diff-back at interior points, and update
`docs/spec/builtins/` + the weekly changelog per `CLAUDE.md`.

---

## 8. Risks & watch-items

- **pmint runs first and is unverified** — extend the diff-back gate to it (§3). Do *not*
  put Cherry logic in pmint; it belongs in the verified recursive engine.
- **Argument-search cost** — Cherry's "is `−kL` a perfect square / does the residue close?"
  tests must reuse existing exact primitives (`Factor`, resultants, `SolveAlways`) and stay
  cap-free per `[[feedback_no_arbitrary_caps_decision_procedures]]`; watch for expression
  swell in the coupled solves.
- **Branch/real-form correctness** — imaginary-argument images (`Si/Ci/Fresnel`) must carry
  the right branch; verify via `rt_verify_antideriv` and real-interior `assert_rm_num`, not
  `PossibleZeroQ` (see `[[project_possiblezeroq_decay_false_positive]]`).
- **ABI skew** — new `special_functions/*.c` + `expr`-adjacent changes: `make clean` between
  incremental builds to avoid phantom segfaults
  (`[[project_matrix_stale_build_abi_skew]]`).

---

### Reference map (for the implementer)

- Recognizers to generalize: `integrate_risch_transcendental.c:279-508`
- Ansatz substrate to reuse: `rt_hyperexp_case:1719-1948`, `rt_solve_rde*:793-977`,
  `rt_field_hyperexp_coupled:3243-3425`
- Deferred field work (Phase 3): `:2726-2728`, `rt_field_integrate:3431-3494`
- Derivative rules to extend: `deriv.c:339-392`
- pmint `getSpecial` (currently Tan/Tanh/LambertW only): `intrischnorman.c:3213-3249`
- Dispatch position: `integrate.c:805` (RischNorman), `:811` (RischTranscendental)
