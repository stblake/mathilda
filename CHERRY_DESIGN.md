# Design: Incorporating Cherry's Special Functions into the Transcendental Risch Integrator

**Date:** 2026-07-14
**Targets:** G. W. Cherry, *Integration in Finite Terms with Special Functions: The
Logarithmic Integral* (SIAM J. Comput. 1986) → `LogIntegral` / `li`; and *An Analysis
of the Rational Exponential Integral* (SIAM J. Comput. 1989) → `ExpIntegralEi` / `ei`
and `Erf`/`Erfi`. Both PDFs in the repo root.
**Subject:** `src/calculus/integrate_risch_transcendental.c` and the `Risch\`` /
`Integrate\`` foundation modules (`risch_canonical`, `risch_structure`, `risch_field`,
`intrat` LRT).
**Status of prerequisites:** P0 (canonical representation) and P1 A+B (structure-theorem
ℚ-decision + `LogarithmicDerivativeOfRadical`) are **on `main`** — these are exactly
Cherry's foundation. This document is the architecture + the refactors to do *now* so the
Cherry engines drop in cleanly later. It is a design, not an implementation.

> **Superseded for execution by [`CHERRY_PLAN.md`](CHERRY_PLAN.md) (2026-07-16)**, which carries the
> code-level algorithm specs, the paper-exact test pins, and the corrected file map. This document
> predates the `risch_*` module split, so some file locations below are stale — the Status section
> (§3) and `CHERRY_PLAN.md` §1.2 give the current locations.

---

## 0. The one idea

Both Cherry papers prove an **extended Liouville theorem** of the identical shape:

```
∫ γ dx  =  v  +  Σ_i  k_i · SF_i(a_i)                         (†)
```

- `v` — the **elementary part** (an element of the differential field, i.e. the ordinary
  Risch answer: polynomial-in-tower + log terms).
- `SF_i` — a **special function** from a fixed menu (`li`, `ei`, `erf`, and later
  dilog/Fresnel/Si/Ci — all reducible to `ei`/`erf`, §7).
- `a_i` — a **generated argument** (a specific field/algebraic element), produced by a
  *finite* candidate generator unique to each `SF`.
- `k_i` — an **undetermined constant** (possibly with a fixed algebraic prefactor such as
  `e^{-α_i}` or `e^{-β_i}`).

Differentiating (†) gives the **matching identity** the algorithm actually solves:

```
γ  =  D_tower[v]  +  Σ_i  k_i · T_i(a_i)                       (‡)
```

where `T_i(a_i) = d/dx[SF_i(a_i)] / k_i` is the **derivative template** — a *rational
function of the tower variables* once `a_i` is fixed. The whole problem is: **generate the
finite set of candidate `a_i`, then solve (‡) as a linear system over the constants
`{coeffs(v), k_i}`.** If no constants work, *no antiderivative of this form exists* — the
decision property.

This is precisely the structure of the existing `rt_tower_solve` (build an undetermined-
coefficient ansatz `Q`, form `D_tower[Q] − F`, clear denominators, `SolveAlways` over the
constants). **Cherry = `rt_tower_solve` with extra basis terms `Σ k_i T_i(a_i)` whose
arguments come from new candidate generators.** That single observation drives the design.

### The two instantiations of (†)/(‡)

| | **1986 — Log Integral** | **1989 — Rational Exponential** |
|---|---|---|
| Integrand class | any transcendental-elementary `γ`, top monomial a **log** | `γ = g·e^f`, `g,f ∈ F(x)` (an **exp** top monomial `θ=e^f`) |
| `SF` menu | `li(u)` | `ei(u)`, `erf(ũ)` |
| Extended-Liouville form | `γ = w₀' + Σcᵢwᵢ'/wᵢ + Σdᵢuᵢ'/vᵢ`, `vᵢ=log(uᵢ)` (Thm 2.2, eq 2.1) | `g = y' + f'y + Σcᵢ(uᵢ'/uᵢ) + Σdᵢũᵢ'` (Thm 2.3, eq 2.2) |
| Argument constraint | `uᵢ = ∏ fⱼ^{αᵢⱼ}` (products of irreducible factors) | `uᵢ = f+αᵢ` (ei); `ũᵢ² = f+βᵢ` (erf) |
| **Argument generator** | **Σ-decomposition** (§4, Thm 4.4) | ei: **`Res_x(g₁, p+αq)`** (RT analogue, §4.2 P1); erf: **completing-square** `β²−Uβ+V` (§3, §4.3 R1/R2) |
| Derivative template `T(a)` | `u'/log(u)` | ei: `f'/(f+αᵢ)` (rel. to `e^f`); erf: `dᵢ ũᵢ'` (rel. to `e^f`) |
| Answer term `k·SF(a)` | `dᵢ li(uᵢ)` | `cᵢ e^{-αᵢ} ei(uᵢ)`, `dᵢ e^{-βᵢ} erf(ũᵢ)` |
| Coefficient solve | Risch Main Thm part (b) = linear over C | undetermined coeffs E1/E2/E3 = linear over C |
| Structural bound | Thm 4.4 termination (degree gate + monotone tails) | **≤ 2 erf args, quadratic** (Thm 3.2); exactly one "q-side" ei term (P2) |
| Recursion | rank-descending tower (Thm 5.4: peel top log/exp) | single exp level; y solved by multiplicity descent |

---

## 1. What we already have (the mapping is near 1:1)

| Cherry sub-procedure | Existing Mathilda facility | File |
|---|---|---|
| Partial-fraction split by top monomial (Lemma 5.1 / eq 4.2) | canonical rep `f = f_p + f_s + f_n` | `risch_canonical.c` (`Risch\`CanonicalRepresentation`) |
| ordinary-log part `Σcᵢwᵢ'/wᵢ`; `∫P/Q` (footnote 3) | Rothstein–Trager / LRT log part | `intrat.c` (`Integrate\`TranscendentalLogPart`), `rt_frac_lrt`, `rt_field_lrt_logpart` |
| **ei-argument resultant** `Res_x(g₁, p+αq)` | the *same* parametric RT resultant, differently parametrized | LRT resultant in `intrat.c` |
| Repeated-factor / multiplicity descent (Thm 5.4a; E1/E2) | Hermite reduction | `risch_hermite.c` |
| `vᵢ = rᵢaₙ + Σrᵢⱼaⱼ + …` (Roca79 structure thm) | **P1 structure oracle** | `risch_structure.c` (`LogReducible`/`ExpReducible`/`LogarithmicDerivativeOfRadical`) |
| final `B, cᵢ, dᵢ` linear solve (Risch69 part b) | undetermined-coeff ansatz + `SolveAlways` | `rt_tower_solve` in `integrate_risch_transcendental.c` |
| rank-descending tower recursion | `RtTower` + `rt_field_integrate` recursion | `integrate_risch_transcendental.c` |
| squarefree / irreducible factorization, gcd, Bézout, PF, resultant | poly subsystem (+ FLINT) | `poly/`, `parfrac.c`, `subresultants.c` |

**Genuinely new machinery Cherry needs** (nothing else):
1. **Σ-decomposition engine** (1986 §4) — the `li`-argument generator.
2. **Completing-square β-finder** (1989 §3/§4.3, Route B) — the `erf`-argument generator.
3. **Polynomial perfect-square test + polynomial square-root** (`q = s²`, `rᵢ² = p+βᵢq → rᵢ`).
4. A thin **ei-argument driver** re-parametrizing the existing RT resultant as `Res_x(g₁,p+αq)`.

The four existing recognizers (`rt_try_erf/ei/li/dilog`) are exactly Cherry's *narrow enumerated
cases* — one Gaussian, one `Ei`, one `li` monomial, one dilog template — each already in the
correct **`template → diff-back gate`** frame. Cherry replaces *template* with *decision
procedure* (finite generator + bound), keeping the gate and **adding the decision property**
(prove non-existence).

---

## 2. Target architecture

Three components, layered on the existing dispatch. `rt_integrate` today is
`rational → transcendental → special`; the special stage becomes registry-driven.

### 2.1 `SpecialFunctionForm` registry (the load-bearing extension seam)

A table of special-function "plugins", each a struct of function pointers:

```c
typedef struct {
    const char* name;                       /* "Erf", "ExpIntegralEi", "LogIntegral", ... */

    /* Applicability: does the top monomial kind admit this SF?  (Li ⇐ log top;
     * Ei/Erf ⇐ exp top.)  Cheap pre-filter before generation. */
    bool (*applicable)(const RtTower* T, long top);

    /* Argument generator: the FINITE candidate set {a_i} for this integrand.
     *   Li  -> Sigma-decomposition (1986 Thm 4.4)
     *   Ei  -> Res_x(g1, p + alpha q)          (1989 §4.2 P1)
     *   Erf -> completing-square beta finder   (1989 §3 / §4.3)
     * Returns count; fills args[] (owned) and per-arg fixed prefactors pref[]
     * (e.g. e^{-alpha}, e^{-beta}, or 1). NULL/0 => this SF contributes nothing. */
    size_t (*gen_arguments)(Expr* gamma, const RtTower* T, long top, Expr* x,
                            Expr*** args_out, Expr*** pref_out);

    /* Derivative template T_i(a): the field expression contributed to gamma per unit
     * coefficient k_i (eq ‡).  Li: u'/log(u).  Ei: f'/(f+alpha) (rel. e^f).
     * Erf: d ua'  (rel. e^f). */
    Expr* (*deriv_template)(Expr* a, const RtTower* T, Expr* x);

    /* Answer constructor: builds  k * pref * SF(a). */
    Expr* (*answer_term)(Expr* k, Expr* pref, Expr* a);

    /* Structural bound for the ansatz (max terms): Erf -> 2; Ei/Li -> from generator. */
    int max_terms;
} SpecialFunctionForm;
```

The registry is a static array; the current recognizers become entries (§4, migration).
Adding a new special function later = adding one entry + its generator, with **zero changes
to the driver**.

### 2.2 `extended_liouville_solve` driver (generalizes `rt_tower_solve`)

```
extended_liouville_solve(gamma, x, T, top):
  1. forms = [ sf in REGISTRY if sf.applicable(T, top) ]
  2. cand  = concat over forms of sf.gen_arguments(gamma, T, top, x)     # finite
     (bail to NULL early if a form's constant-existence pre-test already fails —
      that is Cherry's "no integral of this form" decision, e.g. P1/P2/R1/R2 tests)
  3. Build the ansatz basis:
       - elementary part v: the SAME undetermined-coefficient polynomial/log ansatz
         rt_tower_solve already builds (degrees/multiplicities from the canonical
         rep + Hermite structure — NOT capped);
       - one basis term k_c * T_c(a_c) per candidate argument a_c.
  4. gamma - D_tower[v + Σ k_c T_c(a_c)]  -> Together -> Numerator == 0
     -> SolveAlways over { coeffs(v), k_c } .                # linear over C
  5. if solved: answer = v* + Σ k_c* pref_c * SF_c(a_c),  back-substitute tower->kernels,
                diff-back gate against ORIGINAL gamma (kept — correct by construction).
     else:      return NULL  (decision: no antiderivative of this form).
```

Step 4 is *literally* `rt_tower_solve` with an enlarged basis. The elementary part `v` reuses
the existing polynomial/log ansatz construction; the special terms just add columns to the
same linear system. **This is why the refactor in §3.1 (extract the ansatz-solver) is the
keystone.**

### 2.3 Argument generators (the three new engines)

- **`cherry_sigma_decomp.c`** — Σ-decomposition (1986 Thm 4.3/4.4). Input: `Φ ∈ K(x)`, an
  irreducible sequence `Σ=(f₁,…,fₘ)` (from the denominator's factorization), a degree bound
  `d`. Output: the finite set of `(bᵢ, exponent-vector αᵢ)` giving `li` arguments
  `uᵢ = ∏ fⱼ^{αᵢⱼ}`, or a proof of non-existence. Steps (all reuse existing poly ops):
  multiplicity read-off (Lemma 4.1 = squarefree/irreducible factorization),
  `bᵢ = (p mod f₁)/(q mod f₁)` (modular reduction in `K[x]/(f₁)`), interpolation-based degree
  gate, monotone-tail termination. **Stated for arbitrary degree `d`** — `li` needs `d=1`; the
  dilogarithm/higher SFs reuse the same engine at `d≥2` (Cherry's explicit hook, 1986 p.7).
- **`cherry_completing_square.c`** — the `erf` β-finder (1989 Thm 3.2, Route B). Input: `f=p/q`.
  Gate: `q` a perfect square (`q=s²`) else **0 erf terms**. Then `p² = R·h² + S·q` (poly
  division), `−S = U·p + V·q` (undetermined constants), roots of `β²−Uβ+V` are the ≤2 `β`;
  extract `rᵢ = √(p+βᵢq)`; `ũᵢ = rᵢ/s`. Needs the new poly-square-root primitive.
- **ei generator** — a thin wrapper re-parametrizing the existing LRT resultant as
  `Res_x(g₁, p+αq)` (1989 §4.2 P1): roots `αᵢ` → `uᵢ = f+αᵢ`, weight `e^{-αᵢ}`; PF-match
  `g̃₁/g₁ = Σ p̃ᵢ/pᵢ` against `Σ cᵢ pᵢ'/pᵢ` for the constant-existence test; the single "q-side"
  term `c = −γ_∞/∂q` (P2). Reuses `intrat`'s resultant + the P1 structure oracle for the
  existence test.

### 2.4 New polynomial primitives (small, generally useful)

- `PolynomialSqrt[p]` / `[p, x]` — the polynomial square root when it exists, with the
  perfect-square test folded in (returns `$Failed` otherwise). **Landed** as
  `builtin_polynomialsqrt` in `src/poly/facpoly.c` (static helper `ps_is_numeric`; see R5). A
  separate `PolynomialPerfectSquareQ` / `poly_perfect_square_q` predicate proved unnecessary — the
  square-root builtin subsumes it. Self-contained addition to `poly/`; useful beyond Cherry.

---

## 3. Refactors to do NOW (make the incorporation clean, zero behaviour change)

These are behaviour-preserving substrate moves that turn today's ad-hoc special layer into
the registry/driver seam. Each is independently landable and testable (byte-identical
outputs), the same discipline as the P0 spine.

1. **Extract the ansatz solver.** Pull the "build undetermined-coefficient basis → form
   `D_tower[·] − F` → `SolveAlways` → substitute free params → back-rules" core out of
   `rt_tower_solve` into a reusable `rt_ansatz_solve(basis_terms, unknowns, F, T, x)` that
   accepts **extra basis terms** (the future `Σ k_c T_c(a_c)`). Today's callers pass an empty
   extra-basis list → identical output. *This is the keystone; everything else composes onto it.*
   **(Superseded by R1 below: no extraction was done — `rt_tower_solve` already IS this solver and
   accepts SF basis terms directly, so `rt_ansatz_solve` was never created.)**

2. **Promote `rt_special_case` to a registry.** Wrap `rt_try_erf/ei/li/dilog` as four
   `SpecialFunctionForm` entries behind the four-method interface (each still a narrow
   generator internally). Dispatch becomes a loop over the registry. No behaviour change; it
   establishes the seam so generalized generators swap in one at a time.

3. **Expose the parametric RT resultant.** Factor the LRT resultant
   (`Res_x(den, num − α·den')`) in `intrat.c` into a named primitive
   `rt_parametric_resultant(A, B, α, x)` so the ei generator can call `Res_x(g₁, p+αq)` and
   the erf Route-A fallback can call `Res_x(p+βq, p'+βq')`. Pure extraction.

4. **Give the driver a "lower-field integrate" hook.** Cherry's recursion peels the top
   monomial and integrates the remainder over the lower field (Thm 5.4). `rt_field_integrate`
   already recurses; expose `rt_integrate_lower(remainder, T, top-1, x)` so the SF driver can
   delegate the elementary/lower part cleanly instead of re-deriving it.

5. **Land the two poly primitives** (`poly_perfect_square_q`, `poly_sqrt`) as standalone,
   unit-tested builtins now — they gate the erf path and are independently useful.

None of 1–5 changes any current integral; together they make Cherry a matter of *adding
generators + registry entries*, not surgery.

### Status (2026-07-14) — the five refactors, as landed

Surveying the code sharpened three of these to avoid adding dead abstraction ahead of a
consumer (the project's "hack-free / simplicity" bar):

- **R5 — DONE (new code + tests).** `PolynomialSqrt[p]` / `[p, x]` in `src/poly/facpoly.c`
  (`tests/test_polynomialsqrt.c`). Returns `s` with `s²==p` for a perfect square (even
  multiplicities; numeric content via `Sqrt`), else `$Failed`; accepted only after the exact
  `Expand[s²−p]==0` certificate. This is the erf-argument square-root primitive (`rᵢ=√(p+βᵢq)`).
- **R2 — DONE (registry).** `RtSpecialForm` table in `src/calculus/risch_special.c` — the module
  was split after this doc was written, so the registry and the four recognizers moved out of
  `integrate_risch_transcendental.c` (now a ~710-line driver). The four recognizers are entries and
  `rt_special_case` loops the registry (`risch_special.c:262`). Byte-identical behaviour (verified).
  The struct is the seam that grows `gen_arguments`/`deriv_template`/`answer_term`/`max_terms`.
- **R3 — DONE (thin primitive over existing `Resultant`).** `Integrate\`RothsteinTragerResultant[
  num, den, z, x] = Resultant[num − z D[den,x], den, x]` in `src/calculus/intrat.c`
  (`tests/test_rt_resultant.c`) — roots in `z` are the residues; the Ei/log argument generator.
  (The `Resultant` builtin already exists; the ei case `Res_x(g₁,p+αq)` is a plain `Resultant`
  call. This names the RT form for reuse.)
- **R1 — DONE (contract, no signature change).** `rt_tower_solve` already IS the
  reduce-to-a-linear-system-over-constants solver and already accepts special-function basis
  terms in `Q` (Mathilda's `D[]` differentiates `LogIntegral`/`ExpIntegralEi`/`Erf`
  correctly), so it needed a documented extension-point contract, not an extraction. The
  Cherry driver appends `Σ kᵢ SF(aᵢ)` to `Q` and the `kᵢ` to `syms` and calls it unchanged.
- **R4 — DEFERRED to the Cherry module (intentionally).** The lower-field recursion hook is
  `rt_field_integrate` / `rt_transcendental_case`, both internal and already the recursion
  entry. A standalone wrapper now has no consumer until the Cherry driver exists, so exposing
  it then (with the driver) avoids dead code. No change now.

---

## 4. Migration of the existing recognizers

Each current template becomes the *degenerate case* of its Cherry generator, verified to still
fire on its old inputs:

| Recognizer | Becomes | Cherry generalization |
|---|---|---|
| `rt_try_erf` (`f'/f` deg-1, one Gaussian) | `Erf` form, `gen_arguments` via completing-square | ≤2 quadratic `erf` args over any `f=p/q` with `q=s²` |
| `rt_try_ei` (linear exp / linear denom) | `ExpIntegralEi` form, `gen_arguments` via `Res_x(g₁,p+αq)` | multi-term `ei`, algebraic `αᵢ`, any rational `g` |
| `rt_try_li` (`c wᵖ⁻¹w'/Log[w]`) | `LogIntegral` form, `gen_arguments` via Σ-decomposition | multi-`li`, product arguments `∏fⱼ^{αⱼ}` |
| `rt_try_dilog` (`K Log[1+px]/x`) | `PolyLog[2,·]` form (future) | degree-2 Σ-decomposition (Cherry's `d≥2` hook) |

Keep each old recognizer as a **fast-path** inside its form's `gen_arguments` (cheap template
match first, full generator second) so common cases stay O(1) and nothing regresses.

---

## 5. Phased sequencing

Dependencies: refactors (§3) → generators → driver wiring. Ship each phase green.

- **C0 — Substrate refactors (§3.1–3.5).** Behaviour-preserving. Keystone = ansatz-solver
  extraction. Land the two poly primitives with tests.
- **C1 — `ei` (1989 exponential, no erf).** Cheapest generator (reuses the RT resultant);
  highest hit rate. Registry entry + `Res_x(g₁,p+αq)` generator + P1/P2 constant tests +
  the E1/E2/E3 `y` solve (undetermined coeffs). Pins: Ex 5.1, 5.3, and the p.894 `e^x/(x²−2)`.
- **C2 — `erf` (1989 completing-square).** Adds `cherry_completing_square.c` + poly-sqrt.
  ≤2-arg bound. Pins: Ex 5.2, 5.4.
- **C3 — `li` (1986 Σ-decomposition).** Adds `cherry_sigma_decomp.c` (degree `d=1`) + the
  rank-descending recursion (Thm 5.4). Pins: Ex 1.1, 5.1, the Macsyma d1/d3/d4, and the
  **non-elementary decisions** Ex 5.2 / d12.
- **C4 — reductions & dilog.** `si/ci/Fresnel` via `ei/erf` linear combos (1986 §6, 1989 App);
  dilogarithm via degree-2 Σ-decomposition. Retire the four narrow recognizers.

Each of C1–C3 is *correct by construction* (diff-back gate retained) **and** gains the
**decision property** (the constant-existence tests + structural bounds prove non-existence),
which the current recognizers lack.

---

## 6. Test pins (from the papers — exact I/O)

**1986 (li):**
- `∫ x/log(x)² dx = 2 li(x²) − x²/log(x)` (Ex 1.1 / Macsyma d1)
- `∫ x²/log(x+1) dx = li(x³+3x²+3x+1) − 2 li(x²+2x+1) + li(x+1)` (d3)
- `∫ x³/log(x²−1) dx = ½ li(x⁴−2x²+1) + ½ li(x²−1)` (Ex 5.1)
- `∫ 1/(log(x)+3) dx = e^{−3} li(e³x)` (d2 — transcendental-constant rescaling)
- `∫ (log(x)²+3)/(log(x)²+3log(x)+2) dx = −7e^{−2} li(e²x) + 4e^{−1} li(ex) + x` (d4)
- **decision NO:** `∫ x²/log(x²−1) dx` (Ex 5.2), `∫ (x²+1)/(x²+x+1) eˣ dx` under alg-closed C (d12)

**1989 (ei/erf):**
- `∫ e^{1/x} dx = x e^{1/x} − ei(1/x)` (Ex 5.1)
- `∫ (1/x + 1/x²) e^{1/x²} dx = −½ ei(1/x²) − erf(1/x)` (Ex 5.2)
- `∫ (x⁴−1)/(x⁵+x) e^{(x²+1)/x} dx = ½e^{√2} ei((x²−√2x+1)/x) + ½e^{−√2} ei((x²+√2x+1)/x)` (Ex 5.3)
- `∫ e^{(x⁴+a)/x²} dx = ½e^{−2√a} erf((x²+√a)/x) + ½e^{2√a} erf((x²−√a)/x)` (Ex 5.4)
- `∫ e^x/(x²−2) dx = e^{√2} ei(x+√2) + e^{−√2} ei(x−√2)` (p.894 — new algebraic constant)

Note Mathilda's convention differences (`ei ↔ ExpIntegralEi`, `li ↔ LogIntegral`, Cherry's
`erf ↔ Erfi`); the answer constructors must map to Mathilda's `Erf`/`Erfi` (Cherry
`Erf(x)=(1/i)Erfi(x)`), `ExpIntegralEi`, `LogIntegral`, `PolyLog`.

---

## 7. Non-goals / risks

- **Algebraically-closed constant field.** Both theorems assume `C = C̄`; the arguments carry
  algebraic constants (`√2`, `√a`, quadratic `β`). The engine must represent/solve over
  `F̄` — `RootReduce`/`qqbar` + the algebraic-tower normalizer already exist and are the tool
  (this is the single biggest engineering surface; Cherry's d12 deficiency is exactly a
  not-alg-closed miss).
- **Not the general Risch DE.** 1989 deliberately avoids it (E1/E2/E3 undetermined coeffs);
  do not route the `y` solve through `rde_base`.
- **Soundness is free, completeness is the work.** The diff-back gate keeps every emitted
  answer correct; incompleteness only declines. The decision property (proving NO) is the new
  guarantee and rests on the structural bounds — those must be implemented faithfully
  (≤2 erf, Σ-decomp termination) or the "NO" answers are unsound-as-a-decision (though still
  never a wrong antiderivative).
- **Keep the registry small and total.** Every `SpecialFunctionForm` must declare its bound
  and its constant-existence test; a generator that can loop (unbounded candidate set) breaks
  the decision property.

---

## 8. One-line summary

Cherry is not a rewrite: it is `rt_tower_solve` with special-function basis terms, fed by
three new finite argument-generators (Σ-decomposition, completing-square, α-resultant), behind
a `SpecialFunctionForm` registry — sitting on the canonical-rep + structure-oracle + LRT +
Hermite foundation that P0/P1 already put on `main`. Do the five §3 refactors now; add
generators one paper-example at a time.
