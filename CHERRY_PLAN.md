# CHERRY_PLAN.md — Implementation Plan: Cherry's Special-Function Risch Extensions

**Date:** 2026-07-16
**Targets:** `ExpIntegralEi` (ei), `Erf`/`Erfi` (erf), `LogIntegral` (li) as *decided* outputs of
the recursive transcendental Risch integrator `Integrate\`RischTranscendental`.
**Sources:**
- G. W. Cherry, *An Analysis of the Rational Exponential Integral*, SIAM J. Comput. 18(5), 1989,
  pp. 893–905 → **ei, erf** (`Cherry Integration of Exponentials.pdf`).
- G. W. Cherry, *Integration in Finite Terms with Special Functions: The Logarithmic Integral*,
  SIAM J. Comput. 15(1), 1986, pp. 1–21 → **li** (`Cherry Log Integrals.pdf`).
**Companion:** `CHERRY_DESIGN.md` (2026-07-14) — the architecture; this document is the executable
sequence and the paper-faithful algorithm specs. Where the two differ on file locations, **this
document is authoritative** (the module was split after CHERRY_DESIGN.md; see §1.2).
**Scope (confirmed):** full transcendental-tower recursion (Thm 5.4), all three engines at full
depth (C1 ei → C2 erf → C3 li), algebraically-closed constants (√2, quadratic β) handled from
phase 1.

---

## 1. Context

### 1.1 Why

The current special-function layer is four narrow template matchers in
`src/calculus/risch_special.c`:

| Recognizer | Line | Matches | Emits |
|---|---|---|---|
| `rt_try_erf`   | `:61`  | `K·E^(ax²+bx+c)` (log-derivative `f'/f` is degree-1 poly) | `Erf` |
| `rt_try_ei`    | `:119` | `M·E^(ax+b)/(cx+d)` (linear exponent, linear denom)       | `ExpIntegralEi` |
| `rt_try_li`    | `:166` | `c·w^(p−1)·w'/Log[w]` (loops `p=1..6`)                     | `c·LogIntegral[w^p]` |
| `rt_try_dilog` | `:208` | `K·Log[1+px]/x`                                            | `−K·PolyLog[2,−px]` |

Each is correct where it fires (diff-back gated against a real `D[]` rule) but fires on a
vanishingly small set of shapes and — decisively — **cannot prove non-existence**. Cherry replaces
each template with a **finite argument-generator + structural bound**, which yields:
- *completeness* — multi-term outputs, algebraic arguments (√2), product arguments (`∏ fⱼ^{αⱼ}`);
- *the decision property* — a constant-existence pre-test or a terminating bound proves "no
  antiderivative of this special-function class exists" and declines soundly.

### 1.2 Verified current state (reconciles CHERRY_DESIGN.md against the tree)

CHERRY_DESIGN.md's *claims of existence* are correct, but the module split moved things; several
of its file/name references are stale. Corrected map:

| Facility | CHERRY_DESIGN.md | **Actual (verified 2026-07-16)** |
|---|---|---|
| `RtSpecialForm` registry + 4 recognizers | `integrate_risch_transcendental.c` | **`risch_special.c`** — struct `:249`, `RT_SPECIAL_FORMS[]` `:254`, `rt_special_case` loop `:262`. `integrate_risch_transcendental.c` is now a 710-line driver. |
| poly square-root | `poly_sqrt` / `poly_perfect_square_q` | **`PolynomialSqrt`** = `builtin_polynomialsqrt` (`src/poly/facpoly.c:99`, reg `:188`); static `ps_is_numeric`; test `tests/test_polynomialsqrt.c`. No standalone predicate. |
| RT resultant | `Integrate\`RothsteinTragerResultant` | `builtin_intrat_rt_resultant` (`src/calculus/intrat.c:3925`, reg `:4058`) = `Resultant[num − z·D[den,x], den, x]`; test `tests/test_rt_resultant.c`. |
| ansatz solver | `rt_ansatz_solve` | **`rt_tower_solve`** only (`src/calculus/risch_field_integrate.c:354`). SF-basis capability is **latent** (header contract `:342–353`): accepts `Q = poly_ansatz + Σ kᵢ·SF(aᵢ)` with `kᵢ` in `syms`; `D[]` already differentiates li/ei/erf. The generator wiring that *fills* `aᵢ` does not exist yet. |
| diff-back verifier | (implied) rejects residual `Integrate` | **`rt_verify_antideriv`** (`src/calculus/risch_util.c:317`) does **NOT** reject an unevaluated `Integrate[…,x]` head. |

**Reusable substrate (do not rebuild):** literal tower Hermite `risch_hermite.c` +
partial-log-part residue split (`intrat_log_part_core`, `intrat.c:1685`); `rde_base`/`rde_tower`
(`integrate_risch_rde.c`, 1124 LoC — but see §3, we do **not** route the y-solve through it);
LRT `rt_field_lrt_logpart` → `intrat_log_part_core`; `rt_field_integrate` recursion +
`rt_tower_build_min` (`risch_tower.c`); `rt_hyperexp_case` Laurent split + `rt_exp_kernelize`
commensurate-exponent collapse; `rt_rde_var_bound` degree bound.

**Phase 0 is already complete** — every derivative rule and builtin exists:
`Erf/Erfi/Erfc` (`deriv.c:382/391/400`; builtins `erf.c:489`, `erfi.c:568`, `erfc.c:440`),
`ExpIntegralEi` (`:427`; `expintegralei.c:680`), `LogIntegral` (`:433`; `logintegral.c:146`),
`PolyLog` (`:1239`; `polylog.c:918`), `SinIntegral/CosIntegral/SinhIntegral/CoshIntegral`
(`:340–358`), `FresnelC/FresnelS` (`:365/369`).

### 1.3 Three A4 hazards baked into every engine

From `RISCH_AUDIT_A4.md` (2026-07-15):
1. **Decline → `RT_DEC_UNKNOWN`** (→ unevaluated `Integrate`), **never** `RT_DEC_NONELEMENTARY`
   (→ `False`), unless a genuine Cherry certificate fired. `rt_dec_nonelem`
   (`risch_field_integrate.c:696`) is write-once; call it *only* at authoritative Cherry decline
   points (P1/P2 PF-match failure, Σ-decomp monotone-tail NO), never at a routing/out-of-scope
   decline. (A4 Finding 2.)
2. **Reject residual `Integrate` heads in the Cherry emission gate** (A4 Finding 5). Because
   `D[Integrate[f,x],x]→f`, a Cherry candidate that left a residual integral would diff-back to
   zero and be rubber-stamped. **Design correction (2026-07-16):** the guard is placed in the
   Cherry engine's emission gate — `rt_free_of_head(Q,"Integrate") && rt_verify_antideriv(...)` —
   **not** in the shared `rt_verify_antideriv`. A blanket guard there would regress an *intended*
   feature: `risch_singleext.c:548–563` deliberately assembles `logs + Integrate[remainder, x]`
   (the partial-log-part emission, A4 Finding 4 = not a defect) and relies on the FTC close to pass
   the same verifier. The strictness belongs where an open integral is illegitimate (Cherry answers
   are complete antiderivatives), not in the shared gate.
3. **Certify only with exact `rt_verify_antideriv`** (`Together∘TrigToExp` / `Simplify` zero-test).
   **Never** numeric sampling — A4 Finding 6 deleted the last sampler; do not reintroduce one.

---

## 2. The one idea, and the architecture

### 2.1 The matching identity

Both papers prove an extended-Liouville theorem of one shape:

```
∫ γ dx  =  v  +  Σ_i  k_i · SF_i(a_i)                                (†)
```

- `v` — the **elementary part** (ordinary Risch answer: tower-polynomial + logs).
- `SF_i(a_i)` — a special function (`li`/`ei`/`erf`) at a **generated argument** `a_i`.
- `k_i` — an undetermined constant, possibly with a fixed prefactor `e^{−α_i}` / `e^{−β_i}`.

Differentiating (†) gives the identity the algorithm actually solves:

```
γ  =  D_tower[v]  +  Σ_i  k_i · T_i(a_i),   T_i(a_i) = d/dx[SF_i(a_i)] / k_i          (‡)
```

`T_i(a_i)` is a **rational function of the tower variables once `a_i` is fixed**:
- li: `T = u'/log(u)`
- ei: `T = f'/(f+α_i)` (relative to the `e^f` factor)
- erf: `T = d_i·ũ_i'` (relative to `e^f`)

So **Cherry = `rt_tower_solve` with extra basis terms `Σ k_i T_i(a_i)`**, whose arguments come from
three new *finite* generators. Solve the linear system over `{coeffs(v), k_i}`; if it has no
solution, no antiderivative of this form exists.

### 2.2 Grow the registry into a generator interface

Upgrade `RtSpecialForm` (`risch_special.h`) from `{name, recognize}` to the four-method plugin
struct of CHERRY_DESIGN.md §2.1:

```c
typedef struct {
    const char* name;                    /* "ExpIntegralEi" | "Erf" | "LogIntegral" | ... */
    /* cheap pre-filter: does the top monomial kind admit this SF? (li ⇐ log top; ei/erf ⇐ exp) */
    bool  (*applicable)(const RtTower* T, long top);
    /* the FINITE candidate set {a_i} + per-arg fixed prefactor {pref_i} (e^{-α}, e^{-β}, or 1).
       Returns count; 0 => this SF contributes nothing. Early-declines on the paper's
       constant-existence pre-tests (P1/P2/R1/R2 / Σ-decomp monotone tail). */
    size_t (*gen_arguments)(Expr* gamma, const RtTower* T, long top, Expr* x,
                            Expr*** args_out, Expr*** pref_out);
    Expr* (*deriv_template)(Expr* a, const RtTower* T, Expr* x);   /* T_i(a) of (‡) */
    Expr* (*answer_term)(Expr* k, Expr* pref, Expr* a);           /* k * pref * SF(a) */
    int   max_terms;                                              /* erf -> 2; ei/li -> generator */
} RtSpecialForm;
```

Keep each existing recognizer (`rt_try_erf/ei/li/dilog`) as an **O(1) fast-path** at the top of its
form's `gen_arguments` — cheap template match first, full generator second — so common cases stay
constant-time and nothing regresses.

### 2.3 The `extended_liouville_solve` driver

New file `src/calculus/cherry_driver.c`; generalizes the `rt_special_case` loop:

```
extended_liouville_solve(gamma, x, T, top):
  1. forms = [ sf in RT_SPECIAL_FORMS if sf.applicable(T, top) ]
  2. cand  = concat over forms of sf.gen_arguments(gamma, T, top, x)          # finite
       - if any form's constant-existence pre-test FAILS -> that is Cherry's "no integral of
         this form" verdict; record it (see decision aggregation below), contribute 0 args.
  3. basis = v_ansatz(T, canonical rep, Hermite)          # SAME poly/log ansatz rt_tower_solve
                                                          # already builds — NOT capped
             + [ k_c * sf.deriv_template(a_c) for a_c in cand ]
  4. Q = basis;  syms = coeffs(v_ansatz) ++ [k_c ...]
     solved = rt_tower_solve(Q, syms, ...)               # gamma - D_tower[Q] -> Together
                                                          # -> Numerator == 0 -> SolveAlways
  5. if solved:
        answer = v* + Σ k_c* * pref_c * sf.answer_term(...)
        back-substitute tower -> kernels
        if rt_verify_antideriv(answer, gamma, x): return answer
     return NULL           # + decision aggregation (§ below)
```

Step 4 is **literally `rt_tower_solve` with an enlarged basis** — the special terms add columns to
the same linear system. This is why the latent SF-basis contract in `rt_tower_solve` (`:342–353`)
is the keystone: no new solver.

**Decision aggregation (A4 F2).** The driver distinguishes:
- *no candidate applies / routing decline* → return `NULL`, leave `g_rt_decision = UNKNOWN`;
- *a form's constant-existence certificate fired NO* (e.g. Σ-decomp monotone-tail proves the
  restricted decomposition cannot exist, or P1 PF-match yields non-constant residues over the
  full menu) → call `rt_dec_nonelem()` so `ElementaryIntegralQ` can answer `False`.

---

## 3. Engine C1 — `ExpIntegralEi` (Cherry 1989 §4.2)

**Class:** `γ = g·e^f`, `g, f ∈ F(x)` with `f = p/q`, `gcd(p,q)=1`, `q` monic.
**Integral form (Thm 2.3):**
`∫ g e^f dx = y·e^f + Σ_i c_i e^{−α_i} ei(u_i) + Σ_i d_i e^{−β_i} erf(ũ_i)` with `u_i = f + α_i`.
This engine produces `y` and the `ei` terms; the `erf` terms are C2 (§4). New file
`src/calculus/cherry_ei.c` (may share a compilation unit `cherry_exp.c` with C2).

**Structure theorem constraint (Thm 2.3):** `g = y' + f'y + Σ c_i(u_i'/u_i) + Σ d_i ũ_i'`, with
`u_i = f + α_i` (α_i algebraic over F) and `ũ_i² = f + β_i`.

### 3.1 The y-solve — E1/E2/E3 undetermined coefficients (Dav86), **not** `rde_base`

Cherry deliberately avoids the general Risch DE `y' + f'y = g`; a multiplicity-descent solves `y`
directly. Setup: `g = Ḡ/(g₁g₂)`, `g₂` = product of irreducible factors of `den(g)` dividing `q`,
`g₁` = the rest; PF `g = G + g̃₁/g₁ + g̃₂/g₂`.

- **E1 — fractional part of `y` from `g₁`.** Squarefree `g₁ = h_r^r ⋯ h₁`; PF
  `g̃₁/g₁ = ΣΣ G_{ij}/h_i^j`. Repeatedly substitute `y = ȳ − Y_{nn} h_n'/h_r`, solving the
  congruence `G_{rr} ≡ −n·Y_{nn}·h_n' (mod h_r)` **uniquely** for `Y_{nn}`, until `g₁` is
  squarefree (each irreducible factor of `den(y)` divides `q`). Set `n = 0`.
- **E2 — reduce multiplicities in `g₂`.** With `φ^{n+1} ‖ den(y')`, `φ^{m+1} ‖ den(f')`: find
  squarefree `h` with `h' ‖ g₂`; PF `g̃₂/g₂ = ĝ/h^r + ⋯`, `ỹ₂/y₂ = ŷ/h^{r−(m+1)} + ⋯`,
  `f' = f̂/h^{m+1} + ⋯`; solve the congruence `ĝ ≡ f̂·ŷ (mod h)` **uniquely** for `ŷ`; substitute
  `y = ȳ − ŷ/h^{r−(m+1)}` to reduce multiplicity. After finitely many steps `y` is a polynomial.
- **E3 — polynomial part.** `y = y_n xⁿ + ⋯ + y₀`; identity becomes
  `g = (y_n xⁿ+⋯)' + f'(y_n xⁿ+⋯) + Σ c_i(u_i'/u_i) + Σ d_i ũ_i'`. Cases:
  - (i) `∂p > ∂q` (`lc(f')≠0`): `lc(g) = lc(f')·y_n`, solve `y_n`, substitute `y = ȳ − y_n xⁿ`.
  - (ii) `∂p ≤ ∂q` (`∂f' ≤ −2`): subcase `n>0` ⇒ `lc(g) = n·y_n`; subcase `n=0` uses
    `h^{m+1} ‖ q`, PF `g = ĝ/h^r + ⋯`, `f' = f̂/h^{m+1} + ⋯`, comparison ⇒ `r = m+1`, `y₀ = ĝ/f̂`.

Each substitution strictly reduces a denominator multiplicity (E1/E2) or a degree (E3);
termination is structural, no linear system over C. Implement with Mathilda `Apart`,
`PolynomialGCD`, squarefree factorization, `PolynomialRemainder`, `PolynomialExtendedGCD`.

### 3.2 P1 — the ei argument generator (resolves `g̃₁/g₁ = Σ c_i u_i'/u_i`)

`u_i = f + α_i ⇒ p_i = p + α_i q`. Generate the α's as roots of the parametric resultant:

```
h(α) = Res_x(g₁, p + α q)            # Integrate`RothsteinTragerResultant re-parametrized
```

Factor `g₁ = ∏ p_i` from the roots `α_i`; do the PF `g̃₁/g₁ = Σ p̃_i/p_i` (Mathilda `Apart`).
**Constant-existence test:** each `p̃_i` must be a *constant* `c_i` (independent of x). If any `p̃_i`
is non-constant ⇒ **no integral in a special incomplete Γ-extension exists** — record the NO
certificate (§2.3) and decline.

*Re-parametrization note:* the existing `RothsteinTragerResultant[num,den,z,x] =
Resultant[num − z·D[den,x], den, x]`. Cherry's `Res_x(g₁, p+αq)` is a **plain** `Resultant[g₁,
p+α q, x]` in the parameter α — expose a thin sibling `Integrate\`ExpIntegralEiResultant[g1,p,q,α,x]
= Resultant[g1, p + α q, x]` in `intrat.c` beside the existing one (both are one-line `Resultant`
calls; naming makes the reuse explicit and testable).

### 3.3 P2 — the single "q-side" term (resolves `G + g̃₂/g₂`)

`p_i = p + α_i q ∈ F̄` in exactly two mutually exclusive cases: (i) `∂p = ∂q` and `α_i = lc(p)`;
(ii) `p ∈ F` and `α_i = 0`. Either way there is **at most one** such term `c·(q'/q)`. Transposing
the known `−(q'/q) Σ_{p_i∉F̄} c_i` to the left gives `γ = −c(q'/q) + Σ d_i(r_i/s)'` with `γ ∈ F(x)`
known; then if `q ≠ 1`, comparing proper-fraction coefficients:

```
c = −γ_∞ / ∂q
```

where `_∞` is the proper-fraction leading coefficient (`f_∞ = c_{r−1}` of `f_n/f_d`; `f_∞ = 0` if
`∂f_n < r−1`; `(f')_∞ = 0`; additive). Implement `_∞` as a small helper on the `Apart` form.

### 3.4 C1 test pins (exact I/O; `ei ↔ ExpIntegralEi`)

- **Ex 5.1** — `∫ e^{1/x} dx = x e^{1/x} − ei(1/x)`. (`q=x` not a square ⇒ 0 erf; y=x; P2 c=−1.)
- **Ex 5.3** — `∫ (x⁴−1)/(x⁵+x) e^{(x²+1)/x} dx`
  `= ½ e^{−√2} ei((x²+√2x+1)/x) + ½ e^{√2} ei((x²−√2x+1)/x)`.
  (`h(α)=α⁴−4α²+4 ⇒ α=±√2`; PF ⇒ `c₁=c₂=½`; algebraic constants — needs §7. **NB:** the paper's
  printed exponent signs are OCR-ambiguous on this line; the authoritative rule is weight
  `e^{−α_i}` with `u_i=f+α_i`, so the `+√2x` numerator (α=√2) pairs with `e^{−√2}`.)
- **p.894** — `∫ e^x/(x²−2) dx = (1/(2√2))[ e^{√2} ei(x−√2) − e^{−√2} ei(x+√2) ]`.
  (New quadratic constants √2, e^{±√2}; §7. The paper's printed line is OCR-corrupted — the
  algorithm's rule `u_i=x+α_i`, weight `e^{−α_i}`, PF `1/(x²−2)=1/(2√2)[1/(x−√2)−1/(x+√2)]` is
  the source of truth.)
- **Macsyma d8** — `∫ e^x/(x+1)² dx = e^{−1} ei(x+1) − e^x/(x+1)`.
- **Macsyma d11** — `∫ (x²+3)e^x/(x²+3x+2) dx = −7 e^{−2} ei(x+2) + 4 e^{−1} ei(x+1) + e^x`.

---

## 4. Engine C2 — `Erf`/`Erfi` (Cherry 1989 §3, §4.3)

**Integral form:** `Σ d_i e^{−β_i} erf(ũ_i)`, `ũ_i² = f + β_i`, `ũ_i = r_i/s`.
File `src/calculus/cherry_erf.c` (shares `cherry_exp.c` with C1 — same `γ = g e^f` front-end and
the same E1/E2/E3 y-solve run once for both). **Convention:** Cherry `erf(u) = ∫ u' e^{+u²} = Erfi(u)`
(positive exponent). Map to Mathilda: emit `Erfi` directly, or `Erf` via `Erf(x) = Erfi(x)/i`.

### 4.1 Gate — `q` a perfect square

`ũ_i = r_i/s` with `s` the unique monic square root of `q` (`q = s²`, `gcd(r_i,s)=1`). Test with
`PolynomialSqrt[q, x]`; if it returns `$Failed` ⇒ **0 erf terms**, halt the erf branch cleanly
(this is a decision, not a failure).

### 4.2 The ≤2-erf quadratic bound (Thm 3.2)

There are **at most two** constants `β_i ∈ F̄` with `f + β_i = g_i²`; if `β₁,β₂ ∉ F` they are
**quadratic conjugates** over F (Cor 3.3). This bound (`max_terms = 2`) is what makes erf a
decision. Cherry deliberately **avoids** the degree-`3n−1` resultant `Res_x(p+βq, p'+βq')` and uses
the undetermined-constant route below.

### 4.3 R1 / R2 — the completing-square β-finder (§4.3)

From `g = Σ d_i(r_i/s)'` and `(r_i/s)' = ½(s/r_i)(p/q)'`:
`2g / (s·(p/q)') = Σ d_i/r_i`, with denominator `c·h = ∏ r_i` (`c` an unknown constant, `h` the
monic associate of `∏ r_i`).

- **R1 (one erf):** `c·h = r₁ ⇒ p = c²h² − β₁ q`. Solve for constants `R, S` with `p = R h² + S q`
  (unique if they exist) ⇒ `β₁ = −S`; `r₁ = √(p + β₁ q)` via `PolynomialSqrt`; `d₁` from
  `2g/(s(p/q)') = d₁/r₁`. If `R,S` don't exist or the ratio isn't constant, one-erf fails.
- **R2 (two erf):** `c²h² = (p+β₁q)(p+β₂q)`, so `p² = c²h² − [(β₁+β₂)p + β₁β₂ q]q`. Since
  `∂[(β₁+β₂)p + β₁β₂q]q < ∂(h²)`, compute unique `R, S` (`∂R<∂q`, `∂S<∂h²`) with `p² = R h² + S q`;
  then `−S = (β₁+β₂)p + β₁β₂ q`. Determine `U, V` with `−S = U p + V q`; then **`β₁,β₂` are the
  roots of `β² − Uβ + V`** (→ §7 for the algebraic pair). `r_i = √(p + β_i q)` via `PolynomialSqrt`;
  solve `d₁,d₂` from the numerator comparison `k = d₁ r₂ + d₂ r₁` (`k` known).

Implement with `PolynomialQuotient`/`PolynomialRemainder` (for `p = Rh²+Sq`, `p² = Rh²+Sq`),
`CoefficientList`/`SolveAlways` (for the `U,V` and `R,S` undetermined constants), `PolynomialSqrt`,
and `Roots`/`RootReduce` for `β² − Uβ + V`.

### 4.4 C2 test pins

- **Ex 5.2** — `∫ (1/x + 1/x²) e^{1/x²} dx = −½ ei(1/x²) − erf(1/x)`.
  (`q=x²=s²`, `s=x`; R1: `1 = R + Sq ⇒ S=0, β₁=0, r₁=1, d₁=−1`; ei P2: `c=−½`.)
- **Ex 5.4** — `∫ e^{(x⁴+a)/x²} dx = ½ e^{−2√a} erf((x²+√a)/x) + ½ e^{2√a} erf((x²−√a)/x)`,
  over `F = Q(a)`. (`h=x⁴−a`; `p²=Rh²+Sq ⇒ R=1, S=4ax²`; `−4ax²=Up+Vq ⇒ U=0, V=−4a`; `β²−4a` ⇒
  `β=±2√a`; `r₁=x²+√a, r₂=x²−√a`. Weight `e^{−β_i}`: `β₁=2√a` pairs with the `+√a` numerator. §7.)

---

## 5. Engine C3 — `LogIntegral` (Cherry 1986 §4 + §5)

**li-elementary case (iv):** adjoins `∫ u'/log(u) dx = li(u)`. New file
`src/calculus/cherry_sigma_decomp.c` implements the argument generator; the driver wiring lives in
the tower recursion (§6). **Extended-Liouville form (Thm 2.2, eq 2.1):**
`γ = w₀' + Σ c_i(w_i'/w_i) + Σ d_i(u_i'/v_i)`, `v_i' = u_i'/u_i` (so `∫ d_i u_i'/v_i = d_i li(u_i)`),
with **argument constraint `u_i = ∏_j f_j^{α_{ij}}`** (products of powers of irreducible factors).

### 5.1 The Σ-decomposition engine (Thm 4.3 / Thm 4.4)

**Definition:** `Φ ∈ K(x)` has a Σ-decomposition over `Σ = (f₁,…,f_m)` (distinct irreducibles,
none in K) if `Φ = Σ_{i=1}^{n} b_i ∏_j f_j^{α_{ij}}`, `b_i ∈ K`. Neither existence nor uniqueness
holds unrestricted (Ex 4.1: `x` has none over `Σ=(x²+1)`; Ex 4.2: `0` has infinitely many over
`(x, x+1)`). A **restriction map** `g: T → Zᵐ` (with `α_{i1} ∈ T`, `g(α_{i1}) = (α_{i1},…,α_{im})`)
makes it unique (Thm 4.2) and — with a degree bound — decidable (Thm 4.4).

**Inputs:** `Φ`, `Σ` (from denominator factorization), restriction `g`, degree bound `d`
(li needs `d = 1`; dilog / higher SFs reuse at `d ≥ 2` — Cherry's explicit hook, p.7).
**Output:** the finite `{(b_i, α_i-vector)}` giving `u_i = ∏_j f_j^{α_{ij}}`, or a **proof of
non-existence**.

**Steps (Thm 4.4, all reuse existing poly ops):**
1. `α_{11} = multiplicity of f₁ in Φ` (Lemma 4.1: `Φ = f^a p/q`, `a` = multiplicity). If
   `α_{11} ∉ T`, **STOP: no decomposition.** Else `α_{12},…,α_{1m} = g(α_{11})`.
2. `Φ/∏_j f_j^{α_{1j}} = p/q`; `b₁ = (p mod f₁)/(q mod f₁)` (eq 4.3, via `PolynomialRemainder`).
   If `b₁ ∉ K`, **STOP.**
3. Recurse on `Φ₂ = Φ − b₁ ∏_j f_j^{α_{1j}}`; generate `(b₂,α₂),(b₃,α₃),…`.
4. Interpolation polynomials `p_j(α)` (unique, degree `< n`, `p_j(α_{i1}) = α_{ij}`); the *degree*
   of the decomposition is `max_j deg p_j`. Find `α*` beyond which every `p_j` is monotone (if some
   `p_j` has degree 0, divide `Φ` through by that `f_j^{α_j}` to force strict monotonicity).
5. Iterate until `α_{k1} > α*` for `k ≤ d`; if ever `α_{ij} ≠ p_j(α_{i1})`, **STOP: none.**
   **Monotone-tail termination (the NO decision):**
   - *decreasing* sequence in `j`: `α_{nj} = mult(f_j, Φ_{k+1})`; terminate if for some `i>r`
     `α_{nj} < r` (`r = mult(f_j, ·)`).
   - *increasing* sequences: iterate until all `α_{r·}` positive; if `Φ_r` has degree 0, terminate
     in `K[x]` with **failure**; else terminate when for some `i>r`
     `deg(∏_j f_j^{α_{ij}}) > deg(Φ_r)`.

Expose as an internal `cherry_sigma_decompose(Phi, Sigma, g, d)` plus a debuggable builtin
`Integrate\`SigmaDecomposition[Phi, {f1,…}, x]` (degree-1 default) for tests.

### 5.2 Lemma 5.2 divisibility solve

Reconstruct the `v_i` from the denominator of `A`: to find polynomial `d` (`deg d < n = deg f`)
with `f | (p + d q)`: `d ≡ −p y (mod f)` where `x f + y q = 1` (`PolynomialExtendedGCD`). Unique
when it exists; the mechanism behind `u_i` reconstruction in Thm 5.3 cases (a)/(b2).

### 5.3 C3 test pins (`li ↔ LogIntegral`)

- **Ex 1.1 / d1** — `∫ x/log(x)² dx = 2 li(x²) − x²/log(x)`.
- **d2** — `∫ 1/(log(x)+3) dx = e^{−3} li(e³x)` (transcendental-constant rescale: `ū₁ = λ₁u₁`,
  `d₁ → d₁/λ₁`; §7-adjacent but rational-constant).
- **d3** — `∫ x²/log(x+1) dx = li(x³+3x²+3x+1) − 2 li(x²+2x+1) + li(x+1)`.
- **d4** — `∫ (log(x)²+3)/(log(x)²+3log(x)+2) dx = −7 e^{−2} li(e²x) + 4 e^{−1} li(ex) + x`
  (rational-in-log, always integrable — Thm 5.5 family).
- **Ex 5.1** — `∫ x³/log(x²−1) dx = ½ li(x⁴−2x²+1) + ½ li((x²−1)²)` (two-log tower
  `C(x, log(x−1), log(x+1))`; degree-1 Σ-decomp over `Σ=(x−1, x+1)`).
- **Decision-NO — Ex 5.2** — `∫ x²/log(x²−1) dx` is **not** elementary+li: the required
  Σ-decomposition of `(3x⁴−4x²+1)/(4x)` provably does not exist (monotone-tail termination,
  increasing case). Must decline with `rt_dec_nonelem()` → `ElementaryIntegralQ → False`.
- **d12 distinction (documentation pin, not a NO):** `∫ (x²+1)e^x/(x²+x+1) dx` *is* a genuine
  logarithmic/exp integral; Cherry's Macsyma prototype missed it only for lacking algebraically-
  closed constants (deficiency (ii)). With §7 our engine should **find** it — mark it a positive
  pin, contrasting Ex 5.2's true NO.

---

## 6. Full-tower recursion (Cherry 1986 Thm 5.4 / Thm 5.3) — hook into `rt_field_integrate`

The generators must fire *inside* the tower recursion, not only on the outermost integrand.

**Thm 5.4 (induction on `rank(E) = (m_n,…,m_1,1)`):**
- **Base:** `rank = (1)`, `γ ∈ C(x)` — rational Risch (`Integrate\`BronsteinRational`); no li needed.
- **Inductive step — peel the top monomial `θ` of rank `r`:**
  - **Lemma 5.1 monomial split.** Write `γ = Σ_k A_k θᵏ + A₀ + P(θ)/Q(θ)` (Laurent + proper PF in
    θ). `γ` is li-elementary-integrable **iff each term is** — integrate monomial-by-monomial.
    Reuse `rt_hyperexp_case` Laurent split for exp-θ; PF (`Apart`) for log-θ.
  - **Polynomial / PF part `A₀ + P(θ)/Q(θ)`:** *drop to the lower field* — the existing
    `rt_field_integrate` recursion **is** the induction hypothesis. (Footnote 3: `P(θ)/Q(θ)` must
    be elementary over F, licensing `Rothstein` / `intrat_log_part_core` on `∫P/Q`.)
  - **Essential monomial-coefficient part `A_j θʲ`:** route to the reduced problem **Thm 5.3**:
    - *top `θ` logarithmic (Thm 5.4 case a):* the linear-in-θ pieces (eqs 5.10–5.11) form a
      **degree-1 restricted Σ-decomposition** → **C3** (`cherry_sigma_decompose`, `d=1`). Its
      polynomial part reduces to deciding `A₀ − B̄₁(a'/a)` over F (induction hypothesis).
    - *top `θ` exponential (Thm 5.4 case b / Thm 5.3):* `Aθ` → the **C1/C2** generators
      (`Res_x(g₁,p+αq)` for ei, completing-square for erf) via the reduced problem. If a second
      rank-`r` exponential exists (`m_r > 1`), reduce over `D(ψ)` and recurse (induction
      hypothesis). Commensurate exponents (`E^{2u}→t²`, `E^{x/2},E^{x/3}→E^{x/6}`) already
      normalized by `rt_exp_kernelize`.

**Where to wire it.** Register the Cherry hook at the same site as the (already-landed) tower
Hermite / residue-log work — `rt_field_hyperexp_hermite` in `risch_field_integrate.c` — after the
elementary proper-part attempt and before the field-level decline. Per RISCH_REVIEW §5, "landing
the deferred elementary work and landing Cherry are the same refactor"; that elementary work is now
on `main`, so the hook site exists. The driver call is `extended_liouville_solve(remainder, x, T,
top)` on whatever the elementary passes could not absorb.

---

## 7. Algebraically-closed constant layer (`C = C̄`)

Both theorems assume `C = C̄`; arguments carry algebraic constants (√2 in Ex 5.3 / p.894, √a in
Ex 5.4, the quadratic `β² − Uβ + V`, and complex log-roots in d12). This is the single largest
engineering surface.

- **Roots** of `h(α) = Res_x(g₁,p+αq)` (ei) and of `β² − Uβ + V` (erf): obtain via `RootReduce` /
  the `qqbar` layer; keep them in a field of lowest degree (Cherry guarantees erf constants are
  **quadratic** over F, ei constants come from the resultant's factorization).
- **Representation:** carry `√2`, `√a`, `Root[...]` symbolically through `u_i = f + α_i`,
  `ũ_i = r_i/s`, and the prefactors `e^{−α_i}`, `e^{−β_i}`.
- **Verification / normalization:** `flint_algebraic_field_normalize` (the algebraic-tower
  normalizer) to reduce and compare candidate expressions over `Q(√d)` / `Q(ζ_n)` / towers; this is
  the same machinery `RootReduce`/`Together`-over-kernels already use.
- **Conjugate pairing:** algebraic α's/β's come in conjugate pairs producing conjugate SF terms
  (Ex 5.3's two ei, Ex 5.4's two erf); emit both, with conjugate prefactors.

This layer is what lifts d12 from a miss to a hit and makes Ex 5.3/5.4/p.894 exact.

---

## 8. Phasing (each phase ships green; diff-back gate retained throughout)

> **Status (2026-07-16): C3 LogIntegral engine — LANDED (single-log multi-li, first increment).**
> `cherry_li.c` (`rt_cherry_li`) integrates over a single-log tower `theta = Log[w]` (w polynomial)
> emitting `v(x,theta) + Sum_k d_k li(w^k)` — the degree-1 Σ-decomposition over one generator `w`,
> matched IN THE TOWER (`Log[w^k]=k theta` by construction, `D[li(w^k)]=w^(k-1)w'/theta`), so the
> exact rational identity over {theta,x} is the certificate (plain Simplify can't reduce
> `Log[w^k]=k Log[w]`; emission also PowerExpand-verified). Closes Cherry d1 `x/Log[x]^2`, d3
> `x^2/Log[x+1]` (multi-li), Ex 5.1 `x^3/Log[x^2-1]` (via `w=x^2-1` directly), and a broad range.
> Also **transcendental-constant rescaling** (d2 `1/(Log[x]+3) = e^-3 li(e^3 x)`, d4 rational-in-Log
> `= x + 4e^-1 li(ex) - 7e^-2 li(e^2 x)`): a constant root `rho` of the θ-denominator adds a basis
> term `k w'/(θ-rho)` and answer `k e^rho li(e^-rho w)`. Deferred: multi-log towers (reducible `w`
> needing a product Σ-decomposition, `Log[x^3-x]`), the Σ-decomposition NON-existence decision
> (Ex 5.2). Generalises `rt_try_li` (kept as fast path).
>
> **Status (2026-07-16): C2 Erf/Erfi engine — LANDED (completing-square).**
> `cherry_ei.c` now emits error-function terms too (Cherry 1989 §3): when `q = s²` is a perfect
> square, up to two `Erfi` terms `u~_j = r_j/s` (`r_j² = p + β_j q`) are generated by the
> completing-square β-finder (discriminant roots + `b=0`, validated by `PolynomialSqrt`) and solved
> in the SAME linear system as the ei terms — `s = √q` is a polynomial so `u~ = r/s` is rational.
> Closes Ex 5.2 `(1/x+1/x²)E^(1/x²)`, Ex 5.4 `E^((x⁴+a)/x²)` (symbolic parameter), pure `E^(1/x²)`.
> The coefficient solve moved from `SolveAlways` to explicit `Solve[coeff-eqs, {unknowns}]` so a
> symbolic integrand parameter (the `a` in Ex 5.4) stays a parameter. Cherry's `erf = (√π/2)Erfi`,
> so answers carry `√π`. Complex β and non-square `q` decline cleanly.
>
> **Status (2026-07-16): dilogarithm — LANDED (degree-2 Σ-decomposition, partial).**
> `cherry_dilog.c` (`rt_cherry_dilog`) generalises `rt_try_dilog` to the `Log·Log + PolyLog[2]`
> answer form via TOWER matching: dilog args are the interpolants `(x-r_i)/(r_j-r_i)` between the
> rational roots of the linear factors; candidates whose derivative leaves a `Log` (reversed pair
> → `Log[-x]` = iπ shift) are filtered. Closes `Log[x]/(1+x) = Log[x]Log[1+x]+PolyLog[2,-x]`,
> `Log[x]/(1-x)`, `Log[x]/(x^2-1)`, `Log[2x+1]/(x+1)`. Deferred (decline cleanly): transcendental
> root spacings (`Log[2+x]/x`), degree>1 in the log tower (`Log[x]Log[1+x]/x`). Registered as a
> `PolyLog` form behind the `rt_try_dilog` fast path.
>
> **Status (2026-07-16): C4 Fresnel — LANDED.** `Integrate` closes Gaussian-phase trig
> `K Sin[a x^2+b x+c]` / `K Cos[...]` → `FresnelS`/`FresnelC` by completing the square (the trig
> sibling of the Gaussian→`Erf` recognizer), a new Automatic-cascade stage
> `src/calculus/integrate_fresnel.c` after Goursat. Diff-back verified; `Sin[x^2]`, `Cos[x^2]`,
> `Sin[x^2+x+1]`, etc. (Note: basic `Si/Ci/Shi/Chi` — `Sin[x]/x` etc. — already closed pre-C4.)
>
> **Status (2026-07-16): complex-quadratic ei — LANDED (lone conjugate pair).**
> `E^x/(x^2+1)` etc. now close with a complex-conjugate pair of `ExpIntegralEi` (diff-back exact/fast
> over Q(i)). Gate changed to **degree ≤ 2** (`MinimalPolynomial`; real or complex quadratic;
> symbolic admitted; deg ≥ 3 deferred). A complex root is admitted only as the SOLE conjugate pair
> (`q=1`, `g1` irreducible quadratic) — any mixing (P2/reciprocal, extra factor) or a Q(i√3) field
> (d12 `x^2+x+1`) blows up `Together` and declines cleanly. Erf stays real-β.
>
> **Status (2026-07-16): C2 algebraic-constant layer for ei — LANDED (real algebraic).**
> `cherry_ei.c` now admits real algebraic ei-argument constants (Cherry `C=C-bar`): the P1
> resultant roots may be irrational reals (√2, √5, golden ratio) and the coefficients solve over
> the extension they generate — Cherry p.894 `E^x/(x^2-2)`, `E^x/(x^2-x-1)`, etc. all close and
> diff-back. A **numeric reality gate** selects real roots (symbolic `Im`/`Element[.,Reals]` don't
> reduce radicals); **complex** constants (`E^x/(x^2+1)`→±I, d12) are deferred and decline cleanly.
> The y-denominator uses only the **degree-1 (rational-root) factors** of `den(g)·q` — an
> algebraic-root factor stays out of `y`'s denominator, which both matches Cherry's structure and
> sidesteps a shared `Together`/`SolveAlways`-over-`Q(√d)` stack overflow (d≥5). Also: higher-
> multiplicity poles now close; `test_cherry_ei.c` has an extensive stress battery (~60 cases).
>
> **Status (2026-07-16): C1 (ExpIntegralEi, base field, rational constants) — LANDED.**
> `src/calculus/cherry_ei.c` (`rt_cherry_ei`) + `tests/test_cherry_ei.c`. Registered as a second
> `ExpIntegralEi` entry in the `RtSpecialForm` table behind the `rt_try_ei` fast path.
> Deviation from the C0 plan below (deliberate, matching CHERRY_DESIGN.md's "no dead abstraction
> ahead of a consumer" bar, R1/R4): the registry struct was **not** rewritten to the full 4-method
> interface yet — the existing `{name, recognize}` entry point already accommodates the ei engine,
> and the richer `gen_arguments`/`deriv_template`/`answer_term`/`max_terms` seam will be introduced
> when the erf/li generators (C2/C3) actually need it. Solves Ex 5.1, d8, d11, and nonlinear
> exponents; declines cleanly on algebraic constants (deferred to §7). Two now-obsolete negative
> tests in `test_integrate_risch_transcendental.c` updated to positive ei assertions.


- **C0 — substrate + seam (behaviour-preserving).**
  1. Upgrade `RtSpecialForm` to the 4-method struct (§2.2); wrap the current four recognizers as
     entries with their template as the `gen_arguments` fast-path. Byte-identical behaviour.
  2. **Fix `rt_verify_antideriv`** to reject any result containing an unevaluated `Integrate[…,x]`
     head (A4 F5) — precondition for every later phase.
  3. `cherry_driver.c` skeleton (`extended_liouville_solve`) dispatching the registry; wire it in
     place of the raw `rt_special_case` loop. No new outputs yet.
  4. Expose `Integrate\`ExpIntegralEiResultant` beside `RothsteinTragerResultant` (§3.2).
- **C1 — ei.** `cherry_ei.c`: E1/E2/E3 y-solve + P1 resultant generator + P2 q-side term + PF-match
  constant-existence test. Pins: Ex 5.1, d8, d11 (rational constants); Ex 5.3 / p.894 gated on §7.
- **C2 — erf + algebraic constants (§7).** `cherry_erf.c`: `q=s²` gate, R1/R2 completing-square,
  ≤2 bound, `d_i` solve. Land the RootReduce/qqbar constant layer here (Ex 5.4 needs √a). Pins:
  Ex 5.2, 5.4; back-fill Ex 5.3 / p.894 for C1.
- **C3 — li + tower recursion.** `cherry_sigma_decomp.c` (Thm 4.4, `d=1`) + the Thm 5.4 hook (§6).
  Pins: Ex 1.1/d1, d2, d3, d4, Ex 5.1; **decision-NO** Ex 5.2; positive-pin d12.
- **C4 — reductions & dilog.** `SinIntegral/CosIntegral/SinhIntegral/CoshIntegral/FresnelS/FresnelC`
  as ei/erf linear combos (1986 §6 / 1989 App): `li(u)=ei(log u)`, `si(u)=(1/2i)[ei(iu)−ei(−iu)]`,
  `ci(u)=½[ei(iu)+ei(−iu)]`, Fresnel via erf at `√(±πi/2)u`. Dilog via degree-2 Σ-decomposition
  (Cherry's `d≥2` hook), superseding `rt_try_dilog`. Retire the four narrow recognizers once their
  fast-paths are subsumed.

---

## 9. Files

**New**
- `src/calculus/cherry_driver.{c,h}` — `extended_liouville_solve`, registry dispatch, decision
  aggregation.
- `src/calculus/cherry_ei.c` — ei engine (E1/E2/E3, P1/P2). *(May merge with erf as `cherry_exp.c`
  since both share the `γ=g e^f` front-end + y-solve.)*
- `src/calculus/cherry_erf.c` — erf engine (`q=s²` gate, R1/R2).
- `src/calculus/cherry_sigma_decomp.{c,h}` — Σ-decomposition (Thm 4.4).
- Tests: `tests/test_cherry_ei.c`, `tests/test_cherry_erf.c`, `tests/test_cherry_li.c`,
  `tests/test_cherry_sigma.c`.

**Modified**
- `src/calculus/risch_special.{c,h}` — struct upgrade + fast-path wrapping.
- `src/calculus/risch_field_integrate.c` — Thm 5.4 hook + driver call.
- `src/calculus/risch_util.c` — `rt_verify_antideriv` Integrate-head guard (C0).
- `src/calculus/intrat.c` — `Integrate\`ExpIntegralEiResultant` (thin sibling).
- `tests/CMakeLists.txt` — add new sources to `COMMON_SRC` + new `*_tests` targets.
- `Mathilda_spec.md` overview row + `docs/spec/builtins/` (calculus page) + current-week
  `docs/spec/changelog/<Mon>.md`.

**Optional cleanup**
- Correct the three stale references in `CHERRY_DESIGN.md` (file locations, `poly_sqrt` naming,
  `rt_ansatz_solve` → `rt_tower_solve`).

---

## 10. Verification

- **Per engine, foreground only** (memory: no background Mathilda pollers): `make -j`, then
  `cd tests/build && cmake .. && make cherry_ei_tests && ./cherry_ei_tests` (and the analogous erf/li
  targets). Only run the affected `*_tests` binary.
- **Each paper pin** as **exact I/O** plus a **numeric diff-back** (`assert_rm_num`, the existing
  `test_special_functions` pattern in `tests/test_integrate_risch_transcendental.c`). Every emission
  additionally passes `rt_verify_antideriv` by construction.
- **Decision-NO pins** (Ex 5.2 li; the genuinely non-elementary exp cases `E^(x²)`, `E^x/x`,
  `1/Log[x]` from the A4 battery) must return **unevaluated** (`UNKNOWN`) or `ElementaryIntegralQ →
  False` *only* where truly non-elementary — never a wrong antiderivative, never a spurious `False`.
- **Regression:** the full `test_integrate_risch_transcendental` battery unchanged; C0 is
  byte-identical (verify by diffing outputs before/after). `valgrind --leak-check=full ./Mathilda`
  on the new paths (watch the builtin-ownership contract: NULL-out reused sub-`Expr*`s before the
  evaluator frees the wrapper).

---

## 11. Test-pin appendix (exact I/O, both papers)

Convention mapping: `ei ↔ ExpIntegralEi`, `li ↔ LogIntegral`, Cherry `erf ↔ Erfi`
(`Erf(x) = Erfi(x)/i`).

**1986 — Logarithmic integral (li):**
| # | Integral | Result | Note |
|---|---|---|---|
| Ex 1.1 / d1 | `∫ x/log(x)² dx` | `2 li(x²) − x²/log(x)` | |
| d2 | `∫ 1/(log(x)+3) dx` | `e^{−3} li(e³x)` | transcendental-constant rescale |
| d3 | `∫ x²/log(x+1) dx` | `li(x³+3x²+3x+1) − 2 li(x²+2x+1) + li(x+1)` | multi-li |
| d4 | `∫ (log(x)²+3)/(log(x)²+3log(x)+2) dx` | `−7 e^{−2} li(e²x) + 4 e^{−1} li(ex) + x` | Thm 5.5 family |
| Ex 5.1 | `∫ x³/log(x²−1) dx` | `½ li(x⁴−2x²+1) + ½ li((x²−1)²)` | two-log tower |
| **Ex 5.2** | `∫ x²/log(x²−1) dx` | **NON-ELEMENTARY** (decline) | monotone-tail NO |

**1989 — Rational exponential (ei, erf):**
| # | Integral | Result | Note |
|---|---|---|---|
| Ex 5.1 | `∫ e^{1/x} dx` | `x e^{1/x} − ei(1/x)` | |
| Ex 5.2 | `∫ (1/x + 1/x²) e^{1/x²} dx` | `−½ ei(1/x²) − erf(1/x)` | ei + erf |
| Ex 5.3 | `∫ (x⁴−1)/(x⁵+x) e^{(x²+1)/x} dx` | `½ e^{−√2} ei((x²+√2x+1)/x) + ½ e^{√2} ei((x²−√2x+1)/x)` | √2; §7 |
| Ex 5.4 | `∫ e^{(x⁴+a)/x²} dx` | `½ e^{−2√a} erf((x²+√a)/x) + ½ e^{2√a} erf((x²−√a)/x)` | √a; §7 |
| p.894 | `∫ e^x/(x²−2) dx` | `(1/(2√2))[ e^{√2} ei(x−√2) − e^{−√2} ei(x+√2) ]` | √2; §7 (OCR-reconstructed) |
| d8 | `∫ e^x/(x+1)² dx` | `e^{−1} ei(x+1) − e^x/(x+1)` | |
| d11 | `∫ (x²+3)e^x/(x²+3x+2) dx` | `−7 e^{−2} ei(x+2) + 4 e^{−1} ei(x+1) + e^x` | |
| d12 | `∫ (x²+1)e^x/(x²+x+1) dx` | ei terms with complex-conjugate args | positive pin (needs §7) |

*OCR caveat:* the exponent signs on Ex 5.3's result line and the entire p.894 formula render
ambiguously in the scanned paper. The authoritative rules — weight `e^{−α_i}` with `u_i = f + α_i`
(ei), weight `e^{−β_i}` with `ũ_i² = f + β_i` (erf) — are the source of truth for those signs and
must be re-derived, not copied, when pinning.

---

## 12. One-line summary

Cherry is not a new engine: it is `rt_tower_solve` given three finite argument-generators — the
`Res_x(g₁,p+αq)` resultant (ei), the completing-square β-finder (erf), and the Σ-decomposition
(li) — behind the `RtSpecialForm` registry, hooked into the `rt_field_integrate` tower recursion
(Thm 5.4) and backed by the algebraic-constant layer, on the canonical-rep + structure-oracle +
LRT + tower-Hermite foundation already on `main`. Ship C0 (seam + verifier guard) first, then one
paper example at a time: ei → erf → li.
