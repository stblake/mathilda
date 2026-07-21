# Design: Knowles' Error-Function & Logarithmic-Integral Integration of Transcendental Liouvillian Functions

**Date:** 2026-07-21
**Targets:** P. H. Knowles —
- *Integration of Liouvillian Functions with Special Functions* (SYMSAC '86, Proc. pp. 179–184)
  → the two decision procedures (elementary + `li`; elementary + `erf`) for **transcendental
  Liouvillian** integrands.
- *Integration of a Class of Transcendental Liouvillian Functions with Error-Functions, Part I*
  (J. Symb. Comput. 13, 1992, 525–543) → the rigorous `erf` procedure for the hardest subroutine
  (a quasiquadratic top exponential monomial, Theorem 7.1).
- *… Part II* (J. Symb. Comput. 16, 1993, 227–241) → the general `erf` procedure (Prop 3.4,
  Lemma 3.5 / Theorem 3.6), relaxing the `*`-reduced precondition to a post-hoc branch.

All three PDFs are in the repo root.

**Companions (the layer *below* this one):** [`CHERRY_DESIGN.md`](CHERRY_DESIGN.md),
[`CHERRY_PLAN.md`](CHERRY_PLAN.md), [`CHERRY_BLOCKERS.md`](CHERRY_BLOCKERS.md). Cherry integrates
transcendental *elementary* integrands in terms of `li`/`ei`/`erf`; **Knowles is the strict
extension to transcendental *Liouvillian* integrands** (the integrand may itself already contain
`li`, `erf`, `Ei`, nested logs/exps). Knowles reuses Cherry's engine wholesale and only enlarges
(i) the integrand class and (ii) the argument generators/predicates.

**Subject:** `src/calculus/` Risch/Cherry modules (`risch_tower.c`, `risch_field_integrate.c`,
`risch_special.c`, `cherry_driver.c`, `cherry_sigma_decomp.c`, and new `knowles_*.c`).

> **This document is a design/spec, not an implementation** — the same discipline as
> `CHERRY_DESIGN.md`. The executable tracker is `tasks/todo.md`; the full planning record is
> `/Users/user/.claude/plans/let-s-perform-an-extensive-iterative-bentley.md`.

---

## 0. The one idea (identical to Cherry's)

Both Cherry and Knowles prove an **extended Liouville theorem** of the same shape (Singer–Saunders–
Caviness 1985, via Cherry's corollaries):

```
∫ g dx  =  v  +  Σ_i k_i · SF_i(a_i)                                   (†)
```

- `v` — the **elementary part** (an element of the differential field).
- `SF_i(a_i)` — a **special-function term** from a fixed menu (`erf`, `li`, `ei`), with a
  **generated argument** `a_i` and an **undetermined constant** `k_i`.

Differentiating (†) gives the **matching identity** the algorithm solves as a **linear system over
the constants**:

```
g  =  D_tower[v]  +  Σ_i k_i · T_i(a_i)                                (‡)
```

where `T_i(a_i) = d/dx[SF_i(a_i)]/k_i` is a rational function of the tower variables once `a_i` is
fixed. **The whole problem is: generate the finite set of candidate `a_i`, then solve (‡).** If no
constants work, no antiderivative of this form exists — the **decision property**.

This is exactly `rt_tower_solve` (`risch_field_integrate.c`) — build an undetermined-coefficient
ansatz, form `D_tower[Q] − F`, clear denominators, `SolveAlways` over the constants. **Knowles =
`rt_tower_solve` with special-function basis terms whose arguments come from new generators, over a
tower that now admits non-elementary primitives.**

**What Knowles adds over Cherry** (the only differences):
1. the integrand class widens *elementary → Liouvillian* — the tower must admit `li`/`erf`/`Ei`
   generators (§2.1);
2. the argument generators run *over the tower*, not just `C(x)` (§2.2);
3. the `erf` case adds the **quasiquadratic** predicate and the `*`-reduced post-hoc branch (§2.3);
4. the recursion fires generators on *inner* peeled coefficients (Prop 3.4), not just the outermost
   integrand (§2.4).

---

## 1. Mathematical review (algorithmic viewpoint)

### 1.1 Lineage & problem class
Characteristic zero, algebraically closed constants `C`, `x' = 1`.
- Risch 1969: transcendental *elementary* integrand → *elementary* answer.
- Cherry 1983/85/86: transcendental *elementary* integrand → answer may use `erf` / `li` / `ei`.
- **Knowles 1986/92/93: transcendental *Liouvillian* integrand** (rationals + algebraics +
  **primitives (antiderivatives)** + exponentials) → answer may use elementary + `erf` (Part I/II)
  or elementary + `li` (1986 §2). Strictly contains Cherry: integrates over fields that are *not
  reduced* (e.g. `∫ exp(½ log log x − 1/log x)/(x log²x) dx = 2 erf(1/√log x)`).

### 1.2 Structural foundation (SSC 1985 / Cherry Lemma 5.1)
If `g ∈ F` (Liouvillian) has an **erf-elementary** antiderivative, then for `c_i, d_i ∈ C`,
`w_i ∈ F`, and `u_i, v_i` algebraic over `F` with `v_i'/v_i = −(u_i²)'`:

```
g = w_0' + Σ c_i (w_i'/w_i) + Σ d_i (u_i' v_i)              (erf terms are the u_i'v_i)
```

and the `li` analogue with `d_i (u_i'/v_i)`, `v_i = log(u_i)` (Theorem 2.1, 1986). Once the finite
candidate `{u_i, v_i}` set is known, solving for `c_i, d_i` is linear over `C` (Mack/SSC Thm A1) —
`rt_tower_solve`.

### 1.3 Σ-decompositions (the candidate engine; Cherry §4 / Part I §3–4)
`Φ = Σ_i b_i ∏_j f_j^{α_ij}` (`b_i ∈ K`, `f_j` distinct irreducible, `α_ij ∈ Z`). Lemmas 3.1/3.2 +
Cor 4.1–4.4 give a *terminating* computation (or a non-existence proof) with exponent bounds. This
pins "which monomials can appear inside an erf/li argument" to a finite search. Implemented for the
degree-1 all-equal case as `cherry_sigma_decompose`; Knowles needs the non-all-equal / product
generalization.

### 1.4 The erf procedure (Part I core; Part II general)
For integrand `g·θ_n` with top monomial `θ_n`:
1. **Rothstein–Trager** on the erf argument: `−u_i² = s + β_i`, `β_i = Σ_L r_ij θ_j + Σ_E r_ij a_j +
   ν_i` (`2r_ij ∈ Z`, `ν_i ∈ C`) — the resultant/LRT idea reused to find erf arguments.
2. **Perfect-square gate (Lemma 6.2):** a finite computable set of squarefree `Π_i` with
   `Π_i·(s+β_i) = R_i²`; only `u_i` making this a perfect square can occur (finiteness).
3. **Case analysis (Cases 1–2.3.2)** by `rank(s)` vs top rank `r`, each subcase reducing to a
   divisor/perfect-square test (Lemma 5.3–5.5), a Σ-decomposition (Cor 4.4), or a linear system
   (SSC A1).
4. **Part II general:** Lemma 2.1/2.2 + Cor 2.3 peel the top monomial by partial fractions,
   reducing `g`'s erf-terms to those of a lower coefficient `Ā_0`; **Prop 3.4** inducts on rank;
   **Theorem 3.6** makes `*`-reduced a post-hoc branch (return {not-`*`-reduced | no-antiderivative
   | ∫}).

### 1.5 Predicates the erf case adds
- **Quasiquadratic exponential** (Part I §6 / Part II §2 / 1986 §3): `θ = exp(a)` with
  `a = (ℓ/2)θ* + (aθ²+bθ+c)/(θ+f)` and side-conditions (a)–(f)/(h) on ranks and the square-ness of
  `c − bf + af²`. The shape that *forces* an `erf` into the answer.
- **`*`-reduced field:** no relevant exponential monomial is quasiquadratic. Cherry needs *reduced*;
  Knowles needs the weaker `*`-reduced (Part I), and Part II Theorem 3.6 needs it only as a
  fallback branch — so we implement a **per-monomial quasiquadratic detector**, not a whole-field
  normalizer.

### 1.6 The li procedure (1986 §2, Theorem 2.3)
`li-elementary` extension adds `θ_i = li(u)` (`θ_i' = u'/log u`). Decide li-elementary
integrability: enumerate candidate `li(u)` (Σ-decomposition), then SSC A1 solves the constants.
Example `∫ dx/(log x · log li(x)) = li(li(x))`. No quasiquadratic/`*`-reduced machinery — the
cleaner procedure.

### 1.7 Regression pins (verbatim from the papers)
- `∫ exp(½·loglog x − 1/log x)/(x log²x) dx = 2 erf(1/√log x)`   (Part I Ex 8.1; 1986 Ex 3.2)
- `∫ exp(−x² − erf²x) dx = erf(erf(x))`                          (Part II Ex 4.1; 1986 Ex 3.1)
- `∫ x·exp(−x² − erf²x) dx` → **no** erf-elementary antiderivative   (Part II Ex 4.2; decision NO)
- `∫ erf(x)·exp(−x² − erf²x) dx = −½ exp(−erf²x)`   (Part II Ex 4.3; elementary over the tower)
- `∫ [2 e^{−x²} erf(x) − 3 e^{−1/x²}/x²] dx = erf²(x) + 3 erf(1/x)`   (Part II Ex 4.4)
- `∫ dx/(log x · log li(x)) = li(li(x))`                          (1986 §2)
- Convention: `d/dx ERF = (2/√π)e^{−x²}` (Mathilda's classical `Erf`) vs Knowles `erf(u)=∫e^{−u²}du`;
  differ by `2/√π` — a rescale in the answer constructor (`erf(u) = (√π/2)·Erf(u)`).

---

## 2. Fit with Mathilda (the mapping is near-1:1 with Cherry)

### 2.0 Substrate that Knowles reuses as-is

| Knowles need | Existing facility | File |
|---|---|---|
| SSC/Mack constant solve (Lemma 5.1 / Thm 2.1) | `rt_tower_solve` (SF-aware; `D[]` knows `Erf`/`LogIntegral`/`ExpIntegralEi`) | `risch_field_integrate.c` |
| Σ-decomposition (Part I §3–4; 1986 §2) | `cherry_sigma_decompose`, `Integrate\`SigmaDecomposition` | `cherry_sigma_decomp.c` |
| Rothstein–Trager on erf arg; li residues | `rt_frac_lrt`, `rt_field_lrt_logpart`, `Integrate\`RothsteinTragerResultant` | `risch_singleext.c`, `intrat.c` |
| perfect-square gate, `√(p+βq)` (Lemma 6.2, 5.3–5.5) | `PolynomialSqrt`, completing-square β-finder | `poly/facpoly.c`, `cherry_ei.c` |
| peel top monomial by partial fractions (Part II Lem 2.1/2.2) | `Risch\`CanonicalRepresentation`, `parfrac.c` | `risch_canonical.c` |
| differential tower with per-generator `Dcoef` | `RtTower` | `risch_tower.c` |
| rank-descending recursion (Prop 3.4 / Thm 5.4) | tower depth index in `rt_field_integrate` | `risch_field_integrate.c` |
| diff-back soundness gate | `rt_verify_antideriv` | `risch_util.c` |
| answer heads | `SYM_Erf/Erfi/Erfc`, `SYM_ExpIntegralEi`, `SYM_LogIntegral` | `sym_names.c` |
| algebraic constants (`C=C̄`) | `RootReduce`/`qqbar`, FLINT number-field | poly/FLINT |

### 2.1 New machinery — the gap (verified absent)
1. **Liouvillian-primitive tower generator (load-bearing).** `RtTower` admits only `{LOG,EXP,TAN}`.
   Knowles' integrand may contain `li(u)`/`erf(u)`/`Ei(u)`, each a **primitive monomial** `θ'=a'`:
   `li: θ'=u'/log u`, `Ei: θ'=e^u u'/u`, `erf: θ'=(2/√π)e^{−u²}u'`. Add an `RT_PRIM` kind carrying an
   explicit derivation-coefficient `Expr`. **Subtlety:** an `RT_PRIM`'s `Dcoef` references *lower*
   tower kernels (erf's `e^{−u²}` is itself an `RT_EXP` generator that must sit below the `erf`), so
   the build/order/structure-check must seed those kernels first and express `Dcoef` in tower
   variables — the same triangularity discipline the existing log/exp ordering already enforces.
2. **Deep-tower recursion firing generators on inner coefficients (Part II Lem 2.1/2.2, Prop 3.4).**
   `extended_liouville_solve` fires only on the outermost integrand (`CHERRY_BLOCKERS.md` B2). Wire
   `extended_liouville_solve(remainder, x, top−1)` at the `rt_field_integrate` decline.
3. **Quasiquadratic detector + `*`-reduced post-hoc branch** (`knowles_quasiquadratic.c`).
4. **erf/li argument enumeration lifted to the tower** (Part I RT `−u²=s+β` + Lemma 6.2 gate +
   Cases 1–2.3.2; Σ-decomposition non-all-equal for li).
5. **Explicit rank** — optional; start from the depth index, promote only if peeling needs the
   `F_i/F_{i-1}` stratification explicitly.

### 2.2 Where it plugs in (no cascade surgery)
New `RT_SPECIAL_FORMS[]` entries in `risch_special.c` (`{ "Erf", knowles_erf_liouvillian,
RT_SF_TOP_EXP }` after Cherry's `erf`; `{ "LogIntegral", knowles_li_liouvillian, RT_SF_TOP_LOG }`
after Cherry's `li`), reached through `extended_liouville_solve` — firing on the outer integrand and
(once the K0 hook lands) every peeled coefficient, with zero edits to `builtin_integrate`. New
files: `src/calculus/knowles_erf.c`, `knowles_quasiquadratic.c`, optional `knowles_li.c`.

### 2.3 Language: C (decided)
Beside `cherry_*.c`, not `.m`. Knowles is a **decision procedure** (Part I Cases 1–2.3.2, Prop 3.4
induction), not an identity/table — the project's two-tier rule puts decision procedures in C. It
binds to C-internal APIs (`rt_tower_solve`, `RtTower`, `rt_frac_lrt`, `cherry_sigma_decompose`,
`extended_liouville_solve`) not exposed as builtins; K0 is a struct/enum change; the dispatch seam
is a C function-pointer registry. Algebra-heavy inner steps are already builtins, called through the
evaluator via `rt_eval_call`/`rt_template` — so the C reads as orchestration of Wolfram-level ops.

---

## 3. Phased plan (critical path K0 → K2; each ships green)

**Target: K2 (erf-Liouvillian).** K0 is a hard prerequisite; K1 (li) is an optional warm-up.
Soundness is free (every emission diff-back-gated); completeness (and the decision property) is the
work.

- **K0 — Liouvillian-primitive tower substrate.** `RT_PRIM` kind + `Dcoef`; `rt_tower_build_min`
  admits `Erf/Erfi/LogIntegral/ExpIntegralEi` kernels (seeding their lower exp/log kernels);
  `rt_tower_deriv`/`rt_dt_i` handle `RT_PRIM` (`θ'=Dcoef`); deep-tower recursion hook. Test:
  `test_knowles_tower.c` (round-trip `D`), Cherry battery byte-identical.
- **K1 — li-Liouvillian (1986 §2; OPTIONAL warm-up).** `knowles_li.c` over the K0 tower; non-all-
  equal Σ-decomposition; register after Cherry's `li`. Pin `∫dx/(log x·log li(x))=li(li(x))`;
  extend `Integrate\`LiElementaryQ`.
- **K2 — erf-Liouvillian (Part I + Part II Thm 3.6). ★ TARGET.**
  - K2a `knowles_quasiquadratic.c`: quasiquadratic detector + `Π(s+β)=R²` perfect-square gate.
  - K2b erf-argument enumerator (RT `−u²=s+β` + Cases 1–2.3.2, over the tower).
  - K2c `knowles_erf.c`: Part II peel + Prop 3.4 induction + `rt_tower_solve`; Thm 3.6 post-hoc
    branch; `erf(u)=(√π/2)Erf(u)` convention shim; registry entry.
  - Pins: §1.7 erf battery + the negative decision.
- **K3 — algebraic-constant completion (`C=C̄`).** Shared with Cherry A1 (FLINT number-field solve);
  off the critical path for the rational-constant pins.

**Effort/risk:** K0 medium/low · K1 medium/low-med · K2 large/med (Part I is the 19-page hardest
subroutine; diff-back keeps emissions sound — only the "NO" verdict rests on the structural bounds)
· K3 large/med (infra-gated).

---

## 4. Testing strategy (extensive; stress tests → future docs)

Two tiers, string-I/O through the real evaluator (`eval_is`/`assert_ei` pattern from
`test_cherry_ei.c`), standalone `*_tests` binaries added to `tests/CMakeLists.txt` `COMMON_SRC`:
1. **Unit tests** (component contracts): `test_knowles_tower.c` (K0 round-trip),
   `test_knowles_quasiquadratic.c` (detector + gate), `test_knowles_erf_args.c` (candidate sets).
2. **Stress tests → tutorial corpus** (`test_knowles_erf.c`, `test_knowles_li.c`): every §1.7 pin +
   systematically harder generated cases (deeper nesting, mixed `li`/`erf` towers, varied
   constants), each asserting {result contains the SF head} **and** {diff-back `D[∫]−f ≡ 0`}; plus
   a **decision battery** of negatives asserting sound declines. Authored as worked examples — the
   seed corpus for a `docs/spec/` tutorial ("Integration in terms of error functions").
3. **Memory:** valgrind-clean on the new binaries (builtin ownership; NULL-out-before-free).

---

## 5. Non-goals / risks
- **`C = C̄`** — the arguments carry algebraic constants (`√2`, quadratic `β`); the complex/higher-
  degree constant case shares Cherry A1's infra dependency (FLINT number-field solve). Not on the
  critical path for the rational-constant pins.
- **Soundness free, completeness is the work.** The diff-back gate keeps every *emitted* answer
  correct; a "NO" verdict is only a sound *decision* if the structural bounds (≤ 2 erf terms,
  Σ-decomposition termination, perfect-square gate) are implemented faithfully.
- **Keep generators finite/bounded.** Every new form must declare its bound and constant-existence
  test; an unbounded candidate set breaks the decision property.

---

## 6. One-line summary
Knowles is not a rewrite: it is Cherry's `rt_tower_solve` linear-solve over a tower that now admits
non-elementary primitives (`RT_PRIM`), fed by tower-lifted argument generators plus the
quasiquadratic detector, behind the same `RT_SPECIAL_FORMS[]` seam — extending Cherry from
transcendental-*elementary* to transcendental-*Liouvillian* integrands. Do K0 (the substrate) first;
then the erf procedure (K2), one paper example at a time.
