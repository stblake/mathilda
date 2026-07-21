# CHERRY_BLOCKERS.md — What stands between us and "complete Cherry"

**Date:** 2026-07-17
**Companion:** `CHERRY_PLAN.md` (the executable spec), `CHERRY_DESIGN.md` (architecture).
**Successor layer:** [`KNOWLES_DESIGN.md`](KNOWLES_DESIGN.md) (2026-07-21) — Knowles' extension of
Cherry from transcendental-*elementary* to transcendental-*Liouvillian* integrands (integrand may
itself contain `li`/`erf`/`Ei`), in terms of `erf` (1992/93) and `li` (1986 §2). It reuses this
subsystem wholesale; its B2 deep-tower peel and the `RT_PRIM` Liouvillian-primitive tower generator
are the shared substrate.
**Purpose:** an honest, verified ledger of everything remaining before Cherry's algorithm
(1986 li + 1989 ei/erf) is *complete* — i.e. has the **decision property**, fires **inside the
tower recursion**, and handles **algebraically-closed constants** in full. Each entry names the
gated test pins, the root cause, its dependencies, and a concrete resolution path.

---

## 0. Where we actually are (verified 2026-07-16 / 2026-07-17)

Three narrow-but-correct engines shipped **ahead of** the C0 substrate the plan sequenced first
(`CHERRY_PLAN.md:497–506` acknowledges this reversal). All three test suites are green
(`cherry_ei_tests`, `cherry_li_tests`, `cherry_dilog_tests`), and the following pins are **verified
working in the REPL**:

| Engine | Working pins (exact I/O confirmed) |
|---|---|
| ei (`cherry_ei.c`) | Ex 5.1 `E^(1/x)`, Ex 5.3 `√2`, p.894 `E^x/(x²−2)`, d8, d11, real-algebraic (√5, golden), lone complex-conjugate pair `E^x/(x²+1)` |
| erf (in `cherry_ei.c`) | Ex 5.2 `(1/x+1/x²)E^(1/x²)`, Ex 5.4 `E^((x⁴+a)/x²)` (symbolic `a`), `E^(1/x²)` |
| li (`cherry_li.c`) | d1 `x/Log[x]²`, d2 `1/(Log[x]+3)`, d3 `x²/Log[x+1]` (multi-li), d4 (rational-in-log), Ex 5.1 `x³/Log[x²−1]` |
| dilog (`cherry_dilog.c`) | `Log[x]/(1+x)`, `Log[x]/(1−x)`, `Log[x]/(x²−1)`, transcendental spacing `Log[2+x]/x` |
| Fresnel (`integrate_fresnel.c`) | `Sin[x²]`, `Cos[x²+x+1]`, … (C4) |

**Dispatch reality:** all three engines are reached from a *single* site —
`rt_special_case` (`risch_special.c:268`) called from `rt_integrate_core`
(`integrate_risch_transcendental.c:316`) — so they fire **only on the outermost integrand**, and
every decline is a bare `return NULL` → `RT_DEC_UNKNOWN` → unevaluated `Integrate`. **No Cherry
engine can currently prove non-existence.**

**Update (2026-07-17):** the C0 seam now exists — `RtSpecialForm` carries a top-monomial
applicability mask + `rt_special_case_routed`, `cherry_driver.c`/`extended_liouville_solve` is the
single dispatch point, and `Integrate\`ExpIntegralEiResultant` is registered (B3). The flat
multi-term ei generator (`rt_cherry_exp_multiterm`, Thm 5.4 case b) and the lone complex-conjugate
constant layer (A1) also landed. Not built: the deep-tower peel of Thm 5.4 for depth-≥2 nested
towers, and the general complex/degree-≥3 constant case (see B2 / A1 below).

---

## A. True blockers (stuck pending hard math or missing infrastructure)

### A1 — Complex algebraically-closed constants (`C = C̄` over `Q(i)` and beyond) — ✅ LANDED for the lone pair (2026-07-17)
- **Delivered (lone conjugate pair over `Q(i)` and many `Q(i√d)`):** d12
  `∫ (x²+1)e^x/(x²+x+1) dx` closes with a complex-conjugate `ExpIntegralEi` pair
  (`Q(i√3)`), as do `E^x/(x²+x+1)`, `E^x/(x²+3)`, `E^x/(x²+2x+5)`,
  `(2x+1)e^x/(x²+x+3)` (`Q(i√11)`), … — all diff-back exact.
- **Stress-test finding + FIX (2026-07-17):** the lone conjugate pair now closes over
  **ANY `Q(i√d)`**, uniformly. A stress sweep found that the direct-`Solve` path only
  cracked the fields Mathilda's `Solve` happened to handle (`E^x/(x²+2x+3)` over
  `Q(i√2)` declined). Root cause: `Together` fails to cancel the two complex-linear
  factors back against the real quadratic, over-determining the system. New fallback
  **`rt_cherry_ei_conjpair_nf`** (`cherry_ei.c`) solves it **over Q** by the `{1, chs}`
  number-field basis method (write `a = center ± chs`, `chs² = disc`; carry `c0,c1` and
  `y` with rational unknowns; reduce mod `chs² − disc`; split into the `{1, chs}` basis;
  `Solve` over Q). Fires only when the direct solve failed, so prior closures stay
  byte-identical. Closes `E^x/(x²+2x+3)`, `E^x/(x²+2x+7)`, `E^x/(x²−2x+3)`,
  `(3x+1)E^x/(x²+2x+3)`, … (`center` needs `Simplify`, not `Together`).
- **Constant exponent offset — FIXED (2026-07-17):** `E^(c + h(x))` with `c` x-free
  (e.g. `E^(1/x+2) = E²·E^(1/x)`) now closes — the constant is factored out before the
  P2 recognition (which the inflated `deg(p)` had defeated). Fires only when the
  polynomial part of the exponent is a nonzero constant.
- **Root-cause correction (important):** the coefficient `Solve` over `Q(i√d)` is NOT
  intrinsically the blocker — the *isolated* small system solves fine natively. The
  failure was **system size**: the generous real-case `Y`-degree bound (`Ny≈6`)
  inflated the ansatz so `Solve` over the algebraic field could not reduce the larger
  mixed system. Fix (`cherry_ei.c`): for an admitted complex candidate — where the
  only pole is the irreducible quadratic `g1`, so `y` is a small polynomial — tighten
  `Ny` to the polynomial-part degree + `deg(sden)` + 1. Native `Solve` then handles
  `Q(i√d)`; the exact diff-back gate keeps it sound (a too-small bound only declines).
  So the full FLINT number-field linear solve was **not needed** for the lone-pair
  case; it remains the path for the *still-deferred* general case below.
- **Still deferred (decline cleanly):** a complex pair MIXED with a `P2`/reciprocal
  term (`E^(1/x)/(x²+1)`) or an extra factor, degree-`≥3` constants (`E^x/(x³−2)`),
  and any erf with complex β — these are the cases where the ansatz cannot be kept
  small and generic `Together`/`Solve` over the tower still blows up; the FLINT
  number-field arithmetic route (below) is the resolution for them.
- **Symptom:** `E^x/(x²+x+1)` returns unevaluated. The lone-conjugate-pair path only survives when
  `q=1` and `g₁` is a *single* irreducible quadratic over `Q(i)`; a `Q(i√3)` field (or any mixing
  with a P2/reciprocal term or an extra factor) makes `Together`/`SolveAlways` over the algebraic
  tower blow up (stack overflow / non-termination), so the engine gates it out
  (`cherry_ei.c` real-β / degree-≤2 gate, `CHERRY_PLAN.md:479–495`).
- **Root cause:** this is genuinely an **infrastructure** blocker, not just unbuilt Cherry logic —
  symbolic `Together`/`SolveAlways` over `Q(i√d)` / `Q(ζ_n)` towers is too slow and stack-hungry
  for the coefficient solve. `CHERRY_PLAN.md:411–428` (§7) calls this "the single largest
  engineering surface."
- **Dependencies:** the FLINT algebraic layer — `flint_algebraic_field_normalize`, `RootReduce`/
  `qqbar`, cyclotomic `Together` (see memories *FLINT extension engine*, *RootReduce qqbar*,
  *Cyclotomic extension support*). The coefficient solve must run **inside the number field**, not
  via generic `Together`.
- **Resolution path:** route the ei/erf coefficient linear system through the FLINT number-field
  arithmetic (represent α/β as `qqbar`/`Root[...]`, solve over `Q(α)` with FLINT rather than
  `SolveAlways` over kernels). Pair conjugate roots and emit conjugate SF terms with conjugate
  prefactors (`e^{−α_i}`). Verify with `flint_algebraic_field_normalize` instead of `Simplify`.
- **Effort:** large. **Risk:** medium — correct-by-construction if the field arithmetic is exact,
  but easy to time out. Prove non-regression against the currently-green real-algebraic battery.

### A2 — Σ-decomposition NON-existence decision (the li decision property) — ✅ LANDED (2026-07-17)
- **Delivered as `Integrate`LiElementaryQ[f, x]`** (`cherry_sigma_decomp.c` + `rt_cherry_li_nonelem`
  in `cherry_li.c`): `True` when the Cherry li engine exhibits an elementary+`LogIntegral`
  antiderivative, `False` when the Σ-decomposition non-existence certificate fires
  (Ex 5.2 `x²/Log[x²−1]`, `x²/Log[x³−x]`), unevaluated outside the single-log scope.
- **Important correction to the original plan.** The plan targeted `ElementaryIntegralQ → False` for
  Ex 5.2, but that is **already** satisfied by the existing Bronstein field decision (anything
  needing `li` is not *elementary*, so `Risch`ElementaryIntegralQ` returns `False` for Ex 5.1 AND
  Ex 5.2 already). The genuinely-new decision is **li-elementarity** — Ex 5.1 `x³/Log[x²−1]` is
  li-elementary (`True`) but not elementary (`False`); Ex 5.2 is neither. That distinction is what
  `Integrate`LiElementaryQ` surfaces, and it is the faithful Cherry-1986 decision property.
- **Soundness:** the `False` verdict fires only for the pure essential form `A/Log[w]` (w squarefree)
  via the faithful Thm 4.4 termination on `Φ = A/w'` — never on a routing decline. `rt_dec_nonelem`
  was **not** needed (the ElementaryIntegralQ target was already met by existing machinery).

---

## B. Large unbuilt components (tractable, but substantial — sequence these)

### B1 — General Σ-decomposition engine (Cherry 1986 Thm 4.4) — ✅ LANDED (2026-07-17)
- **Delivered as `cherry_sigma_decompose` + `Integrate`SigmaDecomposition[Phi, {f1,…}, x]`**
  (`src/calculus/cherry_sigma_decomp.c`): the faithful Thm 4.4 decision — multiplicity extraction,
  `b_i = (p mod f₁)/(q mod f₁)`, recursion, the cross-factor consistency check, and the
  increasing-case degree-overshoot terminator that **proves** non-existence (`$Failed`). Verified on
  the paper's exact numbers (Ex 5.1 `x²/2 → {{1/2,{0,0}},{1/2,{1,1}}}`; Ex 5.2 `x/2 → $Failed`).
- **Scope note:** the restriction is the degree-1 **all-equal** `g(r)=(r,…,r)` (all li single-log
  decompositions collapse to this — the args are powers of `rad(w)`). Genuinely *non-proportional*
  product decompositions arise only inside the tower recursion (**B2**), which is still deferred;
  the engine's structure carries the general shape but non-all-equal restriction maps are the B2
  extension.
- **Positive li integration** stays on the existing tower-solve in `cherry_li.c` (correct and
  diff-back-gated); B1 powers the **decision** (via A2) and the debuggable builtin, not a new
  positive path. `∫ x²/log(x³−x)` is (correctly) proven **not** li-elementary rather than integrated.

### B2 — Full-tower recursion (Cherry 1986 Thm 5.4 / Thm 5.3) — ✅ FLAT case b LANDED (2026-07-17); deep-tower peel deferred
- **Delivered (Thm 5.4 case b, flat exponential level):** `rt_cherry_exp_multiterm`
  (`cherry_ei.c`, registered behind `rt_cherry_ei`) integrates a single-kernel `E^w`
  integrand with SEVERAL commensurate Laurent terms `Σ p_i E^(i w)` — which the
  single-shape `rt_cherry_ei` cannot peel (its cofactor keeps a residual exponential)
  — by Laurent-splitting in `t = E^w` and integrating each `p_i E^(i w)` with the ei
  engine, summing (diff-back verified). Closes `∫(E^x+E^(2x))/(x-1) = E ei(x-1) +
  E² ei(2x-2)`, `E^x(E^x+1)/((x-1)(x-2))`, `(E^x+E^(2x))/(x²-2)`, etc.
- **Also landed (B3):** the C0 seam that made this a one-registration change — the
  `RtSpecialForm` top-mask + `extended_liouville_solve` driver + `ExpIntegralEiResultant`.
- **Still deferred:** the DEEP tower peel — a Cherry structure exposed only after
  peeling an outer monomial of a depth-≥2 tower. Investigation (2026-07-17) found the
  reachable nested cases (e.g. `(Log[Log[x]]+1)/(x Log[Log[x]]) = Log[x]+li(Log[x])`)
  are ALREADY closed by the existing primitive-polynomial recursion, and a hook at the
  `rt_field_integrate` decline hit verification fragility (an ei answer's derivative
  reintroduces a shifted exponential `E^(iw+α)` that `rt_tower_deriv` does not collapse
  to a tower monomial). The flat case-b engine covers the demonstrable multi-exponential
  gap; a genuinely-new depth-≥2 pin that the existing recursion misses was not found.
- **Root cause:** engines dispatch only at `integrate_risch_transcendental.c:316` (outermost). The
  §6 induction (generators firing on `A_j θʲ` inside the Lemma 5.1 monomial split) is unwired.
- **Dependencies:** cleanest atop **C0** (the 4-method registry + `extended_liouville_solve`
  driver), so the hook is one `extended_liouville_solve(remainder, x, T, top)` call rather than
  three ad-hoc dispatches. The hook site (`rt_field_hyperexp_hermite`, `risch_field_integrate.c`)
  already exists per RISCH_REVIEW §5.
- **Resolution path:** register the Cherry hook after the elementary proper-part attempt and before
  the field-level decline (`CHERRY_PLAN.md:402–407`); route log-top → C3 (B1), exp-top → C1/C2.
- **Effort:** large. **Risk:** medium (touches the hot Risch recursion — guard with the full
  `test_integrate_risch_transcendental` battery for non-regression).

### B3 — C0 substrate + seam — ✅ LANDED (2026-07-17)
- **Delivered (behaviour-preserving, byte-identical across the whole pin battery):**
  1. `RtSpecialForm` gains a **top-monomial applicability mask** (`RT_SF_TOP_EXP` /
     `RT_SF_TOP_LOG` / `RT_SF_TOP_ANY`) + a routed dispatch `rt_special_case_routed`
     (`risch_special.c`). *Design note:* the full 4-method struct (`gen_arguments`/
     `deriv_template`/`answer_term`/`max_terms`) was **deliberately not** added — the
     Cherry engines are self-contained and re-derive their own structure from kernel-form
     `f`, so those method pointers would have no consumer (the project's "no dead
     abstraction ahead of a consumer" bar). The consumed seam is the mask + driver.
  2. ~~`rt_verify_antideriv` Integrate-head guard~~ — already resolved at the engine gates (A3).
  3. `src/calculus/cherry_driver.c` / `extended_liouville_solve(f, x, top_mask)` — the
     single Cherry dispatch point the outermost integrator and any future tower hook call.
  4. `Integrate\`ExpIntegralEiResultant[g1,p,q,a,x] = Resultant[g1, p+a q, x]` in
     `intrat.c`, registered + docstring + `tests/test_rt_resultant.c` coverage.

---

## A3 — Integrate-head guard in the Cherry emission gates — ✅ RESOLVED (verified 2026-07-17)
- **Status:** **done** where the 2026-07-16 design correction (`CHERRY_PLAN.md:74–83`) says it
  belongs. All three engines guard every emission with `rt_free_of_head(Q,"Integrate")`:
  `cherry_ei.c:506` (`… && rt_verify_antideriv(Q,f,x)`), `cherry_li.c:259`, `cherry_dilog.c:339`.
  A candidate carrying a residual `Integrate[…,x]` is dropped before diff-back, so it cannot be
  rubber-stamped by `D[Integrate[f,x],x]→f` (A4 Finding 5 closed for the Cherry path).
- **By design, the shared `rt_verify_antideriv` (`risch_util.c:317–351`) does NOT carry the guard** —
  a blanket guard there would regress the *intended* partial-log-part emission
  (`risch_singleext.c:548–563`, A4 Finding 4 = not a defect). Leave it as-is.
- **Carry-forward:** B1/B2 assemble multi-term answers — keep the same
  `rt_free_of_head(Q,"Integrate")` guard on any new emission path they add.

---

## C. Narrow deferred gaps (bounded, independent of A/B)

| ID | Pin | Gap | Root cause | Path |
|---|---|---|---|---|
| C-i | `∫ Log[3+2x]/x dx` | non-monic linear dilog kernel — **cherry_dilog-internal only; the full Integrate cascade closes it, so not user-facing** (verified 2026-07-17) | `cherry_dilog.c` builds interpolants from **monic** `x−r`; a non-monic `w=3+2x` isn't normalized, AND rational-root kernels (`Log[x+1/2]/(x+1)`) decline even after monic-normalizing — two stacked issues | extract rational content `Log[a x+b] = Log[a] + Log[x+b/a]` AND add rational-root interpolant support; deprioritized since the cascade already returns a correct answer |
| C-ii | `∫ Log[x]Log[1+x]/x dx` | degree > 1 in the log tower | dilog gate requires `deg ≤ 1` in every `t_i` (`cherry_dilog.c:194`) | genuinely non-elementary in the *dilog* class as written — verify it is `Integrate::nonelem` (the test suite already asserts this decline) or needs a higher-`d` Σ-decomposition (Cherry's `d≥2` hook, ties to B1) |
| C-iii | li output form | `x/Log[x]²` emits `ExpIntegralEi[2 Log[x]]` where the pin reads `2 li(x²)` | `li(u)=ExpIntegralEi[Log u]` — **mathematically equal**, cosmetic only | optional: normalize ei-of-log → `LogIntegral` in the emission for output parity with the paper |

---

## Recommended sequence

```
A3  (verifier guard)                         ✅ resolved at the engine gates
B1  (general Σ-decomposition, Thm 4.4)        ✅ LANDED 2026-07-17 (Integrate`SigmaDecomposition)
 └─ A2 (li decision)                          ✅ LANDED 2026-07-17 (Integrate`LiElementaryQ)  → 1986 li DECISION complete
B3  (C0 seam: mask + driver + resultant)      ✅ LANDED 2026-07-17 (byte-identical)
B2  (Thm 5.4 case b, flat multi-term ei)      ✅ LANDED 2026-07-17 (rt_cherry_exp_multiterm); deep-tower peel deferred
A1  (complex C=C̄, lone conjugate pair)        ✅ LANDED 2026-07-17 (d12 + Q(i√d) family); mixed/deg≥3 deferred
C-i / C-ii / C-iii                            — bounded, land opportunistically
```

**Bottom line (updated 2026-07-17):** the *narrow* Cherry is done and green; the **1986 li
decision property** (B1 + A2) is complete; and the **C0 seam (B3)**, the **flat multi-term ei
generator (Thm 5.4 case b, B2)**, and the **lone complex-conjugate-pair constant layer (A1,
d12 + `Q(i√d)`)** all landed 2026-07-17. What remains for *fully* complete Cherry is bounded and
well-understood: the **deep-tower peel** of Thm 5.4 (a Cherry structure exposed only after peeling
an outer monomial of a depth-≥2 tower — no reachable pin found that the existing recursion misses),
and the **general complex/higher-degree constant case** (a complex pair mixed with a P2 term, or a
degree-≥3 constant tower) — the one item with a real *infrastructure* dependency: an exact FLINT
number-field linear solve replacing generic `Together`/`Solve` (the lone-pair case was closed
without it, by keeping the ansatz small). None is fundamentally impossible.
