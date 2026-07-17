# CHERRY_BLOCKERS.md — What stands between us and "complete Cherry"

**Date:** 2026-07-17
**Companion:** `CHERRY_PLAN.md` (the executable spec), `CHERRY_DESIGN.md` (architecture).
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

The registry is still the 2-field `RtSpecialForm {name, recognize}` (`risch_special.c:252`); the
4-method interface, `cherry_driver.c`/`extended_liouville_solve`, `Integrate\`ExpIntegralEiResultant`,
`cherry_sigma_decomp.c`, the `rt_verify_antideriv` Integrate-guard, and the Thm 5.4 tower hook **do
not exist**.

---

## A. True blockers (stuck pending hard math or missing infrastructure)

### A1 — Complex algebraically-closed constants (`C = C̄` over `Q(i)` and beyond)
- **Gated pins:** d12 `∫ (x²+1)e^x/(x²+x+1) dx` (needs `Q(i√3)`); the general complex case of
  Ex 5.3 / p.894 beyond the single-quadratic-pair fast path; any erf with complex β.
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

### B2 — Full-tower recursion (Cherry 1986 Thm 5.4 / Thm 5.3) — hook into `rt_field_integrate`
- **Gated pins:** any Cherry-integrable integrand whose special-function structure only appears
  *after* peeling an outer monomial (nested towers, e.g. ei/li terms inside `E^(…)` over
  `C(x, Log[x])`). Today the engines never see peeled monomial coefficients.
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

### B3 — C0 substrate + seam (the plan's *prerequisite*, skipped)
- **Not a capability gap** but the architectural debt that makes B1/B2 clean and A2 soundly-gated.
  Four sub-items, all behaviour-preserving (`CHERRY_PLAN.md:509–516`):
  1. Upgrade `RtSpecialForm {name, recognize}` → the 4-method struct (`applicable`,
     `gen_arguments`, `deriv_template`, `answer_term`, `max_terms`); wrap the 7 current recognizers
     as fast-path `gen_arguments`. (`risch_special.c:252`.)
  2. ~~`rt_verify_antideriv` Integrate-head guard~~ — already resolved at the engine-gate level, see
     **A3 below** (✅). Nothing to do here.
  3. `cherry_driver.c` / `extended_liouville_solve` skeleton dispatching the registry, replacing
     the raw `rt_special_case` loop. No new outputs.
  4. `Integrate\`ExpIntegralEiResultant[g1,p,q,α,x] = Resultant[g1, p+α q, x]` beside
     `RothsteinTragerResultant` in `intrat.c` (`CHERRY_PLAN.md:226–230`).
- **Effort:** medium; low risk (byte-identical behaviour — verify by diffing REPL outputs before/
  after). **Recommended to do *before* B1/B2** so those land against a real seam.

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
| C-i | `∫ Log[3+2x]/x dx` | non-monic linear dilog kernel | `cherry_dilog.c` builds interpolants from **monic** `x−r`; a non-monic `w=3+2x` isn't normalized | extract rational content: `Log[a x+b] = Log[a] + Log[x+b/a]`, feed monic factor + fold `Log[a]` constant into the `Log·Log` term |
| C-ii | `∫ Log[x]Log[1+x]/x dx` | degree > 1 in the log tower | dilog gate requires `deg ≤ 1` in every `t_i` (`cherry_dilog.c:194`) | genuinely non-elementary in the *dilog* class as written — verify it is `Integrate::nonelem` (the test suite already asserts this decline) or needs a higher-`d` Σ-decomposition (Cherry's `d≥2` hook, ties to B1) |
| C-iii | li output form | `x/Log[x]²` emits `ExpIntegralEi[2 Log[x]]` where the pin reads `2 li(x²)` | `li(u)=ExpIntegralEi[Log u]` — **mathematically equal**, cosmetic only | optional: normalize ei-of-log → `LogIntegral` in the emission for output parity with the paper |

---

## Recommended sequence

```
A3  (verifier guard)                         ✅ already resolved at the engine gates
B1  (general Σ-decomposition, Thm 4.4)        ✅ LANDED 2026-07-17 (Integrate`SigmaDecomposition)
 └─ A2 (li decision)                          ✅ LANDED 2026-07-17 (Integrate`LiElementaryQ)  → 1986 li DECISION complete
B3  (C0 seam: registry + driver + resultant)  — medium, prerequisite for B2
B2  (Thm 5.4 tower recursion)                 — large   ┐ non-proportional product decomps + nested towers
A1  (complex constant layer over FLINT NF)    — large   ┘ 1989 ei/erf complete
C-i / C-ii / C-iii                            — bounded, land opportunistically
```

**Bottom line (updated 2026-07-17):** the *narrow* Cherry is done and green, and the **1986 li
decision property is now complete** — **B1** (general Σ-decomposition, `Integrate`SigmaDecomposition`)
and **A2** (li decision, `Integrate`LiElementaryQ`) landed. What remains for *complete* Cherry:
**B2** (Thm 5.4 tower recursion — the positive path for non-proportional product decompositions and
nested towers), **B3** (the C0 seam — its clean prerequisite), and **A1** (the complex `C=C̄`
constant layer, the only one with a real *infrastructure* dependency: FLINT number-field arithmetic
replacing generic `Together`/`SolveAlways`). None is fundamentally impossible.
