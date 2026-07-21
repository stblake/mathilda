# CHERRY_GAUSSIAN_PLAN.md — closing `∫ P(x) e^{-c x²} dx` (symbolic constant + numeric-c hang)

**Date:** 2026-07-21
**Trigger:** `Integrate[x^4 Exp[-c x^2] (a x^4 - b), x]` (symbolic `c`) → currently
`Integrate::nonelem`. Mathematica returns a rational·`E^{-c x²}` part + `Sqrt[Pi] Erf[Sqrt[c] x]`.
This is inside Cherry's 1989 error-function extension (`γ = g e^f`, `f` quadratic → `erf`); the
engine already handles the *shape* (polynomial × Gaussian, ≤2 completing-square erf terms) and
trips on exactly two narrow issues, both in `src/calculus/cherry_ei.c`.

**Companions:** `CHERRY_PLAN.md` §4 (C2 erf), `CHERRY_BLOCKERS.md`. Root causes verified by two
independent code traces + scratch-patch prototypes (diff-back = 0, patches reverted).

---

## Root cause (verified)

Current behaviour (isolated REPL runs):

| Integrand | exponent coeff | now |
|---|---|---|
| `x^n Exp[-x^2]` | `c=1` | ✅ solves |
| `x^4 (a x^4 - b) Exp[-2 x^2]` | `c=2`, symbolic poly | ✅ solves |
| `x^4 (a x^4 - b) Exp[-c x^2]` | **symbolic `c`** | ❌ declines ← target |
| `x^4 Exp[-2 x^2]` (pure numeric) | `c=2` | ⚠️ **hangs** |

- **D1 — decline (symbolic `c`).** `PolynomialSqrt[-c x^2, x] → $Failed`: `ps_is_numeric`
  (`src/poly/facpoly.c:87`) admits only Integer/Real/Rational/Complex content, so the x-free
  **symbolic** factor `-c` falls to the "bare irreducible" branch (`facpoly.c:145`). The erf
  β-finder gets 0 candidates (`cherry_ei.c:212-213`); with 0 ei candidates the engine declines
  (`cherry_ei.c:528`). Nothing else gates on symbolic `c` — `real_admissible` admits undecided
  values, and the coefficient solve uses explicit `Solve` (not `SolveAlways`) so `c` stays a param.
- **D2 — hang (numeric `c≠1`).** The erf argument becomes `Sqrt[-c]·x` = `I Sqrt[2] x` for `c=2`;
  the matching identity carries `2/Sqrt[Pi]` (`cherry_ei.c:621-623`). `Together` at
  `cherry_ei.c:634` → generic multivariate subresultant GCD (`rat.c:793`) over ~9 pseudo-generators
  (`x, a0..a4, k0, Sqrt[2], Sqrt[Pi]`) → super-exponential blowup (true GCD = 1). c=1 (no
  irrational) and symbolic c/a,b (different generator ordering) dodge it. Hang, not a wrong answer.

---

## Plan

### Fix A — `PolynomialSqrt` carries x-free (symbolic) content  → unblocks the target
- [ ] **A1.** `src/poly/facpoly.c` `builtin_polynomialsqrt`: capture `psx = (arg_count==2) ? args[1]
  : NULL`. In the factor loop (`facpoly.c:124-148`), treat a factor `t` as **content** (fold into
  `konst`, carried through `Sqrt`) when `ps_is_numeric(t) || (psx && free_of(t, psx))`; same for a
  `Power[base,e]` with x-free `base`. Add a small local `free_of(e,var)` helper. **Restrict to the
  2-arg form** so bare `PolynomialSqrt[p]` is untouched.
  Soundness unchanged: the exact `Expand[s^2 - p] == 0` certificate (`facpoly.c:162-170`) still gates.
- [ ] **A2.** `tests/test_polynomialsqrt.c`: add `PolynomialSqrt[-c x^2, x]` (→ `Sqrt[-c] x`, s²==p),
  `c (x+1)^2`; assert bare 1-arg `PolynomialSqrt[c x^2]` still `$Failed`.

### Fix B — keep `2/Sqrt[Pi]` out of the erf linear system  → kills the D2 hang, hardens D1
- [ ] **B1.** `src/calculus/cherry_ei.c`: fold the transcendental prefactor out of the erf unknown.
  At `:621-623` drop the `twopi = 2 Pi^(-1/2)` factor → `rhs[2+m+j] = K_j · (r_j/s)'` (system stays
  over `Q(algebraic)`, no `Sqrt[Pi]` denominator). At answer assembly `:667-673` multiply the Erfi
  term by `Sqrt[Pi]/2`: `ans = K_j · (Sqrt[Pi]/2) · E^(-beta_j) · Erfi[r_j/s]`. Algebraically
  identical; `Together` at `:634` no longer sees a transcendental-constant denominator. ei terms
  unaffected.
- [ ] **B2.** (Optional, only if a pin still hangs) `rat.c` cancel early-out: when FLINT GCD declines
  and the denominator is a unit in every polynomial generator (constant), return uncancelled —
  mirrors `rat_has_dependent_power_generators` (`rat.c:786`).

### Verification
- [ ] **V1.** `make -j`; `cd tests/build && cmake .. && make cherry_ei_tests polynomialsqrt_tests`.
- [ ] **V2.** New pins in `tests/test_cherry_ei.c`, each **exact I/O + numeric diff-back** with a
  wall-clock ceiling (must not hang):
  - target: `x^4 (a x^4 - b) Exp[-c x^2]`, `x^4 Exp[-c x^2]` (symbolic c)
  - numeric c≠1: `x^2 Exp[-2 x^2]`, `x^4 Exp[-2 x^2]`, `x^2 Exp[-3 x^2]`, `x^2 Exp[-(1/2) x^2]`
  - regression: `x^2 Exp[-x^2]`, `x^4 Exp[-x^2]`, Ex 5.2 `(1/x+1/x²)E^(1/x²)`, Ex 5.4 `E^((x⁴+a)/x²)`.
- [ ] **V3.** Full `cherry_ei_tests` / `cherry_li_tests` / `cherry_dilog_tests` /
  `test_integrate_risch_transcendental` unchanged.
- [ ] **V4.** `valgrind --leak-check=full` on the new paths (NULL-out reused sub-`Expr*` per the
  builtin ownership contract).

### Docs
- [ ] **D.** `docs/spec/builtins/` (calculus + PolynomialSqrt) + current-week
  `docs/spec/changelog/<Mon>.md`; add the closure row to `CHERRY_BLOCKERS.md` §C.

---

## Notes
- Fix A alone makes the **target** (symbolic c) solve — scratch-patch verified, diff-back 0.
  Fix B is required for the numeric-`c≠1` hang and hardens the symbolic path. Do both.
- Output form: `Erfi[Sqrt[-c] x]` (Cherry `erf = (√π/2)Erfi`) equals Mathematica's
  `Erf[Sqrt[c] x]` (`Erfi[i z]=i Erf[z]`). Cosmetic `Erfi→Erf` normalization is an optional
  follow-up, not a correctness item.
- No new caps/magic bounds; soundness rests on the diff-back gate.
