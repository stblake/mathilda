# Definite & Line Integration — Review of Deficiencies

Review date: **2026-07-07**. Scope: `Integrate[f, {x, a, b}]` (Newton–Leibniz /
fundamental-theorem path, `src/calculus/integrate_newton_leibniz.c`) and complex
line/contour integration (`src/calculus/integrate_line.c`), driven against a
46-case corpus of classical and stress examples through the live REPL.

Legend: **FIXED** = resolved in this review (with regression test); **OPEN** =
known gap, not yet addressed. Every OPEN item is a *safe* failure — the integral
is left **unevaluated**, never returned wrong.

---

## 1. Bugs found and FIXED

### 1.1 Risch–Norman stack-overflow crash — *critical* — FIXED
- **Symptom:** hard `SIGSEGV` killing the whole process (not just the integral).
- **Repro:**
  - `Integrate[x*2.71828^(-x), x]`
  - `N[Integrate[x*Exp[-x], {x, 0, Infinity}]]` (N turns `Exp` into `2.71828^…`)
  - `N[Integrate[x^2*Exp[-x], {x, 0, Infinity}]]`
- **Cause:** an inexact-base exponential integrand drove `split_factor`
  (`intrischnorman.c`) into unbounded recursion; content extraction never
  reduced the indeterminate set, so the recursion overflowed the C stack (lldb:
  guard-page hit in a `split_factor → … → split_factor` chain of identical state).
- **Fix:** recursion-depth guard `PMINT_MAX_SPLIT_DEPTH` in `split_factor_rec`;
  on overflow it returns failure and the integral is left unevaluated.
- **Test:** `test_intrischnorman.c :: test_phase1_inexact_base_no_stack_overflow`.

### 1.2 Reversed limits (`a > b`) returned a wrong value — *high* — FIXED
- **Symptom:** a divergent integral silently returned a spurious finite (and
  spuriously complex) value.
- **Repro:**
  - `Integrate[1/z, {z, 1, -1}]` → `I Pi` (correct: divergent — path crosses the
    pole at 0)
  - `Integrate[1/x, {x, 2, -2}]` → `I Pi`
- **Cause:** `nl_classify` assumed `a < b`, so an interior pole with reversed
  bounds was classified *exterior*; the driver then returned `F(b) − F(a)` =
  `Log[-1] − Log[1]`.
- **Fix:** made interior-pole classification orientation-independent (a pole
  strictly between the bounds is detected for both `a < b` and `a > b`).
  Divergence is now caught (`Integrate::idiv`, unevaluated); non-pole reversed
  integrals still negate correctly (`Integrate[x^2, {x, 1, 0}]` → `-1/3`).
- **Test:** `test_integrate_newton_leibniz.c :: test_divergent`,
  `test_reversed_limits`, and reversed cases in `test_detector`.

### 1.3 Line integrator missed on-path singularities of quadratic denominators — FIXED
- **Symptom:** a contour segment passing through a pole returned unevaluated but
  **without** the `Integrate::idiv` divergence warning.
- **Repro:** `Integrate[1/(z^2+1), {z, I-1, I+1}]` (segment passes through the
  pole at `I`).
- **Cause:** the parametrised denominator is a polynomial in the real parameter
  `t` with *complex* coefficients; `Solve[den == 0, t, Reals]` drops its real
  roots (`Solve[…, Reals]` returns `{}` for complex-coefficient equations).
- **Fix:** `line_segment_singularities` now solves over ℂ and keeps the real
  roots via its existing `im ≈ 0` numeric filter.
- **Test:** `test_integrate_line.c :: test_detection`, `test_divergence`.

---

## 2. OPEN deficiencies (documented gaps, safe failures)

All items below currently leave the integral **unevaluated**. None returns a
wrong answer. They are ordered roughly by how classical / high-value the missing
result is. The root causes live in *upstream* engines (Limit, indefinite
integration), not in the definite-integration drivers themselves.

### 2.1 Improper integrals blocked by the `Limit` engine
The Newton–Leibniz path is correct, but the boundary limit it needs is not
resolved, so the whole integral bails.

| Integral | Expected | Blocking limit (returns unevaluated) |
|---|---|---|
| `Integrate[x*Exp[-x], {x, 0, Infinity}]` | `1` | `Limit[(-1-x) E^(-x), x→∞]` = 0 |
| `Integrate[Log[x], {x, 0, 1}]` | `-1` | `Limit[x Log[x], x→0]` = 0 |
| `Integrate[x*Log[x], {x, 0, 1}]` | `-1/4` | `Limit[x Log[x], x→0]` = 0 |
| `Integrate[Exp[-x^2], {x, -∞, ∞}]` | `Sqrt[Pi]` | `Limit[Erf[x], x→∞]` = 1 (also 2.3 below) |
| `Integrate[x^2 Exp[-x^2], {x, -∞, ∞}]` | `Sqrt[Pi]/2` | same |

- **Root cause:** `Limit` does not evaluate polynomial × decaying-exponential at
  ∞, `x Log[x]` at 0, or `Erf` at ∞.
- **Fix direction:** strengthen `src/calculus/limit.c` (dominant-term / L'Hôpital
  handling for `poly·exp` and `x^a Log[x]`; known limits of `Erf`, `Erfc`).

### 2.2 `Integrate[Sin[x]^2, {x, 0, Pi}]` fails on an unsimplified antiderivative
- **Symptom:** returns unevaluated; the *numeric-bound* form additionally emits a
  spurious `Integrate::idiv` (the integrand `Sin[x]^2` is bounded and obviously
  convergent). `Integrate[Sin[x]^2, {x, 0, 2Pi}]` → `Pi` works.
- **Cause:** the **indefinite** integral returns a monstrous rational-trig
  antiderivative instead of `x/2 − Sin[2x]/4`, and its denominator
  `6 − 2Cos^4 + 12Cos^2 + 16Cos + 2Sin^4` **vanishes at x = π**. FTC then hits a
  0/0 form at the upper limit that the Limit engine cannot resolve. `[0, 2π]`
  survives only because plain substitution avoids the interior split.
- **Fix direction:** power-reduction / cleaner antiderivative for
  `Sin^2`/`Cos^2` (and trig-power integrands generally) in the indefinite
  cascade, so the antiderivative has no spurious singularities.

### 2.3 Gaussian antiderivative not produced
- **Symptom:** `Integrate[Exp[-x^2], x]` returns unevaluated (no `Erf`).
- **Cause:** the indefinite integrator emits no `Erf`-based antiderivative for
  `Exp[-x^2]` (and relatives). Compounds 2.1 for the definite Gaussian.
- **Fix direction:** recognise `Exp[quadratic]` → `Erf`/`Erfi` antiderivative.

### 2.4 `Abs[x]` and other piecewise integrands
- **Symptom:** `Integrate[Abs[x], {x, -1, 1}]` (expected `1`) unevaluated.
- **Cause:** `Integrate[Abs[x], x]` unsupported.
- **Fix direction:** piecewise handling of `Abs`, `Sign`, `UnitStep`, `Floor`,
  splitting the domain at the breakpoints.

### 2.5 Entire integrand needing a special-function antiderivative (line path)
- **Symptom:** `Integrate[Sin[z]/z, {z, 1, I, -1, -I, 1}]` (closed loop, entire
  integrand, expected `0`) unevaluated.
- **Cause:** the line integrator needs an antiderivative in `x` first; `Sin[z]/z`
  antidifferentiates to `SinIntegral[z]` (`Si`), which is not produced, so it
  bails. (Correctly bails — the answer is not returned wrong.)
- **Fix direction:** `SinIntegral`/`CosIntegral`/`ExpIntegral` antiderivatives; or
  a residue/Cauchy shortcut for closed contours of meromorphic integrands.

---

## 3. What is solid (verified, no issues)

- **Line/contour integrator:** every stress case correct — path independence,
  winding number 2 (`→ 4 Pi I`), unenclosed poles → 0, double/triple poles → 0,
  integrable endpoint branch points, and residues (`1/(z^2+1)` around `+I` → `Pi`).
  On-path and endpoint singularities are reported divergent.
- **Core FTC:** polynomials, elementary trig, rational integrands, infinite-bound
  `1/(1+x^2)`, full-period Weierstrass rationals (`1/(2+Cos[x])` etc.), and
  integrable endpoint singularities (`1/Sqrt[x]`, `1/Sqrt[1-x^2]`).
- Reversed limits without an interior pole negate correctly.

---

## 4. Corpus / harness

The review harness (`/tmp/intreview/harness.py`) runs each case through the REPL
JSON protocol and verifies via a numeric difference check
`N[Chop[(integral) - (expected), 10^-8]] == 0`, classifying divergent cases by an
`idiv` warning + unevaluated result. After the fixes: **35 / 46 pass**; the 11
non-passes are exactly the OPEN items in §2 (all unevaluated, none wrong).
