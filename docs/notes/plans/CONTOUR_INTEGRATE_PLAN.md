# Definite Integration by Contour Integration — Implementation Plan

## Context

Mathilda evaluates real definite integrals two ways today, both
antiderivative-driven:

- **Newton–Leibniz (FTC)** — find an antiderivative, take one-sided limits at the
  bounds (`src/calculus/integrate_newton_leibniz.c`).
- **Complex line integration** — parametrize a *given* contour γ(t) and integrate
  (`src/calculus/integrate_line.c`).

Both are unreliable, or produce messy limit-of-arctan output, for the large class
of **improper integrals on ℝ** and **periodic trigonometric integrals** that
classical complex analysis dispatches cleanly via the **residue theorem** — e.g.

- ∫₋∞^∞ 1/(1+x⁴) dx = π/√2
- ∫₋∞^∞ cos x/(1+x²) dx = π/e
- ∫₀^{2π} 1/(2+cos θ) dθ = 2π/√3
- ∫₋∞^∞ sin x / x dx = π (principal value / Jordan half-residue)

This plan adds a residue-theorem method that **recognizes** these forms, sums
residues over the poles enclosed by the appropriate contour, and returns the
closed form. It is **correct-by-construction and gated by an `NIntegrate` numeric
crosscheck**: if the symbolic residue sum does not close to a verified scalar, the
method returns `NULL` and the existing Newton–Leibniz path takes over. No existing
integration path is weakened.

**Scope (confirmed with the user):** the four families below **including
principal-value / indented-contour handling** for poles on the real axis. The
public `Residue[]` primitive it builds on is already shipped (see below).

## Prerequisite — SHIPPED

`Residue[f, {z, z0}]` (`src/calculus/residue.c`, commit `250b8cc`) is the
primitive this method leans on: it returns the coefficient of `(z-z0)^-1` in the
Laurent expansion via `Series`, adaptively raising the expansion order for
unknown-numerator poles, and returns unevaluated at branch points (`den != 1`).
The contour method calls `Residue` (or the internal `residue_compute`, to be
exported from `residue.h`) once per enclosed pole. See
[`project_residue_implementation`] in memory and
[`docs/spec/builtins/calculus.md#residue`].

## Confirmed infrastructure to reuse

- **Dispatch:** `integrate_definite(Expr* res)` in `src/calculus/integrate.c`
  routes a single real spec `{x,a,b}` to `integrate_newton_leibniz_try(f,x,a,b,
  method)`; complex/polyline specs go to the line-integral parametrizer. The
  method enum, `method_from_string()`, `definite_parse_method()`, and a switch
  dispatch are the extension points. `integrate_init()` calls sub-module
  `*_init()`s.
- **Bounds classification:** `nl_numeric(e,&val)` → 0=symbolic, 1=finite, 2=+Inf,
  3=−Inf. `Infinity` = `SYM_Infinity`; `-Infinity` = `DirectedInfinity[-1]`.
- **Numeric crosscheck:** `nl_crosscheck_ok` (in integrate_newton_leibniz.c)
  builds `NIntegrate[f,{x,a,b}]`, mutes arithmetic warnings, and compares to a
  `1e-3*(1+|v|)` tolerance. `NIntegrate` already handles infinite ranges.
- **Pole finding:** `Solve[Q==0, x]` over Complexes with radical options;
  `solve_cubic_radical` / `solve_quartic_radical` in `src/poly/solvepoly.c`
  (quartic path fires over Complexes — exactly what we want); numeric `NRoots`
  returns `Complex[re,im]` roots for the UHP/unit-disk sign tests.
- **Closure:** `RootReduce` (FLINT qqbar, `src/rootreduce.c`) canonicalises
  radical / root-of-unity residue sums; `Simplify`, `Together`, `Numerator`,
  `Denominator`, `builtin_im` (`Im[]`), `is_rational_in` (already in integrate.c).
- **`RootSum` is NOT a closure vehicle** — it is Lagrange-collapse-only
  (`src/root.c`), cannot sum arbitrary UHP residues or restrict to a half-plane.
  Never emit it as a definite answer.

## Approach — one method, four recognizers

Each recognizer is a narrow, conjunctive gate. Poles are found with the solver,
classified by `N[]`-evaluated sign of `Im`/`Abs`, and residues summed via the
`Residue` primitive.

**Family A — rational on (−∞,∞).** Gates: `f = P/Q` rational in `x`, `deg Q ≥
deg P + 2`, `Q` has no real root, bounds exactly `(−∞,+∞)`. Value `= 2πi · Σ Res`
over poles with `Im > 0` (upper-half-plane semicircle). Exploit the real-
coefficient conjugate-pair structure: sum only the UHP half.

**Family B — Fourier / Jordan on (−∞,∞).** Gates: `f = R(x)·K`,
`K ∈ {Exp[I a x], Cos[a x], Sin[a x]}`, `R` rational with deg-drop ≥ 1, `R`'s
denominator has no real root, `a` a nonzero real of **known sign**. Compute
`J = 2πi·Σ_UHP Res[R·Exp[I a x]]` for `a>0` (close UHP; for `a<0` close LHP with
`−2πi·Σ_LHP`). Then `∫R cos = Re[J]`, `∫R sin = Im[J]`. Do **not** split
`Cos → (e^{iax}+e^{-iax})/2` — only one exponential decays under Jordan.
Symbolic `a` of unknown sign → NULL.

**Family C — trig over a full period** `(0,2π)` or `(−π,π)`. Recognize by
transformation: `TrigExpand`, substitute `Cos[x]→(z+1/z)/2`,
`Sin[x]→(z−1/z)/(2I)`, multiply by `dz/(I z)`; if the result is rational in `z`
(`is_rational_in`), the family fires. Value `= 2πi·Σ Res` over poles with
`|z| < 1`. This is the **strongest** family — finite-range crosscheck always
available.

**Family PV — real-axis poles (principal value / Jordan half-residue).** Simple
poles on ℝ (or on `|z|=1` for C) contribute a **half-residue** `πi·Res`; combine
with the full residues of off-axis poles. Enables ∫ sin x / x = π and PV forms.
Requires all on-axis poles to be **simple** (higher-order on-axis → NULL). Needs
a modified crosscheck (`NIntegrate` PV / limit of symmetric truncation).

**Half-line `[0,∞)` via even symmetry.** If `PossibleZeroQ[f(x) − f(−x)]`, then
`∫₀^∞ = ½ · ∫₋∞^∞` (e.g. ∫₀^∞ 1/(1+x⁴) = π/(2√2)). Small, high-value add-on.

## Closure pipeline (per family)

1. Prefer **radical roots**: `deg Q ≤ 4`, or `Factor`/`FactorList` first so pieces
   are ≤ quartic — qqbar closes radicals reliably.
2. Sum `Σ Residue[...]` over the enclosed poles; multiply by `2 π I`
   (or `π I` half-residues for on-axis simple poles).
3. `Together → Simplify → RootReduce → Simplify`. Chop a residual imaginary part
   below tol (a real integral must be real — a surviving non-negligible `Im`
   means the recognizer mis-fired → NULL).
4. **Numeric-verify** vs `NIntegrate` at 1e-3 rel. For infinite-range A/B where
   `NIntegrate` cannot validate, fall back to a large-finite-window check with
   loosened tol; if that also fails, NULL.
5. If the result did not collapse to a closed scalar (still contains `Root` /
   `Power[...,1/n]` after `RootReduce`) → **NULL** (fall through to N–L). Never
   emit an unsimplified `Root`/`RootSum` as the definite value.

## Ordering vs Newton–Leibniz

In `integrate_definite`, the residue method runs **before** Newton–Leibniz for a
single real spec under `Automatic` (or explicit `Method->"Residue"`), because for
improper/periodic forms it yields cleaner, verified closed forms. On NULL it falls
through to N–L unchanged. When the user pins `Method->"NewtonLeibniz"`, the residue
method is skipped. Engage only for a single real spec with matching bounds; leave
the complex-spec and multi-spec (iterated) paths untouched.

## File decomposition

**New `src/calculus/integrate_residue.{c,h}`:**

```
/* Master entry, called from integrate_definite before NL. Borrowed args. */
Expr* integrate_residue_try(Expr* f, Expr* x, Expr* a, Expr* b, const char* method);

static Expr* residue_family_rational(Expr* f, Expr* x, Expr* a, Expr* b);  /* A */
static Expr* residue_family_fourier (Expr* f, Expr* x, Expr* a, Expr* b);  /* B */
static Expr* residue_family_trig    (Expr* f, Expr* x, Expr* a, Expr* b);  /* C */

/* shared helpers */
static bool  is_neg_pos_infinity_pair(Expr* a, Expr* b);   /* reuse nl_numeric */
static bool  is_full_period(Expr* a, Expr* b);             /* {0,2Pi} or {-Pi,Pi} */
static int   poly_degree_in(Expr* p, Expr* x);
static Expr* solve_complex_roots(Expr* Q, Expr* x, bool radicals);
static bool  pole_in_uhp(Expr* z);                         /* N[Im[z]] > tol */
static bool  pole_in_unit_disk(Expr* z);                   /* N[Abs[z]] < 1-tol */
static bool  pole_on_real_axis(Expr* z);                   /* |Im| < tol -> PV */
static bool  pole_on_unit_circle(Expr* z);                 /* |Abs-1| < tol -> PV/NULL */
static bool  has_real_root(Expr* Q, Expr* x);
static Expr* sum_residues(Expr* f, Expr* x, Expr** poles, size_t n, bool half);
static bool  contour_crosscheck(Expr* f, Expr* x, Expr* a, Expr* b, Expr* V);

Expr* builtin_integrate_contour_residue(Expr* res);  /* Integrate`ContourResidue */
void  integrate_residue_init(void);
```

**Export from `src/calculus/residue.h`:** `Expr* residue_compute(Expr* f, Expr* x,
Expr* x0)` (currently the body of `builtin_residue`) so the contour method calls
it directly without re-parsing.

**Modify `src/calculus/integrate.c`:** add `METHOD_RESIDUE` to the enum; map
`"Residue"` / `"ContourResidue"` in `method_from_string`; allow it in
`definite_parse_method`; in `integrate_definite`'s real-spec branch, try
`integrate_residue_try(...)` before `integrate_newton_leibniz_try(...)` when the
method permits; call `integrate_residue_init()` from `integrate_init()`;
`#include "integrate_residue.h"`.

**Build & tests:** add `../src/calculus/integrate_residue.c` to
`tests/CMakeLists.txt` COMMON_SRC (with the other `integrate_*.c`); add an
`integrate_residue_tests` executable block (copy `integrate_line_tests`). New
`tests/test_integrate_residue.c` following the `check_eq` +
`symtab_init()/core_init()` pattern; numeric verification via
`N[Abs[Integrate[f,{x,a,b}] - NIntegrate[f,{x,a,b}]]] < 1/10^6`. Top-level
`makefile` auto-globs `src/calculus/*.c` — no edit.

**Docs:** a `### Contour / residue-theorem definite integration` subsection in
`docs/spec/builtins/calculus.md` (after the LineIntegral section), plus a
changelog note under the current week's `docs/spec/changelog/<Monday>.md`.

## Phased delivery

- **Phase 1** — Family C (finite-range crosscheck, strongest); Family A for
  `deg Q ≤ 4` (radical closure); Family B for a positive rational literal `a`.
  All crosscheck-gated, NULL → NL.
- **Phase 1b** — half-line `[0,∞)` via even symmetry; Family B general positive
  real `a`, and `Cos/Sin` via `Re/Im[J]`.
- **Phase 2 (this scope)** — Family PV: indented-contour / half-residue for
  simple real-axis poles (∫ sin x / x = π); Family A `deg Q ≥ 5` when it factors
  to ≤ quartic.
- **Deferred** — keyhole / branch-cut contours (`∫₀^∞ x^{s-1}/(1+x)`): non-
  rational, out of scope for the rational-recognizer architecture.

## Top risks & mitigations

1. **Residue sum won't close** (radicals/`Root` survive `RootReduce`) → prefer
   radical roots (deg ≤ 4 / factor first); chop residual `Im`; if not a scalar,
   NULL → NL. Never ship a messy `Root`/`RootSum`.
2. **Im-sign misclassification** near the real axis → elevated-precision `N[]`,
   magnitude-scaled tol, conjugate-pair UHP-only sum; `|Im| < tol` → treat as
   on-axis (PV branch or NULL).
3. **Recognizer false-positive** (bad degree/convergence gate, hidden real pole)
   → every accepted value passes the `NIntegrate` crosscheck; mismatch → NULL.
4. **`NIntegrate` can't validate** infinite-range/oscillatory A/B → large-finite-
   window fallback with loosened tol; else NULL (lose a correct clean answer
   rather than return an unverified one).
5. **Dispatch/ownership regression** in `integrate_definite` → mirror
   `integrate_newton_leibniz_try`'s borrowed-args/owned-result contract; engage
   only for a single real spec; regression-test the existing NL/line suites.

## Worked test examples (expected closed forms)

Family A: `∫1/(1+x^2)=π`; `∫1/(1+x^4)=π/√2`; `∫1/(1+x^2)^2=π/2` (order-2 pole);
`∫x^2/(1+x^4)=π/√2`; `∫1/((x^2+1)(x^2+4))=π/6`.
Family B: `∫cos(x)/(1+x^2)=π/e`; `∫x sin(x)/(1+x^2)=π/e`;
`∫cos(a x)/(x^2+b^2)=(π/b)e^{-ab}` for `a,b>0` (else NULL).
Family C: `∫₀^{2π} 1/(2+cosθ)=2π/√3`; `∫₀^{2π} 1/(5−4cosθ)=2π/3`.
Family PV: `∫ sin(x)/x = π` (half-residue at 0); PV `∫1/(x^2−1)`.
Half-line: `∫₀^∞ 1/(1+x^4)=π/(2√2)`.

Negative controls (must return NULL → N–L handles / leaves unevaluated):
`∫1/(1+x^3)` (real axis pole); `∫1/Sqrt[1+x^4]` (branch point, not rational);
`∫₀^π 1/(2+cosθ)` (not a full period → Family C must not fire).

## Verification

1. `make -j` (foreground).
2. `cd tests/build && cmake .. && make -j integrate_residue_tests`, run it, grep
   for `FAIL:`.
3. Regression: `integrate_newton_leibniz_tests`, `integrate_line_tests`,
   `integrate_dispatch_tests`, `residue_tests` — confirm the reordered dispatch
   changed nothing for finite/antiderivative-driven integrals.
4. `valgrind --leak-check=full` on the test binary; diff against the macOS
   dyld/objc baseline (~12.8KB) — no Mathilda-source leaks.
5. Confirm `Method->"NewtonLeibniz"` still bypasses the residue path, and that
   the residue answer matches `NIntegrate` on every positive case.
