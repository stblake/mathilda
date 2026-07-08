# Integrate with Assumptions + contour-integration families

Goal: `Integrate[f, {x,a,b}, Assumptions -> {...}]` accepts the option and
evaluates parametric definite integrals via new/existing residue methods.

Flagship examples:
1. `Integrate[Cos[k x]/(x^2+a^2), {x,-Inf,Inf}, Assumptions->{a>0,k>0}]` → `(Pi E^(-a k))/a`  — Fourier family (existing) + symbolic params
2. `Integrate[Exp[a x]/(Exp[x]+1), {x,-Inf,Inf}, Assumptions->0<a<1]` → `Pi/Sin[Pi a]`  — NEW rectangular / quasi-periodic-in-Exp[x] contour
3. `Integrate[1/(1 + x^n), {x, 0, Infinity}, Assumptions->n>1]` → `Pi/(n Sin[Pi/n])`  — NEW sector contour (angle 2Pi/n)
4. `Integrate[x^(1/3)/(x^2 + 1), {x, 0, Infinity}]` → `Pi/Sqrt[3]`  — NEW keyhole / Mellin contour ∫_0^∞ x^(s-1) R(x) dx

## Design (shared)
- Classify poles / kernel-frequency by **sign-consistent numeric instantiation** σ
  of the free parameters (built from the assumptions' sign/interval constraints).
- Keep residue arithmetic symbolic; **cross-check** the closed form against
  `NIntegrate[f/.σ, spec/.σ]` at the instantiation point — this restores the rigor
  lost to symbolic classification and is the correctness gate for every family.
- No global `$Assumptions` mutation — thread the assumptions expr into the engine.

## Tasks
- [x] Phase 0: parse `Assumptions -> ...` in `integrate.c`; thread into `integrate_residue_try`.
- [x] Phase 1: σ-instantiation + σ-aware `res_reim`; PowerExpand poles + I->-I
      conjugation for symbolic Fourier; NO NIntegrate crosscheck (correct-by-
      construction; refuse unconstrained params). → Ex1 DONE.
- [x] Phase 4: `mellin_core` (keyhole, -π Σ z^(s-1)Res/sin(πs)). → Ex4 DONE.
- [x] Phase 2: `residue_family_rectangular` (w=Exp[x] substitution → mellin_core). → Ex2 DONE.
- [x] Phase 3: `residue_family_sector` (x^m/(c+x^n), symbolic n). → Ex3 DONE.
- [~] Ex5 log-keyhole DEFERRED: needs Arg/Log of symbolic on-circle poles in the
      (0,2π) branch (poles return as Cos[a]±Sqrt[..], Residue at E^(I a) = 0);
      no such assumption-aware complex-arg machinery exists. Answer is 0.
- [ ] Phase 5: extensive unit tests (tests/test_integrate_contour_assume.c).
- [ ] Phase 6: docs (docs/spec/builtins/calculus) + changelog + valgrind.
- [ ] Polish: recognise Power[E, linear] as Exp kernel (Ex1 Exp[I k x] form).

## Review (2026-07-08)
Shipped Ex1–Ex4 + rational-symbolic, all correct-by-construction, valgrind-clean,
full residue/dispatch/integrals suites green.

- `integrate.c`: `Assumptions -> …` accepted on definite (multi-option, any order)
  and indefinite forms; threaded into `integrate_residue_try` (new param).
- `integrate_residue.c`: sign-consistent instantiation `g_inst` (+ `g_all_pos`,
  `g_bounds`); σ-aware `res_reim`; `res_powerclean`/`res_conjugate`; `mellin_core`
  (keyhole) + `residue_family_mellin`/`_rectangular`/`_sector`; convergence gates
  on assumption-**guaranteed** bounds via `param_interval` (refuse under-constrained).
- NO NIntegrate crosscheck (per user directive: bugs must surface). Refuse a
  two-sided-unbounded parameter. Family A poles PowerExpand-cleaned → `π/a`.
- Tests: 7 new functions in `tests/test_integrate_residue.c`.
- Docs: `docs/spec/builtins/calculus.md` + changelog `2026-07-06.md`.

### Deferred / gaps
- Ex5 log-keyhole `∫₀^∞ Log[x]/(x²−2x Cos[a]+1)` (= 0): needs assumption-aware
  `Arg`/`Log` `(0,2π)`-branch reasoning for symbolic on-circle poles (poles come
  back as `Cos[a]±Sqrt[..]`; `Residue` at reparametrised `E^(I a)` = 0). No such
  machinery exists — would be a research-grade addition.
- `Exp[I k x]` (`Power[E,…]`) Fourier kernel form not recognised (`Cos`/`Sin` are);
  match_kernel + the family's kernel-arg extraction assume a `Cos/Sin/Exp` head.

---

## Review: Integrate`DiffUnderInt (differentiation under the integral sign) — 2026-07-08

**Delivered.** New definite-integration method `Integrate`DiffUnderInt` (Feynman's
trick), reachable via `Method -> "DiffUnderInt"` / `"DifferentiationUnderIntegral"`
and tried last in the automatic definite cascade (after residue + Newton-Leibniz).
New file `src/calculus/integrate_diffunderint.{c,h}`; wired into `integrate.c`
(method enum/parse, `integrate_definite`, `integrate_init`); test
`tests/test_integrate_diffunderint.c` (+CMake target, +COMMON_SRC); docs in
`docs/spec/builtins/calculus.md` and changelog `docs/spec/changelog/2026-07-06.md`.

**How it works.** Differentiate w.r.t a parameter → evaluate the simpler inner
integral J → integrate J back over the parameter → fix the constant from an exact
base value. Verification is SYMBOLIC (`Simplify[D[I,p]-J]===0` + exact base) — no
NIntegrate (project rule); the conditional-convergence pitfall is caught for free.

**Key discovery / design.** The general `Integrate` HANGS uninterruptibly (from
inside a builtin — `TimeConstrained` does not bound a nested `evaluate`) on the
inner Laplace/Fourier/Gaussian forms Feynman produces, and cannot do Gaussians at
all. So the method is FAMILIES-ONLY: it computes the standard inner integrals with
its own closed-form evaluators (Laplace/Fourier half-line, sinc/Frullani,
even-rational half-line via v=x^2+Apart+Beta) and declines (fast, unevaluated,
never wrong, never hangs) anything outside them.

**Coverage of the 24 target examples: 12 correct**, all fast (<=6 s), zero hangs:
#1,2,3,5,6,7,13,14,17,20,21,22. (#5 returns a `Sqrt[b^2]` form that equals
`ArcTan[b/a]` for b>0.)

**Deferred (return unevaluated, fast) — each needs more machinery:**
- Gaussian #4,10,12,24 — engine can't integrate `e^{-x^2}` and `Limit[Erf[x],x->Inf]`
  is unknown; needs a Gaussian moment-integral evaluator + Erf/Erfc support.
- Trig-period {0,Pi}/{0,Pi/2} #9,15 — needs a Weierstrass-substitution definite
  evaluator (the engine hangs on trig over a period).
- Piecewise #8,19 — result is case-split / `Min`; needs the Stage-C branch engine.
- #23 (sinc with decay) and #16 (log-rational base) — need a complex-Log→ArcTan
  reduction the engine lacks; #18 needs a finite-radical family.

**Verification.** `integrate_diffunderint_tests` green; leak-clean (0 frames in
integrate_diffunderint via macOS `leaks`; the 3 residual leaks are pre-existing
engine allocations); `integrate_newton_leibniz_tests` green (no cascade
regression); strict-C99 `-Wall -Wextra` clean.
