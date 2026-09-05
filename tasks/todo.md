# M12 — Second-order symmetry & integrating-factor DSolve methods

Goal: solve the ~110 nonlinear 2nd-order ODEs in the 12000.org "SymPy-failed"
corpus via algorithmic (not overfit) methods, validated by 10–20 generalised
problems per seed deficiency.

## Tasks

- [ ] Read reference paper physics/9703082 (2nd-order symmetry methods) + Cheb-Terrab
      & Roche integrating-factor equations
- [ ] Study substrate: dsolve_common.{c,h}, dsolve_lie.c, dsolve_liouville.c,
      dsolve_autonomous.c, dsolve_reduce_order (order-reduction template)
- [ ] Implement M12a `SecondOrderIntegratingFactor` (src/calculus/dsolve_muode.c)
      — μ classes _mu_y1, _mu_x_y1, _mu_xy, _mu_y_y1 → first integral → run_implicit
- [ ] Implement M12b `SecondOrderSymmetry` (src/calculus/dsolve_lie2.c)
      — 2nd-prolongation determining PDE, ansatz/NullSpace, canonical-coord reduction
- [ ] Wire cascade in dsolve.c (order-2 group, after liouville/autonomous, before
      frobenius); externs; init chain
- [ ] Register builtins (ATTR_PROTECTED + docstrings); sym_names if needed
- [ ] Tests: tests/test_dsolve_m12_stress.c (forward-generator families) + pinned
      cases in test_dsolve.c; tests/CMakeLists.txt
- [ ] Anti-overfit: N4/u=y² (15), N3/scaling (12), μ(y′)/μ(x,y′) (15),
      μ(y,y′) (12), projective (10); D2 hang regression
- [ ] Docs: DSOLVE_PLAN.md M12 block + §1d entries; docs/spec/builtins/calculus.md;
      docs/spec/changelog/2026-08-31.md; docs/design/dsolve_lie_symmetry.md
- [ ] Verify: make -j, ctest (no regression), make check-c99, valgrind spot-checks,
      corpus re-measure, D2 hang gone

## Review

**Done (M12 end-to-end).** `DSolve`SecondOrderSymmetry` (`src/calculus/dsolve_lie2.c`)
— nonlinear 2nd-order Lie point-symmetry method: 2nd-prolongation determining system
(`NullSpace` ansatz) → canonical-coordinate order reduction → cascade → invert. Wired
into the scalar cascade after `Liouville`, before Frobenius. Registered pinned builtin
+ docstring + `ATTR_PROTECTED`.

**Solves** (numerically verified): N3, N4, N5, projective `x²(x+y)y''=(xy'−y)²`, Eq.3,
plus three arbitrary-coefficient anti-overfit families (projective `x³y''=a(y−xy')²`,
projective+linear, scaling `k x²yy''+y²=x²y'²`) — all 100% in
`tests/test_dsolve_m12_stress.c` + pinned unit in `test_dsolve.c`.

**Two hard lessons (see lessons.md):**
1. `dsolve_run`'s symbolic verify KEEPS an undecidable residual → a heuristic reduction
   producing logs/radicals must NUMERICALLY back-substitute before returning (N1/N2 were
   WRONG answers that passed symbolic verify).
2. The evaluator re-invokes a declining DSolve builtin ~3× per call → an expensive
   declining method needs its own per-toplevel decline memo (the dispatcher fail-memo
   misses it because the equation re-normalizes each pass).

**Bounded, never hangs:** all recursive sub-solves `TimeConstrained`-wrapped + wall-clock
deadline + undefined-fn gate + linearity gate. Declines are clean fall-throughs (~5 s).

**Gates:** all 4 DSolve ctest suites + `make check-c99` green. Per-call leak is the
inherited Integrate/Solve-engine leak (as `AlmostLinear`); lie2's own ownership is clean.

**Out of scope / follow-ups:** pre-existing `DSolve`Separable` inversion hang on
`q'==-2q(1+q)(1+2q)/r`; linear var-coeff timeouts (EQ4/5/9/10/11 — Kovacic/special-fn,
M14 territory). Roadmap M13–M17 in `DSOLVE_PLAN.md`.

## M13/M14 session addendum

- **M13 (Abel/AIR): DEFERRED** — investigated thoroughly; the corpus `[_rational,_Abel]`
  set (≈9/10 first-kind are `y'=f3 y³+f2 y²`) is not the constant-invariant class
  (canonical `G0/G3` non-constant; reciprocal→Chini(n=−1) also declines). Needs the
  full Abel Invariant Rational hierarchy (research-grade; Nasser's own solver + SymPy
  both fail on precisely these). Not hacked (would risk wrong answers/overfit).
- **M14 (ChangeOfVariable): DONE** — `src/calculus/dsolve_changevar.c`. 2nd-order linear
  with transcendental coeffs → `t=φ(x)` (Cos/Sin/Tan) rationalizes → recurse → back-sub.
  Flagship: `y''+Cot[x]y'+k(k+1)y==0 → Legendre`. Numerically verified on the original
  (Legendre-Q has Log → undecidable residual). Bounded (TimeConstrained + deadline +
  decline memo + re-entry guard). Anti-overfit: 2 families 9/9 (`test_dsolve_m14_stress.c`).
  Valgrind leak-flat vs baseline. All 5 DSolve suites + check-c99 green.
