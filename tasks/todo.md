# DSolve M6 (Phase 2) — second-order constant-coefficient linear PDE

Plan: `/Users/user/.claude/plans/mossy-gliding-blum.md`

Delivers `DSolve`PDELinearSecondOrder` — homogeneous, principal-part-only,
constant-coefficient 2nd-order linear PDE via operator factoring (trial
`u=f(v2+λ v1)`, λ-quadratic `A λ²+B λ+C=0`). One method covers hyperbolic
(wave/d'Alembert), elliptic (Laplace, complex chars), parabolic (repeated root).

## Phase 1 — Generalize shared PDE verify (`dsolve_common.c`) ✅ DONE
- [x] `pde_term_order`/`pde_scan_order(e,u)` — max i+j over `Derivative[i,j][u][v1,v2]` nodes
- [x] verify: enumerate all (i,j), i+j=maxord..1 → `D[bodyC,{v1,i},{v2,j}]`, then u→bodyC
- [x] verify: bodyC via ReplaceAll with a *list* {C[1]:>Sin, C[2]:>Cos, C[3]:>Exp, C[4]:>#^2}
- [x] backward compatible: first-order PDE tests unchanged

## Phase 2 — New method `dsolve_pde2.c` ✅ DONE
- [x] `dsolve_pde2_solve`: guard, extract A,B,C,Dc,Ec,Fc,f; gate (homogeneous/principal/const)
- [x] factor via `dsolve_analyze_roots` (A≠0); swap v1↔v2 (A=0,C≠0); `C[1][v1]+C[2][v2]` (mixed)
- [x] distinct roots → F+G; repeated → F(w)+v1 G(w)
- [x] `builtin_dsolve_pde2` + `dsolve_pde2_init` (PDELinearSecondOrder, ATTR_PROTECTED, docstring)

## Phase 3 — Wire-up ✅ DONE
- [x] `dsolve.c`: externs, cascade line after pde1, `dsolve_pde2_init()`
- [x] `tests/CMakeLists.txt`: add `../src/calculus/dsolve_pde2.c` to COMMON_SRC

## Phase 4 — Tests (`test_dsolve.c`) ✅ DONE
- [x] wave/d'Alembert, Laplace (elliptic), u_xy==0, repeated-root, asymmetric-distinct
- [x] pinned `DSolve`PDELinearSecondOrder`; decline a 1st-order PDE + ODE
- [x] existing PDE tests still green (t_pde_transport/forcing/zeroth_order, t_stress_pde)

## Phase 5 — Gates + docs ✅ DONE
- [x] `make -j` clean; `make check-c99` clean; dsolve ctest suites green (all pass)
- [x] valgrind: total at documented baseline 13,440/6,312; pde2 frames only in pre-existing
      zero_test/numericalize uninitialised-value contexts (no new leaks)
- [x] REPL spot-checks (wave/Laplace/mixed/repeated/asymmetric, independent residual=0)
- [x] DSOLVE_PLAN.md §2b + M6; docs/spec/builtins/calculus.md; changelog 2026-08-31.md; memory

## Review — DONE (2026-09-04)

Delivered `DSolve`PDELinearSecondOrder` — the first Phase-2 second-order PDE
method — as a clean, fully back-substitution-verifiable increment.

- **One method, three discriminant types.** Operator factoring via the
  characteristic quadratic `A λ²+B λ+C==0` (trial `f(v2+λ v1)`) covers hyperbolic
  (wave/d'Alembert `C[1][x-c t]+C[2][x+c t]`), elliptic (Laplace, complex chars
  `C[1][y-I x]+C[2][y+I x]` — matches Mathematica), and parabolic (repeated root
  `C[1][w]+v1 C[2][w]`). Reuses `dsolve_analyze_roots`, so complex/repeated λ are
  free. This realizes §2b `PDEHyperbolicGeneral` in full generality.
- **Substrate reuse + one generalization.** Reused `dsolve_run_pde` /
  `dsolve_assemble_pde` / `dsolve_method_builtin_pde` unchanged. Generalized only
  the shared `dsolve_verify_pde` to arbitrary derivative order (scanned from the
  residual — `max_order` is 0 for PDEs) and up to 4 arbitrary functions with
  distinct test functions; backward compatible with the first-order PDE tests.
- **Independently validated:** all residuals (wave, Laplace, repeated, and the
  asymmetric distinct-root `2u_xx+5u_xy+2u_yy`) `PossibleZeroQ`-zero under concrete
  substituted arbitrary functions — not vacuous.

**Gotcha → memory:** `max_order` is never populated for PDEs (`ds_scan` /
`ds_match_funcapp` only recognise single-index `Derivative[m][u][x]`, arg_count 1;
a PDE term `Derivative[i,j][u][v1,v2]` has arg_count 2 and is skipped). Any PDE
method / verifier must scan the residual for its order rather than read
`P->max_order`.

Deferred (documented future in §2b): inhomogeneous forcing + lower-order terms
(telegraph/damped wave — the full symbol must factor into first-order operators);
the wave-IVP d'Alembert formula (initial data, half-line, `Piecewise`); heat
kernel; `PDEClassify` (discriminant); quasilinear / nonlinear first-order (§2a).

---

# DSolve M6 (Phase 2) — separation of variables  ✅ DONE (2026-09-04)

`DSolve`SeparationOfVariables` (`src/calculus/dsolve_pdesep.c`, pinned-only) —
separated product mode `u == X(v1) Y(v2)` for a homogeneous, constant-coefficient,
**no-mixed-term** linear PDE. Divide by `X Y` → two constant-coefficient ODEs in a
separation constant λ (`Σ a_i X^(i) − λ X == 0`, `Σ b_j Y^(j) + (e+λ) Y == 0`),
each solved by recursing into the scalar cascade; the one redundant overall scale
absorbed (fix a first-order side's lone constant to 1). Pinned-only because the
general solution is a superposition over λ.

- [x] `dsolve_pdesep_solve`: order scan + per-(i,j) coefficient extraction; gate
      (linear / constant-coeff / no-mixed / homogeneous / has both x- and y-deriv)
- [x] build + solve the two ODEs via recursion into `DSolve` (applied form,
      `dsolve_extract_applied_bodies`); absorb scale, renumber, λ → C[k]
- [x] pinned builtin + init (ATTR_PROTECTED, docstring); NOT in the auto cascade
- [x] tests `t_pdesep` (heat, heat-with-k, Helmholtz residual=0; decline mixed /
      inhomogeneous / ODE-in-disguise); all 3 DSolve ctest suites green
- [x] `make check-c99` PASS (root); valgrind at baseline 13,440/6,312, **zero**
      `dsolve_pdesep` frames in the report
- [x] docs: DSOLVE_PLAN.md §2b + M6, calculus.md (table + prose), changelog
- [x] independently verified residual = 0 (heat, heat-k, Helmholtz)

Design note: a product mode ≠ general solution, so pinned-only (matches
`FirstOrderPowerSeries` / `EigenvalueProblem`). The `Sqrt[-4 λ]/2` form in some
outputs is constcoeff's quadratic-formula spelling (cosmetic) and the sign of λ is
a free-constant relabeling — the family is identical and back-sub verified.

Still open in §2b: wave-IVP d'Alembert formula; heat kernel / `Erf`; `PDEClassify`;
lower-order/inhomogeneous 2nd-order; quasilinear/nonlinear first-order (§2a).
