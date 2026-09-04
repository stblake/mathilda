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

---

# DSolve M6 (Phase 2) — PDEClassify  ✅ DONE (2026-09-04)

`PDEClassify[eqn, u, {v1, v2}]` (`src/calculus/dsolve_pdeclassify.c`) — a standalone
top-level classifier: the discriminant `Δ = B² − 4 A C` of the principal part
(`A u_{v1 v1} + B u_{v1 v2} + C u_{v2 v2}`) → `"Hyperbolic"` (Δ>0) / `"Parabolic"`
(Δ=0) / `"Elliptic"` (Δ<0). Reuses `dsolve_parse` for the residual + the same
2nd-order coefficient extraction as pde2. Sign decided via a small ladder
(`ds_is_zero` → numeric type → `Sign[]`); an undecidable-sign / parameter-dependent
Δ (Tricomi `y u_xx + u_yy`) leaves the call unevaluated (honest decline).

- [x] `builtin_pdeclassify`: parse, extract A,B,C (2nd-order terms), Δ, sign ladder → String
- [x] only principal part matters (lower-order terms ignored); declines 1st-order / undecidable
- [x] register top-level `PDEClassify` (Protected, docstring); init from dsolve_init
- [x] tests `t_pdeclassify` (wave/heat/Laplace/mixed/lower-order/Tricomi-decline/1st-order-decline)
- [x] `make check-c99` PASS; dsolve_tests green; valgrind baseline 13,440/6,312, **0** pdeclassify frames
- [x] docs: DSOLVE_PLAN.md §2b + M6, calculus.md prose, changelog

Design note: standalone builtin (not a solver → outside the `DSolve\`` namespace),
matching the plan's bare `PDEClassify` naming. First cut: linear in the 2nd-order
terms, two independent variables, constant-signed discriminant.

Still open in §2b: wave-IVP d'Alembert formula; heat kernel / `Erf`; lower-order/
inhomogeneous 2nd-order (telegraph); quasilinear/nonlinear first-order (§2a).

---

# DSolve M6 (Phase 2) — wave d'Alembert IVP  ✅ DONE (2026-09-04)

`DSolve`WaveDAlembert` (`src/calculus/dsolve_wave.c`) — the whole-line wave IVP
`u_tt==c² u_xx`, `u(x,t0)==f(x)`, `u_t(x,t0)==g(x)` → d'Alembert
`½(f(x−cτ)+f(x+cτ)) + 1/(2c)∫_{x−cτ}^{x+cτ} g(K)dK`, `τ=t−t0`. Auto-dispatched
(new `is_pde` cascade route) + pinned.

- [x] parse (ICs arrive as equations → `neq==3`); sort PDE vs 2 ICs; fixed var=time, free=space
- [x] sign-robust IC rhs extraction; wave speed from principal coeffs (reject parabolic/elliptic)
- [x] d'Alembert body build (dummy `K`, collision-guarded; integral kept for undefined `g`)
- [x] own multi-equation verify (rebuild with `f=Cos, g=Sin`; PDE residual + both ICs → 0)
- [x] own assemble `{{u->Function[{x,t},…]}}`; auto route in dsolve.c + pinned builtin
- [x] tests `t_wave_dalembert` (concrete verified; undefined at displacement IC; auto+pinned; declines)
- [x] `make check-c99` PASS; all 3 DSolve ctest suites green
- [x] valgrind: decline paths at baseline 13,440/6,312 (dsolve_wave leak-free); solving paths
      inherit only the documented per-call Integrate-engine leak (LinearFirstOrder/AlmostLinear too)

Key gotcha → memory: a PDE initial condition `u[x,t0]==f[x]` is a 2-arg point
condition that `ds_is_condition` (single-var funcapps only) does NOT split off, so a
wave IVP reaches the solver as `neq==3` (PDE + 2 ICs as equations); the method must
sort them itself. Verifying a d'Alembert solution with undefined `f,g` is impossible
(no Leibniz rule for the unevaluated `Integrate[g[K],…]`), so verify rebuilds with
`f=Cos, g=Sin`.

Still open in §2b: heat kernel / `Erf`; inhomogeneous / half-line / `Piecewise` wave;
lower-order/inhomogeneous 2nd-order (telegraph); quasilinear/nonlinear first-order (§2a).

---

# DSolve M6 (Phase 2) — heat-kernel Cauchy problem  ✅ DONE (2026-09-04)

`DSolve`HeatKernel` (`src/calculus/dsolve_heat.c`) — the whole-line heat Cauchy
problem `u_t==k u_xx`, `u(x,t0)==f(x)` → the heat-kernel convolution
`1/Sqrt[4πkτ] ∫_{-∞}^{∞} f(K) Exp[-(x−K)²/(4kτ)] dK`, `τ=t−t0`. Auto (new `is_pde`
route, `neq==2`) + pinned.

- [x] sort PDE vs 1 IC (reuse wave's approach); fixed var=time, free=space; sign-robust f rhs
- [x] pure-heat gate: require u_t + u_xx only (no u_tt/u_x/u_xy/u/forcing), k>0
- [x] build convolution with Integrate node RAW (no integration attempt → no hang, no engine leak)
- [x] verify at KERNEL level (`G_t − k G_xx == 0`; confirms k); IC by delta-convergence (theory)
- [x] auto route + pinned builtin `DSolve`HeatKernel`
- [x] tests `t_heat_kernel` (auto+pinned solve & carry the convolution; kernel solves PDE; declines)
- [x] `make check-c99` PASS; all 3 DSolve ctest suites green
- [x] valgrind EXACTLY at baseline 13,440/6,312 — dsolve_heat fully leak-clean (raw-Integrate
      avoids the per-call Integrate-engine leak that wave inherits)

Key decisions → memory: the Gaussian convolution is nonelementary (Mathilda can't do
the improper integral; Simplify on it HANGS), so build the Integrate node raw
(unevaluated) and verify the KERNEL not the convolution; the IC (delta-convergence)
rests on theory, not back-substitution. First bug found+fixed: time-derivative
multi-index was swapped (space=v1 → time=v2 → u_t is Derivative[0,1], not [1,0]).

Still open in §2b: Erf-producing step/box heat data; finite-interval Fourier series;
inhomogeneous / half-line / `Piecewise` wave; telegraph; quasilinear/nonlinear 1st-order (§2a).
