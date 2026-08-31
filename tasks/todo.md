# DSolve — Lagrange/d'Alembert (1a), parametric general solution

`y == x·φ(y') + ψ(y')` → parametric `{x=X(t,C), y=X·φ+ψ}`, t=slope. New parametric
substrate path (mirror the implicit path). General solution only; singular +
IVP-fitting deferred.

## Substrate (dsolve_common.{c,h}) — mirror implicit quartet
- [ ] dsolve_verify_parametric(P, X, Y, tname): sub y'->D[Y,t]/D[X,t], y[x]->Y, x->X; permissive + PossibleZeroQ fallback
- [ ] dsolve_assemble_parametric: {{x->Function[{t},X], y->Function[{t},Y]}}, ds_rename_param
- [ ] dsolve_run_parametric(P, fn): unpack DSolve`Param[X,Y], verify, assemble; decline if ncond>0
- [ ] dsolve_method_builtin_parametric(res, fn)
- [ ] parameter symbol: first of {t,s,u} not in equation (ds_contains guard)
- [ ] declarations in dsolve_common.h

## New file src/calculus/dsolve_lagrange.c
- [ ] dsolve_lagrange_try: recognize (Clairaut-style), extract φ=D[Yexpr,x], ψ=Yexpr-x·φ (free of x; ψ free of Y); decline φ≡p; solve linear ODE for X(t) via dsolve_linear_factor_solve; Y=X·φ+ψ; return DSolve`Param[X,Y]
- [ ] builtin_dsolve_lagrange via dsolve_method_builtin_parametric; dsolve_lagrange_init

## Wiring src/calculus/dsolve.c
- [ ] enum DS_LAGRANGE; "Lagrange" map; externs; cascade line (run_parametric after clairaut); pinned case; init
- [ ] tests/CMakeLists.txt: add dsolve_lagrange.c to mathilda_common

## Tests
- [ ] test_dsolve.c: t_method_lagrange, t_lagrange_more, t_lagrange_declines; register
- [ ] test_dsolve_stress.c: lagrange_ok generator + t_stress_lagrange; register

## Docs
- [ ] calculus.md row; changelog section; DSOLVE_PLAN.md §1a flip

## Gates
- [ ] make -j clean; REPL spot-checks (parametric output, residual, Clairaut unchanged, decline)
- [ ] ctest -R dsolve 3/3; make check-c99; valgrind baseline; graph rebuild

## Review

Done. Lagrange/d'Alembert (§1a) implemented with a new **parametric substrate
path** mirroring the implicit path, plus one new method file.

- Substrate (`dsolve_common.{c,h}`): `dsolve_run_parametric` /
  `_verify_parametric` / `_assemble_parametric` / `_method_builtin_parametric`.
  Try-fn returns `DSolve\`Param[X, Y, t]`; verify substitutes `y'=Y'(t)/X'(t)` into
  the residual (permissive + PossibleZeroQ fallback); output
  `{{x->Function[{t},X], y->Function[{t},Y]}}`; IVP declines (`ncond>0`).
- `src/calculus/dsolve_lagrange.c`: Clairaut-style recognition; extract
  `φ=D[Yexpr,x]`, `ψ=Yexpr−xφ`; decline `φ≡p` (Clairaut) and genuinely-linear
  (`φ` const ∧ `ψ` affine → LinearFirstOrder); linear ODE for `X(t)` via
  `dsolve_linear_factor_solve`; `Y=Xφ+ψ`. Cascade slot after Clairaut.
- Parameter symbol: collision-safe bare symbol (first of {t,s,u,w,r,q} not in eqn).
- Tests: `t_method_lagrange` / `t_lagrange_more` / `t_lagrange_declines` (unit,
  parametric verify) + `t_stress_lagrange` (8-case forward generator).

Verified: `y==2xy'+(y')^2 → {x->(C[1]-2/3 t^3)/t^2, y->(2C[1]-1/3 t^3)/t}` (residual
PossibleZeroQ True); transcendental (Log) forms verify via PossibleZeroQ; Clairaut
still explicit lines (unchanged); pinned Lagrange declines Clairaut/linear/IVP.

Gates: `ctest -R dsolve` 3/3 green (19.5+6.5+10.1 s); `make check-c99` exit 0;
valgrind 1×==6× = 13,440 def + 6,312 indir (one-time engine baseline, no per-call
leak); graph rebuilt. No user correction; no lesson added (the linear-decline guard
was found via my own spot-check, not a correction). LSP `ATTR_READPROTECTED` /
missing-include noise is stale — GCC/CMake build is clean.
