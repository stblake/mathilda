# NDSolve: systems of 1-D PDEs + upwind schemes

## Goal
Extend the method-of-lines engine to solve **coupled 1-D PDE systems**
(`{h, u}` shallow-water, coupled reaction-diffusion, coupled waves) and add
two selectable **upwind** spatial schemes for hyperbolic/advective terms.

## Design decision
- Keep the proven scalar path in `nd_mol_solve` **untouched** (152-check suite
  is the regression gate). Add a dedicated `nd_mol_solve_system` for N>1,
  reusing all existing static helpers (they are already parameterized by
  fname/ysym/torder). Route to it when `{...}` has >1 function.

## Phase 1 — N-equation systems, centered stencils
- [ ] Per-function descriptor array (fname, torder, sord, ic, BCs, periodic, layout base)
- [ ] Global block-per-function reduced-state layout
- [ ] Equation classification: interior (owner = fn whose top-time-deriv appears),
      IC, BC, periodic — dispatched per function
- [ ] Solve each evolution eq for its top time-derivative -> G_f
- [ ] Coupled RHS build: substitute literals of ALL functions (values + spatial
      stencils + lower time-derivs) into each G_f
- [ ] Per-function BC elimination + IC sampling
- [ ] Result: one InterpolatingFunction per function -> {{u->IF, v->IF, ...}}
- [ ] Guards: complex systems, MPFR systems, implicit/coupled mass matrix -> warn
- [ ] Tests: coupled reaction-diffusion (manufactured), coupled waves (eigenmode),
      shallow-water small-perturbation (gravity-wave speed), missing-BC guard

## Phase 2 — upwind schemes
- [ ] B1: Lax-Friedrichs artificial viscosity (default for hyperbolic systems)
- [ ] B2: sign-biased donor-cell upwind (scalar, sharp)
- [ ] Scheme selection via Method/DifferenceOrder option
- [ ] Tests: linear-advection translate (centered vs upwind, convergence),
      dam-break bounded & positive under B1

## Review — DONE (2026-07-26)

Both phases complete, verified, no regressions.

**Phase 1 (systems):** added `nd_mol_solve_system` in `ndsolve_mol.c`; scalar
path untouched. `NDSolve[{u,v,...}, ...]` routes there. Block-per-function
reduced state in one global ODE vector; each function's RHS substitutes all
functions' stencils/node-values (coupling). One InterpolatingFunction per
function → `{{u->if,...}}`. Guards: missing BC/IC, coupled mass matrix,
complex, 2-D system → unevaluated.

**Phase 2 (upwind):** donor-cell (scalar, `DifferenceOrder->1`/`"Upwind"->True`,
wind = `-∂G/∂u_x` per node) + Lax-Friedrichs viscosity (systems + scalar,
`"LaxFriedrichs"->True`). Centered default unchanged.

**Verified vs exact references:** decoupled/coupled reaction-diffusion
(manufactured), linearized shallow-water gravity wave (`c=√(gH)`), coupled 2nd-
order wave (normal-mode, exact to semi-discrete eigenvalue), linear advection
(centered exact / donor phase-accurate-diffusive-convergent / reversed-wind
stable), top-hat monotonicity (centered rings, upwind bounded), guards.

**Tests:** `tests/test_ndsolve_pde.c` 151→172 checks, 0 failures. ODE (27) and
classical (87) suites unchanged. Valgrind: only macOS ObjC/dyld baseline noise
(no stacks reference the new code).

**Docs:** `docs/spec/changelog/2026-07-20.md` + `docs/spec/builtins/numerical-
calculus.md` updated; file header scope note refreshed.

**Out of scope (stated):** incompressible Navier–Stokes (needs pressure-Poisson
DAE constraint).
