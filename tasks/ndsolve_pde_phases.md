# NDSolve PDE — remaining phases

Plan: `~/.claude/plans/ndsolve-pde-method-of-lines.md`. Each phase: implement +
extensive unit/stress tests + valgrind clean + docs/changelog + commit.

- [x] **Phase 1** — 1-D MoL front-end, 2nd-order central stencils, Dirichlet,
      temporal order 1&2, machine precision. (`ndsolve_mol.c`, 43 tests) — DONE (ceb4c57)
- [x] **Phase 2a** — Fornberg arbitrary-order stencils + `DifferenceOrder`
      (default 4); spatial-order cap lifted. (48 tests) — DONE (ffa6a5f)
- [x] **Phase 2b** — Neumann / Robin / Periodic boundary conditions.
      (62 tests) — DONE
- [x] **Phase 3** — compiled banded operator (linear PDE → A·U+s), exact
      Jacobian, banded LU solve, stiffness auto-select. ~10x speedup. (68 tests) — DONE
- [x] **Phase 4** — 2 spatial dimensions (tensor grid, 2-D stencils, Dirichlet). (81 tests) — DONE
      Tests: 2-D heat/wave separable eigenmodes.
- [x] **Phase 5** — MPFR MoL (non-stiff, `nd_solve_mpfr_mol`) + nonlinear
      (viscous Burgers) + **complex (Schrödinger, `nd_realify` + two-real-IF
      output)** — DONE (99 tests). All phases complete.
- [x] **Adaptive-implicit stepping** — variable-step BDF (orders 1-2) + Adams
      with predictor-corrector (Milne) local error control and Newton-failure
      step recovery; solves incompatible IC/BC corners that previously diverged
      (ndcf). (`ndsolve_implicit.c`, `ndsolve_adams.c`) — DONE.
- [x] **Higher-order BDF (VSVO, orders 1-5)** — exact nonuniform-mesh Lagrange
      coefficients, order ramp/hold, uniform-mesh stability. (`ndsolve_implicit.c`;
      +test_bdf_high_order) — DONE.
- [x] **Stiff arbitrary precision** — MPFR variable-order BDF (`mpfr_bdf_dir`):
      MPFR state/coeffs/residual/solve, double Jacobian. Stiff ODEs+PDEs at
      WP>machine. (`ndsolve_mpfr.c`; +test_bdf_mpfr_stiff) — DONE.
- [x] **2-D general BCs (Neumann/Robin)** — the tensor-grid solver
      (`nd_mol_solve_2d`) now accepts `a·u + b·u_n + r == 0` on each edge
      (Dirichlet/Neumann/Robin), eliminating boundary nodes into interior
      stencils (`nd_bc_eliminate_2d`); corners resolve via the transverse edge.
      Periodic + mixed derivatives stay deferred but are detected/reported.
      (+test_pde_2d_neumann, +test_pde_2d_robin_steady) — DONE.
- [ ] **2-D periodic + mixed derivatives**; **robustness cluster**
      (BackwardEuler maxsteps, complex ODE realification, x0!=xmin).
