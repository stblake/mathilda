# NDSolve PDE — remaining phases

Plan: `~/.claude/plans/ndsolve-pde-method-of-lines.md`. Each phase: implement +
extensive unit/stress tests + valgrind clean + docs/changelog + commit.

- [x] **Phase 1** — 1-D MoL front-end, 2nd-order central stencils, Dirichlet,
      temporal order 1&2, machine precision. (`ndsolve_mol.c`, 43 tests) — DONE (ceb4c57)
- [ ] **Phase 2a** — Fornberg arbitrary-order stencils + `DifferenceOrder`
      option (default 4); lift spatial-order cap. Tests: order-2/4/6 spatial
      convergence (ratios ~4/16/64), eigenmode exactness pinned to order 2.
- [ ] **Phase 2b** — Neumann / Robin / Periodic boundary conditions.
      Tests: Neumann steady state, periodic traveling wave, Robin.
- [ ] **Phase 3** — compiled banded numeric operator (linear PDE → A·U+s),
      exact banded Jacobian, banded LU solve, stiffness auto-select. Tests:
      fast-path vs symbolic parity, efficiency, auto-method.
- [ ] **Phase 4** — 2 spatial dimensions (tensor grid, 2-D stencils/BCs).
      Tests: 2-D heat/wave separable eigenmodes.
- [ ] **Phase 5** — MPFR MoL (stencils+operator+output), plus nonlinear
      (Burgers) and complex (Schrödinger) stress. Tests: high-precision digits,
      Burgers, Schrödinger norm conservation.
