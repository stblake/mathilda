# NDSolve PDE — remaining phases

Plan: `~/.claude/plans/ndsolve-pde-method-of-lines.md`. Each phase: implement +
extensive unit/stress tests + valgrind clean + docs/changelog + commit.

- [x] **Phase 1** — 1-D MoL front-end, 2nd-order central stencils, Dirichlet,
      temporal order 1&2, machine precision. (`ndsolve_mol.c`, 43 tests) — DONE (ceb4c57)
- [x] **Phase 2a** — Fornberg arbitrary-order stencils + `DifferenceOrder`
      (default 4); spatial-order cap lifted. (48 tests) — DONE (ffa6a5f)
- [x] **Phase 2b** — Neumann / Robin / Periodic boundary conditions.
      (62 tests) — DONE
- [ ] **Phase 3** — compiled banded numeric operator (linear PDE → A·U+s),
      exact banded Jacobian, banded LU solve, stiffness auto-select. Tests:
      fast-path vs symbolic parity, efficiency, auto-method.
- [ ] **Phase 4** — 2 spatial dimensions (tensor grid, 2-D stencils/BCs).
      Tests: 2-D heat/wave separable eigenmodes.
- [ ] **Phase 5** — MPFR MoL (stencils+operator+output), plus nonlinear
      (Burgers) and complex (Schrödinger) stress. Tests: high-precision digits,
      Burgers, Schrödinger norm conservation.
