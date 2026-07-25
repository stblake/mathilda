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
