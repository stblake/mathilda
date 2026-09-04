# DSolve M6 — First-order nonlinear PDEs: PDEQuasilinear (Lagrange) + PDEClairaut

Plan: `/Users/user/.claude/plans/let-s-continue-our-implementation-purring-cherny.md`

## Tasks

- [x] 1. Baseline build + confirmed the three target PDEs currently decline
- [x] 2. Substrate: implicit-PDE path in `dsolve_common.{c,h}`
      (`dsolve_run_pde_implicit` dispatching `PDEImplicit`/`PDEExplicit`/`PDEBranches`;
      `dsolve_verify_pde_implicit`, `dsolve_assemble_pde_implicit`,
      `dsolve_method_builtin_pde_implicit`)
- [x] 3. `src/calculus/dsolve_pdeclairaut.c` — PDEClairaut (complete integral +
      singular envelope; reuses explicit verify/assemble)
- [x] 4. `src/calculus/dsolve_pdequasi.c` — PDEQuasilinear (Lagrange): semilinear
      (explicit) + conservation-law (implicit) classes
- [x] 5. Wired into `dsolve.c` (extern decls, is_pde cascade [quasi before Clairaut],
      init chain)
- [x] 6. Wired new sources into `tests/CMakeLists.txt` COMMON_SRC
- [x] 7. Unit tests in `tests/test_dsolve.c` (`t_pde_quasilinear`, `t_pde_clairaut`)
- [x] 8. Stress families in `tests/test_dsolve_stress.c`
      (`t_stress_pde_quasilinear` 10 members, `t_stress_pde_clairaut` 5 members)
- [x] 9. All 3 dsolve ctest suites green (182 / 35 / 9); `make check-c99` PASS
- [x] 10. valgrind: byte-for-byte identical to baseline (zero new leaks)
- [x] 11. Docs: `docs/spec/builtins/calculus.md` (2 method rows), weekly changelog,
      DSOLVE_PLAN.md §2a `[✓]` + M6 status

## Review

Implemented two first-order nonlinear PDE methods, closing most of DSOLVE_PLAN.md
§2a (Charpit remains the one deferred nonlinear item).

**`DSolve`PDEQuasilinear`** — Lagrange's method of characteristics for
`P u_v1 + Q u_v2 == R`. Semilinear class recurses the scalar ODE cascade for the
characteristic invariant, returns the explicit `x C[1][y/x]`-style general
solution; conservation-law class (Burgers) returns the implicit
`φ1 == C[1][u]` relation. Generalizes the constant-coefficient
`PDELinearFirstOrder`.

**`DSolve`PDEClairaut`** — the Clairaut form → complete integral
`C[1] v1 + C[2] v2 + f(C[1],C[2])` + singular envelope under
`IncludeSingularSolutions`.

**Substrate** — new `dsolve_run_pde_implicit` (implicit relation verified by the
implicit-function rule with concrete test functions; also routes explicit +
multi-branch cases), mirroring the ODE `dsolve_run_implicit`.

Worked examples verified end-to-end:
- `x u_x + y u_y == u` → `x C[1][y/x]`
- `u_x + x u_y == y` → (Y-along-characteristic path)
- `u u_x + u_y == 0` → `(y u - x)/u == C[1][u]` (Burgers, implicit)
- `u == x u_x + y u_y + u_x u_y` → `C[1]C[2] + C[1]x + C[2]y` (+ envelope `-x y`)

All correct, all declines clean, zero new memory leaks.

Deferred (documented): PDE Charpit; genuinely-quasilinear non-conservation
(`P/Q` depends on `u` with `R ≠ 0`); non-elementary characteristic integrals.
