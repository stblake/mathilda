# M8 — general linear systems (defective / singular / arbitrary forcing) + triangular

Fixes: `DSolve[{y'[t]+1==1, x'[t]+y[t]==0}, {y[t],x[t]}, t]` returning unevaluated.
Root cause: coupled constant-coeff system with **defective + singular** matrix
`A={{0,0},{-1,0}}`; the eigen-only `LinearFirstOrderSystem` declines on the zero
eigenvector, and `DecoupleSystem` declines (eq 2 mentions two functions).

## Plan (approved: DSOLVE_PLAN.md M8)

- [x] 1. Reworked `dsolve_linsys.c` around `mat_exp` (Jordan fundamental matrix,
      finite nilpotent series, var-params forcing, Cosh+Sinh realifier).
- [x] 2. Shared `dsolve_renumber_constants` / `dsolve_extract_system_body` in
      `dsolve_common.c`; `dsolve_decouple.c` refactored onto them.
- [x] 3. New `dsolve_triangular.c` — `TriangularSystem` (DAG forward substitution,
      private-head constant parking).
- [x] 4. Wired into `dsolve.c` (extern + init + cascade decouple→triangular→linsys)
      and `tests/CMakeLists.txt` COMMON_SRC.
- [x] 5. Build clean, `make check-c99` clean, `dsolve_tests` green (28 tests, 5 new),
      valgrind: no definitely-lost delta vs scalar-DSolve baseline.
- [x] 6. Docs: `calculus.md`, changelog `2026-08-31.md`, DSOLVE_PLAN M8 done; lessons +
      3 memories.

## Review

Closed the reported `DSolve[{y'[t]+1==1, x'[t]+y[t]==0}, {y[t],x[t]}, t]` decline
(defective + singular A) in full generality, not with a special case:

- **`LinearFirstOrderSystem`** now solves any constant-`A` system via the Jordan
  fundamental matrix — diagonalizable, defective (`x^k e^{λx}`), complex
  (`e^{αx}Cos/Sin`), and singular-with-forcing (variation of parameters). Verified
  by back-substitution across all spectral cases.
- **`TriangularSystem`** (new) closes coupled-triangular systems at variable
  coefficient too; the two together leave only non-triangular variable-`A` (future).
- One real bug found+fixed mid-implementation: forward-substituted constants
  collided with the scalar engine's fresh `C[k]` (`x -> C[2]-C[2]t`); fixed by
  parking solved constants in the private head `DSolve`sysK`. `GeneratedParameters`
  does not fix this (it renames all `C[k]`). Caught by residual verification.

No regressions: scalar DSolve, decoupled systems, and complex-eigenvalue systems
all unchanged; genuinely out-of-scope systems decline gracefully (no hang/crash).
