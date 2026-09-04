# DSolve M11 — Finish Phase 1 ODE gaps

Plan: `/Users/user/.claude/plans/let-s-continue-the-dsolve-parsed-heron.md`

## Phase 1 — Refactor dsolve_linsys.c ✅ DONE
- [x] New `src/calculus/dsolve_linsys.h` (tidy / matexp / extract_Ab / assemble + fwd decls)
- [x] `mat_exp` → `dsolve_linsys_matexp`; `tidy` → `dsolve_linsys_tidy` (skips ComplexExpand on real bodies)
- [x] `dsolve_linsys_extract_Ab` (x-guards moved to caller); `dsolve_linsys_assemble(M,t,xvar,b,b_zero,n)`
- [x] Rewrite `dsolve_linsys_solve`; `dsolve.c` `#include "dsolve_linsys.h"`
- [x] Gate: constant systems (diag/defective/complex/forced/singular/IVP/pinned) verified identical

## Phase 2 — LinearSystemVarCoeff ✅ DONE
- [x] `dsolve_linsys_varcoeff.c` (scalar-factor A=f(x)·B, τ=∫f, reuse assemble); COMMON_SRC; cascade; pinned; init

## Phase 3 — Linear BVP hardening ✅ DONE
- [x] `dsolve_fit_constants`/`dsolve_fit_system` no_solution out-param on empty Solve
- [x] `dsolve_run`/`dsolve_run_system` thread → empty List{}; header sig update
- [x] over-det→{}, under-det→free C kept, IVP/undecided unaffected (verified)

## Phase 4 — DSolve`EigenvalueProblem ✅ DONE
- [x] `dsolve_eigenvalue.c`: detect eigenparameter+homogeneous BCs, analytic family, verify under n∈Integers
- [x] TrigReduce+Simplify verifier (C[1]→bare symbol); pinned; COMMON_SRC; init; DD/NN/mixed all solve

## Phase 5 — Tests ✅ DONE
- [x] varcoeff unit (t_sys_varcoeff_*) + stress (t_stress_linsys_varcoeff)
- [x] BVP unit (t_bvp_overdetermined/underdetermined/system/undecided)
- [x] Eigenvalue unit (t_eig_dirichlet/neumann/mixed/no_misfire)

## Gates
- [x] ctest: dsolve_tests (170), dsolve_stress_tests (33), dsolve_m5_stress_tests (9) all green
- [x] `make check-c99` clean; main build clean; valgrind leak-clean (no new leaks; total at baseline 13,440/6,312)
- [x] SymPy cross-validate (python3.11): c1 spans {x,x³}, c2 oscillatory — agree
- [x] Docs: DSOLVE_PLAN.md (§1e/§1f + M11), docs/spec/builtins/calculus.md, changelog 2026-08-31.md; memory

## Review — DONE (2026-09-04)

Delivered all three Phase-1 ODE gaps, each verified + tested:
1. **Refactor** `dsolve_linsys.{c,h}` — shared `matexp`/`tidy`/`extract_Ab`/`assemble`;
   `tidy` now skips ComplexExpand on real bodies (avoids Log→Abs/Arg leak). LinearFirstOrderSystem unchanged.
2. **LinearSystemVarCoeff** — scalar-factor A=f(x)B via τ=∫f; solves the coupled var-coeff class that declined before.
3. **BVP soundness** — inconsistent BVP → {} (was: silent unfitted general). Empty-Solve signal threaded via no_solution.
4. **EigenvalueProblem** — pinned Sturm-Liouville first cut (DD/NN/mixed), analytic family + back-sub verify under n∈Integers.

Key gotchas (→ memory): Simplify integer assumption needs a bare symbol not C[k];
TrigReduce before Simplify for factored Cos[Pi(n-1/2)]; ComplexExpand splits Log[x] on real vars.

Deferred (documented future): commutative-antiderivative/Floquet systems; non-constant-weight
Sturm-Liouville, Robin/periodic BCs, λ==0 Neumann mode, DEigensystem/DEigenvalues surface.
