# DSolve — OperatorFactor / DSolve`DFactor (§1c linear-operator factoring)

Factor a linear operator L by finding a first-order right factor (D-r), r∈C(x)
(hyperexponential solution Exp[∫r]); peel via operator right-division, recurse
DSolve on the order-(n-1) quotient, close with a trailing first-order solve.
- DSolve`OperatorFactor: cascade method, homogeneous, order ≥ 3 (Kovacic owns 2).
- DSolve`DFactor[eqn,y,x]: returns {Dx - r1, Dx - r2, ...} (inert Dx = d/dx).

Self-contained new file (NO changes to dsolve_kovacic.c → zero regression risk).

Note (learned in design): distinct-pole factors → dilog in reductions (decline);
shifted-Euler (shared pole x-b, b≠0) is elementary AND unique to OperatorFactor
(EulerCauchy only detects pole at 0). Use shifted-Euler + all-constant in tests.

## Implementation
- [x] src/calculus/dsolve_operator_factor.c — of_build_ansatz, of_riccati_residual
      (Bell polys), of_find_factor (branch-iterate + remainder guard), of_divide
      (binomial recurrence), trailing solve (measured C[kk]); dsolve_operfactor_try
      + builtin + init; builtin_dsolve_dfactor + recursive collector
- [x] src/calculus/dsolve.c — extern, DS_OPERFACTOR enum, "OperatorFactor" map,
      cascade after kovacic, pinned case, init (also registers DSolve`DFactor)
- [x] tests/CMakeLists.txt — add dsolve_operator_factor.c

## Tests
- [x] tests/test_dsolve.c — 6 unit tests (method/more/ivp/declines/dfactor/auto)
- [x] tests/test_dsolve_stress.c — operfactor_ok + dfactor_ok generators

## Docs
- [x] docs/spec/builtins/calculus.md — OperatorFactor + DFactor
- [x] docs/spec/changelog/2026-08-31.md — ## DSolve — OperatorFactor / DFactor
- [x] DSOLVE_PLAN.md — §1c [✓] OperatorFactor; M5 note

## Gates
- [x] make -j clean (no warnings; fixed a gcc-16 -Warray-bounds false positive
      by pre-sizing the ansatz arrays instead of realloc)
- [x] REPL spot-checks — shifted-Euler/constant/resonant/order-4 solve; distinct-pole
      /Airy-3/order-2 decline cleanly; DFactor correct; AUTO dispatch fires; regressions ok
- [x] ctest -R dsolve — 3/3 green (Kovacic/M5 intact = zero regression)
- [x] make check-c99 — exit 0
- [x] valgrind — OperatorFactor AND DFactor paths: 1x==6x==13,440+6,312 B (no per-call leak)
- [ ] commit + push (awaiting user go)

## Review — DONE

`DSolve`OperatorFactor` (§1c) solves homogeneous linear ODEs of order ≥ 3 by
factoring the operator: find a first-order right factor `(D-r)`, `r∈C(x)` (Bell-
polynomial Riccati `Σ a_k P_k(r)==0`, undetermined-coefficient search with per-branch
remainder guard), peel via operator right-division (binomial recurrence), recurse
`DSolve` on the order-(n-1) quotient, close with the trailing first-order solve
(constant `C[kk]` measured contiguous). `DSolve`DFactor` returns `{Dx-r1,...}`.
New self-contained file `src/calculus/dsolve_operator_factor.c` (no Kovacic changes).

**Verified:** shifted-Euler flagship (unique niche — EulerCauchy only detects pole 0),
constant order 3/4, resonant repeated, order-4 all solve + back-substitute; distinct-
pole (dilog), Airy-3, order-2 all decline cleanly; DFactor reconstructs the operator;
AUTO dispatch fires; Kovacic/Euler/Airy/constcoeff regressions intact.

**Two bugs caught during testing (both fixed):** (1) gcc-16 -Warray-bounds false
positive on the realloc+macro ansatz builder → rewrote to pre-size exactly; (2) the
DFactor reconstruct test used `Reverse[fs]` which only passed for commuting (constant)
factors — the correct order folds over `fs` directly (innermost factor applied first),
and the reference operator-on-tf must be built by the same composition, not by
`op /. y->Function[tf]` (which mis-evaluates the derivative terms).

**Follow-ups (out of scope):** inhomogeneous forcing; irregular-singular (double-pole)
r; 2nd-order right factors (Beke); distinct-pole compositions (dilog reductions);
Sturm–Liouville EigenvalueProblem (§1f, the last open §1c/1f item).
