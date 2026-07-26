# NDSolve: numeric RHS compiler (evaluator-free nonlinear stepping)

## Goal
Compile the nonlinear reduced RHS f_i(t, Y) into fast numeric bytecode so the
time-stepper evaluates it WITHOUT the symbolic evaluator (no symbol binding, no
expr copy, no numericalize). Big win for explicit nonlinear (shallow-water,
Burgers) and nonlinear stiff (Newton residual + FD Jacobian).

Requirements (user): extensive unit tests, highly efficient, no memory leaks.

## Design
- New module `ndsolve_compile.{c,h}`: a stack-machine compiler + VM.
  - `nd_compile_rhs(P)` -> NdCompiled* or NULL (graceful bail on any unsupported
    construct -> existing symbolic sampler remains the fallback).
  - `nd_compiled_eval(C,t,Y,out)` -> out[0..d-1] via the VM.
  - `nd_compiled_jacobian(C,t,Y,Jout)` -> sparse colored FD (CPR) using the
    per-component variable dependencies the compiler already records.
  - `nd_compiled_free`.
- Opcodes: CONST/VAR/TVAR, ADD/SUB/MUL/DIV/NEG/INV, POW/POWI, and unary
  elementary fns (Sqrt/Exp/Log/Sin/Cos/Tan/.../Abs/Sign/Erf...), MAX/MIN/ATAN2.
- Leaves: numeric (nd_to_double) -> CONST; state sym NDSolve`w<k> -> VAR k;
  tvar -> TVAR; Pi/E/EulerGamma/... -> CONST; anything else -> bail.
- NdProblem gains `NdCompiled* compiled; bool compile_failed;` (zero-init).
  nd_rhs_real lazily compiles on first non-operator call; nd_jacobian_real
  prefers the compiled colored FD when available.

## Stages
- [ ] S1: compiler + VM; standalone unit test vs evaluate() on a big battery.
- [ ] S2: wire nd_rhs_real (lazy) + free at all teardown sites; end-to-end
      nonlinear NDSolve correctness (Burgers, shallow-water, nonlinear RD).
- [ ] S3: colored sparse FD Jacobian; prefer it in nd_jacobian_real.
- [ ] S4: valgrind (no leaks), benchmarks, docs, commit.

## Review — DONE (2026-07-26)

New module `ndsolve_compile.{c,h}`: stack-machine compiler + VM for the nonlinear
reduced RHS. `nd_rhs_real` lazily compiles on first non-operator call and runs
bytecode (no symbol binding / expr copy / numericalize); bails to the symbolic
sampler on any unsupported construct or when an EvaluationMonitor is attached.
`nd_jacobian_real` uses a CPR colored finite-difference Jacobian over the
bytecode (O(bandwidth) evals) via the per-component variable dependencies the
compiler records. Freed at all NdProblem teardown sites (mol ×3, ODE
nd_problem_free); NdProblem gained `compiled`/`compile_failed` (zero-init).

**Tests** (`tests/test_ndsolve_compile.c`, 6 groups): compiled eval vs the
symbolic evaluator on a broad arithmetic/elementary battery + nonlinear-PDE
couplings + a tridiagonal system — matched to **machine precision** (<1e-15);
colored-FD Jacobian vs analytic (symbolic-D) to ~1e-9; banded coloring gives 3
colors for 24 vars; and three graceful-bail cases. All NDSolve suites stay green
(172 PDE + 27 ODE + 87 classical + 6 compile).

**Measured** (compiler ON vs OFF): shallow-water dam-break (explicit) 1.83s ->
0.22s (~8.3x); porous-medium nonlinear diffusion (BDF stiff) >120s -> 0.14s
(>800x, the colored-FD Jacobian replaces the pathological symbolic per-entry
Jacobian). Valgrind: production compiler leak-free (main binary = macOS
baseline; zero compiler stacks); test-harness evaluate()-input leaks fixed with
eval_and_free.
