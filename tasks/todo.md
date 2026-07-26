# Compile[] engine — M0 substrate (scalar core)

Per docs/design/compile.md. This increment: the reusable engine core, scalar
lattice only (Bool/Int/Real/Complex). Register VM designed for control flow.
Extensive unit + stress tests. Highly efficient VM (monomorphic typed opcodes,
no runtime type dispatch, Real hot path, stack-discipline register reuse).

## Scope (this pass)
- [ ] src/compile/compile.{h,c}: typed value slots, IR, register VM, front-end
      (type inference + lowering + bail), arithmetic/comparison/boolean/
      elementary opcodes, coercions.
- [ ] compiled_eval (boxed) + compiled_eval_real (fast all-real path) +
      compiled_arg_deps (for NDSolve reuse later).
- [ ] makefile SRC += src/compile/*.c ; tests/CMakeLists target.
- [ ] tests/test_compile.c: parity vs interpreter across all types/ops/coercions/
      nesting; stress (deep+wide expressions, many args, fuzz); perf smoke.
- [ ] valgrind + ASan clean.

## Deferred to later milestones (design doc)
- Generic KERNEL over shared ndkernels registry (all numeric special fns) — M1/M4.
- Control flow If/Do/For/While/Nest — M2.
- Arrays / NDArray / lists of machine numbers — M3.
- User-facing Compile[] builtin + CompiledFunction object — M1.
- NDSolve migration onto the engine — after core is proven.

## Review — DONE (2026-07-26)

Delivered the reusable engine core: `src/compile/compile.{h,c}`.
- Type lattice Bool/Int/Real/Complex; bottom-up inference; widening coercions
  inserted where operand types differ; bail (NULL) on anything outside the subset.
- Monomorphic typed opcodes (a Real add is one instruction; no runtime tag);
  register machine with stack-discipline temp allocation (O(depth) registers);
  reusable frame per program (no per-call malloc).
- Ops: arithmetic (I/R/C), Mod/Quotient, integer/real/complex Power, all the
  elementary functions (Sqrt/Exp/Log/trig/hyperbolic/inverse/Erf/…), Abs/Sign,
  Floor/Ceiling/Round (→Int), Re/Im/Arg/Conjugate, Max/Min, ArcTan[x,y],
  comparisons (→Bool), And/Or/Not/Xor, named constants.
- API: compile_expr / compiled_eval (boxed) / compiled_eval_real (fast all-real,
  no boxing) / compiled_arg_deps (sparsity, for NDSolve reuse) / compiled_free.

Tests: `tests/test_compile.c` (33 checks) — parity vs the interpreter to machine
precision (max_rel ~1e-16) across all types/ops/coercions/nesting; all-real fast
path == boxed path; arg-dep introspection; graceful bail; stress (400-deep nest,
500-term sum, 8 args); performance **~234× faster** than the interpreter
(~86 ns/call). Valgrind + ASan clean. Existing suites unaffected.

## Next (per design doc)
- M1: generic KERNEL over shared ndkernels (all special fns) + user Compile[]
  builtin + CompiledFunction object.
- M2: control flow (If/Do/For/While/Nest). M3: arrays/NDArray.
- Then migrate NDSolve onto the engine.
