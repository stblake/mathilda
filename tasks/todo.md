# Compile engine M3 — arrays / NDArray (increment M3a: rank-1 vectors)

Goal: machine arrays as first-class Compile values, delegating to the existing
NDArray infra (ndkernels/ndreduce/ndstruct/BLAS). Big milestone → ship in slices.
**M3a = rank-1 machine-real/complex vectors** (proves the architecture:
array-typed registers, NDArray ownership in the VM, EWKERNEL delegation,
broadcast, reduction, array args + result). Matrices/DOT/MATMUL/Part = M3b.

## Type representation
- Extend `CompileType` with array types. Encode (elem, rank) compactly:
  `CT_ARR` base; `CT_IS_ARRAY(t)`, `CT_ELEM(t)`, `CT_RANK(t)`, `CT_ARRAY(elem,rank)`.
  M3a uses rank-1 only (CT_ARRAY(CT_REAL,1), CT_ARRAY(CT_COMPLEX,1)).
- Arg spec `{v, _Real, 1}` → rank-1 real array arg. (extend compiled_function
  argspec parser later; for now the internal compile_expr accepts array arg_types.)

## Value / lifetime model (the crux)
- Array register uses `Slot.p` = owned `NDArrayData*` (args: borrowed).
- Frame is REUSED across calls → NEVER free-on-overwrite (stale ptr). Instead:
  emit explicit `OP_ARR_FREE dst` driven by the compile-time temp-stack
  discipline — when free_if_tmp pops an ARRAY temp, emit ARR_FREE. Args borrowed
  (never freed). Result array's ownership transfers to caller (not freed).
  Invariant: every allocated array temp is paired with ARR_FREE or is the result.

## Opcodes (finalize against ndkernels/ndreduce APIs — research pending)
- Elementwise vec⊕vec: delegate to ndkernels binary kernel if Plus/Times/... have
  one (OP_VKERN2 carrying fn ptr), else dedicated OP_VADD_R/C etc. Shape-checked
  at runtime (bail/NaN on mismatch).
- Broadcast scalar⊕vec: OP_VBCAST_* (or fold into VKERN2 with a splat).
- Unary fn over buffer: OP_VKERN (Sin/Exp/... via ndkernels vectorized path).
- Reduction: OP_VTOTAL_R/C (ndreduce Total) vec→scalar.
- OP_VLEN vec→int. OP_ARR_FREE.

## Public API
- `compiled_eval_array` or extend compiled_eval to accept/return NDArray handles
  via CompileValue (add an array case to the boxed union). Caller passes NDArray*
  args, receives an owned NDArray* result (or scalar).
- User surface (Compile[{{v,_Real,1}}, ...]) wiring = later; M3a proves the engine
  via test_compile.c with direct compile_expr + NDArray build/read.

## Steps
- [ ] Research ND APIs (agent) → finalize opcodes.
- [ ] CompileType array encoding + helpers.
- [ ] Slot/array-register plumbing; ARR_FREE + temp-stack frees.
- [ ] emit: array arg resolution, elementwise, broadcast, unary kernel, Total.
- [ ] VM cases (delegate to ndkernels/ndreduce; allocate/free NDArray buffers).
- [ ] Public eval API for array in/out.
- [ ] test_compile.c: parity vs interpreter (build NDArray, compile v→..., compare
      buffers); ownership (ASan/leaks clean, repeated calls reuse frame).
- [ ] Docs + changelog + memory.

## Deferred (M3b+)
- Rank-2 matrices, DOT/MATMUL (BLAS), Part/Slice, MAKEARR (pack tuple), Int
  arrays, List-of-machine-numbers coercion at the boundary, user Compile[] array
  argspec, Map/Table fusion.
