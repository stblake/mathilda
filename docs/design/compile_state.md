# Compiler build — state & handoff

Snapshot for resuming the `Compile[]` numeric-compiler work with fresh context.
Companion to [`compile.md`](compile.md) (the full design) and the memory files
`project_compile_engine`, `project_autocompile_numeric_builtins`,
`project_compile_engine` (read those too).

_Last updated: 2026-07-27. All work below is committed and pushed to `main`._

---

## 1. What exists and works today

The engine lives in **`src/compile/`** — a typed register-machine bytecode VM that
evaluates numeric expressions over machine numbers with no `Expr` allocation and
no runtime type dispatch (the opcode carries the type).

- **`compile.c`** — the VM + the emitter (`compile_expr`). Scalar type lattice
  `CT_BOOL/CT_INT/CT_REAL/CT_COMPLEX`; bottom-up inference; widening coercions;
  monomorphic typed opcodes; stack-discipline register allocation; reusable
  per-program `frame` (no per-call malloc); computed-goto dispatch on GCC/Clang.
  Coverage: full arithmetic, comparisons, boolean, elementary + special-function
  kernels (via the shared `ndkernels` registry), `If`, `Sum`/`Product`,
  `With`/`Module` locals, `Set`/`AddTo`/…/`Increment`, `CompoundExpression`,
  `Do`/`While`/`For`, `Nest`.
- **`compiled_function.{c,h}`** — user-facing `Compile[argspec, body]` →
  `CompiledFunction` object (new `EXPR_COMPILED` atom, refcounted immutable
  payload). Numeric args run the bytecode; symbolic args / uncompilable bodies
  fall back to the interpreter. `HoldAll | Protected`.
- **`autocompile.{c,h}`** — the adapter that lets the numeric builtins compile a
  held body once and evaluate it over machine numbers, with per-point interpreter
  fallback. Header is deliberately self-contained (NO `<complex.h>`/`compile.h`
  include — `double _Complex` is a builtin type — so it never leaks the `I`
  macro into callers).

### Auto-compile wiring (all shipped this session)

| Builtin | Chokepoint | Scope | Speedup |
|---|---|---|---|
| Plot / Plot3D | `plot_eval_fn`, `plot3d_eval_z/_fn` | real result; non-real → exclude point | ~215× / ~11× |
| Table | `table.c` numeric-range loop | **inexact iterator only**; exact untouched | ~128× |
| NIntegrate 1-D | `ni_eval_at` | finite/half/whole-line machine; complex fallback | ~353× |
| NIntegrate multi-D | `ni_mc_sample` | cubature + Monte-Carlo | ~504× |
| FindRoot scalar | `fr_eval_with_bindings` (pointer-identity `main_f`) | machine real Newton/Secant/Brent + FD | ~19× |
| FindRoot systems | `fr_run_newton_system_real` | per-component + Jacobian programs | ~6.9× |

**Shared correctness rule:** a compiled REAL program returns non-finite exactly
where the interpreter would produce a COMPLEX value (`Sqrt` of a negative).
Callers that exclude non-real points (Plot) drop it; callers needing the
contribution (NIntegrate/FindRoot/real-Table) fall back to the interpreter at
that one point. MPFR paths are UNTOUCHED everywhere. Uncompilable bodies (e.g.
`Zeta`) → `NULL` program → interpreter.

### Tests
`tests/test_compile.c` (engine, 63 checks), `tests/test_compiledfunction.c`
(user `Compile[]`), `tests/test_autocompile.c` (all 6 builtin wirings, parity +
fallback + oscillatory-regression + systems). All pass; `leaks`-clean. The only
red suite is `simplify_tests` (1 pre-existing radical failure, unrelated —
verified against HEAD).

### Milestones: M0, M1a, M1b, M2 (a/b/c), M4 (auto-compile) DONE.

---

## 2. Build & disk discipline (IMPORTANT — read before rebuilding)

- Main binary: `make -j4` (produces `./Mathilda`). ~53 MB objects + 7 MB binary.
- Tests: `tests/CMakeLists.txt` now compiles all `COMMON_SRC` **once** into an
  `OBJECT` library `mathilda_common`, spliced into every test target via
  `$<TARGET_OBJECTS:mathilda_common>` (commit `9ec6ad9`). Build a target:
  `cd tests/build && cmake --build . --target mathilda_common <target> -j4`.
- **Disk gotcha (bit us hard this session):** editing a `COMMON_SRC` file forces a
  `mathilda_common` rebuild + relink of every built test exe (~6 MB each). Dozens
  of build/test cycles = many GB of *write churn* (net small, but fills activity
  monitors). Prefer incremental builds; build only the specific target you need.
- **The real disk killer was `.claude/worktrees/`** — stale git worktrees from
  prior `isolation: worktree` subagents, each holding a full pre-OBJECT-lib
  `tests/build` (~1.3 GB). Reclaimed 7.5 GB by deleting their build artifacts
  (source + uncommitted changes preserved). If disk fills again, check
  `du -sh .claude/worktrees` and `git worktree list` first. Those 6 worktrees
  still exist (~456 MB, uncommitted risch/poly/core edits) pending the user's
  decision to `git worktree remove --force`.

---

## 3. NEXT: M3 — arrays / NDArray in the compiler

The engine is scalar-only. M3 adds machine **arrays** as a first-class value
category, delegating to the existing NDArray infra. Ship in slices; **M3a =
rank-1 machine-real/complex vectors** (proves the architecture end-to-end).

### 3.1 Design decisions (made, ready to implement)

**Type representation.** Extend `CompileType` (currently the enum
`CT_BOOL/INT/REAL/COMPLEX`) to encode array types compactly — `CT_ARR` base with
`CT_IS_ARRAY(t)/CT_ELEM(t)/CT_RANK(t)/CT_ARRAY(elem,rank)` helper macros. This
keeps the single-int flow through `infer_type` (returns `CompileType`), `Val.type`,
`Ctx.scope[].type`, `Ctx.arg_types[]`, `num_common`, `coerce`. M3a uses rank-1
only: `CT_ARRAY(CT_REAL,1)`, `CT_ARRAY(CT_COMPLEX,1)`. Matrices (rank-2) = M3b.

**Value slot.** `Slot` (`compile.c:29`) is `union { long long i; double r; double
_Complex z; const void* p; }` — the `p` field already exists; an array register
holds an owned `NDArrayData*`/`Expr*` there.

**Lifetime — THE crux (frame is reused across calls, so NEVER free-on-overwrite;
that would touch a stale pointer from a prior call).** Instead emit explicit
`OP_ARR_FREE dst` driven by the existing compile-time temp-stack discipline: when
`free_if_tmp` pops an ARRAY-typed temp, emit `ARR_FREE`. Array **args are
borrowed** (never freed). The **result array's** ownership transfers to the
caller (not freed). Invariant: every allocated array temp is paired with an
`ARR_FREE` or is the result. This mirrors malloc/free and needs no runtime
liveness pass. (Design doc §13 mentions an arena/free-list optimization — defer.)

**Result / args API.** Add an array case to `CompileValue` (boxed) or a new
`compiled_eval_array` entry point: caller passes `NDArray*`/`Expr*` array args,
receives an owned array (or scalar). For M3a the *user surface*
(`Compile[{{v,_Real,1}}, …]` argspec parsing in `compiled_function.c`) can come
last — prove the engine first with direct `compile_expr` + NDArray build/read in
`test_compile.c`.

### 3.2 ND delegation API (researched — these are the functions to call)

Kernel structs & helpers live in **`src/ndarray.h`** (there is NO `ndkernels.h`).
Buffer = row-major flat, complex interleaved (re,im), NOT C99 `_Complex`.

- **Elementwise arith / broadcast:** `ndarray_elementwise(Expr** args, size_t n,
  bool is_plus)` (`ndarray.h:127`) — Plus/Times over flat buffers with numpy
  scalar broadcast; new EXPR_NDARRAY or NULL. Power: `ndarray_elementwise_power`,
  `ndarray_scalar_power`, `ndarray_base_scalar_power` (`ndarray.h:134-145`).
- **Unary/binary kernel map:** `ndarray_map_unary(const Expr* a, const
  NDUnaryKernel* k)` / `ndarray_map_binary(a0,a1,k)` (`ndarray.h:185-186`) — loop
  the whole buffer for you, allocate a new array, multithreaded. Look up the
  kernel via `symtab_lookup(head)->ndarray_unary_kernel` /
  `->ndarray_binary_kernel` (`symtab.h:88-89`). Kernel structs `NDUnaryKernel`
  (`ndarray.h:166`, fns cplx/real/real_closed/to_real), `NDBinaryKernel`
  (`ndarray.h:179`).
- **Reductions:** `ndred_total/mean/max/min/...(Expr* res)` (`ndreduce.h:35-48`) —
  reduce the LEADING axis; take the WHOLE call Expr (borrowed). For a bare buffer,
  wrap it in an EXPR_NDARRAY + call node, or use raw helpers in
  `ndarray_internal.h` (`nd_parallel_reduce`, `nd_gather_real`, …).
- **Dot/MatMul:** `ndarray_dot2(const Expr* a, const Expr* b, bool* shape_error)`
  (`ndarray.h:118`); BLAS `dgemm` wrapper `dot2(Expr*,Expr*,bool*)`
  (`src/linalg/dot.c:60`). — M3b.
- **Part/structural:** `ndarray_part(a, indices, n, &degrade)` (`ndarray.h:99`);
  `ndstruct_sort/reverse/transpose/...(Expr* res)` (`ndstruct.h`). — M3b.
- **Build/read buffers:** construct with `expr_new_ndarray(rank, dims, buf,
  dtype)` (adopts `buf`, copies `dims`; `expr.h:164`). Element access with
  `ndt_get/ndt_set(buf,k,dtype,&re,&im)` (`ndarray.h`). `na_load_vector` /
  `na_build_vector` / `na_read_scalar` / `na_scalar` (`src/linalg/numarray.h`) move
  between Expr and raw `double*` buffers. Free a built array with `expr_free`.
  **All ND fast paths allocate a new buffer and never mutate inputs — none are
  in-place.** `NDType`: `NDT_FLOAT64=0` (default), `FLOAT32`, `COMPLEX64`,
  `COMPLEX32` (`expr.h:45`).

### 3.3 M3a opcode sketch (finalize during impl)

- Elementwise vec⊕vec `+ - * /`: either dedicated `OP_VADD_R/C` etc. or a generic
  `OP_VKERN2` carrying the binary-kernel fn ptr; simplest to call
  `ndarray_elementwise`/`ndarray_map_binary` from the VM op. Shape-checked at
  runtime (bail → the whole `compiled_eval` fails, caller falls back).
- Scalar⊕vec broadcast: `ndarray_elementwise` already broadcasts; emit a splat or
  a dedicated `OP_VBCAST_*`.
- Unary fn over buffer `Sin[v]`, `Exp[v]`: `OP_VKERN` carrying the unary-kernel ptr
  → `ndarray_map_unary`.
- Reduction `Total[v]` (vec→scalar): `OP_VTOTAL_R/C` → `ndred_total` (wrap buffer).
- `OP_ARR_FREE dst`, `OP_VLEN` (vec→int).

### 3.4 M3a step list

1. `CompileType` array encoding + helper macros.
2. Array-register plumbing in `Slot`/`Val`; `OP_ARR_FREE` + emit frees on
   `free_if_tmp` of an array temp; array-arg loading (borrowed).
3. `emit`: array-arg symbol resolution; elementwise; broadcast; unary kernel;
   `Total`. `infer_type` array cases.
4. VM cases delegating to `ndarray_elementwise`/`ndarray_map_unary`/`ndred_total`
   (allocate/free NDArray buffers; write handle into `Slot.p`).
5. Public array in/out eval API.
6. `test_compile.c`: build an NDArray, compile `v → v+1` / `2 v` / `Sin[v]` /
   `Total[v]`, compare buffers to interpreter; ownership (leaks + ASan clean,
   repeated calls reuse the frame with no leak/UAF).
7. Docs + changelog (`docs/spec/changelog/2026-07-27.md`) + memory.

**Deferred to M3b+:** rank-2 matrices, `DOT`/`MATMUL` (BLAS), `Part`/`Slice`,
`MAKEARR` (pack a fixed tuple), Int arrays, List-of-machine-numbers coercion at
the boundary, user `Compile[]` array argspec, Map/Table fusion.

### 3.5 Watch-outs

- Frame reuse ⇒ never free-on-overwrite; only the explicit `ARR_FREE`/result model.
- `ndreduce`/`ndstruct` take the whole *call* Expr (borrow-only); the
  map/elementwise/dot helpers take *operand* Exprs. Wrap a bare buffer accordingly.
- Verify the `CompiledFunction` `expr_free` path frees any array-returning
  program's result correctly (ownership handoff), like `InterpolatingFunction`.
- Gate leaks with BOTH `leaks` and an ASan build (ASan caught a latent OOB the
  memcheck missed earlier in this project).
