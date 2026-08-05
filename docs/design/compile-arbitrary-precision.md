# `Compile[]` arbitrary precision — GMP bignum & MPFR real/complex in compiled code

Status: **implemented (P1–P2, P3 substrate); sampler wiring remaining**. Author: (this doc).

Implemented as of 2026-08-05: the managed-scalar engine (GMP bignum, MPFR real,
MPFR complex) behind `Compile[]`'s `WorkingPrecision -> n` and
`"BigIntegers" -> True` options (P1+P2), plus the precision-aware
auto-compilation entry point `autocompile_new_prec` / `autocompiled_eval_mpfr`
(the P3 substrate). Remaining: wiring the individual numeric samplers
(`NIntegrate`/`Plot`/`FindRoot`/`NDSolve`) to call that entry at high
`WorkingPrecision`, managed arrays, managed control flow, and the warm
thread-local container arena (a speed follow-up; §7 notes the v1 uses
correctness-first per-call allocation, which never touches the machine path).

Supersedes the non-goal in [`compile.md`](compile.md) §2 ("Arbitrary precision
(MPFR) inside compiled code — compiled code is machine precision by definition;
MPFR paths stay interpreted"). This document lifts precisely that deferral, and
only that one — every other property of the engine in `compile.md` /
[`compile_state.md`](compile_state.md) is preserved, including the machine path
unchanged to the byte.

---

## 1. Motivation

The engine in `compile.md` compiles a numeric body once into a monomorphic
bytecode program and runs it over raw machine numbers with **no `Expr`
allocation and no runtime type dispatch** — an `M`-milestone lineage that reached
8×–800× on NDSolve and 3.8× fused elementwise. Its scalar lattice is exactly four
machine types: `CT_BOOL`, `CT_INT` (int64), `CT_REAL` (float64), `CT_COMPLEX`
(`double _Complex`).

Everything outside that lattice is turned away at the boundary. In particular,
`literal()` (`src/compile/compile.c:1433-1435`) **downcasts** an `EXPR_BIGINT`
via `mpz_get_d` and an `EXPR_MPFR` via `mpfr_get_d` to a machine `double`, and an
int64 register that overflows at runtime raises `OP_FAIL` and the **whole body**
silently bails to the interpreter. Both are correct — the second is exactly why
the interpreter promotes to GMP — but both leave arbitrary precision *entirely
uncompiled*.

Meanwhile the interpreter's arbitrary-precision paths are the worst case for its
per-op overhead: every `mpz`/`mpfr` operation allocates a fresh `Expr`, runs the
full `evaluate()` fixed-point loop, reads attributes, and negotiates precision
(`numeric_combined_bits`) — then does one `mpz_add` / `mpfr_add` and frees. The
arithmetic is a rounding-error fraction of the surrounding tree-walk. **That
surrounding cost is exactly what a compiled loop removes.**

Two workloads want this:
- **High-precision numerics** (MPFR): tight iteration at tens-to-hundreds of
  digits — Newton/Halley, series partial sums (`exp`, Bessel, ζ), orbit maps
  (Lyapunov/bifurcation), and the sample bodies of `NIntegrate` / `Plot` /
  `FindRoot` / `NDSolve` at high `WorkingPrecision`, which today cannot compile
  at all and run fully interpreted.
- **Exact big-integer loops** (GMP): Fibonacci/factorial/binomial/Catalan,
  modular exponentiation, continued-fraction convergents, exact-integer DP.

## 2. Goals / non-goals

**Goals**
- Compile bignum, MPFR-real, and MPFR-complex **scalar** arithmetic so a tight
  loop runs with no `Expr` allocation and no `evaluate()` tree-walk.
- **Opt-in per `Compile[]`.** Default behaviour is byte-for-byte the current
  machine path.
- **Zero machine-precision regression** (§11) — the non-negotiable constraint,
  proven by an A/B timing gate, not by inspection.
- One fixed working precision per compiled program (§5), which makes the arena
  containers warm and sidesteps dynamic precision-contagion.
- Correctness parity with the interpreter's `N[body, n]` / exact result, to full
  precision, verified by tests. No leaks.

**Non-goals (deferred)**
- Bignum / MPFR **arrays** — no packed dtype exists (`NDType`, `expr.h:71-78`),
  and there is no SIMD/BLAS benefit; Mathematica never packs high-precision
  arrays either. Managed types are scalar-only here.
- **Managed-op CSE / constant-folding.** Managed registers are treated as
  opaque pointers by the optimiser in v1 (§10); value-numbering across managed
  ops is a later optimisation.
- LLVM/native codegen and MPC (there is no MPC linkage in the build; complex is a
  pair of `mpfr_t`, §8).
- **Auto-selecting** bignum. The compiler cannot know at compile time whether an
  int64 register overflows, so bignum is a user-requested type (`"BigIntegers"`),
  never inferred, and is not an auto-compilation feature.

## 3. What the interpreter already has (and we reuse)

- **GMP bignum** — `EXPR_BIGINT` holds an inline `mpz_t` (`expr.h:198`);
  `expr_bigint_normalize` (`expr.c:438`) demotes back to int64 when it fits.
- **MPFR real** — `EXPR_MPFR` (under `USE_MPFR`) carries its **own** precision in
  bits (`expr.h:206`); constructors `expr_new_mpfr_bits/_from_d/_from_mpz` in
  `expr.c`.
- **MPFR complex** — structural `Complex[EXPR_MPFR, EXPR_MPFR]`, and an
  inner-loop value type `ncpx { mpfr_t re, im; }` (`src/numeric_complex.h:137`)
  with kernels `ncpx_add/sub/mul/div/exp/log/sin/cos/sqrt/pow/...` at an explicit
  working precision — **no MPC**. Collapse-to-real when the imaginary part rounds
  to zero via `numeric_mpfr_make_complex`.
- **Precision machinery** — `numeric_digits_to_bits` / `numeric_bits_to_digits`
  (`numeric.c:48-56`).
- **Overflow discipline** — `checked_int.h` (`long long` `ci_*` family, confined
  to `src/compile/`) already detects int64 overflow → `OP_FAIL` → interpreter.
  Bignum compilation *replaces* that bail with a bignum register when the user
  asks; the machine path keeps the bail unchanged.

The design adds a compiled *substrate* over these existing kernels; it does not
re-implement any arithmetic. Each new opcode delegates to an `mpz_*` / `mpfr_*` /
`ncpx_*` routine, so it is correct by construction, exactly as the array opcodes
delegate to the NDArray layer (`compile.md` §8a).

## 4. User surface

Two options on `Compile[]`, parsed in `builtin_compile`
(`src/compile/compiled_function.c`) and threaded into the engine (§9):

```
Compile[{{x, _Real}}, body, WorkingPrecision -> 50]   (* MPFR reals @ 50 digits *)
Compile[{{x, _Real}}, body]                            (* MachinePrecision: unchanged *)
Compile[{{n, _Integer}}, body, "BigIntegers" -> True]  (* exact GMP integers      *)
```

- `WorkingPrecision -> n` — numeric `n` (decimal digits) makes real/complex
  registers MPFR at `numeric_digits_to_bits(n)` bits. `MachinePrecision` (the
  default) selects the current machine-only path, unchanged.
- `"BigIntegers" -> True` — integer-typed registers become GMP bignums (the
  overflow-bail is replaced by exact arithmetic). Default `False`. Orthogonal to
  `WorkingPrecision` (integers have no precision).

Auto-compilation (§12) does not use user syntax: the sampler passes its own
`WorkingPrecision` programmatically to the compile entry point.

`WorkingPrecision -> MachinePrecision` and no `"BigIntegers"` ⇒ the program has
**no managed register**, and every mechanism below is inert (§11).

## 5. The fixed-precision-per-program contract

Mathematica's `Compile` is machine-only, so arbitrary precision here is a defined
Mathilda extension and we choose the simplest contract that is still useful:

> A `CompiledProgram` carries a single `mpfr_prec_t`, fixed when the
> `CompiledFunction` is instantiated. Every MPFR/`ncpx` register and temporary in
> the program is that precision. Bignum registers are exact and carry no
> precision.

This deliberately drops Mathematica's dynamic precision-*contagion* (the
"tighter operand wins", `numeric_min_inexact_bits`) — that is a symbolic
evaluation concern. A compiled kernel is a fixed-precision numeric routine, which
is exactly what `NIntegrate`/`FindRoot`/`NDSolve` want at a chosen
`WorkingPrecision`, and it is what lets the arena preallocate every container to
one precision and keep them warm across calls (§7).

## 6. Type lattice extension

The lattice (`compile.h:37-64`) packs machine scalars into `0..3`, array types
into `[4,36)` (`CT_ARR + 4*(rank-1) + elem`), and associations into `[36,40)`.
Add three scalar kinds **above** that window:

```c
CT_BIGINT     = 40,   /* GMP mpz_t                       */
CT_BIGREAL    = 41,   /* MPFR mpfr_t                     */
CT_BIGCOMPLEX = 42,   /* ncpx = { mpfr_t re, im }        */
```

This is **safe by default**: every existing scalar gate — `t == CT_REAL`,
`t < CT_REAL`, `t <= CT_COMPLEX`, `CT_IS_ARRAY`, `CT_IS_ASSOC` — already *rejects*
a type ≥ 40, so each machine fast path bails on a managed type until explicitly
taught. Managed handling is purely additive; nothing in the machine path changes
meaning.

Widening (in `num_common` / `coerce`, `compile.c`): the join crosses the
machine→managed boundary — `CT_INT ⊔ CT_BIGINT = CT_BIGINT`,
`CT_REAL ⊔ CT_BIGREAL = CT_BIGREAL`, `CT_BIGREAL ⊔ CT_BIGCOMPLEX = CT_BIGCOMPLEX`,
and a mixed body promotes machine operands with an explicit coercion opcode
(`mpz_set_si` / `mpfr_set_d` / build an `ncpx` with zero imaginary). `finite_result`
(`compile.c:8107`) gains a `CT_BIGREAL`/`CT_BIGCOMPLEX` arm so a compiled MPFR
NaN bails like a machine NaN.

Type inference (`infer_type`, `compile.c:2450`) is only a routing hint; when
`WorkingPrecision` is numeric, real/complex argument declarations are typed
`CT_BIGREAL`/`CT_BIGCOMPLEX` at instantiation, and a wrong guess still ends in a
clean bail.

## 7. The managed register bank and the TLS precision-keyed arena (the crux)

The single most important decision. It is derived from a constraint the frame
model already documents (`compile.c:8121-8128`): the register frame *used to be*
program-owned and was therefore "neither reentrant nor thread-safe"; it was moved
to the per-call, naturally-nested **C stack** (`Slot stackframe[512]`,
`compile.c:6936/8167`), with array/tile temporaries living in that same per-call
frame. A managed store attached to the `CompiledProgram` would reintroduce
exactly the eliminated bug — it would race across `nd_parallel_for` workers and
corrupt across reentrant `vm_call`s.

But an `mpz_t`/`mpfr_t` is **not** POD: it owns a heap limb buffer and needs
`init`/`clear`. It cannot live for free in the flat `Slot` frame (a `Slot` is a
16-byte union; the frame has no per-slot setup). So:

**A managed `Slot` holds an index into a thread-local, precision-keyed container
arena** — not a value, not a program-global pointer. The arena:
- is a `VM_TLS` structure mirroring `vm_call_depth` (`compile.c:8157`): each
  thread has its own, so worker threads never share containers;
- is a **stack** whose pointer advances by the program's managed-register count
  on `vm_call` / `compiled_eval` entry and retreats on exit — so it **nests with
  the C call stack** and is reentrant by construction;
- **grows on demand and never shrinks**; containers are `mpz_init` /
  `mpfr_init2(_, prec)` on first use and `mpz_clear`/`mpfr_clear` only at thread
  exit — so after warmup there is **no `malloc` per call**, only an index bump
  and (when a reused slot's precision differs) an `mpfr_set_prec`. This is what
  turns the sampler case (one `compiled_eval` per sample point) from
  malloc-bound into arithmetic-bound.

Register-bank layout gains a fourth boundary. Today the frame is
`[0, arr_base)` scalar, `[arr_base, tile_base)` array, `[tile_base, nreg)` tile
(`CompiledProgram`, `compile_internal.h:498-500`; `patch_reg`, `compile.c:8009`).
Add `managed_base` after `tile_base`; managed registers are `[managed_base, nreg)`
and their `Slot.i` is the arena index. Sizing is static and identical to the
array/tile banks — the register allocator is already stack-disciplined
(O(expression depth)), so the peak count of simultaneously-live managed temps is
known at compile time.

**For a machine-only program `managed_base == nreg`: the arena is never entered,
never sized, never touched.** (§11.)

Invariant, asserted at emit time: **no managed register is read or written inside
an `OP_APAR` (parallel strip-loop) range.** `par_chunk` (`compile.c:6944`)
`memcpy`s the parent frame's slots to each worker; a copied arena index would
alias one thread's container into another. This holds trivially in v1 (managed
types are scalar-only and never appear in a fused/parallel array loop), but the
assertion is what keeps a future array milestone from silently violating it.

## 8. Opcodes and kernel dispatch

The opcode set is generated from one `OPLIST` X-macro
(`compile_internal.h:137-341`) that drives the enum, the VM computed-goto jump
table, the optimiser property table, and the disassembler. Adding opcodes there
is O(1) for every *existing* opcode (the jump table is indexed by opcode number),
so this does not touch the machine hot path (§11).

Keep the count small by reusing the **function-pointer immediate** mechanism
(`Slot.p`) that special functions already use (`try_kernel`, `compile.c:3217`):

| Group | Opcodes | Dispatch |
|-------|---------|----------|
| MPFR real | `OP_R_UN`, `OP_R_BIN` | raw `mpfr_*` fn ptr (`mpfr_add`, `mpfr_sin`, …) writing the preallocated dst container at program precision |
| MPFR complex | `OP_Z_UN`, `OP_Z_BIN` | `ncpx_*` kernels at program precision |
| Bignum | `OP_ZI_ADD/SUB/MUL/POWI/NEG` | `mpz_*` |
| Comparisons → `CT_BOOL` | `OP_R_CMP`, `OP_ZI_CMP` | `mpfr_cmp` / `mpz_cmp` |
| Coercions | `OP_INT_TO_ZI`, `OP_REAL_TO_R`, `OP_R_TO_Z` | `mpz_set_si`, `mpfr_set_d`, build `ncpx` |
| Move | `OP_MANAGED_MOVE` | **deep value copy** `mpfr_set`/`mpz_set`/`ncpx_set` (never a `Slot`/pointer copy — see §10) |
| Literal | materialised into a managed const register at program precision | replaces the `literal()` downcast *only* when a managed type is in play |

MPFR's ABI fits this cleanly: `int mpfr_add(mpfr_t, mpfr_srcptr, mpfr_srcptr,
mpfr_rnd_t)` and `int mpfr_sin(mpfr_t, mpfr_srcptr, mpfr_rnd_t)` — one generic
binary and one generic unary handler cover every registered `mpfr_*`. In-place
aliasing (`dst == a`) is legal for `mpfr_*`/`mpz_*` and is already the shape the
optimiser tolerates (§10).

All MPFR/`ncpx` code is under `#ifdef USE_MPFR`; bignum uses GMP, which is always
present. Watch the `int64_t`-vs-`long long` split (`checked_int.h`): the compiler
already uses the `long long` `ci_*` family, and the new `mpz`/`mpfr` glue must not
cross into the `_i64` family.

## 9. Boxing and calling convention

`CompileValue.v` already carries an `Expr* a` arm (`compile.h:71`), so **no struct
grows**: a boxed `EXPR_BIGINT` / `EXPR_MPFR` / `Complex[MPFR,MPFR]` argument rides
in `v.a`.

- **In** — `load_arg` (`compile.c:8079`, currently `return false` on an unknown
  `v->type`, a safe default) gains arms that box the incoming `Expr` into a
  managed register: `mpz_set`/`mpfr_set` (or, from a machine literal, `set_si`/
  `set_d`) at program precision.
- **Out** — the `compiled_eval` result switch (`compile.c:8220`) gains
  `CT_BIGINT`/`CT_BIGREAL`/`CT_BIGCOMPLEX` arms that unbox the result register's
  arena container into a fresh `Expr`: `expr_bigint_normalize` (so a bignum that
  shrank re-canonicalises to `EXPR_INTEGER`), a fresh `EXPR_MPFR`, or
  `numeric_mpfr_make_complex` (collapsing to real when the imaginary part rounds
  to zero).

`compiled_function_new` (`compiled_function.c:390-408`) and `compile_expr_ex`
(`compile.h:120`) gain the program precision and the `"BigIntegers"` flag
alongside the existing `unsigned flags` (via a small params struct or an
`_ex2` entry — the machine entry points keep their current signatures so no
existing caller changes).

## 10. Optimiser integration

Register-number granularity already absorbs in-place mutation: `op_desc`
(`optimize.c:58`) marks `OP_R_BIN` as reading `a`,`b` and writing `dst`, and that
is honest even when `dst == a` — `mpfr_add(dst,dst,b)` reads `a` before writing,
exactly as a machine `dst = a + b` with `dst == a` does. So "in-place mutation
breaks purity" is a non-issue.

Two real fixes:
1. **Bounded bank predicates.** `is_tile_reg(r)` is the open-ended `r >=
   tile_base` (`optimize.c:172`). A fourth bank above tile makes that
   true for every managed register. Convert it to the bounded
   `r >= tile_base && r < managed_base`, add `is_managed_reg`, and extend
   `compile_optimize`'s signature (`optimize.c:693`) with `managed_base`.
2. **Managed = opaque pointer in v1.** Classify managed registers as
   `is_ptr_reg`, which blocks CSE-to-`MOVE` (the aliasing hazard documented at
   `optimize.c:168-171` — a `MOVE` that copied the pointer would alias two SSA
   values onto one container) and blocks constant-folding (there is no managed
   constant in a `Slot` to fold). They remain **removable when dead** and
   **hoistable when invariant** like tiles — both correct because managed ops are
   side-effect-free and the container is frame-owned. `pass_vn` must not treat a
   managed opcode as foldable.

`OP_MANAGED_MOVE` is therefore always a **deep value copy** (`mpfr_set`/`mpz_set`/
`ncpx_set`). If a later milestone wants managed-op CSE, it flips managed off
`is_ptr_reg` and relies on this value-copy `MOVE` — copy propagation (`alias[]`,
`optimize.c:527-540`) is already correct under value-copy semantics — but that is
out of scope here.

## 11. Zero machine-precision regression (the non-negotiable)

Every mechanism above is either compile-time or gated on the program actually
declaring a managed register. Concretely, for a `WorkingPrecision ->
MachinePrecision` program with no `"BigIntegers"`:

- **Lattice / inference / widening (§6)** run at *emit* time only, and every
  machine gate rejects types ≥ 40 as before — no runtime cost, and the emitted
  bytecode is identical.
- **Opcodes (§8)** are added to `OPLIST`; the computed-goto jump table is indexed
  by opcode number, so existing opcodes dispatch at the same cost regardless of
  how many new ones exist. No machine opcode handler changes.
- **The arena (§7)** has `managed_base == nreg`: never entered, never sized. The
  per-call enter/exit is one `if (prog->nmanaged)` branch at `compiled_eval`
  entry — outside `vm_run`'s inner loop entirely — and it is false.
- **`vm_run`'s inner dispatch loop is not touched at all.** No stop-pc test, no
  new per-instruction branch (the same discipline that kept `OP_APAR` out of the
  hot path, `compile_internal.h:472-475`).
- **Optimiser (§10)** gains a bounded range check that is compile-time; for a
  machine program `managed_base == nreg` makes it identical in effect to today.
- **Boxing (§9)** adds arms reached only for managed argument/result types; the
  machine boxing path is unchanged.

**Proof obligation, not assertion.** This claim is gated by measurement:
- A/B the machine-precision compile benchmarks (the timing half of
  `make check-fastpath-sweep`, plus the compile micro-benchmarks) on `main` vs
  the branch; require parity within noise on every machine workload.
- Byte-compare emitted bytecode for a corpus of machine bodies via
  `CompileDiagnostics` / the disassembler — it must be identical pre/post.

If either fails, the design is wrong, not the measurement.

## 12. Auto-compilation and samplers

`ac_make` (`autocompile.c:61`) currently types all vars `CT_REAL` (or
`CT_COMPLEX`). Add a precision parameter: when a sampler requests precision `n`,
type the vars `CT_BIGREAL`/`CT_BIGCOMPLEX` and instantiate the program at `n`
bits. Add `autocompiled_eval_mpfr` returning an `EXPR_MPFR` /
`Complex[MPFR,MPFR]`. The existing `MATHILDA_NO_AUTOCOMPILE` switch,
`autocompile_enabled()` gate, and the `COMPILE_FOLD_GLOBALS & ~COMPILE_WRAP_INT`
flag discipline are unchanged (the wrap-int guard is moot for MPFR).

Callers already know their `WorkingPrecision` and pass it through:
`src/numerical_calculus/nint.c`, `.../nsum.c`, `src/numerical_roots/findroot.c`,
`src/graphics/plot.c`, and the NDSolve RHS. The uniform pattern
(`autocompile_new` once, `autocompiled_eval_*` per point, `autocompile_free`) is
unchanged; only the precision-aware variant is selected when
`WorkingPrecision > MachinePrecision`. This is where the MPFR win is realised
end-to-end: high-precision sampling that is fully interpreted today.

Bignum is **not** wired into auto-compilation (§2) — samplers are float
workloads, and the compiler cannot infer that an integer body needs GMP.

## 13. Coverage audit

`tools/compile_coverage.py` hardwires "compilable subset = machine number": every
probe shape is `_Real`/`_Integer` (`SHAPES`), and it ratchets a `BASELINE` so a
newly-unlowered head fails the build. Arbitrary precision adds a dimension it
cannot currently express, which creates two silent failure modes: a head that
lowers only at machine precision looks fully covered (false negative), and if an
MPFR kernel ever counts as a "fast path", a machine-only head flips to a spurious
new gap (false positive).

Add a third **precision family** alongside `scalar`/`array` (`FAMILIES`): probe
shapes carrying a `WorkingPrecision` declaration, tracked and ratcheted as a
separate column, with an explicit exemption block (mirroring `EXEMPT`) so partial
precision coverage is a recorded judgement, not a silent signal. **Ship this in
the same change as the first code (Phase 1) or the audit rots.**

## 14. Memory management and leaks

The arena is the primary leak surface. Contract:
- Every arena container is `init`'d once (lazily, first use) and `clear`'d once
  (thread exit). The per-call enter/exit only moves the stack pointer — it never
  `init`s or `clear`s, so there is no per-call `init`/`clear` imbalance to leak.
- The result container is **copied out** into a fresh `Expr` at the boundary (§9)
  and the arena slot is reused, never freed at call scope.
- Worker threads (`par_chunk`) get their own TLS arena; v1 forbids managed
  registers inside `OP_APAR` (§7), so no arena crosses a thread boundary.
- `valgrind --leak-check=full` over the arbitrary-precision tests, on reentrant
  and (once threading touches it) threaded runs, is a release gate. Note the
  macOS valgrind baseline noise (~12.8 KB / 400 blocks) is pre-existing.

## 15. Testing

- `tests/test_compile_arbprec.c`: for a corpus of bodies and argument values,
  assert `Compile[…, WorkingPrecision->n][args]` equals `N[body, n]` to full
  precision, and `Compile[…, "BigIntegers"->True][args]` equals the exact
  interpreter result; assert the **machine-only** compile of each body is
  unchanged (bytecode identity, §11).
- Reentrancy: a compiled managed function that calls a compiled managed callee
  (nested arena frames); recursion depth up to `VM_MAX_CALL_DEPTH`.
- Precision edge cases: `n` at machine boundary; imaginary part rounding to zero
  (collapse to real); bignum that shrinks back to int64 (`expr_bigint_normalize`).
- `USE_MPFR=0` build: MPFR path compiles out, `"BigIntegers"` still works,
  `WorkingPrecision->n` degrades with the existing one-shot warning.

## 16. Milestones

- **P0 — this design doc.** ✅ Reviewed before hot-path code.
- **P1 — substrate + MPFR real.** ✅ Lattice types (`CT_BIGINT/BIGREAL/
  BIGCOMPLEX`); managed register bank + per-call container enter/exit gated on
  `nmgd_slots`; `WorkingPrecision` option + `compile_expr_prec` plumbing;
  MPFR-real opcodes + `mpfr_*` dispatch; boxing via `CompileValue.v.a`; separate
  `emit_mgd` so the optimiser never sees a managed register. **Deviation from
  the plan:** rather than the optimiser `managed_base` fix, a managed program
  forces `COMPILE_NO_OPT | COMPILE_NO_CSE` so `optimize.c` is untouched — simpler
  and equally zero-regression (the optimiser only ever runs on machine
  programs). The coverage precision axis and the warm TLS arena are **not yet
  done** (see below).
- **P2 — bignum + MPFR complex.** ✅ On the same substrate: `"BigIntegers"` +
  `mpz_*` opcodes; `ncpx` complex opcodes (selector-dispatched) + boxing to
  `Complex[MPFR,MPFR]` via `numeric_mpfr_make_complex`.
- **P3 — auto-compilation + samplers.** Substrate ✅ (`autocompile_new_prec` /
  `autocompiled_eval_mpfr`). ✅ **`NIntegrate`, `NSum`, `FindRoot` wired** — each
  compiles the held body in MPFR at high `WorkingPrecision` and evaluates each
  sample point / summand term / Newton step through it, with the same per-point
  faithful-degrade fallback as the machine path (results identical to the
  interpreter; leak-free; existing suites pass). Measured: NIntegrate ~1.5× on
  arithmetic-heavy bodies, FindRoot ~1.2×; NSum benefits less (its
  Euler–Maclaurin contour-derivative path is complex-input). Remaining:
  complex-input autocompile (`autocompile_new_prec_z`) for the NSum contour path
  and complex samplers, `Plot`/`NDSolve`, and end-to-end benchmarks.

**Known-remaining, tracked here so it does not rot:** (a) the coverage-audit
precision family (§13) — not yet added; the existing machine audit is unaffected
because no new machine head was introduced; (b) the warm TLS arena (§7); (c)
managed control flow and managed arrays.

Verification in place: `tests/test_compile_arbprec.c` (14 checks — real, mixed
kernels, divide, constants, precision, bignum exactness, complex, fallback,
machine parity, and the autocompile MPFR substrate); all existing compile suites
(`compile`, `compiledfunction`, `pred_compile`, `autocompile`) pass unchanged;
valgrind shows no leak beyond the macOS baseline; `make check-c99` clean.

Each milestone updates the `Compile` entry in `docs/spec/builtins/` (the two new
options), adds a changelog note under the current week's
`docs/spec/changelog/<Mon>.md`, keeps docstrings current, and passes
`make check-c99`, `make check-compile-coverage`, and the C99 Linux CI.

## 17. Risks and open questions

- **Container lifecycle is the make-or-break (ranked #1).** If the arena is not
  warm across calls, the sampler case degrades from 10–50× to a disappointing
  3–10× (per-call `mpfr_init2`/`clear` dominating a few-op body). The TLS
  stack-nested arena (§7) is the mitigation; P1 must prove it on the sampler
  shape, not just a tight in-body loop.
- **Value proposition honesty (ranked #2).** MPFR 10–50× holds for tens-to-low-
  hundreds of digits with many ops per call; it compresses toward 2–5× at
  thousands of digits (Amdahl: mpfr arithmetic dominates). Bignum is 2–5×,
  collapsing toward ~2× as operands grow. Set expectations accordingly; do not
  promise machine-class speedups.
- **Optimiser (ranked #3).** Mechanical and well-scoped (§10), but the
  `is_ptr_reg`/`is_tile_reg`/`is_managed_reg` classification must be got exactly
  right or a managed register is mis-bound as a tile buffer (`frame_bind_tiles`,
  `compile.c:8137`). Covered by the bounded range check + tests.
- **Coverage audit (ranked #4).** Must ship with P1 or it silently rots.
- Open: whether to key the TLS arena by precision (multiple live precisions per
  thread) or `mpfr_set_prec` on reuse. Default to the latter (cheap when
  precision is stable, the common case); revisit only if a workload mixes
  precisions hot.

## 18. Recommendation

Proceed, MPFR-first within a substrate built to carry all three managed types
(bignum, MPFR real, MPFR complex) from the start, per the approved scope. The
lattice extension is safe-by-default, the boxing needs no struct growth, and the
opcode strategy reuses the existing kernel-pointer machinery — the genuinely new
and load-bearing piece is the per-call TLS precision-keyed arena (§7), which is
also what determines whether the headline speedup is real. Build and measure that
first (P1), hold the machine path to bytecode identity (§11), and only then add
bignum and complex (P2) and the samplers (P3).
