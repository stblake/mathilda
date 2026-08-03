# `Compile[]` — a unified numeric compiler for Mathilda

Status: **design / proposed**. Author: (this doc). Supersedes the bespoke NDSolve
RHS compiler (`src/numerical_calculus/ndsolve_compile.c`), which becomes a client
of the general engine.

---

## 1. Motivation

Numeric-heavy evaluation in Mathilda repeatedly walks and rebuilds `Expr` trees:
every RHS sample in NDSolve, every function sample in `Plot`/`NIntegrate`/
`FindRoot`, every iteration of a numeric `Do`/`For`/`While`/`Nest`. The NDSolve
work proved the payoff of compiling the numeric body once into bytecode and
running it over raw machine numbers: **8×–800×** on real problems, with the
symbolic evaluator kept only as a fallback.

That compiler is a *special case* — straight-line, real-only, fixed variable set.
This document generalizes it into a single `Compile[]` engine that:

1. targets a small **type lattice** (`Boolean`, `MachineInteger`, `Real`,
   `Complex`) rather than one scalar type, keeping a fast real path;
2. supports **control flow** (`If`, `Do`, `For`, `While`, `Nest`, `NestList`,
   `Fold`, `FoldList`, `Which`, `Piecewise`) — unifying the procedural builtins
   with the arithmetic core;
3. is reused **internally** by NDSolve / Plot / NIntegrate / FindRoot / Table and
   also exposed as a user-facing `Compile[]` / `CompiledFunction`.

Decisions already taken (2026-07-26): small lattice; internal engine **and** user
`Compile[]`; control flow designed in from the start.

---

## 2. Goals / non-goals

**Goals**
- One IR + VM shared by all numeric fast paths; the symbolic interpreter is the
  universal fallback (compile-time bail, and optional runtime callback).
- Static type inference → **monomorphic typed opcodes** → no runtime type
  dispatch (the whole point: a `Real` add is one `faddd`, not a tagged branch).
- **First-class arrays.** Lists of machine numbers and **`NDArray` objects** are
  compilable value types (packed `Integer`/`Real`/`Complex` tensors), not a
  deferred add-on. A `List` of machine numbers is coerced to a packed array at
  the boundary; `NDArray` args flow in directly.
- **Complete numeric-function coverage via one machine-precision kernel layer.**
  *Every* numeric function (elementary **and** special — `Gamma`, `Erf`,
  `Bessel*`, `Zeta`, `PolyLog`, …) has a machine `double`/`double complex` kernel
  that avoids building any `Expr`. This is not a new registry: it is the
  **existing `ndkernels` layer** (`symtab_set_ndarray_unary/binary_kernel`,
  real+complex variants, special functions self-registering) — unified so the
  Compile VM (scalar), the NDArray subsystem (elementwise, BLAS), and the
  interpreter's numeric/`Listable` threading over machine numbers all share it.
  "Completeness" = filling the current *degrade sentinels* with real kernels.
- Correctness parity with the interpreter to rounding, verified by tests.
- No `Expr` allocation on the hot path; no memory leaks.

**Non-goals (deferred)**
- LLVM/native codegen (`CompilationTarget->"C"`). The VM is the target; native
  codegen can come later behind the same IR.
- Arbitrary precision (MPFR) inside compiled code — compiled code is machine
  precision by definition; MPFR paths stay interpreted.
- Compiling pattern matching, string ops, or structural list surgery on
  non-numeric elements.
- Ragged / non-rectangular lists (only rectangular → packed tensors compile;
  ragged lists bail).

---

## 3. Prior art: the NDSolve compiler (what we reuse)

`ndsolve_compile.c` already gives us, in miniature:
- a **stack-machine VM** over `double` with an opcode set for arithmetic +
  elementary functions (`run_prog`);
- a **compile-or-bail** front-end (`emit`) that returns `NULL` on any unsupported
  construct;
- **variable resolution** by matching against the problem's state symbols
  (`NameMap`, interned-pointer hash);
- a **dependency/coloring** analysis used for a sparse finite-difference
  Jacobian.

The general engine subsumes the first three. NDSolve keeps the last (colored FD)
as a *client-side* concern layered on top of a compiled function. See §12.

---

## 4. Type system

A minimal **scalar** lattice, ordered by widening, plus **array** types layered
on the same element types:

```
scalars:  Boolean        (1 byte; True/False)
          MachineInteger (int64_t)
          Real           (double)
          Complex        (double complex)

arrays:   Array[elem, rank, dims?]   where elem ∈ {Integer, Real, Complex}
          — backed by an NDArray (packed buffer + shape).
          rank 1 = list/vector, rank 2 = matrix, rank n = tensor.
          A List of machine numbers infers to Array[·, 1]; an NDArray arg
          carries its dtype/rank directly.
```

A register slot holds **either** a scalar (in-line int64/double/complex) **or**
an owned `NDArray*` (for array types). The opcode's static type says which — no
runtime tag.

- **Widening coercions** (implicit): `Boolean ⊄ numeric` (booleans never widen to
  numbers implicitly; guards must be Boolean). `MachineInteger → Real → Complex`.
- **Inference**: bottom-up over the `Expr`. Literals: exact integer → `Integer`,
  real → `Real`, `Complex[...]`/`I` → `Complex`, `True`/`False` → `Boolean`.
  Argument types come from the `Compile` type spec (default `Real`). An op's
  result type is a function of its argument types (a per-op type rule).
- **Overflow**: `MachineInteger` is modular `int64` (WL uses bignums; we do not —
  a compiled `Integer` that overflows is UB in WL Compile too unless
  `"CompileOptions" -> ...`). Documented limitation; `Real` is unaffected.
- **`Power`/`Sqrt`/`Log` domain**: of a `Real` base stays `Real` and yields `NaN`
  out of domain (matches NDSolve; the caller treats non-finite as a failure).
  Declaring an argument `Complex` (or a `Complex` literal) forces the complex
  kernel, which is total. This mirrors WL's real-vs-complex Compile behavior.
- **Inference failure** (can't assign a lattice type, e.g. a subexpression is a
  `List` or an un-inlinable symbolic value) → **bail** (§10).

Runtime storage is **untyped raw slots** (16 bytes: big enough for `double
complex` or an `NDArray*` handle); the *opcode* carries the type, so there is no
runtime tag. `Real` values never pay complex cost. Array-typed slots own their
`NDArray` and are freed when overwritten / at frame teardown (§13).

**Inference over arrays.** Element type widens as for scalars; rank propagates:
scalar⊕array → array (broadcast), array⊕array → array (shape-checked, or bail on
mismatch when shapes aren't known at compile time — a runtime shape guard is
emitted). A numeric function of an array threads elementwise (its result element
type from the scalar rule). Reductions (`Total`, `Dot`, `Max`) lower rank.

---

## 5. IR & bytecode

A **typed, three-address-ish register machine** with structured control flow
lowered to jumps. (A register machine — not the pure stack VM of NDSolve —
because loops mutate named locals and we want basic blocks for jumps.)

### Frame
- `regs[]`: raw 16-byte slots. Registers hold: function arguments, `Compile`
  locals, loop variables (`Do`/`For`/`Nest` accumulators & counters), and
  compiler temporaries. Count fixed at compile time.
- A tiny **expression stack** is still allowed *within a basic block* for
  temporaries, to keep lowering simple; or temporaries can be registers. (Impl
  choice; start with registers-only for uniformity.)

### Instruction
```
{ uint16 op; uint32 dst; uint32 a; uint32 b; union imm; }
```
`op` is **monomorphic and typed** (the type is baked in by inference):

- **Load/const**: `CONST_I/R/C/B dst, imm`, `MOVE dst, a`.
- **Arithmetic** per type: `ADD_I/ADD_R/ADD_C`, `SUB_*`, `MUL_*`, `DIV_R/DIV_C`,
  `NEG_*`, `POW_R/POW_C/POWI` (integer exponent fast path), `INV_R/INV_C`,
  `MOD_I`, `QUOT_I`.
- **Coercions**: `I2R`, `R2C`, `I2C`.
- **Comparisons** (→ Boolean): `LT_I/LT_R`, `LE_*`, `EQ_*`, `GT_*`, `GE_*`,
  `EQ_C` (complex ==), plus `AND/OR/NOT/XOR` on Boolean.
- **Elementary kernels** (unary/binary), per numeric type where defined and
  *hot enough to inline as a dedicated opcode*: `SIN_R/SIN_C`, `EXP_R/EXP_C`,
  `LOG_R/LOG_C`, `SQRT_R/SQRT_C`, trig/hyperbolic/inverse, `ABS_R/ABS_C`
  (→ Real), `SIGN_*`, `FLOOR_R/CEIL_R/ROUND_R`, `RE_C/IM_C/CONJ_C/ARG_C`,
  `MAX_*/MIN_*`, `ATAN2_R`.
- **Generic numeric kernel**: `KERNEL1 head, dst, a` / `KERNEL2 head, dst, a, b`
  — invoke the registered machine kernel for *any* numeric function not worth a
  bespoke opcode (`Gamma`, `Erf`, `Bessel*`, `Zeta`, `PolyLog`, …). The opcode's
  type selects the real vs complex kernel; the pointer is resolved at compile
  time from the shared registry (§8). This is how "all numeric functions" get an
  Expr-free path without an opcode explosion.
- **Array ops**: `EWADD_R/EWADD_C/…` (elementwise arithmetic; SIMD/loop),
  `EWKERNEL head` (apply a numeric kernel over a buffer — the existing
  `ndkernels` vectorized path, with complex-promotion), `BCAST` (scalar⊕array),
  `DOT_R/DOT_C` & `MATMUL_R/MATMUL_C` (BLAS `dgemv`/`dgemm`/`z*`),
  `REDUCE_TOTAL/MAX/MIN`, `PART`/`SLICE`, `MAKEARR` (pack a fixed numeric tuple),
  `LEN`, shape guards. Array temporaries are owned (§13).
- **Control flow**: `JMP target`, `JZ cond, target` (branch if Boolean false),
  `JNZ`. Structured constructs lower to these (§11).
- **Fallback**: `CALLBACK slot, dst` — reconstruct a sub-`Expr` from live
  registers, invoke the interpreter, read the result back (§10). Optional;
  Phase 2+.
- **Return**: `RET a`.

Basic blocks are implicit (targets are instruction indices). No SSA/register
allocation in v1 — the front-end allocates a fresh register per temporary and
reuses named-variable registers; a later pass can coalesce.

---

## 6. Runtime values & VM

- The VM is a `for`-loop over the instruction array with a `switch` on `op`
  (computed-goto dispatch as a later optimization). Each case reads/writes raw
  register slots as the type implied by the opcode — **no tag checks**.
- Real-typed programs touch only the `double` half of slots → identical speed to
  the current NDSolve VM.
- Errors: a domain error yields `NaN`/`Inf`; the runner checks finiteness at
  `RET` (and optionally per-step) and reports failure to the caller, which
  decides (NDSolve: sample failure; user `cf[x]`: return `Indeterminate` +
  message, matching WL's `CompiledFunction::cfn` behavior of falling back).

---

## 7. Compiler front-end

`compile_expr(body, argspec) -> CompiledProgram* | NULL`:

1. **Bind arguments** to registers with their declared types (default `Real`).
2. **Infer + lower** recursively: for each node, infer child types, pick the
   monomorphic opcode, insert coercions, emit into the current basic block.
   Control-flow heads dispatch to §11 lowerers.
3. **Bail** (return `NULL`) the moment a node is not compilable *and* runtime
   callback is disabled/inapplicable — the whole function stays interpreted.
   With callbacks enabled, emit a `CALLBACK` instead of bailing (§10).
4. Record: register count, program, result type, and (for NDSolve reuse) the set
   of argument registers each output depends on.

Variable resolution generalizes NDSolve's `NameMap`: a scope maps symbol →
(register, type). `Compile`'s formal args seed it; `Do`/`For`/`Nest`/`With`/
`Module` loop and local variables extend it in nested scopes.

---

## 8. The unified machine-precision kernel registry

The single most important reuse: Mathilda **already has** a machine-precision,
`Expr`-free numeric kernel layer — `src/ndkernels.c` — that maps each numeric
builtin's libc/`complex.h` math over a buffer, registered on the symbol table via
`symtab_set_ndarray_unary_kernel(name, &kernel)` / `_binary_kernel`, with a
`REXPR` (real) and `CEXPR` (complex-in-`z`) form per function; special functions
(`Gamma`, `Erf`, `Bessel*`, …) self-register their kernels; and *degrade
sentinels* mark functions that currently lack a machine kernel (they fall back to
`Expr`-based List threading).

`Compile[]` does **not** introduce a second registry. It **shares this one**:

- A kernel entry exposes both a **scalar** form (`double→double`,
  `cx→cx`) and a **vectorized** form (buffer→buffer, incl. the BLAS/SIMD fast
  paths NDArray already has). Both are generated from the one `REXPR`/`CEXPR`
  definition, so scalar and elementwise never diverge.
- The **Compile VM** uses the *scalar* form via `KERNEL1/2` (for special
  functions) — elementary hot ops still get bespoke opcodes for speed.
- **Array-typed** Compile ops and the **NDArray** builtins use the *vectorized*
  form via `EWKERNEL`.
- The **interpreter** can route `Listable` numeric functions over lists of
  machine numbers through the same vectorized kernels (an Expr-free fast path in
  `eval.c`), instead of threading and rebuilding `Expr`s per element.

**"All numeric functions" = complete the coverage.** The work is filling the
degrade sentinels with real machine kernels (scalar+vector), drawing on the
existing scalar numeric implementations (many special functions already have
`double complex` kernels — 16+ files under `special_functions/`). Each becomes
callable, Expr-free, from Compile, NDArray, and the interpreter at once.

**Not compilable → bail/callback** (§10): `List` surgery on non-numeric elements,
`Part`/`Map` with symbolic bodies, pattern matching, strings, ragged lists, and
anything with no numeric kernel.

Adding/advancing a function = one kernel definition (scalar+vector from one
macro) + a lowering rule + a parity test. Explicit, local maintenance tax.

Initial compilable heads: `Plus Times Subtract Divide Power Minus Abs Sign Mod
Quotient Sqrt Exp Log Sin Cos Tan Cot Sec Csc Sinh Cosh Tanh ArcSin ArcCos
ArcTan(2) Erf Erfc Gamma LogGamma Beta Zeta PolyGamma Bessel* Max Min Floor
Ceiling Round Re Im Conjugate Arg Chop  ·  Equal Unequal Less LessEqual Greater
GreaterEqual And Or Not Xor  ·  If Which Piecewise Do For While Nest NestList
Fold FoldList Sum Product  ·  Set AddTo SubtractFrom TimesBy Increment Decrement
·  List(→pack) Part Dot Total Length Table(numeric)`.

## 8a. Relationship to the NDArray subsystem

`Compile[]` and `NDArray` become two front-ends over one numeric core:

| Concern | Existing NDArray | Compile reuse |
|---|---|---|
| dtypes | `f64/f32/c64/c32` (`DataType`) | element types of the array lattice |
| elementwise fns | `ndkernels.c` (real+cplx, degrade) | `EWKERNEL` / scalar `KERNEL` share it |
| reductions | `ndreduce.c` (`Total`, …) | `REDUCE_*` opcodes call it |
| structural | `ndstruct.c` (`Part`, `Take`, `Sort`) | `PART/SLICE` opcodes call it |
| linear algebra | BLAS/LAPACK fast paths | `DOT`/`MATMUL` opcodes call it |
| lifetime | NDArray object model | array-typed register slots own handles |

The compiler treats an `NDArray` argument as an `Array[elem, rank]` value and its
buffer as the operand of array opcodes — no copy, no `Expr`. A `List` of machine
numbers is packed to an `NDArray` on entry and unpacked (or returned as an
`NDArray`/`List`) on exit, per the declared result type.

### 8a.1 Delegation, and the two tables that drive it

Most array heads are not re-implemented in the VM. They are **delegated** to the
interpreter's own NDArray entry point, which takes the whole call and promises
the List call's answer — so the compiled subset of these heads *is* the
interpreted one by construction, and the rounding is identical too. Two tables
in `src/compile/compile.c`:

| table | shape | opcode | heads |
|---|---|---|---|
| `ND_FNS` | array → array, plus trailing **integer** arguments | `A_NDFN` | `Reverse` `Sort` `Accumulate` `Flatten` `Transpose` `Take` `Drop` `Differences` `Ratios` `Most` `Rest` `Clip` `RotateLeft` `RotateRight` `MovingAverage` `MovingMedian` `TakeLargest` `TakeSmallest` |
| `ND_REDS` | array → **scalar** | `V_NDRED` | `Mean` `Median` `Variance` `StandardDeviation` `RootMeanSquare` `Max` `Min` |

`Total` keeps its own `V_TOTAL` because an int64 sum must stay exact past
2^53; everything in `ND_REDS` is real-valued, so it needs no such split.
`ND_REDS` entries are **real element type only** — `Mean[{1, 2}]` is `3/2` and
no machine slot holds a Rational — except `Max` and `Min`, which *select* an
element and so carry `int_ok`.

Adding a head is one table row when its entry point already has the shape. The
win is not the operation (it was already a buffer walk) but that a body
*containing* one no longer bails wholesale.

### 8a.2 The declared element type must match the kernel's

An array opcode's declared element type is what every **downstream** opcode
reads the destination slot as. If a kernel writes a buffer of a different
dtype, the consumer does not lose precision — it reinterprets the bits.

The narrowing kernels (`Floor`, `Ceiling`, `Round`, `Sign`, `IntegerPart`,
`UnitStep`: real in, exact `Integer` out) are where this bit. They declared
Real and wrote `NDT_INT64`, so `Total[Floor[v]]` returned `2.96439*10^-323` —
the int64 6 read as a double. `nd_unary_elem` and `nd_binary_elem` now mirror
`ndarray_map_unary` / `ndarray_map_binary` condition for condition; **the pairs
must not drift.**

A missing complex arm on a kernel means three different things and the guards
have to tell them apart:

| `cplx` | `real` | `to_int` | meaning |
|---|---|---|---|
| — | — | — | degrade sentinel: no machine kernel, bail |
| — | — | `to_int_r` | narrowing-only (`UnitStep`, `IntegerPart`): array yes, scalar no |
| — | — | `to_int_i` | integer-only (`MoebiusMu`, `EulerPhi`, `GCD`): integer array yes |
| yes | · | · | ordinary kernel |

Keying a guard on `to_int_r` rather than `to_int` reads the integer-only
kernels as sentinels, which is exactly how `MoebiusMu[v]` stayed uncompilable
after the narrowing case was fixed.

### 8a.3 The coverage audit

`make check-compile-coverage` (`tools/compile_coverage.py`) joins the kernel
registry and `pack.c`'s `AWARE` list against `CompileDiagnostics` over every
typed argument shape. It is a work queue, not a report: the subset is a cliff,
so a head that quietly falls out of it costs every other head in its body. It
ratchets against a checked-in `BASELINE`, so a NEW gap fails the build while
the standing queue is reported.

That queue is written up head by head, grouped by the mechanism that would
close each group, in [`COMPILE_MISSING.md`](../../COMPILE_MISSING.md).

---

## 9. `CompiledFunction` object & calling convention

- New head **`CompiledFunction[argtypes, resultType, <bytecode>]`**, printed
  opaquely (like `InterpolatingFunction`). The bytecode is held as opaque C data
  attached to the `Expr` — reuse whatever mechanism `InterpolatingFunction` /
  `NDArray` use to carry a C payload (verify at impl time; likely a boxed pointer
  node the printer/`expr_free` special-case).
- **Apply**: evaluating `cf[a1,…,an]` dispatches to `builtin_compiledfunction_apply`,
  which coerces the numeric args into the frame, runs the VM, and boxes the
  result back to an `Expr` (`Integer`/`Real`/`Complex`/`True`/`False`). Non-finite
  → `Indeterminate` + `CompiledFunction::cfn`-style message.
- **User syntax**: `Compile[{x, y, …}, body]` (default `Real`), and typed
  `Compile[{{x, _Real}, {n, _Integer}, {z, _Complex}}, body]`. `Compile` is
  `HoldAll` (the body must not pre-evaluate).

---

## 10. Fallback strategy ("MainEvaluate")

Two levels, in order of ambition:

- **v1 — whole-function bail.** If any node is uncompilable, `Compile` returns the
  body uncompiled; internal callers keep interpreting, and user `Compile[]` emits
  `Compile::noinfo`/`::cfse`-style message and returns a `CompiledFunction` that
  simply calls the interpreter (still correct, just not fast). Simple, safe, and
  exactly the NDSolve model.
- **v2 — runtime callback (`CALLBACK`).** For a mostly-compilable body with a few
  uncompilable leaves, emit a callback: at runtime, materialize a sub-`Expr` with
  live register values substituted, call `evaluate()`, read a number back. This
  reintroduces `Expr` building *only for those leaves* (WL's behavior; the source
  of "why is my Compile slow" — surfaced via `RuntimeAttributes`/warnings).
  Deferred until the compilable core is solid.

---

## 11. Control-flow lowering

All structured; lowered to `JZ`/`JMP` over basic blocks. Loop/accumulator
variables are registers, mutated in place (no `Expr`).

- **`If[c, t, f]`** → infer `c:Boolean`; `t`,`f` must unify to a common result
  type (widen). `JZ c, Lf; <t>; JMP Le; Lf: <f>; Le:`. `If[c,t]` → `f = Null`
  only valid where the value is unused (statement position).
- **`Which[c1,e1,…]`**, **`Piecewise`** → chained `If`.
- **`Do[body, {i, i0, i1, di}]`** → integer/real counter register; `body` compiled
  in a scope binding `i`; standard counted loop with `LT`/`ADD`/`JMP`.
- **`For[start, test, incr, body]`** → lower the four parts directly; `start`/
  `incr`/`body` in statement position, `test:Boolean`.
- **`While[test, body]`** → `Lt: <test>; JZ →Le; <body>; JMP Lt; Le:`.
- **`Nest[f, x, n]`** → `x` register, counted loop applying compiled `f` (f is a
  compiled sub-body or a `CompiledFunction`), `n` `MachineInteger`.
  **`NestList`/`FoldList`** produce arrays → require the tensor path (Phase 3) or
  a bail; **`Fold[f, x, list]`** over a *packed numeric list* also Phase 3, but
  `Fold` over a compiled counted range is fine.
- **`Sum`/`Product`** with numeric bounds → counted-loop accumulator (a common,
  high-value case).

Statement vs. value position matters (loops return `Null` unless the body's last
value is used, e.g. `Nest`). The lowerer tracks position.

---

## 12. Integration & reuse

- **NDSolve** (Phase 0/1): replace `ndsolve_compile.c`'s bespoke VM with a
  `CompiledProgram` whose arguments are `(t, Y[])`. The RHS is `d` compiled
  outputs sharing one frame; the **colored-FD Jacobian stays** (it is a client
  algorithm over the compiled RHS + the per-output dependency sets, which the
  general front-end also produces). Must reproduce today's numbers and pass the
  172+27+87 suites unchanged.
- **Plot / Plot3D / ParametricPlot**: compile `f` once, sample the VM in the
  adaptive sampler. (Plotting samples heavily → large win.)
- **NIntegrate / FindRoot / FindMinimum / NDSolve events**: compile the
  integrand/objective.
- **Table / Sum / Product / Nest** with numeric bodies and bounds: compile the
  body + loop when the body is in the subset; else interpret.

Each integration is independently shippable once the engine exists.

---

## 13. Memory management

- A `CompiledProgram` owns: the instruction array, a constant pool, register
  metadata, and (v2) the `Expr` templates for `CALLBACK` nodes. One `*_free`.
- The **VM allocates nothing** for scalar programs (frame is caller-provided or
  arena-reused). **Array programs** do allocate `NDArray` temporaries — these are
  owned by their register slot and freed when the slot is overwritten or at frame
  teardown; a compile-time liveness pass frees eagerly, and an arena/free-list
  keeps allocation churn low across repeated calls (the NDArray object model
  already manages buffer lifetime). Array args are borrowed (not freed) unless
  the program mutates in place.
- `CompiledFunction` `Expr` owns its `CompiledProgram` and frees it in `expr_free`
  (integration point to verify, like `InterpolatingFunction`).
- Follow the builtin ownership contract for `builtin_compiledfunction_apply`.
- Valgrind- and ASan-clean gates (this session already showed ASan catching a
  latent OOB the memcheck missed — ASan is part of the gate).

---

## 14. Errors & numerical semantics

- Compiled code is IEEE-754 machine arithmetic. Parity with the interpreter is
  "to rounding" for values both compute as finite reals/complex; divergence only
  where the interpreter would stay symbolic/exact (documented).
- Domain errors → `NaN`/`Inf` → finiteness check → caller policy (§6/§9).
- Integer overflow: modular `int64` (documented; differs from WL bignums).
- `Real` vs `Complex` of `Sqrt`/`Log`/`Power`: type-driven, as in §4.

---

## 15. Testing strategy

- **Parity battery** (extend `test_ndsolve_compile.c` → `test_compile.c`):
  random inputs, compiled vs interpreter, to machine precision, across every
  opcode and type, including `Complex` and `Integer`, and coercion paths.
- **Control flow**: `If`/`Which`/`Do`/`For`/`While`/`Nest`/`Fold`/`Sum` compiled
  vs interpreted over random inputs; loop-variable mutation; nested loops;
  early-exit conditions.
- **Type inference**: assert inferred result types; assert `Real` stays `Real`
  (no accidental complexification), `Complex` totalizes `Sqrt` of negatives.
- **Bail/fallback**: uncompilable constructs bail (v1) or callback (v2) and still
  produce the correct value.
- **Integration**: NDSolve suites unchanged; Plot/NIntegrate/FindRoot numeric
  results unchanged; performance gates (bench harness, like `bench_assoc`).
- **Memory**: valgrind + ASan clean on the whole battery.

---

## 16. Milestones (control flow designed in from the start)

Even with control flow in the core VM, ship in reviewable increments:

- **M0 — substrate.** Extract a `compile/` module: typed value slots, IR, VM
  skeleton, and a *shared* kernel-registry handle onto `ndkernels`. NDSolve
  migrated onto it, straight-line, real — must reproduce current perf + all
  suites. *No behavior change; pure refactor.*
- **M1 — scalar types + arithmetic core + `Compile[]`.** Full scalar lattice,
  inference, arithmetic/comparison/boolean/elementary opcodes, the generic
  `KERNEL` op over the shared registry (→ **all** scalar numeric functions,
  Expr-free), user `Compile[{args}, body]` for straight-line bodies,
  `CompiledFunction` object + apply. Complex path.
- **M2 — control flow.** `If`/`Which`/`Piecewise`/`Do`/`For`/`While`/`Nest`/
  `Sum`/`Product` (+ loop-var mutation). Delivers the Nest/Do/For/While
  unification.
- **M3 — arrays / NDArray.** `Array[elem, rank]` value type; lists-of-machine-
  numbers coercion; `EWKERNEL`/`EWADD`/`BCAST`/`REDUCE_*`/`DOT`/`MATMUL`/`PART`
  opcodes delegating to `ndkernels`/`ndreduce`/`ndstruct`/BLAS; array temporary
  lifetime. `NDArray` args flow through with no copy. (Core, per the new
  requirement — not deferred.)
- **M4 — complete the kernel coverage + interpreter reuse.** Fill `ndkernels`
  degrade sentinels with real scalar+vector kernels for the remaining special
  functions; route the interpreter's `Listable` numeric threading over machine-
  number lists through the shared vectorized kernels. Wire Plot / NIntegrate /
  FindRoot / Table to auto-compile.
- **M5 — fallback callback (v2)** (partial compilation). Optional, demand-driven.

M0+M1 generalize the NDSolve win to all straight-line scalar numeric callers and
ship `Compile[]`; M2 is the procedural unification; M3+M4 are the array / "all
numeric functions" requirements.

**Status (2026-07-26).** M0 (substrate + NDSolve migration), M1a (generic
`KERNEL` over `ndkernels`), M2 control flow (`If`, `Sum`/`Product`,
`With`/`Module` mutable locals, `Set`/`AddTo`/`SubtractFrom`/`TimesBy`/
`Increment`/`Decrement`, `CompoundExpression`, `Do`/`While`/`For`,
`Nest[Function[u, body], x, n]`), **M1b** (user-facing
`Compile[]`/`CompiledFunction` — `EXPR_COMPILED` opaque atom, numeric bytecode
path + interpreter fallback), and the **M4 auto-compile wiring** of the numeric
builtins are done. The shared adapter `src/compile/autocompile.{c,h}` transparently
compiles a held body once and evaluates it over machine numbers, with interpreter
fallback; it is wired into **Plot/Plot3D** (~215×/~11×), **Table** (real iterator
only, ~128×), **NIntegrate** (1-D machine ~353×, multi-D cubature/Monte-Carlo
~504×), and **FindRoot** (scalar machine real ~19×, multivariate systems ~6.9×).
Remaining: `Which`/`Piecewise`, `Fold`/`NestList`/`FoldList` (need arrays, M3),
M3 arrays/NDArray, M4 full special-function kernel coverage, complex-contour
NIntegrate compilation, and MPFR fast paths.

**Status update (2026-07-27) — M3a done.** Rank-1 machine arrays are now a
first-class value category. `CompileType` packs array types into the same
integer as the scalars (`CT_ARR + 4*(rank-1) + elem`); array registers live in a
dedicated contiguous bank above the scalar registers; operand ownership transfer
is encoded in the consuming instruction's flag bits. Elementwise arithmetic,
broadcast, power, the unary/binary `ndkernel` maps, and the `Total`/`Length`
reductions all delegate to the NDArray subsystem, giving *exact* parity with the
interpreter.

The measurement that should shape M3b: because both paths call the same ND
kernels, delegation alone buys only the per-operation evaluator round-trip —
**2.3× at length 16, 1.0× at length 4096**. Unlike the scalar case (11–353×),
array speed does *not* come from removing interpretation; it comes from removing
**intermediate buffers**. `Sin[v] Exp[-v] + Sqrt[v]` currently makes five
full-length passes and four temporary buffers. The high-value M3b item is
therefore **elementwise fusion** — lowering a chain of elementwise ops to one
buffer pass driven by the existing scalar VM over each element — with `Dot`
(BLAS) and `Part` next. Rank-2, array locals (which need a copy op or handle
refcounting), and the user `Compile[]` array argspec follow.

**Status update (2026-07-27) — coverage gaps closed.** A Newton-fractal
benchmark exposed the *other* way a compiler loses an order of magnitude: not by
being slow, but by silently not compiling. `Do[body, {n}]` — the plainest
counted loop — was outside the subset, so the whole `Compile[]` bailed to the
interpreter, at 37×. Fixed by giving `Do`/`Sum`/`Product` one iterator-spec
parser matching `src/iter.c`, with `Sum`/`Product` still rejecting the bare
count because the interpreter does.

Three follow-ons let a nested `Table` compile end to end: `COMPILE_FOLD_GLOBALS`
(a non-argument symbol holding a machine number folds to that constant — opt-in,
autocompile-only, because those programs die inside one builtin call while a
user `Compile[]` object outlives its scope); `autocompiled_eval_boxed` (an
integer-valued element stays an Integer instead of being flattened to a double);
and inlining of `CompiledFunction` calls inside a compiled body, which removes
an evaluator round-trip that was 72% of per-point cost. Together: **4.45 s →
0.046 s (97×)**, with the un-wrapped form now compiling too.

The design lesson is that the compilable subset is a **cliff, not a slope** — a
construct just outside it costs the entire body, and `compile_expr` returning
NULL is indistinguishable from success at the call site. Any future subset
extension should start by auditing which spellings the interpreter accepts that
the compiler does not.

---

## 17. Risks & open questions

- **Opcode/kernel sprawl** (types × functions). Mitigation: generate the
  per-type opcode tables from a single macro/table; only `Real`+`Complex` need
  full elementary coverage; `Integer`/`Boolean` are small.
- **`CompiledFunction` opaque-payload plumbing**: must confirm how the `Expr`
  layer carries and frees a C payload (`InterpolatingFunction`/`NDArray`
  precedent). *Open — verify before M1.*
- **Type-inference surprises** (real vs complex `Sqrt`, integer overflow). We
  adopt WL's semantics and document them; tests pin them down.
- **Maintenance tax**: the compilable subset must track interpreter semantics.
  Mitigation: parity tests are mandatory per added builtin; bail is always safe.
- **Array shape inference & temporary churn.** Ranks/dtypes are inferred but
  concrete shapes often aren't known until runtime → emit runtime shape guards
  and size checks; array temporaries must be freed eagerly (liveness) and pooled
  to avoid per-call allocation. Mitigation: lean on the NDArray object model +
  an arena; parity/valgrind tests on array programs.
- **Scope creep toward LLVM/native codegen.** Explicitly deferred; the IR is
  designed so a C/LLVM backend can slot in behind it later. Tensors/NDArray are
  now *in* scope (M3), per the requirement.
- **Is M2's control flow worth the VM complexity vs. just compiling the loop
  *body* and iterating in C?** For `Do`/`Sum`/`Nest` over numeric ranges, a
  compiled body + C-level loop (no `JMP`) is simpler and captures most of the
  win; full in-VM control flow matters most for data-dependent `While`/`If`
  inside the loop. *Open: possibly stage "compiled body, C loop" before full
  in-VM control flow.*

---

## 18. Recommendation

Proceed M0 → M1 → M2 → M3 → M4. M0 is a pure, low-risk refactor that de-risks
everything else by putting NDSolve on the shared substrate first. M1 generalizes
the proven win to all straight-line scalar numeric callers, ships user-facing
`Compile[]`, and — via the generic `KERNEL` op over the existing `ndkernels`
registry — already gives an Expr-free path for *every* scalar numeric function.
M2 delivers the `Nest`/`Do`/`For`/`While` unification. **M3 makes lists of
machine numbers and `NDArray` objects first-class**, delegating to the NDArray
subsystem's kernels/BLAS rather than duplicating them. M4 completes numeric-
function coverage (fills the degrade sentinels) and extends the Expr-free path to
the interpreter's list threading.

Two things to nail before building: (1) the `CompiledFunction` opaque-payload
mechanism (§17); (2) exposing a **scalar** entry point from each `ndkernels`
kernel (they are currently vectorized-only) so the Compile VM and NDArray share
one definition — a small refactor of `ndkernels.c`'s `REXPR`/`CEXPR` macros.
