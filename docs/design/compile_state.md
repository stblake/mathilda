# Compiler build — state & handoff

## 0e. M8 — machine integers as a peer of machine reals (2026-07-29)

READ THIS FIRST if you are touching anything integer-shaped.

**The engine had a real type and an integer type in name only.** 18 integer
opcodes against 69 real; no integer buffers; and integer arithmetic WRAPPED,
which is the one thing the engine forbids — a compiled body must answer
identically to the interpreter or not answer at all.

**Overflow is now detected, and the fallback IS the answer.** `ci_add` / `ci_sub`
/ `ci_mul` / `ci_neg` / `ci_abs` / `ci_powi` in `compile_internal.h` return TRUE
ON OVERFLOW (the `__builtin_*_overflow` convention) under a `__GNUC__` guard with
a strict-C99 fallback; `-DCOMPILE_NO_OVERFLOW_BUILTIN` selects the fallback and
the two are asserted equal at the boundaries. A VM body is one line:
`OP(ADD_I): IOP(ci_add(RA.i, RB.i, &RD.i));` where `IOP` branches to `vm_fail`.

- **The optimiser folds too.** `fold_op` uses the same helpers and REFUSES to
  fold an overflowing operation, so the optimiser cannot bake in a wrapped
  constant that the runtime check would have deferred.
- **`OP_LOOP` is deliberately unchecked** while `OP_INC_I` is checked. The loop
  index only overflows if the limit sits within one step of INT64_MAX, which
  takes ~10^18 iterations to reach; `INC_I` is also `Increment`'s opcode, and
  `x = 2^63-1; x++` is ordinary.

**GOTCHA THAT COST A DEBUG CYCLE — do not repeat it.** The build switch was first
called `COMPILE_WRAP_INT`, the same name as the *public flag bit* in `compile.h`.
`#ifdef COMPILE_WRAP_INT` in the VM was therefore true in EVERY build, so the
checks were compiled out and every test silently exercised the wrap path. The
build switch is now `VM_NO_INT_CHECK` (named for the VM, like `VM_NO_THREADED`).
An `#ifdef` on a name that `compile.h` `#define`s as a value is always true.

Also re-learned: **`make | head` kills make with SIGPIPE before the LINK**, so
the binary is older than the objects (see the memory note). And **same-second
mtimes defeat make** — after restoring a file with `cp`, `touch` it.

**The A/B flag is per-INSTRUCTION, and that is why it is free.** `COMPILE_WRAP_INT`
(a `compile_expr_ex` flag, surfaced as `RuntimeOptions ->
{"CatchMachineIntegerOverflow" -> False}`) sets `IF_NOCHK` in the instruction's
`flags`, stamped in `ins_f` — the one place every instruction passes through.
`IOP` reads that bit only AFTER `&&` has confirmed an overflow, so the
no-overflow path is byte-identical either way. **Measured: the option buys 0%;
only `-DVM_NO_INT_CHECK` recovers the 0–4% the detection costs.** Report it as a
semantics switch, not a speed one.

`autocompile.c` masks the flag off explicitly (`COMPILE_FOLD_GLOBALS &
~COMPILE_WRAP_INT`). Auto-compiled Plot/Table/NIntegrate must never wrap: the
user never asked for compilation there.

**Result HEADS are a separate axis from values, and only a sweep finds them.**
A value compares equal whether it is `35` or `35.`, so the parity tests are
structurally blind here. Sweeping all 103 `NumericFunction` heads with `_Integer`
arguments (`Head[h[3]]` vs `Head[Compile[…][3]]`) found 68 divergences; 10 were
integer-closed and are fixed via `INT_CLOSED` in compile.c — ONE table consulted
by both `infer_type` and `emit`, because those two dispatchers drifting apart is
the failure mode. The other 58 are the accepted divergence (symbolic `Sin[3]`,
Rational `3/2`). **`LegendreP` looks integer-closed and is not**: `LegendreP[2, 2]`
is `11/2`.

Every integer-closed head ALSO has a registered real kernel behind it — reaching
the kernel first is exactly how they came back as Reals — so these branches must
sit BEFORE the generic kernel dispatch and must DECLINE (not bail) for
non-integer arguments, or the real fast path disappears silently.

**Integer arrays: `NDT_INT64`, compiler-internal.** Reusing `EXPR_NDARRAY` rather
than a private buffer keeps the whole M3a ownership discipline intact — `Slot.arr`
is still an `Expr*`, `OP_ARR_FREE` / the frame sweep / `AF_FREE_A|B` are
untouched. Containment is what makes it safe: `ndt_from_string` will not produce
it, `NDArray[…]` never infers it, and the `Compile[]` boundary always unpacks to
a List of Integers.

**The lossy-`double` trap, hit twice.** `ndt_get` / `ndt_set` route through a
`double` and are exact only to 2^53. Both the boundary PACK (`flatten_into`) and
`Total` (`ndred_total_all`) silently went through them, so
`Total[{9007199254740993, 1}]` came back one short — with no error, because every
individual step "worked". Exact `ndt_get_i` / `ndt_set_i` now serve every real
read/write of an int64 buffer, and `OP_V_TOTAL` sums in `int64` with the same
overflow rule. **Before delegating an integer array to any ND-layer function,
check what it accumulates in.**

**Still real/complex only:** elementwise FUSION (the tile opcodes `VADD_R`,
`VMUL_C`, … have no integer forms), so an integer chain takes the ordinary
per-element loop. `ndstruct_sort` compares through a double as well and would
need the same treatment. Both bail cleanly, which is the correct default.

`Bit*` is NOT a coverage gap to close: `BitAnd`/`BitOr`/`BitXor`/`BitNot`/
`BitShiftLeft`/`BitShiftRight` are unimplemented in the INTERPRETER (they return
unevaluated), so compiling them would answer where the interpreter declines.
Implement them there first.

**Two tests were passing for the wrong reason** and are worth remembering as a
class: `must_bail_raw("NestList integer body", …, in2, RI2, 1)` passed with
nargs=1, so `n` was a free symbol and the body bailed on THAT, not on the integer
history its comment blamed. A bail test proves nothing unless you know WHICH bail
you got — `CompileDiagnostics` reports it.


Snapshot for resuming the `Compile[]` numeric-compiler work with fresh context.
Companion to [`compile.md`](compile.md) (the full design) and the memory files
`project_compile_engine`, `project_autocompile_numeric_builtins` (read those
too).

_Last updated: 2026-07-29 (M7: compile-time function values, and the functional
heads — Nest/Fold/FixedPoint/NestWhile/Map/Scan/Table.  Earlier: M5 optimiser,
coverage audit, any-rank arrays, strip-mined fusion, Expr-level CSE, per-call
frames, OP_CALL; M6 bail diagnostics + nine more auto-compiled builtins;
M3c indexed arrays)._

---

## 0d. M7 — the functional heads, and the abstraction they needed

`Nest`, `Fold`, `FixedPoint`, `NestWhile`, `Map`, `Scan` and `Table` are in the
subset. What unblocked all seven at once was **one** thing: a compile-time
function value.

**`fn_resolve` / `emit_apply` / `infer_apply` replace `extract_function`.** The
VM still has no runtime function value and should not gain one — a lambda is
inlined. What was missing was the compile-time question "which function is this,
and how do I paste it in?", answered in one place instead of once per head.
Retargeting `Nest` onto it made `Nest[Cos, x, n]`, `Nest[#^2 &, x, n]`,
`Nest[Composition[Sin, Cos], x, n]` and `Nest[compiledFn, x, n]` compile with no
lowering of their own. **That is the test of whether the abstraction is right:
adding a head is now a lowering, not a lowering plus a private parser.**

`FN_HEAD` (a bare `Sin`, `Plus`, …) is the awkward case, because `emit_node`
dispatches on a head name plus an `Expr**`, not on `Val`s. Rather than refactor
~45 lowerings to take `Val`s — a large change with no behaviour delta — it
synthesizes `h[$1, …, $n]` over reserved placeholder symbols bound to the
argument registers, the same scaffolding trick the multi-iterator `Do` uses. A
few nodes per *compile*, never per call, and every head the compiler already
knows becomes usable as a function value for free.

**`Slot[k]` needs a binding mechanism of its own** — it is an `EXPR_FUNCTION`,
not a symbol, and `scope[]` matches interned symbol POINTERS. One flat frame is
exact, not an approximation: `substitute_slots` (`src/purefunc.c:76`)
deliberately does not recurse into a nested `Function`, so only the innermost
frame is ever visible. And slots are hidden under a NAMED lambda, because the
interpreter's named path substitutes names only — `Function[u, # + u]` leaves a
live `Slot[1]` in its answer, so that is not a machine number.

**`Map` has two lowerings and it is a cost split, not a subset split.** When the
body threads elementwise, `Map[Function[u, body], v]` IS `body` with `u` bound to
the whole array, which the existing fusion strip-mines and threads — ten lines,
6.6x over the general per-element loop. The legality gate is `fuse_listable`
(would the interpreter thread this?) plus "exactly one array leaf and it is the
parameter" plus a result-type check. Everything else takes the general loop,
which is correct for any body. Same shape as `Part`'s inline/delegated split.

**Three real bugs fell out, all of the same family — a value whose KIND or
FAILURE was decided in the wrong place:**

1. **A built array took its result kind from an unrelated argument.**
   `Compile[{{v,_Real,1}}, ConstantArray[1., 3]][NDArray[{1., 2.}]]` returned an
   `NDArray` where the interpreter's `ConstantArray` returns a `List`. The
   boundary chose from the arguments alone; a body that CONSTRUCTS its result has
   no kind to inherit. `Val.built` / `CompiledProgram.result_built` now carry it.
   Every array-constructing head would have inherited the bug, so this had to be
   fixed before `Table` landed, not after.
2. **`compiled_eval_real` computed the abort flag and dropped it.** Invisible for
   as long as no opcode an all-Real program could contain was able to fail (array
   ops can; an all-Real program has no array registers). `OP_FAIL` made it
   reachable in one step.
3. **`Do`/`While`/`For` in result position returned `0`, not `Null`.**
   Pre-existing and reachable as `Compile[{n}, Do[…, {n}]][3]`. `Null` is not
   worth a fifth lattice type — it would ripple through `num_common`, `coerce`,
   `finite_result` and `cf_unbox` — so a statement-shaped head is declined in the
   one position where the difference is observable.

**Where these heads decline, they decline because compiling would DIVERGE.**
That distinction is the whole design and each case is worth remembering:
`Table` needs integer iterators (the interpreter walks a real one by repeated
addition against a `1e-14` slack, so a closed form differs in the last bits and
at the endpoint in the element COUNT) and a non-integer body (no integer dtype);
`Map` needs rank 1 (rank ≥ 2 maps over ROWS) and a result element type equal to
the source's (the repack uses the SOURCE dtype, so `Map[Abs, complexvec]` comes
back complex-typed); `Fold` over `{}` and `Nest`/`FixedPoint` with a negative
count fail the call because all three are UNEVALUATED in the interpreter.

**`SameQ` is not `Equal`, and that is what makes an unbounded iteration
terminate.** `expr_eq` calls two NaNs the same (`src/expr.c:622`), so a
`FixedPoint` whose orbit reaches NaN stops there. `OP_SAMEQ_R`/`OP_SAMEQ_C`
implement exactly that (componentwise for complex, mirroring `expr_eq`'s
recursion into `Complex[re, im]` — NOT `creal==creal && cimag==cimag`, which
differs on a NaN component). Both also carry the interpreter's `ITER_SAFETY_CAP`
of 10⁶ and `OP_FAIL` on reaching it, so the compiled path gives up exactly where
the interpreter does.

**A prerequisite nobody would guess: `Fold` left an `NDArray` unevaluated.** An
`NDArray` is atomic, so `Fold`'s element walk looked straight past it while the
identical `List` call folded fine. There was nothing to be parity WITH, so
compiling `Fold` first required giving the *interpreter* a packed path. Expect
the same for `Select`, `Join`, `First`, `Differences`, `RotateLeft/Right`,
`Riffle`, `Partition`, `TakeWhile`, `AllTrue`/`AnyTrue`/`NoneTrue` — all of them
return unevaluated on a packed argument today.

**Measure with a PACKED argument.** Over a plain `List` a compiled `Map` reads
1.0x, for two reasons that have nothing to do with the loop: both sides are then
dominated by packing 200 000 `Expr` nodes at the boundary and unpacking them, and
the interpreter's `Map` over a `List` already has the legacy `numloop` fast path,
so the "interpreted" side is not interpreted. Packed in and packed out at 200k
elements: `Map[u^2 + 1. &]` 277x, `Map[Sin[u] Exp[-u] + Sqrt[u] &]` 109x,
`Map[If[…] &]` (general loop) 47x, `Fold[Plus, 0.]` 19x.

**Update (2026-07-30): the boundary is gone, and it was never the marshalling.**
Automatic packed arrays (`docs/design/packed_arrays.md`) mean an ordinary
`Range[1., 200000.]` *is* a buffer, so the packed-in half now happens without the
caller arranging it, and the boundary returns the result buffer with its
presentation set instead of materialising it. Measured against the same value with
`MATHILDA_NO_PACK=1`: `u^2 + 1.` **50x**, a compiled `Map` body **4.8x**,
`Fold[Plus, 0.]` **4.2x** — and each now matches the *visible*-`NDArray` timing to
within a few percent, which is the real test that no marshalling is left.

What actually cost the 50x was the transparency gate, not the copying. A
`CompiledFunction`'s head is an `EXPR_COMPILED` with no `SymbolDef`, so the gate's
allowlist read it as unaware and materialised the very buffer the boundary exists
to borrow. An allowlist keyed on the symbol table cannot see a head that is not a
symbol.

The list of heads above that "return unevaluated on a packed argument" is also
stale: `Fold`, `Select`, `Join`, `First`, `Differences`, `RotateLeft/Right`,
`Riffle`, `Partition` and `TakeWhile` all have interpreter paths now, and the
remaining gap is recorded in `packed_arrays.md` §6.

---

## 0c. M3c — indexed arrays, and what a stencil costs

`Part` is in the subset now, both directions, which is what makes a hand-written
numeric kernel compilable at all. Before this, one `u[[i, j]]` anywhere put the
*whole* body on the interpreter.

**Two lowerings, chosen by the shape of the subscript list — and together they
cover every spec `Part` accepts on a dense array.** One scalar subscript per
axis lowers inline (`A_AXIS` per axis, then the `A_LOAD` fusion already had);
everything else — Span, All, position lists, partial indexing, and any mixture —
delegates to the interpreter's own `ndarray_part` via `A_PART`. The split is not
a subset restriction, it is a cost distinction: the inline path allocates
nothing and is where a stencil lives, the delegated path allocates its result.

**`A_AXIS` folds four things into one instruction** — multiply by the axis
length, 1-based resolution, negative-from-the-end resolution, and the range
check. Doing them as separate opcodes would have tripled the index cost. **The
check has to be per axis:** `m[[1, ncols + 5]]` is inside the buffer as a linear
offset and reads the next row, which a flat bounds check cannot catch.

**`A_LOAD` had to stop being pure.** It was pure while the only writer was
fusion, which stores solely into a fresh result buffer. Once user code can write
a buffer it also reads, a pure load is CSE'd across the store (`u[[1]] = 1.;
a = u[[1]]; u[[1]] = 2.; b = u[[1]]` gives `b - a == 0`) and LICM hoists a
loop-invariant load out of the loop that mutates it. Both are regression tests
now. Generalises: *the moment a read and a write of the same memory can appear
in one body, the read is not pure, however local the analysis looks.*

**Ownership is decided by the register bank.** Array arguments live below
`nlocals` in the scalar range; owned arrays live in the `ARR_VREG` bank. So
"may I write through this?" is `reg >= ARR_VREG` — one comparison, no extra
field for a construction site to forget. Argument arrays are read-only on
purpose: they are borrowed, and for a `List` argument packed at the boundary a
write would vanish without a trace.

**What it costs.** The 2-D wave-equation stencil worked example is 72 instructions per interior
point and runs a 641x641 grid for 639 steps (2.6e8 updates) in 27 s — **569x**
the same march interpreted, and **1.9x** Wolfram Language 14.0's own `Compile`.
Notably WL's `CompilationTarget -> "C"` (verified to be genuine native code, not
a silent fallback) is *slower* than its bytecode VM on this body past n = 101.
That is a caution for this project's own "native backend next" plan: in a
tensor-heavy kernel the cost is array element access, not dispatch.

**Two interpreter bugs fell out of the parity tests**, both of the same kind —
the compiled path implementing something the interpreter only pretended to do.
`TimesBy` was registered with no implementation (and `*=` / `/=` did not parse),
and `Part` assignment into an `NDArray` silently ignored every non-integer spec
while working fine on the equivalent `List`. Both fixed; `ndarray_part_set` now
shares the per-axis selector with `ndarray_part`, so a write names exactly the
elements a read of the same spec would.

---

## 0b. M6 — bails are no longer silent

`CompileDiagnostics[argspec, expr]` (and `MATHILDA_COMPILE_DIAG=1` for the
internal wirings) reports whether a body compiles and, if not, the **innermost**
subexpression that stopped it. Built as ONE wrapper around `emit` — the lowering
proper is now `emit_node`, and `emit` records the node on the way back up — so
no bail site knows diagnostics exist and a bail added tomorrow is diagnosed the
day it is written. `EmitMark` carries the record so speculative lowering (fusion
probing) rolls it back; a speculative failure is not a bail.

**Codegen: constant operands now live in the instruction** (`K_BINK`, 18
opcodes). `pass_vn` already tracked which registers hold constants, so a binary
op with exactly one constant operand becomes its immediate form and DCE deletes
the `CONST`. The surviving operand always becomes `a` (one register read, no
branch); a comparison with the constant on the left swaps the PREDICATE rather
than adding eight opcodes, which is NaN-safe where negating would not be. Each
form is one arithmetic operation, so no FMA contraction is possible and the
opt/no-opt `memcmp` gate covers it.

**And the lesson from measuring it: a third fewer instructions bought ~9%, not
a third.** Horner deg-40 went 121 → 81 instructions for 181 → 173 ns/call. The
removed `CONST`s were the cheapest instruction in the set and the multiply-add
chain around them is serially dependent — the VM was never instruction-count
bound. Quote the wall-clock, not the instruction count.

Three findings from having the diagnostics:

1. **Two tests had silently rotted.** Both built their "interpreter reference"
   out of `Zeta` because it had no machine kernel. It got one during M5, and both
   became compiled-vs-compiled. **A reference built out of a coverage gap expires
   when the gap closes.** The replacement is a user DownValue (`uncid[t_] := t`)
   — uncompilable by construction, and exactly value-preserving, so parity checks
   became exact instead of approximate.
2. **The complex gap is measurable, not theoretical.** `ComplexPlot[Zeta[z], …]`
   on a 40×40 grid takes 0.7 s because `Zeta` has a real kernel and no complex
   one. That is `NUMERIC_FUNCTION_MISSING.md` Class B showing up as wall-clock.
3. **Wiring a builtin can give zero speedup and still be right.** `NSum` gained
   nothing until its *second* sampler (the Euler–Maclaurin continuous-x path) was
   covered too. Check for a second path to the same body before blaming the
   engine.

---

## 0. M5 — read this first

Four things changed since M4. Full write-up in
[`../spec/changelog/2026-07-27.md`](../spec/changelog/2026-07-27.md).

1. **There is a middle end now** — `src/compile/optimize.c`. CFG + backward
   liveness, per-block value numbering (folding/CSE/copy-prop), liveness-driven
   DCE, and LICM. `compile_internal.h` holds `Slot`/`Instr`/the opcode enum;
   `OPLIST` carries a KIND per opcode and drives the enum, the VM jump table and
   the optimiser's property table from one list. **Add an opcode there or the VM
   crashes.** `COMPILE_NO_OPT` is the A/B switch and the gate is *bitwise*
   identity, not accuracy.
2. **Scalar dispatch is ~1.45× faster** — 295 → 200 ns/call on the Horner
   micro-bench at `-O3` (3.7 → 2.5 ns per arithmetic op), from lazy operand
   addressing in `NEXT()` plus the optimiser and `OP_LOOP`.  (An earlier `-O0`
   measurement made this look like 2.0×; see the measurement trap below.)
3. **Any rank works** and needed no new machinery: the delegated ND path was
   already rank-general, and the compiler's own rank-1 front gate was the only
   blocker. `Total` stays rank-1 on purpose (it reduces the leading axis only).
4. **Coverage is audited by the build, and is now 93 of 103 (90%), from 55** —
   `tests/test_compile_coverage.c` probes every `NumericFunction` head by
   actually compiling `Head[x]` … `Head[x,y,z,w]`. It fails if a listed gap
   starts compiling too, so the exception list cannot go stale.
   **All 10 remaining entries are deliberate exclusions, not pending work**: five
   heads the *interpreter* leaves unevaluated on machine reals (`BesselJZero`,
   `BarnesG`, `Hyperfactorial`, `Factorial2`, `FactorialPower` — a kernel there
   would answer where the interpreter declines), and five that are not one
   machine number (`GCD`, `LCM`, `DigitSum`, `ReIm`, `QuotientRemainder`).
   Heads whose real signature is not a flat scalar list need a `PROBES` entry
   (`Clip`, `Rescale`, `HypergeometricPFQ`) — pFq was reported as a coverage gap
   for a while purely because the probe used the wrong shape.

5. **Elementwise fusion is ON by default and strip-mined** (`COMPILE_NO_FUSE`
   disables it). Each opcode processes a tile of `VBLOCK` = 64 elements in a
   vectorisable C loop, so dispatch amortises 64-fold and temporaries stay in L1.
   1.9-3.4x over the delegated NDArray path at `-O3`.

5b. **Fused MAPs are threaded** above `NDARRAY_THREAD_THRESHOLD` (`OP_APAR`,
   `COMPILE_NO_PAR` for A/B). 3.2x / 5.5x / 6.6x at 1M elements on 16 cores,
   rising with per-element cost as memory bandwidth stops being the limit.
   **Maps only** — a map is bit-identical however it is split, a reduction is
   not (FP addition is not associative), so reductions stay serial.
   The loop is LIFTED into a standalone sub-program at finalize (jump targets
   rebased, own `OP_RET`) so a worker calls the ordinary `vm_run`; teaching the
   VM a stop-pc would have cost a compare on every dispatch forever. Extraction
   MUST run after the optimiser — LICM and compaction move instruction indices.
   Each worker gets its own frame (registers copied, tiles re-pointed), so there
   is no lock anywhere; TSan clean.

**MEASUREMENT TRAPS — read before trusting any number here.** Three found so far,
each of which silently invalidated a whole round of figures:

1. The `tests/` CMake build had **no `CMAKE_BUILD_TYPE`**, i.e. no optimisation
   flags at all. Now `Release`, but anything measured in that tree before
   2026-07-27 was `-O0`. Check `grep CMAKE_BUILD_TYPE tests/build/CMakeCache.txt`
   before quoting a figure.
2. `tests/CMakeLists.txt` **never defined `MATHILDA_THREADS`** (the makefile
   always did), so every `nd_parallel_for` in the test build ran serially and no
   test had ever reached the threaded ND path. Now set via `find_package(Threads)`.
3. `bench_compile.c` timed with **`clock()`, which is CPU time summed over all
   threads** — a perfectly scaling parallel region reports as N times SLOWER.
   The threaded map first measured 0.56-0.83x while actually running 6.4x
   faster. Any benchmark a threaded path can reach must use
   `clock_gettime(CLOCK_MONOTONIC)`, as the rest of `tests/bench_*.c` does.

6. **CSE now fires, at the Expr level.** The bytecode value-numbering CSE was
   defeated structurally (`binop`/`unop` write into an operand's register, so the
   entry is invalidated by the instruction that created it). A repeated subtree
   is instead hoisted to a register reserved BELOW the temp stack by raising
   `nlocals`. 1.48x on a body with repeats. `COMPILE_NO_CSE` is the A/B switch.
7. **Frames come from the C stack, and `OP_CALL` exists.** A program is now
   reentrant and thread-safe (it was neither: one shared register file, and after
   M5b one shared set of tile buffers). A non-inlined compiled callee is CALLed
   rather than bailing the whole body. Inlining stays the default and the choice
   is a size-based cost model.

**Findings that should shape the next round:**

- **Self-recursive `Compile[]` still does not compile**, and `OP_CALL` does not
  change that: `Compile[]` deliberately does not fold globals, so a body cannot
  resolve the symbol it is about to be assigned to. Would need a self-reference
  patch at object construction.
- **`infer_type` is load-bearing now.** It used to be a routing hint — a branch
  that reported a scalar where the value was really an array cost nothing,
  because the ND layer picked its own result dtype. Fusion sizes the output
  buffer from it, so the same sloppiness is now a wrong answer. Five branches had
  to be fixed (unary math, 2-arg `Log`, `ArcTan`, `Power` with a literal
  exponent, the binary-kernel path). **Any new head must propagate array-ness.**
- **The tile bank has its own aliasing rule.** A tile op whose result element
  WIDTH differs from its operand's (complex-to-real: `Abs`, `Re`, `Im`, `Arg`)
  must not reuse the operand's register — it reads `double _Complex*` and writes
  `double*`, and the compiler may assume those cannot overlap. This is invisible
  until the loop vectorises.
- **Fusion may only thread where the interpreter threads.** Listable is necessary
  but not sufficient: with BOTH operands arrays, a head with a registered binary
  kernel has no interpreter path at all (`ArcTan[nd, nd]` comes back
  unevaluated), so fusion must decline it too.
- **The parity tests keep finding INTERPRETER bugs, not kernel bugs.** Three so
  far, each surfaced only because a fast double path was written next to the
  existing one and the two were then compared over a few hundred points:
  `ProductLog` wrong on `[0.35, 1/e]` (seed gap), `Zeta` at 0, and
  `HypergeometricPFQ` at negative real `z` — off by 5e-2 at `z = -40` because
  `machine_sum` summed a series whose largest term is `e^|z|` times its own sum
  in plain doubles. The fix pattern that generalises: **measure the cancellation
  (`max|term| / |sum|`) rather than estimating it from the argument**, and
  re-sum through the MPFR path at `53 + lost + 16` bits, rounding back so the
  argument precision still decides the output type. Costs nothing when nothing
  cancels. When a parity test fails, check which side is wrong before assuming
  it is the new one.

---

## 1. What exists and works today

The engine lives in **`src/compile/`** — a typed register-machine bytecode VM that
evaluates numeric expressions over machine numbers with no `Expr` allocation and
no runtime type dispatch (the opcode carries the type).

- **`compile.c`** — the VM + the emitter (`compile_expr`). Type lattice
  `CT_BOOL/CT_INT/CT_REAL/CT_COMPLEX` plus packed array types (§3); bottom-up
  inference; widening coercions; monomorphic typed opcodes; stack-discipline
  register allocation; reusable per-program `frame` (no per-call malloc);
  computed-goto dispatch on GCC/Clang. Coverage: full arithmetic, comparisons,
  boolean, elementary + special-function kernels (via the shared `ndkernels`
  registry), `If`, `Sum`/`Product`, `With`/`Module` locals, `Set`/`AddTo`/…/
  `Increment`, `CompoundExpression`, `Do`/`While`/`For`, `Nest`, and **rank-1
  machine arrays** (§3).
- **`compiled_function.{c,h}`** — user-facing `Compile[argspec, body]` →
  `CompiledFunction` object (new `EXPR_COMPILED` atom, refcounted immutable
  payload). Numeric args run the bytecode; symbolic args / uncompilable bodies
  fall back to the interpreter. `HoldAll | Protected`. **No array argspec yet.**
- **`autocompile.{c,h}`** — the adapter that lets the numeric builtins compile a
  held body once and evaluate it over machine numbers, with per-point interpreter
  fallback. Header is deliberately self-contained (NO `<complex.h>`/`compile.h`
  include — `double _Complex` is a builtin type — so it never leaks the `I`
  macro into callers).

### Auto-compile wiring

| Builtin | Chokepoint | Scope | Speedup |
|---|---|---|---|
| Plot / Plot3D | `plot_eval_fn`, `plot3d_eval_z/_fn` | real result; non-real → exclude point | ~215× / ~11× |
| Table | `table.c` numeric-range loop | **inexact iterator only**; exact untouched; nested via FOLD_GLOBALS | ~128× (97× on the Newton fractal) |
| NIntegrate 1-D | `ni_eval_at` | finite/half/whole-line machine; complex fallback | ~353× |
| NIntegrate multi-D | `ni_mc_sample` | cubature + Monte-Carlo | ~504× |
| FindRoot scalar | `fr_eval_with_bindings` (pointer-identity `main_f`) | machine real Newton/Secant/Brent + FD | ~19× |
| FindRoot systems | `fr_run_newton_system_real` | per-component + Jacobian programs | ~6.9× |
| ContourPlot | `eval_at` (`GridCtx.ac`) | both grid paths (function + equation form) | 5.2× |
| DensityPlot | `dp_eval` | grid | 2.0× |
| ComplexPlot | `cp_eval` | **complex argument** (`autocompile_new_z`) | 2.5× |
| ParametricPlot / PolarPlot | `param_eval_at` | each coordinate its own program; 1- and 2-iterator | 8.4× / 16.8× |
| ParametricPlot3D | `param3d_eval_at` | three programs | 1.1× |
| VectorPlot | `vp_eval` | both field components | 1.5× |
| StreamPlot | `eval_field` | both components; hottest (several samples per RK step) | 7.1× |
| NSum / NProduct | `ns_term_machine` **and** `ns_eval_complex_machine` | machine precision only; MPFR untouched | 6.4× / 8.0× |

The sub-2× rows are not underperforming fast paths — the diagnostic confirms all
three bodies compile. Sampling is no longer their bottleneck: they build one
`Rectangle`/`Arrow`/`Polygon` `Expr` per cell and that construction now dominates.

**`NSum` had two samplers, not one.** Wiring only the integer-indexed term
chokepoint gave *no speedup*, because the Euler–Maclaurin correction samples the
same summand at CONTINUOUS real x through a different function. Covering both
took `Sin[n]/n^3` from 122 ms to 18 ms. When a wiring shows no gain, look for a
second path to the same body before concluding the fast path is not the
bottleneck.

**Shared correctness rule:** a compiled REAL program returns non-finite exactly
where the interpreter would produce a COMPLEX value (`Sqrt` of a negative).
Callers that exclude non-real points (Plot) drop it; callers needing the
contribution (NIntegrate/FindRoot/real-Table) fall back to the interpreter at
that one point. MPFR paths are UNTOUCHED everywhere. Uncompilable bodies (e.g.
`Zeta`) → `NULL` program → interpreter.

### Tests
`tests/test_compile.c` (engine, scalar + array), `tests/test_compiledfunction.c`
(user `Compile[]`), `tests/test_autocompile.c` (all 6 builtin wirings, parity +
fallback + oscillatory-regression + systems). All pass; `leaks`- and ASan-clean.

### Milestones: M0, M1a, M1b, M2 (a/b/c), M3a (rank-1 arrays), M4 (auto-compile) DONE.

### Coverage-gap fixes (2026-07-27) — read this before optimising anything

A Newton-fractal benchmark ran ~30× slower than expected. **Nothing was slow.**
The body simply did not compile, and a bail is silent: it costs an order of
magnitude and looks exactly like a working fast path. Four fixes:

1. **`Do` accepted only `{i, lo, hi}`.** `Do[body, {n}]` — the plainest possible
   loop — made the *whole* `Compile[]` bail. 37× on its own. `Do`/`Sum`/
   `Product` now share `loop_spec_parse`, matching `src/iter.c`: `n`, `{n}`,
   `{i, hi}`, `{i, lo, hi}`, `{i, lo, hi, di}` (literal nonzero `di`; the loop
   test flips to `OP_GE_I` when `di < 0`). `Sum`/`Product` reject the bare count
   because the interpreter does — **the compiled path must never answer where
   the interpreter declines**, even when the compiler could.
2. **`COMPILE_FOLD_GLOBALS`** (`compile_expr_ex`): a non-argument symbol holding
   a machine-number OwnValue folds to that constant. This is what lets a nested
   `Table` compile — it desugars to nested `Table`s and the inner body sees the
   outer variable only through the symbol table. Opt-in, and used **only** by
   `autocompile.c`, whose programs live and die inside one builtin call. User
   `Compile[]` must NOT fold: its object outlives the defining scope. Staleness
   is impossible because any body that could reassign a global (a `Set` to a
   non-local) already bails.
3. **`autocompiled_eval_boxed`**: `Table` was flattening every element through a
   `double`, so `If[…, 1, 2]` produced `1.` instead of `1`. Element type is
   user-visible in a returned list; `Plot`/`NIntegrate`/`FindRoot` keep the
   unboxed real path.
4. **CompiledFunction inlining**: `newt[x + I y, 25]` inside a `Table` body used
   to bail, costing an evaluator round-trip per point (~2.3 µs of a 3.1 µs call,
   72% of runtime). `emit`/`infer_type` now inline the callee. The subtle part
   is scope: arguments lower in the CALLER's environment, then the body lowers
   with **only** the parameters visible (`c->inlining` makes `arg_find` return
   -1), or a caller argument sharing a name with a callee global would capture
   it. Depth-capped at 8 for self-reference.

Net on `Table[newt[x + I y, 25], {y,-1,1,2./199}, {x,-1,1,2./199}]`:
**4.45 s → 0.046 s (97×)**. Written inline without `Compile[]` it now also fully
compiles (0.042 s) and agrees element-for-element.

**Generalisable lesson:** the compilable subset is a cliff, not a slope. Any
construct just outside it costs the whole body. When a compiled path
underperforms, first check that it *compiled at all* — a `compile_expr`
returning NULL is indistinguishable from success at the call site. The counted
iterator forms are now covered; the same audit is worth doing for the other
constructs the interpreter accepts in more spellings than the compiler does.

---

## 2. Build & disk discipline (IMPORTANT — read before rebuilding)

- Main binary: `make -j4` (produces `./Mathilda`). ~53 MB objects + 7 MB binary.
- Tests: `tests/CMakeLists.txt` compiles all `COMMON_SRC` **once** into an
  `OBJECT` library `mathilda_common`, spliced into every test target via
  `$<TARGET_OBJECTS:mathilda_common>` (commit `9ec6ad9`). Build a target:
  `cmake --build tests/build --target mathilda_common <target> -j4`.
  Use **absolute paths** — the Bash tool's working directory persists between
  calls, so a bare `cd tests/build` can silently run a stale binary.
- ASan gate: configure a throwaway tree in the scratchpad, build only the target
  you need, run, then delete it (~174 MB):
  `cmake -S tests -B <scratch>/asan -DCMAKE_BUILD_TYPE=Debug
   -DCMAKE_C_FLAGS="-fsanitize=address -fno-omit-frame-pointer -g -O1"
   -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address"`, then
  `ASAN_OPTIONS=detect_leaks=0 <scratch>/asan/compile_tests`.
- **Disk gotcha:** editing a `COMMON_SRC` file forces a `mathilda_common` rebuild
  + relink of every built test exe (~6 MB each). Build only the target you need.
- **The `.claude/worktrees/` risk** — stale git worktrees from prior
  `isolation: worktree` subagents each held a full pre-OBJECT-lib `tests/build`
  (~1.3 GB); 7.5 GB was reclaimed once. Six worktrees still exist (~456 MB, with
  uncommitted risch/poly/core edits) pending the user's decision to
  `git worktree remove --force`. If disk fills, check `du -sh .claude/worktrees`
  and `git worktree list` first.

---

## 3. M3a as built — rank-1 machine arrays

### 3.1 Type representation

`CompileType` packs array types into the same integer as the scalars:

```c
typedef enum { CT_ERR = -1, CT_BOOL = 0, CT_INT, CT_REAL, CT_COMPLEX, CT_ARR = 4 } CompileType;
#define CT_ARRAY(elem, rank)  /* CT_ARR + 4*(rank-1) + elem */
#define CT_IS_ARRAY(t) / CT_ELEM(t) / CT_RANK(t)
```

Every field that already carried a `CompileType` — `infer_type`'s result,
`Val.type`, `Ctx.scope[].type`, `Ctx.arg_types[]` — carries array types unchanged.
`CT_ERR` exists so the enum's underlying type is signed and the pre-existing
`(int)t < 0` checks stay well defined. `num_common` was extended: an array
absorbs a scalar (broadcast), two arrays must agree on rank, element types widen
as scalars do.

`CompileValue` gained `Expr* a` — an EXPR_NDARRAY. **Argument arrays are
borrowed; a result array is owned by the caller.**

### 3.2 The three decisions that make it safe

1. **Array registers live in their own bank.** Array temps are allocated into a
   virtual range tagged `ARR_VREG (0x40000000)`; at finalize, `patch_reg` rewrites
   the tags to `maxreg + k`, so every array register sits in one contiguous bank
   *above* the scalar registers. A slot is therefore either always-array or
   never-array, and teardown can never mistake a `double` for a pointer. (Jump
   targets also live in the `b` field but are small integers, so the rewrite skips
   `OP_JMP`/`OP_JZ`.)
2. **Ownership transfer is encoded in the instruction, not in separate frees.**
   `Instr` gained a `uint16_t flags` field in what was padding after `op` (verified:
   `sizeof(Instr)` is 32 either way). `AF_FREE_A`/`AF_FREE_B` tell the op to free
   the operand *after* reading it, so the result may still reuse an operand's
   register exactly as the scalar `binop` does. `free_if_tmp` (value is dead)
   emits an explicit `OP_ARR_FREE`; `pop_tmp` (consumer handles the free) does not
   — that split is what keeps a temp inside a `Do` loop from accumulating one
   buffer per iteration.
3. **Teardown is a range sweep.** Every array register is NULLed before a call
   (the frame is reused across calls, so a handle from a previous call must never
   be touched) and swept afterwards — including on the abort path, which is why
   `vm_run` gained a `bool* failed` out-param.

### 3.3 Opcodes and delegation

`OP_ARR_FREE`, `OP_V_EW` (→ `ndarray_elementwise`, which already folds broadcast
scalars), `OP_V_POW` (→ `ndarray_elementwise_power` / `ndarray_scalar_power` /
`ndarray_base_scalar_power`), `OP_V_KERN` (→ `ndarray_map_unary`), `OP_V_KERN2`
(→ `ndarray_map_binary`), `OP_V_TOTAL` (→ the new `ndred_total_all`), `OP_V_LEN`.

Lowerings with no direct ND counterpart: `a - b` → `a + (-1)b`, `-v` → `(-1)v`,
`a / b` → `a * b^-1` (scalar divisor reciprocated in a register), `Sqrt[v]` →
`v^0.5` (Sqrt has no registered ND unary kernel).

`ndreduce.h` gained **`ndred_total_all(const Expr* a)`** — full reduction taking
the *array* rather than the enclosing call, so the VM needs no wrapper call node.
Shares `ndred_total`'s summation, so rounding is identical.

### 3.4 Guard rails

- **`scalar_only()` in `binop`/`unop`/`kern_unop`/`kern_binop`** is the single
  choke point: every scalar opcode is emitted through one of those four, so any
  head with no array lowering (comparisons, `Max`/`Min`, `Mod`, …) bails the
  moment an array reaches it rather than reinterpreting a handle as a double.
- **Explicit array bails at every `OP_MOVE` site** (`If` branches, `Sum`/`Product`
  accumulator, `With`/`Module` locals, `Set` family, `Nest` state): a MOVE
  duplicates a handle without duplicating ownership, so array values are not in
  the M3a subset there. Emit-time types are ground truth for these checks —
  `infer_type` can under-report arrays inside nested expressions, which is fine
  (it only routes; a wrong "false" still ends in a bail).
- **A borrowed argument array cannot be the whole result** (`compile_expr`
  rejects it), otherwise the caller would free a value it does not own.
- **Element-type promise check at runtime:** if a program promised a real element
  type and the ND layer returned a complex dtype (`Sqrt[v]` with a negative
  entry, `Log[v]` of a negative), the op fails → whole call fails → interpreter.
  This is the scalar "non-finite where the interpreter would go complex" contract
  lifted to buffers.

### 3.5 THE FINDING THAT SHOULD SHAPE M3b

Parity with the interpreter is **exact** (max relative error 0.0 across 34 array
bodies) — because both paths call the same ND kernels. Which is also why the
speedup is small: delegation only removes the per-operation evaluator round-trip.

| body | length | speedup |
|---|---|---|
| `Total[Sin[v] Exp[-v] + Sqrt[v]]` | 16 | **2.3×** |
| same | 4096 | **1.0×** |

Contrast the scalar engine's 11–353×. **Array speed does not come from removing
interpretation; it comes from removing intermediate buffers.** That body
currently makes five full-length passes and allocates four temporary buffers.

So the highest-value M3b item is **elementwise fusion**: lower a chain of
elementwise array ops into ONE buffer pass that runs the existing scalar VM over
each element (the scalar engine is already the right per-element machine — this
is mostly a matter of compiling the chain's scalar body once and looping it over
the buffers, with the reduction folded into the same pass for `Total[...]`).
Expect that to be where the real multiple appears, especially for long vectors.

---

## 4. NEXT: M3b

In rough value order:

1. **Elementwise fusion** (§3.5) — the actual performance win.
2. **`Dot`/`MatMul`** → `ndarray_dot2` (`ndarray.h:118`) / BLAS `dgemm` wrapper
   `dot2` (`src/linalg/dot.c:60`). Needs rank-2 first.
3. **Rank-2 arrays.** The type encoding already supports rank ≤ 8
   (`CT_MAX_RANK`); `compile_expr` currently rejects rank ≠ 1 and `Total` requires
   rank 1 (`ndred_total_all` collapses every axis, so a rank-2 `Total` — which
   reduces only the LEADING axis — needs the `ndred_total` call path or a new
   partial-reduction entry point).
4. **Array locals / `If` branches / `Nest` state** — needs either an
   `OP_ARR_COPY` or handle refcounting (`expr_ref` on EXPR_NDARRAY; verify
   whether `expr_copy` shares or deep-copies the buffer before relying on it).
5. ~~**User `Compile[]` array argspec**~~ — DONE. A `List` argument is packed at
   the boundary and freed after the call; the result kind follows the argument
   kind (Lists in → List out, NDArray in → NDArray out), and a body that BUILDS
   its array with no array argument returns a List, because that is what the
   interpreter running the same body returns.
6. ~~**`Part`/`Slice`**~~ — DONE in M3c (§0c), full spec vocabulary, both
   directions. Still open from this item: `MAKEARR` (pack a fixed tuple), Int
   arrays, and `Table` as an array constructor inside a body (today that needs
   `ConstantArray` plus a `Do` loop).

### ND delegation API reference (still accurate)

Kernel structs & helpers live in **`src/ndarray.h`** (there is NO `ndkernels.h`).
Buffer = row-major flat, complex interleaved (re,im), NOT C99 `_Complex`.

- **Elementwise arith / broadcast:** `ndarray_elementwise(Expr** args, size_t n,
  bool is_plus)` (`ndarray.h:127`) — folds scalar operands into one constant, so
  broadcast is free. Power: `ndarray_elementwise_power`, `ndarray_scalar_power`,
  `ndarray_base_scalar_power` (`ndarray.h:134-145`).
- **Unary/binary kernel map:** `ndarray_map_unary(a, k)` /
  `ndarray_map_binary(a0, a1, k)` (`ndarray.h:185-186`); look the kernel up via
  `symtab_lookup(head)->ndarray_unary_kernel` / `->ndarray_binary_kernel`
  (`symtab.h:88-89`). Kernel structs `NDUnaryKernel` (`ndarray.h:166`),
  `NDBinaryKernel` (`ndarray.h:179`). A NULL `cplx` is the *degrade sentinel*
  (no machine kernel — `Sqrt`, `Zeta`, `Erfi`, … ) → bail to the interpreter.
- **Reductions:** `ndred_total_all(a)` (new, operand-taking) or
  `ndred_total/mean/max/min/…(Expr* res)` (`ndreduce.h`), which take the WHOLE
  call Expr (borrowed) and reduce the LEADING axis.
- **Build/read buffers:** `expr_new_ndarray(rank, dims, buf, dtype)` (adopts
  `buf`, copies `dims`; `expr.h:164`); `ndt_get/ndt_set(buf, k, dtype, &re, &im)`.
  `NDType`: `NDT_FLOAT64=0` (default), `FLOAT32`, `COMPLEX64`, `COMPLEX32`.
  **All ND fast paths allocate a new buffer and never mutate inputs.**

### Watch-outs carried forward

- Frame reuse ⇒ never free-on-overwrite; only the explicit `ARR_FREE`/flag/result
  model.
- Array bodies in tests must be parsed **without** evaluating: with the array
  parameters still free symbols the evaluator rewrites the constructs under test
  (`Total[v]` collapses to `v`, `With[{u=v},…]` inlines, `Length[v]` → `0`).
  For the same reason the interpreter reference substitutes via
  `env_new`/`env_set`/`replace_bindings`, not `ReplaceAll` — `ReplaceAll`
  evaluates its first argument before binding.
- Verify the `CompiledFunction` `expr_free` path frees an array-returning
  program's result correctly once the user argspec lands.

---

## Two-argument `If` (2026-07-30)

`If[test, var = val]` — with no else branch, in statement position — had no
lowering. Because the compilable subset is a **cliff, not a slope**, that one
form cost the whole body: a Sieve of Eratosthenes and a Collatz longest-chain
search, both otherwise entirely within the subset, silently ran interpreted.

It lowers exactly like `While`: emit the guard, `JZ` over the body, answer the
integer 0, and join `stmt_valued_head` so that 0 can never be a program's result
(the interpreter answers `Null` when the test is False). `compile.md` §11 had
specified this — *"`If[c,t]` → `f = Null` only valid where the value is
unused (statement position)"* — it was simply never emitted.

- Collatz longest chain below 10⁶: **240 s → 4.65 s** (52×)
- Sieve of Eratosthenes to 10⁷: did not compile at all → **666 ms**

Both are now *faster than Wolfram's* `Compile` on the same source (1.50× and
1.18×). The general lesson stands and is worth repeating: **audit which
spellings the interpreter accepts that the compiler does not.** `If[c, t]` is
not an exotic construct — it is how anyone writes a conditional store — and it
was missing for as long as control flow has existed.

`CompileDiagnostics` named the cause exactly (`"Subexpression" -> "If[len > bl,
bl = len]"`), which is the whole reason it exists; the failure was in not asking
it about these bodies sooner.
