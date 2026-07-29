# Compile[] M5 — optimising code generation, coverage, any-rank NDArray

Plan: `~/.claude/plans/i-would-like-to-fluffy-whisper.md`
Handoff: `docs/design/compile_state.md` §0
Changelog: `docs/spec/changelog/2026-07-27.md`

Baseline (2026-07-27, before any change): Horner dispatch 649 ns/call for 80
arith ops = 8.0 ns/op; mixed-libm 150 ns/call; array len-4096 1.0x, 61 ns/element.

## Done

- [x] **Benchmark gate** — `tests/bench_compile.c`. Primary assertion is that
      every benchmarked body *compiled*; plus a machine-independent ratio check
      that the optimiser never makes a body slower.
- [x] **Lazy operand addressing** — `NEXT()` no longer computes all three operand
      pointers per instruction. 649 → 427 ns/call (1.53x).
- [x] **Bytecode optimiser** — `src/compile/optimize.c`: CFG + backward liveness,
      per-block value numbering (folding/CSE/copy-prop), DCE, LICM. Plus
      `compile_internal.h` (KIND-carrying `OPLIST` drives enum + jump table +
      property table) and `OP_LOOP`. 427 → 335 ns/call. Nest 1.29x, Newton 1.25x.
- [x] **Optimiser correctness gate** — 18 bodies, opt vs `COMPILE_NO_OPT`, must
      agree *bitwise* (`memcmp`: NaN == NaN, -0.0 != +0.0).
- [x] **Any rank** — lifted the rank-1 front gate; rank 2/3/4 exact parity.
- [x] **Elementwise fusion** — `COMPILE_FUSE`, correct at every rank, real and
      complex, gated on `Listable`. **Off by default: not yet faster.**
- [x] **Coverage audit** — `tests/test_compile_coverage.c`. 103 NumericFunction
      heads, 55 compile, 48 listed gaps; fails in both directions.
- [x] **Closed 3 silent gaps the audit found** — real `Mod`, real `Quotient`,
      real `Arg`. Each verified against the interpreter first; parity 0.0.
- [x] **Call boundary** — no malloc/free per `cf[x]`; all-Real signatures take
      the unboxed `compiled_eval_real` path.
- [x] **Block strip-mining** — each opcode processes a tile of VBLOCK=64 elements
      in a vectorisable C loop. Fusion is now ON by default and 1.9-3.4x over the
      delegated path. Turned up five array-blind `infer_type` branches, a tile
      aliasing rule, and one place fusion had to be made less capable than it
      could be (ArcTan[nd,nd], which the interpreter declines).
- [x] **Fixed the test build's missing CMAKE_BUILD_TYPE** — it was compiling at
      -O0, so every absolute figure measured there had been inflated.
- [x] **Stress tests** — 340 body x length combinations across the tile boundary,
      21 fused-vs-delegated, 120 randomised trees, 200 repeated calls.
      leaks-clean and ASan+UBSan clean.

- [x] **Expr-level CSE** — the bytecode CSE was defeated structurally, so
      repeated subtrees are hoisted to registers reserved below the temp stack.
      1.48x on a body with repeats. Found a real `pass_vn` bug on the way (a
      write cleared aliases pointing AT the register but not its own).
- [x] **C-stack frames** — a program is now reentrant and thread-safe; tile
      storage moved into the frame.
- [x] **`OP_CALL`** — a non-inlined compiled callee is called rather than bailing
      the whole body. Inlining stays the default via a size cost model.
- [x] **User `Compile[]` array argspec** — `Compile[{{v, _Real, r}}, body]`. A
      List argument is packed at the boundary and the result KIND follows the
      argument kind, so the compiled path and the interpreter fallback agree.
- [x] **Exponential-integral kernels** — 9 heads, coverage 55 -> 64 of 103.
- [x] **Plot3D leak** — `build_surface_primitives` freed the primitive array but
      not the primitives, on the no-surface path.

## Next, in value order
- [ ] **Self-recursive `Compile[]`** — still cannot compile. `Compile[]`
      deliberately does not fold globals (the object outlives its scope), so a
      body cannot resolve the symbol it is about to be assigned to. Needs a
      self-reference patch at object construction.
- [x] **Thread the strip-mined loop — DONE. 3.2x / 5.5x / 6.6x at 1M elements.**
      `OP_APAR` fans a fused MAP out over `nd_parallel_for` and jumps past the
      loop, or falls through to run it serially — the fallback IS the serial
      loop, so declining is free. Maps only: a map is bit-identical however it
      is split (asserted by `memcmp` under `COMPILE_NO_PAR`), a reduction is not.
      The loop is lifted into a standalone sub-program at finalize so workers
      call the ordinary `vm_run`. TSan clean, zero leaks.
      Also fixed two things that had been hiding results: `tests/CMakeLists.txt`
      never defined `MATHILDA_THREADS` (so NO test had ever run the threaded ND
      path), and `bench_compile.c` timed with `clock()` — CPU time summed over
      threads, which reports perfect scaling as an N-fold SLOWDOWN.
- [x] **Fill the remaining kernels — DONE. Coverage 93/103 (90%), from 55, and
      no pending numerics remain.** `src/special_functions/sf_machine.c` covers
      the exponential-integral family, Erfi, ProductLog, Fresnel, PolyGamma,
      HurwitzZeta, HarmonicNumber, Zeta, Fibonacci/LucasL, Pochhammer, Binomial,
      LegendreP, Airy x4, and the final twelve — PolyLog, LerchPhi, QPochhammer,
      BesselI/K and the four hypergeometrics (one `sf_machine_pfq`; 0F1/1F1/2F1
      are wrappers, since 1F1 canonicalises to pFq before the evaluator sees it).
      UnitStep/Clip/Rescale stay bespoke lowerings — their result TYPE is the
      difficulty, not the numerics.
      **The 10 remaining listed gaps are all deliberate**, and the audit's header
      now says so rather than calling them pending: `BesselJZero`, `BarnesG`,
      `Hyperfactorial`, `Factorial2`, `FactorialPower` (the *interpreter* leaves
      these unevaluated on machine reals, so a kernel would answer where it
      declines) and `GCD`/`LCM`/`DigitSum`/`ReIm`/`QuotientRemainder`
      (exact-integer, or not a single machine number).
- [x] **Airy's uncovered band — DONE.** `2.5 < |x| < 8` is covered by Taylor
      marching of `y'' = x y`, seeded from whichever expansion is exact at the
      nearer end. Errors 1e-16 to 2.6e-15, at least as good as the two regions
      already accepted. The direction is the whole problem: each solution is
      marched where it DOMINATES (Bi forward from 2.5, Ai backward from 8), since
      marching toward a recessive solution amplifies seed error by the dominant
      one's growth — 2.4e5 across this band.
- [ ] **Complex arguments — 41 heads compile for `_Real` and bail for
      `_Complex`.** See `NUMERIC_FUNCTION_MISSING.md`: this is now the single
      largest source of silent interpreter fallback in the engine.
      **Done: `Sign`, `FractionalPart`, `Rescale`** (52 -> 55 of 103).
      `Floor`/`Ceiling`/`Round`/`IntegerPart` turned out TYPE-BLOCKED, not
      mechanical — they return `Complex[Integer, Integer]` and the lattice has no
      complex-integer type. Remaining: 28 special functions needing genuine
      complex numerics with matching branch cuts.
      **Done: `Gamma`, `LogGamma`** (55 -> 57). No new numerics — both already
      had a `double complex` Lanczos INSIDE the interpreter; exposing and
      sharing it makes compiled and interpreted agree bit for bit. Check for an
      existing `static double complex` implementation before writing a kernel.
      The real work was a branch-cut bug: the `Re < 1/2` reflection used a
      PRINCIPAL log in both the machine and MPFR paths, correct only in the
      strip `-1 < Re < 0`.
      **Done: the exponential-integral family** — ExpIntegralEi, LogIntegral,
      SinIntegral, CosIntegral, SinhIntegral, CoshIntegral (57 -> 63). Again no
      new numerics: each had a `double complex` series dead behind
      `#ifndef USE_MPFR`. The work was a CANCELLATION GATE with the budget set
      from measurement (1e9 -> 1.3e-8 error; 1e3 -> 4e-13 at the SAME coverage).
      The gate lives in the ABI wrapper, not the series, so the interpreter's
      no-MPFR last resort keeps answering where it has no fallback.
      **Next:** extend their coverage with the complex continued fraction for
      E1 (the real path already has one past |x|=40), which would lift the
      whole family's declines; then the complex hypergeometrics.
- [x] **`CompileDiagnostics` — DONE.** `CompileDiagnostics[argspec, expr]`
      reports whether a body compiles and, on failure, the INNERMOST
      subexpression that stopped it; on success the result type, instruction
      count, CSE count, and the instruction count with the optimiser off.
      `MATHILDA_COMPILE_DIAG=1` prints the same whenever an auto-compiled
      builtin falls back. One wrapper around `emit` (the lowering proper is now
      `emit_node`) does it, so no bail site knows diagnostics exist and a bail
      added tomorrow is diagnosed the day it is written.
      It immediately found two tests that had rotted into vacuity when `Zeta`
      gained a machine kernel — **never build an interpreter reference out of a
      coverage gap; build it out of a user DownValue**, which cannot expire and
      is exactly value-preserving.
- [x] **Auto-compile nine more builtins — DONE**, measured against the previous
      build rather than a proxy: PolarPlot 16.8x, ParametricPlot 8.4x,
      NProduct 8.0x, StreamPlot 7.1x, NSum 6.4x, ContourPlot 5.2x,
      ComplexPlot 2.5x, DensityPlot 2.0x, VectorPlot 1.5x,
      ParametricPlot3D 1.1x.
      New `autocompile_new_z` — ComplexPlot needs a complex ARGUMENT, not just a
      complex result, and its subset is genuinely smaller.
      `NSum` gave ZERO speedup until its second sampler (the Euler–Maclaurin
      continuous-x path) was covered too.
      Also fixed a real leak on the way: ParametricPlot/ParametricPlot3D/
      PolarPlot passed `expr_new_real(t)` inline to `symtab_add_own_value`,
      which COPIES — one leaked node per sample point (529 blocks / 34 KB for
      one default ParametricPlot).
- [x] **Constant operands folded into the instruction — DONE.** New `K_BINK`
      kind, 18 opcodes (real/int add/sub/mul/div + the four order comparisons).
      Rewritten in the optimiser using the constant tracking `pass_vn` already
      had; DCE removes the dead `CONST`. Horner deg-40 121 -> 81 instructions,
      `1.5 + 2.5 x` 5 -> 3.
      **But measure the time, not the count: -33% instructions bought ~9%**
      (`Table` of a degree-5 polynomial over 10^6 points, 140 -> 129 ms), and
      1-3% on bodies that are not constant-heavy. The removed `CONST`s were the
      cheapest instruction in the set and the surrounding chain is serially
      dependent — the VM was never instruction-count bound.
      Anti-vacuity guard added: the A/B gate passes whether or not the rewrite
      fires, so a separate check asserts the instruction count really drops.
- [ ] **Next codegen step is a native backend, not more peepholes.** At ~2 ns
      per instruction the VM is memory-port bound (Instr load, operand load,
      result store, jump-table load), which is why removing a third of the
      instructions moved 9%. A multiply-add superinstruction is the one
      remaining cheap idea and it needs floating-point contraction turned off
      to stay bit-identical — `-ffp-contract=off` on one file, or `#pragma STDC
      FP_CONTRACT OFF`, whose GCC support is unreliable.
- [ ] **Plot primitive construction** — now the bottleneck for DensityPlot /
      VectorPlot / ParametricPlot3D: one `Rectangle`/`Arrow`/`Polygon` `Expr`
      per cell. All three bodies compile; sampling is simply no longer where the
      time goes. Not a compiler problem.
- [ ] **Native backend** (`CompilationTarget -> "C"`), behind a build flag.
- [ ] **Pre-existing diffuse leaks in `builtin_parametricplot`'s option
      handling** (~1–9 blocks per option path, unrelated to sampling). Surfaced
      while leak-checking the wiring above; not chased.

## Review — the last 12 numerics

Coverage 84 -> 93 of 103. The work split unevenly: nine of the twelve were
ordinary double implementations, and the twelfth (`HypergeometricPFQ`) cost more
than the other eleven combined, for two reasons worth recording.

**pFq was never actually an uncovered head.** The audit probes with
`Head[x]` … `Head[x,y,z,w]`, and pFq takes *two lists and a scalar*. It had a
lowering the whole time; the probe was measuring the wrong signature and
reporting a gap that was really a defect in the audit. `PROBES` now carries the
shape override alongside `Clip` and `Rescale`. Worth remembering that a coverage
number is only as honest as its probe.

**The parity failure was the interpreter's, not the kernel's.** This is now the
third time (`ProductLog`, `Zeta`, pFq) — writing a second implementation and
diffing it over a few hundred points is turning out to be a better bug-finder for
the existing numerics than it is a risk to them. Do not assume the new side is
wrong.

`machine_sum` summed pFq in plain doubles. For negative real `z` the series
alternates and `max|term|` exceeds the sum by ~`e^|z|`, so `1F1(1;2;-40)` came
back with 5.3e-2 relative error against its closed form `(E^z-1)/z` — the bits
were gone before the loop ended, and no amount of extra terms recovers them. The
fix measures the loss instead of guessing it from `|z|`: track `max|term|/|sum|`,
report `log2` of it, and re-sum through the MPFR path at `53 + lost + 16` bits
when it exceeds 4, rounding back. Exact for every `p`, `q` and parameter set, and
free where nothing cancels (all-positive terms give `lost = 0`). Machine
arguments still give a machine answer; it is now the correctly rounded one.

Also moved the NDSolve bail example onto the structural half of `KNOWN_GAPS`.
`Zeta`, `AiryAi` and `PolyLog` were each used as "a head with no machine kernel"
and each in turn started compiling; `BarnesG` and `QuotientRemainder` cannot move
without breaking interpreter parity, so they will not.

Verified: `compile_tests`, `compile_coverage_tests`, `compiledfunction_tests`,
`autocompile_tests`, `ndsolve_compile_tests`, `hypergeopfq_tests`,
`numeric_tests`, `beta_tests`, `integrate_ramanujan_tests`, `legendre_tests`,
`sum_tests` all pass; `bench_compile` within gate; `leaks` reports 0 bytes on
all four compile suites.

## Review — threading the fused map

3.2x on `Sqrt[v] + v^2`, 5.5x on `Sin[v] Exp[-v] + Sqrt[v]`, 6.6x on
`Gamma[v] + Erf[v]`, at 1M elements on 16 cores. The gain rises with per-element
cost, which is the right shape: the cheap body moves 16 MB for three flops and is
bandwidth-bound long before it is core-bound.

Two design choices worth keeping:

**Fall through, don't branch to a fallback.** `OP_APAR` either fans out and jumps
past the loop, or falls through into the serial loop that was already emitted
immediately after it. There is no second implementation to keep in step, so
declining — too small, no threads, a worker failed — is always safe by
construction rather than by care.

**Lift the loop instead of teaching the VM to stop.** A worker needs to run one
instruction range, and the obvious way is a stop-pc argument to `vm_run` — which
costs a comparison on every dispatch, forever, for one loop. Copying the range
out at finalize with rebased targets and its own `OP_RET` costs nothing at run
time. It has to happen after the optimiser: LICM and compaction both move
instruction indices, so a range recorded at emit time names the wrong
instructions by the time it is used.

**Two measurement bugs, and the more embarrassing one is mine.** The first round
of numbers said threading was a 0.56-0.83x LOSS. It was not: `bench_compile.c`
timed with `clock()`, which sums CPU time over threads, so perfect scaling reads
as an N-fold slowdown. I had already recorded this exact trap in memory for the
NDArray parallel map and walked into it again. Before that I also chased a wrong
explanation (concurrent first-touch page faults) far enough to write a pre-fault
pass, which measurement then showed cost about 5%; it was removed. The lesson
that generalises: instrument the region directly before theorising about why a
number is bad — `clock_gettime` around `nd_parallel_for` gave the answer in one
run and said 6.4x while the benchmark was still calling it a loss.

Separately, `tests/CMakeLists.txt` never defined `MATHILDA_THREADS`, so the whole
threaded ND layer had been dead in the test build — every `nd_parallel_for` there
compiled to its serial fallback and no test had ever exercised it. Now on; the ND
suites pass with it.

Verified: compile, compile_coverage, compiledfunction, autocompile,
ndsolve_compile, ndarray, ndarray_functions, ndarray_reduce, ndarray_linalg,
mapthread, linalg, ndsolve all pass; `bench_compile` within gate (and now gates
that threading never makes a body slower); `leaks` 0 bytes; ThreadSanitizer
reports no races.

---

# Compile M3c — indexed machine arrays (2026-07-27)

Motivated by a request for a `Compile[]` tutorial built on a 2-D finite-difference
wave-equation solver. The solver was not writable: `Part` was outside the
compilable subset, so `u[[i, j]]` put the whole body on the interpreter.

- [x] `Part` reads, full spec vocabulary — inline path for one scalar subscript
      per axis (`A_AXIS` + existing `A_LOAD`), delegated path (`A_PART` →
      `ndarray_part`) for Span / All / position lists / partial indexing / any
      mixture.
- [x] `Part` assignment, same vocabulary, plus `+=` `-=` `*=` `/=` on a scalar
      position. Target must be an owned array (a `Module`/`With` local); writing
      through a borrowed argument is refused at compile time.
- [x] `ConstantArray` at any rank (`A_NEW`); array-typed `Module`/`With` locals,
      including as the result (`A_COPY`, `A_XFER`).
- [x] `Length` at any rank; multi-iterator `Do`.
- [x] `A_LOAD` made impure — a pure load is CSE'd across a store and hoisted out
      of the loop that mutates it. Both regression-tested.
- [x] Interpreter bugs found by the parity tests: `TimesBy` had no
      implementation at all (and `*=` / `/=` did not parse); `Part` assignment
      into an `NDArray` silently ignored every non-integer spec.
- [x] `COMPILE_EXAMPLE.md` — the tutorial, with an exact-discrete-solution
      correctness check and measurements against the interpreter, Wolfram
      Language 14.0 (WVM and native C) and `NDSolve` on both systems.

## Review

**The interesting number is not the 569x.** Compiled-vs-interpreted ratios say as
much about the interpreter as the compiler, and Mathilda's interpreter is ~12x
slower than Wolfram's on this array-heavy code. The number that isolates the
compiler is the cross-system one: Mathilda's compiled stencil is **1.8-2.1x
faster than Wolfram Language 14.0's**, stable across n = 41 to 641.

**WL's native-C target is slower than its own bytecode VM here** (14.9 s vs
12.2 s at n = 401), verified to be genuine `LibraryFunction` code rather than a
silent fallback. Worth weighing before this project's own "native backend next"
plan: in a tensor-heavy kernel the cost is array element access, not dispatch.

**At matched accuracy `NDSolve` beats the hand-written scheme by ~12x** (0.031 s
for 5.4e-5, versus 0.357 s at n = 153 for 1.57e-5). That is the honest framing
for `Compile[]`: it does not make your scheme competitive with a library solver
on the library solver's home ground; it makes the schemes that are *not* in the
library run at machine speed.

## New, from this work

- [ ] **`Module` costs 3x in the interpreter.** The identical interpreted march
      is 4.09 s at top level and 12.28 s inside a `Module` (n = 41). Not a
      compiler issue; worth its own look, and worth knowing when benchmarking —
      quoting the Module form would have inflated the speedup to 1700x.
- [x] **`Table` as an array constructor inside a compiled body — DONE (M7).**
      `A_NEW` sized from `max(0, floor((hi-lo)/di) + 1)` per axis (OP_QUOT_I is
      floor division, which is what the interpreter's `val <= hi` walk amounts
      to for an integer step), then k nested loops sharing one flat store index.
      INTEGER iterators only and a non-integer-valued body only — a real
      iterator is walked by repeated addition against a `1e-14` slack, and a
      packed buffer has no integer dtype.
- [ ] **Array-valued `If` branches / `Sum` accumulators / `Nest` state.** Each
      needs either an `OP_ARR_COPY` at the join or handle refcounting.
      (`OP_A_COPY` exists since M3c, so this is now a matter of using it.)

## M7 — the functional heads (2026-07-29)

- [x] **Compile-time function values** — `fn_resolve` / `emit_apply` /
      `infer_apply` replace `extract_function`. `Function[u,body]`,
      `Function[{u,...},body]`, `Function[body]` with `#`, a bare head, a
      `CompiledFunction` symbol, `Composition`, `Identity`. `Slot[k]` gets its
      own binding frame (it is not a symbol). Retargeting `Nest` gave it every
      spelling with no lowering of its own.
- [x] **`Fold`, `FixedPoint` (n / SameTest), `NestWhile`** — scalar-result
      iteration. `OP_SAMEQ_R`/`_C` (NaN is SameQ to NaN, per `expr_eq`) and
      `OP_FAIL` + `VM_ITER_SAFETY_CAP` matching the interpreter's 10^6.
- [x] **`Table`** (above), **`Map`**, **`Scan`**. `Map` takes the fused route
      when the body threads elementwise (6.6x over the element loop) and the
      general loop otherwise.
- [x] **Delegated structural heads** — `OP_A_NDFN` calls the interpreter's own
      `ndstruct_*` / `ndred_*` entry point, so `Reverse`, `Sort`, `Accumulate`,
      `Flatten`, `Transpose`, `Take[a,n]`, `Drop[a,n]` are parity by
      construction. Their degrade path returns a nested List, which is the
      signal to decline. `Total[Take[Sort[v], 3]]` is now 8 instructions.
- [x] **`NestList` / `FoldList`** — known-length history buffers.
- [x] **Three divergences fixed**: a BUILT array taking its kind from an
      unrelated argument; `compiled_eval_real` ignoring the abort flag;
      `Do`/`While`/`For` answering `0` instead of `Null` in result position.
- [x] **Interpreter prerequisite**: `Fold`/`FoldList` left an `NDArray`
      unevaluated, so there was nothing to be parity with.

### Still open from M7
- [ ] **`FixedPointList` / `NestWhileList`** — data-dependent length. Plan:
      compile only the LITERALLY BOUNDED forms (`FixedPointList[f, x, n]`,
      `NestWhileList[f, x, test, 1, max]`), allocating `n + 1` and truncating
      with a new `OP_A_TRUNC` (shrink `dims[0]` + realloc down, no copy; `k == 0`
      must give a length-0 array so it unpacks to `{}`). The unbounded forms
      should keep bailing: allocating `ITER_SAFETY_CAP` elements per call is not
      an acceptable default, and that is a resource limit, not a semantic one.
- [ ] **The selection heads** (`Select`, `TakeWhile`, `LengthWhile`,
      `SelectFirst`, `AllTrue`/`AnyTrue`/`NoneTrue`, `Join`, `Differences`,
      `First`/`Last`/`Most`/`Rest`, `RotateLeft`/`RotateRight`, `Riffle`,
      `Partition`, `MapThread`) — blocked on the interpreter, which leaves every
      one of them UNEVALUATED on a packed argument, so there is nothing to be
      parity with. Give them NDArray paths first (valuable on its own: packed
      arrays become first-class for a dozen more heads), then the compiled
      lowerings are mechanical and reuse `OP_A_TRUNC` and `OP_A_NDFN`.
- [x] **Retarget `src/numloop.c` onto this engine — INVESTIGATED, and the
      premise is FALSE.  Do not do it, and do not delete numloop.**

      The plan said this engine "now dominates" the legacy double-only stack
      machine, so every call site should move to `autocompile` and numloop
      should go.  Measured on identical bodies (median of 3, `-O3`):

      | body | this engine | numloop | |
      |---|---|---|---|
      | `Nest[(u + 2./u)/2 &, 3., 2e6]` | 0.032 s | 0.081 s | engine **2.5x** |
      | `While`, 1e6 iterations | 0.015 s | 0.023 s | engine **1.6x** |
      | `Do[s += Sin[1. i], {i,1,1e6}]` | 0.026 s | 0.022 s | numloop **1.2x** |
      | `Do[s += 1. i, {i,1,1e6}]` | 0.0106 s | 0.0102 s | even |
      | `Map[u^2+1. &, list]`, 2e5 | 0.030 s | 0.028 s | numloop **1.1x** |

      Two reasons the domination does not hold.  **`Map` over a plain `List`**
      pays this engine's boundary cost — packing 200k `Expr` nodes into a buffer
      and unpacking them again — which numloop skips entirely by walking the
      List.  (Over a PACKED argument this engine is 277x; the comparison is
      about the argument kind, not the loop.)  And **`Do`** was losing on
      instruction count, not on dispatch: see the next entry.

      What numloop is actually worth today, measured with `MATHILDA_NO_NUMLOOP`:
      `Nest` 37x, `While` 40x, `Do` 35x, `Map` 5.2x, `Fold` 1.08x.  Retargeting
      wholesale would have made the two most common loop constructs slower.

      Still worth doing, separately and on evidence: route `Nest` (2.5x, and it
      would gain complex accumulators, bare heads and `Composition`) and `While`
      (1.6x) through `autocompile`, leaving `Do` and `Map` on numloop.  That is a
      split rule across two engines in one file, so it needs to be worth the
      complexity — which is a judgement to make with the numbers above in hand.

- [x] **`Do` now closes with `OP_LOOP` — 8 inner instructions down to 5.**
      The disassembly said it plainly: four of the eight instructions in
      `Do[s = s + 1. i, {i,1,n}]` were pure loop control (test, branch,
      increment, back-edge), while `OP_LOOP` — which does increment, test and
      branch in ONE, and which the fused array loops had been using all along —
      was not used by any counted loop.  A unit step now closes with it (bound
      register holds `hi + 1`, since `OP_LOOP` tests `++i < a`), keeping the
      entry guard because a bottom-tested loop would otherwise run an empty
      range once.  Non-unit and negative steps keep the general form: `OP_LOOP`
      steps by one.  1.21x on the plain loop, 1.23x with a libm body.
      `Sum`/`Product` have the identical shape and should get the same
      treatment.
- [ ] **`MapIndexed`** — the index reaches `f` as a LIST (`{i}`), which the
      lattice cannot hold. Only worth doing under a restriction such as "the
      index appears solely inside `Part[idx, 1]`".
- [ ] **`Array`** — not simply `Table[f[i], {i,1,n}]`: it carries origin and
      range specs (`src/list/array.c:30`) that the rewrite would have to honour.
- [ ] **`x^0` is exact `1` in the interpreter and `1.` compiled.** Found via
      `Table[x^i, {i, 0, 4}]`, but PRE-EXISTING at the scalar level —
      `Compile[{x}, x^0][2.]` shows it alone. Not a one-line fix: `0.^0` is
      `Indeterminate`, so the result type genuinely depends on the value.
      `parity_head` has no zero-exponent case, which is why it survived.
- [ ] **Span assignment on a `List`** goes through `expr_part_assign_rec`, which
      rebuilds the structure per element — O(n) per assignment. The NDArray path
      is now O(selected); the List path is untouched.
- [ ] `Increment`/`Decrement` on a `Part` target in compiled code (only the
      binary forms handle `Part` today).
