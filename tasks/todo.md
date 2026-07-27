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
- [ ] **Thread the strip-mined loop** via `nd_parallel_for` — now unblocked by
      C-stack frames. Note `-DMATHILDA_THREADS` is not set in tests/CMakeLists.
- [ ] **Fill the remaining kernels** — 39 listed gaps (was 48). The
      exponential-integral family is done (`src/special_functions/expint_machine.c`).
      Remaining are genuinely new numerics: `Zeta`, `PolyLog`, `HurwitzZeta`,
      `LerchPhi`, `PolyGamma`, `Erfi`, `Fresnel*`, `Airy*`, `ProductLog`,
      `BesselI/K`, `LegendreP`, the hypergeometrics. Plus `UnitStep`/`Clip`/
      `Rescale`, which need an n-ary kernel form (the registry is unary/binary).
- [ ] **`CompileDiag`** — a bail still reports nothing. The audit covers heads;
      per-body diagnostics would cover the rest.
- [ ] **Native backend** (`CompilationTarget -> "C"`), behind a build flag.
