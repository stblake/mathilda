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
- [ ] **Airy's uncovered band** — `2.5 < |x| < 8`, where neither the ascending
      series nor the asymptotic expansion reaches double precision. Needs a third
      method (Chebyshev fits, or the modified-Bessel route). Declines today.
- [ ] **`CompileDiag`** — a bail still reports nothing. The audit covers heads;
      per-body diagnostics would cover the rest.
- [ ] **Native backend** (`CompilationTarget -> "C"`), behind a build flag.

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
