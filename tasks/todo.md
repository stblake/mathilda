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

## Next, in value order

- [ ] **Block strip-mining for fusion** — each opcode processes a tile of ~64
      elements in a vectorisable C loop; dispatch amortises 64x and temporaries
      stay in L1. This is where the array order of magnitude is. Fusion stays
      off until this lands.
- [ ] **Make CSE actually fire** — it is defeated structurally: `binop`/`unop`
      write into an operand's register, invalidating the value-number entry the
      same instruction created. Either emit in SSA form + linear-scan regalloc,
      or (much cheaper) do CSE at the `Expr` level in `emit`, hoisting
      structurally-equal subtrees into persistent registers like `With` locals.
- [ ] **`OP_CALL`** — compiled-to-compiled calls without the inline depth cap of
      8, which is what recursion needs. Requires the frame stack (below).
- [ ] **Frame stack** — replace the single per-program `p->frame`, which makes a
      `CompiledProgram` non-reentrant and non-thread-safe today. Prerequisite for
      `OP_CALL` and for threading a strip-mined loop.
- [ ] **Fill the remaining kernels** — 48 listed gaps. Most are MPFR-only modules
      (`Zeta`, `PolyLog`, `Airy*`, …) needing genuinely new double implementations,
      each with a parity test against the MPFR path. Trivial tier first
      (`Sinc`, `UnitStep`, `Fibonacci`, `InverseErf`/`Erfc` — the last two already
      have double kernels, just unregistered).
- [ ] **User `Compile[]` array argspec** — `{v, _Real, 1}`; `cf_box`/`cf_unbox`
      still have `default: break` for array types.
- [ ] **`CompileDiag`** — a bail still reports nothing. The audit covers heads;
      per-body diagnostics would cover the rest.
- [ ] **Native backend** (`CompilationTarget -> "C"`), behind a build flag.
