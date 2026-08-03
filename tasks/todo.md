# Todo: Compile lowering for COMPILE_MISSING.md §4 (Fourier) + §5 (matrix producers)

Plan: /Users/user/.claude/plans/let-s-continue-with-out-quirky-eclipse.md

## Part A — §4 Fourier lowering
- [x] `src/fourier.c`: add `force_complex` param to `machine_build_ndarray`; update 4 callers
- [x] `src/fourier.c`: add `fourier_compile_core` + `fourier_compile` / `inverse_fourier_compile`
- [x] `src/fourier.h`: declare the two compile entry points
- [x] `src/compile/compile.c`: `#include ../fourier.h`; `NdFnSpec.complex_result`; `nd_fn_result` elem line; 4 Fourier `ND_FNS` rows

## Part B — §5 producer lowering
- [x] `src/compile/compile.c`: `nd_fn_result` `rank_rule 5`; 4 producer `ND_FNS` rows

## Part C — §5 buffer producers (interpreter efficiency)
- [x] `src/linalg/hankelmat.c`: `ndbuild_open` buffer branch in `hk_build`
- [x] `src/linalg/toeplitzmat.c`: buffer branch in `tz_build`
- [x] `src/linalg/vandermondemat.c`: buffer branch in `vm_build` (checked int power, overflow fallback)

## Part D — audits, docs, changelog
- [x] `tools/compile_coverage.py`: drop 8 heads from BASELINE (+ record pre-existing `Boole` gap)
- [x] `COMPILE_MISSING.md`: close §4/§5, update count (35→32/239), record "no new opcode" insight
- [x] `docs/spec/changelog/2026-08-03.md`: changelog entry
- [x] `docs/spec/builtins/{packed-arrays,fourier-transforms,linear-algebra}.md` + `docs/design/compile.md` §8a.1

## Part E — tests
- [x] `tests/test_compile_transforms.c` (new) + `tests/CMakeLists.txt`
- [x] producer packed-output correctness + Vandermonde overflow fallback

## Verification
- [x] `make` clean (C99); `make check-c99` green
- [x] `make check-compile-coverage` green (8 target heads no longer gaps; `Boole` recorded)
- [x] `make check-packed-aware` green
- [x] 13 affected test suites pass, 0 FAILs (incl. new compile_transforms, existing producer/fourier/compile suites)
- [x] valgrind: leak profile byte-identical to a pre-existing suite (macOS objc baseline) — zero new leaks
- [x] `make check-array-exactness` green (343 probes, **0 MIXED**)

## Review

**What shipped.** The `COMPILE_MISSING.md §4` transforms (`Fourier`,
`InverseFourier`, `FourierDCT`, `FourierDST`) and `§5` matrix producers
(`DiagonalMatrix`, `HankelMatrix`, `ToeplitzMatrix`, `VandermondeMatrix`) now
lower inside `Compile[]` and auto-compilation, closing their cliffs. The three
non-`DiagonalMatrix` producers also gained a direct rank-2 machine-buffer output
(REPL win too).

**Key finding — no new opcode.** `A_NDFN` stores whatever NDArray the delegate
returns; the compiler tracks only `CT_ARRAY(elem, rank)` statically. So §4 is one
flag (`NdFnSpec.complex_result`, plus a Fourier wrapper that always builds
`NDT_COMPLEX64`) and §5 is one rank rule (`rank_rule 5`), both in `nd_fn_result`.
The speculative "producing opcode" / "result-dtype field" the doc called for were
unnecessary.

**Two load-bearing subtleties.**
1. *Fourier commits to complex.* `machine_build_ndarray` collapses complex→real
   at run time by tolerance; a static register can't track that, so the compiled
   wrapper never collapses. Values match; only a tiny-nonzero-imaginary case
   differs (roundoff, not pinned by tests).
2. *Vandermonde does not coerce.* Its `Power[node,e]` cells keep the node's head,
   so a mixed int+real list stays unpacked (two heads). The buffer path is gated
   to a *uniform* node head — stricter than Hankel/Toeplitz, which coerce via
   `hk_cell`/`tz_cell`. Integer powers use `ci_powi_i64` and fall back to the
   exact bignum path on overflow.

**Also fixed (pre-existing).** `Boole` was on `AWARE` (from the bool-dtype
commit) but neither lowered nor in `BASELINE`, so `check-compile-coverage` was
already red before this work; recorded it in `BASELINE`.

**Verification.** Clean `-std=c99 -O3` build; `check-c99`, `check-compile-coverage`
(8 heads no longer gaps), `check-packed-aware` all green. 13 affected test suites
pass (new `test_compile_transforms.c` + existing producer/fourier/compile suites).
Valgrind leak profile byte-identical to a pre-existing suite → zero new leaks.
`check-array-exactness` run separately.

**Deferred (documented).** ~~The two-vector `Hankel`/`Toeplitz[c, r]` forms
(rank-1 × rank-1 → rank-2) want an `A_NDFN2` lowering.~~ **Done** (follow-up).

## Follow-up — two-vector Hankel/Toeplitz A_NDFN2 lowering (2026-08-03)

- `src/compile/compile.c`: new `R2_MATRIX` rank rule (rank 1 × rank 1 → rank 2)
  in `NdFn2Rank` + `nd_fn2_result`; two `ND_FN2S` rows for `HankelMatrix` /
  `ToeplitzMatrix` (`NDF_INT|NDF_REAL`, `NDF2_SAME`), delegating to the builtins
  (which delist the two arrays and write the buffer directly via `hk_build`/
  `tz_build`, so the O(m·n) output is zero-Expr).
- int+int → int64, real+real → real; mixed int/real, complex, or a rank-2
  operand decline cleanly to the interpreter (which coerces exactly).
- Tests: `test_cf_hankelmatrix_2vec` / `test_cf_toeplitzmatrix_2vec` in
  `tests/test_compile_transforms.c` (lowers, real/int parity, rectangular,
  composed `Tr`, clean declines) + a two-vector case in the leak loop.
- Verified: build clean, `check-c99`, `check-compile-coverage` green; 8 compile/
  producer suites pass (0 FAILs); valgrind leak profile byte-identical → zero new
  leaks. Docs updated (`COMPILE_MISSING.md` §5, changelog, `packed-arrays.md`,
  `linear-algebra.md`).
