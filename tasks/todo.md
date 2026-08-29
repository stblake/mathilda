# Task: Implement SchurDecomposition

Plan: `/Users/user/.claude/plans/following-on-from-our-tingly-plum.md`

## Checklist

### 1. LAPACK Schur drivers
- [ ] `lapack.h`: prototypes `dgees_`/`zgees_`/`dgges_`/`zgges_`/`dgebal_`/`zgebal_` (non-Accelerate block)
- [ ] `lapack.h`: `mat_lapack_*` wrapper decls (+ stub-safe)
- [ ] `lapack.c`: wrapper impls (query→alloc→call→free) + no-LAPACK stubs

### 2. Expose MPFR Schur Q
- [ ] `eigen.h`: decl `eigen_schur_real_mpfr` (USE_MPFR)
- [ ] `eigen_direct.c`: impl reusing Hessenberg + Francis QR, return Q + T

### 3. Schur module
- [ ] `schurdecomp_internal.h` (opts struct, dispatch decls)
- [ ] `schurdecomp.h` (builtin + init decls)
- [ ] `schurdecomp.c` (entry, options parse, std/gen detection, precision dispatch, result build, ndla)
- [ ] `schurdecomp_machine.c` (dgees/zgees/dgges/zgges + dgebal pivoting + RBDF->False complexify)
- [ ] `schurdecomp_mpfr.c` (standard real via eigen_schur_real_mpfr)

### 4. Registration
- [ ] `sym_names.{h,c}`: `SYM_RealBlockDiagonalForm`
- [ ] `core.c`: `schurdecomp_init()` call
- [ ] `info.c`: docstring
- [ ] `schurdecomp_init`: builtin + ATTR_PROTECTED

### 5. Surfaces
- [ ] `pack.c`: AWARE entry (not INT64_OK)
- [ ] `tools/compile_coverage.py`: EXEMPT entry

### 6. Tests
- [ ] `tests/test_schurdecomp.c`
- [ ] `tests/CMakeLists.txt`: COMMON_SRC + executable stanza

### 7. Docs / version
- [ ] `docs/spec/builtins/linear-algebra.md`
- [ ] `docs/spec/changelog/2026-08-24.md`
- [ ] `src/version.h`: 0.117 → 0.118

### 8. Verify
- [x] build + check-c99 (both clean)
- [ ] unit tests (schur + regression: jordan/eigen/svd) — building
- [x] golden examples (piped transcript) — all pass, 100x100 recon 4.5e-14
- [x] packed audits: packed-aware OK, nd-surfaces OK, array-exactness OK (0 MIXED);
      compile-coverage: Schur correctly EXEMPT (22 pre-existing unrelated failures)
- [ ] valgrind
- [ ] rebuild graph + self-review

## Notes
- All sections 1-7 complete. LAPACK dgees/zgees/dgges/zgges + dgebal/dgebak added.
- check-compile-coverage has 22 PRE-EXISTING failures (image/fit/interp/Hermite),
  none from this change; Schur is exempt and not flagged.

## Review

**Done (v0.118).** `SchurDecomposition[m]` (standard) and `SchurDecomposition[{m,a}]`
(generalized/QZ) implemented across `src/linalg/schurdecomp{,_internal}.h`,
`schurdecomp.c` (dispatch), `schurdecomp_machine.c` (LAPACK), `schurdecomp_mpfr.c`
(arbitrary-precision standard real).

- **LAPACK layer**: new `mat_lapack_dgees/zgees/dgges/zgges` + `dgebal/zgebal/dgebak/zgebak`
  wrappers in `lapack.{c,h}` (Fortran prototypes guarded out under Accelerate; no-LAPACK stubs).
- **MPFR**: new `eigen_schur_real_mpfr` in `eigen_direct.c` surfaces the orthogonal Q the
  Francis QR already accumulates.
- **Options**: `Pivoting` ({q,t,d}, m.d==d.q.t.q^H via dgebal/dgebak), `RealBlockDiagonalForm`
  (real blocks vs complex triangular), `TargetStructure->"Structured"` (returns dense).
- **Surfaces**: NDArray/packed handled via na_load_matrix (on pack.c AWARE); result matches
  input representation (List→List, NDArray→NDArray) to avoid Plus mis-threading; symbolic
  constants (Pi) numericalised on load miss; Compile-exempt.
- **Registration**: `schurdecomp_init` in core.c, ATTR_PROTECTED, docstring, `SYM_RealBlockDiagonalForm`.
- **Docs**: linear-algebra.md section, changelog 2026-08-24.md, version 0.117→0.118.

**Verification**: 24/24 unit tests pass; golden transcript reconstructs to ~1e-14 (100×100 in
a few ms); check-c99 / check-packed-aware / check-nd-surfaces / check-array-exactness pass;
Schur correctly EXEMPT in check-compile-coverage (22 pre-existing unrelated failures untouched);
valgrind clean (only Accelerate's own LAPACK baseline leak, identical to Eigenvalues/SVD);
JordanDecomposition + eigen regression suites pass.

**Independent code review**: no blockers/leaks/ABI errors. Two warnings applied as fixes —
(1) check `dgebak`/`zgebak` return and gate `ok=0` (was discarded, inconsistent with `dgebal`);
(2) docstring caveat that MPFR is native only for the standard-real default-options case.

**Known limitations (documented)**: complex/generalized/RBDF-False/Pivoting at arbitrary
precision fall back to machine precision (would need a complex MPFR QR / an MPFR QZ).
