# Task: Matrix decompositions stay on packed/NDArray buffers (no Expr boxing)

Plan: `/Users/user/.claude/plans/following-on-from-our-tingly-plum.md`

## Phase 0 — NDArray[] idempotency (prerequisite)
- [ ] `src/ndarray.c` `builtin_ndarray`: idempotent `is_ndarray(arg)` branch (mirror `builtin_tondarray`)
- [ ] verify `NDArray[NDArray[m]] === NDArray[m]`, dtype re-cast works

## Phase 1 — shared helpers
- [ ] `numarray.{c,h}`: `na_result_presentation(input)` + `na_build_matrix_as(...)`
- [ ] refresh stale "NDArray real-only" header prose

## Phase 2 — SchurDecomposition
- [ ] `schur_load_cm`: real load first (fast path), complex fallback
- [ ] `schur_build`: drop pack_unpack; `na_build_matrix_as` + presentation

## Phase 3 — pack complex RESULTS (Eigen/Jordan)
- [ ] `eigen_direct.c`: complex eigenvalue/vector builders → COMPLEX64 packed
- [ ] `jordandecomp.c`: complex spectrum on packed input → COMPLEX64 packed

## Phase 4 — complex NDArray INPUT on buffer (reduce delist)
- [ ] LU → zgetrf; SVD → zgesdd; Eigen → zgeev; QR → cplx=1
- [ ] scope: complex-dtype decline only (options/generalized may keep delisting)

## Phase 5 — presentation consistency (boxed → packed-list)
- [ ] QR boxed-input branch → packed-list stamp
- [ ] audit all heads: boxed input → transparent packed-list, no boxing

## Phase 6 — docs, tests, benchmark
- [ ] docstrings/spec + fix stale "complex boxed" comments
- [ ] tests: packed real+complex I/O, 3 input forms, idempotency
- [ ] regression + audits + valgrind
- [ ] scipy baseline on packed inputs (≈1.0x); formalize benchmarks/30-schur-decomposition/
- [ ] changelog + version bump

## Review

### Milestone committed (v0.119): mechanism + Schur + QR + NDArray[] fix
- **Phase 0 ✓** `NDArray[]` idempotent on an NDArray arg (was malformed `{1}`);
  EXEMPT in check_packed_aware (idempotency guard, not a fast path).
- **Phase 1 ✓** `na_result_presentation` + `na_build_matrix_as`/`_vector_as` in
  numarray. Two subtle correctness fixes discovered here: (a) arm
  `pack_g_any_created` for NDA_HEAD_LIST results so the gate normalizes them in
  reconstruction arithmetic (else `m - q.t.q^H` mis-threads); (b) a complex
  result for a boxed-List caller must be boxed to a Complex[] List, NOT a
  complex packed-list (the gate corrupts complex packed-lists → float64).
- **Phase 2 ✓** Schur: real-first load (fast memcpy path) + packed factors.
  Packed 400×400 ~65ms, on par with / faster than scipy ~87ms.
- **Phase 5 (partial) ✓** QR boxed-input → transparent packed-list (was boxed).
- Verified: all input forms reconstruct ~1e-14; packed-aware/array-exactness
  audits pass; schur/jordan/eigen/packed_list/ndarray_linalg tests pass;
  valgrind clean (Accelerate baseline only).

### Phase 4 ✓ (v0.120): complex NDArray input on the buffer
- LU → zgetrf (+ zlange/zgecon); SVD → zgesdd (u,v complex64, sigma real,
  v=ConjugateTranspose[VT]); Eigenvalues/Eigenvectors → zgeev. All return packed
  complex64; verified vs numpy (evals/svals match, recon ~1e-15); valgrind clean.

### Phase 3 — DECIDED NOT TO DO
- A real matrix's mixed real+complex spectrum stays a boxed List: force-packing
  it to complex64 would turn real eigenvalues into Complex[re,0], diverging from
  Mathematica, and the boxed result already agrees across surfaces. The input is
  never unpacked (na_load reads the buffer); only the inherently-mixed result is
  boxed. Complex INPUT (uniformly complex result) is handled by Phase 4.

### Remaining
- **Phase 6** benchmark folder `benchmarks/30-schur-decomposition/`.
