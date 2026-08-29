# Packed/NDArray fast paths + packed output for all matrix decompositions

Bring every decomposition to the **LUDecomposition standard**: read the packed
buffer directly (`na_load_matrix`, no delist), and return packed output
(`na_build_matrix`/`na_build_vector`) inheriting the input's presentation
(`out->present_as = arg->present_as`). Real results only (na_build declines
complex → nested List, which is correct). Boxed-List input → boxed output.

## Heads
- [ ] **SingularValueDecomposition** — no fast path today (delists). Add
      pre-delist na_load + dgesdd kernel; packed {u, w, v} (all real). Biggest win.
- [ ] **QRDecomposition** — reads packed already, but `mat_qr_mathilda`
      *unpacks* the output (ndarray_to_nested_list). Make it inherit present_as → packed q, r.
- [ ] **Eigenvalues** — delists + re-buffers. Direct na_load read; packed
      eigenvalue vector for a real spectrum (complex → boxed).
- [ ] **Eigenvectors** — same; packed eigenvector matrix for a real spectrum.
- [ ] **JordanDecomposition** — reads packed (done); make s,j packed for a
      real spectrum (inherit present_as); complex spectrum stays boxed.

Not applicable: HermiteDecomposition (integer-only), CharacteristicPolynomial
(returns a polynomial).

## Per-head checklist
Build clean · jordandecomp/eigen/linalg/svd/qr/lu tests pass · packed==unpacked
agreement · valgrind clean · check-packed-aware / nd-surfaces / array-exactness /
c99 green.

## Review — DONE

All heads brought to the LU standard (read packed buffer directly; packed
output inheriting present_as; boxed-List input keeps boxed path):

- **SVD** — `ndla_singularvaluedecomposition` (dgesdd). Packed {u,w,v}. 24.8→16.3ms (200².
- **QR** — `mat_qr_mathilda` no longer unpacks; inherits present_as.
- **Eigenvalues/Eigenvectors** — `ndla_*` with structure dispatch: dsyev
  (symmetric→packed real) / dgeev (general→packed real, boxed complex). Fixed a
  would-be regression (dgeev on symmetric was 5× slower than dsyev).
- **JordanDecomposition** — packed s,j for a real spectrum.

Verified: 18 linalg suites pass; packed==NO_PACK agreement; valgrind-clean
(byte-identical to baseline across sym/general/SVD/QR/Jordan); c99, packed-aware,
nd-surfaces (EXIT 0), array-exactness (0 MIXED) all green; 0 build warnings.

Perf (200×200 packed vs boxed): SVD 16.3 vs 24.8ms, Eigenvalues(sym) 2.2 vs
2.7ms, Eigenvectors(sym) 14.3 vs 16.3ms — packed now faster everywhere, no
regressions.

Limits: packed OUTPUT only for real spectra (na_build has no complex machine
array); complex-spectrum results stay boxed, as before. Generalized forms,
truncation, options, and complex-entry matrices fall back to the boxed path.
