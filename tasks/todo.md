# Task: Robust type-crossing (machine ↔ MPFR ↔ exact) for Det / Eigenvalues / Eigenvectors / CharacteristicPolynomial

## Problem statement
`Det[RandomReal[{-10,10},{200,200}]]` returns `-inf.0`. The true determinant is
≈ −1.076×10³⁴⁰ (verified via exact rationalize→N), which genuinely exceeds IEEE
double range (max ≈1.8×10³⁰⁸). The machine LU path multiplies the U-diagonal in
raw `double`, so it overflows to `±inf` (and underflows to `0.0` mid-product)
even when a correct answer exists. Per the directive, an overflow of the machine
type must **defer to arbitrary precision** rather than return `inf`.

## Investigation results (current state)
- **Det machine real**: overflow→`inf`, underflow→`0.0` (non-singular). `ndla_det`
  `src/linalg/ndlinalg.c:243-256`. Uses LAPACK `dgetrf` ✓ but product is naive.
- **Det machine complex**: same class, `ndlinalg.c:267-278`.
- **Det genuine MPFR matrix (>53 bits)**: NO fast path. `builtin_det` uses FLINT
  (exact) or O(n!) Laplace → **8×8 ok, 16×16+ hangs**. (`det.c`). Note the MPFR
  LU kernel already exists (`ludecomp_mpfr.c`) and is unused by Det.
- **Norm vector 2-norm**: `s += mag*mag` overflows; `Norm[{1e200,1e200}]`→`inf`
  (true ≈1.41e200). `ndlinalg.c:608-614`. Scaled BLAS `dnrm2`/`dznrm2` already
  exist (`blas_bridge.c:78,158`) and are not used here.
- **Norm p-norm**: `pow(mag,p)` unscaled, same class. `ndlinalg.c:615-621`.
- **Normalize**: hand-rolled 2-norm then divide; `Normalize[{1e200,1e200}]`→
  **zero vector**. `ndlinalg.c:642-651`.
- **Cleared (NOT bugs)**: matrix Norm (LAPACK-scaled), MatrixPower/Tr (inherent
  growth, no representable answer lost), Dot (BLAS). No `Permanent` exists.
- **Eigenvalues/Eigenvectors**: machine paths use LAPACK `dgeev`, fast & finite
  for 200×200; MPFR/exact via Faddeev (O(n⁴)) fast at 16×16. Healthy.
- **CharacteristicPolynomial**: Faddeev, fast even at 120×120. Returns symbolic
  `Plus`. Healthy.
- **Compile**: `Det` compiles (ND_FNS, AWARE). `Eigenvalues`/`Eigenvectors` are
  correctly EXEMPT (result type/shape data-dependent — cannot be a statically
  typed compiled register; same limitation Mathematica has). `CharacteristicPoly`
  returns a symbolic Plus → not compilable, but is **unrecorded** in the coverage
  tool's EXEMPT list → contributes to `make check-compile-coverage` being red.
  (Note: the gate is *already* red on main from 22 other unrelated AWARE heads —
  ArrayPlot, Binarize, Image*, Fit, …; those are out of scope.)

## Plan

### Phase 1 — Det machine overflow/underflow → MPFR fallback (the reported bug)
1. Add helper `nd_real_det_result(const double* LU, int n, int sign)` in
   `ndlinalg.c`: direct `double` product; if any diagonal is exactly 0.0 →
   singular → `expr_new_real(0.0)`; if the product is finite & nonzero → return
   `expr_new_real`; otherwise (overflow/underflow) re-accumulate the same
   diagonal in a 53-bit-mantissa `mpfr_t` (wide exponent) → `expr_new_mpfr_move`.
   `#ifndef USE_MPFR` degrades to the current `double` result.
2. Refactor `ndla_det` real branch (both LAPACK & hand-LU) to compute `sign`
   then call the helper.
3. Add complex analogue `nd_complex_det_result` → `Complex[mpfr,mpfr]` on
   overflow; refactor the complex branch to use it.

### Phase 2 — Det on genuine MPFR matrices (kill the O(n!) Laplace hang)
4. Add `Expr* mpfr_det_dispatch(Expr* m, int n)` in `ludecomp_mpfr.c` (declared
   in `ludecomp_internal.h`): reuse `lum_load_matrix`/`lum_factor`; multiply the
   U-diagonal in MPFR; fold the permutation-parity sign; real→`EXPR_MPFR`,
   complex→`Complex[mpfr,mpfr]`. Returns NULL on non-numeric leaf / USE_MPFR off.
5. Wire into `builtin_det` (`det.c`) after the FLINT exact path and before the
   Laplace fallback: for inexact-numeric input, try `mpfr_det_dispatch`.
   Fixes the 16×16+ hang; leaves Laplace only for small symbolic/exact matrices.

### Phase 3 — Norm / Normalize overflow (same class, from survey)
6. `ndla_norm` vector 2-norm → route through scaled `BLAS`dnrm2`/`dznrm2`.
7. `ndla_norm` p-norm → scale by max|xᵢ| before accumulating.
8. `ndla_normalize` → use the scaled norm; fixes the zero-vector result.

### Phase 4 — Compile coverage bookkeeping
9. Record `CharacteristicPolynomial` as EXEMPT in `tools/compile_coverage.py`
   (returns a symbolic `Plus`, no machine-buffer result). Re-confirm Det/Eig/
   Evec entries. Report (but do not fix) the 22 pre-existing unrelated red heads.

### Phase 5 — Verify, document, remember
10. Add C tests (or a `.m` corpus check) for: Det 200×200 overflow → finite MPFR
    matching exact ground truth; Det underflow-diagonal → 1e200; Det MPFR 16×16
    fast & correct; Norm/Normalize 1e200 vectors. Rebuild + run.
11. Update `docs/spec/builtins/` (linalg) + this week's changelog
    (`docs/spec/changelog/2026-08-24.md`).
12. Update harness memory (new: linalg overflow→MPFR-fallback pattern).

## Review (completed 2026-08-26)

All phases done and verified.

- **Det machine overflow/underflow** → `nd_real_det_result` / `nd_complex_det_result`
  (`ndlinalg.c`). `Det[RandomReal[{-10,10},{200,200}]]` = `-1.076e340` (was `inf`,
  rel. err 1.6e-14 vs exact); underflow diagonal → `1e200` (was `0`); complex
  overflow → finite; singular → `0`.
- **Det MPFR matrix** → `mpfr_det_dispatch` (`ludecomp_mpfr.c`), wired into
  `det.c`. 16×16 & 40×40 MPFR now instant (were hanging in Laplace).
- **Det exact-integer NDArray** (bonus, found during final sweep) → int64 buffers
  delist to the exact FLINT/bignum path in `ndla_det`; `Det[NDArray[{{1e9,1},
  {1,1e9}}, "int64"]]` = `999999999999999999` (was `1e18`, precision-lost).
- **Norm/Normalize overflow** → scaled `nd_vec_2norm` (dnrm2 recurrence) +
  scaled p-norm (`ndlinalg.c`); machine vectors routed there from `norm.c`.
  `Norm[{1e200,1e200}]`=`1.414e200` (was `inf`); `Normalize`→unit vector (was 0).
  Exact vectors keep `Sqrt[5]`; MPFR vectors keep full precision.
- **Compile** — `Det` compiles (unchanged); `CharacteristicPolynomial` recorded
  EXEMPT in `tools/compile_coverage.py`; Eigenvalues/Eigenvectors stay EXEMPT
  (result type/shape data-dependent). No forced (wrong) lowerings.
- **Eigenvalues/Eigenvectors** already used LAPACK `dgeev`; verified finite/fast
  at 200×200. No change needed.

Verification: 7 linalg test suites green (added `test_det_overflow` + Norm/
Normalize/int64 cases); `check-array-exactness` 0 MIXED; `check-nd-surfaces`
agreement clean (det/norm return scalars, no mismatch); `check-c99` clean;
valgrind shows no new leaks (byte-identical to a trivial-script baseline).

### Out of scope / flagged to user
- `make check-compile-coverage` is **already red on `main`** from 22 unrelated
  AWARE heads never recorded in EXEMPT/BASELINE (ArrayPlot, Binarize, Chop, Clip,
  Image*, Fit, GaussianFilter, Interpolation, …). Recording CharacteristicPolynomial
  removes it from the failing set but does not turn the gate green — the other 22
  are a pre-existing, separate bookkeeping regression.
- `nd-surfaces --survival` lists pre-existing producers that drop packing
  (qr, svd, eigenvalues/vectors, mask/struct ops) — the open array-substrate
  roadmap, untouched by this change.

## Non-goals / decisions
- **Do NOT force a Compile lowering for Eigenvalues/Eigenvectors/CharPoly** — their
  results are not statically-typed machine values; a lowering would be wrong.
  Record the exemptions instead.
- **Do NOT** re-add or modify BLAS/LAPACK plumbing — it is already used (dgetrf,
  dgeev, dnrm2, dlange/dgesdd).
- The overflow fallback uses a **53-bit mantissa** (input is machine precision);
  only the exponent range is widened. Honest, not falsely precise.
