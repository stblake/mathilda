/*
 * src/linalg/ludecomp_internal.h
 *
 * Internal interface for the LUDecomposition translation units:
 *
 *   ludecomp.c          -- builtin entry + lu_dispatch (kernel router)
 *   ludecomp.c           also hosts lu_symbolic_dispatch + lu_symbolic_core,
 *                        the rationalise -> exact Doolittle -> numericalise
 *                        pipeline used for symbolic, integer, rational,
 *                        complex, machine-real, and MPFR-real inputs.
 *   ludecomp_machine.c  -- LAPACK kernel (real + complex)
 *   ludecomp_mpfr.c     -- MPFR Doolittle kernel
 *
 * Public surface is `Expr* builtin_ludecomposition(Expr*)` in
 * ludecomp.h.
 *
 * Architecture mirrors qrdecomp_internal.h: a thin top-level
 * dispatcher classifies the input by leaf precision and routes to the
 * appropriate kernel; each kernel owns its own load -> work -> store
 * round-trip and returns NULL on soft failure so the dispatcher can
 * fall back to the symbolic path.
 */

#ifndef LUDECOMP_INTERNAL_H
#define LUDECOMP_INTERNAL_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "expr.h"

/* ---------------------------------------------------------------------
 * Top-level dispatcher.
 *
 * Classifies `m` by leaf precision via common_scan_inexact and routes
 * to the appropriate kernel.  Accepts any `rows x cols` shape with
 * `rows > 0` and `cols > 0`; the factorisation has `perm` length
 * `min(rows, cols)` (matching Mathematica's contract for rectangular
 * input).  Returns a freshly-allocated Expr* on success or NULL on
 * cannot-evaluate (caller leaves the call unevaluated; never frees
 * `m`).
 * ------------------------------------------------------------------ */
Expr* lu_dispatch(Expr* m, int rows, int cols);

/* ---------------------------------------------------------------------
 * Symbolic kernel dispatcher.
 *
 * Wraps the Doolittle core with:
 *   - inexact-input rationalisation at min(input bits)
 *   - flatten -> work -> wrap conversion to/from List-of-List
 *   - tidy_matrix (element-wise Together) on lu
 *   - numericalisation back to input precision when rationalisation ran
 *   - assembly of {lu, p, c}
 * ------------------------------------------------------------------ */
Expr* lu_symbolic_dispatch(Expr* m, int rows, int cols);

/* ---------------------------------------------------------------------
 * Symbolic Doolittle core.
 *
 * Pure C function over Expr** buffers; drives every primitive through
 * the Mathilda evaluator so it works on integer / rational / Complex /
 * Sqrt-bearing / free-symbolic entries with no special casing.
 *
 * Inputs:
 *   A_flat[i*cols + k]  =  m[i, k],   i in [0, rows), k in [0, cols)
 *
 * Outputs (all caller-owned, freshly allocated):
 *   *out_LU_flat   rows * cols entries row-major
 *                  (LU_flat[i*cols + k]).  Combined Doolittle L
 *                  (strict-lower, unit diag) and U (upper).
 *   *out_perm      length rows, 1-indexed row permutation
 *                  (perm[i] = original-row that ended up at row i).
 *                  Rows past min(rows, cols) - 1 are never swapped
 *                  during elimination, so for tall input the trailing
 *                  entries remain at their identity values.
 *   *out_singular  true if a zero pivot was encountered at any step
 *                  (factorisation completes regardless, matching
 *                  Mathematica's LUDecomposition::sing behaviour).
 *
 * Returns true on success.  A_flat is NOT consumed.
 * ------------------------------------------------------------------ */
bool lu_symbolic_core(Expr** A_flat, int rows, int cols,
                      Expr*** out_LU_flat, int** out_perm,
                      bool* out_singular);

/* ---------------------------------------------------------------------
 * Machine-precision kernel dispatcher (LAPACK).
 *
 * Invoked by lu_dispatch when common_scan_inexact reports an inexact
 * input at min_bits <= 53.  Loads `m` into a column-major double
 * buffer, calls dgetrf (or zgetrf for complex), then dgecon (zgecon)
 * for the L-infinity condition number.  Non-square input is accepted;
 * the condition slot is set to exact Integer 0 in that case (the
 * estimate is only meaningful for square A).
 *
 * Returns NULL (without consuming `m`) when:
 *   - USE_LAPACK is undefined.
 *   - A matrix leaf isn't a recognised numeric value.
 *   - LAPACK reports a fatal info code.
 * ------------------------------------------------------------------ */
Expr* lu_machine_dispatch(Expr* m, int rows, int cols);

/* Packed fast path for a buffered matrix input (the ndla_* pattern, see
 * ndlinalg.h): reads the NDArray buffer in place instead of delisting it
 * to a nested List and re-evaluating -- the delist dominated the cost.
 * Square real float64 only; anything else (rectangular, complex, no
 * LAPACK) falls back to linalg_delist_and_reeval. */
Expr* ndla_ludecomposition(Expr* res);

/* ---------------------------------------------------------------------
 * Arbitrary-precision MPFR kernel dispatcher.
 *
 * Invoked by lu_dispatch when min_bits > 53.  Runs Doolittle with
 * partial pivoting over row-major MPFR arrays at the input's working
 * precision; for square input also computes the L-infinity condition
 * number from the explicit inverse (back-substitution via L then U).
 * For non-square input the condition slot is set to exact Integer 0.
 *
 * Returns NULL (without consuming `m`) when:
 *   - USE_MPFR is undefined.
 *   - A matrix cell can't be reduced to an MPFR value.
 * ------------------------------------------------------------------ */
Expr* lu_mpfr_dispatch(Expr* m, int rows, int cols);

/* ---------------------------------------------------------------------
 * Arbitrary-precision determinant.
 *
 * det(A) = sign(P) * prod(U[k,k]) for an n x n inexact-numeric matrix,
 * computed via the same MPFR Doolittle factorisation as lu_mpfr_dispatch.
 * Used by builtin_det to replace the O(n!) Laplace fallback for inexact
 * input (which hangs at n >= ~12) with an O(n^3) kernel whose wide MPFR
 * exponent never spuriously overflows.
 *
 * Returns an EXPR_MPFR (real) or Complex[mpfr, mpfr] (complex), or NULL
 * (without consuming `m`) when USE_MPFR is undefined, the input is exact,
 * or a matrix cell can't be reduced to an MPFR value.
 * ------------------------------------------------------------------ */
Expr* mpfr_det_dispatch(Expr* m, int n);

#endif /* LUDECOMP_INTERNAL_H */
