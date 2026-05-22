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
 * to the appropriate kernel.  Returns a freshly-allocated Expr* on
 * success or NULL on cannot-evaluate (caller leaves the call
 * unevaluated; never frees `m`).
 * ------------------------------------------------------------------ */
Expr* lu_dispatch(Expr* m, int n);

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
Expr* lu_symbolic_dispatch(Expr* m, int n);

/* ---------------------------------------------------------------------
 * Symbolic Doolittle core.
 *
 * Pure C function over Expr** buffers; drives every primitive through
 * the Mathilda evaluator so it works on integer / rational / Complex /
 * Sqrt-bearing / free-symbolic entries with no special casing.
 *
 * Inputs:
 *   A_flat[i*n + k]  =  m[i, k],   i in [0, n),  k in [0, n)
 *
 * Outputs (all caller-owned, freshly allocated):
 *   *out_LU_flat   n * n entries row-major   (LU_flat[i*n + k])
 *                  combined Doolittle L (strict-lower, unit diag) and
 *                  U (upper).
 *   *out_perm      length n, 1-indexed row permutation
 *                  (perm[k] = original-row used at step k)
 *   *out_singular  true if a zero pivot was encountered at any step
 *                  (factorisation completes regardless, matching
 *                  Mathematica's LUDecomposition::sing behaviour).
 *
 * Returns true on success.  A_flat is NOT consumed.
 * ------------------------------------------------------------------ */
bool lu_symbolic_core(Expr** A_flat, int n,
                      Expr*** out_LU_flat, int** out_perm,
                      bool* out_singular);

/* ---------------------------------------------------------------------
 * Machine-precision kernel dispatcher (LAPACK).
 *
 * Invoked by lu_dispatch when common_scan_inexact reports an inexact
 * input at min_bits <= 53.  Loads `m` into a column-major double
 * buffer, calls dgetrf (or zgetrf for complex), then dgecon (zgecon)
 * for the L-infinity condition number.
 *
 * Returns NULL (without consuming `m`) when:
 *   - USE_LAPACK is undefined.
 *   - A matrix leaf isn't a recognised numeric value.
 *   - LAPACK reports a fatal info code.
 * ------------------------------------------------------------------ */
Expr* lu_machine_dispatch(Expr* m, int n);

/* ---------------------------------------------------------------------
 * Arbitrary-precision MPFR kernel dispatcher.
 *
 * Invoked by lu_dispatch when min_bits > 53.  Runs Doolittle with
 * partial pivoting over column-major MPFR arrays at the input's
 * working precision; computes the L-infinity condition number from the
 * explicit inverse (back-substitution via L then U).
 *
 * Returns NULL (without consuming `m`) when:
 *   - USE_MPFR is undefined.
 *   - A matrix cell can't be reduced to an MPFR value.
 * ------------------------------------------------------------------ */
Expr* lu_mpfr_dispatch(Expr* m, int n);

#endif /* LUDECOMP_INTERNAL_H */
