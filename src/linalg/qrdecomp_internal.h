/*
 * src/linalg/qrdecomp_internal.h
 *
 * Internal interface for the QRDecomposition translation units.
 *
 *   qrdecomp.c          -- builtin entry + qr_dispatch (kernel router)
 *   qrdecomp.c           also hosts qr_symbolic_dispatch + qr_symbolic_core,
 *                        the rationalise -> exact MGS -> numericalise pipeline
 *                        used for symbolic, integer, rational, complex,
 *                        machine-real, and MPFR-real inputs in Phase 2.
 *   qrdecomp_machine.c  -- LAPACK kernel (added in Phase 3, real + complex)
 *   qrdecomp_mpfr.c     -- MPFR Householder kernel (added in Phase 4)
 *
 * Everything declared here is private to the QR sub-module.  Public
 * surface is `Expr* builtin_qrdecomposition(Expr*)` in qrdecomp.h.
 *
 * The architecture mirrors the eigen module split (see
 * src/linalg/eigen_internal.h): a thin top-level dispatcher classifies
 * the input by leaf precision and routes to the appropriate kernel;
 * each kernel owns its own load -> work -> store round-trip.
 */

#ifndef QRDECOMP_INTERNAL_H
#define QRDECOMP_INTERNAL_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "expr.h"

/* ---------------------------------------------------------------------
 * Parsed options.
 *
 *   Pivoting -> True   :  column-pivoted QR; result is {q, r, p} with
 *                         m . p == ConjugateTranspose[q] . r.
 *   Pivoting -> False  :  result is {q, r}, m == ConjugateTranspose[q] . r.
 *
 *   TargetStructure -> "Dense"       :  Phase 2 only supported form.
 *   TargetStructure -> "Structured"  :  reserved for OrthogonalMatrix /
 *                                       UnitaryMatrix wrappers; currently
 *                                       returns NULL so the call is left
 *                                       unevaluated.
 *
 * Produced by qr_parse_options; consumed by qr_dispatch and the
 * kernel-level dispatchers.
 * ------------------------------------------------------------------ */
typedef enum {
    QR_TS_DENSE = 0,
    QR_TS_STRUCTURED,
    QR_TS_INVALID
} QrTargetStructure;

typedef struct {
    bool              pivoting;
    QrTargetStructure target_structure;
} QrOpts;

/* Parse a Mathilda options sequence (positional matrix at args[0],
 * named options at args[1..n-1]).  Returns false and leaves *opts in
 * an indeterminate state if any option is malformed; the caller emits
 * QRDecomposition::opts and returns NULL. */
bool qr_parse_options(Expr* res, QrOpts* opts);

/* ---------------------------------------------------------------------
 * Top-level dispatcher.
 *
 * Classifies `m` by leaf precision via common_scan_inexact and routes
 * to the appropriate kernel.  Phase 2: every input goes to the
 * symbolic dispatcher (which itself handles inexact input through the
 * existing rationalise / numericalise round-trip).
 *
 * In Phase 3 this function will short-circuit MachinePrecision input
 * (min_bits <= 53) to qr_machine_dispatch, and in Phase 4 it will
 * short-circuit higher-precision MPFR input to qr_mpfr_dispatch.
 *
 * Returns a freshly-allocated Expr* on success, or NULL when the
 * kernel cannot evaluate (in which case the caller leaves the call
 * unevaluated; never frees `m`).
 * ------------------------------------------------------------------ */
Expr* qr_dispatch(Expr* m, int n, int p, const QrOpts* opts);

/* ---------------------------------------------------------------------
 * Symbolic kernel dispatcher.
 *
 * Wraps the Modified Gram-Schmidt core (qr_symbolic_core) with:
 *   - inexact-input rationalisation at min(input bits)
 *   - flatten -> work -> wrap conversion to/from List-of-List
 *   - tidy_matrix (element-wise Expand) on q and r
 *   - numericalisation back to input precision when rationalisation ran
 *   - assembly of {q, r} or {q, r, p}
 *
 * This is the same pipeline used in Phase 1; Phase 2 simply pulls it
 * out of builtin_qrdecomposition so that Phase 3 / 4 can introduce
 * sibling dispatchers without entangling them with option parsing.
 * ------------------------------------------------------------------ */
Expr* qr_symbolic_dispatch(Expr* m, int n, int p, const QrOpts* opts);

/* ---------------------------------------------------------------------
 * Symbolic kernel inner loop (Modified Gram-Schmidt).
 *
 * Pure C function over Expr** buffers; drives every primitive through
 * the Mathilda evaluator so it works on integer / rational / Complex /
 * Sqrt-bearing / free-symbolic entries with no special casing.
 *
 * Inputs:
 *   A_flat[i*p + k]  =  m[i, k],   i in [0, n),  k in [0, p)
 *   with_pivoting    chooses largest-residual-norm column per step
 *   perm             when non-NULL, populated with the column order
 *   use_conj         true for Hermitian inner product, false for
 *                    Euclidean (symbolic-real assumption)
 *
 * Outputs (all caller-owned, freshly allocated):
 *   *out_Q_flat   n * rank entries column-major   (Q_flat[i*rank+j])
 *   *out_R_flat   rank * p entries row-major      (R_flat[j*p+k])
 *   *out_rank     the numerical rank (<= min(n,p))
 *
 * Returns true on success.  A_flat is NOT consumed.
 * ------------------------------------------------------------------ */
bool qr_symbolic_core(Expr** A_flat, int n, int p,
                      bool with_pivoting, int* perm,
                      bool use_conj,
                      Expr*** out_Q_flat, Expr*** out_R_flat, int* out_rank);

/* ---------------------------------------------------------------------
 * Machine-precision kernel dispatcher (Phase 3, src/linalg/qrdecomp_machine.c).
 *
 * Invoked by qr_dispatch when common_scan_inexact reports an inexact
 * input at min_bits <= 53.  Loads `m` into a column-major double buffer,
 * calls LAPACK (dgeqrf / dgeqp3 / dorgqr, or the z* complex variants
 * via interleaved (re, im) pairs), reads back R and Q, drops rows
 * below numerical rank, and wraps the result in Mathilda lists.
 *
 * Returns NULL (without consuming `m`) in any of:
 *
 *   - USE_LAPACK is undefined.  The wrappers in lapack.c then report
 *     `info = -1` and the kernel falls back via qr_dispatch's symbolic
 *     branch.
 *   - The matrix contains a leaf the loader can't reduce to a double
 *     (e.g. a free-symbolic entry promoted into a machine-precision
 *     row).  The symbolic dispatcher handles that case correctly.
 *   - LAPACK reports a non-zero `info`.  Surfaces as a one-shot
 *     warning + symbolic fallback rather than a hard error.
 *
 * On success, returns {q, r} or {q, r, p} -- identical contract to
 * qr_symbolic_dispatch, with the same MachinePrecision Real / Complex
 * entries the symbolic round-trip would have produced.
 * ------------------------------------------------------------------ */
Expr* qr_machine_dispatch(Expr* m, int n, int p, const QrOpts* opts);

/* ---------------------------------------------------------------------
 * Arbitrary-precision MPFR kernel dispatcher (Phase 4,
 * src/linalg/qrdecomp_mpfr.c).
 *
 * Invoked by qr_dispatch when common_scan_inexact reports inexact
 * input above IEEE double (min_bits > 53).  Runs a Householder QR
 * over column-major MPFR arrays at the input's working precision;
 * uses paired re/im arrays for complex (no MPC dependency,
 * matches the eigen MPFR kernels).  Pivoting is implemented via
 * Businger-Golub with a rank-revealing cutoff at 2^(-bits/2) * |R[0,0]|.
 *
 * Returns NULL (without consuming `m`) in any of:
 *
 *   - USE_MPFR is undefined (the stub short-circuits to symbolic).
 *   - A matrix cell can't be reduced to an MPFR value -- e.g. a
 *     free-symbolic entry in a row that's otherwise MPFR-precision.
 *   - The input is rank-deficient AND pivoting is off; the symbolic
 *     dispatcher then handles mid-stream rank deficiency cleanly.
 *
 * On success, returns {q, r} or {q, r, p} with EXPR_MPFR (or
 * Complex[EXPR_MPFR, EXPR_MPFR]) entries at the input's precision --
 * identical public contract to qr_symbolic_dispatch.
 * ------------------------------------------------------------------ */
Expr* qr_mpfr_dispatch(Expr* m, int n, int p, const QrOpts* opts);

#endif /* QRDECOMP_INTERNAL_H */
