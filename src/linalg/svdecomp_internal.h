/*
 * src/linalg/svdecomp_internal.h
 *
 * Internal interface for the SingularValueDecomposition translation
 * units.
 *
 *   svdecomp.c          -- builtin entry, option/positional parser,
 *                          top-level svd_dispatch (kernel router), and
 *                          the symbolic dispatcher + symbolic core
 *                          (eigendecomposition of m^H m, with a 2x2
 *                          closed-form fast path).
 *   svdecomp_machine.c  -- LAPACK kernel (dgesdd / zgesdd for the
 *                          standard form, dggsvd3 / zggsvd3 for the
 *                          generalized {m, a} form).
 *   svdecomp_mpfr.c     -- one-sided Jacobi SVD over MPFR arrays
 *                          (real + complex via paired re/im).
 *
 * Everything declared here is private to the SVD sub-module.  Public
 * surface is `Expr* builtin_singularvaluedecomposition(Expr*)` in
 * svdecomp.h.
 *
 * The architecture mirrors qrdecomp_internal.h / ludecomp_internal.h:
 * a thin top-level dispatcher classifies the input by leaf precision
 * and routes to the appropriate kernel; each kernel owns its own
 * load -> work -> store round-trip and falls back to the symbolic
 * dispatcher (return NULL, never frees the input matrix) on any soft
 * failure.
 */

#ifndef SVDECOMP_INTERNAL_H
#define SVDECOMP_INTERNAL_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "expr.h"

/* ---------------------------------------------------------------------
 * Parsed positional + option form.
 *
 * Produced by svd_parse_args; consumed by every dispatcher.  None of
 * the Expr* fields are owned: they alias the live args[] of the
 * builtin's `res` (which the evaluator owns); they must not be freed
 * by the dispatchers.
 * ------------------------------------------------------------------ */

typedef enum {
    SVD_FORM_FULL = 0,   /* no k argument -- full thin SVD */
    SVD_FORM_K,          /* explicit integer k (signed) */
    SVD_FORM_UPTO        /* UpTo[k] -- k largest, or as many as available */
} SvdKForm;

typedef enum {
    SVD_TS_DENSE = 0,
    SVD_TS_STRUCTURED,
    SVD_TS_INVALID
} SvdTargetStructure;

typedef struct {
    Expr*               m;             /* matrix, or first of the pair */
    Expr*               a;             /* second of the pair when generalized; NULL otherwise */
    bool                generalized;   /* true iff arg[0] was a {m, a} pair */

    SvdKForm            k_form;
    int                 k_value;       /* signed; |k| <= min(n, p); 0 for SVD_FORM_FULL */

    Expr*               tolerance;     /* NULL when not supplied */
    SvdTargetStructure  target_structure;
} SvdArgs;

/* Parse the positional and option arguments of a SingularValueDecomposition
 * call.  Returns false (and leaves *args in an indeterminate state) if
 * any positional or option entry is malformed; the caller emits
 * SingularValueDecomposition::opts or ::sval and returns NULL. */
bool svd_parse_args(Expr* res, SvdArgs* args);

/* ---------------------------------------------------------------------
 * Top-level dispatcher.
 *
 * Classifies `args->m` (and `args->a` for the generalized form) by leaf
 * precision via common_scan_inexact and routes to the appropriate
 * kernel.  Returns a freshly-allocated Expr* on success, or NULL when
 * every kernel declined (in which case the caller leaves the call
 * unevaluated; never frees the matrix arguments).
 * ------------------------------------------------------------------ */
Expr* svd_dispatch(const SvdArgs* args, int n, int p, int n_a);

/* ---------------------------------------------------------------------
 * Symbolic kernel dispatcher.
 *
 * Wraps svd_symbolic_core with:
 *   - inexact-input rationalisation at min(input bits)
 *   - flatten -> work -> wrap conversion to/from List-of-List
 *   - tidy_matrix (element-wise Together / Expand) on u, sigma, v
 *   - numericalisation back to input precision when rationalisation
 *     was applied
 *   - assembly of {u, sigma, v} (or {{u, ua}, {sigma, sigma_a}, v}
 *     for the generalized form)
 *   - shared truncation / tolerance / TargetStructure post-processing
 *     via svd_apply_postprocess
 * ------------------------------------------------------------------ */
Expr* svd_symbolic_dispatch(const SvdArgs* args, int n, int p, int n_a);

/* ---------------------------------------------------------------------
 * Symbolic kernel inner loop.
 *
 * Pure C function over Expr** buffers; drives every primitive through
 * the Mathilda evaluator so it works on integer / rational / Complex /
 * Sqrt-bearing / free-symbolic entries with no special casing.
 *
 * Inputs:
 *   A_flat[i*p + k]  =  m[i, k],   i in [0, n),  k in [0, p)
 *   use_conj         true for Hermitian inner product, false for
 *                    Euclidean (symbolic-real assumption)
 *
 * Outputs (all caller-owned, freshly allocated):
 *   *out_U_flat   n * rank entries row-major   (U_flat[i*rank+j])
 *   *out_S_flat   rank entries (singular values, descending where the
 *                 evaluator can order them)
 *   *out_V_flat   p * rank entries row-major   (V_flat[k*rank+j])
 *   *out_rank     the numerical rank (<= min(n, p))
 *
 * Returns true on success.  A_flat is NOT consumed.  Returns false
 * when the symbolic eigendecomposition could not produce a closed
 * form -- the dispatcher then falls back to N[m, 53] before retrying.
 * ------------------------------------------------------------------ */
bool svd_symbolic_core(Expr** A_flat, int n, int p, bool use_conj,
                       Expr*** out_U_flat,
                       Expr*** out_S_flat,
                       Expr*** out_V_flat,
                       int* out_rank);

/* ---------------------------------------------------------------------
 * Machine-precision (LAPACK) kernel dispatcher.
 *
 * Invoked by svd_dispatch when common_scan_inexact reports an inexact
 * input at min_bits <= 53.  Loads `m` (and `a` for the generalized
 * form) into column-major double buffers, calls dgesdd / zgesdd (or
 * dggsvd3 / zggsvd3 for {m, a}), and wraps the result in Mathilda
 * lists.
 *
 * Returns NULL (without consuming the matrices) in any of:
 *   - USE_LAPACK is undefined.
 *   - The matrix contains a leaf the loader can't reduce to a double.
 *   - LAPACK reports a non-zero info.
 * ------------------------------------------------------------------ */
Expr* svd_machine_dispatch(const SvdArgs* args, int n, int p, int n_a);

/* ---------------------------------------------------------------------
 * Arbitrary-precision MPFR kernel dispatcher.
 *
 * Invoked by svd_dispatch when common_scan_inexact reports inexact
 * input above IEEE double (min_bits > 53).  Runs a one-sided Jacobi
 * SVD over column-major MPFR arrays at the input's working precision.
 * Uses paired re/im arrays for complex (no MPC dependency).
 *
 * Returns NULL (without consuming the matrices) in any of:
 *   - USE_MPFR is undefined.
 *   - A matrix cell can't be reduced to an MPFR value.
 * ------------------------------------------------------------------ */
Expr* svd_mpfr_dispatch(const SvdArgs* args, int n, int p, int n_a);

/* ---------------------------------------------------------------------
 * Shared post-processing applied to a freshly assembled {u, sigma, v}
 * (or the generalized {{u, ua}, {sigma, sigma_a}, v}).
 *
 *   - Applies the Tolerance: any singular value with absolute value
 *     below the threshold is zeroed.
 *   - Applies the k / UpTo[k] truncation: keeps the first k columns
 *     (k > 0) or the last |k| columns (k < 0) of u, sigma, v.
 *   - Wraps sigma as DiagonalMatrix[{...}] when TargetStructure
 *     is "Structured".
 *
 * Consumes `result` and returns a fresh Expr* (the same expression
 * when no post-processing is needed).
 *
 * `rank` is the rank inferred by the kernel; `n` and `p` are the
 * original matrix dimensions.  For the generalized form,
 * svd_apply_postprocess is a no-op: tolerance and structure are still
 * honoured but truncation is rejected at parse time with ::nopart.
 * ------------------------------------------------------------------ */
Expr* svd_apply_postprocess(Expr* result, const SvdArgs* args,
                            int n, int p, int rank);

#endif /* SVDECOMP_INTERNAL_H */
