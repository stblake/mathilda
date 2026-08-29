/*
 * schurdecomp_internal.h -- shared declarations for the SchurDecomposition
 * kernels (schurdecomp.c dispatcher, schurdecomp_machine.c LAPACK kernels,
 * schurdecomp_mpfr.c arbitrary-precision kernel).
 */
#ifndef MATHILDA_LINALG_SCHURDECOMP_INTERNAL_H
#define MATHILDA_LINALG_SCHURDECOMP_INTERNAL_H

#include "expr.h"
#include <stdbool.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* Parsed option bundle. */
typedef struct {
    bool pivoting;                 /* Pivoting -> True (standard form only):
                                    * also return the scaling/permutation d.   */
    bool real_block_diagonal_form; /* RealBlockDiagonalForm (default True):
                                    * True keeps a real quasi-triangular t;
                                    * False makes t complex upper-triangular.  */
    bool structured;               /* TargetStructure -> "Structured": accepted,
                                    * but returns dense matrices (Mathilda has
                                    * no structured-matrix type -- t is already
                                    * triangular and q orthonormal), so the flag
                                    * carries no effect.                        */
} SchurOpts;

/* Parse the option arguments args[1..] of `res` into `opts`.  Returns false
 * (so the caller leaves the whole call unevaluated) on an unknown option head
 * or a malformed / unrecognised value. */
bool schur_parse_options(const Expr* res, SchurOpts* opts);

/* Square order of a matrix argument -- a rank-2 NDArray or a boxed
 * List-of-Lists; -1 for anything that is not a non-empty square matrix. */
int schur_matrix_order(const Expr* e);

/* Machine (LAPACK) kernels.  schur_machine_standard returns {q, t} (or
 * {q, t, d} under Pivoting); schur_machine_generalized returns {q, s, p, t}.
 * Each returns NULL when LAPACK is unavailable, the matrix is not numeric, or
 * the factorization did not converge -- the call is then left unevaluated. */
Expr* schur_machine_standard(const Expr* m, int n, const SchurOpts* opts);
Expr* schur_machine_generalized(const Expr* m, const Expr* a, int n,
                                const SchurOpts* opts);

#ifdef USE_MPFR
/* Arbitrary-precision standard real Schur decomposition (RealBlockDiagonalForm
 * -> True, no Pivoting), computed with the in-house MPFR Hessenberg + Francis
 * QR.  Returns {q, t} at `bits` precision, or NULL to fall back to the machine
 * kernel. */
Expr* schur_mpfr_standard_real(const Expr* m, int n, mpfr_prec_t bits,
                               const SchurOpts* opts);
#endif

#endif /* MATHILDA_LINALG_SCHURDECOMP_INTERNAL_H */
