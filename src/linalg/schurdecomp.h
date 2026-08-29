/*
 * schurdecomp.h -- SchurDecomposition[m] and SchurDecomposition[{m, a}].
 *
 * SchurDecomposition[m] gives the Schur decomposition of a numerical square
 * matrix m as {q, t}, where q is orthonormal (unitary), t is block upper-
 * triangular (the Schur form), and m == q . t . ConjugateTranspose[q].  For a
 * real matrix the default RealBlockDiagonalForm -> True keeps t real with 2x2
 * blocks on the diagonal for complex-conjugate eigenvalue pairs; -> False makes
 * t complex upper-triangular.  Pivoting -> True additionally returns a scaling /
 * permutation matrix d with m . d == d . q . t . ConjugateTranspose[q].
 *
 * SchurDecomposition[{m, a}] gives the generalized (QZ) Schur decomposition as
 * {q, s, p, t} with q, p orthonormal, s, t upper-triangular, m == q.s.p^H and
 * a == q.t.p^H.
 *
 * Numerical only -- there is no closed-form symbolic Schur form of a generic
 * matrix, so a non-numeric matrix leaves the call unevaluated.  Machine input
 * uses LAPACK (dgees / zgees / dgges / zgges, with dgebal for Pivoting); an
 * arbitrary-precision real standard matrix uses the in-house MPFR Hessenberg +
 * Francis QR (eigen_schur_real_mpfr).
 *
 * Memory contract: standard builtin ownership (SPEC.md §4).  Never frees `res`;
 * returns a fresh Expr* the evaluator owns, or NULL to leave the call
 * unevaluated.
 */
#ifndef MATHILDA_LINALG_SCHURDECOMP_H
#define MATHILDA_LINALG_SCHURDECOMP_H

#include "expr.h"

#ifdef __cplusplus
extern "C" {
#endif

Expr* builtin_schurdecomposition(Expr* res);
void  schurdecomp_init(void);

#ifdef __cplusplus
}
#endif

#endif /* MATHILDA_LINALG_SCHURDECOMP_H */
