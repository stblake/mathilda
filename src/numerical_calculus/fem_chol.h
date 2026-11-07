/* Mathilda -- sparse Cholesky factorisation for finite element matrices.
 *
 * A symmetric positive definite matrix A (full CSR, both triangles) whose
 * unknowns sit at points of the plane is factored as P A P^T = L L^T:
 *
 *   ordering  geometric nested dissection -- the points are split at the
 *             median of the longer side of their bounding box, the nodes of
 *             one half adjacent to the other half form the separator, which
 *             is numbered last, and both halves are dissected in turn (leaves
 *             of <= 64 nodes keep their order). For a planar mesh this gives
 *             nnz(L) = O(n log n) and O(n^1.5) work (George 1973), against the
 *             O(n^2) band of a Cuthill-McKee ordering.
 *   symbolic  the elimination tree and the row patterns of L by the
 *             "ereach" traversal (Davis, Direct Methods for Sparse Linear
 *             Systems, 2006, ch. 4), which also gives the column counts.
 *   numeric   up-looking Cholesky, one row of L at a time.
 *
 * Used by the FEM solver when SuiteSparse (CHOLMOD) is not linked, where the
 * band factor of sparse_la_band.c needs O(n^2) storage and stops at ~10^5
 * unknowns. */
#ifndef MATHILDA_FEM_CHOL_H
#define MATHILDA_FEM_CHOL_H

#include <stdbool.h>
#include <stdint.h>

typedef struct FemChol FemChol;

/* Factor the n x n matrix (CSR rp, ci, x; symmetric) whose unknown k is at
 * (xy[2k], xy[2k+1]). NULL when it is not positive definite (or too large),
 * with *why set. */
FemChol* fem_chol_factor(int64_t n, const int64_t* rp, const int64_t* ci, const double* x,
                         const double* xy, const char** why);
/* x = A^-1 b (x may alias b). */
void fem_chol_solve(const FemChol* F, const double* b, double* x);
void fem_chol_free(FemChol* F);
int64_t fem_chol_nnz(const FemChol* F);

#endif /* MATHILDA_FEM_CHOL_H */
