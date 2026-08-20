/*
 * hnf.h -- Hermite Normal Form over Z.
 *
 * Row-style HNF with the unimodular transform recorded, plus the
 * `HermiteDecomposition[A]` builtin.  The row HNF is the reusable integer
 * linear-algebra primitive behind exact linear Diophantine system solving
 * (src/solve/solveint_linear.c) and lattice / module computations.
 */
#ifndef MATHILDA_HNF_H
#define MATHILDA_HNF_H

#include <gmp.h>
#include "expr.h"

/* Row Hermite Normal Form over Z.
 *
 * Given an m x n integer matrix A (flat, row-major, `A[i*n + j]`, borrowed),
 * allocate and compute a unimodular matrix P (m x m, |det| = 1) and R (m x n)
 * with  P * A == R  in row-style HNF: R is upper-triangular in echelon shape
 * (pivots move strictly right as rows descend), every pivot is positive, and
 * each entry ABOVE a pivot is reduced into [0, pivot).
 *
 * `*R_out` (m*n entries) and `*P_out` (m*m entries) are freshly malloc'd flat
 * mpz_t arrays, each entry mpz_init'd; the caller releases them with
 * `linalg_hnf_free(ptr, count)`.  Returns the rank r in [0, min(m, n)].
 * Returns -1 on allocation failure (outputs left NULL).
 */
int  linalg_hnf(const mpz_t* A, int m, int n, mpz_t** R_out, mpz_t** P_out);

/* Clear every mpz_t in `mat` (a flat array of `count` entries) and free it. */
void linalg_hnf_free(mpz_t* mat, int count);

/* HermiteDecomposition[A] -> {U, R} with U unimodular, R the row HNF, and
 * U . A == R.  Integer matrix only; declines (NULL) otherwise. */
Expr* builtin_hermite_decomposition(Expr* res);

#endif /* MATHILDA_HNF_H */
