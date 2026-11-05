/*
 * flint_mat_bridge.h
 * ------------------
 * Boundary between Mathilda's Expr matrices and FLINT's exact dense
 * linear-algebra types (fmpq_mat), for the integer / rational case.
 *
 * The classical linalg kernels compute over generic Expr arithmetic so they
 * also work on symbolic entries; but for a matrix whose entries are all
 * integer or rational, FLINT's multimodular / fraction-free routines are
 * exact and dramatically faster (e.g. Det via O(n!) Laplace expansion in
 * det.c becomes polynomial time). Each entry point recognises the all-rational
 * case and returns NULL otherwise — leaving the caller on its classical path —
 * and is a no-op returning NULL when built without FLINT (USE_FLINT undefined),
 * matching the flint_bridge.c graceful-degrade policy.
 */
#ifndef FLINT_MAT_BRIDGE_H
#define FLINT_MAT_BRIDGE_H

#include "expr.h"
#include <gmp.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Determinant of an n×n matrix given row-major in `flat` (flat[r*n + c]), when
 * every entry is integer-like or Rational[·,·]. Computes exactly via FLINT's
 * fmpq_mat_det and renders the result as an Integer / Rational Expr the caller
 * owns. Returns NULL when an entry is not rational (symbolic matrix), n <= 0,
 * or FLINT is absent — the caller then uses the classical Laplace expansion.
 * Does not take ownership of `flat` or its elements.
 */
Expr* flint_mat_det(Expr** flat, int n);

/*
 * Inverse of an n×n matrix given row-major in `flat`, when every entry is
 * integer-like or Rational. Computes exactly via FLINT's fmpq_mat_inv and
 * renders the result as a List-of-Lists Expr the caller owns. Returns NULL
 * when an entry is not rational, the matrix is singular, n <= 0, or FLINT is
 * absent — the caller then uses the classical division-free path.
 * Does not take ownership of `flat` or its elements.
 */
Expr* flint_mat_inverse(Expr** flat, int n);

/*
 * Reduced row echelon form of an r×c matrix given row-major in `flat`, when
 * every entry is integer-like or Rational. Computes exactly via fmpq_mat_rref
 * and renders a List-of-Lists Expr the caller owns. Returns NULL when an entry
 * is not rational, r/c <= 0, or FLINT is absent.
 * Does not take ownership of `flat` or its elements.
 */
Expr* flint_mat_rref(Expr** flat, int r, int c);

/*
 * Solve the square linear system m·x = b, where m is r×c row-major in `flat_m`
 * and b is r×k row-major in `flat_b`, when every entry is integer-like or
 * Rational. Handles only the square unique-solution case (r == c, m
 * invertible) via fmpq_mat_solve. On success returns a freshly malloc'd
 * Expr** of length c*k (row-major, c rows × k cols), every element owned by
 * the caller; the caller must free the array with free(). Returns NULL for a
 * non-rational entry, a non-square/singular system, r/c/k <= 0, or no FLINT —
 * the caller then uses the classical division-free solver (which also handles
 * rectangular / underdetermined / inconsistent systems).
 * Does not take ownership of the input flats or their elements.
 */
Expr** flint_mat_solve(Expr** flat_m, int r, int c, Expr** flat_b, int k);

/*
 * Rank of an r×c matrix given row-major in `flat`, when every entry is
 * integer-like or Rational, via fmpq_mat_rref. Returns the rank (>= 0), or -1
 * when an entry is not rational, r/c <= 0, or FLINT is absent.
 * Does not take ownership of `flat` or its elements.
 */
int flint_mat_rank(Expr** flat, int r, int c);

/*
 * In-place reduced row echelon form of an r×c rational matrix supplied as a
 * flat row-major mpq_t array `M` (M[i*c + j]; every entry already mpq_init'd).
 * Computed exactly via FLINT's fraction-free fmpq_mat_rref — far faster than a
 * hand-rolled mpq_t Gauss-Jordan on systems with large rational entries, whose
 * numerators blow up under repeated pivoting. On return M holds the (unique)
 * RREF in natural row order, and pivot_for_col[j] (length c, caller-allocated)
 * is the physical row that pivots column j, or -1 if column j has no pivot.
 * Returns 1 on success, 0 when FLINT is absent (caller keeps its mpq path).
 * Does not allocate or free M.
 */
int flint_rref_mpq_inplace(mpq_t* M, int r, int c, int64_t* pivot_for_col);

/*
 * QRDecomposition[m] of an n x p matrix of integers / rationals (row-major
 * `flat`), as the reference answers it: Gram-Schmidt over Q on the columns,
 * dependent columns dropped, so q is rank x n with ConjugateTranspose[q] . r
 * == m.  The orthogonal (unnormalised) vectors v_j and their squared norms
 * N_j stay rational; only the normalisation introduces a square root, so
 * q[j, i] = v_j[i] / Sqrt[N_j] and r[j, k] = <v_j, a_k> / Sqrt[N_j], each
 * built once and evaluated to the canonical radical (1/Sqrt[51],
 * 2 Sqrt[3/17], ...).  The symbolic Gram-Schmidt this replaces simplified
 * nested radicals at every step: a 5 x 5 integer matrix took a minute.
 * Returns {q, r}, or NULL when an entry is not rational (or no FLINT).
 */
Expr* flint_mat_qr(Expr** flat, int n, int p);

/* Registers the FLINT` context matrix builtins (FLINT`Det, FLINT`Inverse,
 * FLINT`LinearSolve, FLINT`RowReduce, FLINT`MatrixRank). Called from
 * core_init(). No-op without FLINT. */
void flint_mat_bridge_init(void);

#ifdef __cplusplus
}
#endif

#endif /* FLINT_MAT_BRIDGE_H */
