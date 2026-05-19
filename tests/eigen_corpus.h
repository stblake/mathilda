/* Shared fixtures and tolerance helpers for the numerical Eigenvalues /
 * Eigenvectors tests (Phase 2 onwards).
 *
 * The corpus is intentionally small and predictable — these helpers
 * deliberately do NOT pull randomness from the global state of any
 * other test binary.  Each test that needs a "random" matrix calls
 * eigen_corpus_*_random with a fixed seed so failures are reproducible.
 *
 * Tolerance discipline (matches Golub & Van Loan Ch 7):
 *   - eigenvalue residual: ||A v - lambda v||_inf  <=  c * n * eps * ||A||_inf
 *   - orthonormality:      ||V^T V - I||_inf       <=  c * n * eps
 *   - eigenvalue sum:      |sum lambda_i - trace A|  <=  c * n * eps * ||A||_inf
 *
 * `c` is a small constant (we use 32) to absorb the inevitable
 * coefficient from a stable Householder + Wilkinson-shift QR
 * implementation; tighter constants are achievable but pointlessly
 * fragile.
 */
#ifndef MATHILDA_TESTS_EIGEN_CORPUS_H
#define MATHILDA_TESTS_EIGEN_CORPUS_H

#include <stddef.h>
#include "expr.h"

/* Evaluate `input`, expect a List of n machine-precision real
 * eigenvalues that agree with the externally-supplied `expected`
 * values (in the same descending |lambda| order) to within `tol`.
 * Prints PASS / FAIL and asserts on mismatch. */
void corpus_check_real_eigenvalues(const char* label,
                                    const char* input,
                                    const double* expected, size_t n,
                                    double tol);

/* Build the Mathilda input string `Eigenvalues[<flat>]` (or
 * `Eigenvectors[...]`) from a flat row-major real n x n matrix.
 * The caller must free the returned heap buffer. */
char* corpus_matrix_to_eigenvalues_input(const double* A, size_t n);
char* corpus_matrix_to_eigenvectors_input(const double* A, size_t n);

/* Evaluate `Eigenvalues[A]` and read the result into the freshly
 * allocated `*lambdas` array of length n (caller frees).  Returns the
 * number of eigenvalues actually returned (may be < n if the call
 * went unevaluated; the test should ASSERT on the count). */
size_t corpus_eval_eigenvalues_real(const double* A, size_t n,
                                     double** lambdas);

/* Evaluate `Eigenvectors[A]` and read the result into freshly
 * allocated `*V` of length n*n (row-major: row i is the i-th
 * eigenvector).  Returns the number of eigenvectors returned. */
size_t corpus_eval_eigenvectors_real(const double* A, size_t n,
                                      double** V);

/* Residual check: ||A v - lambda v||_inf  vs  tol * ||A||_inf.
 * Asserts on failure; returns the residual on success. */
double corpus_assert_residual_real(const double* A, size_t n,
                                    double lambda, const double* v,
                                    double tol);

/* ||V V^T - I||_inf check.  V is row-major n x n with row i = i-th
 * eigenvector (the layout corpus_eval_eigenvectors_real returns).
 * Asserts on failure. */
void corpus_assert_orthonormal_real(const double* V, size_t n, double tol);

/* Infinity-norm of an n x n real matrix. */
double corpus_norm_inf_real(const double* A, size_t n);

#endif /* MATHILDA_TESTS_EIGEN_CORPUS_H */
