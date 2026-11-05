/* Mathilda -- sparse linear algebra: shared internals of src/sparse/sparse_la_*.c
 *
 * The native SparseArray linear-algebra kernels (LinearSolve, Det, Inverse,
 * LinearSolveFunction, Eigenvalues / Eigensystem / SingularValueList with k,
 * MatrixExp[s, v], ...) are built from three layers:
 *
 *   SpMat      a borrowed machine view of a rank-2 SparseArray: CSR with
 *              int64 column indices and a double buffer that is real
 *              (x[k]) or complex interleaved (x[2k], x[2k+1]). Built from an
 *              SPV_F64 / SPV_C128 array (and an SPV_INT64 one promoted to
 *              double) with a zero background.
 *   SpFactor   a direct factorisation of a square SpMat: UMFPACK (LU),
 *              CHOLMOD (Cholesky) when SuiteSparse is linked, else a dense
 *              LAPACK LU of the densified matrix. Solves A x = b and
 *              A^T x = b for many right-hand sides; gives Det.
 *   Krylov     CG / BiCGSTAB / GMRES(m) with None / Jacobi / ILU0 / ILUT /
 *              ILUTP / user-function preconditioning, on SpMV.
 *
 * Index convention: everything here is 0-based; the CSR arrays of an SpMat
 * describe A row by row, so as a column-compressed matrix (what UMFPACK and
 * CHOLMOD read) they describe A^T. The factor code accounts for that.
 */
#ifndef MATHILDA_SPARSE_LA_H
#define MATHILDA_SPARSE_LA_H

#include "sparse.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* ---- machine view -------------------------------------------------------- */

typedef struct {
    int64_t m, n, nnz;
    const int64_t* rp;   /* m + 1 row pointers (borrowed)                  */
    const int64_t* ci;   /* nnz column indices (borrowed or owned_ci)       */
    const double* x;     /* nnz reals, or 2*nnz interleaved complex         */
    bool cplx;
    int64_t* owned_ci;   /* non-NULL when ci was widened from int32         */
    double* owned_x;     /* non-NULL when x was converted (int64 -> double)  */
} SpMat;

/* Fill `A` with a machine view of a rank-2 SparseArray whose background is
 * zero (exact 0, 0., or 0.+0.I) and whose values are SPV_F64 / SPV_C128, or
 * SPV_INT64 when `promote_int` (converted to double). False otherwise; then
 * nothing needs freeing. */
bool spla_view(const SparseArray* sa, bool promote_int, SpMat* A);
void spla_view_free(SpMat* A);

/* Build an SpMat (all storage owned) from a dense machine matrix held as a
 * List / packed List / NDArray of reals or complexes. False when `e` is not a
 * machine matrix. Explicit zeros are dropped. */
bool spla_from_dense(const Expr* e, SpMat* A);

/* y = A x  (or A^T x when trans; for complex A the plain, not conjugate,
 * transpose). x and y are length n / m (or m / n) and may not alias. */
void spla_spmv(const SpMat* A, const double* x, double* y, bool trans);

/* Structural queries. */
bool spla_is_square(const SpMat* A);
/* Symmetric (real) or Hermitian (complex) to the exact bit, and every
 * diagonal entry real and positive: the Cholesky-worth-trying test. */
bool spla_hermitian_posdiag(const SpMat* A);

/* ---- direct factorisation ------------------------------------------------ */

typedef enum {
    SPF_OK = 0,
    SPF_SINGULAR,      /* exactly / numerically singular                    */
    SPF_NOTPD,         /* Cholesky requested and the matrix is not HPD      */
    SPF_NOTSQUARE,
    SPF_TOOBIG,        /* dense fallback would exceed the densify cap       */
    SPF_NOMEM,
    SPF_FAILED         /* library error                                     */
} SpfStatus;

typedef enum {
    SPF_WANT_AUTO = 0, /* Cholesky when Hermitian with a positive diagonal and
                        * it succeeds, else LU                               */
    SPF_WANT_LU,       /* "Multifrontal"                                     */
    SPF_WANT_CHOL,     /* "Cholesky"                                         */
    SPF_WANT_DENSE     /* always the dense LAPACK LU                         */
} SpfWant;

typedef struct SpFactor SpFactor;

/* Factor the square machine matrix A. NULL with *st set on failure. The
 * factor does not borrow A. */
SpFactor* spf_factor(const SpMat* A, SpfWant want, SpfStatus* st);
/* Solve A X = B (trans false) or A^T X = B (trans true; plain transpose for
 * complex) for `nrhs` right-hand sides stored one after another, each of
 * length n (complex: 2n doubles, interleaved). X may alias B. */
bool spf_solve(const SpFactor* f, const double* B, double* X, int64_t nrhs, bool trans);
/* det(A) as re + i im (im 0 for real). False if the backend cannot say. */
bool spf_det(const SpFactor* f, double* re, double* im);
/* 1-norm reciprocal condition estimate, or -1 when unavailable. */
double spf_rcond(const SpFactor* f);
int64_t spf_n(const SpFactor* f);
bool spf_cplx(const SpFactor* f);
const char* spf_backend(const SpFactor* f);   /* "UMFPACK", "CHOLMOD", "LAPACK" */
void spf_free(SpFactor* f);

/* Which optional libraries this build links: for $Version-style reporting and
 * for tests that pin a backend. */
bool spla_have_umfpack(void);
bool spla_have_cholmod(void);
bool spla_have_arpack(void);

/* ---- Krylov -------------------------------------------------------------- */

typedef enum { KRY_AUTO = 0, KRY_CG, KRY_BICGSTAB, KRY_GMRES } KryMethod;
typedef enum { KPC_AUTO = 0, KPC_NONE, KPC_ILU0, KPC_ILUT, KPC_ILUTP, KPC_JACOBI, KPC_USER } KryPrec;

/* A user preconditioner: z = M^-1 r, both length n (complex interleaved).
 * Return false to abort the iteration. */
typedef bool (*KryUserPrec)(void* ctx, const double* r, double* z);

typedef struct {
    KryMethod method;
    KryPrec   prec;
    double    tol;            /* relative residual ||b - A x|| / ||b||      */
    int64_t   max_iter;       /* <= 0: the default, which scales with n      */
    int       restart;        /* GMRES(m) restart length (<= 0: 100)         */
    double    ilut_droptol;   /* ILUT / ILUTP drop tolerance                 */
    int       ilut_fill;      /* ILUT / ILUTP extra fill per row             */
    KryUserPrec user;
    void*     user_ctx;
} KryOpts;

typedef enum { KRY_CONVERGED = 0, KRY_MAXIT, KRY_BREAKDOWN, KRY_PRECFAIL, KRY_NOMEM } KryStatus;

void kry_defaults(KryOpts* o);
/* Solve A x = b; x holds the starting vector on entry (zeros for none) and
 * the last iterate on return. *iters receives the iteration count. */
KryStatus kry_solve(const SpMat* A, const double* b, double* x, const KryOpts* o, int64_t* iters);

/* ---- entry points used outside src/sparse --------------------------------- */

/* Registration of the native heads (called from sparse_init after the gate). */
void sparse_la_init(void);
void sparse_la_func_init(void);   /* MatrixExp (sparse_la_func.c) */

/* Solve A x = B for a machine SparseArray A with the LinearSolve options in
 * `opts` (Method -> ...). Used by the LinearSolve native handler and by
 * LinearSolve[m, b, Method -> "Krylov" | ...] on a dense matrix. Returns the
 * answer, NULL to decline, or expr_copy(call) to leave the call unevaluated
 * after a message. */
Expr* spla_linearsolve_call(const Expr* call, const Expr* A, const Expr* B,
                            const Expr* const* opts, size_t nopts);

#endif /* MATHILDA_SPARSE_LA_H */
