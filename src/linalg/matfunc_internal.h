/*
 * matfunc_internal.h -- shared declarations for the matrix-function family
 * (MatrixExp, MatrixLog, MatrixFunction, and MatrixPower's non-integer
 * exponents):
 *
 *   matfunc.c          dispatchers, argument checks, messages, ND entry points
 *   matfunc_machine.c  machine-precision kernels (double / interleaved complex)
 *   matfunc_mpfr.c     arbitrary-precision kernels (complex MPFR matrices)
 *   matfunc_exact.c    exact / symbolic engine over the Jordan form
 *
 * Machine buffers are COLUMN-MAJOR (the LAPACK layout); a complex buffer
 * interleaves (re, im) per element, so an n x n complex matrix occupies 2 n^2
 * doubles.  Every kernel BORROWS its inputs and writes into caller-owned
 * output storage; nothing here allocates an Expr.
 */
#ifndef MATHILDA_LINALG_MATFUNC_INTERNAL_H
#define MATHILDA_LINALG_MATFUNC_INTERNAL_H

#include "expr.h"
#include <stdbool.h>
#include <stddef.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* Kernel status codes. */
enum {
    MF_OK       = 0,
    MF_FAIL     = 1,   /* generic failure: non-finite data, no LAPACK, ...   */
    MF_ZEROEIG  = 2,   /* a zero eigenvalue where the function is singular   */
    MF_NOSOL    = 3,   /* no primary function value (defective at a branch
                        * point, e.g. the square root of a nilpotent block)  */
    MF_NONNUM   = 4    /* the scalar function did not evaluate to a number   */
};

/* The scalar function behind a Schur-Parlett evaluation: f^(k)(z), the k-th
 * derivative at the complex point z = (zr, zi), written to (*fr, *fi).
 * Returns MF_OK, or a status code (MF_NONNUM when the value is not a number,
 * MF_NOSOL / MF_ZEROEIG at a singularity) that aborts the evaluation. */
typedef int (*mf_scalar_fn)(void* ctx, double zr, double zi, int k,
                            double* fr, double* fi);

/* The closed-form scalar functions the machine kernels share with the 1x1
 * path: exp, the principal log and the principal power z^s.  ctx for the
 * power is a double[2] holding s = (re, im). */
int mf_fn_exp(void* ctx, double zr, double zi, int k, double* fr, double* fi);
int mf_fn_log(void* ctx, double zr, double zi, int k, double* fr, double* fi);
int mf_fn_pow(void* ctx, double zr, double zi, int k, double* fr, double* fi);

/* ---- machine kernels (matfunc_machine.c) ------------------------------- */

/* exp(A) by Pade scaling and squaring with the Al-Mohy--Higham (2009)
 * backward-error bounds.  A is n x n (real, or complex when cx); E receives
 * the result in the same element type.  Returns MF_OK or MF_FAIL. */
int mf_expm(const double* A, int n, bool cx, double* E);

/* exp(A) b without forming exp(A): the Al-Mohy--Higham (2011) truncated
 * Taylor action.  The work is complex when either operand is; `out` has the
 * element type cx || bcx and length n. */
/* The truncated-Taylor degrees m and their theta_m (unit roundoff 2^-53),
 * shared with the sparse action exp(A) b (src/sparse/sparse_la_func.c). */
extern const int    EXPMV_M[];
extern const double EXPMV_THETA[];
extern const size_t EXPMV_NTHETA;

int mf_expmv(const double* A, int n, bool cx, const double* b, bool bcx,
             double* out);

/* The principal f(A) for a scalar function f, through the Schur form:
 *   kind MF_KIND_LOG  -- inverse scaling and squaring (Al-Mohy--Higham 2012)
 *   kind MF_KIND_SQRT -- the Bjorck--Hammarling triangular square root
 *   kind MF_KIND_FUNC -- Schur--Parlett (Davies--Higham 2003) with `fn`
 * The output F is always complex (2 n^2 doubles); *is_real reports that a real
 * A gave a real answer (every eigenvalue's value is consistent with a real
 * matrix function), in which case the imaginary parts are exactly zero and
 * the caller may present the result as a real matrix. */
enum { MF_KIND_FUNC = 0, MF_KIND_LOG = 1, MF_KIND_SQRT = 2 };
int mf_schur_fun(const double* A, int n, bool cx, int kind,
                 mf_scalar_fn fn, void* ctx, double* F, bool* is_real);

/* ---- MPFR kernels (matfunc_mpfr.c) -------------------------------------- */
#ifdef USE_MPFR
/* An n x n complex MPFR matrix, ROW-major, (re, im) in two arrays. */
typedef struct {
    int n;
    mpfr_prec_t prec;
    mpfr_t* re;
    mpfr_t* im;
} MfMat;

/* Load an n x n numeric matrix at precision `prec`; false on a non-numeric
 * leaf.  mf_mat_free releases one; mf_mat_to_expr builds a List of MPFR (or
 * Complex[mpfr, mpfr]) entries rounded to `out_bits`; `force_real` drops the
 * imaginary parts. */
bool  mf_mat_load(const Expr* m, int n, mpfr_prec_t prec, MfMat* out);
void  mf_mat_free(MfMat* a);
Expr* mf_mat_to_expr(const MfMat* a, mpfr_prec_t out_bits, bool force_real);

/* exp(A) in place, by scaling and squaring with a Taylor degree chosen for
 * A's precision.  log(A) and A^s (s complex, at A's precision) in place by
 * inverse scaling and squaring over Denman--Beavers square roots.  Each
 * returns an MF_* status. */
int mf_mpfr_expm(MfMat* a);
int mf_mpfr_logm(MfMat* a);
int mf_mpfr_powm(MfMat* a, const mpfr_t s_re, const mpfr_t s_im);
int mf_mpfr_sqrtm(MfMat* a);
/* f(A) for a scalar function f (any expression f[z] evaluates at a complex
 * MPFR z) by the Parlett recurrence on the complex Schur form, at the input
 * precision `bits`.  NULL (the caller keeps its exact path) on a non-numeric
 * f value or eigenvalues closer than 2^-bits of their scale. */
Expr* mf_mpfr_function(const Expr* m, int n, long bits, const Expr* f);
#endif

/* ---- exact / symbolic engine (matfunc_exact.c) -------------------------- */

/* The function kinds the exact engine and the messages distinguish. */
typedef enum {
    MFX_EXP,       /* MatrixExp                  */
    MFX_LOG,       /* MatrixLog                  */
    MFX_POWER,     /* MatrixPower[m, s]          */
    MFX_FUNCTION   /* MatrixFunction[f, m]       */
} MfxKind;

/* f(m) over the Jordan form of an exact or symbolic square matrix m:
 * f(m) = sum over eigenvalues lambda and k of f^(k)(lambda)/k! C_{lambda,k},
 * with the coefficient matrices C taken from m = S J S^-1.  `f` is a function
 * value (symbol or pure function); `head` names the calling builtin for its
 * messages.  Returns a fresh List-of-Lists, or NULL (after any message) to
 * leave the call unevaluated. */
Expr* mf_exact_apply(const Expr* m, int n, const Expr* f, MfxKind kind,
                     const char* head);

/* Function[x, D[f[x], {x, k}]] (f itself for k = 0), a fresh Expr. */
Expr* mf_derivative_function(const Expr* f, int k);

/* The eigen-decomposition route for a machine matrix whose scalar function does not
 * evaluate to numbers (MatrixFunction[g, m] with an undefined g, or
 * MatrixPower[m, x] with a symbolic x): V . DiagonalMatrix[f /@ vals] .
 * Inverse[V].  NULL when the matrix is not diagonalisable. */
Expr* mf_eigen_apply(const Expr* m, const Expr* f);

#endif /* MATHILDA_LINALG_MATFUNC_INTERNAL_H */
