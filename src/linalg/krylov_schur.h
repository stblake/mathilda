/* Mathilda -- restarted Krylov eigensolvers on a matrix-vector callback.
 *
 * ks_solve finds `nev` eigenpairs (theta, x) of a linear operator OP, given
 * only y = OP x, selected by a criterion (largest magnitude, largest real
 * part, both ends of a Hermitian spectrum, ...). OP is whatever the caller
 * makes it: A itself, a shift-and-invert (A - sigma I)^-1 applied through a
 * factorisation, B^-1 A for a pencil, A^H A for singular values, a dense
 * matrix. The machine-precision eigen and singular-value heads with a count
 * (Eigenvalues[m, k], SingularValueList[s, k], Method -> "Arnoldi") all reach
 * the Krylov layer through this one interface.
 *
 * Two backends:
 *
 *   ARPACK  (USE_ARPACK) the implicitly restarted Lanczos / Arnoldi of
 *           arpack-ng -- dsaupd for real symmetric OP, dnaupd for real
 *           general, znaupd for complex (Hermitian included) -- in regular
 *           mode 1 on the callback. Used when linked, unless the environment
 *           variable MATHILDA_NO_ARPACK is set (to anything but "0") or
 *           ks_set_backend(KS_BACKEND_NATIVE) was called (the tests force
 *           the native path that way), and when ARPACK's own size limits
 *           (nev < n - 1, ncv <= n) allow the problem. A run that does not
 *           converge falls through to the native solver.
 *
 *   native  Krylov-Schur (Stewart 2001) with full re-orthogonalisation
 *           (classical Gram-Schmidt twice, "twice is enough"). The basis
 *           V_m of an Arnoldi decomposition A V_m = V_m H_m + beta v e_m^T
 *           is expanded to m = ncv vectors; the Ritz pairs of the small
 *           projected matrix H_m are computed with LAPACK (dsyev / zheev for
 *           a Hermitian OP, dgeev / zgeev otherwise); the wanted ones are
 *           kept by a thick restart onto an orthonormal basis Q of their
 *           invariant subspace -- the Ritz vectors themselves when OP is
 *           Hermitian (thick-restart Lanczos: H_k = diag(theta)), an
 *           orthonormalised set of Ritz vectors (real and imaginary parts
 *           of a conjugate pair together, so the basis stays real for a
 *           real OP) otherwise -- giving A (V_m Q) = (V_m Q)(Q^H H_m Q) +
 *           beta v (e_m^T Q), a Krylov decomposition again, which the next
 *           expansion continues. A pair is converged when its Ritz estimate
 *           beta |e_m^T y| <= tol max(|theta|, eps^(2/3) max|theta|), and
 *           then confirmed by its explicit residual ||OP x - theta x||
 *           before the run stops. Real OP stays in real arithmetic
 *           (conjugate pairs come out exactly conjugate).
 *
 * Without LAPACK (USE_LAPACK=0, or the runtime probe failing) the native
 * solver reports KS_FAIL and the callers fall back to their dense paths.
 *
 * Vectors are interleaved (re, im) doubles when cplx; real otherwise.
 */
#ifndef MATHILDA_KRYLOV_SCHUR_H
#define MATHILDA_KRYLOV_SCHUR_H

#include <stdbool.h>
#include <stdint.h>

/* y = OP x; x and y have n (real) or 2n (complex) doubles and do not
 * alias. Return false to abort the run (KS_FAIL). */
typedef bool (*KsOpFn)(void* ctx, const double* x, double* y);

typedef enum {
    KS_LM = 0,  /* largest magnitude                                      */
    KS_SM,      /* smallest magnitude (of OP -- slow; shift-invert first) */
    KS_LR,      /* largest real part                                      */
    KS_SR,      /* smallest real part                                     */
    KS_LI,      /* largest |imaginary part|                               */
    KS_SI,      /* smallest |imaginary part|                              */
    KS_BE       /* both ends of a Hermitian spectrum: nev - nev/2 from the
                 * top, nev/2 from the bottom                             */
} KsWhich;

typedef enum { KS_BACKEND_AUTO = 0, KS_BACKEND_NATIVE, KS_BACKEND_ARPACK } KsBackend;

typedef struct {
    int64_t n;
    bool cplx;          /* complex vectors                                  */
    bool hermitian;     /* OP self-adjoint (Euclidean inner product)        */
    KsOpFn op;
    void* ctx;
    int nev;            /* 1 <= nev <= n                                    */
    KsWhich which;
    int ncv;            /* basis size; <= 0: max(2 nev + 1, 20), capped at n */
    int64_t maxit;      /* restarts; <= 0: 1000                             */
    double tol;         /* relative; <= 0: machine epsilon                  */
    const double* v0;   /* starting vector (n or 2n doubles) or NULL         */
    bool want_vectors;
} KsProblem;

typedef struct {
    int nconv;          /* pairs returned, in `which` order (best first)    */
    double* re;         /* eigenvalues of OP                                */
    double* im;
    double* vec;        /* nconv unit vectors, one after another: n doubles
                         * each when !vec_cplx, else 2n interleaved         */
    bool vec_cplx;
    int64_t iters;      /* restarts used                                    */
    const char* backend;/* "ARPACK" or "KrylovSchur"                         */
} KsResult;

typedef enum {
    KS_OK = 0,
    KS_MAXIT,   /* the iteration limit was reached: the nev best current
                 * pairs are returned anyway                              */
    KS_FAIL,    /* no answer (no LAPACK, OP failed, bad arguments)         */
    KS_NOMEM
} KsStatus;

KsStatus ks_solve(const KsProblem* p, KsResult* r);
void ks_result_free(KsResult* r);

/* Force a backend for the rest of the process (tests); AUTO restores the
 * default. Returns the previous setting. */
KsBackend ks_set_backend(KsBackend b);
/* The backend AUTO resolves to now: ARPACK when linked and not disabled. */
bool ks_arpack_active(void);

#endif /* MATHILDA_KRYLOV_SCHUR_H */
