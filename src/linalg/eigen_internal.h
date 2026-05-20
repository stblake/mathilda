#ifndef EIGEN_INTERNAL_H
#define EIGEN_INTERNAL_H

/* Internal shared interface for the eigen_*.c files produced by the
 * Phase 1b split of the former monolithic mateigen.c.
 *
 *   eigen.c          -- builtin entry points + top-level dispatch
 *   eigen_common.c   -- symbolic helpers, MatD / MatM, options, sort/perm
 *   eigen_direct.c   -- Hessenberg + QR (machine + MPFR)
 *   eigen_arnoldi.c  -- Krylov (machine + MPFR)
 *   eigen_banded.c   -- band -> tridiag -> QR (machine + MPFR)
 *   eigen_feast.c    -- contour-integral spectral projector (machine + MPFR)
 *
 * Everything declared here is internal to the eigen translation units
 * and must NOT leak into any other compilation unit.  The public
 * interface lives in eigen.h. */

#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

#include "expr.h"
#include "eigen.h"

/* ---------------------------------------------------------------------
 * Shared option struct: parsed positional + named arguments for the
 * Eigenvalues / Eigenvectors builtins.  Produced by eigen_parse_args
 * (defined in eigen_common.c) and consumed by the dispatchers in
 * eigen.c.
 * ------------------------------------------------------------------ */
typedef struct {
    Expr* arg0;             /* m or {m, a}                          */
    Expr* k_spec;           /* Integer k, or UpTo[k], or NULL       */
    bool  cubics_radical;
    bool  quartics_radical;
    bool  method_given;     /* user supplied Method -> ...          */
    MateigenMethod method;  /* parsed Method, or AUTOMATIC if unset */
    Expr* method_value;     /* original Method RHS (for sub-options) */
} EigenOpts;

/* ---------------------------------------------------------------------
 * Machine-precision dense matrix workspace.  Row-major.  Always
 * heap-allocated; matD_free must be called regardless of how the matrix
 * was populated.
 * ------------------------------------------------------------------ */
typedef struct {
    size_t  n;
    int     is_complex;   /* 0: real only, 1: complex (im[] non-NULL) */
    double* re;           /* length n*n */
    double* im;           /* length n*n, NULL when !is_complex */
} MatD;

#ifdef USE_MPFR
/* Arbitrary-precision dense matrix workspace.  Row-major.  Each
 * mpfr_t in `re` / `im` is pre-initialised to `bits` on alloc and
 * cleared by matM_free.  Always heap-allocated. */
typedef struct {
    size_t      n;
    int         is_complex;   /* 0: real only, 1: complex (im[] non-NULL) */
    mpfr_prec_t bits;
    mpfr_t*     re;           /* length n*n, mpfr_init2'd to bits */
    mpfr_t*     im;           /* length n*n, mpfr_init2'd, NULL when !is_complex */
} MatM;
#endif /* USE_MPFR */

/* =====================================================================
 * eigen_common.c -- shared helpers used across the per-method files.
 * ================================================================== */

/* Symbolic / option helpers used by eigen.c and the per-method files. */
bool        eigen_parse_args(Expr* res, EigenOpts* opts);
bool        eigen_extract_matrix_pair(Expr* arg, Expr** m_out, Expr** a_out,
                                       int64_t* n_out);
bool        eigen_matrix_is_inexact(Expr* m);
Expr*       eigen_chop(Expr* val);
void        eigen_sort_by_abs_desc(Expr** vals, size_t n);
Expr**      eigen_null_space(Expr* M, int n, size_t* count_out);
const char* eigen_lambda_name(void);
Expr*       eigen_build_lambda_matrix(Expr* m, Expr* a_or_null,
                                       const char* lambda_name, int64_t n);
Expr*       eigen_compute_det(Expr* matrix, int n);
Expr*       eigen_char_poly_faddeev(Expr* A, const char* lambda_name, int n);
Expr*       eigen_solve_poly(Expr* poly, const char* lambda_name,
                              bool cubics_radical, bool quartics_radical);
Expr**      eigen_extract_values(Expr* solutions, size_t* count_out);

/* MatD helpers (used by direct/arnoldi/banded/feast). */
void   matD_free(MatD* M);
bool   matD_load(Expr* m, size_t n, MatD* out);
double matD_norm_inf_real(const double* A, size_t n);
double matD_norm_inf_complex(const MatD* A);
bool   matD_is_real_symmetric(const MatD* A, double tol);
bool   matD_is_hermitian(const MatD* A, double tol);
size_t matD_bandwidth_real(const double* A, size_t n, double tol);
size_t matD_bandwidth_complex(const MatD* A, double tol);

/* Build a length-k List of complex eigenvectors from a stacked V_re/V_im
 * column-block.  Used by arnoldi_*_machine which materialises a
 * rectangular eigvec matrix rather than a square Q. */
Expr* mateigen_build_complex_eigvec_list_rect(const double* V_re,
                                               const double* V_im,
                                               size_t k_vectors,
                                               size_t n_components);

#ifdef USE_MPFR
/* MatM helpers (used by direct/arnoldi/banded/feast MPFR kernels). */
void    mpfr_array_free(mpfr_t* a, size_t count);
mpfr_t* mpfr_array_alloc(size_t count, mpfr_prec_t bits);
void    matM_free(MatM* M);
bool    matM_load(Expr* m, size_t n, mpfr_prec_t bits, MatM* out);
void    matM_norm_inf_real(const mpfr_t* A, size_t n,
                            mpfr_prec_t bits, mpfr_t out);
void    matM_norm_inf_complex(const MatM* A, mpfr_t out);
bool    matM_is_real_symmetric(const MatM* A, const mpfr_t tol);
bool    matM_is_hermitian(const MatM* A, const mpfr_t tol);
size_t  matM_bandwidth_real(const mpfr_t* A, size_t n, const mpfr_t tol);
size_t  matM_bandwidth_complex(const MatM* A, const mpfr_t tol);
Expr*   mateigen_build_complex_eigvec_list_rect_M(const mpfr_t* V_re,
                                                   const mpfr_t* V_im,
                                                   size_t k_vectors,
                                                   size_t n_components);
#endif /* USE_MPFR */

/* =====================================================================
 * eigen_direct.c -- Hessenberg + QR.
 *
 * Re-used as a building block by eigen_arnoldi.c, eigen_banded.c and
 * eigen_feast.c, which all funnel their reduced problems through the
 * same QR machinery.  The top-level direct_dispatch is called from
 * eigen.c.
 * ================================================================== */

/* Top-level dispatcher: routes inexact ordinary eigenproblems through
 * the Direct kernel matched to the leaf precision.  Returns NULL when
 * the matrix shape isn't yet wired (caller falls back to symbolic). */
Expr* direct_dispatch(Expr* m, Expr* a, int64_t n,
                       MateigenWant want, Expr* k_spec);

/* Trim a full-length result List against a k-spec (Integer k, -k, or
 * UpTo[k]); used by every per-method kernel that has already produced
 * a complete list. */
Expr* direct_apply_k_spec_list(Expr* full_list, Expr* k_spec);

/* Sort permutation helpers (descending |lambda|), shared by every
 * dispatcher that needs to fix up Ritz-pair order. */
void  direct_sort_perm_desc_abs(const double* vals, size_t n, size_t* perm);
void  direct_sort_perm_desc_abs_complex(const double* re, const double* im,
                                         size_t n, size_t* perm);

/* List builders shared by band / feast / arnoldi paths. */
Expr* direct_build_real_eigenvalue_list(const double* vals, size_t n,
                                         const size_t* perm);
Expr* direct_build_real_eigenvector_list(const double* Q, size_t n,
                                          const size_t* perm);
Expr* direct_build_complex_eigenvalue_list(const double* re, const double* im,
                                            size_t n, const size_t* perm);
Expr* direct_build_complex_hermitian_eigvec_list(const double* V_re,
                                                  const double* V_im,
                                                  size_t n,
                                                  const size_t* perm);

/* Tridiagonalisation + Francis QR re-used by arnoldi / banded / feast. */
void direct_hessenberg_real(double* A, size_t n, double* u, double* Q);
int  direct_qr_real_general(double* H, size_t n,
                             double* eval_re, double* eval_im,
                             double* Q);
int  direct_symtridiag_qr(double* diag, double* sub, size_t n,
                           double* Q, bool want_Q);
void direct_tridiag_real_sym(double* A, size_t n,
                              double* diag, double* subdiag, double* Q,
                              bool want_Q,
                              double* u, double* p, double* q);
void direct_tridiag_complex_hermitian(double* A_re, double* A_im,
                                       size_t n,
                                       double* diag,
                                       double* sub_re, double* sub_im,
                                       double* Q_re, double* Q_im,
                                       bool want_Q,
                                       double* u_re, double* u_im,
                                       double* v_re, double* v_im,
                                       double* q_re, double* q_im);
void direct_phase_correct_tridiag(double* sub_re, double* sub_im,
                                   size_t n,
                                   double* Q_re, double* Q_im,
                                   bool want_Q);
void direct_compose_complex_Q_real_Z(const double* Q_re, const double* Q_im,
                                      const double* Z, size_t n,
                                      double* V_re, double* V_im);

/* Schur back-substitution (real + complex), shared by arnoldi. */
void schur_compute_eigvecs(const double* H, const double* Q, size_t n,
                            const double* eval_re, const double* eval_im,
                            const size_t* perm,
                            double* V_re, double* V_im);

#ifdef USE_MPFR
/* MPFR mirrors of the above. */
void  direct_sort_perm_desc_abs_M(const mpfr_t* vals, size_t n, size_t* perm);
void  direct_sort_perm_desc_abs_complex_M(const mpfr_t* re, const mpfr_t* im,
                                           size_t n, size_t* perm);
Expr* direct_build_real_eigenvalue_list_M(const mpfr_t* vals, size_t n,
                                           const size_t* perm);
Expr* direct_build_real_eigenvector_list_M(const mpfr_t* Q, size_t n,
                                            const size_t* perm);
Expr* direct_build_complex_eigenvalue_list_M(const mpfr_t* re,
                                              const mpfr_t* im,
                                              size_t n, const size_t* perm);
Expr* direct_build_complex_hermitian_eigvec_list_M(const mpfr_t* V_re,
                                                    const mpfr_t* V_im,
                                                    size_t n,
                                                    const size_t* perm);

void direct_hessenberg_real_M(mpfr_t* A, size_t n, mpfr_prec_t bits,
                               mpfr_t* u, mpfr_t* Q, mpfr_t* tmp);
int  direct_qr_real_general_M(mpfr_t* H, size_t n, mpfr_prec_t bits,
                               mpfr_t* eval_re, mpfr_t* eval_im,
                               mpfr_t* Q, mpfr_t* tmp);
int  direct_symtridiag_qr_M(mpfr_t* diag, mpfr_t* sub, size_t n,
                             mpfr_prec_t bits,
                             mpfr_t* Q, bool want_Q, mpfr_t* tmp);
void direct_tridiag_real_sym_M(mpfr_t* A, size_t n, mpfr_prec_t bits,
                                mpfr_t* diag, mpfr_t* sub, mpfr_t* Q,
                                bool want_Q,
                                mpfr_t* u, mpfr_t* p, mpfr_t* q,
                                mpfr_t* tmp);
void direct_tridiag_complex_hermitian_M(mpfr_t* A_re, mpfr_t* A_im,
                                         size_t n, mpfr_prec_t bits,
                                         mpfr_t* diag,
                                         mpfr_t* sub_re, mpfr_t* sub_im,
                                         mpfr_t* Q_re, mpfr_t* Q_im,
                                         bool want_Q,
                                         mpfr_t* u_re, mpfr_t* u_im,
                                         mpfr_t* v_re, mpfr_t* v_im,
                                         mpfr_t* q_re, mpfr_t* q_im);
void direct_phase_correct_tridiag_M(mpfr_t* sub_re, mpfr_t* sub_im,
                                     size_t n, mpfr_prec_t bits,
                                     mpfr_t* Q_re, mpfr_t* Q_im,
                                     bool want_Q);
void direct_compose_complex_Q_real_Z_M(const mpfr_t* Q_re, const mpfr_t* Q_im,
                                        const mpfr_t* Z, size_t n,
                                        mpfr_prec_t bits,
                                        mpfr_t* V_re, mpfr_t* V_im);
void schur_compute_eigvecs_M(const mpfr_t* H, const mpfr_t* Q, size_t n,
                              mpfr_prec_t bits,
                              const mpfr_t* eval_re, const mpfr_t* eval_im,
                              const size_t* perm,
                              mpfr_t* V_re, mpfr_t* V_im);
#endif /* USE_MPFR */

/* =====================================================================
 * eigen_arnoldi.c -- Krylov-subspace iteration.
 * ================================================================== */
Expr* arnoldi_dispatch(Expr* m, Expr* a, int64_t n,
                        MateigenWant want, Expr* k_spec,
                        Expr* method_value);
bool  arnoldi_automatic_prefers(Expr* k_spec, size_t n);

/* =====================================================================
 * eigen_banded.c -- band -> tridiag -> QR (Hermitian only).
 * ================================================================== */
Expr* banded_dispatch(Expr* m, Expr* a, int64_t n,
                       MateigenWant want, Expr* k_spec,
                       Expr* method_value);
bool  banded_automatic_prefers(Expr* m, int64_t n);

/* =====================================================================
 * eigen_feast.c -- spectral-projector contour integral (Hermitian only).
 * ================================================================== */
Expr* feast_dispatch(Expr* m, Expr* a, int64_t n,
                      MateigenWant want, Expr* k_spec,
                      Expr* method_value);
bool  feast_automatic_prefers(Expr* m, int64_t n, Expr* method_value);

#endif /* EIGEN_INTERNAL_H */
