/* Mathilda -- sparse linear algebra: the band backend of SpFactor
 * (sparse_la_band.c).  RCM-ordered LAPACK band LU / Cholesky for the builds
 * without SuiteSparse.  Internal to src/sparse/. */
#ifndef MATHILDA_SPARSE_LA_BAND_H
#define MATHILDA_SPARSE_LA_BAND_H

#include "sparse_la.h"

typedef struct SpBand SpBand;

typedef enum {
    SPB_LU = 0,     /* the band LU                                             */
    SPB_CHOL_LU,    /* A is Hermitian with a positive diagonal (the caller
                     * tested it, spla_hermitian_posdiag): the band Cholesky,
                     * the LU when it is not positive definite                  */
    SPB_CHOL        /* the band Cholesky only (SPF_NOTPD when it fails); A must
                     * be Hermitian                                              */
} SpbMode;

/* Factor A in band form after a symmetric RCM renumbering.  NULL with
 * *st == SPF_TOOBIG when the band would not pay (wider than the cost model
 * allows, or no LAPACK) -- the caller's dense path then answers -- or with
 * SPF_SINGULAR / SPF_NOTPD / SPF_NOMEM / SPF_FAILED. */
SpBand* spb_factor(const SpMat* A, SpbMode mode, SpfStatus* st);
bool spb_solve(const SpBand* b, const double* B, double* X, int64_t nrhs, bool trans);
bool spb_det(const SpBand* b, double* re, double* im);
double spb_rcond(const SpBand* b);
bool spb_is_chol(const SpBand* b);
void spb_free(SpBand* b);

#endif
