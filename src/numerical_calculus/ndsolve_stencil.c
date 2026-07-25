/* Mathilda — finite-difference stencil weights (Fornberg).  See header. */
#include "ndsolve_stencil.h"
#include <stdlib.h>
#include <math.h>

void nd_fornberg_weights(int m, const double* x, int npts, double z, double* w) {
    int N = npts - 1;
    /* c[node*(m+1) + k] = weight of node for the k-th derivative (k = 0..m). */
    double* c = calloc((size_t)npts * (size_t)(m + 1), sizeof(double));
    if (!c) { for (int i = 0; i < npts; i++) w[i] = 0.0; return; }
#define C(i, k) c[(size_t)(i) * (size_t)(m + 1) + (size_t)(k)]
    double c1 = 1.0;
    double c4 = x[0] - z;
    C(0, 0) = 1.0;
    for (int n = 1; n <= N; n++) {
        int mn = n < m ? n : m;
        double c2 = 1.0;
        double c5 = c4;
        c4 = x[n] - z;
        for (int nu = 0; nu <= n - 1; nu++) {
            double c3 = x[n] - x[nu];
            c2 *= c3;
            if (nu == n - 1) {
                for (int k = mn; k >= 1; k--)
                    C(n, k) = c1 * ((double)k * C(n - 1, k - 1) - c5 * C(n - 1, k)) / c2;
                C(n, 0) = -c1 * c5 * C(n - 1, 0) / c2;
            }
            for (int k = mn; k >= 1; k--)
                C(nu, k) = (c4 * C(nu, k) - (double)k * C(nu, k - 1)) / c3;
            C(nu, 0) = c4 * C(nu, 0) / c3;
        }
        c1 = c2;
    }
    for (int i = 0; i < npts; i++) w[i] = C(i, m);
#undef C
    free(c);
}

void nd_stencil_build(int j, int nx, int deriv, int order, double h,
                      int* idx, double* w, int* n) {
    /* Points needed: for deriv 1 and 2 a symmetric (order+1)-point central
     * stencil reaches the requested even order (symmetry supplies the extra
     * order); for higher derivatives fall back to the generic deriv+order. */
    int npts = (deriv <= 2) ? (order + 1) : (deriv + order);
    if (npts < deriv + 1) npts = deriv + 1;   /* enough nodes to form ∂^deriv  */
    if (npts > nx) npts = nx;                  /* clamp on very coarse grids    */

    int half = npts / 2;
    int lo = j - half;                         /* centre on j when it fits      */
    if (lo < 0) lo = 0;
    if (lo + npts > nx) lo = nx - npts;

    double xs[64] = { 0 };
    if (npts > 64) npts = 64;                  /* defensive; buffers are 64     */
    for (int i = 0; i < npts; i++) { idx[i] = lo + i; xs[i] = (double)(lo + i - j); }
    nd_fornberg_weights(deriv, xs, npts, 0.0, w);

    double hd = pow(h, (double)deriv);
    for (int i = 0; i < npts; i++) w[i] /= hd;
    *n = npts;
}
