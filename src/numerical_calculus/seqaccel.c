/*
 * seqaccel.c — shared sequence-acceleration kernels (see seqaccel.h)
 *
 * Extracted from nlimit.c so that NLimit and NSum share one implementation of
 * Richardson/Romberg ("EulerSum") and Wynn's epsilon ("SequenceLimit").  The
 * kernels are pure numeric: they neither build Expr nodes nor apply any
 * acceptance policy — the caller owns the convergence gate.
 */

#include "seqaccel.h"

#include <math.h>
#include <stdlib.h>

#ifdef USE_MPFR
#  include "numeric_complex.h"   /* mpfr_complex_div */
#endif

/* ------------------------------------------------------------------ *
 *  Machine precision                                                  *
 * ------------------------------------------------------------------ */

bool seqaccel_richardson_machine(const double _Complex* S, int terms,
                                 double _Complex* result, double* step) {
    if (terms < 1) return false;
    double _Complex* T = malloc(sizeof(double _Complex) * (size_t)terms * terms);
    if (!T) return false;
    for (int i = 0; i < terms; i++) T[(size_t)i * terms + 0] = S[i];
    for (int j = 1; j < terms; j++) {
        double denom = ldexp(1.0, j) - 1.0;        /* 2^j - 1 */
        for (int i = j; i < terms; i++) {
            double _Complex a = T[(size_t)i * terms + (j - 1)];
            double _Complex b = T[(size_t)(i - 1) * terms + (j - 1)];
            T[(size_t)i * terms + j] = a + (a - b) / denom;
        }
    }
    double _Complex last = T[(size_t)(terms - 1) * terms + (terms - 1)];
    double _Complex prev = (terms >= 2)
        ? T[(size_t)(terms - 2) * terms + (terms - 2)] : last;
    *result = last;
    *step = cabs(last - prev);
    free(T);
    return true;
}

bool seqaccel_wynn_machine(const double _Complex* S, int terms, int degree,
                           double _Complex* result, double* step) {
    int maxdeg = (terms - 1) / 2;
    if (maxdeg < 1) return false;          /* need >= 3 terms for one Shanks */
    if (degree > maxdeg) degree = maxdeg;
    if (degree < 1) degree = 1;

    int cols = terms + 1;
    size_t stride = (size_t)cols;
    double _Complex* eps = malloc(sizeof(double _Complex) * stride * stride);
    if (!eps) return false;
    for (size_t i = 0; i < stride * stride; i++) eps[i] = 0.0;
#define SA_E(c, n) eps[(size_t)(c) * stride + (size_t)(n)]
    for (int n = 0; n < terms; n++) SA_E(1, n) = S[n];   /* ε_0 */

    for (int c = 2; c <= terms; c++) {
        int len = terms + 1 - c;
        for (int n = 0; n < len; n++) {
            double _Complex d = SA_E(c - 1, n + 1) - SA_E(c - 1, n);
            double _Complex recip;
            if (cabs(d) <= 0.0) recip = 1e300;          /* singular bridge */
            else                recip = 1.0 / d;
            SA_E(c, n) = SA_E(c - 2, n + 1) + recip;
        }
    }

    /* Read from the deepest even column ε_{2d} (column 2d+1 shifted); pick the
     * entry that best agrees with its neighbour, where the algorithm has
     * converged, rather than the roundoff-amplified bottom corner. */
    int rc = 2 * degree + 1;               /* result column (ε_{2d}) */
    int rl = terms + 1 - rc;               /* its length */
    double _Complex res = SA_E(rc, rl - 1);
    double best_d = INFINITY;
    bool found = false;
    for (int n = 1; n < rl; n++) {
        double _Complex a = SA_E(rc, n), b = SA_E(rc, n - 1);
        if (!isfinite(creal(a)) || !isfinite(cimag(a))) continue;
        if (!isfinite(creal(b)) || !isfinite(cimag(b))) continue;
        double d = cabs(a - b);
        if (d < best_d) { best_d = d; res = a; found = true; }
    }
    if (!found) {                           /* single-entry column fallback */
        res = SA_E(rc, rl - 1);
        double _Complex prev = (rc >= 3)
            ? SA_E(rc - 2, terms + 1 - (rc - 2) - 1) : SA_E(1, terms - 1);
        best_d = cabs(res - prev);
    }

    *result = res;
    *step = best_d;
#undef SA_E
    free(eps);
    return true;
}

/* ------------------------------------------------------------------ *
 *  MPFR precision                                                     *
 * ------------------------------------------------------------------ */

#ifdef USE_MPFR
/* L1 magnitude |re| + |im| as a double, for the caller's gauges. */
static double sa_l1_d(const mpfr_t re, const mpfr_t im) {
    return fabs(mpfr_get_d(re, MPFR_RNDN)) + fabs(mpfr_get_d(im, MPFR_RNDN));
}

bool seqaccel_richardson_mpfr(const mpfr_t* Sr, const mpfr_t* Si, int terms,
                              long bits, mpfr_t out_re, mpfr_t out_im,
                              double* step, bool* finite) {
    if (terms < 1) return false;
    mpfr_prec_t p = (mpfr_prec_t)bits;
    size_t cells = (size_t)terms * terms;
    mpfr_t* Tr = malloc(sizeof(mpfr_t) * cells);
    mpfr_t* Ti = malloc(sizeof(mpfr_t) * cells);
    if (!Tr || !Ti) { free(Tr); free(Ti); return false; }
    for (size_t c = 0; c < cells; c++) { mpfr_init2(Tr[c], p); mpfr_init2(Ti[c], p); }
    mpfr_t tr, ti, denom;
    mpfr_inits2(p, tr, ti, denom, (mpfr_ptr)0);

    for (int i = 0; i < terms; i++) {
        mpfr_set(Tr[(size_t)i * terms + 0], Sr[i], MPFR_RNDN);
        mpfr_set(Ti[(size_t)i * terms + 0], Si[i], MPFR_RNDN);
    }
    for (int j = 1; j < terms; j++) {
        mpfr_set_ui(denom, 1, MPFR_RNDN);
        mpfr_mul_2ui(denom, denom, (unsigned long)j, MPFR_RNDN);
        mpfr_sub_ui(denom, denom, 1, MPFR_RNDN);
        for (int i = j; i < terms; i++) {
            size_t cur = (size_t)i * terms + j;
            size_t a = (size_t)i * terms + (j - 1);
            size_t b = (size_t)(i - 1) * terms + (j - 1);
            mpfr_sub(tr, Tr[a], Tr[b], MPFR_RNDN);
            mpfr_sub(ti, Ti[a], Ti[b], MPFR_RNDN);
            mpfr_div(tr, tr, denom, MPFR_RNDN);
            mpfr_div(ti, ti, denom, MPFR_RNDN);
            mpfr_add(Tr[cur], Tr[a], tr, MPFR_RNDN);
            mpfr_add(Ti[cur], Ti[a], ti, MPFR_RNDN);
        }
    }

    size_t last = (size_t)(terms - 1) * terms + (terms - 1);
    size_t prev = (terms >= 2)
        ? (size_t)(terms - 2) * terms + (terms - 2) : last;
    mpfr_sub(tr, Tr[last], Tr[prev], MPFR_RNDN);
    mpfr_sub(ti, Ti[last], Ti[prev], MPFR_RNDN);
    *step = sa_l1_d(tr, ti);
    *finite = mpfr_number_p(Tr[last]) && mpfr_number_p(Ti[last]);
    mpfr_set(out_re, Tr[last], MPFR_RNDN);
    mpfr_set(out_im, Ti[last], MPFR_RNDN);

    for (size_t c = 0; c < cells; c++) { mpfr_clear(Tr[c]); mpfr_clear(Ti[c]); }
    free(Tr); free(Ti);
    mpfr_clears(tr, ti, denom, (mpfr_ptr)0);
    return true;
}

bool seqaccel_wynn_mpfr(const mpfr_t* Sr, const mpfr_t* Si, int terms,
                        int degree, long bits, mpfr_t out_re, mpfr_t out_im,
                        double* step, bool* finite) {
    int maxdeg = (terms - 1) / 2;
    if (maxdeg < 1) return false;
    if (degree > maxdeg) degree = maxdeg;
    if (degree < 1) degree = 1;

    mpfr_prec_t p = (mpfr_prec_t)bits;
    int cols = terms + 1;
    size_t stride = (size_t)cols;
    size_t cells = stride * stride;
    mpfr_t* er = malloc(sizeof(mpfr_t) * cells);
    mpfr_t* ei = malloc(sizeof(mpfr_t) * cells);
    if (!er || !ei) { free(er); free(ei); return false; }
    for (size_t i = 0; i < cells; i++) {
        mpfr_init2(er[i], p); mpfr_init2(ei[i], p);
        mpfr_set_ui(er[i], 0, MPFR_RNDN); mpfr_set_ui(ei[i], 0, MPFR_RNDN);
    }
    mpfr_t dr, di, rr, ri, one, zero;
    mpfr_inits2(p, dr, di, rr, ri, one, zero, (mpfr_ptr)0);
    mpfr_set_ui(one, 1, MPFR_RNDN);
    mpfr_set_ui(zero, 0, MPFR_RNDN);
#define SA_ER(c, n) er[(size_t)(c) * stride + (size_t)(n)]
#define SA_EI(c, n) ei[(size_t)(c) * stride + (size_t)(n)]
    for (int n = 0; n < terms; n++) {
        mpfr_set(SA_ER(1, n), Sr[n], MPFR_RNDN);
        mpfr_set(SA_EI(1, n), Si[n], MPFR_RNDN);
    }
    for (int c = 2; c <= terms; c++) {
        int len = terms + 1 - c;
        for (int n = 0; n < len; n++) {
            mpfr_sub(dr, SA_ER(c - 1, n + 1), SA_ER(c - 1, n), MPFR_RNDN);
            mpfr_sub(di, SA_EI(c - 1, n + 1), SA_EI(c - 1, n), MPFR_RNDN);
            if (mpfr_zero_p(dr) && mpfr_zero_p(di)) {
                mpfr_set_d(rr, 1e300, MPFR_RNDN);
                mpfr_set_ui(ri, 0, MPFR_RNDN);
            } else {
                mpfr_complex_div(rr, ri, one, zero, dr, di);
            }
            mpfr_add(SA_ER(c, n), SA_ER(c - 2, n + 1), rr, MPFR_RNDN);
            mpfr_add(SA_EI(c, n), SA_EI(c - 2, n + 1), ri, MPFR_RNDN);
        }
    }

    /* Pick the ε_{2d} entry that best agrees with its neighbour. */
    int rc = 2 * degree + 1;
    int rl = terms + 1 - rc;
    int best_n = rl - 1;
    double best_d = -1.0;
    for (int n = 1; n < rl; n++) {
        if (!mpfr_number_p(SA_ER(rc, n)) || !mpfr_number_p(SA_EI(rc, n))) continue;
        if (!mpfr_number_p(SA_ER(rc, n - 1)) || !mpfr_number_p(SA_EI(rc, n - 1))) continue;
        mpfr_sub(dr, SA_ER(rc, n), SA_ER(rc, n - 1), MPFR_RNDN);
        mpfr_sub(di, SA_EI(rc, n), SA_EI(rc, n - 1), MPFR_RNDN);
        double d = sa_l1_d(dr, di);
        if (best_d < 0.0 || d < best_d) { best_d = d; best_n = n; }
    }
    if (best_d < 0.0) {                      /* single-entry fallback */
        int pc = (rc >= 3) ? rc - 2 : 1;
        int pn = (rc >= 3) ? terms + 1 - pc - 1 : terms - 1;
        mpfr_sub(dr, SA_ER(rc, rl - 1), SA_ER(pc, pn), MPFR_RNDN);
        mpfr_sub(di, SA_EI(rc, rl - 1), SA_EI(pc, pn), MPFR_RNDN);
        best_d = sa_l1_d(dr, di);
        best_n = rl - 1;
    }
    *step = best_d;
    *finite = mpfr_number_p(SA_ER(rc, best_n)) && mpfr_number_p(SA_EI(rc, best_n));
    mpfr_set(out_re, SA_ER(rc, best_n), MPFR_RNDN);
    mpfr_set(out_im, SA_EI(rc, best_n), MPFR_RNDN);

#undef SA_ER
#undef SA_EI
    for (size_t i = 0; i < cells; i++) { mpfr_clear(er[i]); mpfr_clear(ei[i]); }
    free(er); free(ei);
    mpfr_clears(dr, di, rr, ri, one, zero, (mpfr_ptr)0);
    return true;
}
#endif /* USE_MPFR */
