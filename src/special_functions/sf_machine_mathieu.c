/* Mathilda -- machine-precision Mathieu functions (MathieuC / MathieuCPrime,
 * MathieuS / MathieuSPrime).
 *
 * See sf_machine_mathieu.h for the definition and normalization. This file holds
 * the double-precision algorithm shared by the scalar builtin, the packed /
 * NDArray buffer loop, the Compile[] VM, and (as a starting point) the
 * arbitrary-precision kernel in mathieu_arb.c.
 *
 * ALGORITHM -- one "plan" per (a, q), then a cheap sum per z.
 *
 *  1. Characteristic exponent, first estimate: Hill's determinant,
 *         sin^2(pi nu / 2) = Delta(0) sin^2(pi sqrt(a) / 2),
 *     rewritten without the poles Delta(0) has at a = (2k)^2: with the rows of the
 *     Hill matrix rescaled by (4k^2 - a)/(4k^2), the k = 0 row by a,
 *         sin^2(pi nu/2) = (pi^2/4) det(T) R^2,
 *         det(T) = D_1 (a D_1 + q^2/2 D_2)        (T symmetric about k = 0),
 *         D_k = (1 - a/4k^2) D_{k+1} - q^2/(16 k^2 (k+1)^2) D_{k+2},
 *     and R = [sin(pi x)/(pi x)] / prod_{k<=K} (1 - x^2/k^2), x = sqrt(a)/2,
 *     evaluated with the near-zero factor divided out analytically (sinc of the
 *     offset from the nearest integer), so a = (2k)^2 is a regular point. The
 *     truncated determinant converges only polynomially (the tail is
 *     exp(-sum q^2/(16k^4)) to first order), so the first two orders of the tail
 *     are added in closed form (trigamma asymptotics) and K is chosen from the
 *     size of the next one.
 *  2. nu = +-acos(cos(pi nu))/pi + 2m, the representative nearest sqrt(a) (so
 *     q -> 0 gives nu -> sqrt(a) and c_0 -> 1). Every such choice describes the
 *     same function; the choice only fixes which coefficient is c_0.
 *  3. |nu - n| < 2e-6 for an integer n: a is (to machine precision) a
 *     characteristic value -> nu = n exactly, DLMF normalization.
 *  4. Otherwise nu is refined by Newton on the consistency function
 *         g(nu) = (a - (nu+2m)^2) - q (c_{m+1}/c_m + c_{m-1}/c_m),
 *     the two ratios from continued fractions run in from +-infinity. Those
 *     converge super-exponentially (ratio ~ q/(4k^2)), so g carries no truncation
 *     error; every root of g is an equivalent exponent (+-nu + 2j), so Newton
 *     landing on any of them is harmless. m is the dominant coefficient index,
 *     found by two steps of inverse iteration on the truncated tridiagonal
 *     Floquet matrix (stable for any q, unlike a one-sided recurrence).
 *  5. Coefficients: the same two continued fractions, glued at m; normalized by
 *     sum c_k^2 = 1, sign Re(c_0) > 0.
 *  6. Evaluation at z: with w = e^{2iz}, theta = nu + 2 lo,
 *         ce(z)  = ( e^{i theta z} P(w) + e^{-i theta z} P(1/w) ) / 2,
 *         ce'(z) = (i/2) ( e^{i theta z} Pw(w) - e^{-i theta z} Pw(1/w) ),
 *     P(w) = sum_j c_{lo+j} w^j (Horner), Pw with weights (nu+2k) c_k. For a real
 *     plan (real nu and c) and real z the two halves are conjugate, so one Horner
 *     pass gives ce = Re(e^{i theta z} P(w)) -- the hot path.
 *
 * Plans are cached (two slots, keyed by the bit patterns of a and q), so a buffer
 * of z at fixed (a, q) builds the plan once. The cache is thread-local under
 * MATHILDA_THREADS (same guard as the Compile VM) and its storage is reused,
 * never freed (still-reachable, bounded by two plans).
 */
#include "sf_machine_mathieu.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#if defined(MATHILDA_THREADS) && (defined(__GNUC__) || defined(__clang__))
#define MATHIEU_TLS __thread
#else
#define MATHIEU_TLS
#endif

/* ------------------------------------------------------------------ */
/* Minimal complex double arithmetic (C99, no <complex.h>)             */
/* ------------------------------------------------------------------ */

typedef struct { double re, im; } mcx;

static inline mcx mc(double re, double im) { mcx z; z.re = re; z.im = im; return z; }
static inline mcx mc_add(mcx a, mcx b) { return mc(a.re + b.re, a.im + b.im); }
static inline mcx mc_sub(mcx a, mcx b) { return mc(a.re - b.re, a.im - b.im); }
static inline mcx mc_mul(mcx a, mcx b) {
    return mc(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}
static inline mcx mc_scale(mcx a, double s) { return mc(a.re * s, a.im * s); }
static inline double mc_abs(mcx a) { return hypot(a.re, a.im); }
/* Smith's division: no spurious overflow for well-scaled operands. */
static inline mcx mc_div(mcx a, mcx b) {
    if (fabs(b.re) >= fabs(b.im)) {
        double r = b.im / b.re, d = b.re + b.im * r;
        return mc((a.re + a.im * r) / d, (a.im - a.re * r) / d);
    } else {
        double r = b.re / b.im, d = b.re * r + b.im;
        return mc((a.re * r + a.im) / d, (a.im * r - a.re) / d);
    }
}
/* Principal square root. */
static mcx mc_sqrt(mcx a) {
    if (a.im == 0.0) {
        if (a.re >= 0.0) return mc(sqrt(a.re), 0.0);
        return mc(0.0, sqrt(-a.re));
    }
    double m = mc_abs(a);
    double t = sqrt(0.5 * (m + fabs(a.re)));
    if (a.re >= 0.0) return mc(t, a.im / (2.0 * t));
    return mc(fabs(a.im) / (2.0 * t), a.im >= 0.0 ? t : -t);
}
static mcx mc_exp(mcx a) {
    double e = exp(a.re);
    return a.im == 0.0 ? mc(e, 0.0) : mc(e * cos(a.im), e * sin(a.im));
}
static mcx mc_log(mcx a) { return mc(log(mc_abs(a)), atan2(a.im, a.re)); }
/* sin(pi x)/(pi x), entire; a short series near 0. */
static mcx mc_sinc_pi(mcx x) {
    if (mc_abs(x) < 1e-4) {
        mcx y = mc_scale(mc_mul(x, x), M_PI * M_PI);            /* (pi x)^2 */
        mcx s = mc_add(mc(1.0, 0.0), mc_scale(y, -1.0 / 6.0));
        return mc_add(s, mc_scale(mc_mul(y, y), 1.0 / 120.0));
    }
    mcx px = mc_scale(x, M_PI);
    /* sin(u) = sin(ur) cosh(ui) + i cos(ur) sinh(ui) */
    mcx s = mc(sin(px.re) * cosh(px.im), cos(px.re) * sinh(px.im));
    return mc_div(s, px);
}
/* acos without the cancellation of the textbook formula for Re z > 0:
 * acos(z) = pi - acos(-z), and for Re z <= 0, acos(z) = -i log(z + i sqrt(1-z^2)). */
static mcx mc_acos(mcx z) {
    if (z.im == 0.0 && fabs(z.re) <= 1.0) return mc(acos(z.re), 0.0);
    int flip = z.re > 0.0;
    if (flip) z = mc(-z.re, -z.im);
    mcx s = mc_sqrt(mc_sub(mc(1.0, 0.0), mc_mul(z, z)));
    mcx w = mc_add(z, mc(-s.im, s.re));                          /* z + i s */
    mcx l = mc_log(w);
    mcx r = mc(l.im, -l.re);                                     /* -i log w */
    if (flip) r = mc(M_PI - r.re, -r.im);
    return r;
}
static inline int mc_finite(mcx a) { return isfinite(a.re) && isfinite(a.im); }

/* ------------------------------------------------------------------ */
/* Step 1: Hill's determinant -> cos(pi nu)                            */
/* ------------------------------------------------------------------ */

/* Trigamma psi_1(x) for x >= ~30 (asymptotic; the truncation K keeps x large). */
static double mth_psi1_large(double x) {
    double x2 = x * x, ix = 1.0 / x, ix2 = 1.0 / x2;
    (void)x2;
    return ix + 0.5 * ix2 + ix * ix2 * (1.0 / 6.0 - ix2 * (1.0 / 30.0 - ix2 / 42.0));
}

static bool mth_hill_cos(mcx a, mcx q, mcx* C_out) {
    double aa = mc_abs(a), aq = mc_abs(q);
    mcx x = mc_scale(mc_sqrt(a), 0.5);                           /* sqrt(a)/2 */
    if (x.re < 0.0) x = mc(-x.re, -x.im);                        /* the product is even */
    /* K: the first omitted tail order is ~ max(|a|^2 |q|^2, |q|^4) / K^7. */
    double big = fmax(aa * aa * aq * aq, aq * aq * aq * aq);
    double K0 = pow(big * 1e17, 1.0 / 7.0);
    double Kd = 40.0 + 2.0 * mc_abs(x) + 2.0 * sqrt(aq) + (isfinite(K0) ? K0 : 1e9);
    if (Kd > 40000.0) return false;
    int K = (int)Kd;

    mcx q2 = mc_mul(q, q);
    mcx Dn2 = mc(1.0, 0.0), Dn1 = mc(1.0, 0.0);                 /* D_{K+2}, D_{K+1} */
    int sc = 0;                                                  /* D scaled by 1e-150^sc */
    for (int k = K; k >= 1; k--) {
        double kk = (double)k * (double)k;
        mcx t = mc_sub(mc(1.0, 0.0), mc_scale(a, 1.0 / (4.0 * kk)));
        mcx e = mc_scale(q2, 1.0 / (16.0 * kk * ((double)k + 1.0) * ((double)k + 1.0)));
        mcx Dk = mc_sub(mc_mul(t, Dn1), mc_mul(e, Dn2));
        Dn2 = Dn1; Dn1 = Dk;
        /* Large |q|: the D_k grow like (q^2/16)^k/(k!)^4 before they settle;
         * rescale both (the recurrence is linear) and remember the exponent. */
        if (mc_abs(Dn1) > 1e150) {
            Dn1 = mc_scale(Dn1, 1e-150); Dn2 = mc_scale(Dn2, 1e-150); sc++;
        }
    }
    mcx D1 = Dn1, D2 = Dn2;
    mcx detT = mc_mul(D1, mc_add(mc_mul(a, D1), mc_scale(mc_mul(q2, D2), 0.5)));

    /* R = [sin(pi x)/(pi x)] / prod_{k=1..K} (1 - x^2/k^2), stable at x ~ j. */
    long j = lround(x.re);
    mcx x2 = mc_mul(x, x);
    mcx prod = mc(1.0, 0.0);
    for (int k = 1; k <= K; k++) {
        if ((long)k == j) continue;
        prod = mc_mul(prod, mc_sub(mc(1.0, 0.0), mc_scale(x2, 1.0 / ((double)k * (double)k))));
    }
    mcx R;
    if (j >= 1 && j <= K) {
        /* sin(pi x)/(1 - x^2/j^2) = (-1)^(j+1) pi sinc_pi(d) j^2/(x + j), d = x - j */
        mcx d = mc_sub(x, mc((double)j, 0.0));
        mcx f = mc_scale(mc_div(mc_sinc_pi(d), mc_add(x, mc((double)j, 0.0))),
                         M_PI * (double)j * (double)j * ((j & 1) ? 1.0 : -1.0));
        R = mc_div(f, mc_mul(mc_scale(x, M_PI), prod));
    } else {
        R = mc_div(mc_sinc_pi(x), prod);
    }

    /* Tail of the scaled determinant beyond K, both sides (hence the factor 2):
     *   log = -sum e_k - (a/4) sum e_k (1/k^2 + 1/(k+1)^2) - ...,
     *   sum_{k>K} 1/(k^2 (k+1)^2) = psi1(K+1) + psi1(K+2) - 2/(K+1),
     *   sum_{k>K} (1/(k^4 (k+1)^2) + 1/(k^2 (k+1)^4)) ~ 2/(5 (K+1)^5). */
    double Kp = (double)K + 1.0;
    double S1 = mth_psi1_large(Kp) + mth_psi1_large(Kp + 1.0) - 2.0 / Kp;
    double S2 = 2.0 / (5.0 * pow(Kp + 0.5, 5.0));
    mcx lt = mc_add(mc_scale(q2, -S1 / 16.0), mc_scale(mc_mul(a, q2), -S2 / 64.0));
    mcx tail = mc_exp(mc_scale(lt, 2.0));

    mcx s;
    if (sc == 0) {
        s = mc_scale(mc_mul(mc_mul(detT, mc_mul(R, R)), tail), M_PI * M_PI / 4.0);
    } else {
        /* det(T) carries 10^(300 sc): combine in log form (D_1^2 scaling). */
        mcx l = mc_add(mc_log(detT), mc_scale(mc_log(R), 2.0));
        l = mc_add(l, mc_scale(lt, 2.0));
        l.re += 300.0 * (double)sc * log(10.0) + log(M_PI * M_PI / 4.0);
        s = mc_exp(l);
    }
    mcx C = mc_sub(mc(1.0, 0.0), mc_scale(s, 2.0));
    if (!mc_finite(C)) return false;
    *C_out = C;
    return true;
}

/* ------------------------------------------------------------------ */
/* The plan                                                            */
/* ------------------------------------------------------------------ */

#define MATHIEU_MAXN 4096        /* coefficient count the machine plan accepts */

typedef struct {
    int valid;                   /* slot holds a plan for (key) */
    int ok;                      /* the plan was built (else: decline) */
    uint64_t key[4];             /* bit patterns of ar, ai, qr, qi */
    int q0;                      /* q == 0: closed form */
    mcx a, q, nu;
    mcx sa;                      /* sqrt(a) for the q == 0 form */
    int integer;                 /* nu snapped to an integer */
    int m;                       /* dominant coefficient index */
    int lo, n;                   /* coefficients c_{lo .. lo+n-1} (to ~1e-300) */
    int rlo, rhi;                /* the |c_k| >= 1e-20 max|c| window (real z) */
    int real;                    /* nu and every c_k real (real hot path) */
    int realpq;                  /* a and q real */
    int imag_kind;               /* realpq, but the normalized ce is imaginary */
    int imag_kind_s;             /* realpq, but the normalized se is imaginary */
    int parity;                  /* integer nu: +1 c_k = c_{-n-k} (a = a_n: ce
                                  * periodic, se = 0), -1 antisymmetric (a = b_n:
                                  * se periodic, ce = 0); 0 for generic nu */
    mcx y0, y1;                  /* ce(0) = sum c_k, se'(0) = sum (nu+2k) c_k: the
                                  * Taylor seeds for the small-|z| arm */
    int cap;                     /* allocated length of the arrays */
    mcx* c;                      /* normalized coefficients */
    mcx* w;                      /* (nu + 2k) c_k, for the derivative */
    mcx* scratch;                /* work space (ratios / LU) */
    mcx* scratch2;
    mcx* scratch3;
    mcx* scratch4;
    double* lc;                  /* log|c_k| (-HUGE_VAL for an exact zero) */
} MathieuPlan;

static MATHIEU_TLS MathieuPlan g_plans[2];
static MATHIEU_TLS int g_plan_next;

static int plan_reserve(MathieuPlan* p, int n) {
    if (p->cap >= n) return 1;
    int cap = n < 64 ? 64 : n;
    mcx* arr[6] = { p->c, p->w, p->scratch, p->scratch2, p->scratch3, p->scratch4 };
    for (int i = 0; i < 6; i++) {
        mcx* r = (mcx*)realloc(arr[i], sizeof(mcx) * (size_t)cap);
        if (!r) return 0;
        arr[i] = r;
    }
    p->c = arr[0]; p->w = arr[1]; p->scratch = arr[2];
    p->scratch2 = arr[3]; p->scratch3 = arr[4]; p->scratch4 = arr[5];
    double* l = (double*)realloc(p->lc, sizeof(double) * (size_t)cap);
    if (!l) return 0;
    p->lc = l;
    p->cap = cap;
    return 1;
}

static inline mcx diag_d(mcx a, mcx nu, int k) {      /* a - (nu + 2k)^2 */
    mcx t = mc_add(nu, mc(2.0 * (double)k, 0.0));
    return mc_sub(a, mc_mul(t, t));
}

/* Two steps of inverse iteration on the truncated Floquet matrix
 *   T = tridiag(-q, a - (nu+2k)^2, -q),   k = lo .. lo+n-1,
 * by Gaussian elimination with partial pivoting (banded, fill <= 2). Returns the
 * index of the largest |x_k|, or INT32_MIN on failure. */
static int dominant_index(MathieuPlan* p, mcx nu) {
    int n = p->n, lo = p->lo;
    mcx* x = p->scratch;     /* rhs / solution */
    mcx* d0 = p->scratch2;   /* row i: col i, i+1, i+2 after elimination */
    mcx* d1 = p->scratch3;
    mcx* d2 = p->scratch4;
    mcx mq = mc(-p->q.re, -p->q.im);
    double scale = mc_abs(p->a) + 2.0 * mc_abs(p->q) + 1.0;
    for (int i = 0; i < n; i++) x[i] = mc(1.0, 0.0);
    int best = 0;
    for (int it = 0; it < 2; it++) {
        /* (re)build the matrix: rows hold (diag, super, super2) */
        for (int i = 0; i < n; i++) {
            d0[i] = diag_d(p->a, nu, lo + i);
            d1[i] = (i + 1 < n) ? mq : mc(0.0, 0.0);
            d2[i] = mc(0.0, 0.0);
        }
        /* forward elimination; the sub-diagonal of row i+1 is -q (before swaps) */
        mcx sub = mq;
        for (int i = 0; i + 1 < n; i++) {
            /* row i: (d0[i], d1[i], d2[i]) at cols i, i+1, i+2
             * row i+1: (sub, d0[i+1], d1[i+1]) at cols i, i+1, i+2 */
            mcx r1a = sub, r1b = d0[i + 1], r1c = d1[i + 1];
            if (mc_abs(r1a) > mc_abs(d0[i])) {
                mcx ta = d0[i], tb = d1[i], tc = d2[i];
                d0[i] = r1a; d1[i] = r1b; d2[i] = r1c;
                r1a = ta; r1b = tb; r1c = tc;
                mcx tx = x[i]; x[i] = x[i + 1]; x[i + 1] = tx;
            }
            if (d0[i].re == 0.0 && d0[i].im == 0.0) d0[i] = mc(1e-300 * scale, 0.0);
            mcx f = mc_div(r1a, d0[i]);
            d0[i + 1] = mc_sub(r1b, mc_mul(f, d1[i]));
            d1[i + 1] = mc_sub(r1c, mc_mul(f, d2[i]));
            x[i + 1] = mc_sub(x[i + 1], mc_mul(f, x[i]));
            sub = mq;
        }
        if (d0[n - 1].re == 0.0 && d0[n - 1].im == 0.0) d0[n - 1] = mc(1e-300 * scale, 0.0);
        /* back substitution */
        for (int i = n - 1; i >= 0; i--) {
            mcx s = x[i];
            if (i + 1 < n) s = mc_sub(s, mc_mul(d1[i], x[i + 1]));
            if (i + 2 < n) s = mc_sub(s, mc_mul(d2[i], x[i + 2]));
            x[i] = mc_div(s, d0[i]);
        }
        double mx = 0.0; best = 0;
        for (int i = 0; i < n; i++) {
            double v = mc_abs(x[i]);
            if (!(v == v)) return INT32_MIN;
            if (v > mx) { mx = v; best = i; }
        }
        if (!(mx > 0.0) || !isfinite(mx)) return INT32_MIN;
        mcx inv = mc_div(mc(1.0, 0.0), x[best]);
        for (int i = 0; i < n; i++) x[i] = mc_mul(x[i], inv);
    }
    return lo + best;
}

/* g(nu) and g'(nu) at the gluing index m (see the file comment). */
static bool consistency_g(const MathieuPlan* p, mcx nu, int m, mcx* g, mcx* dg) {
    int lo = p->lo, hi = p->lo + p->n - 1;
    mcx q = p->q;
    mcx rp = mc(0.0, 0.0), drp = mc(0.0, 0.0);
    for (int k = hi; k > m; k--) {
        mcx t = mc_add(nu, mc(2.0 * (double)k, 0.0));
        mcx d = diag_d(p->a, nu, k);
        mcx dd = mc_scale(t, -2.0);
        mcx den = mc_sub(d, mc_mul(q, rp));
        mcx dden = mc_sub(dd, mc_mul(q, drp));
        if (den.re == 0.0 && den.im == 0.0) return false;
        rp = mc_div(q, den);
        drp = mc_mul(mc_scale(rp, -1.0), mc_div(dden, den));     /* -q dden/den^2 */
    }
    mcx rm = mc(0.0, 0.0), drm = mc(0.0, 0.0);
    for (int k = lo; k < m; k++) {
        mcx t = mc_add(nu, mc(2.0 * (double)k, 0.0));
        mcx d = diag_d(p->a, nu, k);
        mcx dd = mc_scale(t, -2.0);
        mcx den = mc_sub(d, mc_mul(q, rm));
        mcx dden = mc_sub(dd, mc_mul(q, drm));
        if (den.re == 0.0 && den.im == 0.0) return false;
        rm = mc_div(q, den);
        drm = mc_mul(mc_scale(rm, -1.0), mc_div(dden, den));
    }
    mcx tm = mc_add(nu, mc(2.0 * (double)m, 0.0));
    *g = mc_sub(diag_d(p->a, nu, m), mc_mul(q, mc_add(rp, rm)));
    *dg = mc_sub(mc_scale(tm, -2.0), mc_mul(q, mc_add(drp, drm)));
    return mc_finite(*g) && mc_finite(*dg);
}

/* Coefficients by the two continued fractions glued at m; c_m = 1 (unnormalized). */
static bool glue_coefficients(MathieuPlan* p) {
    int lo = p->lo, n = p->n, hi = lo + n - 1, m = p->m;
    mcx q = p->q, nu = p->nu;
    mcx* c = p->c;
    mcx* r = p->scratch;                     /* r[k-lo]: ratio c_k / c_{k-+1} */
    mcx rr = mc(0.0, 0.0);
    for (int k = hi; k > m; k--) {           /* r_k = c_k / c_{k-1} */
        mcx den = mc_sub(diag_d(p->a, nu, k), mc_mul(q, rr));
        if (den.re == 0.0 && den.im == 0.0) return false;
        rr = mc_div(q, den);
        r[k - lo] = rr;
    }
    rr = mc(0.0, 0.0);
    for (int k = lo; k < m; k++) {           /* r_k = c_k / c_{k+1} */
        mcx den = mc_sub(diag_d(p->a, nu, k), mc_mul(q, rr));
        if (den.re == 0.0 && den.im == 0.0) return false;
        rr = mc_div(q, den);
        r[k - lo] = rr;
    }
    c[m - lo] = mc(1.0, 0.0);
    for (int k = m + 1; k <= hi; k++) c[k - lo] = mc_mul(c[k - 1 - lo], r[k - lo]);
    for (int k = m - 1; k >= lo; k--) c[k - lo] = mc_mul(c[k + 1 - lo], r[k - lo]);
    for (int i = 0; i < n; i++) if (!mc_finite(c[i])) return false;
    return true;
}

/* Integer nu = n (a is a characteristic value). The decaying coefficient
 * vector is then symmetric (c_k = c_{-n-k}: a = a_n, ce periodic) or
 * antisymmetric (a = b_n, se periodic) under the fold k -> -n-k about the centre
 * -n/2, and in the antisymmetric case with n even the centre coefficient is
 * exactly 0 -- a continued fraction of ratios run THROUGH the centre divides by
 * it. So only the upper half is computed: the ratios r_k = c_k / c_{k-1} from
 * +infinity down to the centre, and a forward recurrence from the centre (where
 * the parity fixes the start) up to the glue point g (the dominant index, folded
 * onto the upper half). The parity is the one whose two halves agree at g --
 * the classical a_n / b_n consistency test -- and the lower half is the mirror,
 * exactly +-c, so the wrong-parity sum (se at a_n, ce at b_n) is structurally
 * zero. Sets p->parity; c is unnormalized. */
static bool glue_integer(MathieuPlan* p) {
    int n = (int)p->nu.re, lo = p->lo, hi = p->lo + p->n - 1;
    mcx q = p->q, a = p->a, nu = p->nu;
    mcx* c = p->c;
    mcx* r = p->scratch;                 /* r[k - lo] = c_k / c_{k-1}, k > kc */
    mcx* f = p->scratch2;                /* forward values on [kc, gmax] */
    int even = ((n % 2) + 2) % 2 == 0;
    int kc = even ? -n / 2 : -(n - 1) / 2;   /* nu + 2 kc = 0 (even) or 1 (odd) */
    int g = p->m;
    if (-n - g > g) g = -n - g;
    if (g < kc) g = kc;
    if (kc < lo || g + 2 > hi || -n - lo > hi) return false;
    /* r_{kc+1} = c_{kc+1}/c_kc is infinite in the antisymmetric even case (its
     * denominator vanishes, exactly or nearly); it is only read by the
     * symmetric test with g == kc, which then cannot hold */
    int r_inf = 0;
    mcx rr = mc(0.0, 0.0);
    for (int k = hi; k > kc; k--) {
        mcx den = mc_sub(diag_d(a, nu, k), mc_mul(q, rr));
        if (den.re == 0.0 && den.im == 0.0) {
            if (k == kc + 1 && even) { r_inf = 1; break; }
            return false;
        }
        rr = mc_div(q, den);
        r[k - lo] = rr;
    }
    double best = INFINITY;
    int bt = 0, bg = g;
    for (int t = 1; t >= -1; t -= 2) {
        int g0 = (even && t < 0 && g == kc) ? kc + 1 : g;   /* c_kc = 0 there */
        if (r_inf && g0 + 1 == kc + 1) continue;
        int gmax = g0 > kc + 1 ? g0 : kc + 1;
        /* the start: even & symmetric  d_kc c_kc = 2 q c_{kc+1};
         *            even & antisym    c_kc = 0, c_{kc+1} = 1;
         *            odd               c_{kc-1} = t c_kc */
        if (even && t < 0) { f[0] = mc(0.0, 0.0); f[1] = mc(1.0, 0.0); }
        else if (even) {
            f[0] = mc(1.0, 0.0);
            f[1] = mc_div(diag_d(a, nu, kc), mc_scale(q, 2.0));
        } else {
            f[0] = mc(1.0, 0.0);
            f[1] = mc_div(mc_sub(diag_d(a, nu, kc), mc_scale(q, (double)t)), q);
        }
        for (int k = kc + 1; k < gmax; k++) {             /* f_{k+1} */
            mcx v = mc_div(mc_sub(mc_mul(diag_d(a, nu, k), f[k - kc]),
                                  mc_mul(q, f[k - 1 - kc])), q);
            f[k + 1 - kc] = v;
            if (mc_abs(v) > 1e250)
                for (int j = 0; j <= k + 1 - kc; j++) f[j] = mc_scale(f[j], 1e-250);
        }
        /* consistency at g0: d c_g0 - q c_{g0-1} - q c_{g0+1} = 0, with
         * c_{g0+1} = r_{g0+1} c_g0 and c_{kc-1} the mirror t c_{-n-kc+1} */
        mcx fg = f[g0 - kc], fm;
        if (g0 - 1 >= kc) fm = f[g0 - 1 - kc];
        else fm = mc_scale(f[-n - (g0 - 1) - kc], (double)t);
        if (fg.re == 0.0 && fg.im == 0.0) continue;
        mcx res = mc_sub(mc_sub(diag_d(a, nu, g0), mc_mul(q, mc_div(fm, fg))),
                         mc_mul(q, r[g0 + 1 - lo]));
        double scale = mc_abs(diag_d(a, nu, g0)) + 2.0 * mc_abs(q) + 1.0;
        double rel = mc_abs(res) / scale;
        if (rel < best) { best = rel; bt = t; bg = g0; }
    }
    if (bt == 0 || !(best < 1e-6)) return false;
    /* rebuild the chosen forward half (the loop above may have overwritten f) */
    {
        int t = bt, g0 = bg;
        int gmax = g0 > kc + 1 ? g0 : kc + 1;
        if (even && t < 0) { f[0] = mc(0.0, 0.0); f[1] = mc(1.0, 0.0); }
        else if (even) {
            f[0] = mc(1.0, 0.0);
            f[1] = mc_div(diag_d(a, nu, kc), mc_scale(q, 2.0));
        } else {
            f[0] = mc(1.0, 0.0);
            f[1] = mc_div(mc_sub(diag_d(a, nu, kc), mc_scale(q, (double)t)), q);
        }
        for (int k = kc + 1; k < gmax; k++) {
            mcx v = mc_div(mc_sub(mc_mul(diag_d(a, nu, k), f[k - kc]),
                                  mc_mul(q, f[k - 1 - kc])), q);
            f[k + 1 - kc] = v;
            if (mc_abs(v) > 1e250)
                for (int j = 0; j <= k + 1 - kc; j++) f[j] = mc_scale(f[j], 1e-250);
        }
        for (int k = kc; k <= g0; k++) c[k - lo] = f[k - kc];
        for (int k = g0 + 1; k <= hi; k++) c[k - lo] = mc_mul(c[k - 1 - lo], r[k - lo]);
        for (int k = lo; k < kc; k++) {
            int kp = -n - k;                                   /* > k, <= hi */
            c[k - lo] = mc_scale(c[kp - lo], (double)t);
        }
        if (even && t < 0) c[kc - lo] = mc(0.0, 0.0);
        p->parity = t;
    }
    for (int i = 0; i < p->n; i++) if (!mc_finite(c[i])) return false;
    return true;
}

/* Choose the exponent representative nearest sqrt(a). */
static mcx pick_nu(mcx a, mcx C) {
    mcx nu0 = mc_scale(mc_acos(C), 1.0 / M_PI);
    mcx sa = mc_sqrt(a);
    mcx best = nu0; double bd = INFINITY;
    for (int sg = 1; sg >= -1; sg -= 2) {
        mcx v0 = mc_scale(nu0, (double)sg);
        double m0 = floor((sa.re - v0.re) / 2.0 + 0.5);
        for (int dm = -1; dm <= 1; dm++) {
            mcx v = mc_add(v0, mc(2.0 * (m0 + dm), 0.0));
            double d = mc_abs(mc_sub(v, sa));
            /* a tie (sqrt(a) on the imaginary axis: a < 0, nu real) takes
             * Re(nu) >= 0. MathieuC cannot tell nu from -nu, but MathieuS
             * changes sign with it: se -> sin(sqrt(a) z) as q -> 0 needs the
             * branch continuous from a > 0 (the DLMF se_n sign). */
            if (d < bd - 1e-12 || (fabs(d - bd) <= 1e-12 && v.re > best.re)) {
                bd = d; best = v;
            }
        }
    }
    return best;
}

static bool build_plan(MathieuPlan* p, mcx a, mcx q) {
    p->a = a; p->q = q;
    p->realpq = (a.im == 0.0 && q.im == 0.0);
    p->imag_kind = 0;
    p->q0 = (q.re == 0.0 && q.im == 0.0);
    if (p->q0) {
        p->sa = mc_sqrt(a);
        p->real = (a.im == 0.0 && a.re >= 0.0);
        return true;
    }
    mcx C;
    if (!mth_hill_cos(a, q, &C)) return false;
    mcx nu = pick_nu(a, C);
    double nr = floor(nu.re + 0.5);
    double eps = mc_abs(mc_sub(nu, mc(nr, 0.0)));
    p->integer = eps < 2e-6 * fmax(1.0, fabs(nr));
    /* the inverse iteration below runs at the UNSNAPPED estimate: at an exact
     * characteristic value the Floquet matrix at the integer nu is singular */
    mcx nu_h = nu;
    if (p->integer && eps < 1e-9) nu_h = mc_add(mc(nr, 0.0), mc(1e-7, 0.0));
    if (p->integer) nu = mc(nr, 0.0);

    /* Coefficient range: the oscillatory zone |nu+2k|^2 <~ |a| + 2|q| plus the
     * decay tail, until the ratio product falls below 1e-300: a complex z
     * weights c_k by e^{2|k||Im z|}, so the real-z cutoff (1e-20) is not enough
     * there; the evaluator trims per point (see eval_plan). */
    double aa = mc_abs(a), aq = mc_abs(q);
    int W = (int)(sqrt(aa + 2.0 * aq) / 2.0) + 2;
    double r = 1.0; int Kx = 0;
    while (r > 1e-300 && Kx < MATHIEU_MAXN) {
        double k = (double)(W + 1 + Kx);
        double den = 4.0 * k * k - aa - 2.0 * aq;
        r *= den > aq ? aq / den : 1.0;
        Kx++;
    }
    Kx += 2;
    int kc = (int)floor(-nu.re / 2.0 + 0.5);
    int half = W + Kx;
    if (2 * half + 1 > MATHIEU_MAXN) return false;
    p->lo = kc - half;
    p->n = 2 * half + 1;
    if (!plan_reserve(p, p->n)) return false;

    int m = dominant_index(p, p->integer ? nu_h : nu);
    if (m == INT32_MIN) return false;
    p->m = m;

    if (!p->integer) {
        /* Newton on g. Converges quadratically from the Hill estimate; near a band
         * edge (eps small) the two equivalent roots nu, 2n - nu are close, which
         * slows it, hence the generous iteration cap. */
        bool conv = false;
        double prev = INFINITY;
        for (int it = 0; it < 40; it++) {
            mcx g, dg;
            if (!consistency_g(p, nu, m, &g, &dg)) break;
            if (dg.re == 0.0 && dg.im == 0.0) break;
            mcx step = mc_div(g, dg);
            nu = mc_sub(nu, step);
            if (!mc_finite(nu)) break;
            double sa = mc_abs(step), sc = 1.0 + mc_abs(nu);
            if (sa <= 4e-16 * sc) { conv = true; break; }
            /* The noise floor: near a band edge g' is small (the equivalent
             * roots nu and 2n - nu nearly coincide), so the rounding of g
             * moves nu by ~eps/|g'| -- nu itself is that ill-conditioned in a.
             * Once the steps are tiny and no longer shrinking quadratically,
             * the iteration has converged as far as double allows. */
            if (sa < 1e-9 * sc && sa > 0.25 * prev) { conv = true; break; }
            prev = sa;
        }
        if (!conv) return false;
        /* Real inputs in a stability band: nu is real; drop the rounding residue
         * so the real hot path applies (|C| <= 1 <=> nu real). */
        if (p->realpq && C.im == 0.0 && fabs(C.re) <= 1.0) nu.im = 0.0;
        /* Newton may have moved to an equivalent root (+-nu + 2j); re-centre on
         * the representative nearest sqrt(a), keeping the refined digits. */
        mcx Cn = mc(cos(M_PI * nu.re) * cosh(M_PI * nu.im),
                    -sin(M_PI * nu.re) * sinh(M_PI * nu.im));
        mcx nu2 = pick_nu(a, Cn);
        mcx cand1 = mc_add(nu, mc(2.0 * floor((nu2.re - nu.re) / 2.0 + 0.5), 0.0));
        mcx cand2 = mc_add(mc_scale(nu, -1.0),
                           mc(2.0 * floor((nu2.re + nu.re) / 2.0 + 0.5), 0.0));
        mcx nn = mc_abs(mc_sub(cand1, nu2)) <= mc_abs(mc_sub(cand2, nu2)) ? cand1 : cand2;
        if (mc_abs(mc_sub(nn, nu)) > 1e-9) {
            nu = nn;
            p->lo = (int)floor(-nu.re / 2.0 + 0.5) - half;
            m = dominant_index(p, nu);
            if (m == INT32_MIN) return false;
            p->m = m;
        }
    }
    p->nu = nu;
    p->parity = 0;
    if (p->integer) {
        if (!glue_integer(p)) return false;
    } else if (!glue_coefficients(p)) return false;

    /* normalize: sum c_k^2 = 1, sign Re(c_0) > 0 */
    mcx S = mc(0.0, 0.0);
    double S1 = 0.0;
    for (int i = 0; i < p->n; i++) {
        S = mc_add(S, mc_mul(p->c[i], p->c[i]));
        S1 += p->c[i].re * p->c[i].re + p->c[i].im * p->c[i].im;
    }
    if (!(mc_abs(S) > 1e-12 * S1)) return false;                 /* isotropic */
    mcx s = mc_sqrt(S);
    int i0 = -p->lo;
    if (i0 < 0 || i0 >= p->n) i0 = p->m - p->lo;
    mcx c0 = mc_div(p->c[i0], s);
    if (c0.re < 0.0 || (fabs(c0.re) <= 1e-12 * fabs(c0.im) && c0.im < 0.0))
        s = mc(-s.re, -s.im);
    mcx inv = mc_div(mc(p->integer ? sqrt(0.5) : 1.0, 0.0), s);
    p->real = (nu.im == 0.0);
    for (int i = 0; i < p->n; i++) {
        mcx v = mc_mul(p->c[i], inv);
        p->c[i] = v;
        if (v.im != 0.0) p->real = 0;
        p->w[i] = mc_mul(mc_add(nu, mc(2.0 * (double)(p->lo + i), 0.0)), v);
    }
    double cmax = 0.0;
    for (int i = 0; i < p->n; i++) {
        double av = mc_abs(p->c[i]);
        p->lc[i] = av > 0.0 ? log(av) : -HUGE_VAL;
        if (av > cmax) cmax = av;
    }
    p->rlo = 0; p->rhi = p->n - 1;
    while (p->rlo < -p->lo - 1 && mc_abs(p->c[p->rlo]) < 1e-20 * cmax) p->rlo++;
    while (p->rhi > -p->lo && mc_abs(p->c[p->rhi]) < 1e-20 * cmax) p->rhi--;
    p->y0 = mc(0.0, 0.0); p->y1 = mc(0.0, 0.0);
    for (int i = 0; i < p->n; i++) {
        p->y0 = mc_add(p->y0, p->c[i]);
        p->y1 = mc_add(p->y1, p->w[i]);
    }
    if (p->realpq && !p->real) {
        /* Real (a, q): ce = N_C y_e(z), se = N_S y_o(z) with y_e, y_o real and
         * N_C^2, N_S^2 real (conj(c_k) = +-c_{-n-k} for nu = n + i mu), so each
         * value is real or purely imaginary. Read which from N_C = ce(0) and
         * N_S = se'(0). */
        p->imag_kind = fabs(p->y0.im) > fabs(p->y0.re);
        p->imag_kind_s = fabs(p->y1.im) > fabs(p->y1.re);
    } else {
        p->imag_kind_s = 0;
    }
    return true;
}

static void key_of(double ar, double ai, double qr, double qi, uint64_t key[4]) {
    double v[4] = { ar, ai, qr, qi };
    for (int i = 0; i < 4; i++) {
        double d = v[i] == 0.0 ? 0.0 : v[i];                       /* -0.0 == 0.0 */
        memcpy(&key[i], &d, sizeof(double));
    }
}

static MathieuPlan* get_plan(double ar, double ai, double qr, double qi) {
    uint64_t key[4];
    key_of(ar, ai, qr, qi, key);
    for (int i = 0; i < 2; i++) {
        MathieuPlan* p = &g_plans[i];
        if (p->valid && memcmp(p->key, key, sizeof key) == 0) return p->ok ? p : NULL;
    }
    MathieuPlan* p = &g_plans[g_plan_next];
    g_plan_next ^= 1;
    p->valid = 1;
    memcpy(p->key, key, sizeof key);
    p->ok = build_plan(p, mc(ar, ai), mc(qr, qi)) ? 1 : 0;
    return p->ok ? p : NULL;
}

/* ------------------------------------------------------------------ */
/* Evaluation                                                          */
/* ------------------------------------------------------------------ */

/* Small |z|: the Taylor series of the solution itself,
 *     y'' = (2 q cos 2z - a) y,   y(0) = ce(0) (even) or 0, y'(0) = se'(0) (odd),
 * f(z) = 2 q cos 2z - a = sum_j f_j z^j, (k+2)(k+1) y_{k+2} = sum_j f_j y_{k-j}.
 * Used for the two ODD functions (se, ce'), whose Floquet sums lose relative
 * accuracy as z -> 0 (value ~ z, terms ~ 1), and only for |z| <= tau, where the
 * series converges fast. Returns false when it has not converged (the Floquet
 * sum then answers). */
#define MATHIEU_TAYLOR_N 48
static bool taylor_eval(const MathieuPlan* p, int odd, int deriv, mcx z, mcx* out) {
    mcx f[MATHIEU_TAYLOR_N], y[MATHIEU_TAYLOR_N + 2];
    /* f_j: 2q (-4)^(j/2) / (j/2)! ... for even j; f_0 = 2q - a */
    double fac = 1.0;
    for (int j = 0; j < MATHIEU_TAYLOR_N; j++) {
        if (j & 1) { f[j] = mc(0.0, 0.0); continue; }
        int h = j / 2;
        if (h > 0) fac *= -4.0 / ((double)(j - 1) * (double)j);   /* (-4)^h/(2h)! */
        f[j] = mc_scale(p->q, 2.0 * fac);
    }
    f[0] = mc_sub(f[0], p->a);
    y[0] = odd ? mc(0.0, 0.0) : p->y0;
    y[1] = odd ? p->y1 : mc(0.0, 0.0);
    /* sum in powers of z (and the derivative), stop when two consecutive
     * nonzero terms are below 2^-60 of the running sum */
    mcx zp = mc(1.0, 0.0);                       /* z^k */
    mcx zpm = mc(0.0, 0.0);                      /* z^(k-1) */
    mcx sum = mc(0.0, 0.0);
    int small = 0;
    for (int k = 0; k < MATHIEU_TAYLOR_N; k++) {
        if (k >= 2) {
            mcx s = mc(0.0, 0.0);
            for (int j = 0; j <= k - 2; j += 2) s = mc_add(s, mc_mul(f[j], y[k - 2 - j]));
            y[k] = mc_scale(s, 1.0 / ((double)k * (double)(k - 1)));
        }
        mcx term = deriv ? (k >= 1 ? mc_scale(mc_mul(y[k], zpm), (double)k) : mc(0.0, 0.0))
                         : mc_mul(y[k], zp);
        sum = mc_add(sum, term);
        zpm = zp;
        zp = mc_mul(zp, z);
        if ((y[k].re != 0.0 || y[k].im != 0.0) && k >= 3) {
            if (mc_abs(term) <= 1e-18 * mc_abs(sum)) { if (++small >= 2) { *out = sum; return mc_finite(sum); } }
            else small = 0;
        }
    }
    return false;
}

static bool eval_plan(const MathieuPlan* p, int odd, int deriv, mcx z, mcx* out) {
    if (p->q0) {
        if (p->a.re == 0.0 && p->a.im == 0.0) {
            /* ce_0 = 1/sqrt 2; se = sin(0 z) = 0 */
            *out = (deriv || odd) ? mc(0.0, 0.0) : mc(sqrt(0.5), 0.0);
            return true;
        }
        mcx u = mc_mul(p->sa, z);                                 /* sqrt(a) z */
        mcx cu = mc(cos(u.re) * cosh(u.im), -sin(u.re) * sinh(u.im));
        mcx su = mc(sin(u.re) * cosh(u.im), cos(u.re) * sinh(u.im));
        if (odd) *out = deriv ? mc_mul(p->sa, cu) : su;
        else     *out = deriv ? mc_scale(mc_mul(p->sa, su), -1.0) : cu;
        return mc_finite(*out);
    }
    /* A characteristic value of the other parity (se at a = a_n, ce at b_n):
     * the Floquet-normalized function is identically zero (the limit from both
     * sides; no periodic solution of that parity exists). */
    if ((odd && p->parity > 0) || (!odd && p->parity < 0)) {
        *out = mc(0.0, 0.0);
        return true;
    }
    /* The odd functions (se, ce') vanish at 0 like z: near there the Floquet sum
     * keeps only absolute accuracy, so the Taylor series answers. */
    if ((odd ^ deriv) && z.re == 0.0 && z.im == 0.0) {      /* exactly 0 */
        *out = mc(0.0, 0.0);
        return true;
    }
    if (odd ^ deriv) {
        double tau = 0.25 / sqrt(1.0 + mc_abs(p->a) + 2.0 * mc_abs(p->q));
        if (mc_abs(z) <= tau && taylor_eval(p, odd, deriv, z, out)) return true;
    }
    /* Centred sums: c_k for k >= 0 (P) and k < 0 (M), so the phases are
     * e^{+-i nu z} with nu ~ sqrt(a) rather than nu + 2 lo -- a large offset
     * multiplies the rounding of z into the result. With w = e^{2iz}:
     *   ce(z)  = ( e^{i nu z}[P(w) + M(1/w)] + e^{-i nu z}[P(1/w) + M(w)] ) / 2,
     *   ce'(z) = (i/2)( e^{i nu z}[Pw(w) + Mw(1/w)] - e^{-i nu z}[Pw(1/w) + Mw(w)] ),
     * P(w) = sum_{k>=0} c_k w^k, M(u) = sum_{k>=1} c_{-k} u^k. */
    int i0 = -p->lo;                               /* index of k = 0; lo < 0 < hi */
    const mcx* c = deriv ? p->w : p->c;
    if (z.im == 0.0 && p->real) {
        int top = p->rhi, bot = p->rlo, nm = i0 - bot;
        /* Real hot path: the two halves are conjugate, so
         *   ce = Re(e^{i nu z} [P(w) + M(conj w)]),  ce' = -Im(e^{i nu z} [Pw + Mw]). */
        double cz = cos(2.0 * z.re), sz = sin(2.0 * z.re);
        double Pr = c[top].re, Pi = 0.0;
        for (int j = top - 1; j >= i0; j--) {      /* P(w), w = cz + i sz */
            double t = Pr * cz - Pi * sz;
            Pi = Pr * sz + Pi * cz;
            Pr = t + c[j].re;
        }
        double Mr = 0.0, Mi = 0.0;
        if (nm > 0) {                              /* M(conj w) = sum c_{-k} conj(w)^k */
            Mr = c[bot].re;
            for (int j = bot + 1; j < i0; j++) {   /* Horner, then one more conj(w) */
                double t = Mr * cz + Mi * sz;
                Mi = -Mr * sz + Mi * cz;
                Mr = t + c[j].re;
            }
            double t = Mr * cz + Mi * sz;
            Mi = -Mr * sz + Mi * cz;
            Mr = t;
        }
        double ph = p->nu.re * z.re;
        double er = cos(ph), ei = sin(ph);
        double Sr = Pr + Mr, Si = Pi + Mi;
        double vr = er * Sr - ei * Si, vi = er * Si + ei * Sr;
        /* A = e^{i nu z}[P + M] = sum c_k e^{i (nu+2k) z}: Re A = ce, Im A = se;
         * with the weights, Re A_w = se', -Im A_w = ce'. */
        if (odd) *out = deriv ? mc(vr, 0.0) : mc(vi, 0.0);
        else     *out = deriv ? mc(-vi, 0.0) : mc(vr, 0.0);
        return isfinite(out->re);
    }
    /* General complex path. Term k carries |c_k| e^{-2k Im z} (in w = e^{2iz})
     * or e^{+2k Im z} (in 1/w), times |e^{+-i nu z}|: trim each end to the terms
     * within e^-41 (1e-18) of the largest, and keep the largest log-magnitude for
     * the cancellation / overflow guard below. */
    int n = p->n;
    double y = z.im, enu = fabs(mc_mul(p->nu, z).im);
    /* Scan outward from k = 0: log|c_k| falls super-linearly, so once a term is
     * 41.5 below the running maximum AND outside the real window it only falls
     * further -- stop there. */
    double lmax = p->lc[i0];
    int top = i0, bot = i0;
    for (int i = i0 + 1; i < n; i++) {
        double t = p->lc[i] + 2.0 * fabs((double)(p->lo + i) * y);
        if (t > lmax) lmax = t;
        if (t >= lmax - 41.5) top = i;
        else if (i > p->rhi) break;
    }
    for (int i = i0 - 1; i >= 0; i--) {
        double t = p->lc[i] + 2.0 * fabs((double)(p->lo + i) * y);
        if (t > lmax) lmax = t;
        if (t >= lmax - 41.5) bot = i;
        else if (i < p->rlo) break;
    }
    if (!(lmax + enu < 650.0)) return false;                     /* overflow */
    /* The window was sized for |c_k| down to 1e-300; a large |Im z| weights the
     * tail by e^{2|k| |Im z|}, so a term at the window's END can still be within
     * e^-41.5 of the largest -- the sum would then be truncated, not converged.
     * Decline (the ball kernel widens its window for complex z). */
    if ((top == n - 1 && p->lc[n - 1] + 2.0 * fabs((double)(p->lo + n - 1) * y) >= lmax - 41.5) ||
        (bot == 0 && i0 > 0 && p->lc[0] + 2.0 * fabs((double)p->lo * y) >= lmax - 41.5))
        return false;
    if (bot >= i0) bot = i0;                                      /* no k < 0 terms */
    if (deriv) {
        /* the weights (nu+2k) only grow the terms polynomially */
        double wb = log(mc_abs(p->nu) + 2.0 * (double)(n + 1));
        lmax += wb;
    }
    mcx iz2 = mc(-2.0 * z.im, 2.0 * z.re);                       /* 2 i z */
    mcx w = mc_exp(iz2), wi = mc_exp(mc(-iz2.re, -iz2.im));
    mcx Pw = c[top], Pwi = c[top];
    for (int j = top - 1; j >= i0; j--) {
        Pw = mc_add(mc_mul(Pw, w), c[j]);
        Pwi = mc_add(mc_mul(Pwi, wi), c[j]);
    }
    mcx Mw = mc(0.0, 0.0), Mwi = mc(0.0, 0.0);
    if (bot < i0) {                                              /* k = bot-i0 .. -1 */
        Mw = c[bot]; Mwi = c[bot];
        for (int j = bot + 1; j < i0; j++) {
            Mw = mc_add(mc_mul(Mw, w), c[j]);
            Mwi = mc_add(mc_mul(Mwi, wi), c[j]);
        }
        Mw = mc_mul(Mw, w); Mwi = mc_mul(Mwi, wi);
    }
    mcx inz = mc_mul(mc(0.0, 1.0), mc_mul(p->nu, z));          /* i nu z */
    mcx e1 = mc_exp(inz), e2 = mc_exp(mc(-inz.re, -inz.im));
    mcx A = mc_mul(e1, mc_add(Pw, Mwi)), B = mc_mul(e2, mc_add(Pwi, Mw));
    /* A = F(z), B = F(-z), F = sum c_k e^{i (nu+2k) z}:
     *   ce = (A + B)/2, ce' = (i/2)(A_w - B_w), se = -(i/2)(A - B), se' = (A_w + B_w)/2 */
    mcx v;
    if (odd) v = deriv ? mc_scale(mc_add(A, B), 0.5) : mc_mul(mc(0.0, -0.5), mc_sub(A, B));
    else     v = deriv ? mc_mul(mc(0.0, 0.5), mc_sub(A, B)) : mc_scale(mc_add(A, B), 0.5);
    if (!mc_finite(v)) return false;
    /* Cancellation guard: the largest term vs the result. Losing more than four
     * digits of the sixteen is declined -- the arbitrary-precision kernel (ball
     * arithmetic, precision ladder) then answers. MathieuC[-1, 2, 5. I] is the
     * example: terms ~e^40, value ~1e-3. */
    double av = mc_abs(v);
    if (!(av > 0.0) || lmax + enu - log(av) > 9.3) {
        if (!(av > 0.0) && lmax + enu < -700.0) { *out = v; return true; }
        return false;
    }
    *out = v;
    return true;
}

bool sf_machine_mathieu_fn(int odd, int deriv, double ar, double ai, double qr,
                           double qi, double zr, double zi, double* rr, double* ri) {
    if (!isfinite(ar) || !isfinite(ai) || !isfinite(qr) || !isfinite(qi) ||
        !isfinite(zr) || !isfinite(zi)) return false;
    MathieuPlan* p = get_plan(ar, ai, qr, qi);
    if (!p) return false;
    mcx v;
    odd = odd ? 1 : 0;
    if (!eval_plan(p, odd, deriv ? 1 : 0, mc(zr, zi), &v)) return false;
    if (p->realpq && zi != 0.0 && zr == 0.0 && !p->q0) {
        /* Imaginary axis: the equation is real along it too, so the even
         * functions (ce, se') are N times a real there and the odd ones
         * (ce', se) N times i times a real. */
        int ik = odd ? p->imag_kind_s : p->imag_kind;
        int nreal = p->real || !ik;
        int vreal = (odd ^ (deriv ? 1 : 0)) ? !nreal : nreal;
        if (vreal) v.im = 0.0; else v.re = 0.0;
    }
    if (p->realpq && zi == 0.0) {
        /* real or purely imaginary, exactly */
        int ik = odd ? p->imag_kind_s : p->imag_kind;
        if (p->q0) {
            /* cos(sqrt(a) z), sqrt(a) sin: real; sin(sqrt(a) z), sqrt(a) cos(...)
             * for a < 0: i sinh, i cosh -- purely imaginary */
            if (odd && p->a.re < 0.0) v.re = 0.0; else v.im = 0.0;
        }
        else if (p->real) v.im = 0.0;
        else if (ik) v.re = 0.0;
        else v.im = 0.0;
    }
    *rr = v.re; *ri = v.im;
    return true;
}

bool sf_machine_mathieu(int deriv, double ar, double ai, double qr, double qi,
                        double zr, double zi, double* rr, double* ri) {
    return sf_machine_mathieu_fn(0, deriv, ar, ai, qr, qi, zr, zi, rr, ri);
}

bool sf_machine_mathieu_plan_info(double ar, double ai, double qr, double qi,
                                  double* nur, double* nui, int* integer,
                                  int* parity, int* m, int* lo, int* hi) {
    if (!isfinite(ar) || !isfinite(ai) || !isfinite(qr) || !isfinite(qi)) return false;
    MathieuPlan* p = get_plan(ar, ai, qr, qi);
    if (!p || p->q0) return false;
    *nur = p->nu.re; *nui = p->nu.im; *integer = p->integer; *parity = p->parity;
    *m = p->m; *lo = p->lo; *hi = p->lo + p->n - 1;
    return true;
}
