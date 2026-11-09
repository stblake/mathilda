/* Mathilda -- the Jacobi elliptic functions, in the PARAMETER convention m = k^2
 * (as EllipticK / EllipticF):
 *
 *   JacobiSN, JacobiCN, JacobiDN            sn, cn, dn (u | m)
 *   JacobiCD ... JacobiSD                   the nine ratios pq = p/q (n = 1)
 *   JacobiAmplitude[u, m]                   am(u | m): sin am = sn, cos am = cn,
 *                                           continuous in u (am(u + 4K) =
 *                                           am(u) + 2 Pi), the reference branch
 *   InverseJacobiSN[v, m] = EllipticF[ArcSin[v], m]
 *   InverseJacobiCN[v, m] = EllipticF[ArcCos[v], m]
 *   InverseJacobiDN[v, m] = EllipticF[ArcSin[Sqrt[(1 - v^2)/m]], m]
 *
 * Exact rules (the reference system's): u = 0, m = 0 (circular functions), m = 1
 * (hyperbolic), parity (sn, and every ratio with one s, is odd); am(u|0) = u,
 * am(u|1) = 2 ArcTan[Tanh[u/2]], am(0|m) = 0; the inverses at v = 0, 1 and at
 * m = 0, 1. D in u (DLMF 22.13), inert in m.
 *
 * Numerics. Machine kernel (sf_machine_jacobi), shared by the scalar call, the
 * NDArray / packed map and Compile[]:
 *   - real u, 0 <= m < 1: u reduced modulo 4K (the residual formed in
 *     double-double), then the arithmetic-geometric-mean method (DLMF 22.20):
 *     a_0 = 1, b_0 = sqrt(1 - m), phi_N = 2^N a_N u, phi_{n-1} = (phi_n +
 *     asin(c_n sin(phi_n) / a_n)) / 2; sn = sin phi_0, cn = cos phi_0, dn =
 *     sqrt(1 - m sn^2) or cos phi_0 / cos(phi_1 - phi_0) where that cancels;
 *   - complex u: the addition formulas with (x | m) and (y | 1 - m) (DLMF 22.8);
 *   - real m > 1: the reciprocal-parameter transformation (DLMF 22.17.4);
 *     real m < 0: the imaginary-modulus transformation (DLMF 22.17.6);
 *   - complex m, and every declined point: arb at 53 bits (flint_jacobi_d), so a
 *     buffer stays on its fast path;
 *   - the inverses through EllipticF's machine Carlson kernel where the value is
 *     real, the interpreter's EllipticF[ArcSin[...]] form (inheriting its
 *     branches) elsewhere.
 * Arbitrary precision: arb, through the theta quotients sn = (th3/th2)
 * th1(z)/th4(z) etc. at tau = i K(1-m)/K(m), z = u/(2K), after the same real-m
 * transformations (flint_num_jacobi).
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "jacobi_elliptic.h"

#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "flint_num_bridge.h"
#include "elliptic.h"          /* elliptic_machine_k, elliptic_machine_f */
#include "numeric.h"           /* expr_is_numeric_quantity */
#include "sf_scalar.h"
#include "sf_dd.h"
#include "symtab.h"
#include "../sym_names.h"

static const char* const JAC_NAME[JAC_COUNT] = {
    "JacobiSN", "JacobiCN", "JacobiDN", "JacobiCD", "JacobiCS", "JacobiDC",
    "JacobiDS", "JacobiNC", "JacobiND", "JacobiNS", "JacobiSC", "JacobiSD",
    "JacobiAmplitude", "InverseJacobiSN", "InverseJacobiCN", "InverseJacobiDN"
};
/* numerator / denominator letter of each ratio: s = 0, c = 1, d = 2, n = 3 */
static const int JAC_P[12] = { 0, 1, 2, 1, 1, 2, 2, 3, 3, 3, 0, 0 };
static const int JAC_Q[12] = { 3, 3, 3, 2, 0, 1, 0, 1, 2, 0, 1, 2 };

/* ------------------------------------------------------------------ */
/* Machine kernel                                                      */
/* ------------------------------------------------------------------ */

/* Double-double square root (one Newton step on the double root). */
static sdd jac_sqrt(sdd a) {
    if (a.hi <= 0.0) return sdd_from(a.hi <= 0.0 ? 0.0 : sqrt(a.hi));
    double s = sqrt(a.hi);
    sdd r = sdd_sub(a, sdd_two_prod(s, s));
    return sdd_quick(s, r.hi / (2.0 * s));
}

/* The AGM angle phi_0 for argument x from the tabulated a_n, c_n (DLMF 22.20.1):
 * phi_N = 2^N a_N x, phi_{n-1} = (phi_n + asin(c_n sin(phi_n)/a_n))/2; then
 * sn = sin phi_0, cn = cos phi_0, dn = sqrt(1 - m sn^2), or cos phi_0 /
 * cos(phi_1 - phi_0) where the square root would cancel. */
static void jac_agm_phi(const double* a, const double* cc, int N, double m, double x,
                        double* s, double* c, double* d, double* phi_out) {
    double phi = ldexp(a[N] * x, N), phi1 = phi;
    for (int n = N; n >= 1; n--) {
        phi1 = phi;
        phi = 0.5 * (phi + asin(cc[n] / a[n] * sin(phi)));
    }
    *s = sin(phi);
    *c = cos(phi);
    double ms2 = m * (*s) * (*s);
    *d = (ms2 <= 0.5 || N == 0) ? sqrt(1.0 - ms2) : (*c) / cos(phi1 - phi);
    *phi_out = phi;
}
/* sn, cn, dn and the amplitude for real u = uh + ul and 0 <= m = mh + ml <= 1.
 *
 * The argument and the parameter are carried as double-doubles so that neither
 * the transformations upstream (sqrt(m) u, 1/m, -m/(1-m)) nor the reduction
 * u mod 4K costs digits: the period 4K = 2 Pi / AGM(1, sqrt(1-m)) is computed in
 * double-double, so the residual of u = u' + 4K j is accurate to ~1 ulp of u'
 * for every |u| (a double K would leave j ulp of K behind -- 500 ulp at
 * |u| ~ 40 near a pole of the ratios, measured). */
static bool jac_real(double uh, double ul, double mh, double ml,
                     double* s, double* c, double* d, double* am) {
    if (!isfinite(uh) || !isfinite(mh)) return false;
    if (mh == 0.0 && ml == 0.0) {
        *s = sin(uh); *c = cos(uh); *d = 1.0; *am = uh;
        if (ul != 0.0) { *s += *c * ul; *c -= *s * ul; *am += ul; }
        return true;
    }
    if (mh == 1.0 && ml == 0.0) {
        double u = uh + ul;
        *s = tanh(u); *c = *d = 1.0 / cosh(u);
        *am = 2.0 * atan(tanh(0.5 * u));
        return true;
    }
    /* AGM in double-double: a_0 = 1, b_0 = sqrt(1 - m), c_n = (a_{n-1} - b_{n-1})/2 */
    sdd A = sdd_from(1.0);
    sdd Bm = sdd_add(sdd_from(1.0), sdd_neg(sdd_quick(mh, ml)));
    sdd B = jac_sqrt(Bm), B0 = B;
    double a[16], cc[16];
    int N = 0;
    a[0] = 1.0; cc[0] = sqrt(mh + ml);
    while (N < 15) {
        sdd cn = sdd_ldexp(sdd_sub(A, B), -1);
        if (fabs(cn.hi) <= 1e-17 * A.hi) break;
        sdd A1 = sdd_ldexp(sdd_add(A, B), -1);
        B = jac_sqrt(sdd_mul(A, B));
        A = A1;
        N++;
        a[N] = A.hi; cc[N] = cn.hi;
    }
    /* 4K = 2 Pi / a_N */
    sdd PI2; PI2.hi = 2.0 * SDD_PI_HI; PI2.lo = 2.0 * SDD_PI_LO;
    sdd K4 = sdd_div(PI2, A);
    sdd U = sdd_quick(uh, ul);
    double j = 0.0;
    if (fabs(uh) > 0.5 * K4.hi) {
        j = floor(uh / K4.hi + 0.5);
        if (fabs(j) > 1e15) return false;
        U = sdd_sub(U, sdd_mul_d(K4, j));            /* |U| <= 2K */
    }
    /* the nearest quarter period: U = r + q K, |r| <= K/2, so the values near a
     * zero of sn or cn come from a small argument (sd(r), sn(r)) with full
     * relative accuracy instead of from cos(phi) near pi/2 */
    sdd K = sdd_ldexp(K4, -2);
    double q = floor(U.hi / K.hi + 0.5);
    sdd Rr = sdd_sub(U, sdd_mul_d(K, q));
    double r = Rr.hi + Rr.lo;
    double s0, c0, d0, ph0;
    jac_agm_phi(a, cc, N, mh + ml, r, &s0, &c0, &d0, &ph0);
    int qi = (int)q;
    double kp = B0.hi + B0.lo;                       /* k' = sqrt(1 - m) */
    switch (qi) {
    case 0:  *s = s0;  *c = c0;  *d = d0; break;
    case 1:  *s = c0 / d0;  *c = -kp * s0 / d0; *d = kp / d0; break;
    case -1: *s = -c0 / d0; *c = kp * s0 / d0;  *d = kp / d0; break;
    default: *s = -s0; *c = -c0; *d = d0; break;    /* q = +-2 */
    }
    /* the amplitude, continuous in u: the AGM angle of the unshifted U */
    double sa, ca, da, pha;
    jac_agm_phi(a, cc, N, mh + ml, U.hi + U.lo, &sa, &ca, &da, &pha);
    *am = pha + 2.0 * M_PI * j;
    (void)ph0;
    return isfinite(*s) && isfinite(*c) && isfinite(*d);
}

/* Complex u = x + i y, 0 <= m <= 1: the addition formulas. */
typedef struct { double r, i; } jcx;
static jcx jc(double r, double i) { jcx z; z.r = r; z.i = i; return z; }
static jcx jdiv(jcx a, jcx b) {
    double d = b.r * b.r + b.i * b.i;
    if (d == 0.0) return jc(INFINITY, INFINITY);
    /* Smith's algorithm */
    if (fabs(b.r) >= fabs(b.i)) {
        double t = b.i / b.r, e = b.r + b.i * t;
        return jc((a.r + a.i * t) / e, (a.i - a.r * t) / e);
    }
    double t = b.r / b.i, e = b.r * t + b.i;
    return jc((a.r * t + a.i) / e, (a.i * t - a.r) / e);
}

/* x, y, m as double-doubles. */
static bool jac_core(sdd x, sdd y, sdd m, jcx* S, jcx* C, jcx* D, jcx* AM) {
    double s, c, d, am;
    if (!jac_real(x.hi, x.lo, m.hi, m.lo, &s, &c, &d, &am)) return false;
    if (y.hi == 0.0) {
        *S = jc(s, 0.0); *C = jc(c, 0.0); *D = jc(d, 0.0); *AM = jc(am, 0.0);
        return true;
    }
    sdd m1 = sdd_sub(sdd_from(1.0), m);
    double s1, c1, d1, am1;
    if (!jac_real(y.hi, y.lo, m1.hi, m1.lo, &s1, &c1, &d1, &am1)) return false;
    double mm = m.hi + m.lo;
    double del = c1 * c1 + mm * s * s * s1 * s1;
    if (del == 0.0 || !isfinite(del)) return false;
    *S = jc(s * d1 / del, c * d * s1 * c1 / del);
    *C = jc(c * c1 / del, -s * d * s1 * d1 / del);
    *D = jc(d * c1 * d1 / del, -mm * s * c * s1 / del);
    AM->r = NAN; AM->i = NAN;            /* complex amplitude: from sn, cn below */
    return true;
}

/* sdd product of a double-double r and a double u. */
static sdd jac_scale(sdd r, double u) { return sdd_mul_d(r, u); }

/* sn, cn, dn for complex u, any real m. */
static bool jac_scd(double ur, double ui, double m, jcx* S, jcx* C, jcx* D, jcx* AM) {
    if (m >= 0.0 && m <= 1.0)
        return jac_core(sdd_from(ur), sdd_from(ui), sdd_from(m), S, C, D, AM);
    if (m > 1.0) {
        /* sn(u|m) = sn(r u | 1/m)/r, cn(u|m) = dn(r u|1/m), dn(u|m) = cn(r u|1/m) */
        sdd R = jac_sqrt(sdd_from(m));
        sdd im = sdd_div(sdd_from(1.0), sdd_from(m));
        double r = R.hi + R.lo;
        jcx s, c, d, a;
        if (!jac_core(jac_scale(R, ur), jac_scale(R, ui), im, &s, &c, &d, &a)) return false;
        *S = jc(s.r / r, s.i / r); *C = d; *D = c;
        AM->r = NAN; AM->i = NAN;
        return true;
    }
    /* m < 0: mu = -m/(1-m), v = u sqrt(1-m):
     * sn = sd(v|mu)/sqrt(1-m), cn = cd(v|mu), dn = nd(v|mu) */
    sdd w = sdd_sub(sdd_from(1.0), sdd_from(m));
    sdd R = jac_sqrt(w);
    sdd mu = sdd_div(sdd_from(-m), w);
    double r = R.hi + R.lo;
    jcx s, c, d, a;
    if (!jac_core(jac_scale(R, ur), jac_scale(R, ui), mu, &s, &c, &d, &a)) return false;
    jcx sd = jdiv(s, d);
    *S = jc(sd.r / r, sd.i / r);
    *C = jdiv(c, d);
    *D = jdiv(jc(1.0, 0.0), d);
    AM->r = NAN; AM->i = NAN;
    return true;
}
/* Re K(m) for real m (K(1/m)/sqrt(m) above 1). */
static bool jac_re_k(double m, double* out) {
    if (m < 1.0) return elliptic_machine_k(m, out);
    if (m == 1.0) return false;
    double k;
    if (!elliptic_machine_k(1.0 / m, &k)) return false;
    *out = k / sqrt(m);
    return true;
}

static bool jac_ok(double r, double i, double* rr, double* ri) {
    if (!isfinite(r) || !isfinite(i)) return false;
    double big = fmax(fabs(r), fabs(i));
    if (big != 0.0 && big < DBL_MIN) return false;
    if (fabs(r) < DBL_MIN) r = 0.0;
    if (fabs(i) < DBL_MIN) i = 0.0;
    *rr = r; *ri = i;
    return true;
}

static bool jac_inverse(int which, double v, double m, double* out) {
    if (!(m >= 0.0 && m <= 1.0) || !isfinite(v)) return false;
    double phi;
    if (which == JAC_ISN) {
        if (fabs(v) > 1.0) return false;
        phi = asin(v);
    } else if (which == JAC_ICN) {
        if (fabs(v) > 1.0) return false;
        phi = acos(v);
    } else {
        if (m == 0.0 || fabs(v) > 1.0) return false;
        double s2 = (1.0 - v) * (1.0 + v) / m;
        if (s2 > 1.0 || s2 < 0.0) return false;
        phi = asin(sqrt(s2));
    }
    if (m == 1.0 && fabs(phi) >= 0.5 * M_PI) return false;
    return elliptic_machine_f(phi, m, out) && isfinite(*out);
}

bool sf_machine_jacobi(int which, double ur, double ui, double mr, double mi,
                       double* rr, double* ri) {
    if (which < 0 || which >= JAC_COUNT) return false;
    if (!isfinite(ur) || !isfinite(ui) || !isfinite(mr) || !isfinite(mi)) return false;
    if (which >= JAC_ISN) {
        if (ui != 0.0 || mi != 0.0) return false;
        *ri = 0.0;
        if (which == JAC_IDN && mr > 0.0 && mr < 1.0 && ur >= 0.0 && ur * ur < 1.0 - mr) {
            /* 0 <= v < k': dn(K + i y) = k' cd(y | 1-m) runs from k' down to 0,
             * so the inverse is K + i y with y = F(asin(sqrt((1-m-v^2) /
             * ((1-m)(1-v^2)))) | 1-m) -- the reference's branch (the EllipticF
             * form gives the conjugate there) */
            double K, y, w = 1.0 - mr;
            double s2 = (w - ur * ur) / (w * (1.0 - ur) * (1.0 + ur));
            if (!elliptic_machine_k(mr, &K) || s2 < 0.0 || s2 > 1.0) return false;
            if (!elliptic_machine_f(asin(sqrt(s2)), w, &y)) return false;
            return jac_ok(K, y, rr, ri);
        }
        return jac_inverse(which, ur, mr, rr);
    }
    if (mi == 0.0) {
        jcx S, C, D, AM;
        if (jac_scd(ur, ui, mr, &S, &C, &D, &AM)) {
            if (which == JAC_AM) {
                if (ui == 0.0 && mr >= 0.0 && mr <= 1.0 && isfinite(AM.r)) {
                    *rr = AM.r; *ri = 0.0; return true;
                }
                /* the principal -i log(cn + i sn), shifted by the winding count
                 * 2 Pi round((Pi Re u/(2 Re K) - Re phi)/(2 Pi)) */
                double wr = C.r - S.i, wi = C.i + S.r;
                double pr = atan2(wi, wr), pim = -log(hypot(wr, wi));
                double K;
                if (jac_re_k(mr, &K)) {
                    double k = floor((M_PI * ur / (2.0 * K) - pr) / (2.0 * M_PI) + 0.5);
                    pr += 2.0 * M_PI * k;
                }
                if (ui == 0.0) pim = 0.0;
                return jac_ok(pr, pim, rr, ri);
            }
            jcx v[4] = { S, C, D, jc(1.0, 0.0) };
            jcx p = v[JAC_P[which]], q = v[JAC_Q[which]];
            jcx r = (JAC_Q[which] == 3) ? p : jdiv(p, q);
            if (ui == 0.0) r.i = 0.0;
            if (jac_ok(r.r, r.i, rr, ri)) return true;
        }
    }
    double vr, vi;
    if (!flint_jacobi_d(which, ur, ui, mr, mi, &vr, &vi)) return false;
    return jac_ok(vr, vi, rr, ri);
}

/* ------------------------------------------------------------------ */
/* Builtins                                                            */
/* ------------------------------------------------------------------ */

/* `f` (a pure function source) applied to the copied args and evaluated. */
static Expr* jac_apply(const char* f, Expr* a, Expr* b) {
    Expr* fn = parse_expression(f);
    if (!fn) return NULL;
    Expr* call = b ? expr_new_function(fn, (Expr*[]){ expr_copy(a), expr_copy(b) }, 2)
                   : expr_new_function(fn, (Expr*[]){ expr_copy(a) }, 1);
    return eval_and_free(call);
}

static bool jac_is_int(const Expr* e, int64_t v) {
    return e && e->type == EXPR_INTEGER && e->data.integer == v;
}

/* True when e is visibly negated: a negative number, or Times[c, ...] with a
 * negative rational / integer coefficient c. */
static bool jac_negated(const Expr* e) {
    if (e->type == EXPR_INTEGER) return e->data.integer < 0;
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (h == SYM_Rational && e->data.function.arg_count == 2)
        return jac_negated(e->data.function.args[0]);
    if (h == SYM_Times && e->data.function.arg_count >= 2) {
        const Expr* c = e->data.function.args[0];
        return (c->type == EXPR_INTEGER && c->data.integer < 0) ||
               (c->type == EXPR_FUNCTION && c->data.function.head->type == EXPR_SYMBOL &&
                c->data.function.head->data.symbol.name == SYM_Rational &&
                jac_negated(c));
    }
    return false;
}

/* exact tables: index = function code */
static const char* const JAC_M0[13] = {
    "Sin[#1]&", "Cos[#1]&", "1&", "Cos[#1]&", "Cot[#1]&", "Sec[#1]&", "Csc[#1]&",
    "Sec[#1]&", "1&", "Csc[#1]&", "Tan[#1]&", "Sin[#1]&", "#1&" };
static const char* const JAC_M1[13] = {
    "Tanh[#1]&", "Sech[#1]&", "Sech[#1]&", "1&", "Csch[#1]&", "1&", "Csch[#1]&",
    "Cosh[#1]&", "Cosh[#1]&", "Coth[#1]&", "Sinh[#1]&", "Sinh[#1]&",
    "(2 ArcTan[Tanh[#1/2]])&" };
static const char* const JAC_U0[13] = {
    "0&", "1&", "1&", "1&", "ComplexInfinity&", "1&", "ComplexInfinity&",
    "1&", "1&", "ComplexInfinity&", "0&", "0&", "0&" };

static Expr* jac_builtin(Expr* res, int which) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return builtin_arg_error(JAC_NAME[which], argc, 2, 2);
    Expr* u = res->data.function.args[0];
    Expr* m = res->data.function.args[1];

    SFMachineArgs ma;
    if (sf_machine_scalar_args(res, &ma) && ma.n == 2) {
        double rr, ri;
        if (sf_machine_jacobi(which, ma.re[0], ma.im[0], ma.re[1], ma.im[1], &rr, &ri))
            return (ma.any_complex || ri != 0.0) ? sf_machine_complex_result(rr, ri)
                                                 : expr_new_real(rr);
    }

    if (which <= JAC_AM) {
        /* a ratio with sn below, at a machine zero: the pole */
        if (which < JAC_AM && JAC_Q[which] == 0 && u->type == EXPR_REAL && u->data.real == 0.0)
            return expr_new_symbol("ComplexInfinity");
        if (jac_is_int(m, 0)) return jac_apply(JAC_M0[which], u, NULL);
        if (jac_is_int(m, 1)) return jac_apply(JAC_M1[which], u, NULL);
        if (jac_is_int(u, 0)) return jac_apply(JAC_U0[which], u, NULL);
        if (flint_num_call_inexact(res))
            return flint_num_for_call(flint_num_jacobi(which, u, m), res);
        if (jac_negated(u)) {
            bool odd = which == JAC_AM ||
                       ((JAC_P[which] == 0) != (JAC_Q[which] == 0));
            Expr* nu = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                                     (Expr*[]){ expr_new_integer(-1), expr_copy(u) }, 2));
            Expr* call = expr_new_function(expr_new_symbol(JAC_NAME[which]),
                                           (Expr*[]){ nu, expr_copy(m) }, 2);
            if (!odd) return eval_and_free(call);
            return eval_and_free(expr_new_function(expr_new_symbol("Times"),
                                 (Expr*[]){ expr_new_integer(-1), call }, 2));
        }
        return NULL;
    }

    /* the inverses */
    Expr* v = u;
    if (which == JAC_ISN) {
        if (jac_is_int(v, 0)) return expr_new_integer(0);
        if (jac_is_int(v, 1)) return jac_apply("EllipticK[#2]&", v, m);
        if (jac_is_int(m, 0)) return jac_apply("ArcSin[#1]&", v, NULL);
        if (jac_is_int(m, 1)) return jac_apply("ArcTanh[#1]&", v, NULL);
    } else if (which == JAC_ICN) {
        if (jac_is_int(v, 1)) return expr_new_integer(0);
        if (jac_is_int(v, 0)) return jac_apply("EllipticK[#2]&", v, m);
        if (jac_is_int(m, 0)) return jac_apply("ArcCos[#1]&", v, NULL);
        if (jac_is_int(m, 1)) return jac_apply("ArcSech[#1]&", v, NULL);
    } else {
        if (jac_is_int(v, 1)) return expr_new_integer(0);
        if (jac_is_int(m, 1)) return jac_apply("ArcSech[#1]&", v, NULL);
    }
    if (!flint_num_call_inexact(res)) return NULL;
    static const char* const TPL[3] = {
        "EllipticF[ArcSin[#1], #2]&", "EllipticF[ArcCos[#1], #2]&",
        "EllipticF[ArcSin[Sqrt[(1 - #1^2)/#2]], #2]&" };
    const char* tpl = TPL[which - JAC_ISN];
    if (which == JAC_IDN) {
        /* the K + i y branch for real 0 <= v < Sqrt[1 - m], 0 < m < 1 (see the
         * machine kernel) */
        Expr* vd = jac_apply("N[#1]&", v, NULL);
        Expr* md = jac_apply("N[#1]&", m, NULL);
        if (vd && md && vd->type == EXPR_REAL && md->type == EXPR_REAL) {
            double vv = vd->data.real, mm = md->data.real;
            if (mm > 0.0 && mm < 1.0 && vv >= 0.0 && vv * vv < 1.0 - mm)
                tpl = "(EllipticK[#2] + I EllipticF[ArcSin[Sqrt[(1 - #2 - #1^2)/((1 - #2) (1 - #1^2))]], 1 - #2])&";
        }
        expr_free(vd); expr_free(md);
    }
    Expr* r = jac_apply(tpl, v, m);
    if (r && expr_is_numeric_quantity(r)) return flint_num_for_call(r, res);
    expr_free(r);
    return NULL;
}

#define JAC_WRAP(NAME, CODE) \
    static Expr* builtin_##NAME(Expr* res) { return jac_builtin(res, CODE); }
JAC_WRAP(jacobi_sn, JAC_SN)   JAC_WRAP(jacobi_cn, JAC_CN)   JAC_WRAP(jacobi_dn, JAC_DN)
JAC_WRAP(jacobi_cd, JAC_CD)   JAC_WRAP(jacobi_cs, JAC_CS)   JAC_WRAP(jacobi_dc, JAC_DC)
JAC_WRAP(jacobi_ds, JAC_DS)   JAC_WRAP(jacobi_nc, JAC_NC)   JAC_WRAP(jacobi_nd, JAC_ND)
JAC_WRAP(jacobi_ns, JAC_NS)   JAC_WRAP(jacobi_sc, JAC_SC)   JAC_WRAP(jacobi_sd, JAC_SD)
JAC_WRAP(jacobi_am, JAC_AM)   JAC_WRAP(inv_jacobi_sn, JAC_ISN)
JAC_WRAP(inv_jacobi_cn, JAC_ICN) JAC_WRAP(inv_jacobi_dn, JAC_IDN)

static const char* const JAC_DOC[JAC_COUNT] = {
    "JacobiSN[u, m]\n\tgives the Jacobi elliptic function sn(u | m) (parameter m = k^2).",
    "JacobiCN[u, m]\n\tgives the Jacobi elliptic function cn(u | m).",
    "JacobiDN[u, m]\n\tgives the Jacobi elliptic function dn(u | m).",
    "JacobiCD[u, m]\n\tgives cn(u | m)/dn(u | m).",
    "JacobiCS[u, m]\n\tgives cn(u | m)/sn(u | m).",
    "JacobiDC[u, m]\n\tgives dn(u | m)/cn(u | m).",
    "JacobiDS[u, m]\n\tgives dn(u | m)/sn(u | m).",
    "JacobiNC[u, m]\n\tgives 1/cn(u | m).",
    "JacobiND[u, m]\n\tgives 1/dn(u | m).",
    "JacobiNS[u, m]\n\tgives 1/sn(u | m).",
    "JacobiSC[u, m]\n\tgives sn(u | m)/cn(u | m).",
    "JacobiSD[u, m]\n\tgives sn(u | m)/dn(u | m).",
    "JacobiAmplitude[u, m]\n\tgives the amplitude am(u | m): Sin[am] = JacobiSN[u, m],\n"
    "\tCos[am] = JacobiCN[u, m], continuous in u (am(u + 4K) = am(u) + 2 Pi);\n"
    "\tthe inverse of EllipticF[phi, m] in phi.",
    "InverseJacobiSN[v, m]\n\tgives the inverse of JacobiSN: EllipticF[ArcSin[v], m].",
    "InverseJacobiCN[v, m]\n\tgives the inverse of JacobiCN: EllipticF[ArcCos[v], m].",
    "InverseJacobiDN[v, m]\n\tgives the inverse of JacobiDN:\n"
    "\tEllipticF[ArcSin[Sqrt[(1 - v^2)/m]], m]."
};

void jacobi_elliptic_init(void) {
    Expr* (*fns[JAC_COUNT])(Expr*) = {
        builtin_jacobi_sn, builtin_jacobi_cn, builtin_jacobi_dn, builtin_jacobi_cd,
        builtin_jacobi_cs, builtin_jacobi_dc, builtin_jacobi_ds, builtin_jacobi_nc,
        builtin_jacobi_nd, builtin_jacobi_ns, builtin_jacobi_sc, builtin_jacobi_sd,
        builtin_jacobi_am, builtin_inv_jacobi_sn, builtin_inv_jacobi_cn,
        builtin_inv_jacobi_dn };
    for (int i = 0; i < JAC_COUNT; i++) {
        symtab_add_builtin(JAC_NAME[i], fns[i]);
        symtab_get_def(JAC_NAME[i])->attributes |=
            (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
        char buf[600];
        snprintf(buf, sizeof buf, "%s\n%s", JAC_DOC[i],
                 i <= JAC_AM ? "Exact at u = 0, m = 0, 1; numeric at machine and arbitrary\n"
                               "precision, real and complex; threads over packed arrays and\n"
                               "compiles. D in u: DLMF 22.13; inert in m."
                             : "Real on the principal interval through a machine Carlson kernel;\n"
                               "complex beyond it with EllipticF's branches; compiles.");
        symtab_set_docstring(JAC_NAME[i], buf);
    }
}
