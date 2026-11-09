/* Mathilda -- the regularized incomplete beta function (DLMF 8.17):
 *
 *   BetaRegularized[z, a, b]       = I_z(a, b) = B(z; a, b) / B(a, b),
 *                                    B(z; a, b) = z^a/a 2F1(a, 1 - b; a + 1; z)
 *   BetaRegularized[z0, z1, a, b]  = I_z1(a, b) - I_z0(a, b)
 *
 * Exact rules (the reference system's): a a non-positive integer gives 1, b a
 * non-positive integer 0; a = 1 gives 1 - (1 - z)^b and b = 1 gives z^a - 0^a;
 * an exact numeric z with positive integer a, b the terminating binomial sum
 * Sum[Binomial[a + b - 1, j] z^j (1 - z)^(a + b - 1 - j), {j, a, a + b - 1}]
 * (BetaRegularized[1/2, 2, 3] = 11/16); anything else exact stays unevaluated.
 * D in z (and z0, z1): z^(a-1) (1-z)^(b-1)/Beta[a, b]; in a and b the closed
 * forms through HypergeometricPFQ[{a, a, 1-b}, {a+1, a+1}, z] and PolyGamma.
 *
 * Machine kernel sf_machine_betareg (scalar, NDArray / packed through the n-ary
 * map, Compile[]): for 0 <= z <= 1, a, b > 0, the continued
 * fraction (modified Lentz, DLMF 8.17.22) on whichever of I_z(a, b) and
 * 1 - I_(1-z)(b, a) converges (z against (a + 1)/(a + b + 2)), times the
 * prefactor z^a (1-z)^b / (a B(a, b)) -- from pow (within an ulp whatever the
 * exponent) and the shared machine Beta for small parameters, from the
 * Stirling form with a log1p exponent for a, b >= 10 -- with 1 - z carried as
 * a double-double. A complement that cancels below 1/8, a fraction that
 * does not converge in 600 steps, a prefactor outside the double range, and
 * every complex or out-of-range argument are answered by arb in the kernel
 * (acb_hypgeom_beta_lower, regularized, at 53 bits), so the buffer and
 * Compile[] answers bit-match the scalar one. Arbitrary precision: the same arb
 * function on the shared precision ladder.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "betaregularized.h"

#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "numeric.h"
#include "message.h"
#include "flint_num_bridge.h"
#include "sf_scalar.h"
#include "sf_dd.h"
#include "sf_machine_gamma.h"
#include "symtab.h"
#include "../sym_names.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/* Machine kernel                                                      */
/* ------------------------------------------------------------------ */

/* The continued fraction of I_x(a, b) / prefactor (DLMF 8.17.22), modified
 * Lentz; false when it does not converge. */
static bool br_cf(double a, double b, double x, double* h_out) {
    const double TINY = 1e-300;
    double qab = a + b, qap = a + 1.0, qam = a - 1.0;
    double c = 1.0, d = 1.0 - qab * x / qap;
    if (fabs(d) < TINY) d = TINY;
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m <= 600; m++) {
        double m2 = 2.0 * m;
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d; if (fabs(d) < TINY) d = TINY;
        c = 1.0 + aa / c; if (fabs(c) < TINY) c = TINY;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d; if (fabs(d) < TINY) d = TINY;
        c = 1.0 + aa / c; if (fabs(c) < TINY) c = TINY;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (fabs(del - 1.0) < 3e-16) { *h_out = h; return true; }
    }
    return false;
}

/* Stirling's correction lgamma(t) - ((t - 1/2) log t - t + log(2 pi)/2), t >= 10 */
static double br_stirling(double t) {
    double r = 1.0 / t, r2 = r * r;
    return r * (1.0 / 12.0 + r2 * (-1.0 / 360.0 + r2 * (1.0 / 1260.0 + r2 * (-1.0 / 1680.0
               + r2 * (1.0 / 1188.0 - r2 * (691.0 / 360360.0))))));
}

/* log(1 + u) - u without the cancellation: the series for |u| < 1/4 */
static double br_log1pmx(double u) {
    if (!(fabs(u) < 0.25)) return log1p(u) - u;
    double term = u, sum = 0.0;
    for (int k = 2; k < 80; k++) {
        term *= -u;
        double t = term / k;
        sum += t;
        if (fabs(t) < 1e-18 * fabs(sum)) break;
    }
    return sum;
}

/* s^a t^b / (a B(a, b)) with s + t = 1, both double-double (exact). Small
 * parameters: pow (within an ulp whatever the exponent), the low parts to first
 * order, the shared machine Beta. Large ones (a, b >= 10): the Stirling form
 * (s n/a)^a (t n/b)^b sqrt(a b/(2 pi n)) e^(D(n) - D(a) - D(b)), n = a + b,
 * with s n/a = 1 + (s b - t a)/a through log1p, so the exponent has no
 * cancellation; declined (to arb) when it is large enough that its rounding
 * would cost more than ~64 ulp. */
static bool br_pref(sdd s, double a, sdd t, double b, double* out) {
    double p;
    if (a >= 10.0 && b >= 10.0) {
        double n = a + b;
        sdd num = sdd_sub(sdd_mul_d(s, b), sdd_mul_d(t, a));    /* s n - a */
        double d = num.hi + num.lo;
        /* a log1p(u) + b log1p(v) with a u + b v = 0: the linear parts cancel
         * exactly, so sum a log1pmx(u) + b log1pmx(v) -- terms no larger than
         * the exponent itself */
        double t1 = a * br_log1pmx(d / a), t2 = b * br_log1pmx(-d / b);
        if (!isfinite(t1) || !isfinite(t2) || fabs(t1) + fabs(t2) > 64.0) return false;
        double E = t1 + t2 + br_stirling(n) - br_stirling(a) - br_stirling(b);
        p = sqrt(a * b / (2.0 * M_PI * n)) * exp(E) / a;
    } else {
        double B;
        if (!sfg_beta(a, b, &B) || !(B > 0.0) || !isfinite(B)) return false;
        double ps = pow(s.hi, a), pt = pow(t.hi, b);
        if (s.lo != 0.0) ps *= 1.0 + a * (s.lo / s.hi);
        if (t.lo != 0.0) pt *= 1.0 + b * (t.lo / t.hi);
        p = (ps / a) * (pt / B);
    }
    if (!isfinite(p) || p == 0.0 || fabs(p) < DBL_MIN) return false;
    *out = p;
    return true;
}

/* I_x(a, b), 0 <= x <= 1, a, b > 0 */
static bool br_real(double x, double a, double b, double* out) {
    if (!(a > 0.0) || !(b > 0.0) || !(x >= 0.0) || !(x <= 1.0)) return false;
    if (!(a + b < 1e7) || a < 1e-150 || b < 1e-150) return false;
    if (x == 0.0) { *out = 0.0; return true; }
    if (x == 1.0) { *out = 1.0; return true; }
    sdd X = sdd_from(x), Y = sdd_two_sum(1.0, -x);           /* 1 - x, exact */
    double p, h;
    if (x < (a + 1.0) / (a + b + 2.0)) {
        if (!br_pref(X, a, Y, b, &p) || !br_cf(a, b, x, &h)) return false;
        *out = p * h;
        return isfinite(*out);
    }
    /* 1 - I_(1-x)(b, a): the fraction at 1 - x rounded (a ~1e-16 relative
     * perturbation it carries with condition O(1)), the prefactor exact */
    if (!br_pref(Y, b, X, a, &p) || !br_cf(b, a, Y.hi + Y.lo, &h)) return false;
    double r = 1.0 - p * h;
    if (!(r >= 0.125)) return false;                         /* the complement cancels */
    *out = r;
    return true;
}

static bool br_out(double r, double i, double* rr, double* ri) {
    if (!isfinite(r) || !isfinite(i)) return false;
    double big = fmax(fabs(r), fabs(i));
    if (big != 0.0 && big < DBL_MIN) return false;
    if (fabs(r) < DBL_MIN) r = 0.0;
    if (fabs(i) < DBL_MIN) i = 0.0;
    *rr = r; *ri = i;
    return true;
}

static bool br_nonpos_int(double r, double i) { return i == 0.0 && r <= 0.0 && r == floor(r); }

bool sf_machine_betareg(double zr, double zi, double ar, double ai, double br, double bi,
                        double* rr, double* ri) {
    if (!isfinite(zr) || !isfinite(zi) || !isfinite(ar) || !isfinite(ai) ||
        !isfinite(br) || !isfinite(bi)) return false;
    /* the reference system's limits: a a non-positive integer 1, b one 0 */
    if (br_nonpos_int(ar, ai) && !br_nonpos_int(br, bi)) { *rr = 1.0; *ri = 0.0; return true; }
    if (br_nonpos_int(br, bi) && !br_nonpos_int(ar, ai)) { *rr = 0.0; *ri = 0.0; return true; }
    bool real = zi == 0.0 && ai == 0.0 && bi == 0.0;
    double v;
    if (real && br_real(zr, ar, br, &v)) return br_out(v, 0.0, rr, ri);
    double vr, vi;
    if (!flint_betareg_d(zr, zi, ar, ai, br, bi, &vr, &vi)) return false;
    if (real && zr >= 0.0 && zr <= 1.0) vi = 0.0;
    return br_out(vr, vi, rr, ri);
}

bool sf_machine_betareg4(double z0r, double z0i, double z1r, double z1i, double ar, double ai,
                         double br, double bi, double* rr, double* ri) {
    if (!isfinite(z0r) || !isfinite(z0i) || !isfinite(z1r) || !isfinite(z1i) ||
        !isfinite(ar) || !isfinite(ai) || !isfinite(br) || !isfinite(bi)) return false;
    if ((br_nonpos_int(ar, ai) && !br_nonpos_int(br, bi)) ||
        (br_nonpos_int(br, bi) && !br_nonpos_int(ar, ai))) { *rr = 0.0; *ri = 0.0; return true; }
    bool real = z0i == 0.0 && z1i == 0.0 && ai == 0.0 && bi == 0.0;
    if (z0r == z1r && z0i == z1i) { *rr = 0.0; *ri = 0.0; return true; }
    double v0, v1;
    if (real && br_real(z0r, ar, br, &v0) && br_real(z1r, ar, br, &v1)) {
        double d = v1 - v0;
        if (fabs(d) * 8.0 >= fmax(fabs(v0), fabs(v1))) return br_out(d, 0.0, rr, ri);
    }
    double vr, vi;
    if (!flint_betareg4_d(z0r, z0i, z1r, z1i, ar, ai, br, bi, &vr, &vi)) return false;
    if (real && z0r >= 0.0 && z0r <= 1.0 && z1r >= 0.0 && z1r <= 1.0) vi = 0.0;
    return br_out(vr, vi, rr, ri);
}

/* ------------------------------------------------------------------ */
/* Builtin                                                             */
/* ------------------------------------------------------------------ */

static bool br_int(const Expr* e, int64_t* v) {
    if (e->type != EXPR_INTEGER) return false;
    *v = e->data.integer;
    return true;
}

static bool br_exact_number(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           (e->data.function.head->data.symbol.name == SYM_Rational ||
            (e->data.function.head->data.symbol.name == SYM_Complex &&
             br_exact_number(e->data.function.args[0]) && br_exact_number(e->data.function.args[1])));
}

static Expr* br_apply(const char* tmpl, Expr* const* args, size_t n) {
    Expr* f = parse_expression(tmpl);
    if (!f) return NULL;
    Expr* a[4];
    for (size_t i = 0; i < n; i++) a[i] = expr_copy(args[i]);
    return eval_and_free(expr_new_function(f, a, n));
}

static bool br_leaf(const Expr* e, double* re, double* im, bool* inexact, bool* cplx) {
    *cplx = false;
    if (sf_machine_leaf(e, re, inexact)) { *im = 0.0; return true; }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Complex &&
        e->data.function.arg_count == 2) {
        bool a, b;
        if (!sf_machine_leaf(e->data.function.args[0], re, &a) ||
            !sf_machine_leaf(e->data.function.args[1], im, &b)) return false;
        *inexact = a || b;
        *cplx = true;
        return true;
    }
    return false;
}

/* the exact rules for I_z(a, b) at one endpoint (NULL: none applies) */
static Expr* br_exact(Expr* z, Expr* a, Expr* b) {
    int64_t ia, ib;
    bool ha = br_int(a, &ia), hb = br_int(b, &ib);
    if (ha && ia <= 0 && !(hb && ib <= 0)) return expr_new_integer(1);
    if (hb && ib <= 0 && !(ha && ia <= 0)) return expr_new_integer(0);
    if (ha && ia == 1) return br_apply("(1 - (1 - #1)^#2)&", (Expr*[]){ z, b }, 2);
    if (hb && ib == 1) return br_apply("(#1^#2 - 0^#2)&", (Expr*[]){ z, a }, 2);
    if (ha && hb && ia > 0 && ib > 0 && ia + ib <= 2000 && br_exact_number(z))
        return br_apply("Sum[Binomial[#2 + #3 - 1, j] #1^j (1 - #1)^(#2 + #3 - 1 - j), {j, #2, #2 + #3 - 1}]&",
                        (Expr*[]){ z, a, b }, 3);
    return NULL;
}

static long br_bits(Expr* const* a, int n) {
    long bits = 0;
    for (int i = 0; i < n; i++) {
        long b = numeric_min_inexact_bits(a[i]);
        if (b > 0 && (bits == 0 || b < bits)) bits = b;
    }
    return bits;
}

static Expr* builtin_betaregularized(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 3 && argc != 4) return builtin_arg_error("BetaRegularized", argc, 3, 4);
    Expr** A = res->data.function.args;
    double re[4], im[4];
    bool ix[4], cx[4], any_ix = false, any_cx = false, mach = true;
    for (size_t i = 0; i < argc; i++) {
        if (!br_leaf(A[i], &re[i], &im[i], &ix[i], &cx[i])) { mach = false; break; }
        any_ix |= ix[i]; any_cx |= cx[i];
    }
    if (mach && any_ix) {
        double rr, ri;
        bool ok = argc == 3 ? sf_machine_betareg(re[0], im[0], re[1], im[1], re[2], im[2], &rr, &ri)
                            : sf_machine_betareg4(re[0], im[0], re[1], im[1], re[2], im[2],
                                                  re[3], im[3], &rr, &ri);
        if (ok) return (!any_cx && ri == 0.0) ? expr_new_real(rr) : sf_machine_complex_result(rr, ri);
    }
    bool numeric = true;
    for (size_t i = 0; i < argc; i++) numeric &= expr_is_numeric_quantity(A[i]);
    long bits = br_bits(A, (int)argc);
    if (!numeric || bits == 0) {
        /* exact (or symbolic) arguments */
        if (argc == 3) return br_exact(A[0], A[1], A[2]);
        Expr* hi = br_exact(A[1], A[2], A[3]);
        if (!hi) return NULL;
        Expr* lo = br_exact(A[0], A[2], A[3]);
        if (!lo) { expr_free(hi); return NULL; }
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
               (Expr*[]){ hi, expr_new_function(expr_new_symbol(SYM_Times),
                                                (Expr*[]){ expr_new_integer(-1), lo }, 2) }, 2));
    }
    Expr* r = argc == 3 ? flint_num_betareg(A[0], A[1], A[2])
                        : flint_num_betareg4(A[0], A[1], A[2], A[3]);
    if (r && bits <= 53) r = numeric_lower_mpfr_to_machine(r);
    return r;
}

void betaregularized_init(void) {
    symtab_add_builtin("BetaRegularized", builtin_betaregularized);
    symtab_get_def("BetaRegularized")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_set_docstring("BetaRegularized",
        "BetaRegularized[z, a, b]\n"
        "\tgives the regularized incomplete beta function I_z(a, b) = Beta[z, a, b]/Beta[a, b].\n"
        "BetaRegularized[z0, z1, a, b]\n"
        "\tgives I_z1(a, b) - I_z0(a, b).\n"
        "Numeric at machine and arbitrary precision, real and complex (analytic\n"
        "continuation through z^a/a 2F1(a, 1-b; a+1; z)); threads over packed arrays\n"
        "and compiles. D[BetaRegularized[z, a, b], z] = z^(a-1) (1-z)^(b-1)/Beta[a, b].");
}
