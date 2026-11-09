/* Mathilda -- EllipticTheta[a, u, q] and EllipticThetaPrime[a, u, q], the Jacobi
 * theta functions of NOME q (|q| < 1) and their u-derivatives (DLMF 20.2):
 *
 *   theta_1 = 2 Sum_{n>=0} (-1)^n q^((n+1/2)^2) sin((2n+1) u)
 *   theta_2 = 2 Sum_{n>=0}        q^((n+1/2)^2) cos((2n+1) u)
 *   theta_3 = 1 + 2 Sum_{n>=1}        q^(n^2) cos(2 n u)
 *   theta_4 = 1 + 2 Sum_{n>=1} (-1)^n q^(n^2) cos(2 n u)
 *
 * with q^(1/4) principal; EllipticTheta[a, q] = EllipticTheta[a, 0, q].
 *
 * Exact rules (the reference system's): q = 0 gives 0, 0, 1, 1 (and 0 for every
 * derivative); the vanishing values theta_1(0, q), theta'_2(0, q),
 * theta'_3(0, q), theta'_4(0, q) are ConditionalExpression[0, Abs[q] < 1]; any
 * other exact call stays unevaluated; a first argument other than 1..4 gives
 * EllipticTheta::etype. D in u: EllipticTheta -> EllipticThetaPrime; the other
 * derivatives stay inert.
 *
 * Numerics. Machine kernel sf_machine_theta (scalar, NDArray / packed map,
 * Compile[]): the q-series in double complex arithmetic for |q| <= 0.98 with
 * q^(n(n+1)) and q^(n^2) built by recurrence; it declines when the series
 * cancels by more than 2^8 (theta_4 near |q| = 1 is exponentially small) or the
 * terms grow (large |Im u|), and arb answers in the kernel: acb_modular_theta /
 * acb_modular_theta_jet at tau = log(q)/(i pi), z = u/pi, which reduces tau to
 * the fundamental domain (the modular transformation) and so is uniform up to
 * |q| -> 1. Arbitrary precision: the same arb route.
 *
 * Attributes: Listable, NHoldFirst, NumericFunction, Protected.
 */
#include "elliptic_theta.h"

#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "numeric.h"
#include <stdlib.h>
#include "message.h"
#include "flint_num_bridge.h"
#include "sf_scalar.h"
#include "sf_dd.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "symtab.h"
#include "../sym_names.h"

/* ------------------------------------------------------------------ */
/* Machine kernel                                                      */
/* ------------------------------------------------------------------ */

typedef struct { double r, i; } tcx;
static tcx tc(double r, double i) { tcx z; z.r = r; z.i = i; return z; }
static tcx tmul(tcx a, tcx b) { return tc(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r); }
static tcx tadd(tcx a, tcx b) { return tc(a.r + b.r, a.i + b.i); }
static tcx tscale(tcx a, double s) { return tc(a.r * s, a.i * s); }
static double tabs(tcx a) { return hypot(a.r, a.i); }
static tcx tsin(tcx z) { return tc(sin(z.r) * cosh(z.i), cos(z.r) * sinh(z.i)); }
static tcx tcos(tcx z) { return tc(cos(z.r) * cosh(z.i), -sin(z.r) * sinh(z.i)); }

/* principal q^(1/4) */
static tcx tqrt4(tcx q) {
    double m = pow(tabs(q), 0.25), a = 0.25 * atan2(q.i, q.r);
    return tc(m * cos(a), m * sin(a));
}

static bool theta_series(int a, int prime, tcx u, tcx q, tcx* out) {
    double aq = tabs(q);
    if (!(aq <= 0.9)) return false;       /* beyond: arb (modular reduction) */
    if (aq == 0.0) {
        *out = tc(0.0, 0.0);
        if (!prime && (a == 3 || a == 4)) *out = tc(1.0, 0.0);
        return true;
    }
    if (fabs(u.i) > 30.0) return false;
    tcx sum = tc(0.0, 0.0);
    double peak = 0.0;
    if (a == 1 || a == 2) {
        /* q^(n(n+1)), n = 0, 1, ...: ratio q^(2n) */
        tcx qp = tc(1.0, 0.0), q2n = tc(1.0, 0.0), q2 = tmul(q, q);
        for (int n = 0; n < 200; n++) {
            double k = 2.0 * n + 1.0;
            tcx arg = tscale(u, k);
            tcx t;
            if (a == 1) t = prime ? tscale(tcos(arg), k) : tsin(arg);
            else        t = prime ? tscale(tsin(arg), -k) : tcos(arg);
            if (a == 1 && (n & 1)) t = tscale(t, -1.0);
            t = tmul(t, qp);
            sum = tadd(sum, t);
            double at = tabs(t);
            if (at > peak) peak = at;
            if (at < 1e-18 * tabs(sum) && n > 2) break;
            if (!isfinite(at)) return false;
            q2n = (n == 0) ? q2 : tmul(q2n, q2);        /* q^(2(n+1)) */
            qp = tmul(qp, q2n);                         /* q^((n+1)(n+2)) */
            if (n == 199) return false;
        }
        sum = tscale(tmul(sum, tqrt4(q)), 2.0);
        peak *= 2.0 * pow(aq, 0.25);
    } else {
        /* 1 + 2 Sum q^(n^2) cos(2 n u): q^(n^2) = q^((n-1)^2) q^(2n-1) */
        tcx qn = tc(1.0, 0.0), qodd = q, q2 = tmul(q, q);
        if (!prime) { sum = tc(0.5, 0.0); peak = 0.5; }
        for (int n = 1; n < 200; n++) {
            qn = tmul(qn, qodd);                        /* q^(n^2) */
            qodd = tmul(qodd, q2);
            double k = 2.0 * n;
            tcx arg = tscale(u, k);
            tcx t = prime ? tscale(tsin(arg), -k) : tcos(arg);
            if (a == 4 && (n & 1)) t = tscale(t, -1.0);
            t = tmul(t, qn);
            sum = tadd(sum, t);
            double at = tabs(t);
            if (at > peak) peak = at;
            if (!isfinite(at)) return false;
            if (at < 1e-18 * tabs(sum) && n > 2) break;
            if (n == 199) return false;
        }
        sum = tscale(sum, 2.0);
        peak *= 2.0;
    }
    if (!isfinite(sum.r) || !isfinite(sum.i)) return false;
    if (peak > 8.0 * tabs(sum)) return false;           /* cancels: arb */
    *out = sum;
    return true;
}

static bool th_ok(double r, double i, double* rr, double* ri) {
    if (!isfinite(r) || !isfinite(i)) return false;
    double big = fmax(fabs(r), fabs(i));
    if (big != 0.0 && big < DBL_MIN) return false;
    if (fabs(r) < DBL_MIN) r = 0.0;
    if (fabs(i) < DBL_MIN) i = 0.0;
    *rr = r; *ri = i;
    return true;
}

bool sf_machine_theta(int a, int prime, double ur, double ui, double qr, double qi,
                      double* rr, double* ri) {
    if (a < 1 || a > 4) return false;
    if (!isfinite(ur) || !isfinite(ui) || !isfinite(qr) || !isfinite(qi)) return false;
    if (!(hypot(qr, qi) < 1.0)) return false;
    /* reduce Re u by multiples of pi (theta_1, theta_2 change sign, theta_3,
     * theta_4 do not), the residual formed in double-double, so the series'
     * (2n+1) u never carries a large rounded argument */
    double j = floor(ur / M_PI + 0.5), sg = 1.0, urr = ur;
    if (j != 0.0 && fabs(j) < 1e15) {
        sdd p = sdd_mul_d(sdd_quick(SDD_PI_HI, SDD_PI_LO), j);
        sdd r = sdd_sub(sdd_from(ur), p);
        urr = r.hi + r.lo;
        if (a <= 2 && fmod(fabs(j), 2.0) == 1.0) sg = -1.0;
    }
    tcx v;
    if (theta_series(a, prime, tc(urr, ui), tc(qr, qi), &v)) {
        v.r *= sg; v.i *= sg;
        /* real u and q >= 0 (or q < 0 for theta_3, theta_4): a real value */
        if (ui == 0.0 && qi == 0.0 && (qr >= 0.0 || a >= 3)) v.i = 0.0;
        return th_ok(v.r, v.i, rr, ri);
    }
    double vr, vi;
    if (!flint_theta_d(a, prime, ur, ui, qr, qi, &vr, &vi)) return false;
    if (ui == 0.0 && qi == 0.0 && (qr >= 0.0 || a >= 3)) vi = 0.0;
    return th_ok(vr, vi, rr, ri);
}

/* ------------------------------------------------------------------ */
/* Builtins                                                            */
/* ------------------------------------------------------------------ */

static bool th_is_int(const Expr* e, int64_t v) {
    return e && e->type == EXPR_INTEGER && e->data.integer == v;
}

static Expr* th_cond_zero(const Expr* q) {
    Expr* f = parse_expression("(ConditionalExpression[0, Abs[#] < 1])&");
    return eval_and_free(expr_new_function(f, (Expr*[]){ expr_copy((Expr*)q) }, 1));
}

static Expr* th_builtin(Expr* res, int prime) {
    const char* name = prime ? "EllipticThetaPrime" : "EllipticTheta";
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 2) {
        /* EllipticTheta[a, q] = EllipticTheta[a, 0, q], kept in its own form
         * when that stays unevaluated (the reference prints EllipticTheta[3, q]) */
        Expr* r = eval_and_free(expr_new_function(expr_new_symbol(name),
            (Expr*[]){ expr_copy(res->data.function.args[0]), expr_new_integer(0),
                       expr_copy(res->data.function.args[1]) }, 3));
        if (r && r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL &&
            strcmp(r->data.function.head->data.symbol.name, name) == 0) {
            expr_free(r);
            return NULL;
        }
        return r;
    }
    if (argc != 3) return builtin_arg_error(name, argc, 2, 3);
    Expr* A = res->data.function.args[0];
    Expr* u = res->data.function.args[1];
    Expr* q = res->data.function.args[2];
    if (A->type != EXPR_INTEGER) {
        if (expr_is_numeric_quantity(A)) {
            char* s = expr_to_string(A);
            mth_message(name, "etype", "The first argument %s is expected to be 1, 2, 3, or 4.", s);
            free(s);
        }
        return NULL;
    }
    int64_t a = A->data.integer;
    if (a < 1 || a > 4) {
        mth_message(name, "etype", "The first argument %lld is expected to be 1, 2, 3, or 4.",
                    (long long)a);
        return NULL;
    }

    SFMachineArgs ma;
    if (sf_machine_scalar_args(res, &ma) && ma.n == 3) {
        double rr, ri;
        if (sf_machine_theta((int)a, prime, ma.re[1], ma.im[1], ma.re[2], ma.im[2], &rr, &ri)) {
            bool real_out = !ma.any_complex && ri == 0.0;
            return real_out ? expr_new_real(rr) : sf_machine_complex_result(rr, ri);
        }
    }
    if (th_is_int(q, 0))
        return expr_new_integer(!prime && (a == 3 || a == 4) ? 1 : 0);
    if (th_is_int(u, 0)) {
        if ((!prime && a == 1) || (prime && a != 1)) return th_cond_zero(q);
    }
    if (flint_num_call_inexact(res))
        return flint_num_for_call(flint_num_theta((int)a, prime, u, q), res);
    return NULL;
}

static Expr* builtin_elliptic_theta(Expr* res) { return th_builtin(res, 0); }
static Expr* builtin_elliptic_theta_prime(Expr* res) { return th_builtin(res, 1); }

void elliptic_theta_init(void) {
    symtab_add_builtin("EllipticTheta", builtin_elliptic_theta);
    symtab_add_builtin("EllipticThetaPrime", builtin_elliptic_theta_prime);
    const char* names[2] = { "EllipticTheta", "EllipticThetaPrime" };
    for (int i = 0; i < 2; i++)
        symtab_get_def(names[i])->attributes |=
            (ATTR_LISTABLE | ATTR_NHOLDFIRST | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_set_docstring("EllipticTheta",
        "EllipticTheta[a, u, q]\n"
        "\tgives the Jacobi theta function theta_a(u, q), a = 1, 2, 3, 4, of nome q:\n"
        "\ttheta_1 = 2 Sum (-1)^n q^((n+1/2)^2) Sin[(2n+1) u], theta_2 the same without\n"
        "\t(-1)^n and with Cos, theta_3 = 1 + 2 Sum q^(n^2) Cos[2 n u], theta_4 with (-1)^n.\n"
        "EllipticTheta[a, q]\n\tgives theta_a(0, q).\n"
        "Numeric at machine and arbitrary precision, real and complex (arb's modular\n"
        "reduction near |q| = 1); threads over packed arrays and compiles.\n"
        "D[EllipticTheta[a, u, q], u] = EllipticThetaPrime[a, u, q].");
    symtab_set_docstring("EllipticThetaPrime",
        "EllipticThetaPrime[a, u, q]\n"
        "\tgives the derivative of EllipticTheta[a, u, q] with respect to u.\n"
        "EllipticThetaPrime[a, q]\n\tgives the derivative at u = 0.\n"
        "Numeric at machine and arbitrary precision, real and complex; threads over\n"
        "packed arrays and compiles.");
}
