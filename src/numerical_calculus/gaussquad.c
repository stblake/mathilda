/*
 * gaussquad.c
 *
 * Gauss-Legendre quadrature rules.
 *
 *   GaussianQuadratureWeights[n, a, b]        machine-precision rule
 *   GaussianQuadratureWeights[n, a, b, prec]  rule at `prec` decimal digits
 *       -> {{x1, w1}, ..., {xn, wn}},  x1 < x2 < ... < xn,
 *          so that Integrate[f[x], {x, a, b}] ~ Sum[wi f[xi]] exactly for
 *          every polynomial of degree <= 2n - 1.
 *
 *   GaussianQuadratureError[n, f, a, b]
 *       -> the error term  Integrate[f, {x, a, b}] - Sum[wi f[xi]]
 *          = (b - a)^(2n+1) (n!)^4 / ((2n + 1) ((2n)!)^3) f^(2n)(xi),
 *          returned as  c Derivative[2n][f]  (the derivative is taken at an
 *          unspecified point xi of (a, b)).
 *
 * Method.  The nodes are the zeros of the Legendre polynomial P_n, found by
 * Newton's method on the three-term recurrence
 *     (j + 1) P_{j+1}(x) = (2j + 1) x P_j(x) - j P_{j-1}(x),
 *     P_n'(x) = n (P_{n-1}(x) - x P_n(x)) / (1 - x^2),
 * from Tricomi's asymptotic initial guesses
 *     x_k ~ (1 - 1/(8 n^2) + 1/(8 n^3)) cos(pi (4k - 1) / (4n + 2)),
 * which lie inside the basin of the k-th zero for every n.  The weights are
 *     w_k = 2 / ((1 - x_k^2) P_n'(x_k)^2).
 * Only the zeros in (0, 1) are iterated; the rule is mirrored (x -> -x,
 * w -> w), so the computed rule is exactly symmetric.  Each Newton step costs
 * O(n), so a whole rule costs O(n^2) per sweep: n = 1000 at machine precision
 * is a few milliseconds.
 *
 * Machine precision: Newton in double.  The weight is taken from the last
 * Newton step by a first-order expansion about the pre-update iterate, so it
 * belongs to the unrounded zero rather than to its 53-bit rounding (which near
 * the ends would cost ~1e-11 relative at n = 1000); `1 - x^2` is formed as
 * (1 - x)(1 + x), exact-relative near the ends.
 *
 * Arbitrary precision: the double nodes seed Newton in MPFR at the target
 * precision plus guard bits.  The guard covers (i) the rounding error of the
 * recurrence, which grows like n, and (ii) the amplification of a node error
 * in the weight, |d log w / dx| = 2|x|/(1 - x^2) <= ~n^2 at the outermost
 * node; 32 + 3 log2(n) bits bounds both.  Newton runs until the correction is
 * below 2^-(working bits), then one more step (quadratic convergence leaves
 * that last step at the rounding level), so every returned digit is correct.
 * The rule is then affinely mapped onto [a, b] at the working precision and
 * rounded once to `prec` digits.
 *
 * Neither head is element-wise over numbers: GaussianQuadratureWeights builds
 * a rule from an integer order, GaussianQuadratureError is symbolic.  Both
 * are therefore off the packed/NDArray/Compile surfaces (see the EXEMPT lists
 * of the audit tools); a machine rule is returned as a plain n x 2 real List,
 * which the evaluator packs.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "sym_names.h"
#include "numeric.h"
#include "message.h"
#include "common.h"

#define GQ_HEAD  "GaussianQuadratureWeights"
#define GQE_HEAD "GaussianQuadratureError"

/* Largest order accepted: the rule is O(n^2) to build and O(n) to store. */
#define GQ_MAX_N 100000000

/* ------------------------------------------------------------------ *
 *  Machine-precision rule on [-1, 1]                                  *
 * ------------------------------------------------------------------ */

/* P_n(x) and P_{n-1}(x) by the three-term recurrence, with the coefficients
 * ca[j] = (2j+1)/(j+1), cb[j] = j/(j+1) precomputed (no division per step). */
static void gq_legendre(int n, double x, const double* ca, const double* cb,
                        double* pn, double* pn1) {
    double p0 = 1.0, p1 = x;
    for (int j = 1; j < n; j++) {
        double p2 = ca[j] * x * p1 - cb[j] * p0;
        p0 = p1; p1 = p2;
    }
    *pn = p1; *pn1 = p0;
}

/* Tricomi's estimate of the k-th largest zero of P_n (k = 1 .. n). */
static double gq_initial_guess(int n, int k) {
    double dn = (double)n;
    double th = M_PI * (4.0 * k - 1.0) / (4.0 * dn + 2.0);
    return (1.0 - 1.0 / (8.0 * dn * dn) + 1.0 / (8.0 * dn * dn * dn)) * cos(th);
}

/*
 * gq_rule_machine:
 *   Fill x[0..n) (increasing) and w[0..n) with the n-point Gauss-Legendre
 *   rule on [-1, 1].
 */
static void gq_rule_machine(int n, double* x, double* w) {
    double* ca = malloc(sizeof(double) * (size_t)(n + 1));
    double* cb = malloc(sizeof(double) * (size_t)(n + 1));
    if (!ca || !cb) {
        free(ca); free(cb);
        for (int i = 0; i < n; i++) { x[i] = NAN; w[i] = NAN; }
        return;
    }
    for (int j = 1; j <= n; j++) {
        ca[j] = (2.0 * j + 1.0) / (j + 1.0);
        cb[j] = (double)j / (j + 1.0);
    }
    double dn = (double)n, nn1 = dn * (dn + 1.0);
    int half = n / 2;
    for (int k = 1; k <= half; k++) {
        double z = gq_initial_guess(n, k);
        double wk = 0.0;
        for (int it = 0; it < 100; it++) {
            double pn, pn1;
            gq_legendre(n, z, ca, cb, &pn, &pn1);
            double omz2 = (1.0 - z) * (1.0 + z);
            double dp = dn * (pn1 - z * pn) / omz2;
            double dz = pn / dp;
            /* Weight at the updated node z - dz, by a first-order expansion
             * about z: the double z - dz is the node rounded to 53 bits, and
             * a weight evaluated there inherits a relative error of
             * ulp(z) 2|z|/(1 - z^2) (1e-11 at the ends of a 1000-point rule).
             * Expanding keeps the unrounded correction:
             *   1 - (z - dz)^2 = omz2 + 2 z dz - dz^2,
             *   P_n'(z - dz) = dp - P_n''(z) dz,
             *   P_n'' = (2 z P_n' - n (n + 1) P_n) / (1 - z^2)  (Legendre ODE). */
            double d2p = (2.0 * z * dp - nn1 * pn) / omz2;
            double om_t = omz2 + dz * (2.0 * z - dz);
            double dp_t = dp - d2p * dz;
            wk = 2.0 / (om_t * dp_t * dp_t);
            z -= dz;
            if (fabs(dz) <= 2e-16 * fabs(z)) break;
        }
        /* k-th largest zero: positions n-k (positive) and k-1 (negative) */
        x[n - k] = z;   w[n - k] = wk;
        x[k - 1] = -z;  w[k - 1] = wk;
    }
    if (n % 2 == 1) {
        /* x = 0: w = 2 / P_n'(0)^2, P_n'(0) = n P_{n-1}(0). */
        double pn, pn1;
        gq_legendre(n, 0.0, ca, cb, &pn, &pn1);
        double dp = dn * pn1;
        x[half] = 0.0;
        w[half] = 2.0 / (dp * dp);
    }
    free(ca); free(cb);
}

/* ------------------------------------------------------------------ *
 *  Argument helpers                                                   *
 * ------------------------------------------------------------------ */

/* Positive machine integer n, else false. */
static bool gq_order(const Expr* e, int* out) {
    if (!e || e->type != EXPR_INTEGER) return false;
    if (e->data.integer < 1 || e->data.integer > GQ_MAX_N) return false;
    *out = (int)e->data.integer;
    return true;
}

/* Machine value of a real numeric endpoint (exact, inexact or NumericQ). */
static bool gq_real_machine(const Expr* e, double* out) {
    if (e->type == EXPR_REAL) { *out = e->data.real; return isfinite(*out); }
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    Expr* v = numericalize(e, numeric_machine_spec());
    bool ok = false;
    if (v && v->type == EXPR_REAL && isfinite(v->data.real)) { *out = v->data.real; ok = true; }
    else if (v && v->type == EXPR_INTEGER) { *out = (double)v->data.integer; ok = true; }
    expr_free(v);
    return ok;
}

/*
 * Precision argument: a positive real number of digits, or MachinePrecision.
 * Sets *machine for MachinePrecision (and for a request no larger than machine
 * precision expressed as the machine value itself).
 */
static bool gq_precision(const Expr* e, double* digits, bool* machine) {
    *machine = false;
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_MachinePrecision) {
        *machine = true;
        return true;
    }
    double d;
    if (e->type == EXPR_INTEGER) d = (double)e->data.integer;
    else if (e->type == EXPR_REAL) d = e->data.real;
    else return false;
    if (!(d >= 1.0) || d > 1e6) return false;
    *digits = d;
    return true;
}

static Expr* gq_pair(Expr* a, Expr* b) {
    Expr* v[2] = { a, b };
    return expr_new_function(expr_new_symbol(SYM_List), v, 2);
}

/* ------------------------------------------------------------------ *
 *  Machine rule on [a, b] as an Expr                                  *
 * ------------------------------------------------------------------ */

static Expr* gq_machine_expr(int n, double a, double b) {
    double* x = malloc(sizeof(double) * (size_t)n);
    double* w = malloc(sizeof(double) * (size_t)n);
    Expr** rows = malloc(sizeof(Expr*) * (size_t)n);
    if (!x || !w || !rows) { free(x); free(w); free(rows); return NULL; }
    gq_rule_machine(n, x, w);
    double h = 0.5 * (b - a), c = 0.5 * (a + b);
    for (int i = 0; i < n; i++)
        rows[i] = gq_pair(expr_new_real(c + h * x[i]), expr_new_real(h * w[i]));
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    free(rows); free(x); free(w);
    return out;
}

#ifdef USE_MPFR
/* ------------------------------------------------------------------ *
 *  Arbitrary-precision rule                                           *
 * ------------------------------------------------------------------ */

/* P_n(z), P_{n-1}(z) in MPFR; t0/t1/t2 are caller scratch at the same
 * precision. */
static void gq_legendre_mp(int n, const mpfr_t z, mpfr_t pn, mpfr_t pn1,
                           mpfr_t t0, mpfr_t t1, mpfr_t t2) {
    mpfr_set_ui(t0, 1, MPFR_RNDN);
    mpfr_set(t1, z, MPFR_RNDN);
    for (int j = 1; j < n; j++) {
        /* t2 = ((2j+1) z t1 - j t0) / (j+1) */
        mpfr_mul(t2, z, t1, MPFR_RNDN);
        mpfr_mul_ui(t2, t2, (unsigned long)(2 * j + 1), MPFR_RNDN);
        mpfr_mul_ui(pn, t0, (unsigned long)j, MPFR_RNDN);   /* pn as scratch */
        mpfr_sub(t2, t2, pn, MPFR_RNDN);
        mpfr_div_ui(t2, t2, (unsigned long)(j + 1), MPFR_RNDN);
        mpfr_swap(t0, t1);
        mpfr_swap(t1, t2);
    }
    mpfr_set(pn, t1, MPFR_RNDN);
    mpfr_set(pn1, t0, MPFR_RNDN);
}

/* dp = n (pn1 - z pn) / ((1 - z)(1 + z)); omz2 receives (1 - z)(1 + z). */
static void gq_deriv_mp(int n, const mpfr_t z, const mpfr_t pn, const mpfr_t pn1,
                        mpfr_t dp, mpfr_t omz2, mpfr_t s) {
    mpfr_ui_sub(omz2, 1, z, MPFR_RNDN);
    mpfr_add_ui(s, z, 1, MPFR_RNDN);
    mpfr_mul(omz2, omz2, s, MPFR_RNDN);
    mpfr_mul(s, z, pn, MPFR_RNDN);
    mpfr_sub(dp, pn1, s, MPFR_RNDN);
    mpfr_mul_ui(dp, dp, (unsigned long)n, MPFR_RNDN);
    mpfr_div(dp, dp, omz2, MPFR_RNDN);
}

/*
 * gq_mpfr_expr:
 *   The rule on [a, b] (`ae`, `be` borrowed real numeric Exprs) at `digits`
 *   decimal digits, or NULL when an endpoint is not real numeric.
 */
static Expr* gq_mpfr_expr(int n, const Expr* ae, const Expr* be, double digits) {
    long out_bits = numeric_digits_to_bits(digits);
    long guard = 32 + 3 * (long)ceil(log2((double)n + 1.0));
    mpfr_prec_t wb = (mpfr_prec_t)(out_bits + guard);

    /* endpoints at the working precision */
    mpfr_t a, b, im;
    mpfr_inits2(wb, a, b, im, (mpfr_ptr)0);
    bool inexact;
    Expr* an = numericalize(ae, (NumericSpec){ NUMERIC_MODE_MPFR, (long)wb });
    Expr* bn = numericalize(be, (NumericSpec){ NUMERIC_MODE_MPFR, (long)wb });
    bool ok = an && bn && get_approx_mpfr(an, a, im, &inexact) && mpfr_zero_p(im)
              && get_approx_mpfr(bn, b, im, &inexact) && mpfr_zero_p(im)
              && mpfr_number_p(a) && mpfr_number_p(b);
    expr_free(an); expr_free(bn);
    if (!ok) { mpfr_clears(a, b, im, (mpfr_ptr)0); return NULL; }

    double* xd = malloc(sizeof(double) * (size_t)n);
    double* wd = malloc(sizeof(double) * (size_t)n);
    Expr** rows = malloc(sizeof(Expr*) * (size_t)n);
    if (!xd || !wd || !rows) {
        free(xd); free(wd); free(rows); mpfr_clears(a, b, im, (mpfr_ptr)0); return NULL;
    }
    gq_rule_machine(n, xd, wd);   /* seeds */

    mpfr_t z, pn, pn1, dp, omz2, t0, t1, t2, dz, h, c, xv, wv;
    mpfr_inits2(wb, z, pn, pn1, dp, omz2, t0, t1, t2, dz, h, c, xv, wv, (mpfr_ptr)0);
    mpfr_sub(h, b, a, MPFR_RNDN); mpfr_div_2ui(h, h, 1, MPFR_RNDN);
    mpfr_add(c, a, b, MPFR_RNDN); mpfr_div_2ui(c, c, 1, MPFR_RNDN);

    mpfr_prec_t ob = (mpfr_prec_t)out_bits;
    int half = n / 2;
    for (int k = 1; k <= half; k++) {
        /* k-th largest zero, seeded from the double rule */
        mpfr_set_d(z, xd[n - k], MPFR_RNDN);
        int extra = 0;
        for (int it = 0; it < 200; it++) {
            gq_legendre_mp(n, z, pn, pn1, t0, t1, t2);
            gq_deriv_mp(n, z, pn, pn1, dp, omz2, t0);
            mpfr_div(dz, pn, dp, MPFR_RNDN);
            mpfr_sub(z, z, dz, MPFR_RNDN);
            /* converged when |dz| <= 2^-(wb - 4) |z|; then one more step */
            if (mpfr_zero_p(dz)
                || mpfr_get_exp(dz) <= mpfr_get_exp(z) - (mpfr_exp_t)wb + 4) {
                if (++extra >= 2) break;
            }
        }
        gq_legendre_mp(n, z, pn, pn1, t0, t1, t2);
        gq_deriv_mp(n, z, pn, pn1, dp, omz2, t0);
        /* w = 2 / (omz2 dp^2) */
        mpfr_sqr(wv, dp, MPFR_RNDN);
        mpfr_mul(wv, wv, omz2, MPFR_RNDN);
        mpfr_ui_div(wv, 2, wv, MPFR_RNDN);
        mpfr_mul(wv, wv, h, MPFR_RNDN);

        for (int sgn = 0; sgn < 2; sgn++) {
            int pos = sgn == 0 ? n - k : k - 1;
            /* x = c +/- h z */
            mpfr_mul(xv, h, z, MPFR_RNDN);
            if (sgn == 0) mpfr_add(xv, c, xv, MPFR_RNDN);
            else          mpfr_sub(xv, c, xv, MPFR_RNDN);
            Expr* xe = expr_new_mpfr_bits(ob);
            mpfr_set(xe->data.mpfr, xv, MPFR_RNDN);
            Expr* we = expr_new_mpfr_bits(ob);
            mpfr_set(we->data.mpfr, wv, MPFR_RNDN);
            rows[pos] = gq_pair(xe, we);
        }
    }
    if (n % 2 == 1) {
        mpfr_set_ui(z, 0, MPFR_RNDN);
        gq_legendre_mp(n, z, pn, pn1, t0, t1, t2);
        mpfr_mul_ui(dp, pn1, (unsigned long)n, MPFR_RNDN);
        mpfr_sqr(wv, dp, MPFR_RNDN);
        mpfr_ui_div(wv, 2, wv, MPFR_RNDN);
        mpfr_mul(wv, wv, h, MPFR_RNDN);
        Expr* xe = expr_new_mpfr_bits(ob);
        mpfr_set(xe->data.mpfr, c, MPFR_RNDN);
        Expr* we = expr_new_mpfr_bits(ob);
        mpfr_set(we->data.mpfr, wv, MPFR_RNDN);
        rows[half] = gq_pair(xe, we);
    }

    /* A reversed interval (b < a) maps the increasing [-1, 1] nodes to
     * decreasing ones: restore increasing order. */
    if (mpfr_cmp(b, a) < 0)
        for (int i = 0, j = n - 1; i < j; i++, j--) { Expr* t = rows[i]; rows[i] = rows[j]; rows[j] = t; }

    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    mpfr_clears(z, pn, pn1, dp, omz2, t0, t1, t2, dz, h, c, xv, wv, a, b, im, (mpfr_ptr)0);
    free(rows); free(xd); free(wd);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ *
 *  Builtins                                                           *
 * ------------------------------------------------------------------ */

static Expr* builtin_gaussianquadratureweights(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 3 || argc > 4) return builtin_arg_error(GQ_HEAD, argc, 3, 4);
    Expr** args = res->data.function.args;
    int n;
    if (!gq_order(args[0], &n)) {
        mth_message(GQ_HEAD, "intpm",
                    "Positive machine-sized integer expected at position 1.");
        return NULL;
    }
    bool machine = true;
    double digits = 0.0;
    if (argc == 4) {
        if (!gq_precision(args[3], &digits, &machine)) {
            mth_message(GQ_HEAD, "prec",
                        "The precision at position 4 should be a positive number.");
            return NULL;
        }
    }
    if (machine) {
        double a, b;
        if (!gq_real_machine(args[1], &a) || !gq_real_machine(args[2], &b)) return NULL;
        Expr* r = gq_machine_expr(n, a, b);
        if (r && b < a) {   /* keep increasing order on a reversed interval */
            Expr** rows = r->data.function.args;
            for (int i = 0, j = n - 1; i < j; i++, j--) { Expr* t = rows[i]; rows[i] = rows[j]; rows[j] = t; }
        }
        return r;
    }
#ifdef USE_MPFR
    return gq_mpfr_expr(n, args[1], args[2], digits);
#else
    return NULL;
#endif
}

/* (b - a)^(2n+1) (n!)^4 / ((2n+1) ((2n)!)^3) Derivative[2n][f], built
 * symbolically and evaluated. */
static Expr* builtin_gaussianquadratureerror(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 4) return builtin_arg_error(GQE_HEAD, argc, 4, 4);
    Expr** args = res->data.function.args;
    int n;
    if (!gq_order(args[0], &n) || n > 100000) {
        mth_message(GQE_HEAD, "intpm",
                    "Positive machine-sized integer expected at position 1.");
        return NULL;
    }
    /* coefficient (n!)^4 / ((2n+1) ((2n)!)^3), exact */
    Expr* nn = expr_new_integer(n);
    Expr* n2 = expr_new_integer(2 * (int64_t)n);
    Expr* fact_n  = expr_new_function(expr_new_symbol(SYM_Factorial), (Expr*[]){ nn }, 1);
    Expr* fact_2n = expr_new_function(expr_new_symbol(SYM_Factorial), (Expr*[]){ n2 }, 1);
    Expr* num = expr_new_function(expr_new_symbol(SYM_Power),
                                  (Expr*[]){ fact_n, expr_new_integer(4) }, 2);
    Expr* den = expr_new_function(expr_new_symbol(SYM_Power),
                                  (Expr*[]){ fact_2n, expr_new_integer(-3) }, 2);
    Expr* inv = expr_new_function(expr_new_symbol(SYM_Power),
                                  (Expr*[]){ expr_new_integer(2 * (int64_t)n + 1),
                                             expr_new_integer(-1) }, 2);
    Expr* len = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(args[3]),
                               expr_new_function(expr_new_symbol(SYM_Times),
                                   (Expr*[]){ expr_new_integer(-1), expr_copy(args[2]) }, 2) }, 2);
    Expr* lenp = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ len, expr_new_integer(2 * (int64_t)n + 1) }, 2);
    Expr* dhead = expr_new_function(expr_new_symbol(SYM_Derivative),
                    (Expr*[]){ expr_new_integer(2 * (int64_t)n) }, 1);
    Expr* deriv = expr_new_function(dhead, (Expr*[]){ expr_copy(args[1]) }, 1);
    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ num, den, inv, lenp, deriv }, 5);
    return eval_and_free(prod);
}

void gaussquad_init(void) {
    symtab_add_builtin(GQ_HEAD, builtin_gaussianquadratureweights);
    symtab_get_def(GQ_HEAD)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(GQ_HEAD,
        "GaussianQuadratureWeights[n, a, b] gives the list of {xi, wi} pairs of the "
        "n-point Gauss-Legendre rule on [a, b] at machine precision, nodes increasing; "
        "Sum[wi f[xi]] integrates every polynomial of degree <= 2n - 1 exactly.\n"
        "GaussianQuadratureWeights[n, a, b, prec] gives the rule with every node and "
        "weight correct to prec digits.");

    symtab_add_builtin(GQE_HEAD, builtin_gaussianquadratureerror);
    symtab_get_def(GQE_HEAD)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(GQE_HEAD,
        "GaussianQuadratureError[n, f, a, b] gives the error Integrate[f, {x, a, b}] - "
        "Sum[wi f[xi]] of the n-point Gauss-Legendre rule, "
        "(b - a)^(2n+1) (n!)^4/((2n + 1) ((2n)!)^3) f^(2n), the derivative taken at an "
        "unspecified point of (a, b).");
}
