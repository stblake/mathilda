/* Mathilda -- RiemannSiegelZ and RiemannSiegelTheta (see riemann_siegel.h).
 *
 *   theta(t) = (LogGamma(1/4 + i t/2) - LogGamma(1/4 - i t/2)) / (2 i) - t Log(Pi)/2
 *            = Im LogGamma(1/4 + i t/2) - t Log(Pi)/2           (real t)
 *   Z(t)     = Exp(i theta(t)) Zeta(1/2 + i t)
 *
 * theta is the continuous phase that makes Z real on the real line, so the
 * zeros of Zeta on the critical line are the real zeros of Z (the 10th is at
 * t = 49.7738324776723021819...).  theta is odd, Z even.
 *
 * Evaluation:
 *   t exact 0                      theta = 0, Z = Zeta[1/2]
 *   other exact t                  stays symbolic
 *   machine real / complex t       sf_machine_riemann_siegel_*: arb's
 *                                  acb_dirichlet_hardy_z / _theta at a 53-bit
 *                                  target behind the precision ladder of
 *                                  flint_num_bridge.c, rounded to doubles -- the
 *                                  helper the NDArray kernels and Compile[] call
 *   arbitrary precision            the same arb kernels at the argument's
 *                                  precision (flint_num_hardy)
 *   no FLINT                       the defining expressions over LogGamma and
 *                                  Zeta, evaluated numerically
 * A real t answers a real number (the imaginary part of Z is zero by
 * construction; only the real part of the ball is kept).
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "riemann_siegel.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "arithmetic.h"          /* is_complex, make_complex */
#include "attr.h"
#include "common.h"              /* builtin_arg_error, common_machine_real_value */
#include "eval.h"                /* eval_and_free */
#include "expr.h"
#include "flint_num_bridge.h"    /* flint_num_hardy, flint_hardy_d */
#include "numeric.h"             /* numeric_min_inexact_bits */
#include "symtab.h"
#include "../sym_names.h"

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ------------------------------------------------------------------ */
/* Machine helpers                                                     */
/* ------------------------------------------------------------------ */

bool sf_machine_riemann_siegel_z(double t, double* out) {
    double ri;
    return flint_hardy_d(0, t, 0.0, out, &ri);
}
bool sf_machine_riemann_siegel_z_c(double tr, double ti, double* rr, double* ri) {
    if (ti == 0.0) { *ri = 0.0; return sf_machine_riemann_siegel_z(tr, rr); }
    return flint_hardy_d(0, tr, ti, rr, ri);
}
bool sf_machine_riemann_siegel_theta(double t, double* out) {
    double ri;
    return flint_hardy_d(1, t, 0.0, out, &ri);
}
bool sf_machine_riemann_siegel_theta_c(double tr, double ti, double* rr, double* ri) {
    if (ti == 0.0) { *ri = 0.0; return sf_machine_riemann_siegel_theta(tr, rr); }
    return flint_hardy_d(1, tr, ti, rr, ri);
}

/* ------------------------------------------------------------------ */
/* Expression helpers                                                  */
/* ------------------------------------------------------------------ */

static Expr* mk_fn1(const char* h, Expr* a) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}
static Expr* mk_rat(int64_t p, int64_t q) {
    return mk_fn2("Rational", expr_new_integer(p), expr_new_integer(q));
}

/* theta(t) as an expression (the no-FLINT numeric route). */
static Expr* rs_theta_expr(const Expr* t) {
    /* 1/4 +- i t/2 */
    Expr* up = mk_fn2("Plus", mk_rat(1, 4),
                      expr_new_function(expr_new_symbol("Times"), (Expr*[]){
                          mk_rat(1, 2), expr_new_symbol("I"), expr_copy((Expr*)t) }, 3));
    Expr* dn = mk_fn2("Plus", mk_rat(1, 4),
                      expr_new_function(expr_new_symbol("Times"), (Expr*[]){
                          mk_rat(-1, 2), expr_new_symbol("I"), expr_copy((Expr*)t) }, 3));
    Expr* diff = mk_fn2("Plus", mk_fn1("LogGamma", up),
                        mk_fn2("Times", expr_new_integer(-1), mk_fn1("LogGamma", dn)));
    Expr* part = expr_new_function(expr_new_symbol("Times"), (Expr*[]){
                     mk_rat(-1, 2), expr_new_symbol("I"), diff }, 3);      /* /(2i) */
    Expr* lp = expr_new_function(expr_new_symbol("Times"), (Expr*[]){
                   mk_rat(-1, 2), expr_copy((Expr*)t), mk_fn1("Log", expr_new_symbol("Pi")) }, 3);
    return mk_fn2("Plus", part, lp);
}

/* Z(t) as an expression (the no-FLINT numeric route). */
static Expr* rs_z_expr(const Expr* t) {
    Expr* ph = mk_fn1("Exp", mk_fn2("Times", expr_new_symbol("I"), rs_theta_expr(t)));
    Expr* s = mk_fn2("Plus", mk_rat(1, 2), mk_fn2("Times", expr_new_symbol("I"),
                                                  expr_copy((Expr*)t)));
    return mk_fn2("Times", ph, mk_fn1("Zeta", s));
}

/* True when e is a real numeric leaf (no Complex head). */
static bool rs_is_real_number(const Expr* e) {
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) return false;
    return true;
}

/* For a real t keep only the real part of a numeric result (owned). */
static Expr* rs_real_part(Expr* v) {
    Expr *re, *im;
    if (v && is_complex(v, &re, &im)) {
        Expr* r = expr_copy(re);
        expr_free(v);
        return r;
    }
    return v;
}

/* A numeric result demoted to machine Reals (53-bit MPFR leaves from arb). */
static Expr* rs_demote(Expr* v) {
#ifdef USE_MPFR
    if (!v) return NULL;
    if (v->type == EXPR_MPFR) {
        double d = mpfr_get_d(v->data.mpfr, MPFR_RNDN);
        if (!isfinite(d)) return v;
        expr_free(v);
        return expr_new_real(d);
    }
    Expr *re, *im;
    if (is_complex(v, &re, &im) && re->type == EXPR_MPFR && im->type == EXPR_MPFR) {
        double dr = mpfr_get_d(re->data.mpfr, MPFR_RNDN);
        double di = mpfr_get_d(im->data.mpfr, MPFR_RNDN);
        if (!isfinite(dr) || !isfinite(di)) return v;
        expr_free(v);
        return make_complex(expr_new_real(dr), expr_new_real(di));
    }
#endif
    return v;
}

/* Minimum precision (bits) of the inexact parts of e; 0 when e is exact. */
static long rs_inexact_bits(const Expr* e) {
#ifdef USE_MPFR
    return numeric_min_inexact_bits(e);
#else
    Expr *re, *im;
    if (e && e->type == EXPR_REAL) return 53;
    if (e && is_complex((Expr*)e, &re, &im)
            && (re->type == EXPR_REAL || im->type == EXPR_REAL)) return 53;
    return 0;
#endif
}

/* ------------------------------------------------------------------ */
/* Driver                                                              */
/* ------------------------------------------------------------------ */

static Expr* rs_eval(Expr* res, int theta) {
    const char* name = theta ? "RiemannSiegelTheta" : "RiemannSiegelZ";
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error(name, argc, 1, 1);
    Expr* t = res->data.function.args[0];

    if (t->type == EXPR_INTEGER && t->data.integer == 0)
        return theta ? expr_new_integer(0)
                     : eval_and_free(mk_fn1("Zeta", mk_rat(1, 2)));

    long bits = rs_inexact_bits(t);
    if (bits <= 0 || !expr_is_numeric_quantity(t)) return NULL;   /* exact: symbolic */
    bool real = rs_is_real_number(t);

    if (bits <= 53) {
        double tr, ti = 0.0, rr, ri;
        Expr *re, *im;
        bool machine = real ? common_machine_real_value(t, &tr)
                            : (is_complex(t, &re, &im) && common_machine_real_value(re, &tr)
                               && common_machine_real_value(im, &ti));
        if (machine) {
            bool ok = theta ? sf_machine_riemann_siegel_theta_c(tr, ti, &rr, &ri)
                            : sf_machine_riemann_siegel_z_c(tr, ti, &rr, &ri);
            if (ok) return real ? expr_new_real(rr)
                                : make_complex(expr_new_real(rr), expr_new_real(ri));
        }
    }

    Expr* v = flint_num_hardy(t, theta);
    if (!v) {
        /* No FLINT: the defining expression, evaluated numerically. */
        Expr* e = theta ? rs_theta_expr(t) : rs_z_expr(t);
        v = eval_and_free(mk_fn1("N", e));
        if (bits > 53) {
            expr_free(v);
            e = theta ? rs_theta_expr(t) : rs_z_expr(t);
            v = eval_and_free(mk_fn2("N", e, expr_new_real(numeric_bits_to_digits(bits))));
        }
        if (!v || !expr_is_numeric_quantity(v)) { expr_free(v); return NULL; }
    }
    if (real) v = rs_real_part(v);
    if (bits <= 53) v = rs_demote(v);
    return v;
}

Expr* builtin_riemann_siegel_z(Expr* res)     { return rs_eval(res, 0); }
Expr* builtin_riemann_siegel_theta(Expr* res) { return rs_eval(res, 1); }

void riemann_siegel_init(void) {
    symtab_add_builtin("RiemannSiegelZ", builtin_riemann_siegel_z);
    symtab_get_def("RiemannSiegelZ")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_set_docstring("RiemannSiegelZ",
        "RiemannSiegelZ[t]\n"
        "Hardy's function Z(t) = Exp[I RiemannSiegelTheta[t]] Zeta[1/2 + I t]: real\n"
        "for real t, even, with |Z(t)| = |Zeta[1/2 + I t]|, so its real zeros are the\n"
        "zeros of Zeta on the critical line. RiemannSiegelZ[0] = Zeta[1/2].\n"
        "Machine, arbitrary-precision and complex arguments; fast paths over packed\n"
        "arrays and in Compile[].");
    symtab_add_builtin("RiemannSiegelTheta", builtin_riemann_siegel_theta);
    symtab_get_def("RiemannSiegelTheta")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_set_docstring("RiemannSiegelTheta",
        "RiemannSiegelTheta[t]\n"
        "The Riemann-Siegel theta function: Im LogGamma[1/4 + I t/2] - t Log[Pi]/2\n"
        "for real t (analytically continued to complex t), the phase that makes\n"
        "RiemannSiegelZ[t] = Exp[I theta(t)] Zeta[1/2 + I t] real. Odd;\n"
        "RiemannSiegelTheta[0] = 0. Machine, arbitrary-precision and complex\n"
        "arguments; fast paths over packed arrays and in Compile[].");
}
