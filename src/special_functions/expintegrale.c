/* Mathilda -- the generalized exponential integral E_n.
 *
 *   ExpIntegralE[n, z]   E_n(z) = Integrate[E^(-z t)/t^n, {t, 1, Infinity}]
 *
 * Layered evaluation (each argument kind takes the cheapest correct route):
 *   E_0(z) = e^(-z)/z                      exact, any z
 *   E_n(0) = 1/(n-1),  integer n >= 2      exact
 *   integer n >= 1, inexact real x > 0     machine E_n via CF (x>1) / series (x<1)
 *   everything else                        stays symbolic (return NULL)
 *
 * The derivative rule d/dz E_n(z) = -E_(n-1)(z) (with the base E_0 above) lives
 * in deriv.c; the half-line Mellin transform of E_n is produced by the
 * operational-calculus fallback in integrate_ramanujan.c via that chain.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "expintegrale.h"

#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "sym_names.h"
#include "common.h"

#include <math.h>
#include <float.h>
#include <stdbool.h>

/* Small expr builders (local, matching the deriv.c idiom). */
static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_fn1(const char* name, Expr* a) {
    return expr_new_function(expr_new_symbol(name), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(name), (Expr*[]){ a, b }, 2);
}

/* E_n(x) for integer n >= 0 and real x > 0, by the standard continued-fraction
 * (x > 1) / power-series (0 < x <= 1) split.  Returns false outside that domain
 * (x <= 0, which is on/left of the branch point, is deferred to stay symbolic).
 * Verified against NIntegrate in the unit tests, not by a recalled constant. */
bool expintegrale_machine(long n, double x, double* out) {
    static const double EULER_GAMMA = 0.5772156649015328606065;
    const double FPMIN = DBL_MIN / DBL_EPSILON;
    const int    MAXIT = 400;
    const double EPS   = 10.0 * DBL_EPSILON;
    if (n < 0 || x <= 0.0) return false;
    if (n == 0) { *out = exp(-x) / x; return true; }
    long nm1 = n - 1;
    if (x > 1.0) {                              /* Lentz continued fraction */
        double b = x + (double)n;
        double c = 1.0 / FPMIN;
        double d = 1.0 / b;
        double h = d;
        for (int i = 1; i <= MAXIT; i++) {
            double a = -(double)i * (double)(nm1 + i);
            b += 2.0;
            d = 1.0 / (a * d + b);
            c = b + a / c;
            double del = c * d;
            h *= del;
            if (fabs(del - 1.0) <= EPS) { *out = h * exp(-x); return true; }
        }
        return false;
    } else {                                    /* power series */
        double ans = (nm1 != 0) ? 1.0 / (double)nm1 : -log(x) - EULER_GAMMA;
        double fact = 1.0;
        for (int i = 1; i <= MAXIT; i++) {
            fact *= -x / (double)i;
            double del;
            if (i != nm1) {
                del = -fact / (double)(i - nm1);
            } else {
                double psi = -EULER_GAMMA;
                for (long ii = 1; ii <= nm1; ii++) psi += 1.0 / (double)ii;
                del = fact * (-log(x) + psi);
            }
            ans += del;
            if (fabs(del) < fabs(ans) * EPS) { *out = ans; return true; }
        }
        return false;
    }
}

Expr* builtin_expintegrale(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* n = res->data.function.args[0];
    Expr* z = res->data.function.args[1];

    /* E_0(z) = e^(-z)/z  (exact; the derivative base case). */
    if (n->type == EXPR_INTEGER && n->data.integer == 0)
        return mk_fn2("Times",
            mk_fn1("Exp", mk_fn2("Times", mk_int(-1), expr_copy(z))),
            mk_fn2("Power", expr_copy(z), mk_int(-1)));

    /* E_n(0) = 1/(n-1) for integer n >= 2. */
    if (n->type == EXPR_INTEGER && n->data.integer >= 2 &&
        z->type == EXPR_INTEGER && z->data.integer == 0)
        return mk_fn2("Power", mk_int((int64_t)n->data.integer - 1), mk_int(-1));

    /* Machine numeric: integer order n >= 1, inexact real x > 0. */
    if (n->type == EXPR_INTEGER && n->data.integer >= 1 && z->type == EXPR_REAL) {
        double x = z->data.real, out;
        if (x > 0.0 && expintegrale_machine((long)n->data.integer, x, &out))
            return expr_new_real(out);
    }
    return NULL;
}

void expintegrale_init(void) {
    symtab_add_builtin("ExpIntegralE", builtin_expintegrale);
    symtab_get_def("ExpIntegralE")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
