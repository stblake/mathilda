/* Mathilda -- the generalized exponential integral E_n.
 *
 *   ExpIntegralE[n, z]   E_n(z) = Integrate[E^(-z t)/t^n, {t, 1, Infinity}]
 *
 * Layered evaluation (each argument kind takes the cheapest correct route):
 *   E_0(z) = e^(-z)/z                      exact, any z
 *   E_n(0) = 1/(n-1),  integer n >= 2      exact
 *   integer n >= 0, machine real x > 0     sf_machine_expn (sf_machine_expint.c):
 *                                          recurrence/series (x <= 1) / continued
 *                                          fraction (x > 1)
 *   integer n >= 0, machine complex z      sf_machine_expn_c (and real x < 0,
 *   (or real x < 0)                        whose value is complex)
 *   any other numeric call with an         arb's acb_hypgeom_expint, at the call's
 *   inexact argument (real order, MPFR)    precision
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
#include <stdbool.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

#include "sf_machine_expint.h"   /* sf_machine_expn */
#include "sf_scalar.h"           /* the shared machine-scalar gate */
#include "flint_num_bridge.h"    /* E_s(z) for any other inexact call */

/* Small expr builders (local, matching the deriv.c idiom). */
static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_fn1(const char* name, Expr* a) {
    return expr_new_function(expr_new_symbol(name), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(name), (Expr*[]){ a, b }, 2);
}

/* E_n(x) for an integer n >= 0 and real x > 0: the shared machine kernel
 * (sf_machine_expint.c), the same one the NDArray kernel and Compile[] call.
 * Returns false (defer) outside that domain. */
bool expintegrale_machine(long n, double x, double* out) {
    return sf_machine_expn((double)n, x, out);
}

Expr* builtin_expintegrale(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* n = res->data.function.args[0];
    Expr* z = res->data.function.args[1];

    /* Machine-precision call with an integer order (exact, or a machine real
     * that is integer-valued) and a real x > 0: the double kernel, before any
     * expression building. Anything else falls through unchanged. */
    {
        SFMachineArgs ma;
        if (sf_machine_scalar_args(res, &ma) && !ma.is_complex[0]) {
            double v, vi;
            if (!ma.is_complex[1]) {
                if (sf_machine_expn(ma.re[0], ma.re[1], &v)) return expr_new_real(v);
                /* x < 0: E_n is complex there (the upper side of the cut) */
                if (ma.re[1] < 0.0 &&
                    sf_machine_expn_c(ma.re[0], ma.re[1], 0.0, &v, &vi))
                    return sf_machine_complex_result(v, vi);
            } else if (sf_machine_expn_c(ma.re[0], ma.re[1], ma.im[1], &v, &vi)) {
                return sf_machine_complex_result(v, vi);
            }
        }
    }

    /* E_0(z) = e^(-z)/z  (exact; the derivative base case). */
    if (n->type == EXPR_INTEGER && n->data.integer == 0)
        return mk_fn2("Times",
            mk_fn1("Exp", mk_fn2("Times", mk_int(-1), expr_copy(z))),
            mk_fn2("Power", expr_copy(z), mk_int(-1)));

    /* E_n(0) diverges for integer n <= 1 (the integrand ~ 1/t^n near 0, or
     * E_n = Gamma[1-n] z^(n-1) + ... for n < 1): ComplexInfinity, as the
     * inexact ExpIntegralE[1, 0.] (audit ExpIntegralE#2). */
    if (n->type == EXPR_INTEGER && n->data.integer <= 1 &&
        z->type == EXPR_INTEGER && z->data.integer == 0)
        return expr_new_symbol(SYM_ComplexInfinity);

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

    /* Any other numeric call with an inexact argument -- a non-integer order,
     * arbitrary precision: E_s(z) through arb (principal branch). */
    if (flint_num_call_inexact(res)) {
        /* E_s(0) diverges for Re s <= 1 (the integrand ~ 1/t^s): arb's NaN
         * there would leave the call unevaluated, the reference says
         * ComplexInfinity. */
        double zr = 1.0, sr = 2.0;
        bool z_zero = (z->type == EXPR_REAL && z->data.real == 0.0)
#ifdef USE_MPFR
                   || (z->type == EXPR_MPFR && mpfr_zero_p(z->data.mpfr))
#endif
                   ;
        if (z_zero) zr = 0.0;
        if (n->type == EXPR_REAL) sr = n->data.real;
        else if (n->type == EXPR_INTEGER) sr = (double)n->data.integer;
#ifdef USE_MPFR
        else if (n->type == EXPR_MPFR) sr = mpfr_get_d(n->data.mpfr, MPFR_RNDN);
#endif
        if (zr == 0.0 && sr <= 1.0) return expr_new_symbol(SYM_ComplexInfinity);
        return flint_num_for_call(flint_num_expint(n, z), res);
    }
    return NULL;
}

void expintegrale_init(void) {
    symtab_add_builtin("ExpIntegralE", builtin_expintegrale);
    symtab_get_def("ExpIntegralE")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
