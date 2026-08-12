/* Mathilda — shared AccuracyGoal / PrecisionGoal handling. See nc_accuracy.h
 * for the design and the combined-tolerance contract. */
#include "nc_accuracy.h"

#include <stdio.h>
#include <math.h>
#include <float.h>          /* DBL_MANT_DIG, via NUMERIC_MACHINE_PRECISION_DIGITS */

#include "sym_names.h"
#include "numeric.h"        /* NUMERIC_MACHINE_PRECISION_DIGITS */
#include "arithmetic.h"     /* is_rational */

/* Coerce a real-valued numeric leaf to a double. Mirrors the per-module
 * `*_to_double_real` helpers (e.g. nint.c) so the goal parser accepts the same
 * value forms option values arrive in: Integer, Real, BigInt, MPFR, Rational. */
static bool nc_to_double(const Expr* v, double* out) {
    if (!v) return false;
    int64_t rn, rd;
    switch (v->type) {
        case EXPR_INTEGER: *out = (double)v->data.integer;   return true;
        case EXPR_REAL:    *out = v->data.real;              return true;
        case EXPR_BIGINT:  *out = mpz_get_d(v->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(v->data.mpfr, MPFR_RNDN); return true;
#endif
        case EXPR_FUNCTION:
            if (is_rational(v, &rn, &rd)) { *out = (double)rn / (double)rd; return true; }
            return false;
        default: return false;
    }
}

bool nc_parse_goal(const Expr* v, double* digits_out) {
    if (!v || !digits_out) return false;
    if (v->type == EXPR_SYMBOL) {
        const char* n = v->data.symbol.name;
        if (n == SYM_Automatic)        { *digits_out = NC_GOAL_AUTO;    return true; }
        if (n == SYM_Infinity)         { *digits_out = NC_GOAL_OFF;     return true; }
        if (n == SYM_MachinePrecision) { *digits_out = NC_GOAL_MACHINE; return true; }
        return false;
    }
    double d;
    if (nc_to_double(v, &d) && d > 0.0) { *digits_out = d; return true; }
    return false;
}

double nc_resolve_goal(double goal_digits, double wp_digits) {
    if (goal_digits < 0.0) {                 /* NC_GOAL_AUTO or NC_GOAL_MACHINE */
        /* Automatic seeks near-full working precision: two guard digits below
         * WorkingPrecision. This matches the historical NSum/NResidue default
         * and the intent that "Automatic" maximise accuracy — a looser value
         * (e.g. WorkingPrecision/2) would dominate the combined tolerance and
         * silently cap the result even when an explicit AccuracyGoal is tight. */
        double w = (wp_digits > 0.0) ? wp_digits : NUMERIC_MACHINE_PRECISION_DIGITS;
        double a = w - 2.0;
        return (a > 1.0) ? a : 1.0;
    }
    return goal_digits;                      /* explicit number, or +Inf (Off) */
}

double nc_combined_tol(double acc_digits, double prec_digits,
                       double result_mag, double wp_digits) {
    double a = nc_resolve_goal(acc_digits,  wp_digits);
    double p = nc_resolve_goal(prec_digits, wp_digits);
    double abstol = pow(10.0, -a);           /* a == +Inf -> 0 */
    double reltol = pow(10.0, -p);           /* p == +Inf -> 0 */
    double tol = abstol + fabs(result_mag) * reltol;
    if (!(tol > 0.0)) tol = 0.0;             /* NaN/underflow guard */
    return tol;
}

void nc_warn_goal(const char* head, double achieved_err, double tol) {
    fprintf(stderr,
            "%s::accgl: requested accuracy goal not met (estimated error %.3g "
            "exceeds tolerance %.3g); returning the best approximation obtained.\n",
            head ? head : "N", achieved_err, tol);
}
