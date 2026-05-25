/* Mathilda — Precision / Accuracy introspection and coercion.
 *
 * See precision.h for the module overview.
 *
 * Precision semantics:
 *   - Exact leaves (Integer, BigInt, Rational, symbols with no numeric
 *     value) → Infinity.
 *   - EXPR_REAL                                                 → MachinePrecision.
 *   - EXPR_MPFR (Phase 2)                                       → prec bits / log2(10).
 *   - Complex[a, b]                                             → min(Precision[a], Precision[b]).
 *   - Compound f[...] (generic)                                 → min of component precisions.
 *
 * Accuracy semantics (Mathematica-compatible):
 *     accuracy = precision − log10(|x|)          (for x ≠ 0)
 *   Exact values (Integer, BigInt, Rational, symbolic constants) → Infinity.
 *   Inexact zeros are finite:
 *     - Machine 0.0 → MachinePrecision − log10(DBL_MIN) ≈ 323.607.
 *     - MPFR 0 with bit-precision p → p / log2(10) digits (matches MMA's
 *       Accuracy[0``p] = p).
 */

#include "precision.h"

#include "arithmetic.h"
#include "attr.h"
#include "numeric.h"
#include "symtab.h"
#include "sym_names.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef USE_MPFR
#include <mpfr.h>
#endif

static const double LOG2_10 = 3.3219280948873626;

/* Sentinel doubles for the symbolic "precision" values. */
static Expr* make_infinity(void)         { return expr_new_symbol("Infinity"); }
static Expr* make_machineprecision(void) { return expr_new_symbol("MachinePrecision"); }

/* ------------------------------------------------------------------------
 *  Forward-declared helpers
 * ---------------------------------------------------------------------- */
static Expr* precision_of(const Expr* e);
static Expr* accuracy_of(const Expr* e);

/* True if an Expr is the symbol `Infinity` (not -Infinity, not ComplexInfinity). */
static bool expr_is_positive_infinity(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == SYM_Infinity;
}

/* Given two precision results (each Infinity, MachinePrecision, or a
 * Real), return the smaller. Used for combining component precisions. */
static Expr* precision_min(Expr* a, Expr* b) {
    if (expr_is_positive_infinity(a)) { expr_free(a); return b; }
    if (expr_is_positive_infinity(b)) { expr_free(b); return a; }
    /* MachinePrecision < any explicit MPFR digits > machine: we approximate
     * MachinePrecision as the constant ≈15.95 for the purpose of this
     * comparison, but keep the symbolic form in the return. */
    double da = (a->type == EXPR_SYMBOL && a->data.symbol == SYM_MachinePrecision)
                ? NUMERIC_MACHINE_PRECISION_DIGITS
                : (a->type == EXPR_REAL ? a->data.real
                                        : (a->type == EXPR_INTEGER ? (double)a->data.integer : 0.0));
    double db = (b->type == EXPR_SYMBOL && b->data.symbol == SYM_MachinePrecision)
                ? NUMERIC_MACHINE_PRECISION_DIGITS
                : (b->type == EXPR_REAL ? b->data.real
                                        : (b->type == EXPR_INTEGER ? (double)b->data.integer : 0.0));
    if (da <= db) { expr_free(b); return a; }
    expr_free(a); return b;
}

/* ------------------------------------------------------------------------
 *  precision_of — internal
 * ---------------------------------------------------------------------- */
static Expr* precision_of(const Expr* e) {
    if (!e) return make_infinity();

    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:
        case EXPR_STRING:
            return make_infinity();
        case EXPR_REAL:
            return make_machineprecision();
#ifdef USE_MPFR
        case EXPR_MPFR: {
            double digits = (double)mpfr_get_prec(e->data.mpfr) / LOG2_10;
            return expr_new_real(digits);
        }
#endif
        case EXPR_SYMBOL:
            /* Numeric constants (Pi, E, …) are exact when left symbolic. */
            return make_infinity();
        case EXPR_FUNCTION: {
            /* Rational is exact; Complex reduces to min of parts;
             * anything else — conservatively take the min across parts. */
            if (is_rational((Expr*)e, NULL, NULL)) return make_infinity();
            Expr *re, *im;
            if (is_complex((Expr*)e, &re, &im)) {
                return precision_min(precision_of(re), precision_of(im));
            }
            Expr* best = make_infinity();
            for (size_t i = 0; i < e->data.function.arg_count; ++i) {
                best = precision_min(best, precision_of(e->data.function.args[i]));
            }
            return best;
        }
    }
    return make_infinity();
}

/* ------------------------------------------------------------------------
 *  accuracy_of — internal
 *
 *  accuracy = precision − log10(|x|), computed pointwise. Compound
 *  structures take the min over components.
 * ---------------------------------------------------------------------- */
static double expr_log10_abs(const Expr* e) {
    if (!e) return 0.0;
    if (e->type == EXPR_INTEGER) {
        int64_t v = e->data.integer;
        if (v == 0) return 0.0;
        return log10((double)(v < 0 ? -v : v));
    }
    if (e->type == EXPR_REAL) {
        double v = fabs(e->data.real);
        if (v == 0.0) return 0.0;
        return log10(v);
    }
    if (e->type == EXPR_BIGINT) {
        double v = fabs(mpz_get_d(e->data.bigint));
        if (v == 0.0) return 0.0;
        return log10(v);
    }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        double v = fabs(mpfr_get_d(e->data.mpfr, MPFR_RNDN));
        if (v == 0.0) return 0.0;
        return log10(v);
    }
#endif
    int64_t n, d;
    if (is_rational((Expr*)e, &n, &d)) {
        double v = (double)n / (double)d;
        v = fabs(v);
        if (v == 0.0) return 0.0;
        return log10(v);
    }
    return 0.0;
}

/* Accuracy of an inexact zero. Mathematica treats approximate zeros as
 * having the largest accuracy that the representation can resolve:
 *   - Machine 0.0     → MachinePrecision + (-log10($MinMachineNumber))
 *                       = 53/log2(10) + (-log10(DBL_MIN)) ≈ 323.607.
 *   - MPFR 0 @ p bits → p / log2(10) digits, matching Accuracy[0``p] = p. */
static double accuracy_of_machine_zero(void) {
    return NUMERIC_MACHINE_PRECISION_DIGITS - log10(DBL_MIN);
}

static Expr* accuracy_of(const Expr* e) {
    if (!e) return make_infinity();

    /* Exact zero (Integer / BigInt 0, Rational 0/n) → Infinity.
     * Inexact zero (machine 0.0, MPFR 0) falls through to the per-type
     * finite-accuracy formulas below. */
    if (e->type == EXPR_INTEGER && e->data.integer == 0) return make_infinity();
    if (e->type == EXPR_BIGINT  && mpz_sgn(e->data.bigint) == 0) return make_infinity();

    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:
        case EXPR_STRING:
        case EXPR_SYMBOL:
            return make_infinity();
        case EXPR_REAL:
            if (e->data.real == 0.0) return expr_new_real(accuracy_of_machine_zero());
            return expr_new_real(NUMERIC_MACHINE_PRECISION_DIGITS - expr_log10_abs(e));
#ifdef USE_MPFR
        case EXPR_MPFR: {
            double prec_digits = (double)mpfr_get_prec(e->data.mpfr) / LOG2_10;
            if (mpfr_zero_p(e->data.mpfr)) return expr_new_real(prec_digits);
            return expr_new_real(prec_digits - expr_log10_abs(e));
        }
#endif
        case EXPR_FUNCTION: {
            /* Exact rationals with numerator 0 are exact zero → Infinity;
             * otherwise rationals are exact and also → Infinity. */
            if (is_rational((Expr*)e, NULL, NULL)) return make_infinity();
            Expr *re, *im;
            if (is_complex((Expr*)e, &re, &im)) {
                return precision_min(accuracy_of(re), accuracy_of(im));
            }
            Expr* best = make_infinity();
            for (size_t i = 0; i < e->data.function.arg_count; ++i) {
                best = precision_min(best, accuracy_of(e->data.function.args[i]));
            }
            return best;
        }
    }
    return make_infinity();
}

/* ------------------------------------------------------------------------
 *  Builtins
 * ---------------------------------------------------------------------- */
Expr* builtin_precision(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    return precision_of(res->data.function.args[0]);
}

Expr* builtin_accuracy(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    return accuracy_of(res->data.function.args[0]);
}

/* Parse the precision/accuracy argument n of SetPrecision/SetAccuracy into
 * a NumericSpec (or None for MachinePrecision). Returns true on success. */
static bool parse_prec_arg(const Expr* n, NumericSpec* out_spec) {
    if (n->type == EXPR_SYMBOL && n->data.symbol == SYM_MachinePrecision) {
        *out_spec = numeric_machine_spec();
        return true;
    }
    if (n->type == EXPR_SYMBOL && n->data.symbol == SYM_Infinity) {
        /* Infinity precision: keep values exact. Implemented as a signal
         * to the walker — for now, we simply return a machine spec and
         * rely on numericalize's default exact-preservation for atoms. */
        *out_spec = numeric_machine_spec();
        return true;
    }
    int64_t rn, rd;
    double digits = 0.0;
    if (n->type == EXPR_INTEGER) digits = (double)n->data.integer;
    else if (n->type == EXPR_REAL) digits = n->data.real;
    else if (is_rational((Expr*)n, &rn, &rd)) digits = (double)rn / (double)rd;
    else return false;
    if (digits <= 0.0) return false;

#ifdef USE_MPFR
    out_spec->mode = NUMERIC_MODE_MPFR;
    out_spec->bits = numeric_digits_to_bits(digits);
    return true;
#else
    static bool warned = false;
    if (!warned) {
        fprintf(stderr,
                "SetPrecision::prec: arbitrary precision unavailable "
                "(USE_MPFR=0); using machine precision.\n");
        warned = true;
    }
    *out_spec = numeric_machine_spec();
    return true;
#endif
}

Expr* builtin_set_precision(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    NumericSpec spec;
    if (!parse_prec_arg(res->data.function.args[1], &spec)) return NULL;
    return numericalize(res->data.function.args[0], spec);
}

Expr* builtin_set_accuracy(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;

    /* Convert accuracy argument into a precision spec via the expression's
     * magnitude, then drive numericalize. This isn't Mathematica's *exact*
     * significance-arithmetic semantics but it's the standard "accuracy n
     * means n digits after the decimal point" approximation.
     *
     * The precision in digits is: n + log10(|x|). */
    const Expr* value = res->data.function.args[0];
    const Expr* n     = res->data.function.args[1];

    int64_t rn, rd;
    double target_acc = 0.0;
    if (n->type == EXPR_INTEGER) target_acc = (double)n->data.integer;
    else if (n->type == EXPR_REAL) target_acc = n->data.real;
    else if (is_rational((Expr*)n, &rn, &rd)) target_acc = (double)rn / (double)rd;
    else if (n->type == EXPR_SYMBOL && n->data.symbol == SYM_MachinePrecision) {
        NumericSpec s = numeric_machine_spec();
        return numericalize(value, s);
    }
    else return NULL;

    if (target_acc <= 0.0) return NULL;

    double log10_abs = expr_log10_abs(value);
    double digits = target_acc + log10_abs;
    if (digits < 1.0) digits = 1.0;

    NumericSpec spec;
#ifdef USE_MPFR
    spec.mode = NUMERIC_MODE_MPFR;
    spec.bits = numeric_digits_to_bits(digits);
#else
    spec = numeric_machine_spec();
#endif
    return numericalize(value, spec);
}

/* ------------------------------------------------------------------------
 *  Registration
 * ---------------------------------------------------------------------- */
void precision_init(void) {
    symtab_add_builtin("Precision",    builtin_precision);
    symtab_add_builtin("Accuracy",     builtin_accuracy);
    symtab_add_builtin("SetPrecision", builtin_set_precision);
    symtab_add_builtin("SetAccuracy",  builtin_set_accuracy);

    symtab_get_def("Precision")   ->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_get_def("Accuracy")    ->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_get_def("SetPrecision")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_get_def("SetAccuracy") ->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;

    symtab_set_docstring("Precision",
        "Precision[x]\n"
        "\tReturns the number of decimal digits of precision in x.\n"
        "\tExact numbers return Infinity; machine-precision reals return\n"
        "\tthe symbol MachinePrecision; MPFR values return their declared\n"
        "\tprecision in decimal digits.");
    symtab_set_docstring("Accuracy",
        "Accuracy[x]\n"
        "\tReturns the number of digits of accuracy in x — equal to\n"
        "\tPrecision[x] − Log10[Abs[x]]. Exact numbers (including exact 0)\n"
        "\treturn Infinity. Inexact zeros are finite: machine 0. returns\n"
        "\t≈ 323.607; MPFR 0 of precision p returns p digits.");
    symtab_set_docstring("SetPrecision",
        "SetPrecision[x, n]\n"
        "\tReturns an expression equivalent to x with numeric values\n"
        "\tre-rounded or promoted to n decimal digits of precision.\n"
        "\tRequires a USE_MPFR build for n > MachinePrecision.");
    symtab_set_docstring("SetAccuracy",
        "SetAccuracy[x, n]\n"
        "\tReturns an expression equivalent to x with numeric values\n"
        "\tre-rounded or promoted to n decimal digits of accuracy.\n"
        "\tRequires a USE_MPFR build for high-accuracy outputs.");
}
