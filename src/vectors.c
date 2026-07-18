/* UnitVector — the n-dimensional unit vector in the k-th direction.
 *
 *   UnitVector[k]       the 2-D unit vector in the k-th direction
 *                       (equivalent to UnitVector[2, k]).
 *   UnitVector[n, k]    the n-D unit vector: a length-n list with a 1 in
 *                       position k and 0s in every other position.
 *
 * Components are exact integers by default (WorkingPrecision -> Infinity).
 * The WorkingPrecision option selects the component representation, mirroring
 * HilbertMatrix (src/linalg/hilbertmat.c):
 *
 *   WorkingPrecision -> Infinity          exact integers (default)
 *   WorkingPrecision -> MachinePrecision  machine-precision Reals
 *   WorkingPrecision -> d                  d-digit MPFR Reals (d above machine
 *                                          precision; otherwise machine Reals)
 *
 * Diagnostics mirror Wolfram's surface text:
 *   - zero arguments               -> UnitVector::argt  (1 or 2 expected)
 *   - non-option trailing argument -> UnitVector::nonopt
 * Non-integer or out-of-range (k < 1 or k > n) arguments leave the call
 * unevaluated (return NULL), matching Mathilda's "can't evaluate" convention.
 */

#include "vectors.h"
#include "symtab.h"
#include "attr.h"
#include "common.h"
#include "arithmetic.h"
#include "numeric.h"
#include "sym_names.h"
#include "print.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* Representation chosen for the vector components. */
typedef enum {
    UV_PREC_EXACT,    /* exact integers (WorkingPrecision -> Infinity) */
    UV_PREC_MACHINE,  /* machine-precision Reals */
    UV_PREC_MPFR      /* arbitrary-precision MPFR Reals */
} uv_prec_mode;

/* Recognise an options-style rule `sym -> value` / `sym :> value`. */
static bool uv_is_rule(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           (e->data.function.head->data.symbol.name == SYM_Rule ||
            e->data.function.head->data.symbol.name == SYM_RuleDelayed) &&
           e->data.function.arg_count == 2 &&
           e->data.function.args[0]->type == EXPR_SYMBOL;
}

/* Interpret a WorkingPrecision value.  Returns true and sets the mode/bits
 * outputs for Infinity, MachinePrecision, or a positive numeric digit count;
 * returns false for any other value (the caller keeps the prior setting,
 * matching "last valid setting wins"). */
static bool uv_parse_working_precision(const Expr* val,
                                       uv_prec_mode* mode, long* bits) {
    if (val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_Infinity) {
        *mode = UV_PREC_EXACT;
        *bits = 0;
        return true;
    }
    if (val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_MachinePrecision) {
        *mode = UV_PREC_MACHINE;
        *bits = 0;
        return true;
    }

    double digits = 0.0;
    int64_t rn, rd;
    if (val->type == EXPR_INTEGER)       digits = (double)val->data.integer;
    else if (val->type == EXPR_REAL)     digits = val->data.real;
    else if (is_rational(val, &rn, &rd)) digits = (double)rn / (double)rd;
    else return false;
    if (digits <= 0.0) return false;

#ifdef USE_MPFR
    if (digits <= NUMERIC_MACHINE_PRECISION_DIGITS) {
        /* At or below machine precision the doubles path is exact enough. */
        *mode = UV_PREC_MACHINE;
        *bits = 0;
    } else {
        *mode = UV_PREC_MPFR;
        *bits = numeric_digits_to_bits(digits);
    }
#else
    /* Without MPFR a digit count degrades to machine precision. */
    static bool warned = false;
    if (!warned) {
        fprintf(stderr,
                "UnitVector::wprec: arbitrary precision unavailable "
                "(USE_MPFR=0); using machine precision.\n");
        warned = true;
    }
    *mode = UV_PREC_MACHINE;
    *bits = 0;
#endif
    return true;
}

/* Build a single component: `one` selects the value 1, otherwise 0. */
static Expr* uv_component(bool one, uv_prec_mode mode, long bits) {
    switch (mode) {
        case UV_PREC_MACHINE:
            return expr_new_real(one ? 1.0 : 0.0);
#ifdef USE_MPFR
        case UV_PREC_MPFR: {
            mpfr_t out;
            mpfr_init2(out, bits);
            mpfr_set_ui(out, one ? 1u : 0u, MPFR_RNDN);
            return expr_new_mpfr_move(out);
        }
#else
        case UV_PREC_MPFR: /* unreachable: parser never selects MPFR here */
            return expr_new_real(one ? 1.0 : 0.0);
#endif
        case UV_PREC_EXACT:
        default:
            return expr_new_integer(one ? 1 : 0);
    }
}

Expr* builtin_unit_vector(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 0) {
        return builtin_arg_error("UnitVector", argc, 1, 2);
    }

    Expr** args = res->data.function.args;

    /* Count leading required (non-option) arguments, capped at 2: the required
     * forms are UnitVector[k] and UnitVector[n, k]. */
    size_t nreq = 0;
    while (nreq < argc && nreq < 2 && !uv_is_rule(args[nreq])) nreq++;

    /* Every argument beyond the required ones must be an option rule; report
     * the rightmost offender via ::nonopt (Wolfram semantics). */
    Expr* last_bad = NULL;
    for (size_t i = nreq; i < argc; i++) {
        if (!uv_is_rule(args[i])) last_bad = args[i];
    }
    if (last_bad != NULL) {
        char* bad_str = expr_to_string(last_bad);
        char* call_str = expr_to_string(res);
        fprintf(stderr,
                "UnitVector::nonopt: Options expected (instead of %s) beyond "
                "position %zu in %s. An option must be a rule or a list of "
                "rules.\n",
                bad_str ? bad_str : "?", nreq, call_str ? call_str : "?");
        free(bad_str);
        free(call_str);
        return NULL;
    }

    /* Map required arguments to (n, k). */
    Expr* n_expr;
    Expr* k_expr;
    if (nreq == 1) {
        n_expr = NULL;   /* two-dimensional by default */
        k_expr = args[0];
    } else { /* nreq == 2 */
        n_expr = args[0];
        k_expr = args[1];
    }

    /* n and k must be positive machine integers, with 1 <= k <= n; otherwise
     * leave the call unevaluated (symbolic arguments flow through unchanged). */
    if (k_expr->type != EXPR_INTEGER || k_expr->data.integer < 1) return NULL;
    int64_t k = k_expr->data.integer;

    int64_t n;
    if (n_expr == NULL) {
        n = 2;
    } else {
        if (n_expr->type != EXPR_INTEGER || n_expr->data.integer < 1) return NULL;
        n = n_expr->data.integer;
    }
    if (k > n) return NULL;

    /* Parse trailing options.  Only WorkingPrecision is recognised; last valid
     * setting wins, an unparseable value is ignored. */
    uv_prec_mode mode = UV_PREC_EXACT;
    long bits = 0;
    for (size_t i = nreq; i < argc; i++) {
        Expr* opt = args[i];
        const char* name = opt->data.function.args[0]->data.symbol.name;
        if (name == SYM_WorkingPrecision) {
            uv_parse_working_precision(opt->data.function.args[1], &mode, &bits);
        }
    }

    /* Build the length-n list: 1 at position k, 0 elsewhere. */
    Expr** cells = malloc(sizeof(Expr*) * (size_t)n);
    for (int64_t i = 0; i < n; i++) {
        cells[i] = uv_component(i == k - 1, mode, bits);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), cells, (size_t)n);
    free(cells);
    return result;
}

void vectors_init(void) {
    symtab_add_builtin("UnitVector", builtin_unit_vector);
    symtab_get_def("UnitVector")->attributes |= ATTR_PROTECTED;
}
