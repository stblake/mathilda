/* HilbertMatrix — the m x n Hilbert matrix with entries 1/(i + j - 1).
 *
 *   HilbertMatrix[n]        n x n Hilbert matrix.
 *   HilbertMatrix[{m, n}]   m x n Hilbert matrix.
 *
 * Entries are exact Rationals by default (WorkingPrecision -> Infinity).
 * The WorkingPrecision option selects the entry representation:
 *
 *   WorkingPrecision -> Infinity          exact Rationals (default)
 *   WorkingPrecision -> MachinePrecision  machine-precision Reals
 *   WorkingPrecision -> d                  d-digit MPFR Reals (d above
 *                                          machine precision; otherwise
 *                                          machine Reals, matching the
 *                                          rest of Mathilda's numeric
 *                                          tower).
 *
 * Diagnostics mirror Wolfram's surface text:
 *   - zero arguments               -> HilbertMatrix::argx
 *   - bad dimension specification  -> HilbertMatrix::dims
 *   - non-option trailing argument -> HilbertMatrix::nonopt
 */

#include "linalg.h"
#include "arithmetic.h"
#include "numeric.h"
#include "sym_names.h"
#include "print.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* Representation chosen for the matrix entries. */
typedef enum {
    HM_PREC_EXACT,    /* exact Rationals (WorkingPrecision -> Infinity) */
    HM_PREC_MACHINE,  /* machine-precision Reals */
    HM_PREC_MPFR      /* arbitrary-precision MPFR Reals */
} hm_prec_mode;

/* Recognise a `head[args...]` whose head is the symbol `sym`. */
static bool hm_is_call(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol == sym;
}

/* Parse a single positive machine integer (Integer > 0). */
static bool hm_positive_int(const Expr* e, int64_t* out) {
    if (e->type == EXPR_INTEGER && e->data.integer > 0) {
        *out = e->data.integer;
        return true;
    }
    return false;
}

/* Parse the dimension specification (args[0]): a positive integer n
 * (square) or a {m, n} pair of positive integers.  Returns true on a
 * valid spec, writing the row/column counts to *m and *n. */
static bool hm_parse_dims(const Expr* spec, int64_t* m, int64_t* n) {
    if (hm_positive_int(spec, m)) {
        *n = *m;
        return true;
    }
    if (hm_is_call(spec, SYM_List) && spec->data.function.arg_count == 2 &&
        hm_positive_int(spec->data.function.args[0], m) &&
        hm_positive_int(spec->data.function.args[1], n)) {
        return true;
    }
    return false;
}

/* Interpret a WorkingPrecision value.  Returns true and sets the
 * mode/bits outputs for Infinity, MachinePrecision, or a positive
 * numeric digit count;
 * returns false for any other value (the caller leaves the prior setting
 * in place, matching "last valid setting wins"). */
static bool hm_parse_working_precision(const Expr* val,
                                       hm_prec_mode* mode, long* bits) {
    if (val->type == EXPR_SYMBOL && val->data.symbol == SYM_Infinity) {
        *mode = HM_PREC_EXACT;
        *bits = 0;
        return true;
    }
    if (val->type == EXPR_SYMBOL && val->data.symbol == SYM_MachinePrecision) {
        *mode = HM_PREC_MACHINE;
        *bits = 0;
        return true;
    }

    double digits = 0.0;
    int64_t rn, rd;
    if (val->type == EXPR_INTEGER)      digits = (double)val->data.integer;
    else if (val->type == EXPR_REAL)    digits = val->data.real;
    else if (is_rational(val, &rn, &rd)) digits = (double)rn / (double)rd;
    else return false;
    if (digits <= 0.0) return false;

#ifdef USE_MPFR
    if (digits <= NUMERIC_MACHINE_PRECISION_DIGITS) {
        /* At or below machine precision the doubles path is exact enough. */
        *mode = HM_PREC_MACHINE;
        *bits = 0;
    } else {
        *mode = HM_PREC_MPFR;
        *bits = numeric_digits_to_bits(digits);
    }
#else
    /* Without MPFR a digit count degrades to machine precision. */
    static bool warned = false;
    if (!warned) {
        fprintf(stderr,
                "HilbertMatrix::wprec: arbitrary precision unavailable "
                "(USE_MPFR=0); using machine precision.\n");
        warned = true;
    }
    *mode = HM_PREC_MACHINE;
    *bits = 0;
#endif
    return true;
}

/* Build a single entry 1/k (k = i + j - 1, always >= 1). */
static Expr* hm_entry(int64_t k, hm_prec_mode mode, long bits) {
    switch (mode) {
        case HM_PREC_MACHINE:
            return expr_new_real(1.0 / (double)k);
#ifdef USE_MPFR
        case HM_PREC_MPFR: {
            mpfr_t out;
            mpfr_init2(out, bits);
            mpfr_set_ui(out, 1, MPFR_RNDN);
            mpfr_div_ui(out, out, (unsigned long)k, MPFR_RNDN);
            return expr_new_mpfr_move(out);
        }
#else
        case HM_PREC_MPFR: /* unreachable: parser never selects MPFR here */
            return expr_new_real(1.0 / (double)k);
#endif
        case HM_PREC_EXACT:
        default:
            /* make_rational reduces 1/1 to the integer 1. */
            return make_rational(1, k);
    }
}

Expr* builtin_hilbertmatrix(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 0) {
        fprintf(stderr,
                "HilbertMatrix::argx: HilbertMatrix called with 0 arguments; "
                "1 argument is expected.\n");
        return NULL;
    }

    /* args[0] is the dimension spec; args[1..] are options. */
    int64_t m = 0, n = 0;
    if (!hm_parse_dims(res->data.function.args[0], &m, &n)) {
        char* s = expr_to_string(res->data.function.args[0]);
        fprintf(stderr,
                "HilbertMatrix::dims: Dimension specification %s should be a "
                "positive machine integer or a pair of positive machine "
                "integers.\n",
                s ? s : "?");
        free(s);
        return NULL;
    }

    /* Parse trailing options.  Only WorkingPrecision is recognised; any
     * non-Rule (or unknown-option Rule) trailing argument is reported via
     * ::nonopt, with the rightmost offender named (Wolfram semantics). */
    hm_prec_mode mode = HM_PREC_EXACT;
    long bits = 0;
    Expr* last_bad = NULL;
    for (size_t i = 1; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        bool is_rule = (hm_is_call(opt, SYM_Rule) ||
                        hm_is_call(opt, SYM_RuleDelayed)) &&
                       opt->data.function.arg_count == 2 &&
                       opt->data.function.args[0]->type == EXPR_SYMBOL;
        if (!is_rule) {
            last_bad = opt;
            continue;
        }
        const char* name = opt->data.function.args[0]->data.symbol;
        if (name == SYM_WorkingPrecision) {
            /* Last valid setting wins; an unparseable value is ignored. */
            hm_parse_working_precision(opt->data.function.args[1],
                                       &mode, &bits);
        } else {
            last_bad = opt;
        }
    }
    if (last_bad != NULL) {
        char* bad_str = expr_to_string(last_bad);
        char* call_str = expr_to_string(res);
        fprintf(stderr,
                "HilbertMatrix::nonopt: Options expected (instead of %s) "
                "beyond position 1 in %s. An option must be a rule or a list "
                "of rules.\n",
                bad_str ? bad_str : "?", call_str ? call_str : "?");
        free(bad_str);
        free(call_str);
        return NULL;
    }

    /* Build the m x n matrix as a List of Lists; entry (i, j) = 1/(i+j-1). */
    Expr** rows = malloc(sizeof(Expr*) * (size_t)m);
    for (int64_t i = 0; i < m; i++) {
        Expr** cells = malloc(sizeof(Expr*) * (size_t)n);
        for (int64_t j = 0; j < n; j++) {
            cells[j] = hm_entry(i + j + 1, mode, bits);  /* (i+1)+(j+1)-1 */
        }
        rows[i] = expr_new_function(expr_new_symbol("List"), cells, (size_t)n);
        free(cells);
    }
    Expr* result = expr_new_function(expr_new_symbol("List"), rows, (size_t)m);
    free(rows);
    return result;
}
