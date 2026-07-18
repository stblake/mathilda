/* Mathilda — common preprocessing helpers.  See common.h. */

#include "common.h"

#include "numeric.h"
#include "rationalize.h"

#include <stdio.h>
#include <stdint.h>

#ifdef USE_MPFR
#include <mpfr.h>
#endif

bool head_is(const Expr* e, const char* sym) {
    return e && sym &&
           e->type == EXPR_FUNCTION &&
           e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym;
}

Expr* builtin_arg_error(const char* head, size_t got, size_t min, size_t max) {
    const char* got_s = (got == 1) ? "" : "s"; /* "argument" vs "arguments" */

    if (min == max) {
        /* Fixed arity. Mathematica tags the one-argument case `argx`. */
        const char* tag = (min == 1) ? "argx" : "argrx";
        fprintf(stderr,
                "%s::%s: %s called with %zu argument%s; "
                "%zu argument%s %s expected.\n",
                head, tag, head, got, got_s,
                min, (min == 1) ? "" : "s", (min == 1) ? "is" : "are");
    } else if (max == SIZE_MAX) {
        fprintf(stderr,
                "%s::argm: %s called with %zu argument%s; "
                "%zu or more arguments are expected.\n",
                head, head, got, got_s, min);
    } else if (max == min + 1) {
        fprintf(stderr,
                "%s::argt: %s called with %zu argument%s; "
                "%zu or %zu arguments are expected.\n",
                head, head, got, got_s, min, max);
    } else {
        fprintf(stderr,
                "%s::argb: %s called with %zu argument%s; "
                "between %zu and %zu arguments are expected.\n",
                head, head, got, got_s, min, max);
    }
    return NULL;
}

/* Bit precision of a single inexact leaf.  Real → 53 (IEEE 754 double).
 * MPFR → mpfr_get_prec().  Caller is expected to only invoke this on
 * leaves where is_inexact_leaf would have returned true. */
static long leaf_precision_bits(const Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_REAL) return 53;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return (long)mpfr_get_prec(e->data.mpfr);
#endif
    return 0;
}

static bool is_inexact_leaf(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* Recursive walk: update *info in place with min precision and has_inexact
 * flag.  We track the running minimum so a single pass captures both. */
static void scan_inexact_walk(const Expr* e, CommonInexactInfo* info) {
    if (!e) return;
    if (is_inexact_leaf(e)) {
        long bits = leaf_precision_bits(e);
        if (!info->has_inexact || bits < info->min_bits) {
            info->min_bits = bits;
        }
        info->has_inexact = true;
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    scan_inexact_walk(e->data.function.head, info);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        scan_inexact_walk(e->data.function.args[i], info);
    }
}

CommonInexactInfo common_scan_inexact(const Expr* e) {
    CommonInexactInfo info = { false, 0 };
    scan_inexact_walk(e, &info);
    return info;
}

Expr* common_rationalize_input(const Expr* e, long bits) {
    return internal_force_rationalize_bits(e, bits);
}

Expr* common_numericalize_result(const Expr* e, long bits) {
    if (!e) return NULL;
    NumericSpec spec = numeric_machine_spec();
#ifdef USE_MPFR
    if (bits > 53) {
        spec.mode = NUMERIC_MODE_MPFR;
        spec.bits = bits;
    }
#else
    (void)bits;
#endif
    return numericalize(e, spec);
}
