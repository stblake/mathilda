#include "simp.h"
#include "simp_internal.h"
#include "arithmetic.h"
#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expand.h"
#include "facpoly.h"
#include "numeric.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "expr.h"
#include "rationalize.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "trigrat.h"
#include "qa.h"
#include "qafactor.h"
#include "simp_log.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif


/* ----------------------------------------------------------------------- */
/* Default complexity measure                                              */
/* ----------------------------------------------------------------------- */

static size_t int_digit_count_int64(int64_t v) {
    if (v == 0) return 1;
    if (v < 0) {
        /* INT64_MIN edge case: |v| not representable; count digits of the
         * negated value digit-at-a-time without ever forming -INT64_MIN. */
        size_t n = 1;
        int64_t t = v;
        while (t <= -10) { n++; t /= 10; }
        return n;
    }
    size_t n = 0;
    while (v > 0) { n++; v /= 10; }
    return n;
}

/*
 * simp_default_complexity implements Mathematica's SimplifyCount:
 *
 *   Symbol      -> 1
 *   Integer 0   -> 1
 *   Integer p>0 -> Floor[Log10[p]] + 1            == digits(p)
 *   Integer p<0 -> Floor[Log10[|p|]] + 2          == digits(|p|) + 1
 *   Rational    -> SimplifyCount[num] + SimplifyCount[den] + 1
 *   Complex     -> SimplifyCount[re]  + SimplifyCount[im]  + 1
 *   Real / MPFR -> 2                              (NumberQ but not Integer/Rational)
 *   String      -> 1                              (treated as a leaf, Mathilda extension)
 *   Function    -> SimplifyCount[head] + sum SimplifyCount[args]
 *
 * The negative-integer adjustment matches Mathematica's behaviour where
 * the leading "-" contributes one unit of complexity. The explicit
 * Rational/Complex cases keep e.g. 100 Log[2] (score 6) preferred over
 * Log[2^100] (score 32). */
size_t simp_default_complexity(const Expr* e) {
    if (!e) return 0;
    switch (e->type) {
        case EXPR_INTEGER: {
            int64_t v = e->data.integer;
            if (v == 0) return 1;
            size_t d = int_digit_count_int64(v);
            return v > 0 ? d : d + 1;
        }
        case EXPR_BIGINT: {
            int sgn = mpz_sgn(e->data.bigint);
            if (sgn == 0) return 1;
            size_t digits = mpz_sizeinbase(e->data.bigint, 10);
            return sgn > 0 ? digits : digits + 1;
        }
        case EXPR_REAL:    return 2;
        case EXPR_SYMBOL:  return 1;
        case EXPR_STRING:  return 1;
        case EXPR_FUNCTION: {
            const Expr* head = e->data.function.head;
            size_t argc = e->data.function.arg_count;
            /* Rational[n, d] and Complex[re, im] are Mathematica-specials:
             * SimplifyCount adds 1 for the wrapper, not the head's own
             * SimplifyCount. */
            if (head && head->type == EXPR_SYMBOL && argc == 2) {
                if (head->data.symbol == SYM_Rational ||
                    head->data.symbol == SYM_Complex) {
                    return simp_default_complexity(e->data.function.args[0])
                         + simp_default_complexity(e->data.function.args[1])
                         + 1;
                }
            }
            size_t total = simp_default_complexity(head);
            for (size_t i = 0; i < argc; i++) {
                total += simp_default_complexity(e->data.function.args[i]);
            }
            return total;
        }
#ifdef USE_MPFR
        case EXPR_MPFR: return 2;
#endif
    }
    return 1;
}

/* Builtin SimplifyCount[expr] -- exposes the default complexity to users
 * so they can inspect or use it inside a custom ComplexityFunction.
 * The caller (evaluate_step) frees `res` after we return a non-NULL Expr;
 * we must NOT free it here. */
Expr* builtin_simplify_count(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1) return NULL;
    size_t s = simp_default_complexity(res->data.function.args[0]);
    /* size_t comfortably fits in EXPR_INTEGER for any expression we'd
     * realistically see; on 64-bit size_t = 8 bytes, int64_t = 8 bytes. */
    return expr_new_integer((int64_t)s);
}

