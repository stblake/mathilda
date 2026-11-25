/* factorialpower.c -- FactorialPower[].
 * Split from numbertheory.c; see numbertheory.h and
 * numbertheory_internal.h for the subsystem layout. */

#include "numbertheory.h"
#include "numbertheory_internal.h"
#include "numeric.h"   /* expr_is_numeric_quantity */
#include "arithmetic.h"
#include "eval.h"
#include "sym_names.h"
#include "internal.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <gmp.h>

/* FactorialPower[n, k] -- the falling factorial (Pochhammer-style).
 *
 *   FactorialPower[n, k] = n * (n - 1) * (n - 2) * ... * (n - k + 1)
 *                        = n! / (n - k)!     (for non-negative integer k)
 *                        = 1                  (k == 0)
 *
 * Concrete numeric paths:
 *   - Both arguments integer with k >= 0: compute the product (GMP for
 *     k > 30 to avoid int64 overflow).
 *   - Numeric n with negative integer k: 1/((n + 1) (n + 2) ... (n + |k|)).
 *   - FactorialPower[x, n, h], the step-h form: see fp_three_arg.
 *
 * Symbolic n with concrete non-negative integer k expands to a product
 * of (n - i) factors so downstream Expand / D can act on it -- e.g.
 *   FactorialPower[n, 3] -> n (n - 1) (n - 2)
 * Mathematically equal under the reference's definition and lets the
 * surrounding evaluator collapse it via Expand[] when the user wants.
 *
 * Anything else (k symbolic, or both arguments symbolic): NULL --
 * stays as `FactorialPower[n, k]` until more information is supplied.
 */
/* True for an inexact numeric leaf: a machine/MPFR real, or a Complex[..] with
 * an inexact part.  Gates the analytic (Gamma) path so a non-integer k gets a
 * value under N while exact args stay symbolic. */
static bool fp_is_inexact(const Expr* e) {
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) return fp_is_inexact(re) || fp_is_inexact(im);
    return false;
}

/* Generalised falling factorial for non-integer k, built as an expression and
 * evaluated (reusing Gamma and its complex + MPFR kernels):
 *
 *     FactorialPower[n, k] = Gamma[n + 1] / Gamma[n - k + 1].
 */
static Expr* fp_gamma_quotient(Expr* n, Expr* k) {
    Expr* np1 = expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy(n), expr_new_integer(1) }, 2);
    Expr* gn  = expr_new_function(expr_new_symbol(SYM_Gamma), &np1, 1);

    Expr* negk = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(-1), expr_copy(k) }, 2);
    Expr* nmk1 = expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy(n), negk, expr_new_integer(1) }, 3);   /* n - k + 1 */
    Expr* gnmk = expr_new_function(expr_new_symbol(SYM_Gamma), &nmk1, 1);
    Expr* inv  = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ gnmk, expr_new_integer(-1) }, 2);

    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ gn, inv }, 2);
    return eval_and_free(prod);
}

static Expr* fp_call(const char* head, Expr** a, size_t n) {
    return expr_new_function(expr_new_symbol(head), a, n);
}

/* x + c h, unevaluated, for the product factors. */
static Expr* fp_shift(Expr* x, int64_t c, Expr* h) {
    if (c == 0) return expr_copy(x);
    Expr* ch = fp_call(SYM_Times, (Expr*[]){ expr_new_integer(c), expr_copy(h) }, 2);
    return fp_call(SYM_Plus, (Expr*[]){ expr_copy(x), ch }, 2);
}

/* A number, or a numeric quantity (Pi, Sqrt[2]). */
static bool fp_numeric(const Expr* e) {
    return expr_is_numeric_quantity(e);
}

/* FactorialPower[x, n, h] -- the falling factorial with step h:
 *
 *   x (x - h) (x - 2 h) ... (x - (n - 1) h)            n a positive integer
 *   1 / ((x + h) (x + 2 h) ... (x + |n| h))             n a negative integer
 *   h^n Gamma[x/h + 1] / Gamma[x/h - n + 1]             any other numeric n
 *
 * with h = 1 the two-argument form. Evaluated when x and h are numeric (and,
 * for a non-integer n, every argument is a number or one is inexact); a
 * symbolic x or h stays unevaluated, as in the reference. */
static Expr* fp_step_product(Expr* x, Expr* n, Expr* h) {
    if (n->type == EXPR_INTEGER && n->data.integer > 0 && n->data.integer <= 100000) {
        int64_t k = n->data.integer;
        Expr** f = malloc(sizeof(Expr*) * (size_t)k);
        for (int64_t i = 0; i < k; i++) f[i] = fp_shift(x, -i, h);
        Expr* prod = fp_call(SYM_Times, f, (size_t)k);
        free(f);
        return eval_and_free(prod);
    }
    if (n->type == EXPR_INTEGER && n->data.integer < 0 && n->data.integer >= -100000) {
        int64_t k = -n->data.integer;
        Expr** f = malloc(sizeof(Expr*) * (size_t)k);
        for (int64_t i = 0; i < k; i++) f[i] = fp_shift(x, i + 1, h);
        Expr* prod = fp_call(SYM_Times, f, (size_t)k);
        free(f);
        return eval_and_free(fp_call(SYM_Power, (Expr*[]){ prod, expr_new_integer(-1) }, 2));
    }
    return NULL;
}

static Expr* fp_three_arg(Expr* x, Expr* n, Expr* h) {
    if (n->type == EXPR_INTEGER && n->data.integer == 0) return expr_new_integer(1);
    if (h->type == EXPR_INTEGER && h->data.integer == 1)
        return fp_call("FactorialPower", (Expr*[]){ expr_copy(x), expr_copy(n) }, 2);
    if (!fp_numeric(x) || !fp_numeric(h)) return NULL;
    Expr* p = fp_step_product(x, n, h);
    if (p || !fp_numeric(n)) return p;
    /* h^n Gamma[x/h + 1] / Gamma[x/h - n + 1]. */
    Expr* xh = eval_and_free(fp_call(SYM_Times, (Expr*[]){ expr_copy(x),
        fp_call(SYM_Power, (Expr*[]){ expr_copy(h), expr_new_integer(-1) }, 2) }, 2));
    Expr* q = fp_gamma_quotient(xh, n);
    expr_free(xh);
    Expr* hn = fp_call(SYM_Power, (Expr*[]){ expr_copy(h), expr_copy(n) }, 2);
    Expr* v = eval_and_free(fp_call(SYM_Times, (Expr*[]){ hn, q }, 2));
    if (v && (expr_is_numeric_like(v) || fp_numeric(v))) return v;
    expr_free(v);
    return NULL;
}

Expr* builtin_factorialpower(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count == 3)
        return fp_three_arg(res->data.function.args[0], res->data.function.args[1],
                            res->data.function.args[2]);
    if (res->data.function.arg_count != 2) return NULL;
    Expr* n_arg = res->data.function.args[0];
    Expr* k_arg = res->data.function.args[1];

    /* k must be a concrete non-negative integer for the exact product path.
     * A non-integer k with an inexact operand (under N) takes the generalised
     * Gamma quotient; otherwise the expression stays symbolic. */
    if (k_arg->type != EXPR_INTEGER) {
        if (fp_is_inexact(n_arg) || fp_is_inexact(k_arg)) {
            Expr* val = fp_gamma_quotient(n_arg, k_arg);
            if (val && expr_is_numeric_like(val)) return val;
            expr_free(val);
        }
        return NULL;
    }
    int64_t k = k_arg->data.integer;
    if (k < 0) {
        /* 1/((n + 1) (n + 2) ... (n + |k|)) for a numeric n. */
        if (!fp_numeric(n_arg) || k < -100000) return NULL;
        Expr* one = expr_new_integer(1);
        Expr* r = fp_step_product(n_arg, k_arg, one);
        expr_free(one);
        return r;
    }
    if (k == 0) return expr_new_integer(1);

    /* Both n and k integer: compute exact product in GMP. */
    if (n_arg->type == EXPR_INTEGER || n_arg->type == EXPR_BIGINT) {
        mpz_t n_z, term, product;
        mpz_init(n_z);
        if (n_arg->type == EXPR_INTEGER) mth_mpz_set_i64(n_z, n_arg->data.integer);
        else                              mpz_set(n_z, n_arg->data.bigint);
        mpz_init_set_ui(product, 1);
        mpz_init(term);
        for (int64_t i = 0; i < k; i++) {
            mpz_set(term, n_z);
            mpz_sub_ui(term, term, (unsigned long)i);
            mpz_mul(product, product, term);
        }
        mpz_clear(n_z);
        mpz_clear(term);
        Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(product));
        mpz_clear(product);
        return result;
    }

    /* Symbolic n, k a small concrete integer (<= 32 to keep tree size
     * bounded). Expand to Times[n, n-1, ..., n-k+1]. */
    if (k <= 32) {
        Expr** factors = malloc(sizeof(Expr*) * (size_t)k);
        for (int64_t i = 0; i < k; i++) {
            if (i == 0) {
                factors[i] = expr_copy(n_arg);
            } else {
                Expr* shift = expr_new_integer(-i);
                factors[i] = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(n_arg), shift }, 2);
            }
        }
        Expr* product = expr_new_function(expr_new_symbol(SYM_Times), factors, (size_t)k);
        free(factors);
        return eval_and_free(product);
    }
    /* A numeric n (1/2, 2.5, Pi) with a larger k: the same product, as the
     * negative-k branch forms it (FactorialPower[1/2, 33] is an exact
     * rational, FactorialPower[2.5, 33] a machine real). */
    if (fp_numeric(n_arg)) {
        Expr* one = expr_new_integer(1);
        Expr* r = fp_step_product(n_arg, k_arg, one);
        expr_free(one);
        return r;
    }
    return NULL;
}
