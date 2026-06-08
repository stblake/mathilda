/* Mathilda -- the Pochhammer symbol (rising factorial).
 *
 *   Pochhammer[a, n] = a (a+1) ... (a+n-1) = Gamma[a+n] / Gamma[a]
 *
 * The implementation deliberately holds almost no numeric code of its own:
 * two existing, fully-tested mechanisms do the heavy lifting.
 *
 *   1. Times -- for an exact integer n the symbol expands to the product of
 *      n linear factors, Times[a, a+1, ..., a+n-1]. Evaluating that product
 *      collapses to an exact (BigInt / Rational) value for numeric a,
 *      preserves MPFR precision when a is an EXPR_MPFR, does complex
 *      arithmetic when a is a Complex[..], and stays a symbolic polynomial
 *      product for symbolic a. One path covers every kind of a.
 *
 *   2. Gamma -- for a non-integer (or out-of-range) n the symbol evaluates
 *      Gamma[a+n]/Gamma[a], reusing the Gamma builtin's exact half-integer
 *      reductions (-> rational multiples of Sqrt[Pi]), its libm / MPFR real
 *      paths, and its machine-precision complex Lanczos path -- for free.
 *
 * Mirrors how gamma.c reuses Factorial and how FactorialPower (the falling
 * factorial, numbertheory.c) builds a Times product.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "pochhammer.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "attr.h"
#include "eval.h"          /* eval_and_free */
#include "expr.h"          /* expr_is_numeric_like */
#include "symtab.h"

/* Largest |n| for which Pochhammer[a, n] is expanded into its explicit
 * product of n linear factors. Beyond this the tree is more noise than help;
 * numeric a still resolves through the Gamma ratio, symbolic a stays
 * symbolic. Mirrors GAMMA_INT_EXPAND_CAP in gamma.c. */
#define POCH_PRODUCT_CAP 1000

/* ------------------------------------------------------------------ */
/* Small helpers                                                      */
/* ------------------------------------------------------------------ */

/* True if `e` is exactly the symbol `name`. */
static bool poch_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol, name) == 0;
}

/* True if `e` (or a subexpression) is headed by Gamma -- used to detect a
 * Gamma ratio that did not fully reduce, so we can leave Pochhammer
 * symbolic rather than emit a half-evaluated answer. */
static bool poch_contains_gamma(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (poch_is_symbol(e->data.function.head, "Gamma")) return true;
    if (poch_contains_gamma(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (poch_contains_gamma(e->data.function.args[i])) return true;
    return false;
}

/* ------------------------------------------------------------------ */
/* Product path: Pochhammer[a, n] for exact integer n                 */
/* ------------------------------------------------------------------ */

/* For nonzero integer n with |n| <= POCH_PRODUCT_CAP, build and evaluate
 *   n > 0:  Times[a, a+1, ..., a+n-1]
 *   n < 0:  1 / Times[a-1, a-2, ..., a-|n|]
 * (matching Mathematica's Pochhammer[x,-5] = 1/((x-5)(x-4)(x-3)(x-2)(x-1))).
 * The linear factors are Plus[a, k]; expr_new_function copies the factor
 * array, so we free the container but not the adopted elements. */
static Expr* poch_build_product(Expr* a, int64_t n) {
    bool negative = n < 0;
    int64_t m = negative ? -n : n;           /* number of factors */
    size_t cnt = (size_t)m;

    Expr** factors = (Expr**)malloc(cnt * sizeof(Expr*));
    if (!factors) return NULL;

    for (int64_t i = 0; i < m; i++) {
        /* offsets: 0,1,...,m-1 for n>0;  -1,-2,...,-m for n<0. */
        int64_t off = negative ? -(i + 1) : i;
        if (off == 0) {
            factors[i] = expr_copy(a);
        } else {
            factors[i] = expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ expr_copy(a), expr_new_integer(off) }, 2);
        }
    }

    Expr* product = expr_new_function(expr_new_symbol("Times"), factors, cnt);
    free(factors);

    if (negative) {
        product = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ product, expr_new_integer(-1) }, 2);
    }
    return eval_and_free(product);
}

/* ------------------------------------------------------------------ */
/* Numeric path: Pochhammer[a, n] = Gamma[a+n] / Gamma[a]             */
/* ------------------------------------------------------------------ */

/* Evaluate Gamma[a+n]/Gamma[a]. Returns the reduced value, or NULL if the
 * ratio failed to reduce (a residual Gamma head survived), leaving the
 * Pochhammer call symbolic. Only meaningful for numeric a and n. */
static Expr* poch_via_gamma(Expr* a, Expr* n) {
    Expr* apn = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ expr_copy(a), expr_copy(n) }, 2);
    Expr* g_top = expr_new_function(expr_new_symbol("Gamma"), &apn, 1);

    Expr* ga = expr_copy(a);
    Expr* g_bot = expr_new_function(expr_new_symbol("Gamma"), &ga, 1);
    Expr* g_bot_inv = expr_new_function(expr_new_symbol("Power"),
        (Expr*[]){ g_bot, expr_new_integer(-1) }, 2);

    Expr* ratio = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ g_top, g_bot_inv }, 2);
    Expr* out = eval_and_free(ratio);

    if (out && poch_contains_gamma(out)) { expr_free(out); return NULL; }
    return out;
}

/* ------------------------------------------------------------------ */
/* Dispatch                                                           */
/* ------------------------------------------------------------------ */

static Expr* pochhammer_two_arg(Expr* a, Expr* n) {
    bool n_is_int = (n->type == EXPR_INTEGER);

    /* 1. Pochhammer[a, 0] = 1, for any a (including symbolic / Infinity). */
    if (n_is_int && n->data.integer == 0) return expr_new_integer(1);

    /* 2. Pochhammer[0, n] = 0 for a positive integer n (a zero factor zeros
     *    the product); short-circuits before any large factorial. */
    if (a->type == EXPR_INTEGER && a->data.integer == 0 &&
        n_is_int && n->data.integer > 0) {
        return expr_new_integer(0);
    }

    /* 3. Pochhammer[Infinity, n] = Infinity for a positive integer n. */
    if (poch_is_symbol(a, "Infinity") && n_is_int && n->data.integer > 0) {
        return expr_new_symbol("Infinity");
    }

    /* 4. Exact integer n within the cap: expand to the linear-factor product.
     *    Works for symbolic a (polynomial form) and numeric a alike. */
    if (n_is_int) {
        int64_t nv = n->data.integer;
        int64_t absn = nv < 0 ? -nv : nv;
        if (absn <= POCH_PRODUCT_CAP) return poch_build_product(a, nv);
    }

    /* 5. Numeric a and n with non-integer / out-of-range n: Gamma ratio.
     *    Covers exact half-integers (-> Sqrt[Pi]), machine and arbitrary
     *    precision reals, and machine-precision complex values. */
    if (expr_is_numeric_like(a) && expr_is_numeric_like(n)) {
        Expr* out = poch_via_gamma(a, n);
        if (out) return out;
    }

    /* 6. Otherwise leave symbolic. */
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

/* Mathematica-compatible argrx diagnostic for a wrong argument count;
 * returns NULL so the evaluator leaves the call unevaluated. */
static Expr* poch_emit_argrx(size_t argc) {
    fprintf(stderr,
            "Pochhammer::argrx: Pochhammer called with %zu argument%s; "
            "2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_pochhammer(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return poch_emit_argrx(argc);
    return pochhammer_two_arg(res->data.function.args[0],
                              res->data.function.args[1]);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void pochhammer_init(void) {
    symtab_add_builtin("Pochhammer", builtin_pochhammer);
    symtab_get_def("Pochhammer")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
