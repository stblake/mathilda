/*
 * product_infinite.c -- Product`Infinite (Stage 4): convergent infinite
 * products via the limit of the finite closed form.
 *
 * For imax == Infinity:
 *   1. Convergence gate.  For a rational f = P/Q the product converges iff
 *      deg P == deg Q AND the leading AND next-to-leading coefficients agree
 *      (so the term ~ 1 + O(1/k^2), not 1 + c/k).  This separates the
 *      convergent a_k ~ 1 + c/k^2 from the divergent a_k ~ 1 + c/k.  A failed
 *      gate prints Product::div and leaves the product unevaluated (mirroring
 *      Sum's divergent handling).
 *   2. Run the finite cascade for a symbolic upper bound U (Product[f,{k,lo,U}])
 *      then take Limit[CF, U -> Infinity].
 *
 * Transcendental infinite products whose finite form is not elementary
 * (prod (1 +- a/k^2) = Sinh/Sin, Wallis, ...) are recognised by the rule table
 * in src/internal/product_infinite.m, which fires after this builtin (and the
 * rest of the cascade) falls through.
 */

#include "product_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "arithmetic.h"   /* arith_warnings_muted */
#include "sym_names.h"
#include "poly.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

int g_product_verify_convergence = 1;

static bool is_inf_sym(const Expr* e) {
    return e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Infinity;
}

static bool head_sym_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol.name, name) == 0;
}

static void warn_div(void) {
    if (!arith_warnings_muted())
        fprintf(stderr, "Product::div: Product does not converge.\n");
}

/* Rational convergence classifier: 1 converges, 0 diverges, -1 unknown
 * (non-rational; the .m table / limit route may still resolve it). */
static int rational_converges(Expr* f, Expr* var) {
    if (prod_has_symbolic_power(f, var)) return -1;
    Expr* tg  = prod_eval("Together",    (Expr*[]){ expr_copy(f) }, 1);
    Expr* num = prod_eval("Numerator",   (Expr*[]){ expr_copy(tg) }, 1);
    Expr* den = prod_eval("Denominator", (Expr*[]){ expr_copy(tg) }, 1);
    expr_free(tg);
    Expr* enum_ = prod_eval("Expand", (Expr*[]){ num }, 1);  /* adopts num */
    Expr* eden  = prod_eval("Expand", (Expr*[]){ den }, 1);  /* adopts den */

    Expr* vars[1] = { var };
    int conv = -1;
    if (is_polynomial(enum_, vars, 1) && is_polynomial(eden, vars, 1)) {
        int dp = get_degree_poly(enum_, var);
        int dq = get_degree_poly(eden, var);
        if (dp != dq) {
            conv = 0;
        } else {
            Expr* lp = get_coeff(enum_, var, dp);
            Expr* lq = get_coeff(eden, var, dq);
            bool lead_eq = expr_eq(lp, lq);
            bool next_eq = true;
            if (dp >= 1) {
                Expr* np = get_coeff(enum_, var, dp - 1);
                Expr* nq = get_coeff(eden, var, dq - 1);
                next_eq = expr_eq(np, nq);
                expr_free(np); expr_free(nq);
            }
            expr_free(lp); expr_free(lq);
            conv = (lead_eq && next_eq) ? 1 : 0;
        }
    }
    expr_free(enum_); expr_free(eden);
    return conv;
}

/* Weierstrass family:  prod_{k=1}^Inf (1 + c/k^2) = Sinh[Pi Sqrt[c]]/(Pi Sqrt[c])
 * for any c free of k (Sqrt of a negative c yields the Sin form automatically,
 * e.g. c = -z^2/Pi^2 gives Sin[z]/z).  Recognised only at imin == 1 (the
 * canonical Weierstrass lower bound); other lower bounds go to the limit route.
 * Returns the closed form, or NULL if f is not of this shape. */
static Expr* try_sinh_family(Expr* f, Expr* var, Expr* imin) {
    if (!(imin->type == EXPR_INTEGER && imin->data.integer == 1)) return NULL;

    /* A genuine Weierstrass body 1 + c/k^2 (c free of k) never carries a
     * symbolic power of the index; bodies like Exp[-x/k] = Power[E,-x/k] do.
     * Skipping them here loses no coverage and avoids running Simplify on an
     * expression whose intermediate forms evaluate Power[0,-1] (spurious
     * Power::infy: 1/0 warnings). */
    if (prod_has_symbolic_power(f, var)) return NULL;

    /* c = Simplify[(f - 1) * k^2]; the form holds iff c is free of k. */
    Expr* fm1 = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(f), expr_new_integer(-1) }, 2);
    Expr* k2 = expr_new_function(expr_new_symbol(SYM_Power),
                   (Expr*[]){ expr_copy(var), expr_new_integer(2) }, 2);
    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ fm1, k2 }, 2);
    /* Mute arithmetic warnings across the speculative Simplify (defense in
     * depth: the gate above already excludes symbolic-power bodies). */
    arith_warnings_mute_push();
    Expr* c = prod_eval("Simplify", (Expr*[]){ prod }, 1);   /* adopts prod */
    arith_warnings_mute_pop();

    if (!prod_free_of(c, var)) { expr_free(c); return NULL; }
    if (c->type == EXPR_INTEGER && c->data.integer == 0) { expr_free(c); return NULL; }

    /* Sinh[Pi Sqrt[c]] / (Pi Sqrt[c]). */
    Expr* sq = expr_new_function(expr_new_symbol("Sqrt"), (Expr*[]){ c }, 1);  /* adopts c */
    Expr* pisq = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ expr_new_symbol("Pi"), expr_copy(sq) }, 2);
    Expr* sinh = expr_new_function(expr_new_symbol("Sinh"),
                     (Expr*[]){ expr_copy(pisq) }, 1);
    Expr* out = prod_div(sinh, pisq);
    expr_free(sinh); expr_free(pisq); expr_free(sq);
    return out;
}

Expr* builtin_product_infinite(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!product_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite || !is_inf_sym(imax)) return NULL;

    if (g_product_verify_convergence) {
        int conv = rational_converges(f, var);
        if (conv == 0) { warn_div(); return NULL; }
    }

    /* Transcendental Weierstrass sine/cosine family (closed form in C). */
    Expr* sinh = try_sinh_family(f, var, imin);
    if (sinh) return sinh;

    /* Finite closed form for a fresh symbolic upper bound, then its limit.
     * The temporary bound uses an internal context symbol that the user cannot
     * bind, and imax==U (not Infinity) so this stage does not recurse. */
    Expr* U = expr_new_symbol("Product`Infinite`U");
    Expr* spec = expr_new_function(expr_new_symbol(SYM_List),
                     (Expr*[]){ expr_copy(var), expr_copy(imin), expr_copy(U) }, 3);
    Expr* pj = expr_new_function(expr_new_symbol(SYM_Product),
                   (Expr*[]){ expr_copy(f), spec }, 2);
    Expr* CF = evaluate(pj);
    expr_free(pj);

    /* No elementary finite form -> let the .m identity table try. */
    if (head_sym_is(CF, "Product")) { expr_free(CF); expr_free(U); return NULL; }

    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
                     (Expr*[]){ expr_copy(U), expr_new_symbol(SYM_Infinity) }, 2);
    Expr* lim = expr_new_function(expr_new_symbol("Limit"),
                    (Expr*[]){ CF, rule }, 2);   /* adopts CF */
    Expr* out = evaluate(lim);
    expr_free(lim);
    expr_free(U);

    /* Limit could not resolve -> fall through. */
    if (head_sym_is(out, "Limit")) { expr_free(out); return NULL; }
    /* Divergent limit -> Product::div, stay unevaluated. */
    if ((out->type == EXPR_SYMBOL &&
         (out->data.symbol.name == SYM_Infinity || out->data.symbol.name == SYM_ComplexInfinity))
        || head_sym_is(out, "DirectedInfinity")) {
        warn_div();
        expr_free(out);
        return NULL;
    }
    return out;
}

void product_infinite_init(void) {
    symtab_add_builtin("Product`Infinite", builtin_product_infinite);
    symtab_get_def("Product`Infinite")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Product`Infinite",
        "Product`Infinite[f, i, imin, Infinity] evaluates a convergent infinite "
        "product as the limit of its finite closed form, after a rational "
        "convergence gate (equal degree, leading and next-to-leading "
        "coefficients). A divergent product gives Product::div and stays "
        "unevaluated. Transcendental forms are handled by product_infinite.m.");
}
