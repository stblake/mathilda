/*
 * product_rational.c -- Product`Rational (Stage 2), the workhorse.
 *
 * Closed form for any rational f(k) whose numerator and denominator factor
 * into linear factors over Q.  Output is built from Pochhammer (and the powers
 * of a leading constant).
 *
 * Atomic identity (definite, unit step):
 *   prod_{k=imin}^{imax} (k + a) = Pochhammer[imin + a, imax - imin + 1].
 * A linear factor (k - r) has a = -r, so it contributes
 *   Pochhammer[imin - r, count]^mult     (numerator: positive exponent)
 *   Pochhammer[imin - r, count]^-mult     (denominator: divides)
 * with count = imax - imin + 1.  The leading constant c = lead_num / lead_den
 * contributes c^count.
 *
 * Indefinite Product[f, i] uses the Wolfram convention P(i+1)/P(i) = f(i),
 * i.e. P(i) = prod_{k=1}^{i-1} f(k): replace imin -> 1 and count -> (i - 1), so
 * Product[k, k] = Pochhammer[1, i-1] = (i-1)! = Gamma[i].
 *
 * f is rational in k => no k in any exponent => Factor/root extraction are
 * hang-safe (prod_has_symbolic_power gates the symbolic-exponent case out).
 */

#include "product_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>

/* Pochhammer[copy a, copy b]. */
static Expr* poch(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(SYM_Pochhammer),
               (Expr*[]){ expr_copy(a), expr_copy(b) }, 2);
}

/* Power[copy base, e] with small integer exponent e (e may be negative). */
static Expr* powi(Expr* base, int e) {
    return expr_new_function(expr_new_symbol(SYM_Power),
               (Expr*[]){ expr_copy(base), expr_new_integer(e) }, 2);
}

/* anchor - r, evaluated. */
static Expr* anchor_minus_root(Expr* anchor, Expr* r) {
    Expr* negr = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ expr_new_integer(-1), expr_copy(r) }, 2);
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(anchor), negr }, 2);
    Expr* r2 = evaluate(sum);
    expr_free(sum);
    return r2;
}

static void free_roots(Expr** roots, int* mults, size_t n) {
    for (size_t i = 0; i < n; i++) expr_free(roots[i]);
    free(roots);
    free(mults);
}

/* Rewrite every Pochhammer[1, x] node to Factorial[x] so the canonical
 * n! form appears (Pochhammer[1, n] == Gamma[1+n] == n!).  Returns a new
 * owned tree; the input is not consumed. */
static Expr* poch1_to_fact(const Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    const Expr* head = e->data.function.head;
    size_t argc = e->data.function.arg_count;
    if (head->type == EXPR_SYMBOL && head->data.symbol == SYM_Pochhammer
            && argc == 2
            && e->data.function.args[0]->type == EXPR_INTEGER
            && e->data.function.args[0]->data.integer == 1) {
        Expr* x = poch1_to_fact(e->data.function.args[1]);
        return expr_new_function(expr_new_symbol(SYM_Factorial), &x, 1);
    }
    Expr* h = poch1_to_fact(head);
    Expr** na = malloc(sizeof(Expr*) * (argc ? argc : 1));
    for (size_t i = 0; i < argc; i++) na[i] = poch1_to_fact(e->data.function.args[i]);
    Expr* r = expr_new_function(h, na, argc);
    free(na);
    return r;
}

Expr* builtin_product_rational(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!product_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;

    /* Infinite products are handled by Product`Infinite, not here. */
    if (definite && imax->type == EXPR_SYMBOL && imax->data.symbol == SYM_Infinity)
        return NULL;

    /* A symbolic exponent (a^k) is Geometric territory and would hang Factor. */
    if (prod_has_symbolic_power(f, var)) return NULL;

    /* Split f into numerator / denominator polynomials in var. */
    Expr* tg  = prod_eval("Together",    (Expr*[]){ expr_copy(f) }, 1);
    Expr* num = prod_eval("Numerator",   (Expr*[]){ expr_copy(tg) }, 1);
    Expr* den = prod_eval("Denominator", (Expr*[]){ expr_copy(tg) }, 1);
    expr_free(tg);

    Expr *lnum = NULL, *lden = NULL;
    Expr **rn = NULL, **rd = NULL;
    int *mn = NULL, *md = NULL;
    size_t nn = 0, nd = 0;
    bool aln = false, ald = false;
    bool okn = prod_linear_factors(num, var, &lnum, &rn, &mn, &nn, &aln);
    bool okd = okn && prod_linear_factors(den, var, &lden, &rd, &md, &nd, &ald);
    expr_free(num); expr_free(den);

    if (!okn || !okd || !aln || !ald) {
        if (okn) { free_roots(rn, mn, nn); expr_free(lnum); }
        if (okd) { free_roots(rd, md, nd); expr_free(lden); }
        return NULL;
    }

    /* secondarg: definite -> count = imax-imin+1; indefinite -> var-1. */
    Expr* anchor;
    Expr* secondarg;
    if (definite) {
        anchor = expr_copy(imin);
        Expr* negmin = expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ expr_new_integer(-1), expr_copy(imin) }, 2);
        Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus),
                        (Expr*[]){ expr_copy(imax), negmin, expr_new_integer(1) }, 3);
        secondarg = evaluate(sum);
        expr_free(sum);
    } else {
        anchor = expr_new_integer(1);
        Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus),
                        (Expr*[]){ expr_copy(var), expr_new_integer(-1) }, 2);
        secondarg = evaluate(sum);
        expr_free(sum);
    }

    /* Assemble the factor list: leading constant ^ secondarg, then a Pochhammer
     * per root (numerator positive exponent, denominator negative). */
    size_t cap = nn + nd + 2, nf = 0;
    Expr** factors = malloc(sizeof(Expr*) * cap);

    /* c = lnum / lden ; contributes c^secondarg */
    Expr* c = prod_div(lnum, lden);
    factors[nf++] = expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ c, expr_copy(secondarg) }, 2);  /* adopts c */

    for (size_t i = 0; i < nn; i++) {
        Expr* a1 = anchor_minus_root(anchor, rn[i]);
        Expr* pk = poch(a1, secondarg);
        expr_free(a1);
        factors[nf++] = powi(pk, mn[i]);
        expr_free(pk);
    }
    for (size_t i = 0; i < nd; i++) {
        Expr* a1 = anchor_minus_root(anchor, rd[i]);
        Expr* pk = poch(a1, secondarg);
        expr_free(a1);
        factors[nf++] = powi(pk, -md[i]);
        expr_free(pk);
    }

    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times), factors, nf);
    free(factors);
    Expr* raw = evaluate(prod);
    expr_free(prod);
    Expr* out = poch1_to_fact(raw);   /* canonical n! for Pochhammer[1, n] */
    expr_free(raw);

    expr_free(anchor); expr_free(secondarg);
    expr_free(lnum); expr_free(lden);
    free_roots(rn, mn, nn);
    free_roots(rd, md, nd);
    return out;
}

void product_rational_init(void) {
    symtab_add_builtin("Product`Rational", builtin_product_rational);
    symtab_get_def("Product`Rational")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Product`Rational",
        "Product`Rational[f, i, imin, imax] gives the closed form of a product "
        "of a rational function f of i whose numerator and denominator factor "
        "into linear factors over the rationals, in terms of Pochhammer. "
        "Returns unevaluated to fall through when f is non-rational or has an "
        "irreducible quadratic-or-higher factor.");
}
