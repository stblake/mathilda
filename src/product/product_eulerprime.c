/*
 * product_eulerprime.c -- Product`EulerPrime: Euler products over the primes.
 *
 * Recognises bodies that depend on the index only through Prime[k] and maps the
 * Euler factor to a zeta / Dirichlet-L value:
 *
 *   prod_k 1/(1 - Prime[k]^-s)                 = Zeta[s]
 *       (even s auto-reduces: Zeta[2]=Pi^2/6, Zeta[4]=Pi^4/90)
 *   prod_k 1/(1 - (-1)^((Prime[k]-1)/2)/Prime[k]) = L(1, chi_4) = Pi/4
 *       (chi_4 is the non-principal character mod 4; chi_4(2)=0, so k>=1 and
 *        k>=2 give the same product)
 *
 * Recognition: substitute Prime[k] -> p (a fresh formal symbol), require the
 * result free of k, then match the reduced Euler factor r = Simplify[1 - 1/pf]
 * against the templates p^-s (zeta) and (-1)^((p-1)/2) p^-1 (chi_4).
 *
 * chi_4 is hard-coded as a single character; a general Dirichlet-L engine is a
 * much larger project and is intentionally out of scope -- any other character
 * falls through (NULL).
 *
 * Memory contract: takes ownership of res but must not free it; returns an
 * owned closed form or NULL to fall through.
 */

#include "product_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>

static bool is_inf_sym(const Expr* e) {
    return e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Infinity;
}

static bool is_power(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Power
        && e->data.function.arg_count == 2;
}

/* PossibleZeroQ[a - b] == True. */
static bool same_val(Expr* a, Expr* b) {
    Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(b) }, 2);
    Expr* diff = expr_new_function(expr_new_symbol(SYM_Plus),
                     (Expr*[]){ expr_copy(a), neg }, 2);          /* adopts neg */
    Expr* r = prod_eval("PossibleZeroQ", (Expr*[]){ diff }, 1);   /* adopts diff */
    bool yes = (r && r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_True);
    if (r) expr_free(r);
    return yes;
}

/* Replace every Prime[var] in e by the formal symbol p (ReplaceAll). */
static Expr* subst_prime(Expr* e, Expr* var, Expr* p) {
    Expr* pk = expr_new_function(expr_new_symbol("Prime"),
                   (Expr*[]){ expr_copy(var) }, 1);
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
                     (Expr*[]){ pk, expr_copy(p) }, 2);           /* adopts pk */
    Expr* args[2] = { expr_copy(e), rule };                      /* adopts rule */
    return prod_eval("ReplaceAll", args, 2);
}

Expr* builtin_product_eulerprime(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!product_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite || !is_inf_sym(imax)) return NULL;
    /* imin must range over all primes (k>=1) or all odd primes (k>=2). */
    if (!(imin->type == EXPR_INTEGER
          && (imin->data.integer == 1 || imin->data.integer == 2))) return NULL;

    /* Body must depend on k only through Prime[k]. */
    if (prod_free_of(f, var)) return NULL;
    Expr* p = expr_new_symbol("Product`EulerPrime`p");
    Expr* pf = subst_prime(f, var, p);
    if (!prod_free_of(pf, var)) { expr_free(pf); expr_free(p); return NULL; }

    /* Reduced Euler factor r = Simplify[1 - 1/pf]. */
    Expr* invpf = expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ expr_copy(pf), expr_new_integer(-1) }, 2);
    Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), invpf }, 2);   /* adopts invpf */
    Expr* one_minus = expr_new_function(expr_new_symbol(SYM_Plus),
                          (Expr*[]){ expr_new_integer(1), neg }, 2); /* adopts neg */
    Expr* r = prod_eval("Simplify", (Expr*[]){ one_minus }, 1);      /* adopts */
    expr_free(pf);

    Expr* out = NULL;

    /* Zeta template: r == p^(-s), s free of p. */
    if (is_power(r) && expr_eq(r->data.function.args[0], p)) {
        Expr* e = r->data.function.args[1];
        Expr* s = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ expr_new_integer(-1), expr_copy(e) }, 2);
        Expr* se = evaluate(s);
        expr_free(s);
        if (prod_free_of(se, p)) {
            /* Reject a provably-divergent exponent s <= 1. */
            Expr* gt = expr_new_function(expr_new_symbol("Greater"),
                           (Expr*[]){ expr_copy(se), expr_new_integer(1) }, 2);
            Expr* gte = evaluate(gt);
            expr_free(gt);
            bool diverges = (gte && gte->type == EXPR_SYMBOL
                             && gte->data.symbol.name == SYM_False);
            if (gte) expr_free(gte);
            if (!diverges)
                out = prod_eval("Zeta", (Expr*[]){ expr_copy(se) }, 1);
        }
        expr_free(se);
    }

    /* chi_4 template: r == (-1)^((p-1)/2) * p^(-1). */
    if (!out) {
        Expr* pm1 = expr_new_function(expr_new_symbol(SYM_Plus),
                        (Expr*[]){ expr_new_integer(-1), expr_copy(p) }, 2);
        Expr* halfexp = expr_new_function(expr_new_symbol(SYM_Times),
                            (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Rational),
                                          (Expr*[]){ expr_new_integer(1),
                                                     expr_new_integer(2) }, 2),
                                       pm1 }, 2);                 /* adopts pm1 */
        Expr* chi = expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ expr_new_integer(-1), halfexp }, 2); /* adopts halfexp */
        Expr* invp = expr_new_function(expr_new_symbol(SYM_Power),
                         (Expr*[]){ expr_copy(p), expr_new_integer(-1) }, 2);
        Expr* expected = expr_new_function(expr_new_symbol(SYM_Times),
                             (Expr*[]){ chi, invp }, 2);          /* adopts chi, invp */
        Expr* expected_e = evaluate(expected);
        expr_free(expected);
        if (same_val(r, expected_e)) {
            /* L(1, chi_4) = Pi/4. */
            out = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Rational),
                                     (Expr*[]){ expr_new_integer(1),
                                                expr_new_integer(4) }, 2),
                                 expr_new_symbol("Pi") }, 2);
        }
        expr_free(expected_e);
    }

    expr_free(r);
    expr_free(p);
    return out;
}

void product_eulerprime_init(void) {
    symtab_add_builtin("Product`EulerPrime", builtin_product_eulerprime);
    symtab_get_def("Product`EulerPrime")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Product`EulerPrime",
        "Product`EulerPrime[f, i, imin, Infinity] evaluates Euler products over "
        "the primes: prod 1/(1-Prime[i]^-s) = Zeta[s], and the chi_4 product "
        "prod 1/(1-(-1)^((Prime[i]-1)/2)/Prime[i]) = Pi/4. Returns unevaluated "
        "otherwise.");
}
