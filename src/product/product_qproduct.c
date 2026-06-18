/*
 * product_qproduct.c -- Product`QProduct (Stage 5): q-rational products.
 *
 * Recognises factors linear in q^k (q free of k, exponent exactly k):
 *   prod_{k=imin}^{imax} (alpha + beta q^k)
 *       = alpha^count * QPochhammer[(-beta/alpha) q^imin, q, count],
 *   count = imax - imin + 1,
 * so e.g.
 *   Product[1 - q^k,   {k, 1, n}]   = QPochhammer[q, q, n],
 *   Product[1 - a q^k, {k, 0, n-1}] = QPochhammer[a, q, n].
 * A reciprocal factor (denominator) contributes the reciprocal.
 *
 * CRITICAL: q^k is a symbolic power -- Together/Factor would loop on it
 * (project_together_factor_hang_exponential), so this stage walks the Times
 * factors structurally and never normalises across a q^k term.
 */

#include "product_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>

static bool is_head(const Expr* e, const char* sym) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == sym;
}

/* Analyse one factor `base` as alpha + beta * qbase^k (q-linear).
 * On success returns true and sets owned *alpha, *beta, *qbase; returns false
 * (nothing owned) if base is not q-linear with exponent exactly `var`. */
/* Recursively walk the (possibly nested) multiplicative factors of `e`.
 * Sets *qb to the single Power[qbase,k] base found (qbase free of var, exponent
 * exactly var); multiplies every other (var-free) factor into *bprod.  Returns
 * false if var appears outside such a single q^k, or two q^k factors appear. */
static bool split_qpow(Expr* e, Expr* var, Expr** qb, Expr** bprod) {
    if (is_head(e, SYM_Times)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!split_qpow(e->data.function.args[i], var, qb, bprod)) return false;
        return true;
    }
    if (is_head(e, SYM_Power) && e->data.function.arg_count == 2
            && e->data.function.args[1]->type == EXPR_SYMBOL
            && e->data.function.args[1]->data.symbol == var->data.symbol
            && prod_free_of(e->data.function.args[0], var)) {
        if (*qb) return false;                        /* two q^k factors */
        *qb = e->data.function.args[0];
        return true;
    }
    if (prod_free_of(e, var)) {
        Expr* p = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ *bprod, expr_copy(e) }, 2);   /* adopts *bprod */
        *bprod = evaluate(p); expr_free(p);
        return true;
    }
    return false;                                     /* var appears non-q-linearly */
}

static bool q_linear(Expr* base, Expr* var, Expr** alpha, Expr** beta, Expr** qbase) {
    /* terms of base (Plus args, or the single term) */
    size_t nt = is_head(base, SYM_Plus) ? base->data.function.arg_count : 1;
    Expr* asum = expr_new_integer(0);
    Expr* kterm = NULL;
    for (size_t i = 0; i < nt; i++) {
        Expr* t = is_head(base, SYM_Plus) ? base->data.function.args[i] : base;
        if (prod_free_of(t, var)) {
            Expr* s = expr_new_function(expr_new_symbol(SYM_Plus),
                          (Expr*[]){ asum, expr_copy(t) }, 2);   /* adopts asum */
            asum = evaluate(s); expr_free(s);
        } else if (!kterm) {
            kterm = t;
        } else {
            expr_free(asum); return false;            /* two k-terms: not q-linear */
        }
    }
    if (!kterm) { expr_free(asum); return false; }    /* free of k: not q-linear */

    /* kterm = beta * Power[qbase, k]; pull out the q^k subfactor (any nesting). */
    Expr* qb = NULL; Expr* bprod = expr_new_integer(1);
    if (!split_qpow(kterm, var, &qb, &bprod) || !qb) {
        expr_free(asum); expr_free(bprod); return false;
    }

    *alpha = asum;
    *beta = bprod;
    *qbase = expr_copy(qb);
    return true;
}

/* QPochhammer[(-beta/alpha) qbase^imin, qbase, count] * alpha^count, ^expo. */
static Expr* q_contribution(Expr* alpha, Expr* beta, Expr* qbase,
                            Expr* imin, Expr* count, int expo) {
    /* a = -beta/alpha */
    Expr* nbeta = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ expr_new_integer(-1), expr_copy(beta) }, 2);
    Expr* a = prod_div(nbeta, alpha);
    expr_free(nbeta);
    /* a * qbase^imin */
    Expr* qimin = expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ expr_copy(qbase), expr_copy(imin) }, 2);
    Expr* a0raw = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ a, qimin }, 2);   /* adopts a, qimin */
    Expr* a0 = evaluate(a0raw); expr_free(a0raw);
    /* QPochhammer[a0, qbase, count] */
    Expr* qp = expr_new_function(expr_new_symbol("QPochhammer"),
                   (Expr*[]){ a0, expr_copy(qbase), expr_copy(count) }, 3);  /* adopts a0 */
    /* alpha^count */
    Expr* ac = expr_new_function(expr_new_symbol(SYM_Power),
                   (Expr*[]){ expr_copy(alpha), expr_copy(count) }, 2);
    Expr* full = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ ac, qp }, 2);   /* adopts ac, qp */
    if (expo == 1) return full;
    Expr* p = expr_new_function(expr_new_symbol(SYM_Power),
                  (Expr*[]){ full, expr_new_integer(expo) }, 2);  /* adopts full */
    return p;
}

Expr* builtin_product_qproduct(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!product_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite) return NULL;                       /* indefinite q-products deferred */
    if (imax->type == EXPR_SYMBOL && imax->data.symbol == SYM_Infinity) return NULL;
    if (!prod_has_symbolic_power(f, var)) return NULL; /* no q^k -> not our job */

    /* count = imax - imin + 1 */
    Expr* nmin = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ expr_new_integer(-1), expr_copy(imin) }, 2);
    Expr* csum = expr_new_function(expr_new_symbol(SYM_Plus),
                     (Expr*[]){ expr_copy(imax), nmin, expr_new_integer(1) }, 3);
    Expr* count = evaluate(csum); expr_free(csum);

    bool is_times = is_head(f, SYM_Times);
    size_t fc = is_times ? f->data.function.arg_count : 1;
    size_t cap = fc + 1, nout = 0;
    Expr** out = malloc(sizeof(Expr*) * cap);
    bool bail = false;

    for (size_t i = 0; i < fc && !bail; i++) {
        Expr* fac = is_times ? f->data.function.args[i] : f;
        Expr* base = fac; int expo = 1;
        if (is_head(fac, SYM_Power) && fac->data.function.arg_count == 2
                && fac->data.function.args[1]->type == EXPR_INTEGER) {
            base = fac->data.function.args[0];
            expo = (int)fac->data.function.args[1]->data.integer;
        }
        if (prod_free_of(base, var)) {
            /* constant factor -> base^(expo*count) */
            Expr* ec = expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ expr_new_integer(expo), expr_copy(count) }, 2);
            Expr* p = expr_new_function(expr_new_symbol(SYM_Power),
                          (Expr*[]){ expr_copy(base), evaluate(ec) }, 2);
            expr_free(ec);
            out[nout++] = p;
            continue;
        }
        Expr *alpha = NULL, *beta = NULL, *qbase = NULL;
        if (!q_linear(base, var, &alpha, &beta, &qbase)) { bail = true; break; }
        out[nout++] = q_contribution(alpha, beta, qbase, imin, count, expo);
        expr_free(alpha); expr_free(beta); expr_free(qbase);
    }

    if (bail || nout == 0) {
        for (size_t i = 0; i < nout; i++) expr_free(out[i]);
        free(out); expr_free(count);
        return NULL;
    }
    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times), out, nout);
    free(out);
    Expr* result = evaluate(prod);   /* Times only -- never Factor/Together */
    expr_free(prod); expr_free(count);
    return result;
}

void product_qproduct_init(void) {
    symtab_add_builtin("Product`QProduct", builtin_product_qproduct);
    symtab_get_def("Product`QProduct")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Product`QProduct",
        "Product`QProduct[f, i, imin, imax] evaluates products of factors linear "
        "in q^i (q free of i) in terms of QPochhammer: "
        "prod (1 - a q^i) = QPochhammer[a q^imin, q, imax-imin+1]. Returns "
        "unevaluated otherwise.");
}
