/*
 * product_geometric.c -- Product`Geometric (Stage 3).
 *
 * Handles factors base^(p(k)) with base free of k and p polynomial in k, by
 * routing the EXPONENT through the shipped Sum:
 *   prod_{k} base^(p(k)) = base^( Sum[p(k), {k, imin, imax}] ).
 * The residual rational cofactor is handed to Product`Rational and multiplied
 * back.  Example:  Product[k 2^k, {k,1,n}] = n! * 2^(n(n+1)/2).
 *
 * CRITICAL (memory project_together_factor_hang_exponential): a Power[base,k]
 * factor must NEVER be passed to Together/Factor -- they loop forever on a
 * symbolic exponent.  This stage keeps every base^... factored and only ever
 * multiplies (Times) the assembled result, whose exponents are free of k.
 */

#include "product_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>

static bool is_power(const Expr* e) {
    return e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Power
        && e->data.function.arg_count == 2;
}

/* Sum[expo, {k, imin, imax}] (definite) or Sum[expo, k] (indefinite). */
static Expr* sum_exponent(Expr* expo, Expr* var, Expr* imin, Expr* imax,
                          bool definite) {
    Expr* r;
    if (definite) {
        /* Top-level Sum's definite surface form is Sum[f, {i, imin, imax}]. */
        Expr* spec = expr_new_function(expr_new_symbol(SYM_List),
                         (Expr*[]){ expr_copy(var), expr_copy(imin),
                                    expr_copy(imax) }, 3);
        Expr* args[2] = { expr_copy(expo), spec };
        r = prod_eval("Sum", args, 2);
    } else {
        Expr* args[2] = { expr_copy(expo), expr_copy(var) };
        r = prod_eval("Sum", args, 2);
    }
    /* If Sum could not close the exponent series, fall through. */
    if (r && r->type == EXPR_FUNCTION
          && r->data.function.head->type == EXPR_SYMBOL
          && strcmp(r->data.function.head->data.symbol, "Sum") == 0) {
        expr_free(r);
        return NULL;
    }
    return r;
}

Expr* builtin_product_geometric(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!product_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (definite && imax->type == EXPR_SYMBOL && imax->data.symbol == SYM_Infinity)
        return NULL;

    /* Geometric only engages when there is a genuine base^(involves k) factor;
     * pure rational inputs were already tried by Product`Rational. */
    if (!prod_has_symbolic_power(f, var)) return NULL;

    /* Factor list of f (Times args, or a single factor). */
    bool is_times = (f->type == EXPR_FUNCTION
                     && f->data.function.head->type == EXPR_SYMBOL
                     && f->data.function.head->data.symbol == SYM_Times);
    size_t fc = is_times ? f->data.function.arg_count : 1;

    size_t cap = fc + 1, ngeo = 0;
    Expr** out = malloc(sizeof(Expr*) * (cap + 1));   /* geometric base^E terms */
    Expr* cofactor = expr_new_integer(1);             /* rational remainder */
    bool bail = false;

    for (size_t i = 0; i < fc && !bail; i++) {
        Expr* fac = is_times ? f->data.function.args[i] : f;
        if (is_power(fac)) {
            Expr* base = fac->data.function.args[0];
            Expr* expo = fac->data.function.args[1];
            if (!prod_free_of(expo, var)) {
                /* exponent involves k */
                if (!prod_free_of(base, var)) { bail = true; break; }  /* k^k etc */
                Expr* E = sum_exponent(expo, var, imin, imax, definite);
                if (!E) { bail = true; break; }
                out[ngeo++] = expr_new_function(expr_new_symbol(SYM_Power),
                                  (Expr*[]){ expr_copy(base), E }, 2);  /* adopts E */
                continue;
            }
            /* exponent free of k -> ordinary rational factor (e.g. k^2) */
        }
        /* A non-power factor that still hides a symbolic power (e.g. 1 + 2^k)
         * is neither geometric nor rational -- bail. */
        if (prod_has_symbolic_power(fac, var)) { bail = true; break; }
        Expr* nc = expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ cofactor, expr_copy(fac) }, 2);  /* adopts cofactor */
        cofactor = evaluate(nc);
        expr_free(nc);
    }

    if (bail || ngeo == 0) {
        for (size_t i = 0; i < ngeo; i++) expr_free(out[i]);
        free(out);
        expr_free(cofactor);
        return NULL;
    }

    /* Rational cofactor closed form via Product`Rational (skip the trivial 1). */
    Expr* cof_result;
    bool cof_trivial = (cofactor->type == EXPR_INTEGER && cofactor->data.integer == 1);
    if (cof_trivial) {
        cof_result = cofactor;   /* == 1 */
    } else {
        if (definite) {
            Expr* args[4] = { cofactor, expr_copy(var), expr_copy(imin), expr_copy(imax) };
            cof_result = prod_eval("Product`Rational", args, 4);  /* adopts cofactor */
        } else {
            Expr* args[2] = { cofactor, expr_copy(var) };
            cof_result = prod_eval("Product`Rational", args, 2);  /* adopts cofactor */
        }
        if (cof_result && cof_result->type == EXPR_FUNCTION
              && cof_result->data.function.head->type == EXPR_SYMBOL
              && strcmp(cof_result->data.function.head->data.symbol,
                        "Product`Rational") == 0) {
            /* cofactor was not a rational the stage could close -> bail. */
            expr_free(cof_result);
            for (size_t i = 0; i < ngeo; i++) expr_free(out[i]);
            free(out);
            return NULL;
        }
    }

    out[ngeo++] = cof_result;
    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times), out, ngeo);
    free(out);
    Expr* result = evaluate(prod);   /* Times only -- never Factor/Together */
    expr_free(prod);
    return result;
}

void product_geometric_init(void) {
    symtab_add_builtin("Product`Geometric", builtin_product_geometric);
    symtab_get_def("Product`Geometric")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Product`Geometric",
        "Product`Geometric[f, i, imin, imax] evaluates products containing "
        "base^(p(i)) factors (base free of i) by routing the exponent through "
        "Sum: prod base^p(i) = base^Sum[p(i),{i,imin,imax}], multiplying any "
        "rational cofactor via Product`Rational. Returns unevaluated otherwise.");
}
