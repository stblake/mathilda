/*
 * dsolve_decouple.c — DSolve`DecoupleSystem.
 *
 * Solves a system in which each equation involves only one dependent function
 * (a decoupled / triangular-by-inspection system) by recursing into the scalar
 * engine for each function and renumbering the generated constants so they do
 * not collide.  Handles variable-coefficient components (e.g. y' == x^2 y) that
 * the constant-coefficient system solver cannot.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include "../internal.h"
#include <stdlib.h>

/* Renumber C[1..m] in `body` to C[offset+1..offset+m]; body consumed. */
static Expr* renumber(Expr* body, int m, int* offset) {
    if (m <= 0) return body;
    Expr** rules = malloc((size_t)m * sizeof(Expr*));
    for (int j = 1; j <= m; j++)
        rules[j - 1] = expr_new_function(expr_new_symbol(SYM_Rule),
                           (Expr*[]){ ds_const(j), ds_const(*offset + j) }, 2);
    Expr* rl = expr_new_function(expr_new_symbol(SYM_List), rules, (size_t)m);
    free(rules);
    *offset += m;
    return eval_and_free(internal_replace_all((Expr*[]){ body, rl }, 2));
}

/* Extract the body from a DSolve result {{fname -> Function[{x}, body]}}. */
static Expr* extract_body(Expr* r, const char* fname) {
    if (!head_is(r, SYM_List) || r->data.function.arg_count < 1) return NULL;
    Expr* inner = r->data.function.args[0];
    if (!head_is(inner, SYM_List)) return NULL;
    for (size_t k = 0; k < inner->data.function.arg_count; k++) {
        Expr* rule = inner->data.function.args[k];
        if (head_is(rule, SYM_Rule) && rule->data.function.arg_count == 2) {
            Expr* lhs = rule->data.function.args[0];
            if (lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == fname) {
                Expr* rhs = rule->data.function.args[1];
                if (head_is(rhs, SYM_Function) && rhs->data.function.arg_count == 2)
                    return expr_copy(rhs->data.function.args[1]);
                return expr_copy(rhs);
            }
        }
    }
    return NULL;
}

Expr** dsolve_decouple_solve(DSolveProblem* P) {
    size_t n = P->nfun;
    const char* xvar = P->ind_names[0];

    /* group each equation by the single function it mentions */
    Expr*** groups = calloc(n, sizeof(Expr**));
    size_t* gcount = calloc(n, sizeof(size_t));
    bool decoupled = true;
    for (size_t e = 0; e < P->neq && decoupled; e++) {
        long which = -1; int cnt = 0;
        for (size_t j = 0; j < n; j++)
            if (ds_contains(P->eq_residuals[e], P->fun_names[j])) { which = (long)j; cnt++; }
        if (cnt != 1) { decoupled = false; break; }
        groups[which] = realloc(groups[which], (gcount[which] + 1) * sizeof(Expr*));
        groups[which][gcount[which]++] = P->eq_residuals[e];
    }
    for (size_t i = 0; i < n && decoupled; i++) if (gcount[i] == 0) decoupled = false;
    if (!decoupled) {
        for (size_t i = 0; i < n; i++) free(groups[i]);
        free(groups); free(gcount);
        return NULL;
    }

    Expr** bodies = calloc(n, sizeof(Expr*));
    int offset = 0;
    bool ok = true;
    for (size_t i = 0; i < n && ok; i++) {
        Expr** eqs = malloc(gcount[i] * sizeof(Expr*));
        for (size_t k = 0; k < gcount[i]; k++)
            eqs[k] = expr_new_function(expr_new_symbol(SYM_Equal),
                         (Expr*[]){ expr_copy(groups[i][k]), expr_new_integer(0) }, 2);
        Expr* eqarg = (gcount[i] == 1) ? eqs[0]
                       : expr_new_function(expr_new_symbol(SYM_List), eqs, gcount[i]);
        free(eqs);
        Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                         (Expr*[]){ eqarg, expr_new_symbol(P->fun_names[i]), expr_new_symbol(xvar) }, 3);
        Expr* r = eval_and_free(call);
        Expr* body = extract_body(r, P->fun_names[i]);
        expr_free(r);
        if (!body) { ok = false; break; }
        bodies[i] = renumber(body, P->max_order[i], &offset);
    }

    for (size_t i = 0; i < n; i++) free(groups[i]);
    free(groups); free(gcount);
    if (!ok) {
        for (size_t i = 0; i < n; i++) if (bodies[i]) expr_free(bodies[i]);
        free(bodies);
        return NULL;
    }
    return bodies;
}

void dsolve_decouple_init(void) {
    /* system methods are dispatched directly for nfun>1; no backtick builtin */
}
