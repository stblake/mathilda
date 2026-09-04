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

/* The per-function solve/renumber helpers (extract body, renumber constants)
 * are shared with DSolve`TriangularSystem; they live in dsolve_common.c as
 * dsolve_extract_system_body / dsolve_renumber_constants. */

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
        Expr* body = dsolve_extract_system_body(r, P->fun_names[i]);
        expr_free(r);
        if (!body) { ok = false; break; }
        bodies[i] = dsolve_renumber_constants(body, P->max_order[i], &offset);
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

static Expr* builtin_dsolve_decouple(Expr* res) {
    return dsolve_method_builtin_system(res, dsolve_decouple_solve);
}

void dsolve_decouple_init(void) {
    symtab_add_builtin("DSolve`DecoupleSystem", builtin_dsolve_decouple);
    symtab_get_def("DSolve`DecoupleSystem")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`DecoupleSystem",
        "DSolve`DecoupleSystem[{eqns}, {y1, y2, ...}, x] solves a first-order system in "
        "which each equation involves only ONE dependent function (the dependency graph "
        "has no cross-edges): it solves each function independently by recursing into the "
        "scalar cascade and renumbers the generated constants. Handles variable "
        "coefficients (e.g. y'[x] == x^2 y[x]). Declines a genuinely coupled system. The "
        "cheapest system method, tried first for nfun>1.");
}
