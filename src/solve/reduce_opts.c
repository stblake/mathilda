/*
 * reduce_opts.c  --  see reduce_opts.h.
 */
#include "reduce_opts.h"

#include <stdlib.h>
#include <string.h>

#include "sym_names.h"

void reduce_opts_default(ReduceOpts* o) {
    o->poly.cubics_radical   = false;   /* Cubics           -> False    */
    o->poly.quartics_radical = false;   /* Quartics         -> False    */
    o->modulus               = NULL;    /* Modulus          -> 0        */
    o->param_head            = "C";     /* GeneratedParameters -> C     */
    o->backsub               = false;   /* Backsubstitution -> False    */
    o->working_precision     = NULL;    /* WorkingPrecision -> Infinity */
    o->method                = NULL;    /* Method           -> Automatic*/
}

/* Rule[Symbol(name), val] -- consumes `val`. */
static Expr* mk_rule(const char* name, Expr* val) {
    return expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_new_symbol(name), val }, 2);
}

Expr* reduce_opts_build_solve(Expr** base_args, int nbase, const ReduceOpts* o) {
    Expr* extra[3];
    int ne = 0;
    if (o && o->poly.cubics_radical)
        extra[ne++] = mk_rule(SYM_Cubics, expr_new_symbol(SYM_True));
    if (o && o->poly.quartics_radical)
        extra[ne++] = mk_rule(SYM_Quartics, expr_new_symbol(SYM_True));
    if (o && o->param_head && strcmp(o->param_head, "C") != 0)
        extra[ne++] = mk_rule(SYM_GeneratedParameters, expr_new_symbol(o->param_head));

    int n = nbase + ne;
    Expr** args = malloc((size_t)n * sizeof(Expr*));
    for (int i = 0; i < nbase; i++) args[i] = base_args[i];
    for (int i = 0; i < ne; i++) args[nbase + i] = extra[i];
    Expr* call = expr_new_function(expr_new_symbol(SYM_Solve), args, (size_t)n);
    free(args);
    return call;
}

bool reduce_opts_wp_digits(const ReduceOpts* o, double* digits) {
    if (!o || !o->working_precision) return false;
    const Expr* wp = o->working_precision;
    if (wp->type == EXPR_INTEGER) {
        if (wp->data.integer <= 0) return false;
        *digits = (double)wp->data.integer;
        return true;
    }
    if (wp->type == EXPR_REAL) {
        if (wp->data.real <= 0.0) return false;
        *digits = wp->data.real;
        return true;
    }
    /* Infinity, MachinePrecision, or any symbolic value: exact-first. */
    return false;
}
