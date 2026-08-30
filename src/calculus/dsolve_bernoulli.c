/*
 * dsolve_bernoulli.c — DSolve`Bernoulli.
 *
 * Solves the Bernoulli equation  y'[x] == A(x) y + B(x) y^n  (n != 0, 1) by the
 * substitution v = y^(1-n), which linearises it:
 *
 *     v' == (1-n) A v + (1-n) B.
 *
 * The exponent n is recovered from F(x, Y) = A Y + B Y^n without assuming n is
 * an integer: with Q := F - Y F_Y = B(1-n) Y^n, the quantity Y Q_Y / Q equals n
 * exactly and is constant iff F has that single non-linear power.  Then
 * B = Q / ((1-n) Y^n) and A = (F - B Y^n)/Y; the linearised problem is handed to
 * dsolve_linear_factor_solve and the result raised to the 1/(1-n) power.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

/* head[a,b,c] built and evaluated; a,b,c consumed. */
static Expr* ev3(const char* head, Expr* a, Expr* b, Expr* c) {
    return eval_and_free(expr_new_function(expr_new_symbol(head), (Expr*[]){ a, b, c }, 3));
}
static Expr* powneg1(Expr* base) { /* base^-1; base consumed */
    return expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ base, expr_new_integer(-1) }, 2);
}

Expr** dsolve_bernoulli_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar = P->ind_names[0];

    Expr* F = dsolve_solve_top_derivative(P, 1);
    if (!F) return NULL;
    const char* Yn = intern_symbol("DSolve`Y");
    Expr* FY = ds_subst(F, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Yn));

    /* Q = FY - Y F_Y */
    Expr* Q = eval_and_free(ds_call2(SYM_Subtract, expr_copy(FY),
                 eval_and_free(ds_call2(SYM_Times, expr_new_symbol(Yn),
                                        ds_d(expr_copy(FY), expr_new_symbol(Yn))))));
    if (ds_is_zero(Q)) { expr_free(Q); expr_free(FY); return NULL; }   /* purely linear */

    /* n = Y Q_Y / Q, a constant != 0 */
    Expr* nexp = eval_and_free(ds_call2(SYM_Times,
                    eval_and_free(ds_call2(SYM_Times, expr_new_symbol(Yn),
                                           ds_d(expr_copy(Q), expr_new_symbol(Yn)))),
                    powneg1(expr_copy(Q))));
    if (!ds_free_of(nexp, Yn) || !ds_free_of(nexp, xvar) || ds_is_zero(nexp)) {
        expr_free(nexp); expr_free(Q); expr_free(FY); return NULL;
    }
    nexp = ds_simplify(nexp);   /* reduce the (now known constant) exponent */
    Expr* omn = eval_and_free(ds_call2(SYM_Subtract, expr_new_integer(1), expr_copy(nexp))); /* 1-n */
    if (ds_is_zero(omn)) { expr_free(omn); expr_free(nexp); expr_free(Q); expr_free(FY); return NULL; }

    Expr* Ynpow = eval_and_free(ds_call2(SYM_Power, expr_new_symbol(Yn), expr_copy(nexp))); /* Y^n */
    /* B = Q / ((1-n) Y^n) */
    Expr* B = ev3(SYM_Times, expr_copy(Q), powneg1(expr_copy(omn)), powneg1(expr_copy(Ynpow)));
    expr_free(Q);
    /* A = (FY - B Y^n) / Y */
    Expr* A = eval_and_free(ds_call2(SYM_Times,
                 eval_and_free(ds_call2(SYM_Subtract, expr_copy(FY),
                                        ds_call2(SYM_Times, expr_copy(B), expr_copy(Ynpow)))),
                 powneg1(expr_new_symbol(Yn))));

    bool ok = ds_free_of(A, Yn) && ds_free_of(B, Yn);
    if (ok) {
        Expr* recon = eval_and_free(ds_call2(SYM_Plus,
            ds_call2(SYM_Times, expr_copy(A), expr_new_symbol(Yn)),
            ds_call2(SYM_Times, expr_copy(B), expr_copy(Ynpow))));
        Expr* chk = eval_and_free(ds_call2(SYM_Subtract, expr_copy(FY), recon));
        ok = ds_is_zero(chk);
        expr_free(chk);
    }
    expr_free(FY); expr_free(Ynpow);
    if (!ok) { expr_free(A); expr_free(B); expr_free(omn); expr_free(nexp); return NULL; }

    /* v' + Pcoef v == Qcoef,  Pcoef = -(1-n) A,  Qcoef = (1-n) B */
    Expr* Pcoef = ev3(SYM_Times, expr_new_integer(-1), expr_copy(omn), A);  /* consumes A */
    Expr* Qcoef = eval_and_free(ds_call2(SYM_Times, expr_copy(omn), B));     /* consumes B */
    Expr* v = dsolve_linear_factor_solve(Pcoef, Qcoef, xvar);               /* consumes both */
    if (!v) { expr_free(omn); expr_free(nexp); return NULL; }

    /* y = v^(1/(1-n)) */
    Expr* body = eval_and_free(ds_call2(SYM_Power, v, powneg1(expr_copy(omn))));
    expr_free(omn); expr_free(nexp);

    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_bernoulli(Expr* res) {
    return dsolve_method_builtin(res, dsolve_bernoulli_try);
}

void dsolve_bernoulli_init(void) {
    symtab_add_builtin("DSolve`Bernoulli", builtin_dsolve_bernoulli);
    symtab_get_def("DSolve`Bernoulli")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("DSolve`Bernoulli",
        "DSolve`Bernoulli[eqn, y, x] solves y'[x] == A(x) y + B(x) y^n (n != 0, 1) "
        "via the substitution v = y^(1-n), which linearises the equation.");
}
