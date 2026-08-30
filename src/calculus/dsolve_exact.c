/*
 * dsolve_exact.c — DSolve`Exact.
 *
 * Solves  M(x,y) + N(x,y) y' == 0  when it is exact (M_y == N_x), or can be made
 * exact by an integrating factor mu(x) or mu(y):
 *   - (M_y - N_x)/N free of y  ->  mu(x) = Exp[Integrate[(M_y-N_x)/N, x]]
 *   - (N_x - M_y)/M free of x  ->  mu(y) = Exp[Integrate[(N_x-M_y)/M, y]]
 * The potential F with F_x = M, F_y = N is F = Integrate[M, x] + g(y) where
 * g'(y) = N - d/dy Integrate[M, x]; the solution is F(x,y) == C[1], solved for y.
 *
 * Works on the algebraic residual R(x, Y, p) (p = y'), which must be linear in p:
 * N = R_p, M = R|p=0.  Every returned branch is still back-substitution verified
 * by the substrate, so an imperfect potential construction cannot ship a wrong
 * answer — it is dropped.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

static Expr* pw_inv(Expr* base) {
    return expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ base, expr_new_integer(-1) }, 2);
}

Expr** dsolve_exact_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* xvar = P->ind_names[0];
    const char* Yn = intern_symbol("DSolve`Y");
    const char* Pn = intern_symbol("DSolve`p");

    Expr* R = dsolve_algebraic_residual(P, Yn, Pn);
    if (!R) return NULL;

    /* R = M + N p, N = R_p free of p, M = R|p=0 */
    Expr* N = ds_d(expr_copy(R), expr_new_symbol(Pn));
    if (!ds_free_of(N, Pn)) { expr_free(N); expr_free(R); return NULL; }
    Expr* M = ds_subst(expr_copy(R), expr_new_symbol(Pn), expr_new_integer(0));
    Expr* recon = eval_and_free(ds_call2(SYM_Plus, expr_copy(M),
                        ds_call2(SYM_Times, expr_copy(N), expr_new_symbol(Pn))));
    Expr* lin = eval_and_free(ds_call2(SYM_Subtract, expr_copy(R), recon));
    bool linear = ds_is_zero(lin);
    expr_free(lin); expr_free(R);
    if (!linear) { expr_free(M); expr_free(N); return NULL; }

    /* diff = M_y - N_x */
    Expr* My = ds_d(expr_copy(M), expr_new_symbol(Yn));
    Expr* Nx = ds_d(expr_copy(N), expr_new_symbol(xvar));
    Expr* diff = eval_and_free(ds_call2(SYM_Subtract, My, Nx));

    Expr* Mu = NULL;
    if (ds_is_zero(diff)) {
        Mu = expr_new_integer(1);
    } else {
        /* mu(x): (M_y - N_x)/N free of Y */
        Expr* r1 = eval_and_free(ds_call2(SYM_Times, expr_copy(diff), pw_inv(expr_copy(N))));
        if (ds_free_of(r1, Yn)) {
            Expr* r1int = ds_integrate(r1, expr_new_symbol(xvar));
            if (!ds_has_head(r1int, SYM_Integrate)) Mu = eval_and_free(ds_call1("Exp", r1int));
            else expr_free(r1int);
        } else {
            expr_free(r1);
            /* mu(y): (N_x - M_y)/M free of x */
            Expr* r2 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
                expr_new_integer(-1), expr_copy(diff), pw_inv(expr_copy(M)) }, 3));
            if (ds_free_of(r2, xvar)) {
                Expr* r2int = ds_integrate(r2, expr_new_symbol(Yn));
                if (!ds_has_head(r2int, SYM_Integrate)) Mu = eval_and_free(ds_call1("Exp", r2int));
                else expr_free(r2int);
            } else expr_free(r2);
        }
    }
    expr_free(diff);
    if (!Mu) { expr_free(M); expr_free(N); return NULL; }

    /* Mm = Mu M, Nn = Mu N */
    Expr* Mm = eval_and_free(ds_call2(SYM_Times, expr_copy(Mu), M));   /* consumes M */
    Expr* Nn = eval_and_free(ds_call2(SYM_Times, expr_copy(Mu), N));   /* consumes N */
    expr_free(Mu);

    /* F = Integrate[Mm, x] + g(Y),  g'(Y) = Nn - d/dY Integrate[Mm, x] */
    Expr* Fx = ds_integrate(Mm, expr_new_symbol(xvar));
    if (ds_has_head(Fx, SYM_Integrate)) { expr_free(Fx); expr_free(Nn); return NULL; }
    Expr* gp = eval_and_free(ds_call2(SYM_Subtract, Nn, ds_d(expr_copy(Fx), expr_new_symbol(Yn))));
    Expr* g = ds_integrate(gp, expr_new_symbol(Yn));
    if (ds_has_head(g, SYM_Integrate)) { expr_free(g); expr_free(Fx); return NULL; }
    Expr* Fpot = eval_and_free(ds_call2(SYM_Plus, Fx, g));

    /* solution F(x, Y) == C[1], solved for Y */
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ Fpot, ds_const(1) }, 2);
    Expr* solres = ds_solve(eq, expr_new_symbol(Yn));
    size_t nb = 0;
    Expr** bodies = dsolve_extract_solutions(solres, Yn, &nb);
    if (solres) expr_free(solres);
    if (!bodies) return NULL;                 /* only explicit solutions in M1 */
    *nbranch = nb;
    return bodies;
}

static Expr* builtin_dsolve_exact(Expr* res) {
    return dsolve_method_builtin(res, dsolve_exact_try);
}

void dsolve_exact_init(void) {
    symtab_add_builtin("DSolve`Exact", builtin_dsolve_exact);
    symtab_get_def("DSolve`Exact")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("DSolve`Exact",
        "DSolve`Exact[eqn, y, x] solves M(x,y) + N(x,y) y' == 0 when exact "
        "(M_y == N_x), or made exact by an integrating factor mu(x) or mu(y).");
}
