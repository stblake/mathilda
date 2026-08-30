/*
 * dsolve_pde1.c — first-order linear PDE via the method of characteristics.
 *
 * Solves  a u_{v1} + b u_{v2} + c(v1,v2) u == f(v1,v2)  with a, b CONSTANT.  The
 * characteristic invariant is xi = a v2 - b v1 (a (a v2 - b v1)_{v1} .. = 0), and
 * along a characteristic (parametrised by v1, with v2 = (xi + b v1)/a) the PDE
 * becomes the linear ODE  u_{v1} + (c/a) u = f/a, solved by the integrating
 * factor.  The homogeneous solution carries the arbitrary function C[1][xi]:
 *
 *     u = Exp[-Integrate[c/a, v1]] ( C[1][xi] + Integrate[Exp[..] f/a, v1] ),
 *
 * with xi -> a v2 - b v1 substituted back at the end.  Covers the transport
 * equation, `3 u_x + 5 u_y == x`, `u_x + 3 u_y + u == 1`, etc.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>

static Expr* ev1(const char* h, Expr* a) { return eval_and_free(ds_call1(h, a)); }
static Expr* mul(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Times, a, b)); }
static Expr* add(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Plus, a, b)); }
static Expr* inv(Expr* a) { return expr_new_function(expr_new_symbol(SYM_Power),
                                (Expr*[]){ a, expr_new_integer(-1) }, 2); }

/* Derivative[o1,o2][u][v1,v2] */
static Expr* pde_deriv(const char* u, int o1, int o2, const char* v1, const char* v2) {
    Expr* d = expr_new_function(expr_new_symbol(SYM_Derivative),
                  (Expr*[]){ expr_new_integer(o1), expr_new_integer(o2) }, 2);
    Expr* du = expr_new_function(d, (Expr*[]){ expr_new_symbol(u) }, 1);
    return expr_new_function(du, (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2);
}

Expr** dsolve_pde1_solve(DSolveProblem* P) {
    if (P->nfun != 1 || P->nind != 2 || P->neq != 1) return NULL;
    const char* uname = P->fun_names[0];
    const char* v1 = P->ind_names[0];
    const char* v2 = P->ind_names[1];

    /* algebraic residual: u_{v1} -> sUx, u_{v2} -> sUy, u -> sU */
    const char* sUx = intern_symbol("DSolve`pdeUx");
    const char* sUy = intern_symbol("DSolve`pdeUy");
    const char* sU  = intern_symbol("DSolve`pdeU");
    Expr* R = expr_copy(P->eq_residuals[0]);
    R = ds_subst(R, pde_deriv(uname, 1, 0, v1, v2), expr_new_symbol(sUx));
    R = ds_subst(R, pde_deriv(uname, 0, 1, v1, v2), expr_new_symbol(sUy));
    R = ds_subst(R, expr_new_function(expr_new_symbol(uname),
                     (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2), expr_new_symbol(sU));

    Expr* a = ds_d(expr_copy(R), expr_new_symbol(sUx));
    Expr* b = ds_d(expr_copy(R), expr_new_symbol(sUy));
    Expr* c = ds_d(expr_copy(R), expr_new_symbol(sU));
    /* linear + a,b constant */
    bool ok = ds_free_of(a, sUx) && ds_free_of(a, sUy) && ds_free_of(a, sU)
           && ds_free_of(b, sUx) && ds_free_of(b, sUy) && ds_free_of(b, sU)
           && ds_free_of(c, sUx) && ds_free_of(c, sUy) && ds_free_of(c, sU)
           && ds_free_of(a, v1) && ds_free_of(a, v2)
           && ds_free_of(b, v1) && ds_free_of(b, v2);
    if (!ok) { expr_free(R); expr_free(a); expr_free(b); expr_free(c); return NULL; }
    /* f = -(R with all s-symbols -> 0) */
    Expr* R0 = expr_copy(R);
    R0 = ds_subst(R0, expr_new_symbol(sUx), expr_new_integer(0));
    R0 = ds_subst(R0, expr_new_symbol(sUy), expr_new_integer(0));
    R0 = ds_subst(R0, expr_new_symbol(sU), expr_new_integer(0));
    Expr* f = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), R0));
    expr_free(R);

    /* need a != 0; if a == 0, swap the two variables (roles of a and b) */
    if (ds_is_zero(a)) {
        const char* t = v1; v1 = v2; v2 = t;
        Expr* tt = a; a = b; b = tt;
    }
    if (ds_is_zero(a)) { expr_free(a); expr_free(b); expr_free(c); expr_free(f); return NULL; }

    const char* xi = intern_symbol("DSolve`pdexi");
    /* v2 along a characteristic: (xi + b v1)/a */
    Expr* v2_along = mul(add(expr_new_symbol(xi), mul(expr_copy(b), expr_new_symbol(v1))),
                         inv(expr_copy(a)));
    Expr* cchar = ds_subst(expr_copy(c), expr_new_symbol(v2), expr_copy(v2_along));
    Expr* fchar = ds_subst(expr_copy(f), expr_new_symbol(v2), v2_along);

    /* integrating factor mu = Exp[Integrate[c/a, v1]] */
    Expr* Pint = ds_integrate(mul(cchar, inv(expr_copy(a))), expr_new_symbol(v1));
    if (ds_has_head(Pint, SYM_Integrate)) {
        expr_free(Pint); expr_free(fchar); expr_free(a); expr_free(b); expr_free(c); expr_free(f);
        return NULL;
    }
    Expr* mu = ev1("Exp", Pint);

    /* arbitrary function C[1][xi] */
    Expr* arbfn = expr_new_function(ds_const(1), (Expr*[]){ expr_new_symbol(xi) }, 1);
    Expr* inner = arbfn;
    bool have_f = !ds_is_zero(f);
    if (have_f) {
        Expr* integrand = mul(mul(expr_copy(mu), fchar), inv(expr_copy(a)));
        Expr* Qint = ds_integrate(integrand, expr_new_symbol(v1));
        if (ds_has_head(Qint, SYM_Integrate)) {
            expr_free(Qint); expr_free(inner); expr_free(mu);
            expr_free(a); expr_free(b); expr_free(c); expr_free(f);
            return NULL;
        }
        inner = add(inner, Qint);
    } else expr_free(fchar);

    /* body = mu^-1 inner, then xi -> a v2 - b v1 */
    Expr* body_xi = mul(inv(expr_copy(mu)), inner);
    expr_free(mu);
    Expr* xi_char = eval_and_free(ds_call2(SYM_Subtract,
                        ds_call2(SYM_Times, expr_copy(a), expr_new_symbol(v2)),
                        ds_call2(SYM_Times, expr_copy(b), expr_new_symbol(v1))));
    Expr* body = ds_subst(body_xi, expr_new_symbol(xi), xi_char);
    expr_free(a); expr_free(b); expr_free(c); expr_free(f);

    Expr** bodies = malloc(sizeof(Expr*));
    bodies[0] = body;
    (void)sU;
    return bodies;
}

void dsolve_pde1_init(void) {
    /* dispatched directly for PDEs; no backtick builtin yet */
}
