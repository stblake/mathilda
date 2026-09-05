/*
 * dsolve_autosys.c — DSolve`AutonomousSystem.
 *
 * Two-dimensional AUTONOMOUS systems  x' == f(x,y),  y' == g(x,y)  (f, g free of
 * the independent variable t), including NONLINEAR ones, by the phase-plane
 * reduction:
 *
 *   1. Eliminate t: the orbit satisfies  dy/dx == g(x,y)/f(x,y), a single scalar
 *      first-order ODE.  Solve it for the first integral / y = Y(x, C1).
 *   2. Reconstruct t: along the orbit  x' == f(x, Y(x)), a scalar autonomous
 *      (separable) ODE, giving x(t, C2); then y(t) = Y(x(t)).
 *
 * Each stage reuses the full scalar DSolve cascade, so any orbit the scalar
 * engine can integrate (separable, linear, Bernoulli, linearizable, …) yields
 * the system solution.  Examples:
 *   {x'=y, y'=y^2/x}      -> orbit y = C x        -> x = C2 e^{C t}
 *   {x'=-1/y, y'=1/x}     -> orbit x y = C        -> x = C2 e^{-t/C}
 *   {x'=x/y, y'=y/x}      -> orbit 1/x-1/y = C    -> x = 1/C + C2 e^{-C t}
 *   {x'=1/y, y'=1/x}      -> orbit y = C x        -> x = Sqrt[2 t/C + C2]
 *
 * The two arbitrary constants (C1 the orbit, C2 the time-translation) are the
 * right count for a 2-D first-order system.  The orbit constant is renumbered to
 * C[2] up front so the reconstruction's fresh C[1] cannot collide with it.
 *
 * Runs last in the system cascade (after the linear-system solvers), gated to a
 * genuinely autonomous 2-function system; verification is the substrate's.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../eval.h"
#include <stdlib.h>

/* Solve the 2-equation system for the two leading derivatives; returns f,g (the
 * RHS of u1' , u2') as owned exprs in the plain symbols X, Y (u_i[t] -> Xi), or
 * false.  X, Y are the interned plain-symbol names for the two functions. */
static bool autosys_extract_fg(DSolveProblem* P, const char* X, const char* Y,
                               Expr** fout, Expr** gout) {
    const char* xvar = P->ind_names[0];
    const char* u1 = P->fun_names[0];
    const char* u2 = P->fun_names[1];
    Expr* D1 = ds_make_funcapp(u1, 1, xvar);
    Expr* D2 = ds_make_funcapp(u2, 1, xvar);
    /* Solve[{R1==0, R2==0}, {u1'[t], u2'[t]}] */
    Expr* eqs = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){
        expr_new_function(expr_new_symbol(SYM_Equal),
            (Expr*[]){ expr_copy(P->eq_residuals[0]), expr_new_integer(0) }, 2),
        expr_new_function(expr_new_symbol(SYM_Equal),
            (Expr*[]){ expr_copy(P->eq_residuals[1]), expr_new_integer(0) }, 2) }, 2);
    Expr* vars = expr_new_function(expr_new_symbol(SYM_List),
                     (Expr*[]){ expr_copy(D1), expr_copy(D2) }, 2);
    Expr* sol = eval_and_free(ds_call2("Solve", eqs, vars));
    bool ok = false;
    Expr* f = NULL; Expr* g = NULL;
    if (sol && sol->type == EXPR_FUNCTION && sol->data.function.arg_count >= 1
        && sol->data.function.head->type == EXPR_SYMBOL
        && sol->data.function.head->data.symbol.name == SYM_List) {
        Expr* br = sol->data.function.args[0];   /* first branch: list of rules */
        if (br->type == EXPR_FUNCTION && br->data.function.arg_count == 2) {
            /* each rule lhs -> rhs; match D1/D2 */
            for (int pass = 0; pass < 2; pass++) {
                Expr* rule = br->data.function.args[pass];
                if (rule->type != EXPR_FUNCTION || rule->data.function.arg_count != 2) continue;
                Expr* lhs = rule->data.function.args[0];
                Expr* rhs = rule->data.function.args[1];
                if (expr_eq(lhs, D1)) f = expr_copy(rhs);
                else if (expr_eq(lhs, D2)) g = expr_copy(rhs);
            }
            ok = (f && g);
        }
    }
    if (sol) expr_free(sol);
    expr_free(D1); expr_free(D2);
    if (!ok) { if (f) expr_free(f); if (g) expr_free(g); return false; }
    /* to symbols: u1[t] -> X, u2[t] -> Y */
    f = ds_subst(f, ds_make_funcapp(u1, 0, xvar), expr_new_symbol(X));
    f = ds_subst(f, ds_make_funcapp(u2, 0, xvar), expr_new_symbol(Y));
    g = ds_subst(g, ds_make_funcapp(u1, 0, xvar), expr_new_symbol(X));
    g = ds_subst(g, ds_make_funcapp(u2, 0, xvar), expr_new_symbol(Y));
    *fout = f; *gout = g;
    return true;
}

/* Does `e` contain a radical — a Power with a non-integer exponent?  The
 * reconstruction ODE x' == f(x, Y(x)) built from a radical orbit Y (e.g.
 * Y = Sqrt[x^2 + C] from x'=y/(x-y), y'=x/(x-y)) leads to an implicit /
 * nonelementary quadrature the scalar solver cannot invert to explicit x(t)
 * (and churns on to the deadline).  Restricting reconstruction to rational
 * orbits keeps the closed-form cases (y = C x, C/x, x/(1-Cx), …) and declines
 * the radical ones cheaply.  This is a speed gate, not a safety one — the
 * risch_squarefree_t NULL-deref such an integrand once hit is now fixed. */
static bool autosys_has_radical(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Power
        && e->data.function.arg_count == 2
        && e->data.function.args[1]->type != EXPR_INTEGER)
        return true;
    if (autosys_has_radical(h)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (autosys_has_radical(e->data.function.args[i])) return true;
    return false;
}

Expr** dsolve_autosys_solve(DSolveProblem* P) {
    if (P->nfun != 2 || P->neq != 2) return NULL;
    const char* xvar = P->ind_names[0];
    const char* u1 = P->fun_names[0];
    const char* u2 = P->fun_names[1];

    const char* Xs = intern_symbol("DSolve`asX");
    const char* Ys = intern_symbol("DSolve`asY");
    const char* Yf = intern_symbol("DSolve`asYf");   /* orbit dependent funcapp */

    Expr* f = NULL; Expr* g = NULL;
    if (!autosys_extract_fg(P, Xs, Ys, &f, &g)) return NULL;

    /* autonomous: f, g free of the independent variable */
    if (!ds_free_of(f, xvar) || !ds_free_of(g, xvar)) { expr_free(f); expr_free(g); return NULL; }
    /* need f != 0 to form dy/dx */
    if (ds_is_zero(f)) { expr_free(f); expr_free(g); return NULL; }

    Expr** bodies = NULL;

    /* orbit: dY/dX == g/f, with X the independent, Y the dependent */
    Expr* gf = ds_simplify(eval_and_free(ds_call2(SYM_Times, expr_copy(g),
                   ds_call2(SYM_Power, expr_copy(f), expr_new_integer(-1)))));
    /* substitute the plain Y-symbol by the funcapp Yf[X] for the scalar solve */
    Expr* gfY = ds_subst(gf, expr_new_symbol(Ys), ds_make_funcapp(Yf, 0, Xs));
    Expr* orbeq = expr_new_function(expr_new_symbol(SYM_Equal),
                      (Expr*[]){ ds_make_funcapp(Yf, 1, Xs), gfY }, 2);
    Expr* orbcall = expr_new_function(expr_new_symbol(SYM_DSolve),
                        (Expr*[]){ orbeq, ds_make_funcapp(Yf, 0, Xs),
                                   expr_new_symbol(Xs) }, 3);
    Expr* orbres = eval_and_free(orbcall);
    size_t nb = 0;
    Expr** Yexprs = dsolve_extract_applied_bodies(orbres, Yf, &nb);
    if (orbres) expr_free(orbres);

    if (Yexprs && nb >= 1) {
        /* take the first orbit branch; renumber its constant C[1] -> C[2] */
        int off = 1;
        Expr* Yexpr = dsolve_renumber_constants(expr_copy(Yexprs[0]), 1, &off); /* C[k]->C[k+1] */
        /* reconstruct x(t): x1'[t] == f( x1[t], Y(x1[t]) ) */
        Expr* fX = ds_subst(expr_copy(f), expr_new_symbol(Ys), expr_copy(Yexpr)); /* f in X only */
        Expr* rhs = ds_subst(fX, expr_new_symbol(Xs), ds_make_funcapp(u1, 0, xvar));
        /* Only attempt the reconstruction for a rational orbit + rational field.
         * A radical orbit gives an implicit / nonelementary reconstruction that
         * the scalar solver either cannot invert to explicit x(t) or churns on
         * to the deadline (measured: #66 times out, #71 declines after ~7 s).
         * Decline cheaply up front instead — this is a speed gate, not a safety
         * one: the underlying risch_squarefree_t NULL-deref that a radical
         * reconstruction used to hit is fixed (see risch_canonical.c), so
         * removing this gate is now safe but merely slow. */
        if (autosys_has_radical(Yexpr) || autosys_has_radical(rhs)) {
            expr_free(rhs); expr_free(Yexpr);
            for (size_t i = 0; i < nb; i++) expr_free(Yexprs[i]);
            free(Yexprs);
            expr_free(f); expr_free(g);
            return NULL;
        }
        Expr* xeq = expr_new_function(expr_new_symbol(SYM_Equal),
                        (Expr*[]){ ds_make_funcapp(u1, 1, xvar), rhs }, 2);
        Expr* xcall = expr_new_function(expr_new_symbol(SYM_DSolve),
                          (Expr*[]){ xeq, ds_make_funcapp(u1, 0, xvar),
                                     expr_new_symbol(xvar) }, 3);
        Expr* xres = eval_and_free(xcall);
        size_t nx = 0;
        Expr** xbodies = dsolve_extract_applied_bodies(xres, u1, &nx);
        if (xres) expr_free(xres);

        if (xbodies && nx >= 1) {
            Expr* x1 = expr_copy(xbodies[0]);   /* x(t), carries its own C[1] */
            /* y(t) = Yexpr(X -> x1) */
            Expr* y1 = ds_simplify(ds_subst(expr_copy(Yexpr), expr_new_symbol(Xs), expr_copy(x1)));
            bodies = malloc(2 * sizeof(Expr*));
            bodies[0] = x1;
            bodies[1] = y1;
        }
        if (xbodies) { for (size_t i = 0; i < nx; i++) expr_free(xbodies[i]); free(xbodies); }
        expr_free(Yexpr);
    }
    if (Yexprs) { for (size_t i = 0; i < nb; i++) expr_free(Yexprs[i]); free(Yexprs); }

    expr_free(f); expr_free(g);
    (void)u2;
    return bodies;
}

static Expr* builtin_dsolve_autosys(Expr* res) {
    return dsolve_method_builtin_system(res, dsolve_autosys_solve);
}

void dsolve_autosys_init(void) {
    symtab_add_builtin("DSolve`AutonomousSystem", builtin_dsolve_autosys);
    symtab_get_def("DSolve`AutonomousSystem")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`AutonomousSystem",
        "DSolve`AutonomousSystem[eqns, {x,y}, t] solves a 2-D autonomous system "
        "x'==f(x,y), y'==g(x,y) (including nonlinear) by the phase-plane reduction: "
        "solve the orbit dy/dx==g/f, then reconstruct x(t) from x'==f along it.");
}
