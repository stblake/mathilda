/*
 * dsolve_pdesep.c — DSolve`SeparationOfVariables: separated product solutions of
 * a linear PDE by the classical separation ansatz  u(v1, v2) == X(v1) Y(v2).
 *
 * Target class (first cut): a homogeneous linear PDE with constant coefficients
 * and NO mixed derivative term (no `Derivative[i,j][u]` with i>0 AND j>0) —
 *
 *     Σ_{i>=1} a_i u_{v1..v1 (i)} + Σ_{j>=1} b_j u_{v2..v2 (j)} + e u == 0
 *
 * Substituting u = X(v1) Y(v2) and dividing by X Y gives
 *
 *     [ Σ_i a_i X^(i)/X ]  +  [ Σ_j b_j Y^(j)/Y + e ]  ==  0,
 *
 * whose first bracket depends only on v1 and second only on v2, so each is the
 * separation constant λ / −λ.  This yields two constant-coefficient linear ODEs
 * in the separation parameter λ, solved by recursing into the scalar DSolve
 * cascade:
 *
 *     x-ODE:  Σ_i a_i X^(i) − λ X == 0
 *     y-ODE:  Σ_j b_j Y^(j) + (e + λ) Y == 0
 *
 * The product X(v1) Y(v2) is a solution for every λ; the general solution is a
 * superposition over λ (an eigenfunction expansion under boundary conditions),
 * so this returns the *representative product mode* with λ as a generated
 * constant — hence a **pinned-only** method (opt-in, never auto-dispatched;
 * matching DSolve`FirstOrderPowerSeries / DSolve`EigenvalueProblem).  Every
 * returned body back-substitutes to zero.
 *
 * Example: the heat equation  u_t == u_xx  ({x, t}) →
 *   X'' − λ X == 0,  −Y' + λ Y == 0  →  u = E^(λ t)(C[1] E^(−√λ x) + C[2] E^(√λ x)).
 *
 * First cut: constant coefficients, no mixed term, homogeneous.  Variable
 * (product-separable) coefficients and BC-driven eigenfunction expansions are
 * future work.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>
#include <stdio.h>

static Expr* mul(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Times, a, b)); }
static Expr* add(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Plus, a, b)); }

/* Interned fresh symbol "DSolve`sepS<i>_<j>" for the derivative-term placeholder. */
static const char* sepsym(int i, int j) {
    char buf[40];
    snprintf(buf, sizeof buf, "DSolve`sepS%d_%d", i, j);
    return intern_symbol(buf);
}

/* Derivative[o1,o2][u][v1,v2] */
static Expr* pdesep_deriv(const char* u, int o1, int o2, const char* v1, const char* v2) {
    Expr* d = expr_new_function(expr_new_symbol(SYM_Derivative),
                  (Expr*[]){ expr_new_integer(o1), expr_new_integer(o2) }, 2);
    Expr* du = expr_new_function(d, (Expr*[]){ expr_new_symbol(u) }, 1);
    return expr_new_function(du, (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2);
}

/* If `e` is Derivative[i,j][u][v1,v2], return i+j; else -1. */
static int pdesep_term_order(const Expr* e, const char* u) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return -1;
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_FUNCTION || h->data.function.arg_count != 1
        || h->data.function.args[0]->type != EXPR_SYMBOL
        || h->data.function.args[0]->data.symbol.name != u) return -1;
    const Expr* d = h->data.function.head;
    if (d->type != EXPR_FUNCTION || d->data.function.arg_count != 2
        || d->data.function.head->type != EXPR_SYMBOL
        || d->data.function.head->data.symbol.name != SYM_Derivative
        || d->data.function.args[0]->type != EXPR_INTEGER
        || d->data.function.args[1]->type != EXPR_INTEGER) return -1;
    return (int)(d->data.function.args[0]->data.integer
               + d->data.function.args[1]->data.integer);
}

static int pdesep_scan_order(const Expr* e, const char* u) {
    int best = pdesep_term_order(e, u);
    if (best < 0) best = 0;
    if (e && e->type == EXPR_FUNCTION) {
        int ho = pdesep_scan_order(e->data.function.head, u);
        if (ho > best) best = ho;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            int a = pdesep_scan_order(e->data.function.args[i], u);
            if (a > best) best = a;
        }
    }
    return best;
}

/* Solve the ODE  lhs == 0  for the applied function fname[var] by recursing into
 * the scalar DSolve cascade; returns the (first-branch) body, or NULL. */
static Expr* pdesep_solve_ode(Expr* lhs, const char* fname, const char* var) {
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                   (Expr*[]){ lhs, expr_new_integer(0) }, 2);
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ eq, ds_make_funcapp(fname, 0, var), expr_new_symbol(var) }, 3);
    Expr* r = eval_and_free(call);
    size_t nb = 0;
    Expr** bodies = dsolve_extract_applied_bodies(r, fname, &nb);
    expr_free(r);
    if (!bodies || nb == 0) { free(bodies); return NULL; }
    Expr* body = bodies[0];
    for (size_t k = 1; k < nb; k++) expr_free(bodies[k]);
    free(bodies);
    /* require a closed form (no unevaluated DSolve / Integrate) */
    if (ds_has_head(body, "DSolve") || ds_has_head(body, SYM_Integrate)) {
        expr_free(body); return NULL;
    }
    return body;
}

Expr** dsolve_pdesep_solve(DSolveProblem* P) {
    if (P->nfun != 1 || P->nind != 2 || P->neq != 1) return NULL;
    const char* u  = P->fun_names[0];
    const char* v1 = P->ind_names[0];
    const char* v2 = P->ind_names[1];

    int maxord = pdesep_scan_order(P->eq_residuals[0], u);
    if (maxord < 1 || maxord > 6) return NULL;

    /* collect placeholder symbols and substitute every derivative term + u */
    const char* syms[64]; int nsym = 0;
    for (int i = 0; i <= maxord; i++)
        for (int j = 0; i + j <= maxord; j++) syms[nsym++] = sepsym(i, j);

    Expr* R = expr_copy(P->eq_residuals[0]);
    for (int s = maxord; s >= 1; s--)
        for (int i = s; i >= 0; i--)
            R = ds_subst(R, pdesep_deriv(u, i, s - i, v1, v2), expr_new_symbol(sepsym(i, s - i)));
    R = ds_subst(R, expr_new_function(expr_new_symbol(u),
                     (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2),
                 expr_new_symbol(sepsym(0, 0)));

    /* validate: linear, constant coefficients, no mixed term, homogeneous */
    bool ok = true;
    bool has_x = false, has_y = false;
    for (int i = 0; i <= maxord && ok; i++)
        for (int j = 0; i + j <= maxord && ok; j++) {
            Expr* c = ds_d(expr_copy(R), expr_new_symbol(sepsym(i, j)));
            for (int k = 0; k < nsym && ok; k++)
                if (!ds_free_of(c, syms[k])) ok = false;                 /* linear */
            if (ok && (!ds_free_of(c, v1) || !ds_free_of(c, v2))) ok = false; /* constant */
            bool zero = ds_is_zero(c);
            if (ok && i > 0 && j > 0 && !zero) ok = false;               /* no mixed */
            if (i >= 1 && j == 0 && !zero) has_x = true;
            if (i == 0 && j >= 1 && !zero) has_y = true;
            expr_free(c);
        }
    if (ok) {
        Expr* R0 = expr_copy(R);
        for (int k = 0; k < nsym; k++)
            R0 = ds_subst(R0, expr_new_symbol(syms[k]), expr_new_integer(0));
        if (!ds_is_zero(R0)) ok = false;                                 /* homogeneous */
        expr_free(R0);
    }
    if (!ok || !has_x || !has_y) { expr_free(R); return NULL; }

    const char* lam  = intern_symbol("DSolve`sepLam");
    const char* sepX = intern_symbol("DSolve`sepX");
    const char* sepY = intern_symbol("DSolve`sepY");

    /* x-ODE:  Σ_{i>=1} a_i X^(i) − λ X == 0 */
    Expr* xlhs = NULL; int p = 0;
    for (int i = 1; i <= maxord; i++) {
        Expr* ai = ds_d(expr_copy(R), expr_new_symbol(sepsym(i, 0)));
        if (ds_is_zero(ai)) { expr_free(ai); continue; }
        Expr* term = mul(ai, ds_make_funcapp(sepX, i, v1));
        xlhs = xlhs ? add(xlhs, term) : term;
        if (i > p) p = i;
    }
    xlhs = add(xlhs, mul(expr_new_integer(-1),
                         mul(expr_new_symbol(lam), ds_make_funcapp(sepX, 0, v1))));

    /* y-ODE:  Σ_{j>=1} b_j Y^(j) + (e + λ) Y == 0 */
    Expr* ylhs = NULL; int q = 0;
    for (int j = 1; j <= maxord; j++) {
        Expr* bj = ds_d(expr_copy(R), expr_new_symbol(sepsym(0, j)));
        if (ds_is_zero(bj)) { expr_free(bj); continue; }
        Expr* term = mul(bj, ds_make_funcapp(sepY, j, v2));
        ylhs = ylhs ? add(ylhs, term) : term;
        if (j > q) q = j;
    }
    Expr* ecoef = ds_d(expr_copy(R), expr_new_symbol(sepsym(0, 0)));     /* e */
    ylhs = add(ylhs, mul(add(ecoef, expr_new_symbol(lam)), ds_make_funcapp(sepY, 0, v2)));
    expr_free(R);

    Expr* Xbody = pdesep_solve_ode(xlhs, sepX, v1);
    if (!Xbody) { expr_free(ylhs); return NULL; }
    Expr* Ybody = pdesep_solve_ode(ylhs, sepY, v2);
    if (!Ybody) { expr_free(Xbody); return NULL; }

    /* Absorb the single redundant overall scale (the product of two homogeneous
     * families is over-parameterized by one).  Fixing a FIRST-order side's lone
     * constant to 1 loses no generality (the other side already free-scales);
     * renumber the rest to a contiguous C[1..] block and place λ last. */
    int offset = 0, lam_index;
    if (q == 1) {
        Ybody = ds_subst(Ybody, ds_const(1), expr_new_integer(1));
        Xbody = dsolve_renumber_constants(Xbody, p, &offset);
        lam_index = offset + 1;
    } else if (p == 1) {
        Xbody = ds_subst(Xbody, ds_const(1), expr_new_integer(1));
        Ybody = dsolve_renumber_constants(Ybody, q, &offset);
        lam_index = offset + 1;
    } else {
        Xbody = dsolve_renumber_constants(Xbody, p, &offset);
        Ybody = dsolve_renumber_constants(Ybody, q, &offset);
        lam_index = offset + 1;
    }

    Expr* body = mul(Xbody, Ybody);
    body = ds_subst(body, expr_new_symbol(lam), ds_const(lam_index));

    Expr** bodies = malloc(sizeof(Expr*));
    bodies[0] = body;
    return bodies;
}

static Expr* builtin_dsolve_pdesep(Expr* res) {
    return dsolve_method_builtin_pde(res, dsolve_pdesep_solve);
}

void dsolve_pdesep_init(void) {
    symtab_add_builtin("DSolve`SeparationOfVariables", builtin_dsolve_pdesep);
    symtab_get_def("DSolve`SeparationOfVariables")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`SeparationOfVariables",
        "DSolve`SeparationOfVariables[eqn, u, {v1, v2}] finds a separated product "
        "solution u == X(v1) Y(v2) of a homogeneous, constant-coefficient linear "
        "PDE with no mixed derivative term. Dividing by X Y separates it into two "
        "constant-coefficient ODEs in a separation constant λ — Σ a_i X^(i) − λ X == 0 "
        "and Σ b_j Y^(j) + (e + λ) Y == 0 — each solved by recursing into the scalar "
        "cascade; λ becomes a generated constant. Example: the heat equation "
        "u_t == u_xx gives E^(λ t)(C[1] E^(−√λ x) + C[2] E^(√λ x)). Returns the "
        "representative product mode (the general solution is a superposition over λ), "
        "so it is pinned-only — not in the automatic cascade. Declines mixed-derivative, "
        "inhomogeneous, and non-constant-coefficient equations.");
}
