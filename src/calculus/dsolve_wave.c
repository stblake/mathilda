/*
 * dsolve_wave.c — DSolve`WaveDAlembert: the initial-value problem for the
 * one-dimensional wave equation on the whole line, by d'Alembert's formula.
 *
 * Solves    u_{tt} == c^2 u_{xx}    with    u(x, t0) == f(x),  u_t(x, t0) == g(x)
 *
 *     u(x, t) = 1/2 ( f(x - c τ) + f(x + c τ) )
 *               + 1/(2 c) Integrate[g(s), {s, x - c τ, x + c τ}],   τ = t - t0.
 *
 * The dependent function's initial conditions are two-argument point conditions
 * (u[x, t0] == …, Derivative[·,1 in time][u][x, t0] == …) that the shared parser
 * does not recognise (ds_is_condition matches only single-variable funcapps), so
 * they arrive as ordinary equations in P->eq_residuals.  This method therefore
 * receives neq == 3 and sorts them itself: the equation carrying a second-order
 * derivative is the PDE, the two carrying a u-funcapp evaluated at a fixed value
 * of one variable are the initial conditions.  The fixed variable is "time", the
 * free one "space".
 *
 * Verification: because Mathilda cannot differentiate the unevaluated
 * Integrate[g(s), …] of an undefined g, the solution is verified by REBUILDING it
 * with concrete test data f == Cos, g == Sin — a decidable check of the
 * construction (which is a fixed linear formula, correct for any f, g exactly
 * when correct for the test pair) — and requiring the PDE residual and both
 * initial conditions to reduce to zero.
 *
 * First cut: whole-line 1-D wave with constant wave speed, principal part only
 * (no lower-order terms), homogeneous.  Inhomogeneous forcing, the half-line /
 * boundary problems, and Piecewise data are future work.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

static Expr* mul(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Times, a, b)); }
static Expr* add(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Plus, a, b)); }
static Expr* sub(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Subtract, a, b)); }
static Expr* powi(Expr* a, int n)  { return eval_and_free(ds_call2(SYM_Power, a, expr_new_integer(n))); }
static Expr* half(Expr* a) { return mul(eval_and_free(ds_call2(SYM_Power, expr_new_integer(2), expr_new_integer(-1))), a); }

/* Derivative[o1,o2][u][a,b], or u[a,b] when o1==o2==0.  a, b consumed. */
static Expr* u_at(const char* u, int o1, int o2, Expr* a, Expr* b) {
    Expr* head;
    if (o1 == 0 && o2 == 0) head = expr_new_symbol(u);
    else {
        Expr* d = expr_new_function(expr_new_symbol(SYM_Derivative),
                      (Expr*[]){ expr_new_integer(o1), expr_new_integer(o2) }, 2);
        head = expr_new_function(d, (Expr*[]){ expr_new_symbol(u) }, 1);
    }
    return expr_new_function(head, (Expr*[]){ a, b }, 2);
}

/* If node is a u-funcapp evaluated at a point (one arg a fixed value, the other
 * the matching independent-variable symbol), return true and fill the order
 * (oi,oj), whether v1 is the free variable, and t0 (borrowed) = the fixed value. */
static bool match_u_point(const Expr* node, const char* u, const char* v1, const char* v2,
                          int* oi, int* oj, bool* free_is_v1, const Expr** t0) {
    if (!node || node->type != EXPR_FUNCTION || node->data.function.arg_count != 2) return false;
    const Expr* h = node->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        if (h->data.symbol.name != u) return false;
        *oi = 0; *oj = 0;
    } else if (h->type == EXPR_FUNCTION && h->data.function.arg_count == 1
               && h->data.function.args[0]->type == EXPR_SYMBOL
               && h->data.function.args[0]->data.symbol.name == u) {
        const Expr* d = h->data.function.head;          /* Derivative[i,j] */
        if (d->type != EXPR_FUNCTION || d->data.function.arg_count != 2
            || d->data.function.head->type != EXPR_SYMBOL
            || d->data.function.head->data.symbol.name != SYM_Derivative
            || d->data.function.args[0]->type != EXPR_INTEGER
            || d->data.function.args[1]->type != EXPR_INTEGER) return false;
        *oi = (int)d->data.function.args[0]->data.integer;
        *oj = (int)d->data.function.args[1]->data.integer;
    } else return false;

    const Expr* a = node->data.function.args[0];
    const Expr* b = node->data.function.args[1];
    bool a_is_v1 = a->type == EXPR_SYMBOL && a->data.symbol.name == v1;
    bool b_is_v2 = b->type == EXPR_SYMBOL && b->data.symbol.name == v2;
    if (a_is_v1 && b_is_v2) return false;               /* both free: a PDE term */
    if (a_is_v1 && ds_free_of(b, v1) && ds_free_of(b, v2)) { *free_is_v1 = true;  *t0 = b; return true; }
    if (b_is_v2 && ds_free_of(a, v1) && ds_free_of(a, v2)) { *free_is_v1 = false; *t0 = a; return true; }
    return false;
}

/* Depth-first search for the first u-at-a-point node in `e` (borrowed). */
static const Expr* find_u_point(const Expr* e, const char* u, const char* v1, const char* v2,
                                int* oi, int* oj, bool* free_is_v1, const Expr** t0) {
    if (!e) return NULL;
    if (match_u_point(e, u, v1, v2, oi, oj, free_is_v1, t0)) return e;
    if (e->type == EXPR_FUNCTION) {
        const Expr* r = find_u_point(e->data.function.head, u, v1, v2, oi, oj, free_is_v1, t0);
        if (r) return r;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            r = find_u_point(e->data.function.args[i], u, v1, v2, oi, oj, free_is_v1, t0);
            if (r) return r;
        }
    }
    return NULL;
}

/* Build the d'Alembert body for the given data (all borrowed; result owned). */
static Expr* dalembert_body(const Expr* f_expr, const Expr* g_expr, const char* space,
                            const char* time, const Expr* t0, const Expr* c) {
    Expr* tau = sub(expr_new_symbol(time), expr_copy((Expr*)t0));
    Expr* ctau = mul(expr_copy((Expr*)c), tau);
    Expr* lo = sub(expr_new_symbol(space), expr_copy(ctau));
    Expr* hi = add(expr_new_symbol(space), ctau);        /* consumes ctau */

    Expr* f_lo = ds_subst(expr_copy((Expr*)f_expr), expr_new_symbol(space), expr_copy(lo));
    Expr* f_hi = ds_subst(expr_copy((Expr*)f_expr), expr_new_symbol(space), expr_copy(hi));
    Expr* disp = half(add(f_lo, f_hi));

    /* Integration dummy: a clean `K` when it cannot capture (absent from the data
     * and distinct from the variables), else a collision-proof internal symbol. */
    const char* K = intern_symbol("K");
    const char* s = (K != space && K != time && !ds_contains(g_expr, K) && !ds_contains(c, K))
                    ? K : intern_symbol("DSolve`waveS");
    Expr* integrand = ds_subst(expr_copy((Expr*)g_expr), expr_new_symbol(space), expr_new_symbol(s));
    Expr* iter = expr_new_function(expr_new_symbol(SYM_List),
                     (Expr*[]){ expr_new_symbol(s), lo, hi }, 3);   /* consumes lo, hi */
    Expr* integral = eval_and_free(expr_new_function(expr_new_symbol(SYM_Integrate),
                         (Expr*[]){ integrand, iter }, 2));
    Expr* vel = mul(powi(mul(expr_new_integer(2), expr_copy((Expr*)c)), -1), integral);

    return add(disp, vel);
}

/* Verify the construction on concrete test data f == Cos, g == Sin. */
static bool dalembert_verify(const char* space, const char* time, const Expr* t0,
                             const Expr* c, const Expr* c_squared) {
    Expr* fcos = ds_call1("Cos", expr_new_symbol(space));
    Expr* gsin = ds_call1("Sin", expr_new_symbol(space));
    Expr* b = dalembert_body(fcos, gsin, space, time, t0, c);
    expr_free(fcos); expr_free(gsin);

    /* PDE: u_tt - c^2 u_xx == 0 */
    Expr* utt = ds_d(ds_d(expr_copy(b), expr_new_symbol(time)), expr_new_symbol(time));
    Expr* uxx = ds_d(ds_d(expr_copy(b), expr_new_symbol(space)), expr_new_symbol(space));
    Expr* pde = sub(utt, mul(expr_copy((Expr*)c_squared), uxx));
    bool ok = ds_is_zero(pde);
    expr_free(pde);

    /* IC1: b|time=t0 == Cos[space] */
    if (ok) {
        Expr* at0 = ds_subst(expr_copy(b), expr_new_symbol(time), expr_copy((Expr*)t0));
        Expr* d1 = sub(at0, ds_call1("Cos", expr_new_symbol(space)));
        ok = ds_is_zero(d1);
        expr_free(d1);
    }
    /* IC2: D[b,time]|time=t0 == Sin[space] */
    if (ok) {
        Expr* bt = ds_d(expr_copy(b), expr_new_symbol(time));
        bt = ds_subst(bt, expr_new_symbol(time), expr_copy((Expr*)t0));
        Expr* d2 = sub(bt, ds_call1("Sin", expr_new_symbol(space)));
        ok = ds_is_zero(d2);
        expr_free(d2);
    }
    expr_free(b);
    return ok;
}

/* Extract the coefficient of the term `t` (a u-funcapp literal, consumed) in the
 * residual R (borrowed): substitute t -> s, return dR/ds. */
static Expr* coeff_of(const Expr* R, Expr* t) {
    const char* s = intern_symbol("DSolve`waveCoef");
    Expr* Rs = ds_subst(expr_copy((Expr*)R), t, expr_new_symbol(s));
    return ds_d(Rs, expr_new_symbol(s));
}

/* Run the wave IVP end to end; returns {{u -> ...}} or NULL (decline). */
Expr* dsolve_wave_ivp_run(DSolveProblem* P) {
    if (!P->is_pde || P->nfun != 1 || P->nind != 2 || P->neq != 3) return NULL;
    const char* u  = P->fun_names[0];
    const char* v1 = P->ind_names[0];
    const char* v2 = P->ind_names[1];

    /* sort the three equations into two ICs and one PDE */
    int    ic_oi[2], ic_oj[2]; bool ic_free_v1[2]; Expr* ic_t0[2]; Expr* ic_rhs[2];
    int nic = 0, pde_idx = -1;
    for (size_t e = 0; e < 3; e++) {
        int oi, oj; bool fv1; const Expr* t0;
        const Expr* T = find_u_point(P->eq_residuals[e], u, v1, v2, &oi, &oj, &fv1, &t0);
        if (!T) { if (pde_idx < 0) pde_idx = (int)e; else pde_idx = -2; continue; }
        if (nic >= 2) { nic = 99; break; }
        /* rhs = -(R|T->0) / coeff, sign-robust for u on either side */
        Expr* cf   = coeff_of(P->eq_residuals[e], expr_copy((Expr*)T));
        Expr* rest = ds_subst(expr_copy(P->eq_residuals[e]), expr_copy((Expr*)T), expr_new_integer(0));
        Expr* rhs  = ds_simplify(mul(mul(expr_new_integer(-1), rest), powi(cf, -1)));
        ic_oi[nic] = oi; ic_oj[nic] = oj; ic_free_v1[nic] = fv1;
        ic_t0[nic] = expr_copy((Expr*)t0); ic_rhs[nic] = rhs; nic++;
    }
    if (nic != 2 || pde_idx < 0) {
        for (int i = 0; i < nic && i < 2; i++) { expr_free(ic_t0[i]); expr_free(ic_rhs[i]); }
        return NULL;
    }

    /* both ICs must fix the same variable at the same value */
    Expr* dt0 = sub(expr_copy(ic_t0[0]), expr_copy(ic_t0[1]));
    bool same_t0 = ds_is_zero(dt0);
    expr_free(dt0);
    if (ic_free_v1[0] != ic_free_v1[1] || !same_t0) {
        for (int i = 0; i < 2; i++) { expr_free(ic_t0[i]); expr_free(ic_rhs[i]); }
        return NULL;
    }
    bool space_is_v1 = ic_free_v1[0];
    const char* space = space_is_v1 ? v1 : v2;
    const char* time  = space_is_v1 ? v2 : v1;
    /* order of a time-derivative, in the (v1,v2) multi-index, is the time slot */
    int disp = -1, vel = -1;
    for (int i = 0; i < 2; i++) {
        if (ic_oi[i] == 0 && ic_oj[i] == 0) disp = i;
        else {
            int time_ord = space_is_v1 ? ic_oj[i] : ic_oi[i];
            int space_ord = space_is_v1 ? ic_oi[i] : ic_oj[i];
            if (time_ord == 1 && space_ord == 0) vel = i;
        }
    }
    Expr* f_expr = NULL; Expr* g_expr = NULL; Expr* t0 = NULL;
    bool ok = (disp >= 0 && vel >= 0 && disp != vel);
    if (ok) { f_expr = ic_rhs[disp]; g_expr = ic_rhs[vel]; t0 = ic_t0[disp]; }

    /* PDE: A_tt u_{tt} + A_ss u_{ss} == 0 (no mixed, no lower order, homogeneous) */
    Expr* c = NULL; Expr* c_squared = NULL;
    if (ok) {
        const Expr* R = P->eq_residuals[pde_idx];
        int tt1 = space_is_v1 ? 0 : 2, tt2 = space_is_v1 ? 2 : 0;   /* Derivative[·] for u_tt */
        int ss1 = space_is_v1 ? 2 : 0, ss2 = space_is_v1 ? 0 : 2;   /* … for u_ss */
        Expr* Att = coeff_of(R, u_at(u, tt1, tt2, expr_new_symbol(v1), expr_new_symbol(v2)));
        Expr* Ass = coeff_of(R, u_at(u, ss1, ss2, expr_new_symbol(v1), expr_new_symbol(v2)));
        Expr* Axy = coeff_of(R, u_at(u, 1, 1, expr_new_symbol(v1), expr_new_symbol(v2)));
        Expr* Ax  = coeff_of(R, u_at(u, 1, 0, expr_new_symbol(v1), expr_new_symbol(v2)));
        Expr* Ay  = coeff_of(R, u_at(u, 0, 1, expr_new_symbol(v1), expr_new_symbol(v2)));
        Expr* A0  = coeff_of(R, u_at(u, 0, 0, expr_new_symbol(v1), expr_new_symbol(v2)));
        ok = ds_free_of(Att, v1) && ds_free_of(Att, v2)
          && ds_free_of(Ass, v1) && ds_free_of(Ass, v2)
          && !ds_is_zero(Att) && ds_is_zero(Axy)
          && ds_is_zero(Ax) && ds_is_zero(Ay) && ds_is_zero(A0);
        if (ok) {
            c_squared = ds_simplify(mul(expr_new_integer(-1), mul(expr_copy(Ass), powi(expr_copy(Att), -1))));
            /* d'Alembert needs a hyperbolic (positive) speed; reject Δ<=0 (elliptic/parabolic) */
            Expr* sg = eval_and_free(ds_call1("Sign", expr_copy(c_squared)));
            if (sg->type == EXPR_INTEGER && sg->data.integer <= 0) ok = false;
            expr_free(sg);
            if (ds_is_zero(c_squared)) ok = false;
            if (ok) c = ds_call1("Sqrt", expr_copy(c_squared));
        }
        expr_free(Att); expr_free(Ass); expr_free(Axy); expr_free(Ax); expr_free(Ay); expr_free(A0);
    }

    Expr* result = NULL;
    if (ok && dalembert_verify(space, time, t0, c, c_squared)) {
        Expr* body = dalembert_body(f_expr, g_expr, space, time, t0, c);
        /* assemble {{u -> Function[{v1,v2}, body]}} (or applied u[v1,v2] -> body) */
        Expr* lhs; Expr* rhs;
        if (P->applied) {
            lhs = expr_new_function(expr_new_symbol(u),
                      (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2);
            rhs = body;
        } else {
            lhs = expr_new_symbol(u);
            Expr* vars = expr_new_function(expr_new_symbol(SYM_List),
                             (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2);
            rhs = expr_new_function(expr_new_symbol(SYM_Function), (Expr*[]){ vars, body }, 2);
        }
        Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule), (Expr*[]){ lhs, rhs }, 2);
        Expr* inner = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ rule }, 1);
        result = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ inner }, 1);
    }

    if (c) expr_free(c);
    if (c_squared) expr_free(c_squared);
    for (int i = 0; i < 2; i++) { expr_free(ic_t0[i]); expr_free(ic_rhs[i]); }
    return result;
}

static Expr* builtin_dsolve_wave(Expr* res) {
    DSolveProblem P;
    if (!dsolve_parse(res, &P)) { dsolve_problem_free(&P); return NULL; }
    Expr* r = dsolve_wave_ivp_run(&P);
    dsolve_problem_free(&P);
    return r;
}

void dsolve_wave_init(void) {
    symtab_add_builtin("DSolve`WaveDAlembert", builtin_dsolve_wave);
    symtab_get_def("DSolve`WaveDAlembert")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`WaveDAlembert",
        "DSolve`WaveDAlembert[{u_tt == c^2 u_xx, u[x,t0] == f[x], "
        "Derivative[0,1][u][x,t0] == g[x]}, u, {x,t}] solves the initial-value "
        "problem for the one-dimensional wave equation on the whole line by "
        "d'Alembert's formula u(x,t) = (f(x - c τ) + f(x + c τ))/2 + "
        "1/(2 c) Integrate[g(s), {s, x - c τ, x + c τ}], τ = t - t0.  The two "
        "conditions fix one variable (the time) and leave the other free (the "
        "space); the initial velocity integral is kept unevaluated for an undefined "
        "g.  Solved automatically by DSolve as well.  First cut: constant wave "
        "speed, principal part only, homogeneous — inhomogeneous / half-line / "
        "Piecewise data are future.");
}
