/*
 * dsolve_heat.c — DSolve`HeatKernel: the Cauchy (initial-value) problem for the
 * one-dimensional heat equation on the whole line, by the heat kernel.
 *
 * Solves    u_t == k u_xx    with    u(x, t0) == f(x)    (whole line, t > t0):
 *
 *     u(x, t) = 1/Sqrt[4 π k τ] Integrate[f(K) Exp[-(x - K)^2/(4 k τ)], {K, -∞, ∞}],
 *               τ = t - t0.
 *
 * The single initial condition is a two-argument point condition the shared
 * parser does not split off (ds_is_condition matches only single-variable
 * funcapps), so the problem arrives as two equations (neq == 2); this method
 * sorts the PDE from the condition itself (the fixed variable of the condition is
 * "time", the free one "space").
 *
 * The Gaussian convolution has no elementary antiderivative in Mathilda (and the
 * improper integral does not evaluate for a symbolic f), so the solution is
 * returned as the unevaluated heat-kernel integral (matching Mathematica for a
 * general f), and the integral node is built WITHOUT evaluation (the Function
 * body holds it).  Verification is at the KERNEL level: the heat kernel
 * G = 1/Sqrt[4 π k τ] Exp[-(x - K)^2/(4 k τ)] is a decidable Gaussian, and the
 * convolution solves the PDE because G does (differentiation under the integral,
 * constant infinite limits).  We therefore require D[G,t] - k D[G,{x,2}] == 0
 * (which also checks that k was read correctly); the initial condition holds by
 * the heat kernel's convergence to δ(x - K) as τ → 0+ (a distributional limit,
 * not a back-substitution).
 *
 * First cut: pure heat (no advection u_x, no reaction u term), constant positive
 * diffusivity, homogeneous.  The Erf-producing step/box data (where the integral
 * evaluates), advection–diffusion, reaction, and finite-interval Fourier-series
 * problems are future work.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

static Expr* mul(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Times, a, b)); }
static Expr* sub(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Subtract, a, b)); }
static Expr* powi(Expr* a, int n)  { return eval_and_free(ds_call2(SYM_Power, a, expr_new_integer(n))); }

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

/* If node is a u-funcapp evaluated at a point (one arg fixed, the other the
 * matching independent-variable symbol), fill order/free-var/t0 and return true. */
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
        const Expr* d = h->data.function.head;
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
    if (a_is_v1 && b_is_v2) return false;
    if (a_is_v1 && ds_free_of(b, v1) && ds_free_of(b, v2)) { *free_is_v1 = true;  *t0 = b; return true; }
    if (b_is_v2 && ds_free_of(a, v1) && ds_free_of(a, v2)) { *free_is_v1 = false; *t0 = a; return true; }
    return false;
}

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

/* Coefficient of the term `t` (a u-funcapp literal, consumed) in residual R. */
static Expr* coeff_of(const Expr* R, Expr* t) {
    const char* s = intern_symbol("DSolve`heatCoef");
    Expr* Rs = ds_subst(expr_copy((Expr*)R), t, expr_new_symbol(s));
    return ds_d(Rs, expr_new_symbol(s));
}

/* 1/Sqrt[4 π k τ] (owned; k, tau borrowed). */
static Expr* heat_prefactor(const Expr* k, const Expr* tau) {
    Expr* arg = mul(mul(mul(expr_new_integer(4), expr_new_symbol(SYM_Pi)),
                        expr_copy((Expr*)k)), expr_copy((Expr*)tau));
    return powi(eval_and_free(ds_call1("Sqrt", arg)), -1);
}

/* Exp[-(space - pt)^2/(4 k τ)] (owned; pt consumed; k, tau borrowed). */
static Expr* heat_gaussian(const char* space, Expr* pt, const Expr* k, const Expr* tau) {
    Expr* diff2 = powi(sub(expr_new_symbol(space), pt), 2);
    Expr* denom = powi(mul(mul(expr_new_integer(4), expr_copy((Expr*)k)), expr_copy((Expr*)tau)), -1);
    Expr* arg = mul(mul(expr_new_integer(-1), diff2), denom);
    return eval_and_free(ds_call1("Exp", arg));
}

Expr* dsolve_heat_run(DSolveProblem* P) {
    if (!P->is_pde || P->nfun != 1 || P->nind != 2 || P->neq != 2) return NULL;
    const char* u  = P->fun_names[0];
    const char* v1 = P->ind_names[0];
    const char* v2 = P->ind_names[1];

    /* sort: one initial condition u[space,t0]==f (order 0) + one PDE */
    int ic_idx = -1, pde_idx = -1, oi = 0, oj = 0; bool free_v1 = false; const Expr* t0b = NULL;
    for (size_t e = 0; e < 2; e++) {
        int i, j; bool fv1; const Expr* tt;
        const Expr* T = find_u_point(P->eq_residuals[e], u, v1, v2, &i, &j, &fv1, &tt);
        if (T && i == 0 && j == 0) { ic_idx = (int)e; oi = i; oj = j; free_v1 = fv1; t0b = tt; }
        else if (!T) pde_idx = (int)e;
    }
    (void)oi; (void)oj;
    if (ic_idx < 0 || pde_idx < 0) return NULL;

    const char* space = free_v1 ? v1 : v2;
    const char* time  = free_v1 ? v2 : v1;
    Expr* t0 = expr_copy((Expr*)t0b);

    /* f(space) = the condition's right-hand side (sign-robust) */
    const Expr* Ric = P->eq_residuals[ic_idx];
    Expr* T = u_at(u, 0, 0, free_v1 ? expr_new_symbol(v1) : expr_copy(t0),
                             free_v1 ? expr_copy(t0) : expr_new_symbol(v2));
    Expr* cf   = coeff_of(Ric, expr_copy(T));
    Expr* rest = ds_subst(expr_copy((Expr*)Ric), T, expr_new_integer(0));
    Expr* f_expr = ds_simplify(mul(mul(expr_new_integer(-1), rest), powi(cf, -1)));

    /* PDE: require pure heat  a u_t + c u_xx == 0 (no u_tt, u_x, u_xy, u, forcing) */
    const Expr* R = P->eq_residuals[pde_idx];
    /* space is v1 when free_v1, so time is v2: u_t sits in the v2 slot. */
    int t1 = free_v1 ? 0 : 1, t2 = free_v1 ? 1 : 0;      /* u_t  (order 1 in time) */
    int s1 = free_v1 ? 2 : 0, s2 = free_v1 ? 0 : 2;      /* u_xx (order 2 in space) */
    Expr* Aut = coeff_of(R, u_at(u, t1, t2, expr_new_symbol(v1), expr_new_symbol(v2)));
    Expr* Axx = coeff_of(R, u_at(u, s1, s2, expr_new_symbol(v1), expr_new_symbol(v2)));
    Expr* Att = coeff_of(R, u_at(u, 2*t1, 2*t2, expr_new_symbol(v1), expr_new_symbol(v2)));
    Expr* Axy = coeff_of(R, u_at(u, 1, 1, expr_new_symbol(v1), expr_new_symbol(v2)));
    Expr* Ax  = coeff_of(R, u_at(u, free_v1?1:0, free_v1?0:1, expr_new_symbol(v1), expr_new_symbol(v2)));
    Expr* A0  = coeff_of(R, u_at(u, 0, 0, expr_new_symbol(v1), expr_new_symbol(v2)));
    bool ok = ds_free_of(Aut, v1) && ds_free_of(Aut, v2)
           && ds_free_of(Axx, v1) && ds_free_of(Axx, v2)
           && !ds_is_zero(Aut) && !ds_is_zero(Axx)
           && ds_is_zero(Att) && ds_is_zero(Axy) && ds_is_zero(Ax) && ds_is_zero(A0);

    Expr* k = NULL;
    if (ok) {
        k = ds_simplify(mul(expr_new_integer(-1), mul(expr_copy(Axx), powi(expr_copy(Aut), -1))));
        Expr* sg = eval_and_free(ds_call1("Sign", expr_copy(k)));
        if (sg->type == EXPR_INTEGER && sg->data.integer <= 0) ok = false;  /* need k > 0 */
        expr_free(sg);
        if (ds_is_zero(k)) ok = false;
    }
    expr_free(Aut); expr_free(Axx); expr_free(Att); expr_free(Axy); expr_free(Ax); expr_free(A0);

    Expr* result = NULL;
    if (ok) {
        Expr* tau = sub(expr_new_symbol(time), expr_copy(t0));

        /* verify the kernel: D[G,time] - k D[G,{space,2}] == 0 */
        const char* yp = intern_symbol("DSolve`heatY");
        Expr* G = mul(heat_prefactor(k, tau), heat_gaussian(space, expr_new_symbol(yp), k, tau));
        Expr* Gt  = ds_d(expr_copy(G), expr_new_symbol(time));
        Expr* Gxx = ds_d(ds_d(expr_copy(G), expr_new_symbol(space)), expr_new_symbol(space));
        Expr* res = sub(Gt, mul(expr_copy(k), Gxx));
        bool kernel_ok = ds_is_zero(res);
        expr_free(res); expr_free(G);

        if (kernel_ok) {
            /* body = prefactor * Integrate[f(K) Gaussian(space-K), {K, -∞, ∞}]  (unevaluated) */
            const char* K = intern_symbol("K");
            if (K == space || K == time || ds_contains(f_expr, K) || ds_contains(k, K))
                K = intern_symbol("DSolve`heatK");
            Expr* fk = ds_subst(expr_copy(f_expr), expr_new_symbol(space), expr_new_symbol(K));
            Expr* gauss = heat_gaussian(space, expr_new_symbol(K), k, tau);
            Expr* integrand = mul(fk, gauss);
            Expr* iter = expr_new_function(expr_new_symbol(SYM_List),
                             (Expr*[]){ expr_new_symbol(K),
                                        expr_new_function(expr_new_symbol("Times"),
                                            (Expr*[]){ expr_new_integer(-1), expr_new_symbol("Infinity") }, 2),
                                        expr_new_symbol("Infinity") }, 3);
            /* build Integrate[...] RAW so the improper Gaussian is not attempted here */
            Expr* integral = expr_new_function(expr_new_symbol(SYM_Integrate),
                                 (Expr*[]){ integrand, iter }, 2);
            Expr* body = expr_new_function(expr_new_symbol(SYM_Times),
                             (Expr*[]){ heat_prefactor(k, tau), integral }, 2);

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
        expr_free(tau);
    }

    if (k) expr_free(k);
    expr_free(t0); expr_free(f_expr);
    return result;
}

static Expr* builtin_dsolve_heat(Expr* res) {
    DSolveProblem P;
    if (!dsolve_parse(res, &P)) { dsolve_problem_free(&P); return NULL; }
    Expr* r = dsolve_heat_run(&P);
    dsolve_problem_free(&P);
    return r;
}

void dsolve_heat_init(void) {
    symtab_add_builtin("DSolve`HeatKernel", builtin_dsolve_heat);
    symtab_get_def("DSolve`HeatKernel")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`HeatKernel",
        "DSolve`HeatKernel[{u_t == k u_xx, u[x,t0] == f[x]}, u, {x,t}] solves the "
        "Cauchy problem for the one-dimensional heat equation on the whole line by "
        "the heat kernel: u(x,t) = 1/Sqrt[4 π k τ] Integrate[f[K] Exp[-(x-K)^2/"
        "(4 k τ)], {K, -∞, ∞}], τ = t - t0.  The Gaussian convolution is kept "
        "unevaluated (it is nonelementary for a general f).  Verified at the kernel "
        "level (the heat kernel solves the PDE, so the convolution does); the "
        "initial condition holds by the kernel's convergence to a delta as τ -> 0+.  "
        "Solved automatically by DSolve as well.  First cut: pure heat (no advection "
        "or reaction), constant positive diffusivity, homogeneous.");
}
