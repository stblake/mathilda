/*
 * dsolve_ifactor.c — DSolve`ReducibleIntegratingFactor.
 *
 * Solves a nonlinear second-order ODE  y''[x] == Phi(x, y, y')  by finding an
 * INTEGRATING FACTOR mu of a restricted two-variable form and reducing the order
 * with it (Cheb-Terrab & Roche, "Integrating Factors for Second-order ODEs",
 * J. Symb. Comput. 27 (1999) 501-519; Maple's `_reducible, _mu_*` classes).
 *
 * A function mu(x, y, y') is an integrating factor of  y'' = Phi  when
 * mu (y'' - Phi) is a total x-derivative  dR/dx  of a first integral R(x,y,y').
 * From (2.8)-(2.10):
 *     mu = R_{y'} ,  R = G(x,y) + Integrate[mu, y'] ,
 *     R_x + y' R_y + Phi R_{y'} = 0   (fixes G by quadrature).
 * R(x,y,y') == C[1] is then a FIRST-ORDER ODE, solved by recursing into the
 * scalar cascade; the second constant comes from that solve, giving the full
 * y(x, C[1], C[2]).
 *
 * mu is searched in three restricted forms, cheapest-first:
 *   - mu(x, y)   [_mu_xy]   : Section 2.1 — Phi a degree-<=2 polynomial in y';
 *                             a closed-form mu (Case A) or a linear-ODE mu (Case B).
 *   - mu(x, y')  [_mu_x_y1] : Section 2.2, Lemma 3 — algebraic (Stage 2).
 *   - mu(y, y')  [_mu_y_y1] : Section 2.3 — point-swap y<->x, reuse mu(x,y') (Stage 3).
 *
 * Correctness rests on TWO gates, not on the mu formulas being byte-perfect:
 *   (1) the symbolic first-integral test  A(R) = R_x + y' R_y + Phi R_{y'} == 0
 *       (the paper's exactness condition (2.3)) rejects a wrong mu before any solve;
 *   (2) a numeric back-substitution of the final explicit branch into the original
 *       residual (l2_num_ok-style) rejects a bad inversion.
 * So a wrong mu can only ever produce a clean decline, never a wrong answer.
 *
 * Runs in the 2nd-order-nonlinear cascade tail immediately before
 * SecondOrderSymmetry (the heavier polynomial-ansatz symmetry search): the
 * integrating-factor search is algebraic/cheaper and reaches ODEs with no point
 * symmetry (paper Section 3).
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include "../print.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

/* ---- wall-clock budget + per-toplevel decline memo (mirrors dsolve_lie2.c) ---- */

static time_t g_if_deadline;
static bool if_expired(void) { return time(NULL) >= g_if_deadline; }

/* Optional decline tracing: MATHILDA_IF_DEBUG=1 prints where the search stops. */
static void if_dbg(const char* where, const Expr* e) {
    if (!getenv("MATHILDA_IF_DEBUG")) return;
    fprintf(stderr, "[ifactor] %s", where);
    if (e) { char* s = expr_to_string((Expr*)e); fprintf(stderr, ": %s", s ? s : "?"); free(s); }
    fprintf(stderr, "\n");
}

#define IF_MEMO_SLOTS 32
static uint64_t if_memo_epoch = 0;
static int      if_memo_n = 0;
static uint64_t if_memo[IF_MEMO_SLOTS];
static void if_memo_sync(uint64_t tid) { if (tid != if_memo_epoch) { if_memo_epoch = tid; if_memo_n = 0; } }
static bool if_memo_seen(uint64_t h) {
    for (int i = 0; i < if_memo_n; i++) if (if_memo[i] == h) return true;
    return false;
}
static void if_memo_add(uint64_t h) {
    if (if_memo_n < IF_MEMO_SLOTS && !if_memo_seen(h)) if_memo[if_memo_n++] = h;
}

/* ---- small helpers ---- */

static Expr* if_powi(Expr* base, int n) {
    return eval_and_free(ds_call2(SYM_Power, base, expr_new_integer(n)));
}
static Expr* if_d(Expr* e, const char* v) { return ds_d(e, expr_new_symbol(v)); }
static Expr* if_d2(Expr* e, const char* a, const char* b) { return if_d(if_d(e, a), b); }

static void if_count(const Expr* e, long* acc, long budget) {
    if (*acc > budget || !e) return;
    (*acc)++;
    if (e->type == EXPR_FUNCTION) {
        if_count(e->data.function.head, acc, budget);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if_count(e->data.function.args[i], acc, budget);
    }
}
static bool if_too_big(const Expr* e) {
    long acc = 0; if_count(e, &acc, 60000); return acc > 60000;
}

/* undefined-function / inert-Derivative gate (identical to lie2's) */
static bool if_has_undef_fn(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        if (h->data.symbol.name == SYM_Derivative) return true;
        SymbolDef* d = symtab_lookup(h->data.symbol.name);
        if (d && !d->builtin_func && !d->down_values) return true;
    } else if (if_has_undef_fn(h)) {
        return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (if_has_undef_fn(e->data.function.args[i])) return true;
    return false;
}

/* Evaluate `e` under TimeConstrained (returns $Aborted on timeout). e consumed. */
static Expr* if_beval(Expr* e, int secs) {
    Expr* g = expr_new_function(expr_new_symbol("TimeConstrained"),
                  (Expr*[]){ e, expr_new_integer(secs),
                             expr_new_symbol(intern_symbol("$Aborted")) }, 3);
    return eval_and_free(g);
}
static bool if_aborted(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == intern_symbol("$Aborted");
}
/* Bounded Integrate[e, v]; e consumed, result owned. */
static Expr* if_integrate_b(Expr* e, const char* v) {
    return if_beval(expr_new_function(expr_new_symbol(SYM_Integrate),
                        (Expr*[]){ e, expr_new_symbol(v) }, 2), 2);
}

/* Collapse mixed hyperbolic/trig/exp forms.  A candidate mu recovered from a
 * sub-DSolve may be in Exp form while Phi carries Coth/Sinh etc.; plain
 * evaluation then leaves an identically-zero coefficient uncollapsed
 * (Coth(2x)(E^4x-1) - (E^4x+1) != 0 to the evaluator), which would leak a
 * spurious term into R.  TrigToExp puts everything in Exp form; a (bounded)
 * Simplify then combines E^(ax) E^(-ax) = 1 across the mixed-sign exponentials
 * that Together/Cancel treat as independent atoms.  Runs only on a found mu (not
 * a hot path).  `e` consumed, result owned (falls back to the TrigToExp form on
 * a Simplify timeout). */
static Expr* if_hypsimp(Expr* e) {
    Expr* t = eval_and_free(ds_call1("TrigToExp", e));
    Expr* keep = expr_copy(t);
    Expr* r = if_beval(ds_call1(SYM_Simplify, t), 2);
    if (!r || if_aborted(r)) { expr_free(r); return keep; }
    expr_free(keep);
    return r;
}

/* DSolve[eqn, fname[varname], varname] -> first applied-form RHS body (owned) or
 * NULL, under a TimeConstrained bound.  `eqn` consumed. (mirrors lie2's run_applied) */
static Expr* if_run_applied(Expr* eqn, const char* fname, const char* varname) {
    Expr* lhs  = ds_call1(fname, expr_new_symbol(varname));
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ eqn, lhs, expr_new_symbol(varname) }, 3);
    Expr* guarded = expr_new_function(expr_new_symbol("TimeConstrained"),
                     (Expr*[]){ call, expr_new_integer(3),
                                expr_new_symbol(intern_symbol("$Aborted")) }, 3);
    Expr* r = eval_and_free(guarded);
    Expr* body = NULL;
    if (head_is(r, SYM_List) && r->data.function.arg_count >= 1) {
        Expr* inner = r->data.function.args[0];
        if (head_is(inner, SYM_List)) {
            for (size_t k = 0; k < inner->data.function.arg_count && !body; k++) {
                Expr* rule = inner->data.function.args[k];
                if (head_is(rule, SYM_Rule) && rule->data.function.arg_count == 2) {
                    Expr* rl = rule->data.function.args[0];
                    if (rl->type == EXPR_FUNCTION && rl->data.function.head->type == EXPR_SYMBOL
                        && rl->data.function.head->data.symbol.name == fname)
                        body = expr_copy(rule->data.function.args[1]);
                }
            }
        }
    }
    expr_free(r);
    return body;
}

/* Numeric back-substitution guard (identical policy to lie2's l2_num_ok). */
static double if_abs_at(const Expr* R, const char* xv, double xval, double c1, double c2) {
    Expr* e = expr_copy((Expr*)R);
    e = ds_subst(e, ds_const(1), expr_new_real(c1));
    e = ds_subst(e, ds_const(2), expr_new_real(c2));
    e = ds_subst(e, expr_new_symbol(xv), expr_new_real(xval));
    e = eval_and_free(ds_call1("Abs", e));
    double m = (e && e->type == EXPR_REAL)    ? e->data.real
             : (e && e->type == EXPR_INTEGER) ? (double)e->data.integer
             : NAN;
    expr_free(e);
    return m;
}
/* Replace every free scalar PARAMETER of `R` (a bare symbol other than the
 * independent variable `xv`, the generated-constant head `C`, and the protected
 * constants E/Pi/I/…) by a distinct generic real, so a residual carrying symbolic
 * ODE parameters (a, b, c, k, n, …) still numericizes.  Without this, a correct
 * solution of a symbolic-parameter ODE (the common corpus form) samples to a
 * symbolic — hence non-finite — residual and is wrongly rejected (the M16 lesson).
 * `R` consumed, result owned. */
static Expr* if_instantiate_params(Expr* R, const char* xv) {
    /* Collect ARGUMENT-position symbols only (no Heads->True): a genuine parameter
     * like a/b/k appears as a multiplier, whereas function heads (Cos, AiryAi, …)
     * must NOT be instantiated (that would corrupt the residual). */
    Expr* patt = expr_new_function(expr_new_symbol("Blank"),
                     (Expr*[]){ expr_new_symbol("Symbol") }, 1);
    Expr* lev = expr_new_function(expr_new_symbol(SYM_List),
                     (Expr*[]){ expr_new_integer(0), expr_new_symbol("Infinity") }, 2);
    Expr* cases = expr_new_function(expr_new_symbol("Cases"),
                     (Expr*[]){ expr_copy(R), patt, lev }, 3);
    Expr* syms = eval_and_free(ds_call1("DeleteDuplicates", cases));
    static const char* prot[] = { "E","Pi","I","Infinity","ComplexInfinity",
        "Indeterminate","True","False","Degree","EulerGamma","GoldenRatio","Catalan","C", NULL };
    if (head_is(syms, SYM_List)) {
        int idx = 0;
        for (size_t i = 0; i < syms->data.function.arg_count; i++) {
            Expr* s = syms->data.function.args[i];
            if (!s || s->type != EXPR_SYMBOL) continue;
            const char* nm = s->data.symbol.name;
            if (nm == xv) continue;
            bool skip = false;
            for (int k = 0; prot[k]; k++) if (strcmp(nm, prot[k]) == 0) { skip = true; break; }
            if (skip) continue;
            /* distinct generic real, avoiding small integers / collisions */
            double val = 1.3 + (double)(idx++) * 4.0 / 17.0 + 0.11;
            R = ds_subst(R, expr_new_symbol(nm), expr_new_real(val));
        }
    }
    expr_free(syms);
    return R;
}

static bool if_num_ok(const DSolveProblem* P, const Expr* body,
                      const char* xv, const char* yname) {
    Expr* R = expr_copy(P->eq_residuals[0]);
    Expr* b0 = expr_copy((Expr*)body);
    Expr* b1 = ds_d(expr_copy((Expr*)body), expr_new_symbol(xv));
    Expr* b2 = if_d2(expr_copy((Expr*)body), xv, xv);
    R = ds_subst(R, ds_make_funcapp(yname, 2, xv), b2);
    R = ds_subst(R, ds_make_funcapp(yname, 1, xv), b1);
    R = ds_subst(R, ds_make_funcapp(yname, 0, xv), b0);
    R = if_instantiate_params(R, xv);   /* symbolic ODE parameters -> generic reals */
    const double xs[]  = { 1.7, 2.3, 1.15, 0.6, 3.1 };
    const double c1s[] = { 1.181, 1.4, 2.25, 0.7, 1.9 };
    const double c2s[] = { 0.714, 0.75, 1.6, 1.3, 0.4 };
    int small = 0, big = 0;
    for (int i = 0; i < 5; i++) {
        double m = if_abs_at(R, xv, xs[i], c1s[i], c2s[i]);
        if (isnan(m) || !isfinite(m)) continue;
        if (m < 1e-6) small++;
        else if (m > 1e-3) big++;
    }
    expr_free(R);
    return small >= 2 && big == 0;
}

/* ===================================================================== *
 *  Shared pipeline:  mu -> first integral R -> solve -> verify -> emit  *
 * ===================================================================== */

/* Build the first integral  R = Integrate[mu, p] + G(x,y)  from an integrating
 * factor mu(x,y,p), fixing G by the first-integral condition
 *   R_x + p R_y + Phi R_p == 0   (2.10).
 * With R_p = mu, write S = M0_x + p M0_y + Phi mu (M0 = Integrate[mu,p]); then
 * G_x + p G_y = -S, so (S affine in p)  G_y = -S_p,  G_x = -S|_{p=0}, and G is
 * recovered by the exact-equation quadrature.  Returns R (owned) or NULL.
 * mu, Phi borrowed. */
static Expr* ifactor_build_R(const Expr* Phi, const Expr* mu,
                             const char* xv, const char* Yn, const char* Pn) {
    Expr* M0 = if_integrate_b(expr_copy((Expr*)mu), Pn);   /* Integrate[mu, p] */
    if (!M0 || if_aborted(M0) || ds_has_head(M0, SYM_Integrate) || if_too_big(M0)) {
        expr_free(M0); return NULL;
    }
    /* S = M0_x + p M0_y + Phi mu */
    Expr* M0x = if_d(expr_copy(M0), xv);
    Expr* M0y = if_d(expr_copy(M0), Yn);
    Expr* S = ds_call2(SYM_Plus, M0x,
                ds_call2(SYM_Plus,
                    ds_call2(SYM_Times, expr_new_symbol(Pn), M0y),
                    ds_call2(SYM_Times, expr_copy((Expr*)Phi), expr_copy((Expr*)mu))));
    S = if_hypsimp(eval_and_free(S));    /* collapse mixed hyperbolic/exp identities */
    if_dbg("S", S);
    if (if_too_big(S)) { expr_free(S); expr_free(M0); return NULL; }

    /* G_y = -S_p (must be free of p for a genuine integrating factor of this form) */
    Expr* Sp = if_d(expr_copy(S), Pn);
    if (!ds_free_of(Sp, Pn)) { expr_free(Sp); expr_free(S); expr_free(M0); return NULL; }
    Expr* Gy = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), Sp));
    /* G_x = -(S at p=0) */
    Expr* S0 = ds_subst(expr_copy(S), expr_new_symbol(Pn), expr_new_integer(0));
    expr_free(S);
    Expr* Gx = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), S0));

    /* exact-equation quadrature: G = Integrate[Gx,x] + Integrate[Gy - d/dy G1, y] */
    Expr* G1 = if_integrate_b(expr_copy(Gx), xv);
    expr_free(Gx);
    if (!G1 || if_aborted(G1) || ds_has_head(G1, SYM_Integrate)) {
        expr_free(G1); expr_free(Gy); expr_free(M0); return NULL;
    }
    Expr* corr = eval_and_free(ds_call2(SYM_Subtract, Gy, if_d(expr_copy(G1), Yn)));
    Expr* G2 = if_integrate_b(corr, Yn);
    if (!G2 || if_aborted(G2) || ds_has_head(G2, SYM_Integrate)) {
        expr_free(G1); expr_free(G2); expr_free(M0); return NULL;
    }
    Expr* R = ds_call2(SYM_Plus, M0, ds_call2(SYM_Plus, G1, G2));
    R = eval_and_free(R);
    return R;
}

/* A(R) = R_x + p R_y + Phi R_p == 0  and  R genuinely depends on p (a real
 * first integral of a 2nd-order ODE, i.e. involves y').  R, Phi borrowed. */
static bool ifactor_R_ok(const Expr* R, const Expr* Phi,
                         const char* xv, const char* Yn, const char* Pn) {
    if (!R || if_too_big(R) || !ds_contains(R, Pn)) return false;
    Expr* Rx = if_d(expr_copy((Expr*)R), xv);
    Expr* Ry = if_d(expr_copy((Expr*)R), Yn);
    Expr* Rp = if_d(expr_copy((Expr*)R), Pn);
    Expr* A = ds_call2(SYM_Plus, Rx,
                ds_call2(SYM_Plus,
                    ds_call2(SYM_Times, expr_new_symbol(Pn), Ry),
                    ds_call2(SYM_Times, expr_copy((Expr*)Phi), Rp)));
    A = eval_and_free(A);
    bool z = ds_is_zero(A);
    expr_free(A);
    return z;
}

/* Given a validated first integral R(x,y,p), solve  R(x,y[x],y'[x]) == C[2]  as a
 * first-order ODE (recursing into the cascade; the sub-solve introduces C[1]),
 * numerically verify the explicit branch, and return the body y(x,C[1],C[2]) or
 * NULL.  R borrowed. */
static Expr* ifactor_solve_from_R(const DSolveProblem* P, const Expr* R,
                                  const char* xv, const char* Yn, const char* Pn,
                                  const char* yname) {
    const char* ifK = intern_symbol("DSolve`ifK");   /* first-integral constant */
    Expr* Rf = ds_subst(expr_copy((Expr*)R), expr_new_symbol(Pn), ds_make_funcapp(yname, 1, xv));
    Rf = ds_subst(Rf, expr_new_symbol(Yn), ds_make_funcapp(yname, 0, xv));
    Expr* eqn = expr_new_function(expr_new_symbol(SYM_Equal),
                    (Expr*[]){ Rf, expr_new_symbol(ifK) }, 2);
    if (if_expired()) { expr_free(eqn); return NULL; }
    Expr* body = if_run_applied(eqn, yname, xv);       /* y(x, C[1], ifK) */
    if (!body) return NULL;
    body = ds_subst(body, expr_new_symbol(ifK), ds_const(2));   /* ifK -> C[2] */
    if (ds_free_of(body, xv) || ds_has_head(body, SYM_Solve)
        || ds_has_head(body, SYM_Integrate) || if_too_big(body)
        || !if_num_ok(P, body, xv, yname)) {
        expr_free(body); return NULL;
    }
    return body;
}

/* ===================================================================== *
 *  Stage 1:  mu(x, y)   [_mu_xy]   (Section 2.1)                        *
 * ===================================================================== */

/* Extract a,b,c with Phi == a p^2 + b p + c (each free of p), or return false. */
static bool ifactor_abc(const Expr* Phi, const char* Pn,
                        Expr** a, Expr** b, Expr** c) {
    Expr* pq = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
                   (Expr*[]){ expr_copy((Expr*)Phi), expr_new_symbol(Pn) }, 2));
    bool poly = pq && pq->type == EXPR_SYMBOL && pq->data.symbol.name == SYM_True;
    expr_free(pq);
    if (!poly) return false;
    Expr* cl = eval_and_free(expr_new_function(expr_new_symbol(SYM_CoefficientList),
                   (Expr*[]){ expr_copy((Expr*)Phi), expr_new_symbol(Pn) }, 2));
    if (!head_is(cl, SYM_List) || cl->data.function.arg_count > 3) { expr_free(cl); return false; }
    size_t n = cl->data.function.arg_count;
    *c = (n >= 1) ? expr_copy(cl->data.function.args[0]) : expr_new_integer(0);
    *b = (n >= 2) ? expr_copy(cl->data.function.args[1]) : expr_new_integer(0);
    *a = (n >= 3) ? expr_copy(cl->data.function.args[2]) : expr_new_integer(0);
    expr_free(cl);
    return true;
}

/* mu(x,y) integrating factor via Section 2.1, or NULL.  Phi borrowed. */
static Expr* ifactor_mu_xy(const Expr* Phi, const char* xv, const char* Yn, const char* Pn) {
    Expr *a = NULL, *b = NULL, *c = NULL;
    if (!ifactor_abc(Phi, Pn, &a, &b, &c)) return NULL;

    Expr* ax  = if_d(expr_copy(a), xv);
    Expr* by  = if_d(expr_copy(b), Yn);
    /* disc = 2 a_x - b_y */
    Expr* disc = eval_and_free(ds_call2(SYM_Subtract,
                    ds_call2(SYM_Times, expr_new_integer(2), expr_copy(ax)), expr_copy(by)));
    bool caseB = ds_is_zero(disc);
    expr_free(disc); expr_free(by);

    Expr* Ia = if_integrate_b(expr_copy(a), Yn);          /* Integrate[a, y] */
    if (!Ia || if_aborted(Ia) || ds_has_head(Ia, SYM_Integrate)) {
        expr_free(a); expr_free(c); expr_free(ax); expr_free(Ia); return NULL;
    }
    Expr* I = if_d(expr_copy(Ia), xv);                    /* I = d/dx Integrate[a,y] */
    Expr* mu = NULL;

    if (!caseB) {
        /* Case A (2.15)-(2.17), closed form. */
        Expr* cy = if_d(expr_copy(c), Yn);
        Expr* bx = if_d(expr_copy(b), xv);
        Expr* ac = eval_and_free(ds_call2(SYM_Times, expr_copy(a), expr_copy(c)));
        Expr* phi = eval_and_free(ds_call2(SYM_Subtract,
                        ds_call2(SYM_Subtract, cy, ac), bx));           /* c_y - a c - b_x */
        Expr* axx = if_d(expr_copy(ax), xv);
        Expr* phy = if_d(expr_copy(phi), Yn);
        Expr* Ups = eval_and_free(ds_call2(SYM_Plus, axx,               /* a_xx + a_x b + phi_y */
                        ds_call2(SYM_Plus,
                            ds_call2(SYM_Times, expr_copy(ax), expr_copy(b)), phy)));
        /* existence (2.16): Ups_y - a_x == 0  and  Ups_x + phi + b Ups - Ups^2 == 0 */
        Expr* e1 = eval_and_free(ds_call2(SYM_Subtract, if_d(expr_copy(Ups), Yn), expr_copy(ax)));
        Expr* e2 = ds_call2(SYM_Plus, if_d(expr_copy(Ups), xv),
                     ds_call2(SYM_Plus, expr_copy(phi),
                       ds_call2(SYM_Subtract,
                         ds_call2(SYM_Times, expr_copy(b), expr_copy(Ups)),
                         if_powi(expr_copy(Ups), 2))));
        e2 = eval_and_free(e2);
        bool exists = ds_is_zero(e1) && ds_is_zero(e2);
        expr_free(e1); expr_free(e2);
        if (exists) {
            /* mu = Exp[ Integrate[-Ups + d/dx Integrate[a,y], x] - Integrate[a,y] ] */
            Expr* inner = eval_and_free(ds_call2(SYM_Plus,
                              ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(Ups)),
                              expr_copy(I)));
            Expr* Iinner = if_integrate_b(inner, xv);
            if (Iinner && !if_aborted(Iinner) && !ds_has_head(Iinner, SYM_Integrate)) {
                Expr* arg = eval_and_free(ds_call2(SYM_Subtract, Iinner, expr_copy(Ia)));
                mu = eval_and_free(ds_call1(SYM_Exp, arg));
            } else expr_free(Iinner);
        }
        expr_free(phi); expr_free(Ups);
    } else {
        /* Case B (2.18)-(2.21): mu = nu(x) Exp[-Integrate[a,y]], nu'' = A nu' + B nu. */
        Expr* cy = if_d(expr_copy(c), Yn);
        Expr* ac = eval_and_free(ds_call2(SYM_Times, expr_copy(a), expr_copy(c)));
        Expr* phi = eval_and_free(ds_call2(SYM_Subtract, cy, ac));     /* c_y - a c */
        Expr* axx = if_d(expr_copy(ax), xv);
        Expr* phy = if_d(expr_copy(phi), Yn);
        /* existence (2.18): a_xx - a_x b - phi_y == 0 */
        Expr* ex = eval_and_free(ds_call2(SYM_Subtract,
                       ds_call2(SYM_Subtract, axx, ds_call2(SYM_Times, expr_copy(ax), expr_copy(b))),
                       phy));
        bool exists = ds_is_zero(ex);
        expr_free(ex);
        if (exists) {
            /* A = 2 I - b ,  B = phi + (I - d/dx)(b - I) = phi + I (b-I) - d/dx(b-I) */
            Expr* Acoef = eval_and_free(ds_call2(SYM_Subtract,
                              ds_call2(SYM_Times, expr_new_integer(2), expr_copy(I)), expr_copy(b)));
            Expr* bmI = ds_call2(SYM_Subtract, expr_copy(b), expr_copy(I));
            bmI = eval_and_free(bmI);
            Expr* Bcoef = eval_and_free(ds_call2(SYM_Plus, expr_copy(phi),
                              ds_call2(SYM_Subtract,
                                  ds_call2(SYM_Times, expr_copy(I), expr_copy(bmI)),
                                  if_d(expr_copy(bmI), xv))));
            expr_free(bmI);
            /* nu''(x) - A nu' - B nu == 0, solved via cascade; take one basis solution. */
            if (ds_free_of(Acoef, Yn) && ds_free_of(Bcoef, Yn) && !if_expired()) {
                const char* nun = intern_symbol("DSolve`ifNu");
                Expr* nupp = ds_make_funcapp(nun, 2, xv);
                Expr* nup  = ds_make_funcapp(nun, 1, xv);
                Expr* nu0  = ds_make_funcapp(nun, 0, xv);
                Expr* lhs = ds_call2(SYM_Subtract,
                                ds_call2(SYM_Subtract, nupp, ds_call2(SYM_Times, Acoef, nup)),
                                ds_call2(SYM_Times, Bcoef, nu0));
                Expr* nueq = expr_new_function(expr_new_symbol(SYM_Equal),
                                 (Expr*[]){ eval_and_free(lhs), expr_new_integer(0) }, 2);
                Expr* nubody = if_run_applied(nueq, nun, xv);   /* nu(x, C[1], C[2]) */
                if (nubody) {
                    /* one basis solution: C[1]->1, C[2]->0 */
                    Expr* nu1 = ds_subst(ds_subst(nubody, ds_const(1), expr_new_integer(1)),
                                         ds_const(2), expr_new_integer(0));
                    if (!ds_is_zero(nu1) && !ds_has_head(nu1, SYM_Integrate)) {
                        Expr* fac = eval_and_free(ds_call1(SYM_Exp,
                                        ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(Ia))));
                        mu = eval_and_free(ds_call2(SYM_Times, nu1, fac));
                    } else expr_free(nu1);
                }
            } else { expr_free(Acoef); expr_free(Bcoef); }
        }
        expr_free(phi);
    }

    expr_free(a); expr_free(b); expr_free(c);
    expr_free(ax); expr_free(Ia); expr_free(I);
    if (mu && (if_too_big(mu) || ds_is_zero(mu))) { expr_free(mu); mu = NULL; }
    return mu;
}

/* ===================================================================== *
 *  Stage 2/3 (mu(x,y'), mu(y,y')) — BLOCKED on the reduced-ODE solver.     *
 *                                                                         *
 *  mu(x,y')=𝓕·μ̃(x) with 𝓕 from Lemma 3 (Υ=Φ_y (2.35)) and μ̃ from Lemma 2  *
 *  (2.30)-(2.34); mu(y,y') is the point-swap y<->x of that (2.95).  This   *
 *  session implemented and VERIFIED the μ-search for Cases A/C/D (Kamke    *
 *  226 → μ=y'; Kamke 136/66 → valid first integrals R with A(R)=0), but it *
 *  yields 0 new corpus solves: the reduced first integrals R==C[1] are     *
 *  NON-ELEMENTARY first-order ODEs (y'=Sqrt[x²y²+2C], y'=Tan[C+Log[x-y]])  *
 *  that neither our cascade nor — verified — Maple/Mathematica close in    *
 *  elementary explicit form; those CAS return them implicitly.  So Stage 2 *
 *  is BLOCKED on (a) a non-elementary/implicit first-order solver, or (b) a *
 *  policy decision to emit the reduced first integral R(x,y[x],y'[x])==C[1] *
 *  as an implicit answer.  The full Cases A–F + Lemma-2 μ̃ recovery are     *
 *  transcribed in DSOLVE_PLAN.md M18 for that follow-up.                   *
 * ===================================================================== */
static Expr* ifactor_mu_x_y1(const Expr* Phi, const char* xv, const char* Yn, const char* Pn) {
    (void)Phi; (void)xv; (void)Yn; (void)Pn; return NULL;
}
static Expr* ifactor_mu_y_y1(const Expr* Phi, const char* xv, const char* Yn, const char* Pn) {
    (void)Phi; (void)xv; (void)Yn; (void)Pn; return NULL;
}

/* ===================================================================== *
 *  try / builtin / init                                                 *
 * ===================================================================== */

/* Turn a candidate mu into an explicit solution body, or NULL.  mu consumed. */
static Expr* ifactor_from_mu(const DSolveProblem* P, const Expr* Phi, Expr* mu,
                             const char* xv, const char* Yn, const char* Pn,
                             const char* yname) {
    if (!mu) return NULL;
    if_dbg("mu", mu);
    Expr* R = ifactor_build_R(Phi, mu, xv, Yn, Pn);
    expr_free(mu);
    if (!R) { if_dbg("build_R FAILED", NULL); return NULL; }
    if (!ifactor_R_ok(R, Phi, xv, Yn, Pn)) { if_dbg("A(R)!=0", R); expr_free(R); return NULL; }
    if_dbg("R", R);
    Expr* body = ifactor_solve_from_R(P, R, xv, Yn, Pn, yname);
    if (!body) if_dbg("solve_from_R FAILED", NULL);
    expr_free(R);
    return body;
}

Expr** dsolve_ifactor_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 2) return NULL;
    const char* xv    = P->ind_names[0];
    const char* yname = P->fun_names[0];
    const char* Yn = intern_symbol("DSolve`ifY");
    const char* Pn = intern_symbol("DSolve`ifP");

    uint64_t memo_h = expr_hash(P->eq_residuals[0]);
    if_memo_sync(eval_toplevel_id());
    if (if_memo_seen(memo_h)) return NULL;

    Expr* F = dsolve_solve_top_derivative(P, 2);      /* y'' == F(x, y, y') */
    if (!F) return NULL;
    Expr* Phi = ds_subst(expr_copy(F), ds_make_funcapp(yname, 1, xv), expr_new_symbol(Pn));
    Phi = ds_subst(Phi, ds_make_funcapp(yname, 0, xv), expr_new_symbol(Yn));
    expr_free(F);
    if (if_too_big(Phi) || if_has_undef_fn(Phi)) { expr_free(Phi); return NULL; }

    /* The reducible-mu method targets NONLINEAR 2nd-order ODEs.  A linear ODE
     * y'' == b(x)y' + c(x)y (Phi affine in y, y') is the domain of the dedicated
     * linear methods (Euler/Kovacic/SpecialFunctionForm earlier; Frobenius after),
     * and must be declined here for TWO reasons: (1) correctness — this is not the
     * method's class; (2) termination — Case B of mu(x,y) reduces to a LINEAR
     * nu-ODE which it re-solves through the cascade; without this gate that nu-ODE
     * would re-enter ifactor, spawn its own nu-ODE, and recurse without bound
     * (hitting $RecursionLimit and poisoning any caller that recurses through
     * DSolve, e.g. SpecialFunctionForm's normal-form pre-pass).  The gate makes the
     * nu-ODE decline immediately so the linear specialists finish it. */
    {
        Expr* dY = if_d(expr_copy(Phi), Yn);
        Expr* dP = if_d(expr_copy(Phi), Pn);
        bool linear = ds_free_of(dY, Yn) && ds_free_of(dY, Pn)
                   && ds_free_of(dP, Yn) && ds_free_of(dP, Pn);
        expr_free(dY); expr_free(dP);
        if (linear) { expr_free(Phi); return NULL; }
    }

    g_if_deadline = time(NULL) + 6;
    Expr* body = NULL;

    /* try the three mu-forms cheapest-first; each found mu funnels through the
     * shared build-R / verify / solve pipeline. */
    if (!body && !if_expired()) body = ifactor_from_mu(P, Phi, ifactor_mu_xy(Phi, xv, Yn, Pn), xv, Yn, Pn, yname);
    if (!body && !if_expired()) body = ifactor_from_mu(P, Phi, ifactor_mu_x_y1(Phi, xv, Yn, Pn), xv, Yn, Pn, yname);
    if (!body && !if_expired()) body = ifactor_from_mu(P, Phi, ifactor_mu_y_y1(Phi, xv, Yn, Pn), xv, Yn, Pn, yname);

    expr_free(Phi);
    if (!body) { if_memo_add(memo_h); return NULL; }
    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_ifactor(Expr* res) {
    return dsolve_method_builtin(res, dsolve_ifactor_try);
}

void dsolve_ifactor_init(void) {
    symtab_add_builtin("DSolve`ReducibleIntegratingFactor", builtin_dsolve_ifactor);
    symtab_get_def("DSolve`ReducibleIntegratingFactor")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`ReducibleIntegratingFactor",
        "DSolve`ReducibleIntegratingFactor[eqn, y, x] solves a nonlinear "
        "second-order ODE y'' == Phi(x, y, y') by finding an integrating factor "
        "mu of a restricted form (mu(x,y), mu(x,y'), or mu(y,y')), reconstructing "
        "the first integral R(x,y,y') == C[1], and solving that first-order ODE "
        "(Cheb-Terrab & Roche 1999).");
}
