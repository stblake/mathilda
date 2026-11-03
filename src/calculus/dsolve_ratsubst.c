/*
 * dsolve_ratsubst.c — DSolve`RationalSubstitution.
 *
 * A homogeneous second-order linear ODE  y'' + P y' + Q y == 0  with RATIONAL
 * coefficients that no direct rational method closes may become a solvable
 * equation after a rational change of the independent variable that the direct
 * methods never try.  Two families, each a classical reduction:
 *
 *  (1) Power argument t = x^k.  When x P and x^2 Q are rational functions of x^k
 *      (k symbolic -- read off a power x^k in the coefficients -- or k = 2, 3, 4),
 *      y(x) = Y(x^k) with
 *          Y'' + (k - 1 + u)/(k t) Y' + v/(k^2 t^2) Y == 0,
 *          u = x P, v = x^2 Q  (expressed in t).
 *      An equation symmetric under x -> -x (P odd, Q even) is the k = 2 case: its
 *      four singular points 0, +-x1, Infinity fold onto three (Riemann P in t), and
 *      a coefficient polynomial in x^k (a1 + b1 x^k)/(a0 + b0 x^k) becomes the
 *      Gauss / Kummer / Bessel form in t.
 *  (2) Inversion x = x0 + 1/t at a finite IRREGULAR singular point x0 (a pole of
 *      P of order >= 2 or of Q of order >= 3, read off the linear factors of the
 *      denominators): the irregular point moves to Infinity, where the confluent
 *      rows (Kummer / Laplace / Whittaker / Bessel) live,
 *          Y'' + (2/t - P(x0 + 1/t)/t^2) Y' + Q(x0 + 1/t)/t^4 Y == 0,
 *      e.g. x^3 y'' + 2 x y' - y == 0  ->  t Y'' + (2 - 2 t) Y' - Y == 0 (Kummer),
 *      (b x + a)^4 y'' + y == 0 -> constant coefficients after y = (b x + a) w.
 *
 * The transformed equation (denominators cleared) is handed to the scalar cascade
 * (TimeConstrained, self-recursion blocked); the answer is composed with t = g(x)
 * -- PowerExpand-ed first, the raw composition as a fallback -- and gated by
 * dsolve_num_ok_complex on the ORIGINAL equation (complex sample points, the
 * parameters instantiated), so a mis-read transform can only decline.
 *
 * BRIEF rules for a recursive method: top-level only (g_dsolve_depth <= 1), a
 * re-entry guard registered with tc_register_unwind_counter, a wall-clock deadline
 * over all candidates, a per-top-level decline memo, and the killswitch
 * MATHILDA_NO_RATSUBST.  It runs LATE (after DSolve`RiemannP), so it only turns an
 * UNEVAL into an answer.
 */
#include "dsolve_common.h"
#include "dsolve.h"               /* g_dsolve_depth */
#include "../core.h"              /* tc_register_unwind_counter, leaf_count_internal */
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>
#include <time.h>

static int g_rs_active = 0;          /* re-entry guard */
static time_t g_rs_deadline;
static bool rs_expired(void) { return time(NULL) >= g_rs_deadline; }

#define RS_MEMO_SLOTS 32
static uint64_t rs_epoch = 0;
static int rs_memo_n = 0;
static uint64_t rs_memo[RS_MEMO_SLOTS];
static void rs_memo_sync(uint64_t tid) { if (tid != rs_epoch) { rs_epoch = tid; rs_memo_n = 0; } }
static bool rs_memo_seen(uint64_t h) { for (int i = 0; i < rs_memo_n; i++) if (rs_memo[i] == h) return true; return false; }
static void rs_memo_add(uint64_t h) { if (rs_memo_n < RS_MEMO_SLOTS && !rs_memo_seen(h)) rs_memo[rs_memo_n++] = h; }

static Expr* rs_ct(Expr* e) { return eval_and_free(ds_call1("Cancel", ds_call1("Together", e))); }
static Expr* rs_inv(Expr* e) { return ds_call2(SYM_Power, e, expr_new_integer(-1)); }

static bool rs_true(Expr* e) {   /* e consumed */
    bool r = e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_True;
    expr_free(e);
    return r;
}

/* DSolve[eqn, Y[t], t] (Method -> method unless NULL) under TimeConstrained -> the
 * applied-form body, or NULL.  eqn borrowed. */
static Expr* rs_run_applied(const Expr* eqn, const char* fname, const char* tv, int secs,
                            const char* method) {
    Expr* call = method
        ? expr_new_function(expr_new_symbol(SYM_DSolve),
              (Expr*[]){ expr_copy((Expr*)eqn), ds_call1(fname, expr_new_symbol(tv)), expr_new_symbol(tv),
                         expr_new_function(expr_new_symbol(SYM_Rule),
                             (Expr*[]){ expr_new_symbol(intern_symbol("Method")),
                                        expr_new_string(method) }, 2) }, 4)
        : expr_new_function(expr_new_symbol(SYM_DSolve),
              (Expr*[]){ expr_copy((Expr*)eqn), ds_call1(fname, expr_new_symbol(tv)), expr_new_symbol(tv) }, 3);
    Expr* r = eval_and_free(expr_new_function(expr_new_symbol("TimeConstrained"),
                  (Expr*[]){ call, expr_new_integer(secs), expr_new_symbol(intern_symbol("$Aborted")) }, 3));
    Expr* body = NULL;
    if (head_is(r, SYM_List) && r->data.function.arg_count >= 1) {
        Expr* inner = r->data.function.args[0];
        if (head_is(inner, SYM_List))
            for (size_t k = 0; k < inner->data.function.arg_count && !body; k++) {
                Expr* rule = inner->data.function.args[k];
                if (head_is(rule, SYM_Rule) && rule->data.function.arg_count == 2) {
                    Expr* rl = rule->data.function.args[0];
                    if (rl->type == EXPR_FUNCTION && rl->data.function.head->type == EXPR_SYMBOL &&
                        rl->data.function.head->data.symbol.name == fname)
                        body = expr_copy(rule->data.function.args[1]);
                }
            }
    }
    expr_free(r);
    return body;
}

/* Y'' + At Y' + Bt Y == 0 in tv (At, Bt consumed): denominators cleared (a nested
 * DSolve does not re-normalise), solved, composed with t = g (borrowed) and gated on
 * the original y'' + Pc y' + Qc y == 0.  Returns the candidate in x or NULL. */
static Expr* rs_solve(const DSolveProblem* P, const Expr* Pc, const Expr* Qc,
                      Expr* At, Expr* Bt, const Expr* g, const char* tv) {
    const char* xv = P->ind_names[0];
    const char* Yf = intern_symbol("DSolve`rsY");
    At = rs_ct(At); Bt = rs_ct(Bt);
    if (!ds_free_of(At, xv) || !ds_free_of(Bt, xv) ||
        !ds_is_rational_in(At, tv) || !ds_is_rational_in(Bt, tv) ||
        leaf_count_internal(At, true) + leaf_count_internal(Bt, true) > 160) {
        expr_free(At); expr_free(Bt); return NULL;
    }
    Expr* m = eval_and_free(ds_call2(SYM_PolynomialLCM,
                  eval_and_free(ds_call1(SYM_Denominator, expr_copy(At))),
                  eval_and_free(ds_call1(SYM_Denominator, expr_copy(Bt)))));
    At = eval_and_free(ds_call1(SYM_Cancel, ds_call2(SYM_Times, At, expr_copy(m))));
    Bt = eval_and_free(ds_call1(SYM_Cancel, ds_call2(SYM_Times, Bt, expr_copy(m))));
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
        (Expr*[]){ ds_call2(SYM_Plus, ds_call2(SYM_Times, m, ds_make_funcapp(Yf, 2, tv)),
                     ds_call2(SYM_Plus,
                       ds_call2(SYM_Times, At, ds_make_funcapp(Yf, 1, tv)),
                       ds_call2(SYM_Times, Bt, ds_make_funcapp(Yf, 0, tv)))),
                   expr_new_integer(0) }, 2);
    /* The direct special-function rows first (cheap, non-recursive), then the
     * whole cascade with what is left of the budget: an earlier cascade method can
     * spend the sub-budget on a Riemann P-equation the direct row closes at once. */
    static const char* const direct[] = { "SpecialFunctionForm", "RiemannP",
                                          "SpecialFunctionTransform", NULL };
    Expr* Yt = NULL;
    g_rs_active++;
    for (int i = 0; i <= 3 && !Yt; i++) {
        int secs = (int)(g_rs_deadline - time(NULL));
        if (secs > 3) secs = 3;
        if (secs < 1) break;
        Yt = rs_run_applied(eq, Yf, tv, secs, direct[i]);
    }
    g_rs_active--;
    expr_free(eq);
    if (!Yt) return NULL;
    if (ds_has_head(Yt, SYM_DSolve) || leaf_count_internal(Yt, true) > 800) { expr_free(Yt); return NULL; }
    Expr* raw = ds_subst(Yt, expr_new_symbol(tv), expr_copy((Expr*)g));     /* y = Y(g) */
    Expr* pe = eval_and_free(ds_call1(SYM_PowerExpand, expr_copy(raw)));
    if (dsolve_num_ok_complex(pe, Pc, Qc, xv)) { expr_free(raw); return pe; }
    expr_free(pe);
    if (dsolve_num_ok_complex(raw, Pc, Qc, xv)) return raw;
    expr_free(raw);
    return NULL;
}

/* (1) t = x^k: u = x P, v = x^2 Q must be rational in t.  k borrowed. */
static Expr* rs_try_power(const DSolveProblem* P, const Expr* Pc, const Expr* Qc,
                          const Expr* k, const char* tv) {
    const char* xv = P->ind_names[0];
    Expr* root = ds_call2(SYM_Power, expr_new_symbol(tv), rs_inv(expr_copy((Expr*)k)));  /* t^(1/k) */
    Expr* uv[2];
    uv[0] = ds_call2(SYM_Times, expr_new_symbol(xv), expr_copy((Expr*)Pc));
    uv[1] = ds_call2(SYM_Times, ds_call2(SYM_Power, expr_new_symbol(xv), expr_new_integer(2)),
                     expr_copy((Expr*)Qc));
    bool ok = true, varies = false;
    for (int i = 0; i < 2; i++) {
        uv[i] = ds_subst(uv[i], expr_new_symbol(xv), expr_copy(root));
        uv[i] = rs_ct(eval_and_free(ds_call1(SYM_PowerExpand, uv[i])));
        if (!ds_free_of(uv[i], xv) || !ds_is_rational_in(uv[i], tv)) ok = false;
        else if (!ds_free_of(uv[i], tv)) varies = true;
    }
    expr_free(root);
    Expr* out = NULL;
    if (ok && varies) {   /* constant u, v is Euler's equation: owned by DSolve`Euler */
        Expr* kt = ds_call2(SYM_Times, expr_copy((Expr*)k), expr_new_symbol(tv));
        Expr* At = ds_call2(SYM_Times,
                       ds_call2(SYM_Plus, ds_call2(SYM_Subtract, expr_copy((Expr*)k), expr_new_integer(1)),
                                uv[0]),
                       rs_inv(kt));
        Expr* Bt = ds_call2(SYM_Times, uv[1],
                       rs_inv(ds_call2(SYM_Power,
                           ds_call2(SYM_Times, expr_copy((Expr*)k), expr_new_symbol(tv)), expr_new_integer(2))));
        Expr* g = eval_and_free(ds_call2(SYM_Power, expr_new_symbol(xv), expr_copy((Expr*)k)));
        out = rs_solve(P, Pc, Qc, At, Bt, g, tv);
        expr_free(g);
    } else { expr_free(uv[0]); expr_free(uv[1]); }
    return out;
}

/* (2) x = x0 + 1/t.  x0 borrowed. */
static Expr* rs_try_invert(const DSolveProblem* P, const Expr* Pc, const Expr* Qc,
                           const Expr* x0, const char* tv) {
    const char* xv = P->ind_names[0];
    Expr* xs = ds_call2(SYM_Plus, expr_copy((Expr*)x0), rs_inv(expr_new_symbol(tv)));   /* x0 + 1/t */
    Expr* Ps = ds_subst(expr_copy((Expr*)Pc), expr_new_symbol(xv), expr_copy(xs));
    Expr* Qs = ds_subst(expr_copy((Expr*)Qc), expr_new_symbol(xv), xs);
    Expr* At = ds_call2(SYM_Subtract,
                   ds_call2(SYM_Times, expr_new_integer(2), rs_inv(expr_new_symbol(tv))),
                   ds_call2(SYM_Times, Ps,
                       rs_inv(ds_call2(SYM_Power, expr_new_symbol(tv), expr_new_integer(2)))));
    Expr* Bt = ds_call2(SYM_Times, Qs,
                   rs_inv(ds_call2(SYM_Power, expr_new_symbol(tv), expr_new_integer(4))));
    Expr* g = eval_and_free(rs_inv(ds_call2(SYM_Subtract, expr_new_symbol(xv), expr_copy((Expr*)x0))));
    Expr* out = rs_solve(P, Pc, Qc, At, Bt, g, tv);
    expr_free(g);
    return out;
}

/* Roots of the linear factors of den (in x) whose multiplicity is >= minmult;
 * appended to xs (owned) up to cap. */
static void rs_irregular_points(const Expr* den, const char* xv, int minmult,
                                Expr** xs, int* n, int cap) {
    Expr* fl = eval_and_free(ds_call1("FactorList", expr_copy((Expr*)den)));
    if (head_is(fl, SYM_List)) {
        for (size_t i = 0; i < fl->data.function.arg_count && *n < cap; i++) {
            Expr* pair = fl->data.function.args[i];
            if (!head_is(pair, SYM_List) || pair->data.function.arg_count != 2) continue;
            Expr* fac = pair->data.function.args[0];
            Expr* me = pair->data.function.args[1];
            if (me->type != EXPR_INTEGER || me->data.integer < minmult) continue;
            Expr* d = ds_d(expr_copy(fac), expr_new_symbol(xv));
            bool linear = !ds_free_of(fac, xv) && ds_free_of(d, xv);
            if (linear) {
                Expr* r = rs_ct(ds_call2(SYM_Times, expr_new_integer(-1),
                              ds_call2(SYM_Times,
                                  ds_subst(expr_copy(fac), expr_new_symbol(xv), expr_new_integer(0)),
                                  rs_inv(expr_copy(d)))));
                bool dup = false;
                for (int j = 0; j < *n; j++) if (expr_eq(xs[j], r)) dup = true;
                if (dup) expr_free(r); else xs[(*n)++] = r;
            }
            expr_free(d);
        }
    }
    expr_free(fl);
}

/* Symbolic exponents k of x^k in e (k not a number), deduplicated, into ks. */
static void rs_symbolic_exponents(const Expr* e, const char* xv, Expr** ks, int* n, int cap) {
    if (!e || e->type != EXPR_FUNCTION || *n >= cap) return;
    if (head_is((Expr*)e, SYM_Power) && e->data.function.arg_count == 2 &&
        e->data.function.args[0]->type == EXPR_SYMBOL &&
        e->data.function.args[0]->data.symbol.name == xv &&
        ds_free_of(e->data.function.args[1], xv)) {
        Expr* k = e->data.function.args[1];
        Expr* nq = eval_and_free(ds_call1("NumberQ", expr_copy(k)));
        if (!rs_true(nq)) {
            /* the positive form: strip a leading -1 (x^-k comes from 1/x^k) */
            Expr* kp = expr_copy(k);
            if (head_is(kp, SYM_Times) && kp->data.function.arg_count >= 2 &&
                kp->data.function.args[0]->type == EXPR_INTEGER && kp->data.function.args[0]->data.integer < 0)
                kp = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), kp));
            bool dup = false;
            for (int j = 0; j < *n; j++) if (expr_eq(ks[j], kp)) dup = true;
            if (dup) expr_free(kp); else ks[(*n)++] = kp;
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        rs_symbolic_exponents(e->data.function.args[i], xv, ks, n, cap);
}

Expr** dsolve_ratsubst_try(DSolveProblem* P, size_t* nbranch) {
    if (getenv("MATHILDA_NO_RATSUBST")) return NULL;
    if (g_dsolve_depth > 1 || g_rs_active) return NULL;
    if (P->nfun != 1 || P->neq != 1 || P->max_order[0] != 2) return NULL;
    Expr* Pc; Expr* Qc;
    if (!dsolve_second_order_PQ(P, &Pc, &Qc)) return NULL;
    const char* xv = P->ind_names[0];
    /* rational coefficients, or rational in a symbolic power x^k (the power path) */
    bool rat = ds_is_rational_in(Pc, xv) && ds_is_rational_in(Qc, xv);
    Expr* ks[4]; int nk = 0;
    if (!rat) {
        rs_symbolic_exponents(Pc, xv, ks, &nk, 4);
        rs_symbolic_exponents(Qc, xv, ks, &nk, 4);
    }
    if ((!rat && nk == 0) || leaf_count_internal(Pc, true) + leaf_count_internal(Qc, true) > 160) {
        for (int i = 0; i < nk; i++) expr_free(ks[i]);
        expr_free(Pc); expr_free(Qc); return NULL;
    }
    uint64_t h = expr_hash(P->eq_residuals[0]);
    rs_memo_sync(eval_toplevel_id());
    if (rs_memo_seen(h)) {
        for (int i = 0; i < nk; i++) expr_free(ks[i]);
        expr_free(Pc); expr_free(Qc); return NULL;
    }

    g_rs_deadline = time(NULL) + 6;
    const char* tv = intern_symbol("DSolve`rst");
    Expr* body = NULL;

    /* (2) inversion at a finite irregular singular point */
    if (rat) {
        Expr* x0s[4]; int nx = 0;
        Expr* dP = eval_and_free(ds_call1(SYM_Denominator, ds_call1(SYM_Together, expr_copy(Pc))));
        Expr* dQ = eval_and_free(ds_call1(SYM_Denominator, ds_call1(SYM_Together, expr_copy(Qc))));
        rs_irregular_points(dP, xv, 2, x0s, &nx, 4);
        rs_irregular_points(dQ, xv, 3, x0s, &nx, 4);
        expr_free(dP); expr_free(dQ);
        for (int i = 0; i < nx; i++) {
            if (!body && !rs_expired()) body = rs_try_invert(P, Pc, Qc, x0s[i], tv);
            expr_free(x0s[i]);
        }
    }
    /* (1) power argument t = x^k: symbolic exponents first, then k = 2, 3, 4 */
    for (int i = 0; i < nk; i++) {
        if (!body && !rs_expired()) body = rs_try_power(P, Pc, Qc, ks[i], tv);
        expr_free(ks[i]);
    }
    if (!body && !rs_expired() && rat) {
        for (int k = 2; k <= 4 && !body && !rs_expired(); k++) {
            Expr* ke = expr_new_integer(k);
            body = rs_try_power(P, Pc, Qc, ke, tv);
            expr_free(ke);
        }
    }
    expr_free(Pc); expr_free(Qc);
    if (!body) { rs_memo_add(h); return NULL; }
    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_ratsubst(Expr* res) {
    return dsolve_method_builtin(res, dsolve_ratsubst_try);
}

void dsolve_ratsubst_init(void) {
    tc_register_unwind_counter(&g_rs_active);   /* restored on a timeout unwind */
    symtab_add_builtin("DSolve`RationalSubstitution", builtin_dsolve_ratsubst);
    symtab_get_def("DSolve`RationalSubstitution")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`RationalSubstitution",
        "DSolve`RationalSubstitution[eqn, y, x] solves a homogeneous second-order linear "
        "ODE with rational coefficients through a rational change of the independent "
        "variable: the power argument t = x^k (k symbolic, from a power x^k in the "
        "coefficients, or 2, 3, 4) when x P and x^2 Q are rational in x^k, and the "
        "inversion x = x0 + 1/t at a finite irregular singular point x0 (moving it to "
        "Infinity, where the confluent Kummer/Bessel/Whittaker forms live). The "
        "transformed equation is solved by the cascade and the composed answer is "
        "numerically back-substituted at complex points with the parameters sampled.");
}
