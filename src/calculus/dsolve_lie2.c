/*
 * dsolve_lie2.c — DSolve`SecondOrderSymmetry.
 *
 * Solves a second-order ODE  y''[x] == Phi(x, y, y')  by finding a Lie point
 * symmetry  X = xi(x,y) d_x + eta(x,y) d_y  and reducing the order with it.
 *
 * 1. SYMMETRY.  A point symmetry satisfies the second-prolongation determining
 *    equation (Cheb-Terrab, Duarte & da Mota, physics/9703082, Eq. 2; here in the
 *    standard prolongation form with p = y'):
 *      S2 = eta_xx + (2 eta_xy - xi_xx) p + (eta_yy - 2 xi_xy) p^2 - xi_yy p^3
 *           + (eta_y - 2 xi_x - 3 xi_y p) Phi
 *           - xi Phi_x - eta Phi_y - (eta_x + (eta_y - xi_x) p - xi_y p^2) Phi_p
 *           == 0.
 *    With xi, eta general bivariate polynomials in (x,y) of total degree <= d and
 *    Phi rational in (x,y,p), S2 clears to a polynomial identity in (x,y,p) whose
 *    coefficients are linear/homogeneous in the ansatz coefficients; the
 *    determining system's NullSpace is a basis of admissible (xi, eta) — the same
 *    machinery dsolve_lie.c uses at first order, lifted to the second prolongation
 *    and split over {x, y, p}.  (SymPy/Maple `symgen` way=3.)
 *
 * 2. REDUCTION.  Given one symmetry, build canonical coordinates (r, s) with
 *    X r == 0, X s == 1, so the ODE becomes autonomous in s and reduces to a
 *    first-order ODE  dq/dr == F(r, q)  in  q = ds/dr.  That first-order ODE is
 *    solved by recursing into the scalar cascade; then s = Integrate[q, r] + C[2]
 *    and  s(x,y) == S(r(x,y)) + C[2]  is solved for y, giving the explicit general
 *    solution y(x, C[1], C[2]).  Every branch is back-substitution verified by
 *    dsolve_run (no inert heads), so a wrong symmetry or bad inversion is dropped.
 *
 * Runs after the nonlinear-2nd-order specialists (ReductionOfOrder,
 * AutonomousReduction, Liouville) and before the series fallback.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

/* Wall-clock (CPU) budget: the reduction chains recursive DSolve/Solve/Integrate
 * calls whose cost is data-dependent; without a deadline a pathological rational
 * Phi can spin for tens of seconds.  Bound one top-level SecondOrderSymmetry call
 * to ~6 s of CPU, checked between steps, and decline on expiry (never a wrong
 * answer — the cascade falls through to the series fallback). */
static time_t g_l2_deadline;
static bool l2_expired(void) { return time(NULL) >= g_l2_deadline; }

/* Per-toplevel decline memo.  The evaluator's fixed-point loop re-invokes a
 * declining builtin several times per DSolve call (the equation re-normalizes each
 * pass, so the dispatcher's own fail-memo misses it); without this, a ~6 s lie2
 * decline is paid 3x.  Record the equation hash on decline and short-circuit a
 * repeat within the same top-level evaluation. */
#define L2_MEMO_SLOTS 32
static uint64_t l2_memo_epoch = 0;
static int      l2_memo_n = 0;
static uint64_t l2_memo[L2_MEMO_SLOTS];
static void l2_memo_sync(uint64_t tid) { if (tid != l2_memo_epoch) { l2_memo_epoch = tid; l2_memo_n = 0; } }
static bool l2_memo_seen(uint64_t h) {
    for (int i = 0; i < l2_memo_n; i++) if (l2_memo[i] == h) return true;
    return false;
}
static void l2_memo_add(uint64_t h) {
    if (l2_memo_n < L2_MEMO_SLOTS && !l2_memo_seen(h)) l2_memo[l2_memo_n++] = h;
}

/* ---- small local helpers ---- */

static Expr* l2_powi(Expr* base, int n) {
    return eval_and_free(ds_call2(SYM_Power, base, expr_new_integer(n)));
}

static void l2_count(const Expr* e, long* acc, long budget) {
    if (*acc > budget || !e) return;
    (*acc)++;
    if (e->type == EXPR_FUNCTION) {
        l2_count(e->data.function.head, acc, budget);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            l2_count(e->data.function.args[i], acc, budget);
    }
}
static bool l2_too_big(const Expr* e) {
    long acc = 0; l2_count(e, &acc, 60000); return acc > 60000;
}

/* True if `e` carries an undefined function (head is a symbol with no builtin and
 * no down-values, or an inert Derivative): its symmetry search / reduction is
 * non-elementary and must be declined rather than spun on. */
static bool l2_has_undef_fn(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        if (h->data.symbol.name == SYM_Derivative) return true;
        SymbolDef* d = symtab_lookup(h->data.symbol.name);
        if (d && !d->builtin_func && !d->down_values) return true;
    } else if (l2_has_undef_fn(h)) {
        return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (l2_has_undef_fn(e->data.function.args[i])) return true;
    return false;
}

/* Solve DSolve[eqn, fname[varname], varname] -> applied-form RHS body (owned) or
 * NULL.  `eqn` consumed.  The sub-solve is wrapped in TimeConstrained: the reduced
 * / characteristic ODE can hit a pre-existing slow path in the first-order cascade
 * (e.g. the explicit inversion inside DSolve`Separable), and lie2 must decline that
 * symmetry rather than spin.  A single top-level DSolve is not itself wrapped in
 * TimeConstrained, so this is the outermost one in the common case; the bound is
 * generous enough (5 s) that no genuine sub-solve is cut. */
static Expr* l2_run_applied(Expr* eqn, const char* fname, const char* varname) {
    Expr* lhs  = ds_call1(fname, expr_new_symbol(varname));
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ eqn, lhs, expr_new_symbol(varname) }, 3);
    Expr* guarded = expr_new_function(expr_new_symbol("TimeConstrained"),
                     (Expr*[]){ call, expr_new_integer(2),
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

/* Evaluate `e` under a TimeConstrained bound (returns $Aborted on timeout), so a
 * pre-existing slow path in Solve/Integrate/NullSpace makes lie2 decline the
 * symmetry rather than spin.  `e` consumed; result owned. */
static Expr* l2_beval(Expr* e, int secs) {
    Expr* g = expr_new_function(expr_new_symbol("TimeConstrained"),
                  (Expr*[]){ e, expr_new_integer(secs),
                             expr_new_symbol(intern_symbol("$Aborted")) }, 3);
    return eval_and_free(g);
}

/* Bounded Solve[eq, var] and Integrate[e, v]. */
static Expr* l2_solve_b(Expr* eq, const char* var) {
    return l2_beval(expr_new_function(expr_new_symbol(SYM_Solve),
                        (Expr*[]){ eq, expr_new_symbol(var) }, 2), 2);
}
static Expr* l2_integrate_b(Expr* e, const char* v) {
    return l2_beval(expr_new_function(expr_new_symbol(SYM_Integrate),
                        (Expr*[]){ e, expr_new_symbol(v) }, 2), 2);
}
/* Bounded Simplify: returns a $Aborted-guarded Simplify, falling back to the
 * unsimplified input on timeout (still correct, just less tidy).  e consumed. */
static Expr* l2_simplify_b(Expr* e) {
    Expr* keep = expr_copy(e);
    Expr* r = l2_beval(expr_new_function(expr_new_symbol(SYM_Simplify), (Expr*[]){ e }, 1), 2);
    if (r && r->type == EXPR_SYMBOL && r->data.symbol.name == intern_symbol("$Aborted")) {
        expr_free(r); return keep;
    }
    expr_free(keep); return r;
}

/* Solve[expr == rhs, var] and return the first solution value for `var` (owned)
 * or NULL.  `expr` and `rhs` consumed. */
static Expr* l2_solve_eq(Expr* expr, Expr* rhs, const char* var) {
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ expr, rhs }, 2);
    Expr* sol = l2_solve_b(eq, var);
    size_t n = 0;
    Expr** vals = dsolve_extract_solutions(sol, var, &n);
    expr_free(sol);
    Expr* out = NULL;
    if (vals && n >= 1) { out = vals[0]; for (size_t i=1;i<n;i++) expr_free(vals[i]); }
    free(vals);
    return out;
}

/* Solve[expr == rhs, var] and return ALL solution values (owned array of *n), or
 * NULL.  `expr` and `rhs` consumed. */
static Expr** l2_solve_all(Expr* expr, Expr* rhs, const char* var, size_t* n) {
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ expr, rhs }, 2);
    Expr* sol = l2_solve_b(eq, var);
    Expr** vals = dsolve_extract_solutions(sol, var, n);
    expr_free(sol);
    return vals;
}

static Expr* l2_d(Expr* e, const char* v) { return ds_d(e, expr_new_symbol(v)); }
static Expr* l2_d2(Expr* e, const char* a, const char* b) { return l2_d(l2_d(e, a), b); }

/* Numeric back-substitution guard: substitute y -> Function[{x}, body] into the
 * original residual and sample it at several (x, C[1], C[2]) points.  The symbolic
 * verify in dsolve_run KEEPS an undecidable residual (Solve policy), so a wrong
 * reduction branch with a symbolically-intractable residual would slip through;
 * this rejects one that is clearly non-zero.  A candidate passes iff >= 2 sample
 * points give a finite |residual| < 1e-6 and NO finite point gives |residual| >
 * 1e-3 (a singular/non-finite point is skipped, not held against it). */
static double l2_abs_at(const Expr* R, const char* xv, double xval, double c1, double c2) {
    Expr* e = expr_copy((Expr*)R);
    e = ds_subst(e, ds_const(1), expr_new_real(c1));
    e = ds_subst(e, ds_const(2), expr_new_real(c2));
    e = ds_subst(e, expr_new_symbol(xv), expr_new_real(xval));
    e = eval_and_free(ds_call1("Abs", e));
    double m = (e && e->type == EXPR_REAL)    ? e->data.real
             : (e && e->type == EXPR_INTEGER) ? (double)e->data.integer
             : NAN;                                     /* non-numeric -> skip */
    expr_free(e);
    return m;
}
static bool l2_num_ok(const DSolveProblem* P, const Expr* body,
                      const char* xv, const char* yname) {
    Expr* R = expr_copy(P->eq_residuals[0]);
    Expr* b0 = expr_copy((Expr*)body);
    Expr* b1 = ds_d(expr_copy((Expr*)body), expr_new_symbol(xv));
    Expr* b2 = l2_d2(expr_copy((Expr*)body), xv, xv);
    R = ds_subst(R, ds_make_funcapp(yname, 2, xv), b2);
    R = ds_subst(R, ds_make_funcapp(yname, 1, xv), b1);
    R = ds_subst(R, ds_make_funcapp(yname, 0, xv), b0);
    const double xs[] = { 1.7, 2.3, 1.15, 0.6, 3.1 };
    const double c1s[] = { 1.181, 1.4, 2.25, 0.7, 1.9 };
    const double c2s[] = { 0.714, 0.75, 1.6, 1.3, 0.4 };
    int small = 0, big = 0;
    for (int i = 0; i < 5; i++) {
        double m = l2_abs_at(R, xv, xs[i], c1s[i], c2s[i]);
        if (isnan(m) || !isfinite(m)) continue;
        if (m < 1e-6) small++;
        else if (m > 1e-3) big++;
    }
    expr_free(R);
    return small >= 2 && big == 0;
}

/* ---- the second-order determining expression S2(xi, eta) ---- */
static Expr* l2_S2(const Expr* xi, const Expr* eta, const Expr* Phi,
                   const char* xv, const char* Yn, const char* Pn) {
    Expr* p = expr_new_symbol(Pn);
    Expr* eta_xx = l2_d2(expr_copy((Expr*)eta), xv, xv);
    Expr* eta_xy = l2_d2(expr_copy((Expr*)eta), xv, Yn);
    Expr* eta_yy = l2_d2(expr_copy((Expr*)eta), Yn, Yn);
    Expr* eta_x  = l2_d(expr_copy((Expr*)eta), xv);
    Expr* eta_y  = l2_d(expr_copy((Expr*)eta), Yn);
    Expr* xi_xx = l2_d2(expr_copy((Expr*)xi), xv, xv);
    Expr* xi_xy = l2_d2(expr_copy((Expr*)xi), xv, Yn);
    Expr* xi_yy = l2_d2(expr_copy((Expr*)xi), Yn, Yn);
    Expr* xi_x  = l2_d(expr_copy((Expr*)xi), xv);
    Expr* xi_y  = l2_d(expr_copy((Expr*)xi), Yn);
    Expr* Phi_x = l2_d(expr_copy((Expr*)Phi), xv);
    Expr* Phi_y = l2_d(expr_copy((Expr*)Phi), Yn);
    Expr* Phi_p = l2_d(expr_copy((Expr*)Phi), Pn);

    /* A = eta_xx + (2 eta_xy - xi_xx) p + (eta_yy - 2 xi_xy) p^2 - xi_yy p^3 */
    Expr* c1 = ds_call2(SYM_Subtract, ds_call2(SYM_Times, expr_new_integer(2), eta_xy), xi_xx);
    Expr* c2 = ds_call2(SYM_Subtract, eta_yy, ds_call2(SYM_Times, expr_new_integer(2), xi_xy));
    Expr* A = ds_call2(SYM_Plus, eta_xx,
                ds_call2(SYM_Plus, ds_call2(SYM_Times, c1, expr_copy(p)),
                  ds_call2(SYM_Plus,
                    ds_call2(SYM_Times, c2, l2_powi(expr_copy(p), 2)),
                    ds_call2(SYM_Times, expr_new_integer(-1),
                      ds_call2(SYM_Times, xi_yy, l2_powi(expr_copy(p), 3))))));
    /* B = (eta_y - 2 xi_x - 3 xi_y p) Phi */
    Expr* bcoef = ds_call2(SYM_Subtract,
                    ds_call2(SYM_Subtract, expr_copy(eta_y),
                        ds_call2(SYM_Times, expr_new_integer(2), expr_copy(xi_x))),
                    ds_call2(SYM_Times, expr_new_integer(3),
                        ds_call2(SYM_Times, expr_copy(xi_y), expr_copy(p))));
    Expr* B = ds_call2(SYM_Times, bcoef, expr_copy((Expr*)Phi));
    /* C = - xi Phi_x - eta Phi_y */
    Expr* Cc = ds_call2(SYM_Subtract,
                ds_call2(SYM_Times, expr_new_integer(-1),
                    ds_call2(SYM_Times, expr_copy((Expr*)xi), Phi_x)),
                ds_call2(SYM_Times, expr_copy((Expr*)eta), Phi_y));
    /* eta1 = eta_x + (eta_y - xi_x) p - xi_y p^2 ; Dd = - eta1 Phi_p */
    Expr* eta1 = ds_call2(SYM_Plus, eta_x,
                   ds_call2(SYM_Subtract,
                     ds_call2(SYM_Times, ds_call2(SYM_Subtract, eta_y, xi_x), expr_copy(p)),
                     ds_call2(SYM_Times, xi_y, l2_powi(expr_copy(p), 2))));
    Expr* Dd = ds_call2(SYM_Times, expr_new_integer(-1), ds_call2(SYM_Times, eta1, Phi_p));
    expr_free(p);
    Expr* S = ds_call2(SYM_Plus, A, ds_call2(SYM_Plus, B, ds_call2(SYM_Plus, Cc, Dd)));
    return eval_and_free(S);
}

/* ---- polynomial-ansatz symmetry search ---- */

static const char* l2_cname(int k) {
    char buf[32]; snprintf(buf, sizeof buf, "DSolve`l2B%d", k); return intern_symbol(buf);
}

static Expr* l2_poly(Expr** coeffs, int degree, const char* xv, const char* Yn) {
    Expr* sum = expr_new_integer(0);
    int k = 0;
    for (int total = 0; total <= degree; total++)
        for (int i = total; i >= 0; i--) {
            int j = total - i;
            Expr* term = expr_copy(coeffs[k++]);
            if (i > 0) term = ds_call2(SYM_Times, term, l2_powi(expr_new_symbol(xv), i));
            if (j > 0) term = ds_call2(SYM_Times, term, l2_powi(expr_new_symbol(Yn), j));
            sum = ds_call2(SYM_Plus, sum, term);
        }
    return eval_and_free(sum);
}

static Expr* l2_xyp(const char* xv, const char* Yn, const char* Pn) {
    return expr_new_function(expr_new_symbol(SYM_List),
        (Expr*[]){ expr_new_symbol(xv), expr_new_symbol(Yn), expr_new_symbol(Pn) }, 3);
}

static bool l2_trueQ(Expr* e) {
    bool t = (e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_True);
    expr_free(e); return t;
}

/* NullSpace basis of point symmetries at ansatz `degree`.  Returns malloc'd array
 * of *count pairs laid out [xi0,eta0,xi1,eta1,...] (2*count exprs), or NULL. */
static Expr** l2_find_symmetries(const Expr* Phi, int degree,
                                 const char* xv, const char* Yn, const char* Pn,
                                 size_t* count) {
    *count = 0;
    Expr* tg  = eval_and_free(ds_call1(SYM_Together, expr_copy((Expr*)Phi)));
    Expr* num = eval_and_free(ds_call1(SYM_Numerator, expr_copy(tg)));
    Expr* den = eval_and_free(ds_call1(SYM_Denominator, tg));
    bool r1 = l2_trueQ(eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
                  (Expr*[]){ num, l2_xyp(xv,Yn,Pn) }, 2)));
    bool r2 = l2_trueQ(eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
                  (Expr*[]){ den, l2_xyp(xv,Yn,Pn) }, 2)));
    if (!(r1 && r2)) return NULL;

    int nmon = (degree + 1) * (degree + 2) / 2;
    int nc   = 2 * nmon;
    Expr** csym = malloc((size_t)nc * sizeof(Expr*));
    for (int k = 0; k < nc; k++) csym[k] = expr_new_symbol(l2_cname(k + 1));
    Expr* xi  = l2_poly(&csym[0],    degree, xv, Yn);
    Expr* eta = l2_poly(&csym[nmon], degree, xv, Yn);
    Expr* S = l2_S2(xi, eta, Phi, xv, Yn, Pn);
    expr_free(xi); expr_free(eta);
    if (l2_too_big(S)) { expr_free(S); for (int k=0;k<nc;k++) expr_free(csym[k]); free(csym); return NULL; }

    Expr* tS = l2_beval(ds_call1(SYM_Together, S), 2);
    Expr* nS = eval_and_free(ds_call1(SYM_Numerator, tS));
    if (l2_too_big(nS)) { expr_free(nS); for (int k=0;k<nc;k++) expr_free(csym[k]); free(csym); return NULL; }
    Expr* cl = l2_beval(expr_new_function(expr_new_symbol(SYM_CoefficientList),
                   (Expr*[]){ nS, l2_xyp(xv,Yn,Pn) }, 2), 2);
    if (l2_expired()) { expr_free(cl); for (int k=0;k<nc;k++) expr_free(csym[k]); free(csym); return NULL; }
    Expr* forms = eval_and_free(ds_call1(SYM_Flatten, cl));

    Expr** vargs = malloc((size_t)nc * sizeof(Expr*));
    for (int k = 0; k < nc; k++) vargs[k] = expr_copy(csym[k]);
    Expr* vars = expr_new_function(expr_new_symbol(SYM_List), vargs, (size_t)nc);
    free(vargs);
    for (int k = 0; k < nc; k++) expr_free(csym[k]);
    free(csym);

    Expr* M = l2_beval(expr_new_function(expr_new_symbol("Outer"),
                  (Expr*[]){ expr_new_symbol("Coefficient"), forms, vars }, 3), 2);
    if (l2_expired() || !ds_free_of(M, xv) || !ds_free_of(M, Yn) || !ds_free_of(M, Pn)) { expr_free(M); return NULL; }
    Expr* ns = l2_beval(ds_call1("NullSpace", M), 2);

    Expr** out = NULL; size_t n = 0;
    if (ns && ns->type == EXPR_FUNCTION) {
        size_t nv = ns->data.function.arg_count;
        out = malloc((nv ? nv : 1) * 2 * sizeof(Expr*));
        for (size_t rr = 0; rr < nv; rr++) {
            Expr* v = ns->data.function.args[rr];
            if (!v || v->type != EXPR_FUNCTION || v->data.function.arg_count != (size_t)nc) continue;
            Expr** c = v->data.function.args;
            out[2*n]   = l2_poly(&c[0],    degree, xv, Yn);
            out[2*n+1] = l2_poly(&c[nmon], degree, xv, Yn);
            n++;
        }
        if (n == 0) { free(out); out = NULL; }
    }
    expr_free(ns);
    *count = n;
    return out;
}

/* ---- order reduction by canonical coordinates ---- */

static Expr* l2_reduce(const DSolveProblem* P, const Expr* xi, const Expr* eta,
                       const Expr* Phi, const char* xv, const char* Yn,
                       const char* Pn, const char* yname) {
    const char* rsym = intern_symbol("DSolve`l2r");
    const char* qsym = intern_symbol("DSolve`l2q");
    const char* qfun = intern_symbol("DSolve`l2qf");
    const char* ctmp = intern_symbol("DSolve`l2c");

    if (l2_expired()) return NULL;
    bool xi0  = ds_is_zero(xi);
    bool eta0 = ds_is_zero(eta);
    if (xi0 && eta0) return NULL;

    /* 1. canonical coordinates r(x,y), s(x,y). */
    Expr *r_expr = NULL, *s_expr = NULL;
    bool r_has_y;
    if (eta0) {                       /* X = xi d_x : r = y, s = Integrate[1/xi, x] */
        r_expr = expr_new_symbol(Yn); r_has_y = true;
        s_expr = l2_integrate_b(l2_powi(expr_copy((Expr*)xi), -1), xv);
    } else if (xi0) {                 /* X = eta d_y : r = x, s = Integrate[1/eta, y] */
        r_expr = expr_new_symbol(xv); r_has_y = false;
        s_expr = l2_integrate_b(l2_powi(expr_copy((Expr*)eta), -1), Yn);
    } else {
        /* r: first integral of y' == eta/xi (constant expressed in x,y). */
        Expr* crhs = eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)eta),
                         l2_powi(expr_copy((Expr*)xi), -1)));
        Expr* crhs_yx = ds_subst(crhs, expr_new_symbol(Yn), ds_make_funcapp(yname, 0, xv));
        Expr* ceq = expr_new_function(expr_new_symbol(SYM_Equal),
                        (Expr*[]){ ds_make_funcapp(yname, 1, xv), crhs_yx }, 2);
        if (l2_expired()) { expr_free(ceq); return NULL; }
        Expr* cbody = l2_run_applied(ceq, yname, xv);      /* y = B(x, C[1]) */
        if (!cbody) return NULL;
        Expr* cbody_c = ds_subst(cbody, ds_const(1), expr_new_symbol(ctmp));  /* B(x, ctmp) */
        r_expr = l2_solve_eq(cbody_c, expr_new_symbol(Yn), ctmp);  /* ctmp = r(x,y) */
        if (!r_expr) return NULL;
        r_has_y = ds_contains(r_expr, Yn);
        Expr* Yofr = l2_solve_eq(expr_copy(r_expr), expr_new_symbol(rsym), Yn); /* y = Y(x,r) */
        Expr* invxi = l2_powi(expr_copy((Expr*)xi), -1);
        if (Yofr) invxi = ds_subst(invxi, expr_new_symbol(Yn), Yofr);
        Expr* Sr = l2_integrate_b(invxi, xv);
        s_expr = ds_subst(Sr, expr_new_symbol(rsym), expr_copy(r_expr));
    }
    if (!s_expr || ds_has_head(s_expr, SYM_Integrate) || l2_too_big(s_expr)) {
        expr_free(r_expr); expr_free(s_expr); return NULL;
    }

    /* 2. q = ds/dr and dq/dr (with y'' -> Phi). */
    Expr* rx = l2_d(expr_copy(r_expr), xv);
    Expr* ry = l2_d(expr_copy(r_expr), Yn);
    Expr* sx = l2_d(expr_copy(s_expr), xv);
    Expr* sy = l2_d(expr_copy(s_expr), Yn);
    Expr* drdx = ds_call2(SYM_Plus, expr_copy(rx), ds_call2(SYM_Times, expr_copy(ry), expr_new_symbol(Pn)));
    Expr* dsdx = ds_call2(SYM_Plus, sx, ds_call2(SYM_Times, sy, expr_new_symbol(Pn)));
    Expr* q = eval_and_free(ds_call2(SYM_Times, dsdx, l2_powi(expr_copy(drdx), -1)));
    expr_free(rx); expr_free(ry);

    Expr* qx = l2_d(expr_copy(q), xv);
    Expr* qy = l2_d(expr_copy(q), Yn);
    Expr* qp = l2_d(expr_copy(q), Pn);
    Expr* numr = ds_call2(SYM_Plus, qx,
                   ds_call2(SYM_Plus, ds_call2(SYM_Times, qy, expr_new_symbol(Pn)),
                     ds_call2(SYM_Times, qp, expr_copy((Expr*)Phi))));
    Expr* dqdr = eval_and_free(ds_call2(SYM_Times, numr, l2_powi(expr_copy(drdx), -1)));
    expr_free(drdx);
    if (l2_too_big(dqdr) || l2_too_big(q)) { expr_free(dqdr); expr_free(q); expr_free(r_expr); expr_free(s_expr); return NULL; }

    /* 3. eliminate p (via q) then the leftover variable (via r) -> F(r,q). */
    Expr* Pofq = l2_solve_eq(expr_copy(q), expr_new_symbol(qsym), Pn);
    expr_free(q);
    if (!Pofq) { expr_free(dqdr); expr_free(r_expr); expr_free(s_expr); return NULL; }
    Expr* dqdr_p = ds_subst(dqdr, expr_new_symbol(Pn), Pofq);

    Expr* Frq;
    if (r_has_y) {
        Expr* Yofr = l2_solve_eq(expr_copy(r_expr), expr_new_symbol(rsym), Yn);
        if (!Yofr) { expr_free(dqdr_p); expr_free(r_expr); expr_free(s_expr); return NULL; }
        Frq = l2_simplify_b(ds_subst(dqdr_p, expr_new_symbol(Yn), Yofr));
        if (!ds_free_of(Frq, xv)) { expr_free(Frq); expr_free(r_expr); expr_free(s_expr); return NULL; }
    } else {
        Frq = l2_simplify_b(ds_subst(dqdr_p, expr_new_symbol(xv), expr_new_symbol(rsym)));
        if (!ds_free_of(Frq, Yn)) { expr_free(Frq); expr_free(r_expr); expr_free(s_expr); return NULL; }
    }

    /* 4. solve dq/dr == F(r,q) for q(r, C[1]). */
    if (l2_expired()) { expr_free(Frq); expr_free(r_expr); expr_free(s_expr); return NULL; }
    Expr* Frq_qf = ds_subst(Frq, expr_new_symbol(qsym), ds_make_funcapp(qfun, 0, rsym));
    Expr* qeq = expr_new_function(expr_new_symbol(SYM_Equal),
                    (Expr*[]){ ds_make_funcapp(qfun, 1, rsym), Frq_qf }, 2);
    Expr* qgen = l2_run_applied(qeq, qfun, rsym);
    if (!qgen) { expr_free(r_expr); expr_free(s_expr); return NULL; }

    /* 5. s = Integrate[Q, r] + C[2] ; relation s(x,y) == that ; solve for y. */
    Expr* Sint = l2_integrate_b(qgen, rsym);
    if (ds_has_head(Sint, SYM_Integrate)) { expr_free(Sint); expr_free(r_expr); expr_free(s_expr); return NULL; }
    Expr* Srhs = ds_call2(SYM_Plus, Sint, ds_const(2));
    Srhs = ds_subst(Srhs, expr_new_symbol(rsym), r_expr);   /* consumes r_expr */
    /* Try EVERY inversion branch and numerically verify it — the symbolic verify
     * in dsolve_run keeps an undecidable residual, so a wrong branch must be
     * rejected here. */
    if (l2_expired()) { expr_free(s_expr); expr_free(Srhs); return NULL; }
    size_t nb = 0;
    Expr** cands = l2_solve_all(s_expr, Srhs, Yn, &nb);     /* consumes s_expr, Srhs */
    Expr* body = NULL;
    for (size_t i = 0; i < nb; i++) {
        Expr* c = cands[i];
        if (!body && !ds_free_of(c, xv) && !ds_has_head(c, SYM_Solve)
            && !ds_has_head(c, SYM_Integrate) && !l2_too_big(c)
            && l2_num_ok(P, c, xv, yname))
            body = expr_copy(c);
        expr_free(c);
    }
    free(cands);
    return body;
}

Expr** dsolve_lie2_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 2) return NULL;
    const char* xv    = P->ind_names[0];
    const char* yname = P->fun_names[0];
    const char* Yn = intern_symbol("DSolve`l2Y");
    const char* Pn = intern_symbol("DSolve`l2P");

    uint64_t memo_h = expr_hash(P->eq_residuals[0]);
    l2_memo_sync(eval_toplevel_id());
    if (l2_memo_seen(memo_h)) return NULL;   /* declined this equation already */

    Expr* F = dsolve_solve_top_derivative(P, 2);      /* y'' == F(x, y, y') */
    if (!F) return NULL;
    Expr* Phi = ds_subst(expr_copy(F), ds_make_funcapp(yname, 1, xv), expr_new_symbol(Pn));
    Phi = ds_subst(Phi, ds_make_funcapp(yname, 0, xv), expr_new_symbol(Yn));
    expr_free(F);
    if (l2_too_big(Phi) || l2_has_undef_fn(Phi)) { expr_free(Phi); return NULL; }

    /* lie2 targets NONLINEAR 2nd-order ODEs.  A linear ODE y'' == A(x)y + B(x)y'
     * (Phi affine in y, y') is the domain of the dedicated linear methods
     * (Euler/Kovacic/SpecialFunction, earlier; Frobenius series, after): its
     * symmetry algebra is 8-dimensional and the determining-system search over a
     * high-degree rational coefficient is expensive with no benefit here.  Decline
     * so a Kovacic-declined linear ODE falls straight through to the series
     * fallback (its pre-existing behaviour). */
    {
        Expr* dY = l2_d(expr_copy(Phi), Yn);
        Expr* dP = l2_d(expr_copy(Phi), Pn);
        bool linear = ds_free_of(dY, Yn) && ds_free_of(dY, Pn)
                   && ds_free_of(dP, Yn) && ds_free_of(dP, Pn);
        expr_free(dY); expr_free(dP);
        if (linear) { expr_free(Phi); return NULL; }
    }

    g_l2_deadline = time(NULL) + 6;   /* ~6 s wall-clock budget per call */
    Expr* body = NULL;
    for (int degree = 1; degree <= 2 && !body && !l2_expired(); degree++) {
        size_t nsym = 0;
        Expr** syms = l2_find_symmetries(Phi, degree, xv, Yn, Pn, &nsym);
        for (size_t i = 0; i < nsym && !body && !l2_expired(); i++)
            body = l2_reduce(P, syms[2*i], syms[2*i+1], Phi, xv, Yn, Pn, yname);
        if (syms) { for (size_t i = 0; i < 2*nsym; i++) expr_free(syms[i]); free(syms); }
    }
    expr_free(Phi);
    if (!body) { l2_memo_add(memo_h); return NULL; }   /* remember the costly decline */
    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_lie2(Expr* res) {
    return dsolve_method_builtin(res, dsolve_lie2_try);
}

void dsolve_lie2_init(void) {
    symtab_add_builtin("DSolve`SecondOrderSymmetry", builtin_dsolve_lie2);
    symtab_get_def("DSolve`SecondOrderSymmetry")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`SecondOrderSymmetry",
        "DSolve`SecondOrderSymmetry[eqn, y, x] solves a second-order ODE "
        "y'' == Phi(x, y, y') by finding a Lie point symmetry (polynomial-ansatz "
        "determining system) and reducing the order via canonical coordinates.");
}
