/*
 * dsolve_autonomous.c — DSolve`AutonomousReduction (any order n >= 2).
 *
 * Solves the autonomous ODE (independent variable x absent)
 *     y^(n)[x] == f(y, y', ..., y^(n-1))
 * by the order reduction p = y' regarded as a function of y.  The chain
 * d/dx = p d/dy (through y) expresses each derivative in p and its y-derivatives:
 *     D[1] = p,  D[k+1] = p d/dy(D[k])
 *     (y'' = p p_y,  y''' = p^2 p_yy + p p_y^2,  ...),
 * so substituting y^(k) -> D[k] turns the equation into an order-(n-1) ODE in p(y),
 * solved by recursing into the scalar cascade.  With p = P(y, C[1..n-1]) the
 * remaining equation y' == P(y) is itself autonomous (hence separable) and is solved
 * by a second recursion, giving y(x, C[1..n]).  For n = 2 this is the classical
 * p p_y == f(y, p).
 *
 * The n integration constants must stay distinct across the two recursions (each
 * DSolve call independently names its constant C[1]), so the stage-1 constants
 * C[1..n-1] are frozen to C[2..n] (via dsolve_renumber_constants) before the stage-2
 * solve mints its fresh C[1].  A degenerate y = const (which trivially satisfies any
 * autonomous equation and would pass back-substitution — see the M4 lesson) is
 * rejected by requiring the final body to depend on x.
 *
 * Runs after ReductionOfOrder / LowerDerivativeReduction (missing-y) and the linear
 * methods in the cascade, so genuinely linear autonomous equations (y'' == -y, …) and
 * missing-y forms are claimed by the cleaner method first; this one catches the
 * nonlinear remainder such as y y'' == (y')^2  →  y = C[2] E^(C[1] x), and the
 * higher-order nonlinear autonomous families (y y''' == y' y'', 2 y y''' == y', ...).
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

/* Recursion bounding (needed only once the method reduces order 3+: a stage-2
 * quadrature such as Integrate[1/Sqrt[y Log[y]+...]] is non-elementary and makes
 * the separable sub-solve spin; the order-2 path is fast and never trips these).
 * A per-top-level decline memo collapses the evaluator's ~3x re-invocation. */
static time_t g_ar_deadline;
static bool ar_expired(void) { return time(NULL) >= g_ar_deadline; }

#define AR_MEMO_SLOTS 32
static uint64_t ar_epoch = 0;
static int ar_memo_n = 0;
static uint64_t ar_memo[AR_MEMO_SLOTS];
static void ar_memo_sync(uint64_t tid){ if(tid!=ar_epoch){ar_epoch=tid;ar_memo_n=0;} }
static bool ar_memo_seen(uint64_t h){ for(int i=0;i<ar_memo_n;i++) if(ar_memo[i]==h) return true; return false; }
static void ar_memo_add(uint64_t h){ if(ar_memo_n<AR_MEMO_SLOTS && !ar_memo_seen(h)) ar_memo[ar_memo_n++]=h; }

/* Solve DSolve[eqn, fname[varname], varname] (TimeConstrained to `secs`) and return
 * the applied-form RHS body (owned), or NULL if the sub-solve declined / timed out.
 * `eqn` is consumed. */
static Expr* run_dsolve_applied(Expr* eqn, const char* fname, const char* varname, int secs) {
    Expr* lhs  = ds_call1(fname, expr_new_symbol(varname));       /* fname[varname] */
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ eqn, lhs, expr_new_symbol(varname) }, 3);
    Expr* g = expr_new_function(expr_new_symbol("TimeConstrained"),
                  (Expr*[]){ call, expr_new_integer(secs < 1 ? 1 : secs),
                             expr_new_symbol(intern_symbol("$Aborted")) }, 3);
    Expr* r = eval_and_free(g);
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

/* Collect distinct interned symbol names in ARGUMENT position (free vars/params);
 * a function head (Log, Plus, ...) is NOT a parameter and must not be instantiated. */
static void ar_collect_params(const Expr* e, const char** names, int* n, int cap) {
    if (!e || *n >= cap) return;
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        for (int i = 0; i < *n; i++) if (names[i] == nm) return;
        names[(*n)++] = nm; return;
    }
    if (e->type == EXPR_FUNCTION) {
        if (e->data.function.head && e->data.function.head->type != EXPR_SYMBOL)
            ar_collect_params(e->data.function.head, names, n, cap);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            ar_collect_params(e->data.function.args[i], names, n, cap);
    }
}

/* Numeric self-verify of the implicit first integral.  dsolve_run_implicit's verify
 * only substitutes y'[x] and so passes VACUOUSLY for order n>=2 (y''..y^(n) stay
 * free).  So reconstruct EVERY derivative from the reduction chain (y'=p, y''=p p_y,
 * ... = D[k] with pfun -> Function[{y}, pbody]), substitute y^(k)[x] -> D[k] and
 * y[x] -> y into the ORIGINAL autonomous residual, and require it ~0 at sample y and
 * generic constants -- so a wrong reduction is dropped (0 FAIL by construction). */
static bool ar_num_ok(const DSolveProblem* P, const Expr* pbody, int n) {
    const char* xvar = P->ind_names[0];
    const char* yname = P->fun_names[0];
    const char* Ysym = intern_symbol("DSolve`arY");
    const char* pfun = intern_symbol("DSolve`arp");
    /* pfunc = Function[{Ysym}, pbody] */
    Expr* pfunc = expr_new_function(expr_new_symbol(SYM_Function), (Expr*[]){
        expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ expr_new_symbol(Ysym) }, 1),
        expr_copy((Expr*)pbody) }, 2);
    /* chain D[1..n] (abstract pfun) */
    Expr** D = malloc((size_t)(n + 1) * sizeof(Expr*));
    D[0] = NULL;
    D[1] = ds_make_funcapp(pfun, 0, Ysym);
    for (int k = 1; k < n; k++)
        D[k + 1] = eval_and_free(ds_call2(SYM_Times, ds_make_funcapp(pfun, 0, Ysym),
                       ds_d(expr_copy(D[k]), expr_new_symbol(Ysym))));
    /* residual with y^(k)[x] -> (D[k] /. pfun -> pfunc), y[x] -> Ysym */
    Expr* R = expr_copy(P->eq_residuals[0]);
    for (int k = n; k >= 1; k--) {
        Expr* yk = ds_subst(expr_copy(D[k]), expr_new_symbol(pfun), expr_copy(pfunc));
        yk = eval_and_free(yk);
        R = ds_subst(R, ds_make_funcapp(yname, k, xvar), yk);
    }
    R = ds_subst(R, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Ysym));
    for (int k = 1; k <= n; k++) expr_free(D[k]);
    free(D); expr_free(pfunc);
    /* instantiate the constants C[1..n] and every remaining free parameter at
     * distinct generic reals; sample Ysym; require |R| ~ 0. */
    for (int k = 1; k <= n; k++)
        R = ds_subst(R, ds_const(k), expr_new_real(0.37 + 0.11 * k));
    {
        const char* skip[] = { Ysym, xvar, intern_symbol("E"), intern_symbol("Pi"),
            intern_symbol("I"), intern_symbol("EulerGamma"), intern_symbol("Degree") };
        const int nskip = (int)(sizeof(skip) / sizeof(skip[0]));
        const char* syms[64]; int ns = 0;
        ar_collect_params(R, syms, &ns, 64);
        int pi = 0;
        for (int i = 0; i < ns; i++) {
            bool sk = false;
            for (int j = 0; j < nskip; j++) if (syms[i] == skip[j]) { sk = true; break; }
            if (sk) continue;
            R = ds_subst(R, expr_new_symbol(syms[i]), expr_new_real(0.43 + 0.19 * (pi++)));
        }
    }
    const double ys[] = { 1.3, 0.7, 2.1, 1.7 };
    int small = 0, big = 0;
    for (int i = 0; i < 4; i++) {
        Expr* e = ds_subst(expr_copy(R), expr_new_symbol(Ysym), expr_new_real(ys[i]));
        e = eval_and_free(ds_call1("Abs", eval_and_free(ds_call1("N", e))));
        double m = (e->type == EXPR_REAL) ? e->data.real
                 : (e->type == EXPR_INTEGER) ? (double)e->data.integer : NAN;
        expr_free(e);
        if (isnan(m) || !isfinite(m)) continue;
        if (m < 1e-6) small++; else if (m > 1e-3) big++;
    }
    expr_free(R);
    return small >= 2 && big == 0;
}

/* Shared stage-1 reduction (both try-fns): gate + missing-x pre-gate + solve-for-top
 * + autonomous check + the p=y'(y) derivative chain + solve the reduced order-(n-1)
 * ODE + freeze its constants C[1..n-1] -> C[2..n].  Returns pbody = p(y, C[2..n])
 * (owned) with *n_out set, or NULL (not autonomous / reduced ODE did not close).
 * Sets g_ar_deadline.  No memo here -- each try-fn memoizes its own decline. */
static Expr* ar_reduce(DSolveProblem* P, int* n_out) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    int n = P->max_order[0];
    if (n < 2) return NULL;
    const char* xvar  = P->ind_names[0];
    const char* yname = P->fun_names[0];
    const char* Ysym = intern_symbol("DSolve`arY");
    const char* Mask = intern_symbol("DSolve`arMask");

    /* missing-x pre-gate: mask y^(k)[x] and require no bare x survives (else the
     * coefficients depend on x -> not autonomous; solving for y^(n) could then hang). */
    {
        Expr* Rm = expr_copy(P->eq_residuals[0]);
        for (int k = n; k >= 1; k--)
            Rm = ds_subst(Rm, ds_make_funcapp(yname, k, xvar), expr_new_symbol(Mask));
        Rm = ds_subst(Rm, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Ysym));
        bool has_x = ds_contains(Rm, xvar);
        expr_free(Rm);
        if (has_x) return NULL;
    }

    g_ar_deadline = time(NULL) + 5;

    Expr* F = dsolve_solve_top_derivative(P, n);          /* y^(n) == F(x, y, ..., y^(n-1)) */
    if (!F) return NULL;

    /* autonomous check on F: mask y^(1..n-1), y; require no explicit x and genuine y
     * dependence (else the missing-y case ReductionOfOrder/LowerDerivative owns it). */
    {
        Expr* Ft = expr_copy(F);
        for (int k = n - 1; k >= 1; k--)
            Ft = ds_subst(Ft, ds_make_funcapp(yname, k, xvar), expr_new_symbol(Mask));
        Ft = ds_subst(Ft, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Ysym));
        bool bad = !ds_free_of(Ft, xvar) || !ds_contains(Ft, Ysym);
        expr_free(Ft);
        if (bad) { expr_free(F); return NULL; }
    }

    /* Reduction p = y'(y): D[1]=p, D[k+1]=p d/dy(D[k]) (y''=p p_y, y'''=p^2 p_yy+p p_y^2,
     * ...).  D[n] involves p^(n-1) -> substituting y^(k)->D[k] gives an order-(n-1) ODE
     * in p(y).  For n=2 this is the classical p p_y == f(y,p). */
    const char* pfun = intern_symbol("DSolve`arp");
    Expr** D = malloc((size_t)(n + 1) * sizeof(Expr*));
    D[0] = NULL;
    D[1] = ds_make_funcapp(pfun, 0, Ysym);
    for (int k = 1; k < n; k++)
        D[k + 1] = eval_and_free(ds_call2(SYM_Times, ds_make_funcapp(pfun, 0, Ysym),
                       ds_d(expr_copy(D[k]), expr_new_symbol(Ysym))));

    Expr* rhs1 = expr_copy(F);
    expr_free(F);
    for (int k = n - 1; k >= 1; k--)
        rhs1 = ds_subst(rhs1, ds_make_funcapp(yname, k, xvar), expr_copy(D[k]));
    rhs1 = ds_subst(rhs1, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Ysym));
    Expr* lhs1 = expr_copy(D[n]);
    for (int k = 1; k <= n; k++) expr_free(D[k]);
    free(D);
    Expr* eq1  = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ lhs1, rhs1 }, 2);
    Expr* pbody = run_dsolve_applied(eq1, pfun, Ysym, (int)(g_ar_deadline - time(NULL)));
    if (!pbody) return NULL;

    /* Freeze stage-1 constants C[1..n-1] -> C[2..n] (leaving C[1] for stage 2 / the
     * implicit quadrature constant); C[k] stays a recognised constant (the fast path
     * the elliptic stage-2 integrand needs -- see the M58 note). */
    int off = 1;
    pbody = dsolve_renumber_constants(pbody, n - 1, &off);
    if (ar_expired()) { expr_free(pbody); return NULL; }
    *n_out = n;
    return pbody;
}

/* Explicit: reduce, then (elementary-quadrature only) solve y'==p(y) by separation. */
Expr** dsolve_autonomous_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    uint64_t memo_h = expr_hash(P->eq_residuals[0]);
    ar_memo_sync(eval_toplevel_id());
    if (ar_memo_seen(memo_h)) return NULL;

    int n = 0;
    Expr* pbody = ar_reduce(P, &n);
    if (!pbody) { ar_memo_add(memo_h); return NULL; }

    const char* xvar  = P->ind_names[0];
    const char* yname = P->fun_names[0];
    const char* Ysym  = intern_symbol("DSolve`arY");

    /* Stage-2 quadrature-spin guard: y'==p(y) is separable with quadrature
     * Integrate[1/p,y], which SPINS uninterruptibly for a non-elementary integrand (a
     * Log under a radical, or a y-denominator radical).  Decline here (fast) -- the
     * IMPLICIT companion below picks these up and returns the inert first integral.
     * The order-2 forms that solve are untouched (rational p; constant-denominator
     * radical; the elliptic quartic passes and fast-declines in stage 2 as before). */
    {
        Expr* den = eval_and_free(ds_call1(SYM_Denominator,
                        ds_call1(SYM_Together, expr_copy(pbody))));
        bool nonconst_denom = !ds_free_of(den, Ysym);
        expr_free(den);
        if (nonconst_denom || ds_has_head(pbody, intern_symbol("Log"))) {
            expr_free(pbody); ar_memo_add(memo_h); return NULL;
        }
    }

    Expr* pOfY = ds_subst(pbody, expr_new_symbol(Ysym), ds_make_funcapp(yname, 0, xvar));
    Expr* eq2  = expr_new_function(expr_new_symbol(SYM_Equal),
                     (Expr*[]){ ds_make_funcapp(yname, 1, xvar), pOfY }, 2);
    Expr* ybody = run_dsolve_applied(eq2, yname, xvar, (int)(g_ar_deadline - time(NULL)));
    if (!ybody) { ar_memo_add(memo_h); return NULL; }
    if (ds_free_of(ybody, xvar)) { expr_free(ybody); ar_memo_add(memo_h); return NULL; }

    Expr** out = malloc(sizeof(Expr*));
    out[0] = ybody;
    *nbranch = 1;
    return out;
}

/* Implicit companion (M59): where the explicit stage-2 quadrature is non-elementary,
 * return the first integral Integrate[1/p, y] - x == C[1] with the integral kept INERT
 * (Inactive[Integrate] -- dodges the uninterruptible Integrate cascade; verified by the
 * FTC/implicit-function rule).  A numeric self-verify guarantees correctness. */
Expr** dsolve_autonomous_implicit_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    uint64_t memo_h = expr_hash(P->eq_residuals[0]) ^ 0x9E3779B97F4A7C15ull;
    ar_memo_sync(eval_toplevel_id());
    if (ar_memo_seen(memo_h)) return NULL;

    int n = 0;
    Expr* pbody = ar_reduce(P, &n);
    if (!pbody) { ar_memo_add(memo_h); return NULL; }

    if (!ar_num_ok(P, pbody, n)) { expr_free(pbody); ar_memo_add(memo_h); return NULL; }

    const char* xvar  = P->ind_names[0];
    const char* yname = P->fun_names[0];
    const char* Ysym  = intern_symbol("DSolve`arY");

    /* G = Inactive[Integrate][1/pbody, Ysym] /. Ysym -> y[x]  -  x */
    Expr* invp = eval_and_free(ds_call2(SYM_Power, pbody, expr_new_integer(-1)));  /* consumes pbody */
    Expr* inactHead = expr_new_function(expr_new_symbol(SYM_Inactive),
                          (Expr*[]){ expr_new_symbol(SYM_Integrate) }, 1);
    Expr* integ = expr_new_function(inactHead, (Expr*[]){ invp, expr_new_symbol(Ysym) }, 2);
    integ = ds_subst(integ, expr_new_symbol(Ysym), ds_make_funcapp(yname, 0, xvar));
    Expr* G = eval_and_free(ds_call2(SYM_Subtract, integ, expr_new_symbol(xvar)));

    Expr** out = malloc(sizeof(Expr*));
    out[0] = G;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_autonomous(Expr* res) {
    Expr* r = dsolve_method_builtin(res, dsolve_autonomous_try);
    if (!r) r = dsolve_method_builtin_implicit(res, dsolve_autonomous_implicit_try);
    return r;
}

void dsolve_autonomous_init(void) {
    symtab_add_builtin("DSolve`AutonomousReduction", builtin_dsolve_autonomous);
    symtab_get_def("DSolve`AutonomousReduction")->attributes
        |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`AutonomousReduction",
        "DSolve`AutonomousReduction[eqn, y, x] solves an autonomous ODE of order n >= 2, "
        "y^(n) == f(y, y', ..., y^(n-1)) (no explicit x), by the reduction p = y'(y): the "
        "derivative chain gives an order-(n-1) ODE in p(y), solved by recursion, then "
        "y' == p(y) by separation.");
}
