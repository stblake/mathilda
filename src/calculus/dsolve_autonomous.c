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

Expr** dsolve_autonomous_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    int n = P->max_order[0];
    if (n < 2) return NULL;
    const char* xvar  = P->ind_names[0];
    const char* yname = P->fun_names[0];
    const char* Ysym = intern_symbol("DSolve`arY");   /* the reduced independent var (= y) */
    const char* Mask = intern_symbol("DSolve`arMask"); /* a marker for any y^(k), k>=1     */

    /* Cheap missing-x pre-gate, BEFORE the potentially expensive solve-for-top.
     * A genuine autonomous y^(n)=f(y,...,y^(n-1)) has no explicit x once the funcapps
     * y^(k)[x] are masked; if a bare x survives, the coefficients depend on x and this
     * is not autonomous.  Without this, solving an x-dependent linear ODE with nested
     * rational coefficients (an exact-ODE order reduction) for y^(n) blows up the
     * rational arithmetic and hangs. */
    {
        Expr* Rm = expr_copy(P->eq_residuals[0]);
        for (int k = n; k >= 1; k--)
            Rm = ds_subst(Rm, ds_make_funcapp(yname, k, xvar), expr_new_symbol(Mask));
        Rm = ds_subst(Rm, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Ysym));
        bool has_x = ds_contains(Rm, xvar);
        expr_free(Rm);
        if (has_x) return NULL;
    }

    /* decline memo + wall-clock deadline (the recursive sub-solves below can spin on
     * a non-elementary stage-2 quadrature at order 3+). */
    uint64_t memo_h = expr_hash(P->eq_residuals[0]);
    ar_memo_sync(eval_toplevel_id());
    if (ar_memo_seen(memo_h)) return NULL;
    g_ar_deadline = time(NULL) + 5;

    Expr* F = dsolve_solve_top_derivative(P, n);          /* y^(n) == F(x, y, ..., y^(n-1)) */
    if (!F) { ar_memo_add(memo_h); return NULL; }

    /* autonomous check: mask y^(1..n-1) and y, require no explicit x and genuine y
     * dependence (else it is the missing-y case ReductionOfOrder/LowerDerivative own). */
    {
        Expr* Ft = expr_copy(F);
        for (int k = n - 1; k >= 1; k--)
            Ft = ds_subst(Ft, ds_make_funcapp(yname, k, xvar), expr_new_symbol(Mask));
        Ft = ds_subst(Ft, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Ysym));
        bool bad = !ds_free_of(Ft, xvar) || !ds_contains(Ft, Ysym);
        expr_free(Ft);
        if (bad) { expr_free(F); ar_memo_add(memo_h); return NULL; }
    }

    /* Reduction p = y'(y): the chain d/dx = p d/dy through y expresses each y^(k) as
     *   D[1] = p,   D[k+1] = p * d/dy(D[k])
     *   (y'' = p p_y,  y''' = p^2 p_yy + p p_y^2,  ...),
     * so D[n] involves p^(n-1): substituting y^(k) -> D[k] gives an order-(n-1) ODE in
     * p(y).  For n = 2 this is exactly the classical p p_y == f(y,p). */
    const char* pfun = intern_symbol("DSolve`arp");
    Expr** D = malloc((size_t)(n + 1) * sizeof(Expr*));
    D[0] = NULL;
    D[1] = ds_make_funcapp(pfun, 0, Ysym);
    for (int k = 1; k < n; k++)
        D[k + 1] = eval_and_free(ds_call2(SYM_Times, ds_make_funcapp(pfun, 0, Ysym),
                       ds_d(expr_copy(D[k]), expr_new_symbol(Ysym))));

    /* Stage 1: D[n] == F with y^(k)[x] -> D[k] (k=1..n-1), y[x] -> Ysym. */
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
    if (!pbody) { ar_memo_add(memo_h); return NULL; }

    /* Freeze the stage-1 constants C[1..n-1] as the GENERATED constants C[2..n],
     * distinct from the stage-2 solve's fresh C[1] so they cannot collide.  They must
     * NOT be frozen to plain symbols: a plain symbolic parameter inside the stage-2
     * integrand sends the first-order cascade down a ~27x slower path
     * (DSolve[y'==Sqrt[a+y^4]] ~ 11s vs DSolve[y'==Sqrt[C[2]+y^4]] ~ 0.4s), and stage 2
     * is exactly the elliptic / quadrature form autonomous reduction produces.  A C[k]
     * is recognised as a constant (fast decline) and is already the final name (no
     * rename).  dsolve_renumber_constants applies the shift as one simultaneous pass. */
    int off = 1;
    pbody = dsolve_renumber_constants(pbody, n - 1, &off);   /* C[1..n-1] -> C[2..n] */
    if (ar_expired()) { expr_free(pbody); ar_memo_add(memo_h); return NULL; }

    /* Stage-2 quadrature-spin guard.  y' == p(y) is separable with quadrature
     * Integrate[1/p, y]; for a non-elementary integrand (a Log under a radical, or a
     * radical of a rational with a y-dependent denominator) Integrate SPINS
     * uninterruptibly (TimeConstrained cannot preempt it — the 3rd-order reductions
     * routinely produce such p, e.g. Sqrt[y Log y + ...]).  Decline BEFORE stage 2 in
     * those cases (fast, no wrong answer).  The order-2 forms that solve are untouched:
     * a rational p (y y''==(y')^2 -> p = C y) and a constant-denominator radical
     * (a+b(y')^2 -> p = Sqrt[(C E^(2 b y)-a)/b], the Tan/Tanh case) both pass; the
     * elliptic order-2 radical (Sqrt of a quartic) passes and fast-declines in stage 2
     * as before. */
    {
        Expr* den = eval_and_free(ds_call1(SYM_Denominator,
                        ds_call1(SYM_Together, expr_copy(pbody))));
        bool nonconst_denom = !ds_free_of(den, Ysym);
        expr_free(den);
        if (nonconst_denom || ds_has_head(pbody, intern_symbol("Log"))) {
            expr_free(pbody); ar_memo_add(memo_h); return NULL;
        }
    }

    /* Stage 2: y'[x] == P(y[x]) (autonomous, separable). */
    Expr* pOfY = ds_subst(pbody, expr_new_symbol(Ysym), ds_make_funcapp(yname, 0, xvar));
    Expr* eq2  = expr_new_function(expr_new_symbol(SYM_Equal),
                     (Expr*[]){ ds_make_funcapp(yname, 1, xvar), pOfY }, 2);
    Expr* ybody = run_dsolve_applied(eq2, yname, xvar, (int)(g_ar_deadline - time(NULL)));
    if (!ybody) { ar_memo_add(memo_h); return NULL; }

    /* Reject a degenerate x-independent body. */
    if (ds_free_of(ybody, xvar)) { expr_free(ybody); ar_memo_add(memo_h); return NULL; }

    Expr** out = malloc(sizeof(Expr*));
    out[0] = ybody;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_autonomous(Expr* res) {
    return dsolve_method_builtin(res, dsolve_autonomous_try);
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
