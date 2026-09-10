/*
 * dsolve_piecewise.c — DSolve`PiecewiseForcing (M35).
 *
 * Linear ODE initial-value problems with PIECEWISE / STEP / HEAVISIDE forcing:
 *     L[y] == g(t),   g carrying Piecewise[...] / UnitStep[...],   y^(k)(t0) given.
 * These are the Boyce & DiPrima "Laplace transform" problems (§2.2.15, 1492-1500).
 * Mathilda has Piecewise/UnitStep but no LaplaceTransform, so we solve them by
 * INTERVAL CONTINUATION, which reuses the whole scalar cascade and yields a
 * back-substitution-verifiable Piecewise closed form (Mathematica's own form):
 *
 *   1. Detect a scalar linear IVP whose forcing carries Piecewise/UnitStep with
 *      finite ordered breakpoints t0 < b1 < ... < bm.
 *   2. On each interval [b_{j-1}, b_j) the forcing is a single smooth expression
 *      piece_j(t) (resolve Piecewise/UnitStep at an interior point).
 *   3. Solve interval 0 with the given ICs (recurse DSolve); at b1 read off the
 *      continuity data y,y',...,y^(n-1) (a bounded linear ODE with piecewise-
 *      continuous forcing has a C^(n-1) solution) as ICs for interval 1; repeat.
 *   4. Assemble Piecewise[{{body_0,t<b1},{body_1,b1<=t<b2},...,{body_m,t>=bm}}].
 *
 * The assembled body carries NO generated constant, so it routes through the
 * standard dsolve_run: dsolve_verify_body numeric-checks the ODE residual and
 * dsolve_fit_constants is a no-op (npar==0 -> FIT_OK).  Because dsolve_run then
 * never re-checks the ICs for a constant-free body, and its residual probe samples
 * only the first interval for Pi-breakpoint problems, this method carries its OWN
 * per-interval numeric guard pw_num_ok (residual at every interval midpoint AND
 * every IC) before returning — a mis-stitched answer can never ship.
 *
 * Runs BEFORE UndeterminedCoefficients / LinearConstantCoefficients (whose
 * variation-of-parameters leaves Integrate[UnitStep[...]*Cos[...],t] inert); a
 * gate (forcing must carry Piecewise/UnitStep AND be an IVP) keeps it from firing
 * on smooth forcing, and a re-entry guard + TimeConstrained sub-solves + wall-clock
 * deadline + per-top-level decline memo bound it (the M12/M14 pattern).
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

#define PW_MAX_BREAKS 8            /* cap intervals so work stays bounded          */

static time_t g_pw_deadline;
static bool pw_expired(void) { return time(NULL) >= g_pw_deadline; }
static int  g_pw_active = 0;       /* re-entry guard (no self-recursion)           */

/* decline memo (the evaluator re-invokes a declining builtin ~3x/call) */
#define PW_MEMO_SLOTS 32
static uint64_t pw_epoch = 0;
static int pw_memo_n = 0;
static uint64_t pw_memo[PW_MEMO_SLOTS];
static void pw_memo_sync(uint64_t tid){ if(tid!=pw_epoch){pw_epoch=tid;pw_memo_n=0;} }
static bool pw_memo_seen(uint64_t h){ for(int i=0;i<pw_memo_n;i++) if(pw_memo[i]==h) return true; return false; }
static void pw_memo_add(uint64_t h){ if(pw_memo_n<PW_MEMO_SLOTS && !pw_memo_seen(h)) pw_memo[pw_memo_n++]=h; }

/* TimeConstrained[e, secs, $Aborted]; e consumed. */
static Expr* pw_beval(Expr* e, int secs){
    Expr* g = expr_new_function(expr_new_symbol(SYM_TimeConstrained),
                  (Expr*[]){ e, expr_new_integer(secs), expr_new_symbol(intern_symbol("$Aborted")) }, 3);
    return eval_and_free(g);
}

/* N[e] as a C double, or NaN if it does not numericize (e.g. Infinity). e borrowed. */
static double pw_numval(const Expr* e){
    Expr* v = eval_and_free(ds_call1(SYM_N, expr_copy((Expr*)e)));
    double d = (v && v->type==EXPR_REAL) ? v->data.real
             : (v && v->type==EXPR_INTEGER) ? (double)v->data.integer : NAN;
    expr_free(v);
    return d;
}

/* --- breakpoint collection -------------------------------------------------- */

static void pw_add_bound(Expr*** buf, int* n, int* cap, Expr* b){
    if (*n >= *cap){ *cap = *cap ? *cap*2 : 8; *buf = realloc(*buf, (size_t)(*cap)*sizeof(Expr*)); }
    (*buf)[(*n)++] = b;   /* takes ownership */
}

/* Collect candidate breakpoints (owned) from the forcing: the root of every
 * UnitStep[arg] (arg==0), and every t-free operand of a Piecewise condition. */
static void pw_collect_bounds(const Expr* e, const char* tv, Expr*** buf, int* n, int* cap){
    if (!e || e->type != EXPR_FUNCTION) return;
    const Expr* h = e->data.function.head;
    if (h && h->type==EXPR_SYMBOL){
        const char* hn = h->data.symbol.name;
        if (hn==SYM_UnitStep && e->data.function.arg_count==1){
            Expr* sol = ds_solve(ds_call2(SYM_Equal, expr_copy(e->data.function.args[0]),
                                          expr_new_integer(0)), expr_new_symbol(tv));
            size_t ns=0; Expr** vals = dsolve_extract_solutions(sol, tv, &ns);
            for (size_t i=0;i<ns;i++) pw_add_bound(buf,n,cap, vals[i]);
            free(vals); expr_free(sol);
        } else if (hn==SYM_Inequality || hn==SYM_Less || hn==SYM_LessEqual
                   || hn==SYM_Greater || hn==SYM_GreaterEqual){
            for (size_t i=0;i<e->data.function.arg_count;i++){
                Expr* op = e->data.function.args[i];
                if (op->type==EXPR_SYMBOL){
                    const char* on=op->data.symbol.name;   /* skip relation symbols */
                    if (on==SYM_Less||on==SYM_LessEqual||on==SYM_Greater
                        ||on==SYM_GreaterEqual||on==SYM_Equal) continue;
                }
                if (ds_free_of(op, tv)) pw_add_bound(buf,n,cap, expr_copy(op));
            }
        }
    }
    for (size_t i=0;i<e->data.function.arg_count;i++)
        pw_collect_bounds(e->data.function.args[i], tv, buf,n,cap);
}

/* Numeric-filter (> t0, finite), dedupe, and sort ascending; returns the count of
 * kept breakpoints and fills breaks[] (owned) + bnum[] (their numeric keys). */
static int pw_order_breaks(Expr** cand, int nc, double t0num, Expr** breaks, double* bnum){
    int m=0;
    for (int i=0;i<nc;i++){
        double v = pw_numval(cand[i]);
        if (isnan(v) || !isfinite(v) || v <= t0num + 1e-9) continue;
        int dup=0; for (int j=0;j<m;j++) if (fabs(bnum[j]-v) < 1e-7) { dup=1; break; }
        if (dup) continue;
        if (m >= PW_MAX_BREAKS) continue;
        /* insertion sort by numeric key */
        int p=m; while (p>0 && bnum[p-1] > v){ bnum[p]=bnum[p-1]; breaks[p]=breaks[p-1]; p--; }
        bnum[p]=v; breaks[p]=expr_copy(cand[i]); m++;
    }
    return m;
}

/* --- per-interval forcing --------------------------------------------------- */

/* Resolve Piecewise/UnitStep in `e` at the interior sample value `mid` (a real),
 * returning the smooth expression active there (still a function of tv).  Non-
 * conditional structure is rebuilt verbatim. */
static Expr* pw_resolve(const Expr* e, const char* tv, double mid){
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    const Expr* h = e->data.function.head;
    if (h && h->type==EXPR_SYMBOL){
        const char* hn = h->data.symbol.name;
        if (hn==SYM_Piecewise && e->data.function.arg_count>=1){
            Expr* pairs = e->data.function.args[0];
            Expr* deflt = (e->data.function.arg_count>=2)? e->data.function.args[1] : NULL;
            if (head_is(pairs, SYM_List)){
                for (size_t i=0;i<pairs->data.function.arg_count;i++){
                    Expr* pr = pairs->data.function.args[i];
                    if (head_is(pr,SYM_List) && pr->data.function.arg_count==2){
                        Expr* c = pr->data.function.args[1];
                        Expr* ct = eval_and_free(ds_subst(expr_copy(c),
                                       expr_new_symbol(tv), expr_new_real(mid)));
                        bool istrue = (ct && ct->type==EXPR_SYMBOL && ct->data.symbol.name==SYM_True);
                        expr_free(ct);
                        if (istrue) return pw_resolve(pr->data.function.args[0], tv, mid);
                    }
                }
            }
            return deflt ? pw_resolve(deflt, tv, mid) : expr_new_integer(0);
        }
        if (hn==SYM_UnitStep && e->data.function.arg_count==1){
            double av = pw_numval(eval_and_free(ds_subst(expr_copy(e->data.function.args[0]),
                            expr_new_symbol(tv), expr_new_real(mid))));
            return expr_new_integer((!isnan(av) && av>=0.0) ? 1 : 0);
        }
    }
    Expr* nh = pw_resolve(e->data.function.head, tv, mid);
    size_t ac = e->data.function.arg_count;
    Expr** na = malloc(ac*sizeof(Expr*));
    for (size_t i=0;i<ac;i++) na[i] = pw_resolve(e->data.function.args[i], tv, mid);
    Expr* r = expr_new_function(nh, na, ac);
    free(na);
    return r;
}

/* --- recursion into the cascade -------------------------------------------- */

/* y[k]-derivative funcapp evaluated at `point`: y[point] or Derivative[k][y][point]. */
static Expr* pw_icapp(const char* yname, int k, const char* tv, const Expr* point){
    return ds_subst(ds_make_funcapp(yname, k, tv), expr_new_symbol(tv), expr_copy((Expr*)point));
}

/* Solve DSolve[{eqn, ics...}, y[t], t] (TimeConstrained) -> applied-form body or NULL.
 * `eqns` is a List (owned) of the equation + IC equations. */
static Expr* pw_run_applied(Expr* eqns, const char* yname, const char* tv){
    Expr* lhs = ds_call1(yname, expr_new_symbol(tv));
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ eqns, lhs, expr_new_symbol(tv) }, 3);
    Expr* r = pw_beval(call, 5);
    Expr* body = NULL;
    if (head_is(r, SYM_List) && r->data.function.arg_count >= 1){
        Expr* inner = r->data.function.args[0];
        if (head_is(inner, SYM_List))
            for (size_t k=0;k<inner->data.function.arg_count && !body;k++){
                Expr* rule = inner->data.function.args[k];
                if (head_is(rule, SYM_Rule) && rule->data.function.arg_count==2){
                    Expr* rl = rule->data.function.args[0];
                    if (rl->type==EXPR_FUNCTION && rl->data.function.head->type==EXPR_SYMBOL
                        && rl->data.function.head->data.symbol.name==yname)
                        body = expr_copy(rule->data.function.args[1]);
                }
            }
    }
    expr_free(r);
    /* reject an unsolved / half-evaluated body */
    if (body && (ds_has_head(body, SYM_Integrate) || ds_has_head(body, SYM_Solve)
                 || ds_has_head(body, SYM_DSolve)
                 || ds_contains(body, intern_symbol("$Aborted")))){
        expr_free(body); body = NULL;
    }
    return body;
}

/* --- verification (every interval + every IC) ------------------------------ */

static double pw_abs_at(const Expr* R, const char* tv, double tval){
    Expr* e = ds_subst(expr_copy((Expr*)R), expr_new_symbol(tv), expr_new_real(tval));
    e = eval_and_free(ds_call1("Abs", eval_and_free(ds_call1(SYM_N, e))));
    double m = (e&&e->type==EXPR_REAL)?e->data.real:(e&&e->type==EXPR_INTEGER)?(double)e->data.integer:NAN;
    expr_free(e);
    return m;
}

/* Numeric back-substitution: the assembled `body` must zero the ODE residual at an
 * interior point of EVERY interval and satisfy EVERY initial condition. */
static bool pw_num_ok(const DSolveProblem* P, const Expr* body, const char* tv,
                      const char* yname, int order, const double* mids, int nmid){
    /* ODE residual with y^(k)[t] -> D[body,{t,k}] */
    Expr* R = expr_copy(P->eq_residuals[0]);
    for (int k=order; k>=0; k--){
        Expr* d = expr_copy((Expr*)body);
        for (int i=0;i<k;i++) d = ds_d(d, expr_new_symbol(tv));
        R = ds_subst(R, ds_make_funcapp(yname, k, tv), d);
    }
    int small=0, big=0;
    for (int i=0;i<nmid;i++){
        double m = pw_abs_at(R, tv, mids[i]);
        if (isnan(m) || !isfinite(m)) continue;
        if (m < 1e-6) small++; else if (m > 1e-3) big++;
    }
    expr_free(R);
    if (big>0 || small < nmid) return false;
    /* every initial/boundary condition */
    for (size_t c=0;c<P->ncond;c++){
        if (P->conds[c].fi != 0) continue;
        Expr* d = expr_copy((Expr*)body);
        for (int i=0;i<P->conds[c].order;i++) d = ds_d(d, expr_new_symbol(tv));
        Expr* lhs = ds_subst(d, expr_new_symbol(tv), expr_copy(P->conds[c].point));
        Expr* diff = ds_call2(SYM_Subtract, lhs, expr_copy(P->conds[c].value));
        Expr* v = eval_and_free(ds_call1("Abs", eval_and_free(ds_call1(SYM_N, diff))));
        double m = (v&&v->type==EXPR_REAL)?v->data.real:(v&&v->type==EXPR_INTEGER)?(double)v->data.integer:NAN;
        expr_free(v);
        if (!isnan(m) && isfinite(m) && m > 1e-6) return false;
    }
    return true;
}

/* --- the method ------------------------------------------------------------ */

/* Free the linear-coefficient vector and forcing, then decline. */
static Expr** pw_free_cg(Expr** coeffs, int order, Expr* g){
    if (coeffs){ for (int i=0;i<=order;i++) expr_free(coeffs[i]); free(coeffs); }
    expr_free(g);
    return NULL;
}

Expr** dsolve_piecewise_try(DSolveProblem* P, size_t* nbranch){
    if (P->nfun!=1 || P->neq!=1) return NULL;
    if (P->ncond==0) return NULL;                 /* stitching needs left-endpoint ICs */
    if (g_pw_active) return NULL;

    const char* tv = P->ind_names[0];
    const char* yname = P->fun_names[0];

    /* linear + forcing */
    Expr** coeffs=NULL; Expr* g=NULL; int order=0;
    if (!dsolve_linear_coeffs(P, &coeffs, &g, &order)) return NULL;
    if (order < 1) return pw_free_cg(coeffs, order, g);

    /* forcing must carry Piecewise/UnitStep (and no DiracDelta impulse) */
    bool step = ds_has_head(g, SYM_Piecewise) || ds_has_head(g, SYM_UnitStep);
    if (!step || ds_contains(g, intern_symbol("DiracDelta")))
        return pw_free_cg(coeffs, order, g);

    /* the ICs must be a single-point IVP with exactly `order` conditions at t0 */
    if ((int)P->ncond != order) return pw_free_cg(coeffs, order, g);
    Expr* t0 = P->conds[0].point;
    for (size_t c=1;c<P->ncond;c++)
        if (!expr_eq(P->conds[c].point, t0)) return pw_free_cg(coeffs, order, g);
    double t0num = pw_numval(t0);
    if (isnan(t0num)) return pw_free_cg(coeffs, order, g);

    uint64_t h = expr_hash(P->eq_residuals[0]);
    pw_memo_sync(eval_toplevel_id());
    if (pw_memo_seen(h)) return pw_free_cg(coeffs, order, g);

    /* breakpoints */
    Expr** cand=NULL; int nc=0, cap=0;
    pw_collect_bounds(g, tv, &cand, &nc, &cap);
    Expr* breaks[PW_MAX_BREAKS]; double bnum[PW_MAX_BREAKS];
    int m = pw_order_breaks(cand, nc, t0num, breaks, bnum);
    for (int i=0;i<nc;i++) expr_free(cand[i]);
    free(cand);
    if (m == 0) return pw_free_cg(coeffs, order, g);   /* no genuine breakpoint > t0 */

    g_pw_deadline = time(NULL) + 8;

    /* interval midpoints (numeric, interior) for clause selection + verification.
     * interval 0 = (-inf, b0): midpoint below b0 (>= t0); j = 1..m-1: (b_{j-1},b_j);
     * interval m = [b_{m-1}, inf): b_{m-1}+1. */
    int nint = m + 1;
    double mids[PW_MAX_BREAKS+1];
    mids[0] = (t0num < bnum[0]) ? 0.5*(t0num + bnum[0]) : bnum[0]-0.5;
    for (int j=1;j<m;j++) mids[j] = 0.5*(bnum[j-1] + bnum[j]);
    mids[m] = bnum[m-1] + 1.0;

    /* solve interval by interval, threading continuity data forward */
    Expr** bodies = calloc((size_t)nint, sizeof(Expr*));
    Expr** cont = NULL; int ncont = 0;             /* continuity conds for next interval */
    Expr** cont_pt = NULL;                          /* the point b_{j-1} for those        */
    bool ok = true;

    for (int j=0;j<nint && ok && !pw_expired(); j++){
        /* smooth forcing on this interval */
        Expr* piece = ds_simplify(pw_resolve(g, tv, mids[j]));

        /* equation: Sum_{r=0}^{order} coeffs[r]*y^(r)[t] == piece */
        Expr* lhs=NULL;
        for (int r=0;r<=order;r++){
            Expr* term = ds_call2(SYM_Times, expr_copy(coeffs[r]), ds_make_funcapp(yname, r, tv));
            lhs = lhs ? ds_call2(SYM_Plus, lhs, term) : term;
        }
        Expr* eqn = ds_call2(SYM_Equal, lhs, piece);

        /* IC list */
        int nic = (j==0) ? (int)P->ncond : ncont;
        Expr** items = malloc((size_t)(nic+1)*sizeof(Expr*));
        int ni=0;
        items[ni++] = eqn;
        if (j==0){
            for (size_t c=0;c<P->ncond;c++)
                items[ni++] = ds_call2(SYM_Equal,
                                  pw_icapp(yname, P->conds[c].order, tv, P->conds[c].point),
                                  expr_copy(P->conds[c].value));
        } else {
            for (int d=0;d<ncont;d++)
                items[ni++] = ds_call2(SYM_Equal,
                                  pw_icapp(yname, d, tv, cont_pt[0]), expr_copy(cont[d]));
        }
        Expr* eqns = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)ni);
        free(items);

        g_pw_active++;
        Expr* body_j = pw_run_applied(eqns, yname, tv);
        g_pw_active--;
        if (!body_j){ ok=false; break; }
        bodies[j] = body_j;

        /* continuity data at b_j for the next interval (skip after the last) */
        if (cont){ for(int d=0;d<ncont;d++){ expr_free(cont[d]); } free(cont); cont=NULL; }
        if (cont_pt){ expr_free(cont_pt[0]); free(cont_pt); cont_pt=NULL; }
        if (j < m){
            ncont = order;
            cont = malloc((size_t)ncont*sizeof(Expr*));
            for (int d=0;d<ncont;d++){
                Expr* dd = expr_copy(body_j);
                for (int i=0;i<d;i++) dd = ds_d(dd, expr_new_symbol(tv));
                cont[d] = ds_simplify(ds_subst(dd, expr_new_symbol(tv), expr_copy(breaks[j])));
            }
            cont_pt = malloc(sizeof(Expr*));
            cont_pt[0] = expr_copy(breaks[j]);
        }
    }
    if (cont){ for(int d=0;d<ncont;d++){ expr_free(cont[d]); } free(cont); }
    if (cont_pt){ expr_free(cont_pt[0]); free(cont_pt); }
    for(int i=0;i<=order;i++){ expr_free(coeffs[i]); }
    free(coeffs); expr_free(g);

    if (!ok){
        for (int j=0;j<nint;j++){ if (bodies[j]) expr_free(bodies[j]); }
        free(bodies);
        for (int i=0;i<m;i++){ expr_free(breaks[i]); }
        pw_memo_add(h); return NULL;
    }

    /* assemble Piecewise[{{body_0, t<b0}, {body_j, b_{j-1}<=t<b_j}, {body_m, t>=b_{m-1}}}, 0] */
    Expr** pairs = malloc((size_t)nint*sizeof(Expr*));
    for (int j=0;j<nint;j++){
        Expr* cond;
        if (j==0)
            cond = ds_call2(SYM_Less, expr_new_symbol(tv), expr_copy(breaks[0]));
        else if (j<m)
            cond = expr_new_function(expr_new_symbol(SYM_Inequality),
                       (Expr*[]){ expr_copy(breaks[j-1]), expr_new_symbol(SYM_LessEqual),
                                  expr_new_symbol(tv), expr_new_symbol(SYM_Less),
                                  expr_copy(breaks[j]) }, 5);
        else
            cond = ds_call2(SYM_GreaterEqual, expr_new_symbol(tv), expr_copy(breaks[m-1]));
        pairs[j] = expr_new_function(expr_new_symbol(SYM_List),
                       (Expr*[]){ bodies[j], cond }, 2);   /* consumes bodies[j] */
    }
    free(bodies);
    Expr* pairList = expr_new_function(expr_new_symbol(SYM_List), pairs, (size_t)nint);
    free(pairs);
    Expr* body = expr_new_function(expr_new_symbol(SYM_Piecewise),
                     (Expr*[]){ pairList, expr_new_integer(0) }, 2);
    for (int i=0;i<m;i++){ expr_free(breaks[i]); }

    /* thorough per-interval + per-IC numeric verify (dsolve_run cannot) */
    if (!pw_num_ok(P, body, tv, yname, order, mids, nint)){
        expr_free(body); pw_memo_add(h); return NULL;
    }

    Expr** out = malloc(sizeof(Expr*));
    out[0] = body; *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_piecewise(Expr* res){ return dsolve_method_builtin(res, dsolve_piecewise_try); }

void dsolve_piecewise_init(void){
    symtab_add_builtin("DSolve`PiecewiseForcing", builtin_dsolve_piecewise);
    symtab_get_def("DSolve`PiecewiseForcing")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`PiecewiseForcing",
        "DSolve`PiecewiseForcing[{eqn, ics}, y, t] solves a linear ODE initial-value "
        "problem whose forcing is piecewise / step / Heaviside (Piecewise[...] or "
        "UnitStep[...]) by interval continuation: it solves on each interval between "
        "the forcing's breakpoints, matches y,y',...,y^(n-1) across each breakpoint, "
        "and returns a verified Piecewise solution.");
}
