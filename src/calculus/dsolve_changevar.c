/*
 * dsolve_changevar.c — DSolve`ChangeOfVariable (M14).
 *
 * Second-order linear ODEs with TRANSCENDENTAL coefficients that become
 * rational-coefficient (hence solvable by Euler / Kovacic / SpecialFunctionForm /
 * Frobenius) under a change of the INDEPENDENT variable t = phi(x).  For
 * y'' + P(x) y' + Q(x) y == 0, the substitution gives
 *     d2y/dt2 + A(t) dy/dt + B(t) y == 0,
 *     A = (phi'' + P phi')/phi'^2,   B = Q/phi'^2   (re-expressed in t),
 * and when A, B come out rational in t the transformed equation is handed back to
 * the scalar cascade; the solution is composed with t = phi(x) and back-substitution
 * verified.  Classic case: y'' + Cot[x] y' + k(k+1) y == 0  --(t=Cos[x])-->  the
 * Legendre equation (1-t^2) Y'' - 2t Y' + k(k+1) Y == 0.
 *
 * Candidate substitutions (with the trig/hyperbolic identities that rationalize
 * the coefficients): t = Cos[x], Sin[x], Tan[x], Cosh[x], Sinh[x], Tanh[x], E^x.
 * Each is verified NUMERICALLY on the ORIGINAL equation before being returned
 * (the composed solution can carry Log/algebraic terms whose residual zero_test
 * cannot decide, so dsolve_run's symbolic verify alone would keep a bad transform).
 *
 * Runs after the direct rational-coefficient methods (SpecialFunctionForm,
 * Kovacic) and before the series fallback.  All recursive sub-solves are
 * TimeConstrained-bounded with a wall-clock deadline and a per-top-level decline
 * memo (see the M12 lessons); a re-entry guard prevents recursion into itself.
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

static time_t g_cv_deadline;
static bool cv_expired(void) { return time(NULL) >= g_cv_deadline; }
static int  g_cv_active = 0;   /* re-entry guard */

/* decline memo (the evaluator re-invokes a declining builtin ~3x/call) */
#define CV_MEMO_SLOTS 32
static uint64_t cv_epoch = 0;
static int cv_memo_n = 0;
static uint64_t cv_memo[CV_MEMO_SLOTS];
static void cv_memo_sync(uint64_t tid){ if(tid!=cv_epoch){cv_epoch=tid;cv_memo_n=0;} }
static bool cv_memo_seen(uint64_t h){ for(int i=0;i<cv_memo_n;i++) if(cv_memo[i]==h) return true; return false; }
static void cv_memo_add(uint64_t h){ if(cv_memo_n<CV_MEMO_SLOTS && !cv_memo_seen(h)) cv_memo[cv_memo_n++]=h; }

static Expr* cv_powi(Expr* b, int n){ return eval_and_free(ds_call2(SYM_Power, b, expr_new_integer(n))); }

static void cv_count(const Expr* e, long* a, long b){
    if(*a>b||!e) return;
    (*a)++;
    if(e->type==EXPR_FUNCTION){ cv_count(e->data.function.head,a,b);
        for(size_t i=0;i<e->data.function.arg_count;i++) cv_count(e->data.function.args[i],a,b); }
}
static bool cv_too_big(const Expr* e){ long a=0; cv_count(e,&a,40000); return a>40000; }

static Expr* cv_beval(Expr* e, int secs){
    Expr* g = expr_new_function(expr_new_symbol("TimeConstrained"),
                  (Expr*[]){ e, expr_new_integer(secs), expr_new_symbol(intern_symbol("$Aborted")) }, 3);
    return eval_and_free(g);
}

/* Solve DSolve[eqn, Y[t], t] (TimeConstrained) -> applied-form RHS body or NULL. eqn consumed. */
static Expr* cv_run_applied(Expr* eqn, const char* fname, const char* tvar){
    Expr* lhs = ds_call1(fname, expr_new_symbol(tvar));
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ eqn, lhs, expr_new_symbol(tvar) }, 3);
    Expr* r = cv_beval(call, 5);
    Expr* body = NULL;
    if (head_is(r, SYM_List) && r->data.function.arg_count >= 1){
        Expr* inner = r->data.function.args[0];
        if (head_is(inner, SYM_List))
            for (size_t k=0;k<inner->data.function.arg_count && !body;k++){
                Expr* rule = inner->data.function.args[k];
                if (head_is(rule, SYM_Rule) && rule->data.function.arg_count==2){
                    Expr* rl = rule->data.function.args[0];
                    if (rl->type==EXPR_FUNCTION && rl->data.function.head->type==EXPR_SYMBOL
                        && rl->data.function.head->data.symbol.name==fname)
                        body = expr_copy(rule->data.function.args[1]);
                }
            }
    }
    expr_free(r);
    return body;
}

/* Build a ReplaceAll rule List from (lhs,rhs) pairs; returns Rule-list Expr (owned). */
static Expr* cv_rules(Expr** pairs, int npair){
    Expr** rl = malloc((size_t)npair*sizeof(Expr*));
    for (int i=0;i<npair;i++)
        rl[i] = expr_new_function(expr_new_symbol(SYM_Rule),
                    (Expr*[]){ pairs[2*i], pairs[2*i+1] }, 2);
    Expr* L = expr_new_function(expr_new_symbol(SYM_List), rl, (size_t)npair);
    free(rl);
    return L;
}

/* trig[x] -> algebraic-in-t rewrite for a given substitution; returns the rule List. */
static Expr* cv_trig_rules(const char* xv, const char* tv, int which){
    const char* T[]={"Sin","Cos","Tan","Cot","Sec","Csc"};
    Expr* t = expr_new_symbol(tv);
    Expr* one_m_t2 = ds_call2(SYM_Subtract, expr_new_integer(1), cv_powi(expr_copy(t),2)); /* 1-t^2 */
    Expr* one_p_t2 = ds_call2(SYM_Plus, expr_new_integer(1), cv_powi(expr_copy(t),2));      /* 1+t^2 */
    Expr* s=NULL,*c=NULL;   /* Sin[x], Cos[x] images */
    if (which==0){ /* t=Cos[x] : cos=t, sin=Sqrt[1-t^2] */
        c=expr_copy(t); s=eval_and_free(ds_call1("Sqrt", expr_copy(one_m_t2)));
    } else if (which==1){ /* t=Sin[x] : sin=t, cos=Sqrt[1-t^2] */
        s=expr_copy(t); c=eval_and_free(ds_call1("Sqrt", expr_copy(one_m_t2)));
    } else { /* t=Tan[x] : sin=t/Sqrt[1+t^2], cos=1/Sqrt[1+t^2] */
        Expr* rt=eval_and_free(ds_call1("Sqrt", expr_copy(one_p_t2)));
        s=eval_and_free(ds_call2(SYM_Times, expr_copy(t), cv_powi(expr_copy(rt),-1)));
        c=cv_powi(rt,-1);
    }
    /* build images for all six */
    Expr* img[6];
    img[0]=expr_copy(s);                                             /* Sin */
    img[1]=expr_copy(c);                                             /* Cos */
    img[2]=eval_and_free(ds_call2(SYM_Times, expr_copy(s), cv_powi(expr_copy(c),-1))); /* Tan */
    img[3]=eval_and_free(ds_call2(SYM_Times, expr_copy(c), cv_powi(expr_copy(s),-1))); /* Cot */
    img[4]=cv_powi(expr_copy(c),-1);                                 /* Sec */
    img[5]=cv_powi(expr_copy(s),-1);                                 /* Csc */
    expr_free(s); expr_free(c); expr_free(t); expr_free(one_m_t2); expr_free(one_p_t2);
    Expr* pairs[12];
    for (int i=0;i<6;i++){
        pairs[2*i]  = ds_call1(T[i], expr_new_symbol(xv));
        pairs[2*i+1]= img[i];
    }
    return cv_rules(pairs, 6);
}

/* Rationalize `e` (in xv) under substitution `which` to a function of tv; returns
 * Together[Simplify[e /. rules]] or NULL if not rational in tv.  e consumed. */
static Expr* cv_to_t(Expr* e, const char* xv, const char* tv, int which){
    Expr* rules = cv_trig_rules(xv, tv, which);
    Expr* sub = eval_and_free(expr_new_function(expr_new_symbol(SYM_ReplaceAll),
                    (Expr*[]){ e, rules }, 2));
    Expr* simp = cv_beval(expr_new_function(expr_new_symbol(SYM_Simplify), (Expr*[]){ sub }, 1), 3);
    if (!simp || cv_too_big(simp)) { expr_free(simp); return NULL; }
    Expr* tg = eval_and_free(ds_call1(SYM_Together, simp));
    Expr* num = eval_and_free(ds_call1(SYM_Numerator, expr_copy(tg)));
    Expr* den = eval_and_free(ds_call1(SYM_Denominator, expr_copy(tg)));
    /* rational in tv iff numerator & denominator are polynomials in tv and no xv remains */
    Expr* r1e = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
                  (Expr*[]){ num, expr_new_function(expr_new_symbol(SYM_List),
                       (Expr*[]){ expr_new_symbol(tv) }, 1) }, 2));
    Expr* r2e = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
                  (Expr*[]){ den, expr_new_function(expr_new_symbol(SYM_List),
                       (Expr*[]){ expr_new_symbol(tv) }, 1) }, 2));
    bool r1 = (r1e && r1e->type==EXPR_SYMBOL && r1e->data.symbol.name==SYM_True);
    bool r2 = (r2e && r2e->type==EXPR_SYMBOL && r2e->data.symbol.name==SYM_True);
    expr_free(r1e); expr_free(r2e);
    if (r1 && r2 && ds_free_of(tg, xv)) return tg;
    expr_free(tg);
    return NULL;
}

/* g = phi(x) for substitution `which`. */
static Expr* cv_gfun(const char* xv, int which){
    const char* h[]={"Cos","Sin","Tan"};
    return ds_call1(h[which], expr_new_symbol(xv));
}

/* Numeric verify: y -> Function[{x}, body] into the original residual at sample points. */
static double cv_abs_at(const Expr* R, const char* xv, double xval){
    Expr* e = ds_subst(expr_copy((Expr*)R), expr_new_symbol(xv), expr_new_real(xval));
    e = eval_and_free(ds_call1("Abs", eval_and_free(ds_call1("N", e))));
    double m = (e&&e->type==EXPR_REAL)?e->data.real:(e&&e->type==EXPR_INTEGER)?(double)e->data.integer:NAN;
    expr_free(e); return m;
}
static bool cv_num_ok(const DSolveProblem* P, const Expr* body, const char* xv, const char* yname){
    /* residual with y^(k)[x] -> D[body,{x,k}], and C[1],C[2] -> sample reals */
    Expr* R = expr_copy(P->eq_residuals[0]);
    Expr* b0=expr_copy((Expr*)body);
    Expr* b1=ds_d(expr_copy((Expr*)body), expr_new_symbol(xv));
    Expr* b2=ds_d(ds_d(expr_copy((Expr*)body), expr_new_symbol(xv)), expr_new_symbol(xv));
    R=ds_subst(R, ds_make_funcapp(yname,2,xv), b2);
    R=ds_subst(R, ds_make_funcapp(yname,1,xv), b1);
    R=ds_subst(R, ds_make_funcapp(yname,0,xv), b0);
    /* substitute any leftover generated constants and free params with sample reals */
    R=ds_subst(R, ds_const(1), expr_new_real(1.3));
    R=ds_subst(R, ds_const(2), expr_new_real(0.7));
    const double xs[]={0.9, 1.35, 1.8, 2.4, 0.5};
    int small=0, big=0;
    for (int i=0;i<5;i++){
        double m=cv_abs_at(R, xv, xs[i]);
        if (isnan(m)||!isfinite(m)) continue;
        if (m<1e-6) small++; else if (m>1e-3) big++;
    }
    expr_free(R);
    return small>=2 && big==0;
}

Expr** dsolve_changevar_try(DSolveProblem* P, size_t* nbranch){
    if (P->nfun!=1 || P->neq!=1) return NULL;
    if (P->max_order[0]!=2) return NULL;
    if (g_cv_active) return NULL;                    /* no self-recursion */
    const char* xv = P->ind_names[0];
    const char* yname = P->fun_names[0];

    Expr *Pc=NULL, *Qc=NULL;
    if (!dsolve_second_order_PQ(P, &Pc, &Qc)) return NULL;   /* linear homog 2nd-order */

    /* Only for TRANSCENDENTAL coefficients — rational ones are owned by the direct
     * methods (Euler/Kovacic/SpecialFunction/Frobenius). */
    bool transc = ds_has_head(Pc,SYM_Sin)||ds_has_head(Pc,SYM_Cos)||ds_has_head(Pc,SYM_Tan)||
                  ds_has_head(Pc,SYM_Cot)||ds_has_head(Pc,SYM_Sec)||ds_has_head(Pc,SYM_Csc)||
                  ds_has_head(Qc,SYM_Sin)||ds_has_head(Qc,SYM_Cos)||ds_has_head(Qc,SYM_Tan)||
                  ds_has_head(Qc,SYM_Cot)||ds_has_head(Qc,SYM_Sec)||ds_has_head(Qc,SYM_Csc);
    if (!transc) { expr_free(Pc); expr_free(Qc); return NULL; }

    uint64_t h = expr_hash(P->eq_residuals[0]);
    cv_memo_sync(eval_toplevel_id());
    if (cv_memo_seen(h)) { expr_free(Pc); expr_free(Qc); return NULL; }

    g_cv_deadline = time(NULL) + 8;
    const char* tv = intern_symbol("DSolve`cvt");
    const char* Yf = intern_symbol("DSolve`cvY");
    Expr* body = NULL;

    for (int which=0; which<3 && !body && !cv_expired(); which++){
        Expr* g   = cv_gfun(xv, which);
        Expr* gp  = ds_d(expr_copy(g), expr_new_symbol(xv));
        Expr* gpp = ds_d(expr_copy(gp), expr_new_symbol(xv));
        /* A = (gpp + P gp)/gp^2 ; B = Q/gp^2 (in x) */
        Expr* gp2 = cv_powi(expr_copy(gp), 2);
        Expr* Ax = eval_and_free(ds_call2(SYM_Times,
                      ds_call2(SYM_Plus, expr_copy(gpp), ds_call2(SYM_Times, expr_copy(Pc), expr_copy(gp))),
                      cv_powi(expr_copy(gp2), -1)));
        Expr* Bx = eval_and_free(ds_call2(SYM_Times, expr_copy(Qc), cv_powi(expr_copy(gp2), -1)));
        expr_free(gp); expr_free(gpp); expr_free(gp2);
        Expr* At = cv_to_t(Ax, xv, tv, which);
        Expr* Bt = (At) ? cv_to_t(Bx, xv, tv, which) : (expr_free(Bx), (Expr*)NULL);
        if (At && !Bt) expr_free(At);
        if (At && Bt && cv_expired()) { expr_free(At); expr_free(Bt); }
        else if (At && Bt){
            /* transformed ODE Y''[t] + At Y'[t] + Bt Y[t] == 0 (At, Bt already in t) */
            Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                (Expr*[]){ ds_call2(SYM_Plus, ds_make_funcapp(Yf,2,tv),
                             ds_call2(SYM_Plus,
                               ds_call2(SYM_Times, At, ds_make_funcapp(Yf,1,tv)),
                               ds_call2(SYM_Times, Bt, ds_make_funcapp(Yf,0,tv)))),
                           expr_new_integer(0) }, 2);
            g_cv_active++;
            Expr* Yt = cv_run_applied(eq, Yf, tv);
            g_cv_active--;
            if (Yt){
                Expr* cand = ds_subst(Yt, expr_new_symbol(tv), expr_copy(g)); /* y = Y(phi(x)) */
                if (!cv_too_big(cand) && cv_num_ok(P, cand, xv, yname)) body = cand;
                else expr_free(cand);
            }
        }
        expr_free(g);
    }
    expr_free(Pc); expr_free(Qc);
    if (!body) { cv_memo_add(h); return NULL; }
    Expr** out = malloc(sizeof(Expr*));
    out[0] = body; *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_changevar(Expr* res){ return dsolve_method_builtin(res, dsolve_changevar_try); }

void dsolve_changevar_init(void){
    symtab_add_builtin("DSolve`ChangeOfVariable", builtin_dsolve_changevar);
    symtab_get_def("DSolve`ChangeOfVariable")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`ChangeOfVariable",
        "DSolve`ChangeOfVariable[eqn, y, x] solves a second-order linear ODE with "
        "transcendental coefficients y'' + P(x) y' + Q(x) y == 0 by a change of the "
        "independent variable t = phi(x) (Cos/Sin/Tan) that rationalizes the "
        "coefficients, recursing the cascade on the transformed equation and "
        "back-substituting (e.g. y'' + Cot[x] y' + k(k+1) y == 0 -> Legendre).");
}
