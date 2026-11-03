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
 * Candidate substitutions: t = Cos[x], Sin[x], Tan[x] for trig coefficients of the
 * bare x; and, since v0.6148 (see the generalisation block below), Cos/Sin/Tan/Sin^2
 * of s = k x for scaled trig arguments and t = E^(k x), E^(2 k x) for exponential
 * and hyperbolic coefficients.
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
#include "dsolve.h"              /* dsolve_cv_active */
#include "../core.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include "../parse.h"           /* parse_expression: the Sin/Cos rewrite rule */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

static time_t g_cv_deadline;
static bool cv_expired(void) { return time(NULL) >= g_cv_deadline; }
int dsolve_cv_active = 0;       /* re-entry guard; also read by dsolve_specialform.c (dsolve.h) */
#define g_cv_active dsolve_cv_active

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
    } else if (which==3){ /* t=Sin[x]^2 : sin=Sqrt[t], cos=Sqrt[1-t] (v0.6148) */
        Expr* omt = ds_call2(SYM_Subtract, expr_new_integer(1), expr_copy(t));
        s=eval_and_free(ds_call1("Sqrt", expr_copy(t))); c=eval_and_free(ds_call1("Sqrt", omt));
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
    if (which==3) return eval_and_free(cv_powi(ds_call1(SYM_Sin, expr_new_symbol(xv)), 2));
    return ds_call1(h[which], expr_new_symbol(xv));
}

/* Numeric verify: y -> Function[{x}, body] into the original residual at sample points. */
static double cv_abs_at(const Expr* R, const char* xv, double xval){
    Expr* e = ds_subst(expr_copy((Expr*)R), expr_new_symbol(xv), expr_new_real(xval));
    e = eval_and_free(ds_call1("Abs", eval_and_free(ds_call1("N", e))));
    double m = (e&&e->type==EXPR_REAL)?e->data.real:(e&&e->type==EXPR_INTEGER)?(double)e->data.integer:NAN;
    expr_free(e); return m;
}
/* Collect distinct interned symbol names appearing in ARGUMENT position (i.e. the
 * free variables/parameters) — a plain-symbol function head (LegendreP, Cot, Plus,
 * …) is NOT a parameter and must never be substituted, or the expression is
 * destroyed. */
static void cv_collect_params(const Expr* e, const char** names, int* n, int cap){
    if (!e || *n >= cap) return;
    if (e->type == EXPR_SYMBOL){
        const char* nm = e->data.symbol.name;
        for (int i = 0; i < *n; i++) if (names[i] == nm) return;
        names[(*n)++] = nm; return;
    }
    if (e->type == EXPR_FUNCTION){
        /* recurse into the head only if it is itself compound (e.g.
         * Derivative[2][y]); a plain-symbol head is the function name, not a param */
        if (e->data.function.head && e->data.function.head->type != EXPR_SYMBOL)
            cv_collect_params(e->data.function.head, names, n, cap);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            cv_collect_params(e->data.function.args[i], names, n, cap);
    }
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
    /* substitute any leftover generated constants with sample reals */
    R=ds_subst(R, ds_const(1), expr_new_real(1.3));
    R=ds_subst(R, ds_const(2), expr_new_real(0.7));
    /* Instantiate every remaining FREE PARAMETER (symbolic degree k, coefficients
     * a,b,n,...) with a distinct sample real — otherwise a solution carrying a
     * symbolic parameter (e.g. LegendreP[k, Cos[x]] for symbolic k) never
     * numericizes and the whole check spuriously rejects a correct transform. */
    {
        const char* skip[] = { xv, intern_symbol("E"), intern_symbol("Pi"),
            intern_symbol("I"), intern_symbol("C"), intern_symbol("EulerGamma"),
            intern_symbol("Degree"), intern_symbol("GoldenRatio"),
            intern_symbol("Catalan"), intern_symbol("Infinity") };
        const int nskip = (int)(sizeof(skip)/sizeof(skip[0]));
        const char* syms[64]; int ns = 0;
        cv_collect_params(R, syms, &ns, 64);
        int pi = 0;
        for (int i = 0; i < ns; i++){
            bool sk = false;
            for (int j = 0; j < nskip; j++) if (syms[i] == skip[j]) { sk = true; break; }
            if (sk) continue;
            /* distinct generic reals, kept small so Legendre/Bessel args stay in range */
            double v = 0.31 + 0.17 * (double)pi; pi++;
            R = ds_subst(R, expr_new_symbol(syms[i]), expr_new_real(v));
        }
    }
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

/* Transformed ODE Y''[t] + At Y'[t] + Bt Y[t] == 0 (At, Bt consumed, already in t):
 * solve it with the cascade (self-recursion blocked), compose with t = g and verify
 * numerically on the ORIGINAL equation.  `s_of_x` (borrowed, may be NULL) maps the
 * working variable `wv` of g back to the original x (wv = k x); the candidate is
 * returned in x, or NULL. */
static Expr* cv_solve_transformed(const DSolveProblem* P, Expr* At, Expr* Bt,
                                  const Expr* g, const char* wv, const Expr* s_of_x,
                                  const char* tv, bool clear){
    const char* xv = P->ind_names[0];
    const char* Yf = intern_symbol("DSolve`cvY");
    /* `clear`: hand the cascade the denominator-cleared polynomial form
     * m Y'' + (m A) Y' + (m B) Y == 0, m = PolynomialLCM of the denominators.  A
     * nested DSolve skips dsolve_linear_normalize (depth > 1 is assumed canonical),
     * so the uncleared Y'' + a/t^2 Y == 0 never reaches the Euler form t^2 Y'' + a Y. */
    Expr* m = expr_new_integer(1);
    if (clear) {
        m = eval_and_free(ds_call2(SYM_PolynomialLCM,
                eval_and_free(ds_call1(SYM_Denominator, expr_copy(At))),
                eval_and_free(ds_call1(SYM_Denominator, expr_copy(Bt)))));
        At = eval_and_free(ds_call1(SYM_Cancel, ds_call2(SYM_Times, At, expr_copy(m))));
        Bt = eval_and_free(ds_call1(SYM_Cancel, ds_call2(SYM_Times, Bt, expr_copy(m))));
    }
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
        (Expr*[]){ ds_call2(SYM_Plus, ds_call2(SYM_Times, m, ds_make_funcapp(Yf,2,tv)),
                     ds_call2(SYM_Plus,
                       ds_call2(SYM_Times, At, ds_make_funcapp(Yf,1,tv)),
                       ds_call2(SYM_Times, Bt, ds_make_funcapp(Yf,0,tv)))),
                   expr_new_integer(0) }, 2);
    g_cv_active++;
    Expr* Yt = cv_run_applied(eq, Yf, tv);
    g_cv_active--;
    if (!Yt) return NULL;
    Expr* cand = ds_subst(Yt, expr_new_symbol(tv), expr_copy((Expr*)g)); /* y = Y(phi) */
    if (s_of_x) cand = eval_and_free(ds_subst(cand, expr_new_symbol(wv), expr_copy((Expr*)s_of_x)));
    if (!cv_too_big(cand) && cv_num_ok(P, cand, xv, P->fun_names[0])) return cand;
    expr_free(cand);
    return NULL;
}

/* A = (g'' + P g')/g'^2 and B = Q/g'^2 in the working variable wv (Pc, Qc borrowed). */
static void cv_AB(const Expr* g, const Expr* Pc, const Expr* Qc, const char* wv,
                  Expr** Ax, Expr** Bx){
    Expr* gp  = ds_d(expr_copy((Expr*)g), expr_new_symbol(wv));
    Expr* gpp = ds_d(expr_copy(gp), expr_new_symbol(wv));
    Expr* gp2 = cv_powi(expr_copy(gp), 2);
    *Ax = eval_and_free(ds_call2(SYM_Times,
              ds_call2(SYM_Plus, gpp, ds_call2(SYM_Times, expr_copy((Expr*)Pc), gp)),
              cv_powi(expr_copy(gp2), -1)));
    *Bx = eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)Qc), cv_powi(gp2, -1)));
}

/* The trigonometric family t = Cos/Sin/Tan[wv] (and, when `clear`, Sin[wv]^2 -- a
 * coefficient set even in Sin and Cos puts the singular points at 0, 1, oo, the
 * Gauss form) on Y'' + Pc Y' + Qc Y == 0 in wv. */
static Expr* cv_try_trig(const DSolveProblem* P, const Expr* Pc, const Expr* Qc,
                         const char* wv, const Expr* s_of_x, const char* tv, bool clear){
    Expr* body = NULL;
    int nwhich = clear ? 4 : 3;              /* the HEAD bare-x pass keeps its three */
    for (int which=0; which<nwhich && !body && !cv_expired(); which++){
        Expr* g = cv_gfun(wv, which);
        Expr *Ax, *Bx;
        cv_AB(g, Pc, Qc, wv, &Ax, &Bx);
        Expr* At = cv_to_t(Ax, wv, tv, which);
        Expr* Bt = (At) ? cv_to_t(Bx, wv, tv, which) : (expr_free(Bx), (Expr*)NULL);
        if (At && !Bt) expr_free(At);
        if (At && Bt && cv_expired()) { expr_free(At); expr_free(Bt); }
        else if (At && Bt) body = cv_solve_transformed(P, At, Bt, g, wv, s_of_x, tv, clear);
        expr_free(g);
    }
    return body;
}

/* ---- Generalisation (v0.6148): scaled trig arguments, hyperbolic and exponential
 * coefficients.  The original family above only rationalises trig functions of the
 * bare variable x.  Coefficients such as Cot[2 x], Tan[x/2]^2, Cot[b x] (one base
 * frequency k with every argument an integer multiple of k x), and exponentials
 * E^(c x) / hyperbolic functions (rewritten to exponentials) with every c an
 * integer multiple of one k, are brought to that shape:
 *   trig:  s = k x, the ODE in s is Y'' + (P/k) Y' + (Q/k^2) Y == 0 with every
 *          trig function written in Sin/Cos and TrigExpand'ed to Sin[s], Cos[s];
 *   exp:   t = E^(k x), x -> Log[t]/k, PowerExpand (t > 0) leaves powers of t.
 * Both run only after the original family declines, so every equation it closed is
 * unchanged. */

/* A rational NUMBER (Integer / BigInt / Rational). */
static bool cv_is_ratnum(const Expr* e){
    if (!e) return false;
    if (e->type==EXPR_INTEGER || e->type==EXPR_BIGINT) return true;
    return e->type==EXPR_FUNCTION && head_is(e, SYM_Rational);
}

/* Collect (borrowed) the x-dependent arguments of the trig heads (trig != 0) or the
 * exponents of E^(...) (trig == 0) in e. Returns false on overflow. */
static bool cv_collect_args(const Expr* e, const char* xv, bool trig,
                            const Expr** out, int* n, int cap){
    if (!e || e->type != EXPR_FUNCTION) return true;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && e->data.function.arg_count >= 1){
        const char* hn = h->data.symbol.name;
        const Expr* a = NULL;
        if (trig && e->data.function.arg_count == 1 &&
            (hn==SYM_Sin||hn==SYM_Cos||hn==SYM_Tan||hn==SYM_Cot||hn==SYM_Sec||hn==SYM_Csc))
            a = e->data.function.args[0];
        else if (!trig && hn==SYM_Power && e->data.function.arg_count==2 &&
                 e->data.function.args[0]->type==EXPR_SYMBOL &&
                 e->data.function.args[0]->data.symbol.name==intern_symbol("E"))
            a = e->data.function.args[1];
        if (a && !ds_free_of(a, xv)){
            if (*n >= cap) return false;
            out[(*n)++] = a;
            return true;            /* the argument itself is checked for linearity */
        }
    }
    if (!cv_collect_args(h, xv, trig, out, n, cap)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!cv_collect_args(e->data.function.args[i], xv, trig, out, n, cap)) return false;
    return true;
}

/* Base frequency k: every argument must be c_i x with c_i free of x and c_i/c_0 a
 * rational number; k = c_0 GCD[c_i/c_0] (so every c_i/k is an integer).  NULL if
 * not of that shape. */
static Expr* cv_base_freq(const Expr* const* args, int n, const char* xv){
    if (n <= 0) return NULL;
    Expr* c0 = NULL;
    Expr** r = malloc((size_t)n * sizeof(Expr*));
    int nr = 0;
    bool ok = true;
    for (int i = 0; i < n && ok; i++){
        Expr* c = eval_and_free(ds_call1(SYM_Together,
                      ds_call2(SYM_Times, expr_copy((Expr*)args[i]),
                               cv_powi(expr_new_symbol(xv), -1))));
        if (!ds_free_of(c, xv) || (c->type==EXPR_INTEGER && c->data.integer==0)){
            expr_free(c); ok = false; break;
        }
        if (!c0) { c0 = c; r[nr++] = expr_new_integer(1); continue; }
        Expr* q = eval_and_free(ds_call1(SYM_Together,
                      ds_call2(SYM_Times, c, cv_powi(expr_copy(c0), -1))));
        if (!cv_is_ratnum(q)) { expr_free(q); ok = false; break; }
        r[nr++] = q;
    }
    Expr* k = NULL;
    if (ok){
        Expr* gcd = eval_and_free(expr_new_function(expr_new_symbol(SYM_GCD), r, (size_t)nr));
        if (cv_is_ratnum(gcd))
            k = eval_and_free(ds_call2(SYM_Times, expr_copy(c0), gcd));
        else expr_free(gcd);
    } else {
        for (int i = 0; i < nr; i++) expr_free(r[i]);
    }
    free(r);
    expr_free(c0);
    return k;
}

/* Every trig head written in Sin/Cos, then TrigExpand'ed. e consumed. */
static Expr* cv_trig_expand(Expr* e){
    static const char* const src =
        "Function[{e}, TrigExpand[e /. {Tan[DSolve`cvz_] :> Sin[DSolve`cvz]/Cos[DSolve`cvz], "
        "Cot[DSolve`cvz_] :> Cos[DSolve`cvz]/Sin[DSolve`cvz], Sec[DSolve`cvz_] :> 1/Cos[DSolve`cvz], "
        "Csc[DSolve`cvz_] :> 1/Sin[DSolve`cvz]}]]";
    Expr* f = parse_expression(src);
    if (!f) return e;
    return cv_beval(expr_new_function(f, (Expr*[]){ e }, 1), 3);
}

/* Scaled-trig family: s = k x. */
static Expr* cv_try_trig_scaled(const DSolveProblem* P, const Expr* Pc, const Expr* Qc,
                                const char* tv){
    const char* xv = P->ind_names[0];
    const Expr* args[64]; int n = 0;
    if (!cv_collect_args(Pc, xv, true, args, &n, 64) ||
        !cv_collect_args(Qc, xv, true, args, &n, 64) || n == 0) return NULL;
    Expr* k = cv_base_freq(args, n, xv);
    if (!k) return NULL;
    const char* sv = intern_symbol("DSolve`cvs");
    /* x -> s/k ; Ps = P/k, Qs = Q/k^2 */
    Expr* x_of_s = eval_and_free(ds_call2(SYM_Times, expr_new_symbol(sv), cv_powi(expr_copy(k), -1)));
    Expr* Ps = eval_and_free(ds_call2(SYM_Times,
                   ds_subst(expr_copy((Expr*)Pc), expr_new_symbol(xv), expr_copy(x_of_s)),
                   cv_powi(expr_copy(k), -1)));
    Expr* Qs = eval_and_free(ds_call2(SYM_Times,
                   ds_subst(expr_copy((Expr*)Qc), expr_new_symbol(xv), expr_copy(x_of_s)),
                   cv_powi(expr_copy(k), -2)));
    expr_free(x_of_s);
    Ps = cv_trig_expand(Ps);
    Qs = cv_trig_expand(Qs);
    Expr* s_of_x = eval_and_free(ds_call2(SYM_Times, k, expr_new_symbol(xv)));
    Expr* body = NULL;
    if (ds_free_of(Ps, xv) && ds_free_of(Qs, xv) && !cv_too_big(Ps) && !cv_too_big(Qs))
        body = cv_try_trig(P, Ps, Qs, sv, s_of_x, tv, true);
    expr_free(Ps); expr_free(Qs); expr_free(s_of_x);
    return body;
}

/* Rationalise e (in x) under t = E^(k x): x -> Log[t]/k, PowerExpand, Together;
 * returns the rational function of t or NULL.  e consumed. */
static Expr* cv_exp_to_t(Expr* e, const char* xv, const char* tv, const Expr* k){
    Expr* x_of_t = eval_and_free(ds_call2(SYM_Times, ds_call1(SYM_Log, expr_new_symbol(tv)),
                                          cv_powi(expr_copy((Expr*)k), -1)));
    e = ds_subst(e, expr_new_symbol(xv), x_of_t);
    e = cv_beval(ds_call1("PowerExpand", e), 3);
    if (!e || cv_too_big(e)) { expr_free(e); return NULL; }
    Expr* tg = eval_and_free(ds_call1(SYM_Together, e));
    if (!ds_free_of(tg, xv) || ds_has_head(tg, SYM_Log) || !ds_is_rational_in(tg, tv)){
        expr_free(tg); return NULL;
    }
    return tg;
}

/* Exponential family (hyperbolic coefficients first rewritten by TrigToExp):
 * t = E^(k x). */
static Expr* cv_try_exp(const DSolveProblem* P, const Expr* Pc0, const Expr* Qc0,
                        const char* tv){
    const char* xv = P->ind_names[0];
    Expr* Pc = eval_and_free(ds_call1("TrigToExp", expr_copy((Expr*)Pc0)));
    Expr* Qc = eval_and_free(ds_call1("TrigToExp", expr_copy((Expr*)Qc0)));
    const Expr* args[64]; int n = 0;
    Expr* body = NULL;
    if (cv_collect_args(Pc, xv, false, args, &n, 64) &&
        cv_collect_args(Qc, xv, false, args, &n, 64) && n > 0){
        Expr* k0 = cv_base_freq(args, n, xv);
        /* t = E^(k x), then t = E^(2 k x): a coefficient set even under
         * x -> -x (Tanh, Sech^2, Coth) has singular points at t = +-I under the
         * first and the Gauss points 0, -1, oo under the second. */
        for (int pass = 0; k0 && pass < 2 && !body && !cv_expired(); pass++){
            Expr* k = pass ? eval_and_free(ds_call2(SYM_Times, expr_new_integer(2), expr_copy(k0)))
                           : expr_copy(k0);
            Expr* g = eval_and_free(ds_call2(SYM_Power, expr_new_symbol(intern_symbol("E")),
                          ds_call2(SYM_Times, expr_copy(k), expr_new_symbol(xv))));
            Expr *Ax, *Bx;
            cv_AB(g, Pc, Qc, xv, &Ax, &Bx);
            Expr* At = cv_exp_to_t(Ax, xv, tv, k);
            Expr* Bt = At ? cv_exp_to_t(Bx, xv, tv, k) : (expr_free(Bx), (Expr*)NULL);
            if (At && !Bt) expr_free(At);
            if (At && Bt && !cv_expired()) body = cv_solve_transformed(P, At, Bt, g, xv, NULL, tv, true);
            else if (At && Bt) { expr_free(At); expr_free(Bt); }
            expr_free(g); expr_free(k);
        }
        expr_free(k0);
    }
    expr_free(Pc); expr_free(Qc);
    if (body) {
        /* Compose t^p as E^(p k x) (Sqrt[E^(2 x)] -> E^x): on each branch the two
         * differ by a constant factor, so this is still a basis -- re-verified. */
        static const char* const src =
            "Function[{e}, e /. Power[Power[E, DSolve`cvq_], DSolve`cvp_] :> "
            "Power[E, Expand[DSolve`cvp DSolve`cvq]]]";
        Expr* f = parse_expression(src);
        if (f) {
            Expr* c2 = cv_beval(expr_new_function(f, (Expr*[]){ expr_copy(body) }, 1), 3);
            if (c2 && !cv_too_big(c2) && cv_num_ok(P, c2, xv, P->fun_names[0])) {
                expr_free(body); body = c2;
            } else expr_free(c2);
        }
    }
    return body;
}

static bool cv_has_trig(const Expr* e){
    return ds_has_head(e,SYM_Sin)||ds_has_head(e,SYM_Cos)||ds_has_head(e,SYM_Tan)||
           ds_has_head(e,SYM_Cot)||ds_has_head(e,SYM_Sec)||ds_has_head(e,SYM_Csc);
}
static bool cv_has_hyp_or_exp(const Expr* e, const char* xv){
    if (ds_has_head(e,SYM_Sinh)||ds_has_head(e,SYM_Cosh)||ds_has_head(e,SYM_Tanh)||
        ds_has_head(e,SYM_Coth)||ds_has_head(e,SYM_Sech)||ds_has_head(e,SYM_Csch)) return true;
    const Expr* args[1]; int n = 0;
    cv_collect_args(e, xv, false, args, &n, 1);   /* stops at the first E^(..x..) */
    return n > 0;
}

Expr** dsolve_changevar_try(DSolveProblem* P, size_t* nbranch){
    if (P->nfun!=1 || P->neq!=1) return NULL;
    if (P->max_order[0]!=2) return NULL;
    if (g_cv_active) return NULL;                    /* no self-recursion */
    const char* xv = P->ind_names[0];

    Expr *Pc=NULL, *Qc=NULL;
    if (!dsolve_second_order_PQ(P, &Pc, &Qc)) return NULL;   /* linear homog 2nd-order */

    /* Only for TRANSCENDENTAL coefficients — rational ones are owned by the direct
     * methods (Euler/Kovacic/SpecialFunction/Frobenius). */
    bool trig = cv_has_trig(Pc) || cv_has_trig(Qc);
    bool expo = !trig && (cv_has_hyp_or_exp(Pc, xv) || cv_has_hyp_or_exp(Qc, xv));
    if (!trig && !expo) { expr_free(Pc); expr_free(Qc); return NULL; }

    uint64_t h = expr_hash(P->eq_residuals[0]);
    cv_memo_sync(eval_toplevel_id());
    if (cv_memo_seen(h)) { expr_free(Pc); expr_free(Qc); return NULL; }

    g_cv_deadline = time(NULL) + 8;
    const char* tv = intern_symbol("DSolve`cvt");
    Expr* body = NULL;
    if (trig){
        body = cv_try_trig(P, Pc, Qc, xv, NULL, tv, false);    /* t = phi(x), bare x (HEAD form) */
        if (!body && !cv_expired()) body = cv_try_trig_scaled(P, Pc, Qc, tv);
    } else {
        body = cv_try_exp(P, Pc, Qc, tv);
    }
    expr_free(Pc); expr_free(Qc);
    if (!body) { cv_memo_add(h); return NULL; }
    Expr** out = malloc(sizeof(Expr*));
    out[0] = body; *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_changevar(Expr* res){ return dsolve_method_builtin(res, dsolve_changevar_try); }

void dsolve_changevar_init(void){
    tc_register_unwind_counter(&g_cv_active);   /* restored on a timeout unwind */
    symtab_add_builtin("DSolve`ChangeOfVariable", builtin_dsolve_changevar);
    symtab_get_def("DSolve`ChangeOfVariable")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`ChangeOfVariable",
        "DSolve`ChangeOfVariable[eqn, y, x] solves a second-order linear ODE with "
        "transcendental coefficients y'' + P(x) y' + Q(x) y == 0 by a change of the "
        "independent variable t = phi(x) that rationalizes the coefficients: "
        "Cos/Sin/Tan/Sin^2 of k x for trig coefficients whose arguments are integer "
        "multiples of one k x, and E^(k x) or E^(2 k x) for exponential/hyperbolic "
        "ones.  Recurses the cascade on the transformed equation and back-substitutes "
        "(y'' + Cot[x] y' + k(k+1) y == 0 -> Legendre; y'' + c E^(k x) y == 0 -> Bessel).");
}
