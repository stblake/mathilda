/*
 * dsolve_scaling2.c — DSolve`ScalingReduction / DSolve`ScalingReductionImplicit.
 *
 * Second-order ODEs F(x, y, y', y'') == 0 invariant under a SCALING or a TRANSLATION, the
 * second-order counterparts of the first-order homogeneous / class-G methods.  Three
 * symmetries, tried in this order:
 *
 *  1. y -> l y (x fixed): F(x, l Y, l P, l Q) == l^w F(x, Y, P, Q).  The Riccati
 *     substitution u = y'/y (y' = u y, y'' = (u' + u^2) y) gives
 *
 *          F(x, 1, u, u' + u^2) == 0,
 *
 *     a FIRST-order ODE in u(x) solved by the cascade; then y = C[2] E^Integrate[u].
 *     (y y'' == g y^2 + f y y' + y'^2 -> u' == g + f u, linear.)
 *
 *  2. x -> l x, y -> l^k y (isobaric, "equidimensional in x"): with Y = x^k U,
 *     P = x^(k-1) (k U + V), Q = x^(k-2) (W + (2k-1) V + k(k-1) U) and t = Log[x],
 *     F becomes x^w F~(u, u'(t), u''(t)) -- AUTONOMOUS in u(t), which the cascade's
 *     AutonomousReduction owns (often as an inert first integral when the quadrature
 *     is elliptic, as for the Lane-Emden index-5 equation x y'' + 2 y' + x y^5 == 0,
 *     k = -1/2, u'' == u/4 - u^5).  Mapped back by u -> y/x^k, t -> Log[x].
 *
 *  3. translation x -> x + e, y -> y - rho e (R_x == rho R_Y, rho constant): v =
 *     y + rho x makes F autonomous in v(x) -- (x - y) y'' + 2 y' (y' + 1) == 0,
 *     (a x + b y) y'^3 + y'' == 0.
 *
 * Linear equations are left alone (a homogeneous linear one is always y-scale
 * invariant, but its Riccati is no easier than the linear methods that declined it).
 *
 * Weight detection.  The Euler identity of the invariance at l = 1 is linear in the
 * weights: (1) Y F_Y + P F_P + Q F_Q == w F is checked directly (w must come out free
 * of x, Y, P, Q); (2) x F_x - P F_P - 2 Q F_Q + k (Y F_Y + P F_P + Q F_Q) == w F is
 * solved for (k, w) from exact rational samples and then CHECKED identically (Together
 * + zero test) -- a sampling artefact can only decline.  (2) is gated to equations
 * ALGEBRAIC in x, Y, P, Q, which is the only way an isobaric equation can be built.
 *
 * Results.  Explicit u-bodies give explicit y-bodies, verified by back-substitution
 * (dsolve_run).  Otherwise (2) maps each implicit u-relation G(t, u[t]) == C[1] back;
 * an inert indefinite integral over u[t] is first made definite (lower limit a
 * generic point), since u[t] -> y[x] x^-k would rewrite its variable slot.  Each
 * mapped relation is verified NUMERICALLY against the original second-order equation
 * by implicit differentiation (y' = -G_x/G_y and its total derivative), with the
 * constants and every parameter sampled -- the shared implicit verifier only checks
 * y' and so is vacuous at order two.  The verified relations are stashed for the
 * implicit pass.
 *
 * Runs right after Liouville, ahead of the integrating-factor and point-symmetry
 * searches (heuristic, and slow where this is a direct reduction), top level only,
 * re-entry guarded (the sub-solves must not rescale), per-top-level decline memo;
 * MATHILDA_NO_SCALING2 disables both passes for A/B.
 */
#include "dsolve_common.h"
#include "../core.h"
#include "dsolve.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include "../message.h"
#include "integrate.h"          /* g_integrate_quiet */
#include "../parse.h"           /* parse_expression: the symbolic-verify abstraction */
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#define SC_STASH_TAG    7
#define SC_SUBSOLVE_SEC 4
#define SC_QUAD_SEC     2

static int g_sc_active = 0;        /* re-entry guard: no rescaling inside a sub-solve */

#define SC_MEMO_SLOTS 32
static uint64_t sc_epoch = 0;
static int sc_memo_n = 0;
static uint64_t sc_memo[SC_MEMO_SLOTS];
static void sc_memo_sync(uint64_t tid){ if(tid!=sc_epoch){sc_epoch=tid;sc_memo_n=0;} }
static bool sc_memo_seen(uint64_t h){ for(int i=0;i<sc_memo_n;i++) if(sc_memo[i]==h) return true; return false; }
static void sc_memo_add(uint64_t h){ if(sc_memo_n<SC_MEMO_SLOTS && !sc_memo_seen(h)) sc_memo[sc_memo_n++]=h; }

static Expr* sc_tc(Expr* call, int secs) {
    return expr_new_function(expr_new_symbol(intern_symbol("TimeConstrained")),
               (Expr*[]){ call, expr_new_integer(secs),
                          expr_new_symbol(intern_symbol("$Aborted")) }, 3);
}

/* Free of all four jet variables (robust, D-based). */
static bool sc_free4(const Expr* e, const char* x, const char* Y, const char* P, const char* Q) {
    return ds_free_of(e, x) && ds_free_of(e, Y) && ds_free_of(e, P) && ds_free_of(e, Q);
}

/* Algebraic in the jet variables: every subexpression involving one is a Plus, Times
 * or Power with an exponent free of them (an isobaric equation can only be built so). */
static bool sc_algebraic(const Expr* e, const char* x, const char* Y, const char* P, const char* Q) {
    if (!e || e->type != EXPR_FUNCTION) return true;
    if (!ds_contains(e, x) && !ds_contains(e, Y) && !ds_contains(e, P) && !ds_contains(e, Q))
        return true;
    if (head_is(e, SYM_Plus) || head_is(e, SYM_Times)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!sc_algebraic(e->data.function.args[i], x, Y, P, Q)) return false;
        return true;
    }
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* ex = e->data.function.args[1];
        if (ds_contains(ex, x) || ds_contains(ex, Y) || ds_contains(ex, P) || ds_contains(ex, Q))
            return false;
        return sc_algebraic(e->data.function.args[0], x, Y, P, Q);
    }
    return false;
}

/* The residual as R(x, Y, P, Q), or NULL if y still occurs. */
static Expr* sc_residual(const DSolveProblem* P, const char* Yn, const char* Pn, const char* Qn) {
    const char* yname = P->fun_names[0];
    const char* xv    = P->ind_names[0];
    Expr* R = expr_copy(P->eq_residuals[0]);
    R = ds_subst(R, ds_make_funcapp(yname, 2, xv), expr_new_symbol(Qn));
    R = ds_subst(R, ds_make_funcapp(yname, 1, xv), expr_new_symbol(Pn));
    R = ds_subst(R, ds_make_funcapp(yname, 0, xv), expr_new_symbol(Yn));
    if (ds_contains(R, yname)) { expr_free(R); return NULL; }
    return R;
}

/* C[k] <-> DSolve`scK<k>, k = 1..4 (Integrate treats a symbol, not C[k], as a
 * parameter).  `e` consumed. */
static Expr* sc_sym_consts(Expr* e, bool forward) {
    for (int k = 1; k <= 4; k++) {
        char nm[32];
        snprintf(nm, sizeof nm, "DSolve`scK%d", k);
        Expr* s = expr_new_symbol(intern_symbol(nm));
        e = forward ? ds_subst(e, ds_const(k), s) : ds_subst(e, s, ds_const(k));
    }
    return eval_and_free(e);
}

/* Integrate[u, x] when it closes in budget (constants symbolized), else the inert
 * Inactive[Integrate][u, x].  `u` borrowed. */
static Expr* sc_quadrature(const Expr* u, const char* xv) {
    Expr* us = sc_sym_consts(expr_copy((Expr*)u), true);
    g_integrate_quiet++;
    Expr* F = eval_and_free(sc_tc(ds_call2(SYM_Integrate, us, expr_new_symbol(xv)), SC_QUAD_SEC));
    g_integrate_quiet--;
    if (F && !ds_has_head(F, SYM_Integrate) && !ds_contains(F, intern_symbol("$Aborted")))
        return sc_sym_consts(F, false);
    if (F) expr_free(F);
    Expr* inact = expr_new_function(expr_new_symbol(SYM_Inactive),
                      (Expr*[]){ expr_new_symbol(SYM_Integrate) }, 1);
    return expr_new_function(inact, (Expr*[]){ expr_copy((Expr*)u), expr_new_symbol(xv) }, 2);
}

/* The sub-solve DSolve[ode, f[v], v] under the re-entry guard, messages muted.
 * `ode` consumed. */
static Expr* sc_subsolve(Expr* ode, const char* fname, const char* v, int secs) {
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ ode, ds_make_funcapp(fname, 0, v), expr_new_symbol(v) }, 3);
    g_sc_active++;
    mth_msg_suppress_push();
    Expr* r = eval_and_free(sc_tc(call, secs));
    mth_msg_suppress_pop();
    g_sc_active--;
    return r;
}

/* ---- (1) y -> l y : u = y'/y ---- */
static Expr** sc_try_yscale(const DSolveProblem* P, const Expr* R, const char* Yn,
                            const char* Pn, const char* Qn, size_t* nbranch) {
    const char* xv = P->ind_names[0];
    Expr* num = eval_and_free(ds_call2(SYM_Plus, ds_call2(SYM_Plus,
                    ds_call2(SYM_Times, expr_new_symbol(Yn), ds_d(expr_copy((Expr*)R), expr_new_symbol(Yn))),
                    ds_call2(SYM_Times, expr_new_symbol(Pn), ds_d(expr_copy((Expr*)R), expr_new_symbol(Pn)))),
                    ds_call2(SYM_Times, expr_new_symbol(Qn), ds_d(expr_copy((Expr*)R), expr_new_symbol(Qn)))));
    Expr* w = eval_and_free(ds_call1(SYM_Together, ds_call2(SYM_Times, num,
                  ds_call2(SYM_Power, expr_copy((Expr*)R), expr_new_integer(-1)))));
    bool homog = sc_free4(w, xv, Yn, Pn, Qn)
              && !ds_contains(w, intern_symbol("ComplexInfinity"))
              && !ds_contains(w, intern_symbol("Indeterminate"));
    expr_free(w);
    if (!homog) return NULL;

    const char* un = intern_symbol("DSolve`scu");
    Expr* Ru = expr_copy((Expr*)R);
    Ru = ds_subst(Ru, expr_new_symbol(Qn), ds_call2(SYM_Plus, ds_make_funcapp(un, 1, xv),
                      ds_call2(SYM_Power, ds_make_funcapp(un, 0, xv), expr_new_integer(2))));
    Ru = ds_subst(Ru, expr_new_symbol(Pn), ds_make_funcapp(un, 0, xv));
    Ru = ds_subst(Ru, expr_new_symbol(Yn), expr_new_integer(1));
    Ru = eval_and_free(ds_call1(SYM_Numerator, ds_call1(SYM_Together, Ru)));
    if (!ds_contains(Ru, un)) { expr_free(Ru); return NULL; }
    Expr* ode = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ Ru, expr_new_integer(0) }, 2);
    Expr* r = sc_subsolve(ode, un, xv, SC_SUBSOLVE_SEC);

    size_t nb = 0;
    Expr** ub = dsolve_extract_applied_bodies(r, un, &nb);
    expr_free(r);
    if (!ub) return NULL;
    Expr** out = malloc(nb * sizeof(Expr*));
    size_t no = 0;
    for (size_t i = 0; i < nb; i++) {
        /* a u-body without the sub-solve's constant is a singular solution of the
         * Riccati equation: it would give y short of two constants */
        if (!ds_contains(ub[i], intern_symbol("C"))) { expr_free(ub[i]); continue; }
        Expr* I = sc_quadrature(ub[i], xv);
        expr_free(ub[i]);
        out[no++] = eval_and_free(ds_call2(SYM_Times, ds_const(2), ds_call1("Exp", I)));
    }
    free(ub);
    if (no == 0) { free(out); return NULL; }
    *nbranch = no;
    return out;
}

/* ---- (2) isobaric : y = x^k u(Log[x]) ---- */

static Expr* sc_q(int64_t a, int64_t b) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Rational),
               (Expr*[]){ expr_new_integer(a), expr_new_integer(b) }, 2));
}

static Expr* sc_at(const Expr* e, const char* const* vars, const int64_t* pt) {
    Expr* r = expr_copy((Expr*)e);
    for (int i = 0; i < 4; i++)
        r = ds_subst(r, expr_new_symbol(vars[i]), sc_q(pt[2 * i], pt[2 * i + 1]));
    return eval_and_free(ds_call1(SYM_Together, r));
}

static bool sc_finite(const Expr* e) {
    return !ds_contains(e, intern_symbol("ComplexInfinity"))
        && !ds_contains(e, intern_symbol("Indeterminate"))
        && !ds_contains(e, intern_symbol("DirectedInfinity"));
}

/* Solve A + k B - w R == 0 for (k, w) from samples, then check it identically. */
static bool sc_weights(const Expr* R, const char* x, const char* Y, const char* P,
                       const char* Q, Expr** kout, Expr** wout) {
    Expr* Rx = ds_d(expr_copy((Expr*)R), expr_new_symbol(x));
    Expr* RY = ds_d(expr_copy((Expr*)R), expr_new_symbol(Y));
    Expr* RP = ds_d(expr_copy((Expr*)R), expr_new_symbol(P));
    Expr* RQ = ds_d(expr_copy((Expr*)R), expr_new_symbol(Q));
    Expr* PRP = eval_and_free(ds_call2(SYM_Times, expr_new_symbol(P), RP));
    Expr* QRQ = eval_and_free(ds_call2(SYM_Times, expr_new_symbol(Q), RQ));
    /* A = x R_x - P R_P - 2 Q R_Q ;  B = Y R_Y + P R_P + Q R_Q */
    Expr* A = eval_and_free(ds_call2(SYM_Subtract,
                  ds_call2(SYM_Times, expr_new_symbol(x), Rx),
                  ds_call2(SYM_Plus, expr_copy(PRP), ds_call2(SYM_Times, expr_new_integer(2), expr_copy(QRQ)))));
    Expr* B = eval_and_free(ds_call2(SYM_Plus, ds_call2(SYM_Times, expr_new_symbol(Y), RY),
                  ds_call2(SYM_Plus, PRP, QRQ)));
    const char* vars[4] = { x, Y, P, Q };
    static const int64_t pts[2][3][8] = {
        { {1,1, 1,1, 2,13, 3,7}, {1,1, 1,1, 11,9, 5,3}, {1,1, 1,1, 4,17, 7,11} },
        { {3,7, 5,11, 2,13, 7,5}, {7,5, 3,4, 11,9, 2,9}, {5,3, 13,7, 4,17, 9,4} } };
    bool found = false;
    Expr *k = NULL, *w = NULL;
    for (int set = 0; set < 2 && !found; set++) {
        Expr *Ai[3], *Bi[3], *Ci[3];
        bool ok = true;
        for (int i = 0; i < 3; i++) {
            Ai[i] = sc_at(A, vars, pts[set][i]);
            Bi[i] = sc_at(B, vars, pts[set][i]);
            Ci[i] = sc_at(R, vars, pts[set][i]);
            if (!sc_finite(Ai[i]) || !sc_finite(Bi[i]) || !sc_finite(Ci[i])) ok = false;
        }
        for (int j = 1; j < 3 && ok && !k; j++) {
            Expr* det = eval_and_free(ds_call1(SYM_Together, ds_call2(SYM_Subtract,
                            ds_call2(SYM_Times, expr_copy(Ci[0]), expr_copy(Bi[j])),
                            ds_call2(SYM_Times, expr_copy(Bi[0]), expr_copy(Ci[j])))));
            if (ds_is_zero(det)) { expr_free(det); continue; }
            Expr* kn = ds_call2(SYM_Subtract, ds_call2(SYM_Times, expr_copy(Ai[0]), expr_copy(Ci[j])),
                                              ds_call2(SYM_Times, expr_copy(Ci[0]), expr_copy(Ai[j])));
            Expr* wn = ds_call2(SYM_Subtract, ds_call2(SYM_Times, expr_copy(Ai[0]), expr_copy(Bi[j])),
                                              ds_call2(SYM_Times, expr_copy(Ai[j]), expr_copy(Bi[0])));
            k = ds_simplify(eval_and_free(ds_call1(SYM_Together,
                    ds_call2(SYM_Times, kn, ds_call2(SYM_Power, expr_copy(det), expr_new_integer(-1))))));
            w = ds_simplify(eval_and_free(ds_call1(SYM_Together,
                    ds_call2(SYM_Times, wn, ds_call2(SYM_Power, det, expr_new_integer(-1))))));
        }
        for (int i = 0; i < 3; i++) { expr_free(Ai[i]); expr_free(Bi[i]); expr_free(Ci[i]); }
        if (!k) continue;
        if (!sc_finite(k) || !sc_finite(w) || !sc_free4(k, x, Y, P, Q) || !sc_free4(w, x, Y, P, Q)) {
            expr_free(k); expr_free(w); k = w = NULL; continue;
        }
        Expr* Z = eval_and_free(ds_call1(SYM_Together, ds_call2(SYM_Plus, expr_copy(A),
                      ds_call2(SYM_Subtract, ds_call2(SYM_Times, expr_copy(k), expr_copy(B)),
                                             ds_call2(SYM_Times, expr_copy(w), expr_copy((Expr*)R))))));
        found = ds_is_zero(Z);
        expr_free(Z);
        if (!found) { expr_free(k); expr_free(w); k = w = NULL; }
    }
    expr_free(A); expr_free(B);
    if (!found) return false;
    *kout = k; *wout = w;
    return true;
}

/* Collect Inactive[Integrate][f, v] nodes (indefinite). */
static void sc_collect_inert(const Expr* e, Expr** out, int* n, int cap) {
    if (!e || e->type != EXPR_FUNCTION || *n >= cap) return;
    const Expr* h = e->data.function.head;
    if (head_is(h, SYM_Inactive) && h->data.function.arg_count == 1
        && h->data.function.args[0]->type == EXPR_SYMBOL
        && h->data.function.args[0]->data.symbol.name == SYM_Integrate
        && e->data.function.arg_count == 2 && !head_is(e->data.function.args[1], SYM_List)) {
        out[(*n)++] = expr_copy((Expr*)e);
        return;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        sc_collect_inert(e->data.function.args[i], out, n, cap);
}

/* Rewrite each indefinite Inactive[Integrate][f(v), v] (v == u[t]) as the definite
 * Inactive[Integrate][f(s), {s, s0, v}], s0 the first of 1, 2, 3, 1/2 at which f is
 * finite for sampled constants.  Returns NULL when some inert integral is over a
 * different variable or no base point is found.  `G` borrowed. */
static Expr* sc_definite_inert(const Expr* G, const Expr* v) {
    Expr* zs[8]; int nz = 0;
    sc_collect_inert(G, zs, &nz, 8);
    Expr* out = expr_copy((Expr*)G);
    bool ok = true;
    for (int i = 0; i < nz && ok; i++) {
        if (!expr_eq(zs[i]->data.function.args[1], (Expr*)v)) { ok = false; break; }
        const char* sv = intern_symbol("DSolve`scS");
        Expr* f = ds_subst(expr_copy(zs[i]->data.function.args[0]), expr_copy((Expr*)v),
                           expr_new_symbol(sv));
        static const int64_t base[4][2] = { {1,1}, {2,1}, {3,1}, {1,2} };
        Expr* s0 = NULL;
        for (int b = 0; b < 4 && !s0; b++) {
            Expr* val = ds_subst(sc_sym_consts(expr_copy(f), true), expr_new_symbol(sv),
                                 sc_q(base[b][0], base[b][1]));
            val = eval_and_free(val);
            if (sc_finite(val) && !ds_is_zero(val)) s0 = sc_q(base[b][0], base[b][1]);
            expr_free(val);
        }
        if (!s0) { expr_free(f); ok = false; break; }
        Expr* spec = expr_new_function(expr_new_symbol(SYM_List),
                         (Expr*[]){ expr_new_symbol(sv), s0, expr_copy((Expr*)v) }, 3);
        Expr* inact = expr_new_function(expr_new_symbol(SYM_Inactive),
                          (Expr*[]){ expr_new_symbol(SYM_Integrate) }, 1);
        Expr* def = expr_new_function(inact, (Expr*[]){ f, spec }, 2);
        out = ds_subst(out, expr_copy(zs[i]), def);
    }
    for (int i = 0; i < nz; i++) expr_free(zs[i]);
    if (!ok) { expr_free(out); return NULL; }
    return out;
}

/* Free symbols in argument position (parameters), deduplicated. */
static void sc_collect_syms(const Expr* e, const char** names, int* n, int cap) {
    if (!e || *n >= cap) return;
    if (e->type == EXPR_SYMBOL) {
        for (int i = 0; i < *n; i++) if (names[i] == e->data.symbol.name) return;
        names[(*n)++] = e->data.symbol.name;
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    if (e->data.function.head && e->data.function.head->type != EXPR_SYMBOL)
        sc_collect_syms(e->data.function.head, names, n, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        sc_collect_syms(e->data.function.args[i], names, n, cap);
}

/* Every free parameter of R (not x, not Yv, not a named constant) -> a generic real,
 * the values shifted by `off` between draws.  `R` consumed. */
static Expr* sc_subst_params(Expr* R, const char* xv, const char* Yv, double off) {
    const char* skip[] = { xv, Yv, intern_symbol("E"), intern_symbol("Pi"), intern_symbol("I"),
        intern_symbol("EulerGamma"), intern_symbol("Degree"), intern_symbol("GoldenRatio"),
        intern_symbol("Catalan"), intern_symbol("Infinity"), intern_symbol("DSolve`scS") };
    const int nskip = (int)(sizeof(skip) / sizeof(skip[0]));
    const char* syms[64]; int ns = 0;
    sc_collect_syms(R, syms, &ns, 64);
    int pi = 0;
    for (int i = 0; i < ns; i++) {
        bool sk = false;
        for (int j = 0; j < nskip; j++) if (syms[i] == skip[j]) { sk = true; break; }
        if (sk) continue;
        R = ds_subst(R, expr_new_symbol(syms[i]),
                     expr_new_real(0.43 + off + 0.19 * (double)(pi++)));
    }
    return R;
}

/* |N[e]| as a double (NAN when not numeric).  `e` consumed. */
static double sc_absval(Expr* e) {
    e = eval_and_free(ds_call1("Abs", eval_and_free(ds_call1("N", e))));
    double m = (e->type == EXPR_REAL) ? e->data.real
             : (e->type == EXPR_INTEGER) ? (double)e->data.integer : NAN;
    expr_free(e);
    return m;
}

/* Zero test of a residual that cannot be sampled as it stands (an inert quadrature
 * or an arbitrary function survives the differentiation): every inert integral, every
 * undefined-function application and Derivative of one, and every C[k] is abstracted
 * to its own fresh symbol (outermost first), and the result is sampled at MODERATE
 * values of all of them (three draws, all decided samples must be ~0).  Sound in the
 * direction that matters: an identity that holds with those values free holds for
 * their actual values.  Sampling rather than the symbolic zero test, which reads a
 * decaying residual such as E^(-a) (...) as zero.  `R` consumed. */
bool dsolve_abstract_zero(Expr* R, const char* xv, const char* Yv) {
    static const char* const src =
        "Function[{DSolve`scR}, Module[{DSolve`scAt}, "
        "DSolve`scAt = Union[Cases[DSolve`scR, Inactive[Integrate][__] | Derivative[__][_][__] | C[_Integer] | "
        "  (DSolve`scH_Symbol[__] /; FreeQ[Attributes[DSolve`scH], Protected]), "
        "  {0, Infinity}]]; "
        "DSolve`scR /. Thread[DSolve`scAt -> "
        "  Table[Symbol[\"DSolve`scA\" <> ToString[DSolve`scI]], {DSolve`scI, Length[DSolve`scAt]}]]]]";
    Expr* fn = parse_expression(src);
    if (!fn) { expr_free(R); return false; }
    Expr* z = eval_and_free(expr_new_function(fn, (Expr*[]){ R }, 1));
    int small = 0, big = 0;
    for (int d = 0; d < 3; d++) {
        Expr* e = sc_subst_params(expr_copy(z), xv, Yv, 0.071 * d);
        e = ds_subst(e, expr_new_symbol(xv), expr_new_real(1.3 + 0.23 * d));
        e = ds_subst(e, expr_new_symbol(Yv), expr_new_real(0.7 + 0.17 * d));
        double m = sc_absval(e);
        if (isnan(m) || !isfinite(m)) continue;
        if (m < 1e-8) small++; else big++;
    }
    expr_free(z);
    return small >= 2 && big == 0;
}

/* Numeric verification of an implicit relation G(x, y[x]) == C of the SECOND-order
 * equation: along a level curve y' = -G_x/G_y =: s(x, Y) and y'' = s_x + s_Y s.  The
 * residual with these must vanish at every (x, Y) (every point lies on some level
 * curve).  Constants and parameters sampled at generic reals (two draws, as in
 * AutonomousReduction's check: a squared-up relation holds on an open set of the
 * constants); accepted when >= 3 of the 4 points of a draw are ~0 and none is
 * decidedly nonzero. */
bool dsolve_verify_implicit2(const DSolveProblem* P, const Expr* G) {
    const char* yname = P->fun_names[0];
    const char* xv    = P->ind_names[0];
    const char* Yv    = intern_symbol("DSolve`scYv");
    Expr* Gs = ds_subst(expr_copy((Expr*)G), ds_make_funcapp(yname, 0, xv), expr_new_symbol(Yv));
    Expr* Gx = ds_d(expr_copy(Gs), expr_new_symbol(xv));
    Expr* Gy = ds_d(Gs, expr_new_symbol(Yv));
    Expr* s  = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
                   ds_call2(SYM_Times, Gx, ds_call2(SYM_Power, Gy, expr_new_integer(-1)))));
    Expr* s2 = eval_and_free(ds_call2(SYM_Plus, ds_d(expr_copy(s), expr_new_symbol(xv)),
                   ds_call2(SYM_Times, ds_d(expr_copy(s), expr_new_symbol(Yv)), expr_copy(s))));
    Expr* R = expr_copy(P->eq_residuals[0]);
    R = ds_subst(R, ds_make_funcapp(yname, 2, xv), s2);
    R = ds_subst(R, ds_make_funcapp(yname, 1, xv), s);
    R = ds_subst(R, ds_make_funcapp(yname, 0, xv), expr_new_symbol(Yv));
    if (ds_contains(R, yname)) { expr_free(R); return false; }
    /* Not samplable (an inert quadrature or an arbitrary function survives the
     * differentiation): decide symbolically instead. */
    if (ds_has_inactive_integrate(R) || ds_has_undefined_function(R))
        return dsolve_abstract_zero(R, xv, Yv);

    bool pass = false;
    for (int draw = 0; draw < 2 && !pass; draw++) {
        Expr* Rd = expr_copy(R);
        for (int k = 1; k <= 4; k++)
            Rd = ds_subst(Rd, ds_const(k), expr_new_real(draw == 0 ? 0.37 + 0.11 * k : -1.9 - 0.23 * k));
        Rd = sc_subst_params(Rd, xv, Yv, 0.0);       /* the remaining free parameters */
        static const double px[4] = { 1.3, 1.7, 2.2, 1.1 };
        static const double pY[4] = { 0.6, 0.9, 1.4, 0.45 };
        int small = 0, big = 0;
        for (int i = 0; i < 4; i++) {
            Expr* e = ds_subst(expr_copy(Rd), expr_new_symbol(xv), expr_new_real(px[i]));
            e = ds_subst(e, expr_new_symbol(Yv), expr_new_real(pY[i]));
            double m = sc_absval(e);
            if (isnan(m) || !isfinite(m)) continue;
            if (m < 1e-7) small++; else if (m > 1e-4) big++;
        }
        expr_free(Rd);
        if (small >= 3 && big == 0) pass = true;
    }
    expr_free(R);
    return pass;
}

/* Every branch of a second-order sub-solve result as a first-integral LHS G with
 * G == C[1]: an implicit {G == C[1]} (G free of C[1], other constants allowed -- the
 * shared dsolve_result_relations requires a constant-free G and so rejects every
 * second-order relation).  All or nothing.  `r` borrowed; *rel owned. */
static bool sc_relations(const Expr* r, const char* un, Expr*** rel, size_t* nrel) {
    if (!head_is(r, SYM_List) || r->data.function.arg_count == 0) return false;
    size_t nb = r->data.function.arg_count;
    Expr** tmp = malloc(nb * sizeof(Expr*));
    size_t nt = 0;
    Expr* c1 = ds_const(1);
    bool ok = true;
    for (size_t b = 0; b < nb && ok; b++) {
        const Expr* br = r->data.function.args[b];
        ok = false;
        if (!head_is(br, SYM_List) || br->data.function.arg_count != 1) break;
        const Expr* it = br->data.function.args[0];
        if (!head_is(it, SYM_Equal) || it->data.function.arg_count != 2) break;
        const Expr* l = it->data.function.args[0];
        const Expr* rr = it->data.function.args[1];
        if (!expr_eq((Expr*)rr, c1) || !ds_contains(l, un)) break;
        Expr* chk = ds_subst(expr_copy((Expr*)l), expr_copy(c1), expr_new_integer(0));
        bool free1 = expr_eq(chk, (Expr*)l);
        expr_free(chk);
        if (!free1) break;
        tmp[nt++] = expr_copy((Expr*)l);
        ok = true;
    }
    expr_free(c1);
    if (!ok) { for (size_t i = 0; i < nt; i++) expr_free(tmp[i]); free(tmp); return false; }
    *rel = tmp; *nrel = nt;
    return true;
}

/* Map a sub-solve result r (function un of tn) back to y(x) for the change of
 * variables  y = A u(T) + B  with T = tmap(x) (A, B, tmap borrowed; NULL means 1, 0,
 * and tn == x respectively).  All-explicit: returns the y-bodies (verified later by
 * dsolve_run).  Otherwise each implicit relation G(T, u) == C[1] -- an inert
 * quadrature over u made definite first -- is mapped by u -> (y[x] - B)/A, T ->
 * tmap, verified at order two (dsolve_verify_implicit2) and stashed for the
 * implicit pass; returns NULL. */
static Expr** sc_map_back(const DSolveProblem* P, const Expr* r, const char* un,
                          const char* tn, const Expr* A, const Expr* B, const Expr* tmap,
                          uint64_t h, size_t* nbranch) {
    const char* yname = P->fun_names[0];
    const char* xv    = P->ind_names[0];
    size_t nbr = head_is(r, SYM_List) ? r->data.function.arg_count : 0;
    if (nbr == 0) return NULL;
    size_t nb = 0;
    Expr** ub = dsolve_extract_applied_bodies((Expr*)r, un, &nb);
    if (ub && nb == nbr) {
        Expr** out = malloc(nb * sizeof(Expr*));
        for (size_t i = 0; i < nb; i++) {
            Expr* b = ub[i];
            if (tmap) b = ds_subst(b, expr_new_symbol(tn), expr_copy((Expr*)tmap));
            if (A) b = ds_call2(SYM_Times, expr_copy((Expr*)A), b);
            if (B) b = ds_call2(SYM_Plus, b, expr_copy((Expr*)B));
            out[i] = eval_and_free(b);
        }
        free(ub);
        *nbranch = nb;
        return out;
    }
    if (ub) { for (size_t i = 0; i < nb; i++) expr_free(ub[i]); free(ub); }
    Expr** rel = NULL; size_t nrel = 0;
    if (ds_has_active_integrate(r) || !sc_relations(r, un, &rel, &nrel)) return NULL;
    /* u -> (y[x] - B)/A */
    Expr* umap = ds_make_funcapp(yname, 0, xv);
    if (B) umap = ds_call2(SYM_Subtract, umap, expr_copy((Expr*)B));
    if (A) umap = ds_call2(SYM_Times, umap, ds_call2(SYM_Power, expr_copy((Expr*)A),
                                                     expr_new_integer(-1)));
    umap = eval_and_free(umap);
    Expr* uoft = ds_make_funcapp(un, 0, tn);
    Expr** keep = malloc(nrel * sizeof(Expr*));
    size_t nk = 0;
    for (size_t i = 0; i < nrel; i++) {
        Expr* G = rel[i];
        if (ds_has_inactive_integrate(G)) {
            Expr* d = sc_definite_inert(G, uoft);
            expr_free(G);
            G = d;
        }
        if (!G) continue;
        G = ds_subst(G, expr_copy(uoft), expr_copy(umap));
        if (tmap) G = ds_subst(G, expr_new_symbol(tn), expr_copy((Expr*)tmap));
        if (ds_contains(G, un) || (tmap && ds_contains(G, tn)) || !dsolve_verify_implicit2(P, G)) {
            expr_free(G); continue;
        }
        keep[nk++] = G;
    }
    free(rel);
    expr_free(uoft); expr_free(umap);
    if (nk > 0) dsolve_stash_put(SC_STASH_TAG, h, keep, nk);
    else free(keep);
    return NULL;
}

static Expr** sc_try_isobaric(const DSolveProblem* P, const Expr* R, const char* Yn,
                              const char* Pn, const char* Qn, uint64_t h, size_t* nbranch) {
    const char* xv    = P->ind_names[0];
    /* an autonomous equation (x absent) is AutonomousReduction's directly; its own
     * x-scaling would only re-derive that reduction (and the weight solve over a
     * radical such as (1 + (b - P)^2)^(3/2) cost 18 s) */
    if (!ds_contains(R, xv) || !sc_algebraic(R, xv, Yn, Pn, Qn)) return NULL;
    Expr *k = NULL, *w = NULL;
    if (!sc_weights(R, xv, Yn, Pn, Qn, &k, &w)) return NULL;

    const char* Un = intern_symbol("DSolve`scU");
    const char* Vn = intern_symbol("DSolve`scV");
    const char* Wn = intern_symbol("DSolve`scW");
    /* Y = x^k U, P = x^(k-1) (k U + V), Q = x^(k-2) (W + (2k-1) V + k (k-1) U) */
    Expr* xk  = ds_call2(SYM_Power, expr_new_symbol(xv), expr_copy(k));
    Expr* xk1 = ds_call2(SYM_Power, expr_new_symbol(xv), ds_call2(SYM_Subtract, expr_copy(k), expr_new_integer(1)));
    Expr* xk2 = ds_call2(SYM_Power, expr_new_symbol(xv), ds_call2(SYM_Subtract, expr_copy(k), expr_new_integer(2)));
    Expr* Qs = ds_call2(SYM_Times, xk2, ds_call2(SYM_Plus, expr_new_symbol(Wn), ds_call2(SYM_Plus,
                   ds_call2(SYM_Times, ds_call2(SYM_Subtract, ds_call2(SYM_Times, expr_new_integer(2), expr_copy(k)),
                                                expr_new_integer(1)), expr_new_symbol(Vn)),
                   ds_call2(SYM_Times, ds_call2(SYM_Times, expr_copy(k),
                                                ds_call2(SYM_Subtract, expr_copy(k), expr_new_integer(1))),
                            expr_new_symbol(Un)))));
    Expr* Ps = ds_call2(SYM_Times, xk1, ds_call2(SYM_Plus,
                   ds_call2(SYM_Times, expr_copy(k), expr_new_symbol(Un)), expr_new_symbol(Vn)));
    Expr* S = expr_copy((Expr*)R);
    S = ds_subst(S, expr_new_symbol(Qn), Qs);
    S = ds_subst(S, expr_new_symbol(Pn), Ps);
    S = ds_subst(S, expr_new_symbol(Yn), ds_call2(SYM_Times, xk, expr_new_symbol(Un)));
    S = ds_call2(SYM_Times, S, ds_call2(SYM_Power, expr_new_symbol(xv),
                                        ds_call2(SYM_Times, expr_new_integer(-1), w)));
    Expr* Ft = eval_and_free(ds_call1(SYM_Numerator, eval_and_free(ds_call1(SYM_Together,
                   eval_and_free(ds_call1("PowerExpand", eval_and_free(S)))))));
    if (!ds_free_of(Ft, xv) || ds_contains(Ft, xv) || !ds_contains(Ft, Wn)) {
        expr_free(Ft); expr_free(k); return NULL;
    }
    const char* un = intern_symbol("DSolve`scu");
    const char* tn = intern_symbol("DSolve`sct");
    Ft = ds_subst(Ft, expr_new_symbol(Wn), ds_make_funcapp(un, 2, tn));
    Ft = ds_subst(Ft, expr_new_symbol(Vn), ds_make_funcapp(un, 1, tn));
    Ft = ds_subst(Ft, expr_new_symbol(Un), ds_make_funcapp(un, 0, tn));
    Expr* ode = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ Ft, expr_new_integer(0) }, 2);
    Expr* r = sc_subsolve(ode, un, tn, SC_SUBSOLVE_SEC + 1);
    /* y = x^k u(Log[x]) */
    Expr* A = ds_call2(SYM_Power, expr_new_symbol(xv), expr_copy(k));
    Expr* tmap = ds_call1(SYM_Log, expr_new_symbol(xv));
    Expr** out = sc_map_back(P, r, un, tn, A, NULL, tmap, h, nbranch);
    expr_free(r); expr_free(k); expr_free(A); expr_free(tmap);
    return out;
}

/* ---- (3) translation x -> x + e, y -> y - rho e : v = y + rho x ---- *
 * R depends on x and y only through y + rho x exactly when R_x == rho R_Y with rho a
 * constant (free of the jet variables).  Then y = v - rho x, y' = v' - rho, y'' = v''
 * makes the equation AUTONOMOUS in v(x): (x - y) y'' + 2 y' (y' + 1) == 0 (rho = -1),
 * (a x + b y) y'^3 + y'' == 0 (rho = a/b). */
static Expr** sc_try_translation(const DSolveProblem* P, const Expr* R, const char* Yn,
                                 const char* Pn, const char* Qn, uint64_t h, size_t* nbranch) {
    const char* xv = P->ind_names[0];
    Expr* Rx = ds_d(expr_copy((Expr*)R), expr_new_symbol(xv));
    Expr* RY = ds_d(expr_copy((Expr*)R), expr_new_symbol(Yn));
    if (ds_is_zero(Rx) || ds_is_zero(RY)) { expr_free(Rx); expr_free(RY); return NULL; }
    Expr* rho = eval_and_free(ds_call1(SYM_Together, ds_call2(SYM_Times, Rx,
                    ds_call2(SYM_Power, RY, expr_new_integer(-1)))));
    if (!sc_finite(rho) || !sc_free4(rho, xv, Yn, Pn, Qn) || ds_is_zero(rho)) {
        expr_free(rho); return NULL;
    }
    rho = ds_simplify(rho);
    const char* un = intern_symbol("DSolve`scu");
    Expr* S = expr_copy((Expr*)R);
    Expr* rx = eval_and_free(ds_call2(SYM_Times, expr_copy(rho), expr_new_symbol(xv)));
    S = ds_subst(S, expr_new_symbol(Qn), ds_make_funcapp(un, 2, xv));
    S = ds_subst(S, expr_new_symbol(Pn), ds_call2(SYM_Subtract, ds_make_funcapp(un, 1, xv),
                                                  expr_copy(rho)));
    S = ds_subst(S, expr_new_symbol(Yn), ds_call2(SYM_Subtract, ds_make_funcapp(un, 0, xv),
                                                  expr_copy(rx)));
    S = eval_and_free(ds_call1(SYM_Numerator, ds_call1(SYM_Together, S)));
    /* the substituted equation must be x-free (mask the v-jet, then test) */
    Expr* Sm = expr_copy(S);
    const char* Mk = intern_symbol("DSolve`scMask");
    for (int o = 2; o >= 0; o--) Sm = ds_subst(Sm, ds_make_funcapp(un, o, xv), expr_new_symbol(Mk));
    bool xfree = ds_free_of(Sm, xv);
    expr_free(Sm);
    if (!xfree || !ds_contains(S, un)) { expr_free(S); expr_free(rho); expr_free(rx); return NULL; }
    Expr* ode = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ S, expr_new_integer(0) }, 2);
    Expr* r = sc_subsolve(ode, un, xv, SC_SUBSOLVE_SEC + 1);
    Expr* B = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), rx));   /* y = v - rho x */
    Expr** out = sc_map_back(P, r, un, xv, NULL, B, NULL, h, nbranch);
    expr_free(r); expr_free(B); expr_free(rho);
    return out;
}

/* A LINEAR equation (R affine in Y, P, Q: every second partial vanishes -- not the
 * Euler-identity test, which a nonlinear degree-1 homogeneous R such as
 * x^2 Q - Sqrt[b Y^2 + a x^2 P^2] also passes) is left alone: a homogeneous one is always
 * y-scale invariant, but its Riccati u = y'/y is no easier than the linear methods
 * that already declined it, and an isobaric linear one is an Euler equation. */
static bool sc_linear(const Expr* R, const char* Yn, const char* Pn, const char* Qn) {
    static const char* const src =
        "Function[{DSolve`scR, DSolve`scY1, DSolve`scP1, DSolve`scQ1}, "
        "Module[{DSolve`scV = {DSolve`scY1, DSolve`scP1, DSolve`scQ1}}, "
        "And @@ Flatten[Table[Expand[Together[D[DSolve`scR, DSolve`scV[[DSolve`scI]], "
        "DSolve`scV[[DSolve`scJ]]]]] === 0, {DSolve`scI, 3}, {DSolve`scJ, DSolve`scI, 3}]]]]";
    Expr* fn = parse_expression(src);
    if (!fn) return false;
    Expr* r = eval_and_free(expr_new_function(fn, (Expr*[]){ expr_copy((Expr*)R),
                  expr_new_symbol(Yn), expr_new_symbol(Pn), expr_new_symbol(Qn) }, 4));
    bool lin = (r->type == EXPR_SYMBOL && r->data.symbol.name == intern_symbol("True"));
    expr_free(r);
    return lin;
}

Expr** dsolve_scaling2_try(DSolveProblem* P, size_t* nbranch) {
    if (getenv("MATHILDA_NO_SCALING2")) return NULL;
    if (g_dsolve_depth > 1 || g_sc_active) return NULL;
    if (P->nfun != 1 || P->neq != 1 || P->max_order[0] != 2) return NULL;
    uint64_t h = expr_hash(P->eq_residuals[0]);
    sc_memo_sync(eval_toplevel_id());
    if (sc_memo_seen(h)) return NULL;
    sc_memo_add(h);                    /* one attempt per top-level evaluation */

    const char* Yn = intern_symbol("DSolve`scY");
    const char* Pn = intern_symbol("DSolve`scP");
    const char* Qn = intern_symbol("DSolve`scQ");
    Expr* R = sc_residual(P, Yn, Pn, Qn);
    if (!R) return NULL;
    if (!ds_contains(R, Yn) || !ds_contains(R, Qn) || sc_linear(R, Yn, Pn, Qn)) {
        expr_free(R); return NULL;
    }
    Expr** out = sc_try_yscale(P, R, Yn, Pn, Qn, nbranch);
    if (!out) out = sc_try_isobaric(P, R, Yn, Pn, Qn, h, nbranch);
    if (!out && !dsolve_stash_has(SC_STASH_TAG, h))
        out = sc_try_translation(P, R, Yn, Pn, Qn, h, nbranch);
    expr_free(R);
    return out;
}

/* The implicit pass: the relations the explicit pass verified and stashed (pinned
 * on its own, it runs the explicit pass first to populate the stash). */
static Expr** sc_implicit_impl(DSolveProblem* P, size_t* nbranch, bool may_solve) {
    if (P->nfun != 1 || P->neq != 1 || P->max_order[0] != 2) return NULL;
    uint64_t h = expr_hash(P->eq_residuals[0]);
    size_t n = 0;
    Expr** G = dsolve_stash_take(SC_STASH_TAG, h, &n);
    if (!G && may_solve) {
        size_t ne = 0;
        Expr** ex = dsolve_scaling2_try(P, &ne);
        if (ex) { for (size_t i = 0; i < ne; i++) expr_free(ex[i]); free(ex); return NULL; }
        G = dsolve_stash_take(SC_STASH_TAG, h, &n);
    }
    if (!G) return NULL;
    *nbranch = n;
    return G;
}
Expr** dsolve_scaling2_implicit_try(DSolveProblem* P, size_t* nbranch) { return sc_implicit_impl(P, nbranch, false); }
static Expr** sc_implicit_pinned_try(DSolveProblem* P, size_t* nbranch) { return sc_implicit_impl(P, nbranch, true); }

static Expr* builtin_dsolve_scaling2(Expr* res) { return dsolve_method_builtin(res, dsolve_scaling2_try); }
static Expr* builtin_dsolve_scaling2_implicit(Expr* res) {
    return dsolve_method_builtin_implicit(res, sc_implicit_pinned_try);
}

void dsolve_scaling2_init(void) {
    tc_register_unwind_counter(&g_sc_active);   /* restored on a timeout unwind */
    symtab_add_builtin("DSolve`ScalingReduction", builtin_dsolve_scaling2);
    symtab_get_def("DSolve`ScalingReduction")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`ScalingReduction",
        "DSolve`ScalingReduction[eqn, y, x] solves a second-order ODE invariant under a "
        "scaling: homogeneous in y (y -> l y) by u = y'/y, a first-order ODE, then "
        "y = C[2] Exp[Integrate[u, x]]; or isobaric (x -> l x, y -> l^k y) by "
        "y = x^k u(Log[x]); or translation-invariant (y depends on x only through y + rho x) "
        "by v = y + rho x -- the last two autonomous in u, v.");
    symtab_add_builtin("DSolve`ScalingReductionImplicit", builtin_dsolve_scaling2_implicit);
    symtab_get_def("DSolve`ScalingReductionImplicit")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`ScalingReductionImplicit",
        "DSolve`ScalingReductionImplicit[eqn, y, x] is DSolve`ScalingReduction when the "
        "isobaric reduced equation solves only implicitly: each first integral "
        "G[t, u[t]] == C[1] mapped back by u -> y/x^k, t -> Log[x] (an inert quadrature "
        "made definite), verified numerically against the second-order equation.");
}
