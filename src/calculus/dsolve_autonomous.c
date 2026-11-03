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
#include "dsolve.h"          /* g_dsolve_depth */
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include "../parse.h"           /* parse_expression: stage-1 normalisation rules */
#include "../message.h"         /* mth_msg_suppress_push: the stage-1 Solve probe */
#include "integrate.h"          /* g_integrate_quiet: silence the speculative quadrature */
#include <stdlib.h>
#include <stdio.h>
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

/* Does `sub` occur in `e` (structurally)?  Both borrowed. */
static bool ar_stage1_fitted(const Expr* pbody, int n);
static bool ar_occurs(const Expr* e, const Expr* sub) {
    if (!e) return false;
    if (expr_eq((Expr*)e, (Expr*)sub)) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (ar_occurs(e->data.function.head, sub)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (ar_occurs(e->data.function.args[i], sub)) return true;
    return false;
}

/* C[lo..hi] <-> the plain symbols DSolve`arK<k> (forward: to symbols).  Solve and
 * Integrate treat a symbol as a parameter but not an applied C[k].  `e` consumed. */
static Expr* ar_sym_consts(Expr* e, int lo, int hi, bool forward) {
    for (int k = lo; k <= hi; k++) {
        char nm[32];
        snprintf(nm, sizeof nm, "DSolve`arK%d", k);
        Expr* s = expr_new_symbol(intern_symbol(nm));
        e = forward ? ds_subst(e, ds_const(k), s) : ds_subst(e, s, ds_const(k));
    }
    return eval_and_free(e);
}

/* The stage-1 sub-solve, every branch: DSolve[eqn, fname[varname], varname] under
 * TimeConstrained, returning the explicit bodies (owned array, *nb set) or NULL.
 * A branch the sub-solve could only give IMPLICITLY, Phi(Y, p[Y]) == C[1], is
 * solved for p[Y] (bounded Solve, messages muted); its explicit roots join the
 * list -- p(y) need not be explicit to the reduced cascade, but it must be for the
 * stage-2 quadrature (the separable reduced ODE p p' == a (1 + p^2)^(3/2) closes
 * only implicitly, -1/Sqrt[1 + p^2] - a y == C[1]).  Every root is re-verified
 * later (ar_num_ok against the original equation), so a spurious root of the
 * squared-up relation is dropped, not shipped.  `eqn` consumed. */
static Expr** run_dsolve_branches(Expr* eqn, const char* fname, const char* varname,
                                  int secs, size_t* nb) {
    *nb = 0;
    Expr* fv   = ds_call1(fname, expr_new_symbol(varname));
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ eqn, expr_copy(fv), expr_new_symbol(varname) }, 3);
    Expr* g = expr_new_function(expr_new_symbol("TimeConstrained"),
                  (Expr*[]){ call, expr_new_integer(secs < 1 ? 1 : secs),
                             expr_new_symbol(intern_symbol("$Aborted")) }, 3);
    Expr* r = eval_and_free(g);
    if (!head_is(r, SYM_List) || r->data.function.arg_count == 0) {
        expr_free(r); expr_free(fv); return NULL;
    }
    size_t m = r->data.function.arg_count, cap = m + 8, c = 0;
    Expr** out = malloc(cap * sizeof(Expr*));
    for (size_t i = 0; i < m; i++) {
        Expr* inner = r->data.function.args[i];
        if (!head_is(inner, SYM_List) || inner->data.function.arg_count != 1) continue;
        Expr* it = inner->data.function.args[0];
        if (head_is(it, SYM_Rule) && it->data.function.arg_count == 2) {
            /* a Root-object body is no use to the stage-2 quadrature (and the
             * normalisation below cannot reach into its held Function) */
            if (expr_eq(it->data.function.args[0], fv) && c < cap
                && !ds_has_head(it->data.function.args[1], intern_symbol("Root")))
                out[c++] = expr_copy(it->data.function.args[1]);
            continue;
        }
        if (!head_is(it, SYM_Equal) || it->data.function.arg_count != 2) continue;
        if (ds_has_head(it, SYM_Integrate) || ds_has_inactive_integrate(it)) continue;
        /* implicit: Solve[Phi == C, p[Y]] -- simplified first, since the separable
         * relation often arrives in a half-angle substituted form that Solve cannot
         * invert but whose simplified spelling is a plain radical. */
        /* Solve on plain symbols: an applied C[k] or p[Y] is not a parameter or an
         * unknown to Solve's radical isolation (Solve[-a Y^2/2 - 1/Sqrt[1 + p^2] == C[1],
         * p] stays unevaluated; with a symbol K it returns both roots). */
        const char* Pv = intern_symbol("DSolve`arPv");
        Expr* rel = ds_subst(expr_copy(it), expr_copy(fv), expr_new_symbol(Pv));
        rel = ar_sym_consts(rel, 1, 4, true);
        rel = eval_and_free(ds_call1("Simplify", rel));
        Expr* sv = ds_call2("Solve", rel, expr_new_symbol(Pv));
        Expr* sg = expr_new_function(expr_new_symbol("TimeConstrained"),
                       (Expr*[]){ sv, expr_new_integer(2),
                                  expr_new_symbol(intern_symbol("$Aborted")) }, 3);
        mth_msg_suppress_push();
        Expr* sol = eval_and_free(sg);
        mth_msg_suppress_pop();
        if (head_is(sol, SYM_List)) {
            for (size_t k = 0; k < sol->data.function.arg_count; k++) {
                Expr* br = sol->data.function.args[k];
                if (!head_is(br, SYM_List) || br->data.function.arg_count != 1) continue;
                Expr* rule = br->data.function.args[0];
                if (!head_is(rule, SYM_Rule) || rule->data.function.arg_count != 2) continue;
                Expr* lhs = rule->data.function.args[0];
                if (lhs->type != EXPR_SYMBOL || lhs->data.symbol.name != Pv) continue;
                /* Solve's roots of a squared-up radical arrive expanded; a bounded
                 * Simplify recovers the readable radical (value-preserving). */
                Expr* v = expr_copy(rule->data.function.args[1]);
                v = eval_and_free(expr_new_function(expr_new_symbol("TimeConstrained"),
                        (Expr*[]){ ds_call1("Simplify", expr_copy(v)), expr_new_integer(1),
                                   v }, 3));
                v = ar_sym_consts(v, 1, 4, false);
                if (ds_contains(v, Pv) || ds_contains(v, fname)
                    || ds_has_head(v, intern_symbol("Root"))
                    || ds_contains(v, intern_symbol("ConditionalExpression"))) {
                    expr_free(v); continue;
                }
                if (c >= cap) { cap *= 2; out = realloc(out, cap * sizeof(Expr*)); }
                out[c++] = v;
            }
        }
        expr_free(sol);
    }
    expr_free(r); expr_free(fv);
    if (c == 0) { free(out); return NULL; }
    *nb = c;
    return out;
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
    Expr* R0 = expr_copy(R);              /* kept for the abstracted fallback below */
    /* instantiate the constants C[1..n] and every remaining free parameter at
     * distinct generic reals; sample Ysym; require |R| ~ 0.  Two constant draws: a
     * reduced first integral obtained by SQUARING up a radical relation (stage-1 Solve
     * on -a y^2/2 - 1/Sqrt[1 + p^2] == C) is a solution exactly where the radical's
     * sign condition holds -- here C + a y^2/2 < 0 -- so a branch is accepted when it
     * verifies on the whole sample set of EITHER draw (an open set of constants), and
     * never when any sample of that draw is decidedly nonzero. */
    bool decided_any = false, pass = false;
    for (int draw = 0; draw < 2 && !pass; draw++) {
        Expr* Rd = expr_copy(R);
        for (int k = 1; k <= n; k++)
            Rd = ds_subst(Rd, ds_const(k),
                          expr_new_real(draw == 0 ? 0.37 + 0.11 * k : -1.9 - 0.23 * k));
        const char* skip[] = { Ysym, xvar, intern_symbol("E"), intern_symbol("Pi"),
            intern_symbol("I"), intern_symbol("EulerGamma"), intern_symbol("Degree") };
        const int nskip = (int)(sizeof(skip) / sizeof(skip[0]));
        const char* syms[64]; int ns = 0;
        ar_collect_params(Rd, syms, &ns, 64);
        int pi = 0;
        for (int i = 0; i < ns; i++) {
            bool sk = false;
            for (int j = 0; j < nskip; j++) if (syms[i] == skip[j]) { sk = true; break; }
            if (sk) continue;
            Rd = ds_subst(Rd, expr_new_symbol(syms[i]), expr_new_real(0.43 + 0.19 * (pi++)));
        }
        const double ys[] = { 1.3, 0.7, 2.1, 1.7 };
        int small = 0, big = 0;
        for (int i = 0; i < 4; i++) {
            Expr* e = ds_subst(expr_copy(Rd), expr_new_symbol(Ysym), expr_new_real(ys[i]));
            e = eval_and_free(ds_call1("Abs", eval_and_free(ds_call1("N", e))));
            double m = (e->type == EXPR_REAL) ? e->data.real
                     : (e->type == EXPR_INTEGER) ? (double)e->data.integer : NAN;
            expr_free(e);
            if (isnan(m) || !isfinite(m)) continue;
            if (m < 1e-6) small++; else if (m > 1e-3) big++;
        }
        expr_free(Rd);
        if (small + big >= 2) decided_any = true;
        if (small >= 2 && big == 0) pass = true;
    }
    expr_free(R);
    /* Numerically undecidable -- the reduced p carries an ARBITRARY function f[y] or
     * an inert quadrature of one (g[y] + f[y] y'^2 + y'' == 0 gives p^2 == E^(-2 Int f)
     * (C - 2 Int E^(2 Int f) g)), which N cannot evaluate at any sample.  Then decide
     * the residual symbolically: Together telescopes the FTC derivative of each inert
     * integral against its integrand, and only a decided ZERO is accepted --
     * decided by sampling with each inert integral / arbitrary function abstracted
     * to a free symbol (dsolve_abstract_zero), never by the symbolic zero test,
     * which reads a decaying E^(-a) (...) residual as zero. */
    if (!decided_any) return dsolve_abstract_zero(R0, Ysym, Ysym);
    expr_free(R0);
    return pass;
}

/* Stage-1 normalisation.  The reduced first-order ODE in p(y) is typically solved
 * through an exponential integrating factor, so p comes back spelled with
 * E^(c (C[k] + d Log[y])) factors -- correct, but the Log under the exponent trips the
 * stage-2 spin guard (so a quadrature that is elementary never gets tried) and leaves
 * the printed answer unreadable.  Two exact rewrites:
 *   1. split each exponent that carries a Log, E^(u1 + u2 + ...) -> E^u1 E^u2 ...,
 *      after which E^(d Log[y]) evaluates to y^d; then combine each additive symbolic
 *      exponent over a common denominator (y^(-8a^2/(4a+b) - 2ab/(4a+b)) -> y^(-2a));
 *   2. absorb a stage-1 constant that occurs only as E^(alpha_i C[k]), with every
 *      alpha_i a rational multiple of one alpha_0, into a renamed C[k]: the
 *      substitution C[k] -> Log[C[k]]/alpha_0 (so E^(alpha_i C[k]) -> C[k]^(alpha_i/
 *      alpha_0)), a reparametrisation of an arbitrary constant, so the solution
 *      family is unchanged.
 * `pbody` consumed; result owned. */
static Expr* ar_normalise_stage1(Expr* pbody, int n) {
    static const char* const nrm_src =
        "Function[{e}, (e /. Power[E, s_] /; !FreeQ[s, Log] :> "
        "With[{u = Expand[s]}, If[Head[u] === Plus, Times @@ (Power[E, #] & /@ List @@ u), "
        "Power[E, u]]]) /. Power[b_, s_Plus] :> Power[b, Cancel[Together[s]]]]";
    static const char* const absorb_src =
        "Function[{t, c}, Module[{al = Union[Cases[t, Power[E, a_. c] :> a, {0, Infinity}]], a0}, "
        "If[al === {} || !FreeQ[t /. Power[E, a_. c] :> 1, c], t, "
        "a0 = First[al]; "
        "If[And @@ (MatchQ[Together[#/a0], _Integer | _Rational] & /@ al), "
        "t /. Power[E, a_. c] :> c^Together[a/a0], t]]]]";
    if (!ds_contains(pbody, intern_symbol("E"))) return pbody;
    Expr* nrm = parse_expression(nrm_src);
    if (!nrm) return pbody;
    pbody = eval_and_free(expr_new_function(nrm, (Expr*[]){ pbody }, 1));
    for (int k = 2; k <= n; k++) {
        Expr* absorb = parse_expression(absorb_src);
        if (!absorb) break;
        pbody = eval_and_free(expr_new_function(absorb, (Expr*[]){ pbody, ds_const(k) }, 2));
    }
    return pbody;
}

/* True when the stage-2 quadrature Integrate[1/p, Y] is PROVABLY non-elementary:
 * p == c Sqrt[P(Y)] (c free of Y) with P a polynomial of degree >= 3 that is
 * squarefree in Y.  By Chebyshev/Abel, Integrate[P^(-1/2), Y] for such P is an
 * elliptic (deg 3, 4) or hyperelliptic integral, never elementary -- this is the
 * whole y'' == f(y) polynomial family (p^2 == 2 Integrate[f] + C).  Attempting the
 * quadrature (2 s budget) and the stage-2 separable sub-solve (5 s) costs the full
 * budget only to come back inert, which pushed these past the caller's time limit;
 * recognising the shape lets the explicit pass decline at once and the implicit
 * companion emit the inert first integral directly.  A symbolic-coefficient P is
 * judged generically (PolynomialGCD over the rational functions of the
 * parameters), which is the general solution's own genericity.  `pbody` borrowed. */
static bool ar_elliptic_quadrature(const Expr* pbody, const char* Ysym) {
    static const char* const src =
        "Function[{DSolve`arF, DSolve`arV}, Module[{DSolve`arG, DSolve`arP}, "
        "DSolve`arG = Together[DSolve`arF^2]; DSolve`arP = Numerator[DSolve`arG]; "
        "MatchQ[DSolve`arF, Power[_, 1/2] | "
        "  Times[DSolve`arc_ /; FreeQ[DSolve`arc, DSolve`arV], Power[_, 1/2]]] && "
        "FreeQ[Denominator[DSolve`arG], DSolve`arV] && "
        "PolynomialQ[DSolve`arP, DSolve`arV] && Exponent[DSolve`arP, DSolve`arV] >= 3 && "
        "FreeQ[PolynomialGCD[DSolve`arP, D[DSolve`arP, DSolve`arV]], DSolve`arV]]]";
    if (!ds_has_radical_power(pbody)) return false;
    Expr* fn = parse_expression(src);
    if (!fn) return false;
    Expr* r = eval_and_free(expr_new_function(fn,
                  (Expr*[]){ expr_copy((Expr*)pbody), expr_new_symbol(Ysym) }, 2));
    bool yes = (r->type == EXPR_SYMBOL && r->data.symbol.name == intern_symbol("True"));
    expr_free(r);
    return yes;
}

/* Stage-2 quadrature-spin guard: y'==p(y) is separable with quadrature
 * Integrate[1/p,y], which SPINS uninterruptibly for a non-elementary integrand (a
 * Log under a radical, or a y-denominator radical).  True when Integrate must not be
 * attempted.  The order-2 forms that solve are untouched (rational p; constant-
 * denominator radical; the elliptic quartic passes and fast-declines in stage 2). */
static bool ar_quadrature_may_spin(const Expr* pbody, const char* Ysym) {
    if (ar_elliptic_quadrature(pbody, Ysym)) return true;
    Expr* den = eval_and_free(ds_call1(SYM_Denominator,
                    ds_call1(SYM_Together, expr_copy((Expr*)pbody))));
    bool nonconst_denom = !ds_free_of(den, Ysym);
    expr_free(den);
    bool transc_radical = false;
    if (!nonconst_denom && ds_has_radical_power(pbody)) {
        /* A RADICAL over a transcendental function of y is the other non-elementary
         * quadrature shape, and it slips past the denominator test (measured:
         * Denominator[Together[Sqrt[E^-y (C + E^y (2y-2))]]] == 1, yet the stage-2
         * sub-solve costs 90 s and returns an inert relation anyway).  The elliptic
         * quartic the guard above protects has a POLYNOMIAL radicand, so it is
         * untouched. */
        static const char* const tf[] = { "Log", "Exp", "Sin", "Cos", "Tan",
                                          "Sinh", "Cosh", "Tanh", "ArcTan" };
        for (size_t i = 0; i < sizeof(tf) / sizeof(tf[0]) && !transc_radical; i++)
            transc_radical = ds_has_head(pbody, intern_symbol(tf[i]));
        /* E^y is Power[E, y], not an Exp head, so test the symbol itself. */
        if (!transc_radical) transc_radical = ds_contains(pbody, intern_symbol("E"));
    }
    return nonconst_denom || transc_radical || ds_has_head(pbody, intern_symbol("Log"));
}

/* Shared stage-1 reduction (both try-fns): gate + missing-x pre-gate + solve-for-top
 * + autonomous check + the p=y'(y) derivative chain + solve the reduced order-(n-1)
 * ODE + freeze its constants C[1..n-1] -> C[2..n].  Returns pbody = p(y, C[2..n])
 * (owned) with *n_out set, or NULL (not autonomous / reduced ODE did not close).
 * Sets g_ar_deadline.  No memo here -- each try-fn memoizes its own decline. */
static Expr** ar_reduce_compute(DSolveProblem* P, int* n_out, size_t* nb_out) {
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
    /* D[n] == p d/dy(D[n-1]) carries the factor p, and so does every right side that
     * vanishes on y' == 0.  Divide it out: p' == F/p rather than p p' == F.  The two
     * are the same ODE away from the excluded constant solution p == 0, but the
     * divided spelling is the one the first-order cascade recognises directly --
     * p p' + f(y) p + g(y) p^2 == 0 is p' == -f - g p, LINEAR, which the undivided form
     * never reached in 5 s (Factorable sees the p factor but the solve then spun),
     * and p p' + f p^2 + g == 0 is Bernoulli in 0.08 s instead of 4 s. */
    Expr* lhs1 = eval_and_free(ds_call1(SYM_Together,
                     ds_call2(SYM_Times, expr_copy(D[n]),
                              ds_call2(SYM_Power, ds_make_funcapp(pfun, 0, Ysym),
                                       expr_new_integer(-1)))));
    rhs1 = eval_and_free(ds_call2(SYM_Times, rhs1, ds_call2(SYM_Power,
               ds_make_funcapp(pfun, 0, Ysym), expr_new_integer(-1))));
    for (int k = 1; k <= n; k++) expr_free(D[k]);
    free(D);
    /* Two spellings of the same right side, Expand'ed first (Bernoulli/linear are
     * recognised term by term: p' == -f p - g/p) and, only if that gives no usable
     * branch, Together'ed (the separable 2 (p + 1) (p + 2)/(y p) closes in radicals
     * that way, where the expanded spelling ends in degree-12 Root objects).  The
     * Together'ed-first order was measured: several stage-1 solves then spin (20 s). */
    Expr* rhsT = eval_and_free(ds_call1(SYM_Together, expr_copy(rhs1)));
    rhs1 = eval_and_free(ds_call1(SYM_Expand, rhs1));
    size_t nb = 0;
    Expr* eq1  = expr_new_function(expr_new_symbol(SYM_Equal),
                     (Expr*[]){ expr_copy(lhs1), expr_copy(rhs1) }, 2);
    Expr** bodies = run_dsolve_branches(eq1, pfun, Ysym, (int)(g_ar_deadline - time(NULL)), &nb);
    if (!bodies && !expr_eq(rhsT, rhs1) && !ar_expired()) {
        eq1 = expr_new_function(expr_new_symbol(SYM_Equal),
                  (Expr*[]){ expr_copy(lhs1), expr_copy(rhsT) }, 2);
        bodies = run_dsolve_branches(eq1, pfun, Ysym, (int)(g_ar_deadline - time(NULL)), &nb);
    }
    expr_free(lhs1); expr_free(rhs1); expr_free(rhsT);
    if (!bodies) return NULL;

    /* Freeze stage-1 constants C[1..n-1] -> C[2..n] (leaving C[1] for stage 2 / the
     * implicit quadrature constant); C[k] stays a recognised constant (the fast path
     * the elliptic stage-2 integrand needs -- see the M58 note). */
    /* A branch missing any stage-1 constant is a SINGULAR solution of the reduced ODE
     * (p == b +- I for p' == a y (1 + (b - p)^2)^(3/2)/p, where the radical vanishes),
     * not the general family: it would yield an answer short of n constants. */
    size_t kept = 0;
    for (size_t b = 0; b < nb; b++) {
        int off = 1;
        bodies[b] = dsolve_renumber_constants(bodies[b], n - 1, &off);
        bool general = true;
        for (int k = 2; k <= n && general; k++) {
            Expr* ck = ds_const(k);
            general = ar_occurs(bodies[b], ck);
            expr_free(ck);
        }
        if (!general) { expr_free(bodies[b]); continue; }
        bodies[kept++] = ar_normalise_stage1(bodies[b], n);
    }
    nb = kept;
    if (nb == 0) { free(bodies); return NULL; }
    if (ar_expired()) {
        for (size_t b = 0; b < nb; b++) expr_free(bodies[b]);
        free(bodies);
        return NULL;
    }
    *n_out = n;
    *nb_out = nb;
    return bodies;
}

/* Per-top-level cache of the stage-1 reduction, shared by the explicit pass, the
 * elliptic early implicit pass and the late implicit companion, so the reduced
 * first-order ODE (which can cost seconds) is solved once per equation.  A NULL
 * result (not autonomous / reduced ODE did not close) is cached too. */
#define AR_CACHE_SLOTS 8
typedef struct { uint64_t h; int n; size_t nb; Expr** b; } ARCache;
static ARCache ar_cache[AR_CACHE_SLOTS];
static int ar_cache_n = 0;
static uint64_t ar_cache_epoch = 0;
static void ar_cache_sync(void) {
    uint64_t tid = eval_toplevel_id();
    if (tid == ar_cache_epoch) return;
    for (int i = 0; i < ar_cache_n; i++) {
        for (size_t k = 0; k < ar_cache[i].nb; k++) expr_free(ar_cache[i].b[k]);
        free(ar_cache[i].b);
    }
    ar_cache_n = 0; ar_cache_epoch = tid;
}
/* Copies of the cached (or freshly computed) stage-1 branches; NULL when none. */
static Expr** ar_reduce(DSolveProblem* P, int* n_out, size_t* nb_out) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    ar_cache_sync();
    uint64_t h = expr_hash(P->eq_residuals[0]);
    int slot = -1;
    for (int i = 0; i < ar_cache_n; i++) if (ar_cache[i].h == h) { slot = i; break; }
    if (slot < 0) {
        int n = 0; size_t nb = 0;
        Expr** b = ar_reduce_compute(P, &n, &nb);
        if (ar_cache_n < AR_CACHE_SLOTS) {
            slot = ar_cache_n++;
            ar_cache[slot].h = h; ar_cache[slot].n = n;
            ar_cache[slot].nb = b ? nb : 0; ar_cache[slot].b = b;
        } else {
            *n_out = n; *nb_out = nb;
            return b;
        }
    }
    if (ar_cache[slot].nb == 0) return NULL;
    size_t nb = ar_cache[slot].nb;
    Expr** out = malloc(nb * sizeof(Expr*));
    for (size_t k = 0; k < nb; k++) out[k] = expr_copy(ar_cache[slot].b[k]);
    *n_out = ar_cache[slot].n;
    *nb_out = nb;
    return out;
}

/* M61 — fit the STAGE-1 constants C[2..n] from the point conditions, before stage 2.
 *
 * The reduction p = y'(y) carries the stage-1 constants, and for an initial-value
 * problem the conditions determine them EXACTLY and independently of the quadrature:
 * y(x0) = y0 fixes where to evaluate, and y^(k)(x0) = v_k for k = 1..n-1 gives n-1
 * equations D_k(y0) == v_k in the n-1 unknowns, where D_1 = p, D_{k+1} = p d/dy(D_k)
 * is the same chain the reduction is built from.
 *
 * Doing this first is what turns a whole class from an unfitted answer into a closed
 * form: y'' + 2 y y' == 0, y(0) == 0, y'(0) == 1 reduces to p == C[2] - y^2, whose
 * quadrature Integrate[1/(C[2] - y^2), y] is a symbolic-parameter ArcTanh the stage-2
 * guard turns away -- while with C[2] == 1 fitted it is just ArcTanh[y] and the answer
 * is Tanh[x].  Previously the case fell through to the inert first integral with C[2]
 * still free, which the harness rightly scores as unsolved.
 *
 * Narrow by construction: with no conditions, with a condition missing at some order
 * below n, with conditions at more than one point (a BVP), or when `Solve` does not
 * return a usable value, `pbody` is returned UNCHANGED -- so no currently-passing
 * general solution can move.  `pbody` consumed; result owned. */
static Expr* ar_fit_stage1(const DSolveProblem* P, Expr* pbody, int n, const char* Ysym) {
    if (n < 2 || P->ncond == 0) return pbody;
    /* locate y(x0) and y^(k)(x0), k = 1..n-1, all at the SAME point */
    const Expr* pt   = NULL;
    const Expr* y0   = NULL;
    const Expr* vals[16];
    if (n - 1 > 15) return pbody;
    for (int k = 0; k < n; k++) vals[k] = NULL;
    for (size_t i = 0; i < P->ncond; i++) {
        const DSolveCond* c = &P->conds[i];
        if (c->fi != 0 || c->order < 0 || c->order >= n) return pbody;   /* over-determined */
        if (!pt) pt = c->point;
        else if (!expr_eq((Expr*)pt, (Expr*)c->point)) return pbody;     /* BVP: two points */
        if (c->order == 0) { if (y0) return pbody; y0 = c->value; }
        else { if (vals[c->order]) return pbody; vals[c->order] = c->value; }
    }
    if (!y0) return pbody;
    for (int k = 1; k <= n - 1; k++) if (!vals[k]) return pbody;

    /* the derivative chain in terms of pbody: D_1 = p, D_{k+1} = p d/dy(D_k) */
    Expr** D = malloc((size_t)n * sizeof(Expr*));
    D[0] = expr_copy(pbody);
    for (int k = 1; k < n - 1; k++)
        D[k] = eval_and_free(ds_call2(SYM_Times, expr_copy(pbody),
                   ds_d(expr_copy(D[k - 1]), expr_new_symbol(Ysym))));

    /* the system {D_k(y0) == v_k} and the unknowns {C[2], ..., C[n]} */
    Expr** eqs  = malloc((size_t)(n - 1) * sizeof(Expr*));
    Expr** unks = malloc((size_t)(n - 1) * sizeof(Expr*));
    for (int k = 1; k <= n - 1; k++) {
        Expr* at = ds_subst(expr_copy(D[k - 1]), expr_new_symbol(Ysym), expr_copy((Expr*)y0));
        eqs[k - 1]  = expr_new_function(expr_new_symbol(SYM_Equal),
                          (Expr*[]){ at, expr_copy((Expr*)vals[k]) }, 2);
        unks[k - 1] = ds_const(k + 1);
    }
    for (int k = 0; k < n - 1; k++) expr_free(D[k]);
    free(D);
    /* One unknown -> the SCALAR Solve form.  Not cosmetic: Solve inverts the radical
     * Sqrt[E^-y (C[2] + E^y (2y - 2))] == 1 for C[2] in the scalar spelling and
     * bubbles unevaluated in the one-element list spelling. */
    Expr* sol;
    mth_msg_suppress_push();             /* an internal fit probe: no Solve::ifun */
    if (n - 1 == 1) {
        sol = ds_solve(eqs[0], unks[0]);
        free(eqs); free(unks);
    } else {
        sol = ds_solve(expr_new_function(expr_new_symbol(SYM_List), eqs, (size_t)(n - 1)),
                       expr_new_function(expr_new_symbol(SYM_List), unks, (size_t)(n - 1)));
        /* eqs/unks arrays are consumed by expr_new_function */
    }
    mth_msg_suppress_pop();

    /* Use the first branch only when it assigns EVERY stage-1 constant a value free of
     * the remaining unknowns; otherwise leave pbody alone. */
    Expr* out = pbody;
    if (head_is(sol, SYM_List) && sol->data.function.arg_count >= 1) {
        Expr* br = sol->data.function.args[0];
        if (head_is(br, SYM_List) && br->data.function.arg_count == (size_t)(n - 1)) {
            Expr* cand = expr_copy(pbody);
            bool ok = true;
            for (size_t i = 0; i < br->data.function.arg_count && ok; i++) {
                Expr* rule = br->data.function.args[i];
                if (!head_is(rule, SYM_Rule) || rule->data.function.arg_count != 2) ok = false;
                else cand = ds_subst(cand, expr_copy(rule->data.function.args[0]),
                                           expr_copy(rule->data.function.args[1]));
            }
            if (ok && ds_free_of(cand, intern_symbol("C"))
                   && !ds_contains(cand, intern_symbol("Indeterminate"))
                   && !ds_contains(cand, intern_symbol("ComplexInfinity"))) {
                expr_free(pbody); out = cand;
            } else expr_free(cand);
        }
    }
    expr_free(sol);
    return out;
}

/* The stage-2 quadrature Integrate[invp, Ysym], attempted under a short budget with
 * nonelem messages muted (a decline here is expected, not news).  Returns the
 * antiderivative (owned) only when it closed (no Integrate left) and its derivative
 * re-checks against the integrand; otherwise NULL.  `invp` borrowed.
 *
 * The stage-1 constants C[2..n] are replaced by plain symbols for the whole attempt
 * and restored in the result: an applied C[k] is not a symbol to Integrate's
 * parameter handling (Integrate[1/Sqrt[(K - 8/3 Sqrt[y]) y^(5/2)], y] declines in
 * 0.07 s with a symbol K and spends the whole budget with C[2]), and the zero test
 * of the re-check samples symbols but not C[k] applications.
 *
 * Results (closed or not) are cached per top-level evaluation, keyed by the
 * integrand, so the explicit pass, the early and the late implicit passes pay for
 * each quadrature once. */
#define AR_QUAD_SECS 2
#define AR_QCACHE_SLOTS 16
typedef struct { uint64_t h; Expr* F; bool closed; } ARQCache;  /* F NULL: not verified */
static ARQCache ar_qcache[AR_QCACHE_SLOTS];
static int ar_qcache_n = 0;
static uint64_t ar_qcache_epoch = 0;
static void ar_qcache_sync(void) {
    uint64_t tid = eval_toplevel_id();
    if (tid == ar_qcache_epoch) return;
    for (int i = 0; i < ar_qcache_n; i++) if (ar_qcache[i].F) expr_free(ar_qcache[i].F);
    ar_qcache_n = 0; ar_qcache_epoch = tid;
}
static Expr* ar_symbolize_constants(Expr* e, int n, bool forward) {
    return ar_sym_consts(e, 2, n, forward);
}
static Expr* ar_try_quadrature_raw(const Expr* invp, const char* Ysym, bool* closed) {
    *closed = false;
    Expr* call = ds_call2(SYM_Integrate, expr_copy((Expr*)invp), expr_new_symbol(Ysym));
    Expr* g = expr_new_function(expr_new_symbol("TimeConstrained"),
                  (Expr*[]){ call, expr_new_integer(AR_QUAD_SECS),
                             expr_new_symbol(intern_symbol("$Aborted")) }, 3);
    g_integrate_quiet++;
    Expr* F = eval_and_free(g);
    g_integrate_quiet--;
    if (!F || ds_has_head(F, SYM_Integrate)
        || ds_contains(F, intern_symbol("$Aborted"))) {
        if (F) expr_free(F);
        return NULL;
    }
    *closed = true;      /* Integrate answered; the re-check below may still be undecided */
    Expr* chk = eval_and_free(ds_call2(SYM_Subtract,
                    ds_d(expr_copy(F), expr_new_symbol(Ysym)), expr_copy((Expr*)invp)));
    bool ok = ds_is_zero(chk);
    expr_free(chk);
    if (!ok) { expr_free(F); return NULL; }
    return F;
}
/* `closed` (may be NULL): whether Integrate returned a closed form at all -- a
 * closed form whose derivative re-check is merely undecided is NOT evidence that the
 * quadrature is non-elementary, so the callers that skip stage 2 or claim the inert
 * relation early act only on !closed. */
static Expr* ar_try_quadrature(const Expr* invp, const char* Ysym, int n, bool* closed) {
    Expr* sym = ar_symbolize_constants(expr_copy((Expr*)invp), n, true);
    uint64_t h = expr_hash(sym);
    ar_qcache_sync();
    for (int i = 0; i < ar_qcache_n; i++) if (ar_qcache[i].h == h) {
        expr_free(sym);
        if (closed) *closed = ar_qcache[i].closed;
        return ar_qcache[i].F ? ar_symbolize_constants(expr_copy(ar_qcache[i].F), n, false)
                              : NULL;
    }
    bool cl = false;
    Expr* F = ar_try_quadrature_raw(sym, Ysym, &cl);
    expr_free(sym);
    if (closed) *closed = cl;
    if (ar_qcache_n < AR_QCACHE_SLOTS) {
        ar_qcache[ar_qcache_n].h = h;
        ar_qcache[ar_qcache_n].F = F ? expr_copy(F) : NULL;
        ar_qcache[ar_qcache_n].closed = cl;
        ar_qcache_n++;
    }
    return F ? ar_symbolize_constants(F, n, false) : NULL;
}

/* Stage 2 for ONE stage-1 branch: solve y' == p(y) (fitted first; declines on a
 * quadrature that would spin).  Returns the y-body (owned) or NULL.  `pbody`
 * consumed. */
static Expr* ar_stage2(DSolveProblem* P, Expr* pbody, int n) {
    const char* xvar  = P->ind_names[0];
    const char* yname = P->fun_names[0];
    const char* Ysym  = intern_symbol("DSolve`arY");

    /* Fit the stage-1 constants from the conditions FIRST: a numeric p turns several
     * symbolic-parameter integrands the guard below would reject into elementary ones. */
    pbody = ar_fit_stage1(P, pbody, n, Ysym);

    /* Decline fast on a quadrature that would spin -- the IMPLICIT companion below
     * picks these up and returns the inert first integral. */
    if (ar_quadrature_may_spin(pbody, Ysym)) { expr_free(pbody); return NULL; }

    /* The separable stage-2 sub-solve can only close if Integrate[1/p, y] does: try
     * the quadrature first (bounded, cached for the implicit passes) and skip the
     * sub-solve -- whose own separable attempt spends its whole 5 s budget on the
     * same integral -- when it does not. */
    {
        Expr* invp = eval_and_free(ds_call2(SYM_Power, expr_copy(pbody), expr_new_integer(-1)));
        bool closed = false;
        Expr* F = ar_try_quadrature(invp, Ysym, n, &closed);
        expr_free(invp);
        if (F) expr_free(F);
        if (!closed) { expr_free(pbody); return NULL; }
    }

    bool fitted = (P->ncond > 0) && ds_free_of(pbody, intern_symbol("C"));
    Expr* pOfY = ds_subst(pbody, expr_new_symbol(Ysym), ds_make_funcapp(yname, 0, xvar));
    Expr* eq2  = expr_new_function(expr_new_symbol(SYM_Equal),
                     (Expr*[]){ ds_make_funcapp(yname, 1, xvar), pOfY }, 2);
    /* When the stage-1 constants are fitted, hand the value condition y(x0)==y0 to the
     * stage-2 sub-solve as well, so the body comes back with NO generated constant at
     * all.  Otherwise the substrate is left to fit the remaining C[1] out of the full
     * condition set, and for y'' + 2 y y' == 0 that means asking Solve to invert
     * {Tanh[C[1]] == 0, Sech[C[1]]^2 == 1} -- which it cannot decide, so the whole
     * (correct) branch is declined and the case falls through to the inert companion.
     * The first-order sub-solve inverts Tanh[x + c] against one condition directly. */
    Expr* stage2 = eq2;
    if (fitted) {
        const DSolveCond* c0 = NULL;
        for (size_t i = 0; i < P->ncond; i++)
            if (P->conds[i].fi == 0 && P->conds[i].order == 0) { c0 = &P->conds[i]; break; }
        if (c0) {
            Expr* icl = expr_new_function(expr_new_symbol(yname),
                            (Expr*[]){ expr_copy(c0->point) }, 1);
            Expr* ic  = expr_new_function(expr_new_symbol(SYM_Equal),
                            (Expr*[]){ icl, expr_copy(c0->value) }, 2);
            stage2 = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ eq2, ic }, 2);
        }
    }
    Expr* ybody = run_dsolve_applied(stage2, yname, xvar, (int)(g_ar_deadline - time(NULL)));
    if (ybody && ds_free_of(ybody, xvar)) { expr_free(ybody); ybody = NULL; }
    return ybody;
}

/* Explicit: reduce, then (elementary-quadrature only) solve y'==p(y) by separation,
 * on the first stage-1 branch that closes. */
Expr** dsolve_autonomous_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    uint64_t memo_h = expr_hash(P->eq_residuals[0]);
    ar_memo_sync(eval_toplevel_id());
    if (ar_memo_seen(memo_h)) return NULL;

    int n = 0; size_t nb = 0;
    Expr** pb = ar_reduce(P, &n, &nb);
    if (!pb) { ar_memo_add(memo_h); return NULL; }
    g_ar_deadline = time(NULL) + 5;

    Expr** out = malloc(nb * sizeof(Expr*));
    size_t no = 0;
    /* IVP: fit every branch first; once one fits, the branches the conditions cannot
     * fit (the wrong sign of p(y0)) are dropped rather than sent through stage 2 and
     * a doomed constant fit. */
    if (P->ncond > 0) {
        const char* Ysym = intern_symbol("DSolve`arY");
        bool any_fit = false;
        for (size_t b = 0; b < nb; b++) {
            pb[b] = ar_fit_stage1(P, pb[b], n, Ysym);
            any_fit = any_fit || ar_stage1_fitted(pb[b], n);
        }
        if (any_fit)
            for (size_t b = 0; b < nb; b++)
                if (pb[b] && !ar_stage1_fitted(pb[b], n)) { expr_free(pb[b]); pb[b] = NULL; }
    }
    /* The first stage-1 branch whose stage 2 closes is the answer (as before the
     * branches were collected): for an explicit solution the other sign is the same
     * family up to the constants (y y''' == y' y'' would come back twice).  The other
     * branches matter only when this one fails. */
    for (size_t b = 0; b < nb; b++) {
        if (!pb[b]) continue;
        if (no > 0 || ar_expired()) { expr_free(pb[b]); continue; }
        Expr* yb = ar_stage2(P, pb[b], n);
        if (yb) out[no++] = yb;
    }
    free(pb);
    if (no == 0) { free(out); ar_memo_add(memo_h); return NULL; }
    *nbranch = no;
    return out;
}

/* True when the (fitted) stage-1 body carries none of the constants C[2..n]. */
static bool ar_stage1_fitted(const Expr* pbody, int n) {
    for (int k = 2; k <= n; k++) {
        Expr* ck = ds_const(k);
        bool has = ar_occurs(pbody, ck);
        expr_free(ck);
        if (has) return false;
    }
    return true;
}

/* Implicit companion (M59): where the explicit stage-2 quadrature is non-elementary,
 * return the first integral Integrate[1/p, y] - x == C[1] with the integral kept INERT
 * (Inactive[Integrate] -- dodges the uninterruptible Integrate cascade; verified by the
 * FTC/implicit-function rule).  One relation per stage-1 branch (p == +-Sqrt[...]).
 * A numeric self-verify per branch guarantees correctness.
 *
 * `elliptic_only`: the EARLY pass (right after the explicit method) claims the
 * equation only when every branch's quadrature is provably non-elementary
 * (ar_elliptic_quadrature): there the inert relation IS the answer, and the
 * integrating-factor and point-symmetry searches between the two slots only burn
 * the budget (6 s on y'' == a0 + a1 y + a2 y^2 + a3 y^3) before arriving at nothing
 * better. */
static Expr** ar_implicit_impl(DSolveProblem* P, size_t* nbranch, bool elliptic_only) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    uint64_t memo_h = expr_hash(P->eq_residuals[0])
                    ^ (elliptic_only ? 0x7F4A7C159E3779B9ull : 0x9E3779B97F4A7C15ull);
    ar_memo_sync(eval_toplevel_id());
    if (ar_memo_seen(memo_h)) return NULL;
    int n = 0; size_t nb = 0;
    Expr** pb = ar_reduce(P, &n, &nb);
    if (!pb) { ar_memo_add(memo_h); return NULL; }

    const char* xvar  = P->ind_names[0];
    const char* yname = P->fun_names[0];
    const char* Ysym  = intern_symbol("DSolve`arY");

    bool* fitted = malloc(nb * sizeof(bool));
    Expr** out = malloc(nb * sizeof(Expr*));
    size_t no = 0;
    bool all_elliptic = true;
    for (size_t b = 0; b < nb; b++) {
        Expr* pbody = pb[b];
        if (!all_elliptic || !ar_num_ok(P, pbody, n)) { expr_free(pbody); continue; }
        /* Fit the stage-1 constants here too: an inert first integral that still
         * carries a free C[k] is the general solution with the conditions unapplied --
         * i.e. unsolved.  (ar_num_ok ran on the unfitted p, the more general
         * statement, so the fitted p inherits it.) */
        pbody = ar_fit_stage1(P, pbody, n, Ysym);
        bool fit_ok = P->ncond == 0 || ar_stage1_fitted(pbody, n);
        /* G = Integrate[1/pbody, Ysym] /. Ysym -> y[x]  -  x, with the quadrature done
         * when it is elementary and safe to attempt, otherwise kept inert. */
        bool elliptic = ar_elliptic_quadrature(pbody, Ysym);
        bool may_spin = elliptic || ar_quadrature_may_spin(pbody, Ysym);
        Expr* invp = eval_and_free(ds_call2(SYM_Power, pbody, expr_new_integer(-1)));  /* consumes pbody */
        Expr* integ = NULL;
        bool closed = false;
        if (!may_spin) integ = ar_try_quadrature(invp, Ysym, n, &closed);
        /* The early pass claims a branch whose quadrature is provably elliptic, or
         * one Integrate DECLINED (returned unevaluated in budget -- the same
         * quadrature the explicit pass needed, so nothing later can close it); a
         * spin-guarded or closed quadrature is left to the late slot -- except
         * inside a sub-solve (depth > 1: ScalingReduction's autonomous u-equation),
         * where the inert relation is what the caller maps back anyway and the
         * point-symmetry search in between cost 4 s per call. */
        if (elliptic_only && !(elliptic || (!may_spin && !closed)
                               || (g_dsolve_depth > 1 && !integ))) all_elliptic = false;
        if (!all_elliptic) { expr_free(invp); if (integ) expr_free(integ); continue; }
        if (!integ) {
            Expr* inactHead = expr_new_function(expr_new_symbol(SYM_Inactive),
                                  (Expr*[]){ expr_new_symbol(SYM_Integrate) }, 1);
            integ = expr_new_function(inactHead, (Expr*[]){ invp, expr_new_symbol(Ysym) }, 2);
        } else expr_free(invp);
        integ = ds_subst(integ, expr_new_symbol(Ysym), ds_make_funcapp(yname, 0, xvar));
        fitted[no] = fit_ok;
        out[no++] = eval_and_free(ds_call2(SYM_Subtract, integ, expr_new_symbol(xvar)));
    }
    free(pb);
    /* For an IVP, a branch whose stage-1 constants the conditions could not fix is
     * the wrong-sign branch (p(y0) == -v0 has no solution) once a sibling DID fit:
     * shipping it would be the general relation with the conditions unapplied. */
    {
        bool any_fit = false;
        for (size_t i = 0; i < no; i++) any_fit = any_fit || fitted[i];
        if (any_fit) {
            size_t k = 0;
            for (size_t i = 0; i < no; i++) {
                if (fitted[i]) out[k++] = out[i];
                else expr_free(out[i]);
            }
            no = k;
        }
    }
    free(fitted);
    if (no == 0 || !all_elliptic) {
        for (size_t i = 0; i < no; i++) expr_free(out[i]);
        free(out);
        ar_memo_add(memo_h);
        return NULL;
    }
    *nbranch = no;
    return out;
}

Expr** dsolve_autonomous_implicit_try(DSolveProblem* P, size_t* nbranch) {
    return ar_implicit_impl(P, nbranch, false);
}
Expr** dsolve_autonomous_elliptic_try(DSolveProblem* P, size_t* nbranch) {
    return ar_implicit_impl(P, nbranch, true);
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
