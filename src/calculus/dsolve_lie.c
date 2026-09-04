/*
 * dsolve_lie.c — DSolve`LieSymmetry (heuristic Lie point-symmetry method).
 *
 * Integrates a first-order ODE  y' = omega(x, y)  by finding a one-parameter Lie
 * group of point symmetries and reducing the equation to a quadrature.  It is the
 * general first-order backstop in the cascade, run after the deterministic
 * specialists and before the implicit / series fallbacks.
 *
 * For a first-order ODE the linearized symmetry (determining) condition for the
 * generator  V = xi d/dx + eta d/dy  is
 *
 *   S(xi, eta) = eta_x + (eta_y - xi_x) omega - xi_y omega^2
 *                                     - xi omega_x - eta omega_y = 0.        (dag)
 *
 * This one PDE in two unknowns is underdetermined, so — as in SymPy/Maple — a
 * fixed table of ansatze (heuristics) for (xi, eta) is tried; each collapses (dag)
 * to something solvable.  A candidate is accepted only when (dag) is provably zero
 * (the checkinfsol gate, lie_check).  Given an accepted (xi, eta), Lie's theorem
 * gives the integrating factor  mu = 1/(eta - omega xi)  of the 1-form
 * -omega dx + dy = 0, from which a first integral  F(x, y) == C[1]  follows by the
 * exact-equation quadrature; it is returned through dsolve_run_implicit (which
 * verifies by implicit differentiation y' == -F_x/F_y and fits an IVP constant).
 * The first integral is the faithful output of symmetry integration (the source
 * literature returns exactly such implicit quadratures).
 *
 * Staging (see docs/design/dsolve_lie_symmetry.md): L1 ships the substrate and the
 * `abaco1_simple` heuristic (one-variable ansatze); L2 adds `linear` (affine
 * ansatz — the linear-coefficients class, via the determining-system NullSpace),
 * `abaco1_product` (product ansatz [F(x)G(y),0]), `abaco2_similar` ([F(x),H(x)]), and
 * `function_sum` (additive ansatz [F(x)+G(y),0]); L3 adds `bivariate` (general
 * degree-2/3 polynomial ansatz — the same NullSpace determining system at higher
 * degree) and `abaco2_unique_unknown` ([F(x),G(y)]/[G(y),F(x)] from functions or
 * non-integer powers of both variables, incl. the §4.4.1 order-zero candidates
 * [-R,1]/[1,-R]/[1,-1/R]), and `chi` (CPC 101, 5th algorithm: the eta = xi omega + chi
 * reformulation with a rich transcendental-atom basis, reaching genuinely
 * transcendental chi).  The formal §4.4.2 Case I/II (`abaco2_unique_general`) is a
 * documented exemption (impractical per CT-Roche).  All chain in through the same
 * lie_check + lie_first_integral pipeline.
 *
 * Robustness: an omega carrying an undefined function of both variables makes the
 * quadrature classifiers balloon and the general zero_test / Integrate engine hang, so
 * the shared hot-path helpers bail past a node budget (LIE_EXPR_BUDGET), decide zero by
 * a structural polynomial test, and the rational/algebraic ansatze are skipped there.
 *
 * References: Cheb-Terrab & Roche, CPC 113 (1998) 239; Cheb-Terrab, Duarte & da
 * Mota, CPC 101 (1997) 254; Cheb-Terrab & Kolokolnikov, math-ph/0007023.
 */
#include "dsolve_common.h"
#include "integrate.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../expr.h"
#include <stdlib.h>
#include <stdio.h>   /* snprintf */

/* base^n as an evaluated Power; base consumed. */
static Expr* powi(Expr* base, int n) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                         (Expr*[]){ base, expr_new_integer(n) }, 2));
}

/* Node count of e, accumulating into *acc and short-circuiting the recursion the
 * instant *acc exceeds budget — so this is O(budget), not O(size), and stays cheap
 * even on a million-leaf expression.  Borrowed. */
static void lie_count_nodes(const Expr* e, long* acc, long budget) {
    if (!e || *acc > budget) return;
    (*acc)++;
    if (e->type == EXPR_FUNCTION) {
        lie_count_nodes(e->data.function.head, acc, budget);
        for (size_t i = 0; i < e->data.function.arg_count && *acc <= budget; i++)
            lie_count_nodes(e->data.function.args[i], acc, budget);
    }
}

/* True iff e has more than LIE_EXPR_BUDGET nodes.  Every heavy symbolic op on the Lie
 * hot path (Together / Cancel / Expand / the polynomial free-of and zero tests) is
 * super-linear, and on a transcendental omega carrying an UNDEFINED function (e.g.
 * y' == Tan[ArcTan[y] + F[x^2+y^2]]) the symbolic derivatives the quadrature
 * heuristics build balloon — one measured pre-Expand numerator hit 217k nodes and
 * Expand blew it past a million, hanging function_sum before abaco2_unique_* could
 * run — while a genuine rational/algebraic target (the class these heuristics solve)
 * keeps every intermediate to a few hundred nodes.  So the shared helpers below bail
 * the instant an INPUT crosses this budget, before running the expensive op.  This is
 * safe: Lie is the heuristic backstop; a bailed helper only makes its heuristic
 * decline (never a wrong answer — every returned first integral is still
 * back-substitution verified by dsolve_run_implicit), and the cheaper abaco2_unique_*
 * kernels still run.  A deterministic node budget is the machine-independent analogue
 * of the wall-clock timeout SymPy and Maple use on exactly these inputs. */
#define LIE_EXPR_BUDGET 6000L
static bool lie_too_big(const Expr* e) {
    long acc = 0;
    lie_count_nodes(e, &acc, LIE_EXPR_BUDGET);
    return acc > LIE_EXPR_BUDGET;
}

/* True iff e is the literal integer 0.  After Expand, a polynomial that is identically
 * zero in its (possibly transcendental) atoms collapses to this via Plus like-term
 * combination — so this structural test IS the polynomial zero-test the fast helpers
 * below want.  It deliberately replaces the general ds_is_zero / zero_test there: the
 * latter runs Simplify/sampling that HANGS on a transcendental carrying an undefined
 * function (the abaco2_unique_* target class), and its extra power — catching
 * trig-identity zeros Expand misses — is not wanted for a heuristic necessary
 * condition (missing one only makes the heuristic decline; the final lie_check still
 * validates the actual symmetry with the full zero test). */
static bool lie_lit_zero(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* True iff e contains an application of an UNDEFINED function — a user function
 * carrying no derivative / evaluation rule (F[x^2+y^2], H[...]).  Two signatures: a
 * raw application H[...] whose head symbol has neither a builtin nor down-values; and
 * the inert Derivative[n][H][...] node that differentiating such an H leaves behind
 * (a known analytic head never leaves a Derivative — its rule fires).  On such an
 * omega the general zero_test and the Integrate engine hang, and the Lie
 * integrating-factor quadrature is non-elementary anyway, so callers switch to the
 * bounded polynomial zero-test and decline the integration (matching the documented
 * "arbitrary function -> non-elementary -> decline" behavior). */
static bool lie_has_undefined_function(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        if (h->data.symbol.name == SYM_Derivative) return true;   /* inert derivative */
        SymbolDef* d = symtab_lookup(h->data.symbol.name);
        if (d && !d->builtin_func && !d->down_values) return true;
    } else if (lie_has_undefined_function(h)) {
        return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (lie_has_undefined_function(e->data.function.args[i])) return true;
    return false;
}

/* Cancel[Together[e]] — the cheap rational-function normalizer (over a common
 * denominator, then cancel the gcd).  Used in place of the far costlier general
 * Simplify on the heuristics' hot path: it is all these rational (xi, eta) forms
 * need, and it keeps the per-decline cost down to a few ms.  e consumed.  Returns e
 * unchanged when it exceeds the node budget (see lie_too_big) — the downstream
 * consumers are budget-guarded too, so the heuristic simply declines. */
static Expr* lie_ratsimp(Expr* e) {
    if (lie_too_big(e)) return e;
    return eval_and_free(ds_call1("Cancel", eval_and_free(ds_call1(SYM_Together, e))));
}

/* Fast "e is free of var" test, e consumed.  For a rational e = P/Q, e is free of
 * var iff its derivative numerator  P_var Q - P Q_var  is identically zero, which is
 * a POLYNOMIAL zero-test after Expand + lie_lit_zero — far cheaper than the general
 * zero_test on the rational derivative D[e, var] (the dominant per-decline cost of
 * the quadrature heuristics: ~2ms on a moderate rational).  Bails (returns "not free
 * of", the conservative decline) the instant the input, or the pre-Expand numerator,
 * exceeds the node budget; see lie_too_big. */
static bool lie_free_of_var(Expr* e, const char* var) {
    if (lie_too_big(e)) { expr_free(e); return false; }
    Expr* t = eval_and_free(ds_call1(SYM_Together, e));         /* one coprime fraction */
    Expr* P = eval_and_free(ds_call1(SYM_Numerator, expr_copy(t)));
    Expr* Q = eval_and_free(ds_call1(SYM_Denominator, t));      /* consumes t */
    Expr* Pv = ds_d(expr_copy(P), expr_new_symbol(var));
    Expr* Qv = ds_d(expr_copy(Q), expr_new_symbol(var));
    Expr* num = eval_and_free(ds_call2(SYM_Subtract,            /* P_var Q - P Q_var */
                    ds_call2(SYM_Times, Pv, Q),
                    ds_call2(SYM_Times, P, Qv)));               /* consumes Pv,Q,P,Qv */
    if (lie_too_big(num)) { expr_free(num); return false; }
    Expr* ex = eval_and_free(ds_call1("Expand", num));          /* consumes num */
    bool z = lie_lit_zero(ex);
    expr_free(ex);
    return z;
}

/* Fast zero-test, e consumed.  For a rational e, e == 0 iff Expand[Numerator[
 * Together[e]]] collapses to the literal 0 — a polynomial zero-test that sidesteps
 * the general zero_test's sampling/Simplify (which hangs on the transcendental
 * undefined-function derivatives the recursive callers feed the quadrature
 * heuristics).  Bails (returns "not zero") when the input, or the pre-Expand
 * numerator, exceeds the node budget — a numerator too large for Together to have
 * collapsed is a nonzero polynomial in its transcendental atoms, so the conservative
 * answer is also the correct one here; see lie_too_big. */
static bool lie_is_zero(Expr* e) {
    if (lie_too_big(e)) { expr_free(e); return false; }
    Expr* n = eval_and_free(ds_call1(SYM_Numerator, eval_and_free(ds_call1(SYM_Together, e))));
    if (lie_too_big(n)) { expr_free(n); return false; }
    Expr* ex = eval_and_free(ds_call1("Expand", n));
    bool z = lie_lit_zero(ex);
    expr_free(ex);
    return z;
}

/* True iff e is a numeric literal < 0, or a product Times[c, ...] whose leading
 * factor is one.  Used only to sign-normalize an integrand (see below). */
static bool lie_negative_leading(const Expr* e) {
    if (e->type == EXPR_INTEGER) return e->data.integer < 0;
    if (e->type == EXPR_REAL)    return e->data.real < 0.0;
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count > 0 &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Times) {
        const Expr* c = e->data.function.args[0];
        if (c->type == EXPR_INTEGER) return c->data.integer < 0;
        if (c->type == EXPR_REAL)    return c->data.real < 0.0;
    }
    return false;
}

/* Integrate[e, var] working around an Integrate-engine asymmetry: on a NON-
 * elementary integrand a leading negative sign takes a far slower path (measured:
 * Integrate[-1/Sqrt[y^4+C[1]], y] ~ 3.7s vs Integrate[1/Sqrt[...], y] ~ 0.4s and
 * -Integrate[1/Sqrt[...], y] ~ 0.36s).  When e is negative-leading, integrate -e and
 * negate the result — pure linearity, always valid.  e consumed; quiet by caller. */
static Expr* lie_integrate_signfix(Expr* e, const char* var) {
    if (!lie_negative_leading(e)) return ds_integrate(e, expr_new_symbol(var));
    Expr* pe = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), e));   /* -e */
    Expr* r  = ds_integrate(pe, expr_new_symbol(var));                        /* consumes pe */
    if (ds_has_head(r, SYM_Integrate)) return r;   /* still inert: sign is irrelevant */
    return eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), r));       /* -Integrate[-e] */
}

/* ---- the linearized symmetry residual S(xi, eta) of (dag), in the symbols
 *      xv (independent) and Yn (the plain symbol standing for y). xi, eta, omega
 *      are borrowed; result owned (evaluated). ---- */
static Expr* lie_S_expr(const Expr* xi, const Expr* eta, const Expr* omega,
                        const char* xv, const char* Yn) {
    Expr* eta_x = ds_d(expr_copy((Expr*)eta),   expr_new_symbol(xv));
    Expr* eta_y = ds_d(expr_copy((Expr*)eta),   expr_new_symbol(Yn));
    Expr* xi_x  = ds_d(expr_copy((Expr*)xi),    expr_new_symbol(xv));
    Expr* xi_y  = ds_d(expr_copy((Expr*)xi),    expr_new_symbol(Yn));
    Expr* om_x  = ds_d(expr_copy((Expr*)omega), expr_new_symbol(xv));
    Expr* om_y  = ds_d(expr_copy((Expr*)omega), expr_new_symbol(Yn));
    Expr* om2   = powi(expr_copy((Expr*)omega), 2);

    /* eta_x + (eta_y - xi_x) omega - xi_y omega^2 - xi omega_x - eta omega_y */
    Expr* t2 = ds_call2(SYM_Times, ds_call2(SYM_Subtract, eta_y, xi_x),
                        expr_copy((Expr*)omega));
    Expr* t3 = ds_call2(SYM_Times, expr_new_integer(-1), ds_call2(SYM_Times, xi_y, om2));
    Expr* t4 = ds_call2(SYM_Times, expr_new_integer(-1),
                        ds_call2(SYM_Times, expr_copy((Expr*)xi), om_x));
    Expr* t5 = ds_call2(SYM_Times, expr_new_integer(-1),
                        ds_call2(SYM_Times, expr_copy((Expr*)eta), om_y));
    Expr* S = ds_call2(SYM_Plus, eta_x,
                  ds_call2(SYM_Plus, t2,
                      ds_call2(SYM_Plus, t3, ds_call2(SYM_Plus, t4, t5))));
    return eval_and_free(S);
}

/* checkinfsol: is (xi, eta) a genuine symmetry of y' = omega?  Together-normalize
 * the residual, then decide it is zero.  The fast polynomial test (lie_is_zero)
 * settles every rational/algebraic case cheaply; only when that says "not zero" AND
 * the residual is free of undefined functions do we fall back to the full zero_test
 * (which catches identity-based zeros a plain Expand misses — but HANGS on a
 * transcendental carrying an undefined function, the abaco2_unique_* target class, so
 * it is skipped there).  Missing a symmetry only makes the heuristic decline; a
 * false accept is caught downstream by dsolve_run_implicit's back-substitution. */
static bool lie_check(const Expr* xi, const Expr* eta, const Expr* omega,
                      const char* xv, const char* Yn) {
    Expr* S  = lie_S_expr(xi, eta, omega, xv, Yn);
    Expr* St = eval_and_free(ds_call1("Together", S));   /* consumes S */
    if (lie_too_big(St)) { expr_free(St); return false; }
    bool z;
    if (lie_is_zero(expr_copy(St)))          z = true;            /* fast polynomial 0 */
    else if (lie_has_undefined_function(St)) z = false;          /* full test would hang */
    else                                     z = ds_is_zero(St);  /* identity-based 0 */
    expr_free(St);
    return z;
}

/* Given a validated symmetry (xi, eta), build the first integral G(x, y[x]) of
 * y' = omega via the Lie integrating factor mu = 1/(eta - xi omega):
 *   P = -mu omega (coeff of dx), Q = mu (coeff of dy),
 *   Phi = Integrate[P, x],   G = Phi + Integrate[Q - d(Phi)/dY, Y]  with Y -> y[x].
 * Returns G (meaning G == C[1]) or NULL if the symmetry is tangent (eta == xi
 * omega), a quadrature is non-elementary, the exactness residual retains x, or G
 * does not actually depend on y.  xi, eta, omega borrowed. */
static Expr* lie_first_integral(const Expr* xi, const Expr* eta, const Expr* omega,
                                const char* xv, const char* Yn, const char* yname) {
    /* An omega carrying an undefined function makes the Lie integrating-factor
     * quadrature non-elementary and hangs the Integrate engine; decline up front
     * (matching the documented "arbitrary function -> non-elementary -> decline"). */
    if (lie_has_undefined_function(omega)) return NULL;
    Expr* den = eval_and_free(ds_call2(SYM_Subtract, expr_copy((Expr*)eta),
                    ds_call2(SYM_Times, expr_copy((Expr*)xi), expr_copy((Expr*)omega))));
    if (ds_is_zero(den)) { expr_free(den); return NULL; }
    Expr* mu = powi(den, -1);                                     /* consumes den */

    /* P = -mu omega */
    Expr* Pc = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
                   ds_call2(SYM_Times, expr_copy(mu), expr_copy((Expr*)omega))));
    g_integrate_quiet++;
    Expr* Phi = lie_integrate_signfix(Pc, xv);                   /* consumes Pc */
    g_integrate_quiet--;
    if (ds_has_head(Phi, SYM_Integrate)) { expr_free(Phi); expr_free(mu); return NULL; }
    /* corr = Q - d(Phi)/dY = mu - d(Phi)/dY  (exactness => free of x) */
    Expr* dPhiY = ds_d(expr_copy(Phi), expr_new_symbol(Yn));
    Expr* corr = ds_simplify(eval_and_free(ds_call2(SYM_Subtract, expr_copy(mu), dPhiY)));
    expr_free(mu);
    if (!ds_free_of(corr, xv)) { expr_free(corr); expr_free(Phi); return NULL; }
    g_integrate_quiet++;
    Expr* Fy = lie_integrate_signfix(corr, Yn);                  /* consumes corr */
    g_integrate_quiet--;
    if (ds_has_head(Fy, SYM_Integrate)) { expr_free(Fy); expr_free(Phi); return NULL; }

    Expr* F = eval_and_free(ds_call2(SYM_Plus, Phi, Fy));        /* consumes Phi, Fy */
    Expr* G = ds_subst(F, expr_new_symbol(Yn), ds_make_funcapp(yname, 0, xv)); /* consumes F */
    if (!ds_contains(G, yname)) { expr_free(G); return NULL; }   /* trivial (x only) */
    return G;
}

/* Gate the candidate (xi, eta) through lie_check and, if it passes, integrate it to
 * the first integral.  Consumes xi and eta; returns G (owned) or NULL.  The shared
 * accept-and-integrate step for every heuristic that produces explicit (xi, eta). */
static Expr* lie_try_symmetry(Expr* xi, Expr* eta, const Expr* omega,
                              const char* xv, const char* Yn, const char* yname) {
    Expr* G = NULL;
    if (lie_check(xi, eta, omega, xv, Yn))
        G = lie_first_integral(xi, eta, omega, xv, Yn, yname);
    expr_free(xi); expr_free(eta);
    return G;
}

/* Exp[Integrate[R, var]], or NULL if that integral is non-elementary.  R consumed. */
static Expr* lie_exp_integral(Expr* R, const char* var) {
    g_integrate_quiet++;
    Expr* I = lie_integrate_signfix(R, var);                     /* consumes R */
    g_integrate_quiet--;
    if (ds_has_head(I, SYM_Integrate)) { expr_free(I); return NULL; }
    return eval_and_free(ds_call1("Exp", I));                    /* consumes I */
}

/* The product-separable x-factor of L: Exp[Integrate[L_x/L, x]].  Returns NULL when
 * L == 0, when L_x/L is not free of y (so L does not separate as X(x) Y(y): note
 * L_x/L free of y  <=>  d^2/dx dy log L = 0  <=>  L separable), or when the integral
 * is non-elementary.  Used by the quadrature heuristics to extract F(x) / 1/F_xx
 * from the necessary-condition expression.  L borrowed; result owned. */
static Expr* lie_sep_xfactor(const Expr* L, const char* xv, const char* Yn) {
    if (lie_too_big(L)) return NULL;                   /* transcendental blowup: decline */
    if (lie_is_zero(expr_copy((Expr*)L))) return NULL;
    Expr* Lx = ds_d(expr_copy((Expr*)L), expr_new_symbol(xv));
    /* ratio = Together[L_x / L].  Together (value-preserving, so it cannot change
     * the free-of answer) collapses the L^4-denominator mess into one fraction. */
    if (lie_too_big(Lx)) { expr_free(Lx); return NULL; }
    Expr* ratio = eval_and_free(ds_call1(SYM_Together,
                      eval_and_free(ds_call2(SYM_Times, Lx, powi(expr_copy((Expr*)L), -1)))));
    if (!lie_free_of_var(expr_copy(ratio), Yn)) { expr_free(ratio); return NULL; }
    return lie_exp_integral(lie_ratsimp(ratio), xv);           /* consumes ratio */
}

/* Simultaneously swap the symbols xv and Yn throughout e (via a temporary), so an
 * expression written in (x, y) becomes the same expression in (y, x).  e consumed. */
static Expr* lie_swap_xy(Expr* e, const char* xv, const char* Yn) {
    const char* tmp = intern_symbol("DSolve`lieSwap");
    e = ds_subst(e, expr_new_symbol(xv),  expr_new_symbol(tmp));
    e = ds_subst(e, expr_new_symbol(Yn),  expr_new_symbol(xv));
    e = ds_subst(e, expr_new_symbol(tmp), expr_new_symbol(Yn));
    return e;
}

/* The RHS of the inverse ODE (Cheb-Terrab&Roche Def 1): y' = 1/omega(y, x), i.e.
 * omega with x<->y swapped, then reciprocated.  A [xi,eta] symmetry of the inverse
 * ODE maps back (Def 3) to [swap(eta), swap(xi)] of the original, letting one
 * pattern-A extractor also cover its inverse pattern.  omega borrowed; result owned. */
static Expr* lie_inverse_omega(const Expr* omega, const char* xv, const char* Yn) {
    Expr* sw = lie_swap_xy(expr_copy((Expr*)omega), xv, Yn);
    /* Together so the reciprocal is one clean fraction: its derivatives (taken by
     * the extractor) are then as cheap as the direct branch's, not 5x costlier. */
    return eval_and_free(ds_call1(SYM_Together, powi(sw, -1)));  /* consumes sw */
}

/* Heuristic `abaco1_simple`: one infinitesimal zero, the other a function of a
 * single variable.  The three productive sub-ansatze (the fourth, xi=g(y) eta=0,
 * forces omega==0) reduce (dag) to "some ratio is free of the other variable":
 *   A. xi=0, eta=g(x):   omega_y free of y            => g = Exp[Integrate[omega_y, x]]
 *   B. xi=g(x), eta=0:   omega_x/omega free of y      => g = Exp[-Integrate[omega_x/omega, x]]
 *   C. xi=0, eta=g(y):   omega_y/omega free of x      => g = Exp[Integrate[omega_y/omega, y]]
 * Each candidate is gated through lie_check before integration.  Returns the first
 * first integral found, or NULL.  omega, wx, wy borrowed. */
static Expr* lie_abaco1_simple(const Expr* omega, const Expr* wx, const Expr* wy,
                               const char* xv, const char* Yn, const char* yname) {
    Expr* G = NULL;
    Expr* zero = expr_new_integer(0);

    /* A: xi=0, eta=g(x) */
    if (ds_free_of(wy, Yn)) {
        Expr* g = lie_exp_integral(expr_copy((Expr*)wy), xv);
        if (g) {
            if (lie_check(zero, g, omega, xv, Yn))
                G = lie_first_integral(zero, g, omega, xv, Yn, yname);
            expr_free(g);
        }
    }

    /* B: xi=g(x), eta=0 */
    if (!G && !ds_is_zero(omega)) {
        Expr* ratio = lie_ratsimp(eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)wx),
                          powi(expr_copy((Expr*)omega), -1))));
        if (ds_free_of(ratio, Yn)) {
            Expr* neg = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(ratio)));
            Expr* g = lie_exp_integral(neg, xv);                 /* consumes neg */
            if (g) {
                bool _ck = lie_check(g, zero, omega, xv, Yn);
                if (_ck)
                    G = lie_first_integral(g, zero, omega, xv, Yn, yname);
                expr_free(g);
            }
        }
        expr_free(ratio);
    }

    /* C: xi=0, eta=g(y) */
    if (!G && !ds_is_zero(omega)) {
        Expr* ratio = lie_ratsimp(eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)wy),
                          powi(expr_copy((Expr*)omega), -1))));
        if (ds_free_of(ratio, xv)) {
            Expr* g = lie_exp_integral(expr_copy(ratio), Yn);
            if (g) {
                if (lie_check(zero, g, omega, xv, Yn))
                    G = lie_first_integral(zero, g, omega, xv, Yn, yname);
                expr_free(g);
            }
        }
        expr_free(ratio);
    }

    expr_free(zero);
    return G;
}

/* {x, y} as a List, for CoefficientList / PolynomialQ over both variables. */
static Expr* lie_xy_list(const char* xv, const char* Yn) {
    return expr_new_function(expr_new_symbol(SYM_List),
        (Expr*[]){ expr_new_symbol(xv), expr_new_symbol(Yn) }, 2);
}

/* Evaluate e and test whether it is the symbol True; e consumed. */
static bool lie_ev_trueQ(Expr* e) {
    Expr* v = eval_and_free(e);
    Expr* T = expr_new_symbol("True");
    bool t = expr_eq(v, T);
    expr_free(v); expr_free(T);
    return t;
}

/* c1*x + c2*y + c3, evaluated; c1, c2, c3 consumed. */
static Expr* lie_affine(Expr* c1, Expr* c2, Expr* c3, const char* xv, const char* Yn) {
    return eval_and_free(ds_call2(SYM_Plus,
        ds_call2(SYM_Times, c1, expr_new_symbol(xv)),
        ds_call2(SYM_Plus,
            ds_call2(SYM_Times, c2, expr_new_symbol(Yn)),
            c3)));
}

/* Heuristic `linear`: an affine symmetry xi = a1 x + a2 y + a3,
 * eta = a4 x + a5 y + a6.  For rational omega the symmetry condition (dag), after
 * clearing denominators, is a polynomial identity in x, y whose coefficients are
 * linear and homogeneous in the six constants; the determining system is the
 * coefficient matrix, and its NullSpace is a basis of admissible (xi, eta) — no
 * Solve (hence no free-variable message).  Each basis symmetry is gated by
 * lie_check and handed to lie_first_integral.  Captures the linear-coefficients
 * class y' = (a1 x + b1 y + c1)/(a2 x + b2 y + c2) (an affine scaling about the
 * lines' intersection) that the earlier deterministic methods miss.  omega
 * borrowed; returns the first first integral found, or NULL. */
static Expr* lie_linear(const Expr* omega, const char* xv, const char* Yn, const char* yname) {
    /* rationality guard: omega must be a ratio of polynomials in x, y */
    Expr* tg   = eval_and_free(ds_call1(SYM_Together, expr_copy((Expr*)omega)));
    Expr* num0 = eval_and_free(ds_call1(SYM_Numerator, expr_copy(tg)));
    Expr* den0 = eval_and_free(ds_call1(SYM_Denominator, tg));         /* consumes tg */
    bool r1 = lie_ev_trueQ(eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
                  (Expr*[]){ num0, lie_xy_list(xv,Yn) }, 2)));         /* consumes num0 */
    bool r2 = lie_ev_trueQ(eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
                  (Expr*[]){ den0, lie_xy_list(xv,Yn) }, 2)));         /* consumes den0 */
    if (!(r1 && r2)) return NULL;

    const char* a[6] = {
        intern_symbol("DSolve`lieA1"), intern_symbol("DSolve`lieA2"),
        intern_symbol("DSolve`lieA3"), intern_symbol("DSolve`lieA4"),
        intern_symbol("DSolve`lieA5"), intern_symbol("DSolve`lieA6") };
    Expr* xi  = lie_affine(expr_new_symbol(a[0]), expr_new_symbol(a[1]), expr_new_symbol(a[2]), xv, Yn);
    Expr* eta = lie_affine(expr_new_symbol(a[3]), expr_new_symbol(a[4]), expr_new_symbol(a[5]), xv, Yn);
    Expr* S = lie_S_expr(xi, eta, omega, xv, Yn);                      /* borrows xi, eta */
    expr_free(xi); expr_free(eta);

    /* forms = Flatten[CoefficientList[Numerator[Together[S]], {x, y}]] */
    Expr* tS = eval_and_free(ds_call1(SYM_Together, S));              /* consumes S */
    Expr* nS = eval_and_free(ds_call1(SYM_Numerator, tS));           /* consumes tS */
    Expr* cl = eval_and_free(expr_new_function(expr_new_symbol(SYM_CoefficientList),
                   (Expr*[]){ nS, lie_xy_list(xv,Yn) }, 2));         /* consumes nS */
    Expr* forms = eval_and_free(ds_call1(SYM_Flatten, cl));          /* consumes cl */

    /* M = Outer[Coefficient, forms, {a1..a6}]; ns = NullSpace[M] */
    Expr* vars = expr_new_function(expr_new_symbol(SYM_List),
        (Expr*[]){ expr_new_symbol(a[0]),expr_new_symbol(a[1]),expr_new_symbol(a[2]),
                   expr_new_symbol(a[3]),expr_new_symbol(a[4]),expr_new_symbol(a[5]) }, 6);
    Expr* M = eval_and_free(expr_new_function(expr_new_symbol("Outer"),
                  (Expr*[]){ expr_new_symbol("Coefficient"), forms, vars }, 3)); /* consumes forms, vars */
    if (!ds_free_of(M, xv) || !ds_free_of(M, Yn)) { expr_free(M); return NULL; }  /* not a clean system */
    Expr* ns = eval_and_free(ds_call1("NullSpace", M));              /* consumes M */

    Expr* G = NULL;
    if (ns && ns->type == EXPR_FUNCTION) {
        size_t nv = ns->data.function.arg_count;
        for (size_t k = 0; k < nv && !G; k++) {
            Expr* v = ns->data.function.args[k];
            if (!v || v->type != EXPR_FUNCTION || v->data.function.arg_count != 6) continue;
            Expr** c = v->data.function.args;
            Expr* xiv  = lie_affine(expr_copy(c[0]), expr_copy(c[1]), expr_copy(c[2]), xv, Yn);
            Expr* etav = lie_affine(expr_copy(c[3]), expr_copy(c[4]), expr_copy(c[5]), xv, Yn);
            if (lie_check(xiv, etav, omega, xv, Yn))
                G = lie_first_integral(xiv, etav, omega, xv, Yn, yname);
            expr_free(xiv); expr_free(etav);
        }
    }
    expr_free(ns);
    return G;
}

/* Interned coefficient symbol "DSolve`lieB<k>" (k >= 1) for the bivariate
 * polynomial ansatz — a namespace distinct from lie_linear's lieA*.  The interner
 * copies the name, so returning a pointer into the reused static buffer is safe. */
static const char* lie_coeff_name(int k) {
    char buf[32];
    snprintf(buf, sizeof buf, "DSolve`lieB%d", k);
    return intern_symbol(buf);
}

/* Build the general bivariate polynomial  Sum_{i+j<=degree} coeffs[k] x^i y^j  in a
 * fixed monomial order (total degree ascending; within a degree, x-power
 * descending), consuming the first nmon = (degree+1)(degree+2)/2 entries of
 * coeffs.  coeffs[] are borrowed (copied in); result owned and evaluated. */
static Expr* lie_poly_build(Expr** coeffs, int degree, const char* xv, const char* Yn) {
    Expr* sum = expr_new_integer(0);
    int k = 0;
    for (int total = 0; total <= degree; total++) {
        for (int i = total; i >= 0; i--) {
            int j = total - i;
            Expr* term = expr_copy(coeffs[k++]);
            if (i > 0) term = ds_call2(SYM_Times, term, powi(expr_new_symbol(xv), i));
            if (j > 0) term = ds_call2(SYM_Times, term, powi(expr_new_symbol(Yn), j));
            sum = ds_call2(SYM_Plus, sum, term);
        }
    }
    return eval_and_free(sum);
}

/* Heuristic `bivariate` at a fixed degree: xi, eta are general bivariate
 * polynomials of total degree <= `degree` in x, y (the affine `linear` heuristic
 * is the degree-1 case).  For rational omega the symmetry condition (dag) clears
 * to a polynomial identity whose coefficients are linear and homogeneous in the
 * unknown polynomial coefficients; the determining system's NullSpace is a basis
 * of admissible (xi, eta).  Each basis symmetry is gated by lie_check and handed
 * to lie_first_integral.  Catches genuinely quadratic / projective symmetries
 * (xi=x^2, eta=x y; the projective group) that the affine ansatz misses.  omega
 * borrowed; returns the first first integral found, or NULL. */
static Expr* lie_poly_symmetry(const Expr* omega, int degree,
                               const char* xv, const char* Yn, const char* yname) {
    /* rationality guard: omega must be a ratio of polynomials in x, y */
    Expr* tg   = eval_and_free(ds_call1(SYM_Together, expr_copy((Expr*)omega)));
    Expr* num0 = eval_and_free(ds_call1(SYM_Numerator, expr_copy(tg)));
    Expr* den0 = eval_and_free(ds_call1(SYM_Denominator, tg));         /* consumes tg */
    bool r1 = lie_ev_trueQ(eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
                  (Expr*[]){ num0, lie_xy_list(xv,Yn) }, 2)));         /* consumes num0 */
    bool r2 = lie_ev_trueQ(eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
                  (Expr*[]){ den0, lie_xy_list(xv,Yn) }, 2)));         /* consumes den0 */
    if (!(r1 && r2)) return NULL;

    int nmon   = (degree + 1) * (degree + 2) / 2;
    int ncoeff = 2 * nmon;                                            /* xi + eta */
    Expr** csym = malloc((size_t)ncoeff * sizeof(Expr*));
    for (int k = 0; k < ncoeff; k++) csym[k] = expr_new_symbol(lie_coeff_name(k + 1));

    Expr* xi  = lie_poly_build(&csym[0],    degree, xv, Yn);
    Expr* eta = lie_poly_build(&csym[nmon], degree, xv, Yn);
    Expr* S = lie_S_expr(xi, eta, omega, xv, Yn);                      /* borrows xi, eta */
    expr_free(xi); expr_free(eta);

    /* forms = Flatten[CoefficientList[Numerator[Together[S]], {x, y}]] */
    Expr* tS = eval_and_free(ds_call1(SYM_Together, S));              /* consumes S */
    Expr* nS = eval_and_free(ds_call1(SYM_Numerator, tS));           /* consumes tS */
    Expr* cl = eval_and_free(expr_new_function(expr_new_symbol(SYM_CoefficientList),
                   (Expr*[]){ nS, lie_xy_list(xv,Yn) }, 2));         /* consumes nS */
    Expr* forms = eval_and_free(ds_call1(SYM_Flatten, cl));          /* consumes cl */

    /* M = Outer[Coefficient, forms, {c1..c_ncoeff}]; ns = NullSpace[M] */
    Expr** vargs = malloc((size_t)ncoeff * sizeof(Expr*));
    for (int k = 0; k < ncoeff; k++) vargs[k] = expr_copy(csym[k]);
    Expr* vars = expr_new_function(expr_new_symbol(SYM_List), vargs, (size_t)ncoeff);
    free(vargs);                                                      /* elements adopted */
    for (int k = 0; k < ncoeff; k++) expr_free(csym[k]);
    free(csym);

    Expr* M = eval_and_free(expr_new_function(expr_new_symbol("Outer"),
                  (Expr*[]){ expr_new_symbol("Coefficient"), forms, vars }, 3)); /* consumes forms, vars */
    if (!ds_free_of(M, xv) || !ds_free_of(M, Yn)) { expr_free(M); return NULL; }  /* not a clean system */
    Expr* ns = eval_and_free(ds_call1("NullSpace", M));              /* consumes M */

    Expr* G = NULL;
    if (ns && ns->type == EXPR_FUNCTION) {
        size_t nv = ns->data.function.arg_count;
        for (size_t r = 0; r < nv && !G; r++) {
            Expr* v = ns->data.function.args[r];
            if (!v || v->type != EXPR_FUNCTION || v->data.function.arg_count != (size_t)ncoeff)
                continue;
            Expr** c = v->data.function.args;                        /* borrowed */
            Expr* xiv  = lie_poly_build(&c[0],    degree, xv, Yn);
            Expr* etav = lie_poly_build(&c[nmon], degree, xv, Yn);
            if (lie_check(xiv, etav, omega, xv, Yn))
                G = lie_first_integral(xiv, etav, omega, xv, Yn, yname);
            expr_free(xiv); expr_free(etav);
        }
    }
    expr_free(ns);
    return G;
}

/* Heuristic `bivariate`: try degree-2 then degree-3 polynomial symmetries (the
 * degree-1 case is the earlier `linear` heuristic, which has already declined). */
static Expr* lie_bivariate(const Expr* omega, const char* xv, const char* Yn, const char* yname) {
    Expr* G = lie_poly_symmetry(omega, 2, xv, Yn, yname);
    if (!G) G = lie_poly_symmetry(omega, 3, xv, Yn, yname);
    return G;
}

/* abaco1_product candidate [xi = F(x) G(y), eta = 0] (Cheb-Terrab & Roche §4.1).
 * The necessary condition (their Eq 19) is that
 *   L = (omega_xy omega - omega_x omega_y) / omega^4
 * be product-separable in x, y; then F is its x-factor (Prop 1) and
 *   g(y) = F d/dx( 1/(F omega) )   (Eq 20)
 * must be free of x, giving G = Exp[Integrate[g, y]].  Fills *xi = F G, *eta = 0 and
 * returns true on success; on any failure sets both outputs NULL and returns false.
 * omega borrowed; the two outputs are owned by the caller on success. */
static bool lie_abaco1_product_cand(const Expr* omega, const char* xv, const char* Yn,
                                    Expr** xi, Expr** eta) {
    *xi = NULL; *eta = NULL;
    if (ds_is_zero(omega)) return false;

    Expr* omx  = ds_d(expr_copy((Expr*)omega), expr_new_symbol(xv));
    Expr* omy  = ds_d(expr_copy((Expr*)omega), expr_new_symbol(Yn));
    Expr* omxy = ds_d(expr_copy(omx), expr_new_symbol(Yn));       /* copy: omx reused below */
    /* L = (omega_xy omega - omega_x omega_y) / omega^4 */
    Expr* num = eval_and_free(ds_call2(SYM_Subtract,
                    ds_call2(SYM_Times, omxy, expr_copy((Expr*)omega)),
                    ds_call2(SYM_Times, omx, omy)));             /* consumes omxy, omx, omy */
    /* L = num / omega^4, left UNSIMPLIFIED (lie_sep_xfactor's free-of gate does not
     * need it; Simplify here would be paid on every declining omega). */
    Expr* L = eval_and_free(ds_call2(SYM_Times, num,
                  powi(expr_copy((Expr*)omega), -4)));           /* consumes num */
    Expr* F = lie_sep_xfactor(L, xv, Yn);
    expr_free(L);
    if (!F) return false;

    /* g(y) = F d/dx( 1/(F omega) ) — must be free of x */
    Expr* Finv = powi(eval_and_free(ds_call2(SYM_Times, expr_copy(F),
                     expr_copy((Expr*)omega))), -1);             /* 1/(F omega) */
    Expr* dinner = ds_d(Finv, expr_new_symbol(xv));              /* consumes Finv */
    Expr* g = lie_ratsimp(eval_and_free(ds_call2(SYM_Times, expr_copy(F), dinner)));
    if (!ds_free_of(g, xv)) { expr_free(g); expr_free(F); return false; }
    Expr* G = lie_exp_integral(g, Yn);                           /* consumes g */
    if (!G) { expr_free(F); return false; }

    *xi  = lie_ratsimp(eval_and_free(ds_call2(SYM_Times, F, G))); /* consumes F, G */
    *eta = expr_new_integer(0);
    return true;
}

/* A pattern-A candidate extractor: fill *xi, *eta (owned) and return true, or set
 * both NULL and return false.  omega borrowed. */
typedef bool (*lie_cand_fn)(const Expr* omega, const char* xv, const char* Yn,
                            Expr** xi, Expr** eta);

/* Run a pattern-A extractor on omega (the direct pattern) and, failing that, on the
 * inverse ODE 1/omega(y,x), mapping that symmetry back to the original as
 * [swap(eta), swap(xi)] (Def 3) so one extractor covers a pattern AND its inverse.
 * Every candidate is gated by lie_check and integrated by lie_first_integral against
 * the ORIGINAL omega.  Returns the first first integral found, or NULL. */
static Expr* lie_run_with_inverse(lie_cand_fn cand, const Expr* omega,
                                  const char* xv, const char* Yn, const char* yname) {
    Expr *xi = NULL, *eta = NULL, *G = NULL;
    if (cand(omega, xv, Yn, &xi, &eta)) {
        if (lie_check(xi, eta, omega, xv, Yn))
            G = lie_first_integral(xi, eta, omega, xv, Yn, yname);
        expr_free(xi); expr_free(eta);
    }
    if (G) return G;

    Expr* inv = lie_inverse_omega(omega, xv, Yn);
    if (cand(inv, xv, Yn, &xi, &eta)) {
        Expr* xio  = lie_swap_xy(eta, xv, Yn);
        Expr* etao = lie_swap_xy(xi,  xv, Yn);
        if (lie_check(xio, etao, omega, xv, Yn))
            G = lie_first_integral(xio, etao, omega, xv, Yn, yname);
        expr_free(xio); expr_free(etao);
    }
    expr_free(inv);
    return G;
}

/* Heuristic `abaco1_product` (§4.1): the symmetry [F(x) G(y), 0], or its inverse
 * [0, F(x) G(y)].  Reaches rational-but-non-polynomial symmetries (xi = y/x, ...)
 * that the polynomial `linear`/`bivariate` ansatze miss.  omega borrowed. */
static Expr* lie_abaco1_product(const Expr* omega, const char* xv, const char* Yn,
                                const char* yname) {
    return lie_run_with_inverse(lie_abaco1_product_cand, omega, xv, Yn, yname);
}

/* function_sum candidate [xi = F(x) + G(y), eta = 0] (Cheb-Terrab & Roche §4.2).
 * The sum-analogue of abaco1_product.  With the "factor" F = omega . d^2/dx^2(1/omega),
 * Eq (27) gives the *rational* identity
 *   F = F''(x) / (F(x) + G(y))                          (the leading omega cancels the
 *                                                        transcendental part of 1/omega),
 * so d/dy(1/F) = G'(y)/F''(x) is product-separable in x, y (Eq 28); its x-factor is
 * 1/F''(x), and then xi = F(x)+G(y) = F''/F = 1/(xfactor . F).  When F == 0 the ODE is
 * invert-linear (Eq 27, footnote 5) and this pattern declines.  Fills *xi, *eta (owned)
 * and returns true, else sets both NULL and returns false.  omega borrowed. */
static bool lie_function_sum_cand(const Expr* omega, const char* xv, const char* Yn,
                                  Expr** xi, Expr** eta) {
    *xi = NULL; *eta = NULL;
    if (ds_is_zero(omega)) return false;

    /* factor = omega . d^2/dx^2(1/omega) = F''/(F+G), rational (Eq 27) */
    Expr* d2 = ds_d(ds_d(powi(expr_copy((Expr*)omega), -1), expr_new_symbol(xv)),
                    expr_new_symbol(xv));
    Expr* factor = lie_ratsimp(eval_and_free(
                       ds_call2(SYM_Times, expr_copy((Expr*)omega), d2)));
    if (lie_is_zero(expr_copy(factor))) { expr_free(factor); return false; } /* invert-linear */

    /* x-factor of d/dy(1/factor) is 1/F''(x) (up to a constant) */
    Expr* dyM  = ds_d(powi(expr_copy(factor), -1), expr_new_symbol(Yn));
    Expr* xfac = lie_sep_xfactor(dyM, xv, Yn);
    expr_free(dyM);
    if (!xfac) { expr_free(factor); return false; }

    /* xi = F + G = F''/factor = 1/(xfac . factor) — the constant scale is irrelevant
     * (lie_check and lie_first_integral are homogeneous in the generator). */
    *xi  = lie_ratsimp(powi(eval_and_free(ds_call2(SYM_Times, xfac, factor)), -1));
    *eta = expr_new_integer(0);
    return true;
}

/* Heuristic `function_sum` (§4.2): the symmetry [F(x) + G(y), 0], or its inverse
 * [0, F(x) + G(y)] (via the inverse ODE).  The additive counterpart of
 * `abaco1_product`.  omega borrowed. */
static Expr* lie_function_sum(const Expr* omega, const char* xv, const char* Yn,
                              const char* yname) {
    return lie_run_with_inverse(lie_function_sum_cand, omega, xv, Yn, yname);
}

/* abaco2_similar candidate [xi = F(x), eta = H(x)] (both single-variable functions
 * of x) — Cheb-Terrab & Roche §4.3.  Q = omega_y / omega_yy (Eq 39); in the Q_y != 0
 * branch T = Q_x / Q_y = -H/F (Eq 40, free of y), then F = Exp[Integrate[(T omega_y
 * - T_x - omega_x)/(omega + T), x]] (Eq 43) provided the integrand is free of y (Eq
 * 44) and H = -T F.  omega borrowed. */
static bool lie_abaco2_similar_cand(const Expr* omega, const char* xv, const char* Yn,
                                    Expr** xi, Expr** eta) {
    *xi = NULL; *eta = NULL;
    Expr* omy  = ds_d(expr_copy((Expr*)omega), expr_new_symbol(Yn));
    Expr* omyy = ds_d(expr_copy(omy), expr_new_symbol(Yn));
    if (lie_is_zero(expr_copy(omyy))) { expr_free(omy); expr_free(omyy); return false; }  /* linear in y */
    Expr* Q  = eval_and_free(ds_call2(SYM_Times, omy, powi(omyy, -1)));
    Expr* Qy = ds_d(expr_copy(Q), expr_new_symbol(Yn));
    if (lie_is_zero(expr_copy(Qy))) { expr_free(Qy); expr_free(Q); return false; }  /* Q_y == 0: future */
    Expr* Qx = ds_d(expr_copy(Q), expr_new_symbol(xv));
    expr_free(Q);
    Expr* T = eval_and_free(ds_call2(SYM_Times, Qx, powi(Qy, -1)));           /* Q_x/Q_y */
    if (!lie_free_of_var(expr_copy(T), Yn)) { expr_free(T); return false; }
    /* T == 0 is the autonomous [F(x), 0] case abaco1_simple / Separable own; declining
     * it avoids re-attempting their (possibly elliptic) quadrature. */
    if (lie_is_zero(expr_copy(T))) { expr_free(T); return false; }
    T = lie_ratsimp(T);

    Expr* omx  = ds_d(expr_copy((Expr*)omega), expr_new_symbol(xv));
    Expr* omy2 = ds_d(expr_copy((Expr*)omega), expr_new_symbol(Yn));
    Expr* Tx   = ds_d(expr_copy(T), expr_new_symbol(xv));
    Expr* numi = eval_and_free(ds_call2(SYM_Subtract,
                     ds_call2(SYM_Subtract, ds_call2(SYM_Times, expr_copy(T), omy2), Tx),
                     omx));
    Expr* deni = eval_and_free(ds_call2(SYM_Plus, expr_copy((Expr*)omega), expr_copy(T)));
    Expr* integ = lie_ratsimp(eval_and_free(ds_call2(SYM_Times, numi, powi(deni, -1))));
    if (!lie_free_of_var(expr_copy(integ), Yn)) { expr_free(integ); expr_free(T); return false; }
    Expr* Ff = lie_exp_integral(integ, xv);
    if (!Ff) { expr_free(T); return false; }

    *xi  = Ff;
    *eta = lie_ratsimp(eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
               ds_call2(SYM_Times, T, expr_copy(Ff)))));
    return true;
}

/* Heuristic `abaco2_similar` (§4.3): [F(x), H(x)] and its inverse [F(y), H(y)]. */
static Expr* lie_abaco2_similar(const Expr* omega, const char* xv, const char* Yn,
                                const char* yname) {
    return lie_run_with_inverse(lie_abaco2_similar_cand, omega, xv, Yn, yname);
}

/* e depends on BOTH xv and Yn. */
static bool lie_has_both(const Expr* e, const char* xv, const char* Yn) {
    return ds_contains((Expr*)e, xv) && ds_contains((Expr*)e, Yn);
}

/* Recursively collect into *out the distinct subexpressions of e that are
 * "functions of both variables": a non-arithmetic function application (Sin, Log,
 * ArcTan, an undefined head, ...), or a power with a non-integer / variable
 * exponent (u^(1/2), a^(x+y), Exp[x+y]), each containing both xv and Yn.  These are
 * the mappings M of Cheb-Terrab & Roche §4.4.1 (Prop 7): any such M of an ODE with
 * a [F(x),G(y)] / [G(y),F(x)] symmetry must be a function of f(x)+g(y).  Rational
 * omega yields none → instant decline.  Deduplicated by expr_eq; *out grown by
 * realloc, elements are owned copies (caller frees). */
static void lie_collect_kernels(const Expr* e, const char* xv, const char* Yn,
                                Expr*** out, size_t* n, size_t* cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    const Expr* head = e->data.function.head;
    bool sym_head = head->type == EXPR_SYMBOL;
    const char* hn = sym_head ? head->data.symbol.name : NULL;
    bool arith = sym_head && (hn == SYM_Plus || hn == SYM_Times);
    bool ispow = sym_head && hn == SYM_Power;

    bool kernel = false;
    if (sym_head && !arith && !ispow) {
        kernel = lie_has_both(e, xv, Yn);                 /* function application */
    } else if (ispow && e->data.function.arg_count == 2) {
        const Expr* ex = e->data.function.args[1];
        if (ex->type != EXPR_INTEGER && lie_has_both(e, xv, Yn)) kernel = true;
    }
    if (kernel)
        for (size_t i = 0; i < *n; i++)
            if (expr_eq((*out)[i], (Expr*)e)) { kernel = false; break; }
    if (kernel) {
        if (*n == *cap) { *cap = *cap ? *cap * 2 : 4;
                          *out = realloc(*out, *cap * sizeof(Expr*)); }
        (*out)[(*n)++] = expr_copy((Expr*)e);
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        lie_collect_kernels(e->data.function.args[i], xv, Yn, out, n, cap);
}

/* -R/X, evaluated (X may be NULL; caller guards).  R, X borrowed; result owned. */
static Expr* lie_neg_R_over_X(const Expr* R, const Expr* X) {
    return lie_ratsimp(eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
               ds_call2(SYM_Times, expr_copy((Expr*)R), powi(expr_copy((Expr*)X), -1)))));
}

/* Heuristic `abaco2_unique_unknown` (§4.4.1): symmetries [F(x), G(y)] and
 * [G(y), F(x)], found from the functions / non-integer powers of both variables
 * present in omega.  For each such mapping M, R = M_y/M_x (Eq 63):
 *   - when R separates by product with x-factor X (= 1/f_x), scheme (iib) gives the
 *     two candidates  [X, -X/R]  (type [F(x),G(y)])  and  [-R/X, 1/X]  (type
 *     [G(y),F(x)]);
 *   - §4.4.1 "more symmetries from the differential invariant of order zero"
 *     (Eqs 73-81) additionally admits, WITHOUT a separability test, the candidates
 *     [-R, 1] (Eq 75, pattern [f(x)g(y),1]) and [1, -R] / [1, -1/R] (the family-77
 *     patterns, e.g. Kamke ODE 433 -> [1, -x/(2x+y)]).
 * Each candidate is gated by lie_check and integrated, first success wins.  Reaches
 * ODEs carrying an arbitrary function or non-integer power of a combined argument,
 * which every rational ansatz structurally declines.  omega borrowed. */
static Expr* lie_abaco2_unique_unknown(const Expr* omega, const char* xv,
                                       const char* Yn, const char* yname) {
    Expr** ker = NULL; size_t nk = 0, cap = 0;
    lie_collect_kernels(omega, xv, Yn, &ker, &nk, &cap);
    Expr* G = NULL;
    for (size_t i = 0; i < nk && !G; i++) {
        Expr* Mx = ds_d(expr_copy(ker[i]), expr_new_symbol(xv));
        Expr* My = ds_d(expr_copy(ker[i]), expr_new_symbol(Yn));
        if (lie_is_zero(expr_copy(Mx)) || lie_is_zero(expr_copy(My))) {
            expr_free(Mx); expr_free(My); continue;       /* M free of x or of y */
        }
        Expr* R = lie_ratsimp(eval_and_free(ds_call2(SYM_Times, My, powi(Mx, -1))));
        Expr* Rinv = lie_ratsimp(powi(expr_copy(R), -1));  /* 1/R */
        Expr* X = lie_sep_xfactor(R, xv, Yn);              /* x-factor if separable */
        if (X) {
            /* scheme (iib): [X, -X/R] then [-R/X, 1/X] */
            if (!G) G = lie_try_symmetry(expr_copy(X),
                            lie_ratsimp(eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
                                ds_call2(SYM_Times, expr_copy(X), expr_copy(Rinv))))),
                            omega, xv, Yn, yname);
            if (!G) G = lie_try_symmetry(lie_neg_R_over_X(R, X), lie_ratsimp(powi(expr_copy(X), -1)),
                            omega, xv, Yn, yname);
            expr_free(X);
        }
        /* §4.4.1-general order-zero candidates (no separability needed) */
        if (!G) G = lie_try_symmetry(lie_ratsimp(eval_and_free(ds_call2(SYM_Times,
                        expr_new_integer(-1), expr_copy(R)))), expr_new_integer(1),
                        omega, xv, Yn, yname);                            /* [-R, 1] */
        if (!G) G = lie_try_symmetry(expr_new_integer(1),
                        lie_ratsimp(eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
                            expr_copy(R)))), omega, xv, Yn, yname);       /* [1, -R] */
        if (!G) G = lie_try_symmetry(expr_new_integer(1),
                        lie_ratsimp(eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
                            expr_copy(Rinv)))), omega, xv, Yn, yname);    /* [1, -1/R] */
        expr_free(R); expr_free(Rinv);
    }
    for (size_t i = 0; i < nk; i++) expr_free(ker[i]);
    free(ker);
    return G;
}

/* ============================ chi (CPC 101, 5th algorithm) ==================== *
 * The reformulation eta = xi omega + chi.  Because S(xi, xi omega + chi) reduces to
 * S(0, chi), a chi solving the linear PDE
 *     chi_x + omega chi_y - omega_y chi == 0                                (Eq 10)
 * is, for ANY xi, a symmetry [0, chi].  chi is sought by a RICH-BASIS ansatz: a
 * polynomial (undetermined constant coefficients, optionally times a low-degree x,y
 * factor) in the transcendental ATOMS of omega — the Sin/Cos and Sinh/Cosh base pairs,
 * Log, the inverse-trig functions — divided by the product Dtrans of those atoms (which
 * supplies the reciprocal / denominator structure a plain polynomial cannot reach).
 * Substituting into Eq 10, clearing denominators (multiplying by a high power of every
 * atom, which also folds Csc/Sec/Cot/Tan back to Sin/Cos), replacing each atom by an
 * independent generator, reducing modulo the Pythagorean relations of each pair, and
 * splitting the coefficients of the resulting polynomial identity gives a homogeneous
 * linear determining system; its NullSpace is a basis of chi's.  Each nonzero chi is
 * gated by lie_check and integrated by lie_first_integral.  This is the one heuristic
 * whose chi may be a genuine transcendental beyond the polynomial reach of `bivariate`
 * — e.g. Kamke ODE 357  x y' ln(x) sin(y) + cos(y)(1 - x cos(y)) == 0  ->
 * chi = Cos[y]^2/(Log[x] Sin[y]).  First cut: elementary-transcendental atoms only (an
 * undefined function -> decline, the abaco2_unique_* domain; Exp / special functions
 * are not folded, so the PolynomialQ gate declines them).  Source: Cheb-Terrab, Duarte
 * & da Mota, CPC 101 (1997), Eqs 9-10 and the 5th algorithm. */

/* Atom-head classification: 1 circular trig (base pair Sin/Cos), 2 hyperbolic
 * (Sinh/Cosh), 3 other elementary transcendental kept as-is (Log + inverse trig,
 * whose derivatives are rational or self and spawn no new atom), 0 not an atom head. */
static int lie_chi_kind(const char* h) {
    if (h == SYM_Sin || h == SYM_Cos || h == SYM_Tan ||
        h == SYM_Csc || h == SYM_Sec || h == SYM_Cot) return 1;
    if (h == SYM_Sinh || h == SYM_Cosh || h == SYM_Tanh ||
        h == SYM_Csch || h == SYM_Sech || h == SYM_Coth) return 2;
    if (h == SYM_Log || h == intern_symbol("ArcSin") || h == intern_symbol("ArcCos") ||
        h == intern_symbol("ArcTan") || h == intern_symbol("ArcCot") ||
        h == intern_symbol("ArcSinh") || h == intern_symbol("ArcCosh") ||
        h == intern_symbol("ArcTanh")) return 3;
    return 0;
}

/* True iff e contains a circular or hyperbolic trig application (Sin/Cos/.../Sinh/...).
 * On such an omega the rational/algebraic classifier ansatze (abaco1_product,
 * function_sum, abaco2_similar) cannot succeed and their Simplify/Cancel/Integrate
 * construction is very slow (Csc/Cot normalization blows up), so the cascade skips them
 * and lets `chi` — which is built for exactly the trig case — handle it. */
static bool lie_has_trig(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) { int k = lie_chi_kind(h->data.symbol.name); if (k == 1 || k == 2) return true; }
    if (lie_has_trig(h)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (lie_has_trig(e->data.function.args[i])) return true;
    return false;
}

/* Append a copy of item to *arr (grown by realloc) if no expr-equal element is there. */
static void lie_add_unique(Expr*** arr, size_t* n, size_t* cap, const Expr* item) {
    for (size_t i = 0; i < *n; i++) if (expr_eq((*arr)[i], (Expr*)item)) return;
    if (*n == *cap) { *cap = *cap ? *cap * 2 : 4; *arr = realloc(*arr, *cap * sizeof(Expr*)); }
    (*arr)[(*n)++] = expr_copy((Expr*)item);
}

/* Walk e, gathering the distinct circular-trig arguments (circ), hyperbolic arguments
 * (hyp) and other-transcendental applications (other) that contain xv or Yn.  A found
 * atom is maximal: recursion does not descend into it. */
static void lie_chi_collect(const Expr* e, const char* xv, const char* Yn,
                            Expr*** circ, size_t* nc, size_t* cc,
                            Expr*** hyp,  size_t* nh, size_t* ch,
                            Expr*** other, size_t* no, size_t* co) {
    if (!e || e->type != EXPR_FUNCTION) return;
    const Expr* head = e->data.function.head;
    if (head->type == EXPR_SYMBOL && e->data.function.arg_count == 1) {
        int k = lie_chi_kind(head->data.symbol.name);
        const Expr* arg = e->data.function.args[0];
        bool has = ds_contains(arg, xv) || ds_contains(arg, Yn);
        if (k == 1 && has) { lie_add_unique(circ, nc, cc, arg); return; }
        if (k == 2 && has) { lie_add_unique(hyp,  nh, ch, arg); return; }
        if (k == 3 && has) { lie_add_unique(other, no, co, e); return; }
    }
    lie_chi_collect(head, xv, Yn, circ, nc, cc, hyp, nh, ch, other, no, co);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        lie_chi_collect(e->data.function.args[i], xv, Yn, circ, nc, cc, hyp, nh, ch, other, no, co);
}

/* Fresh generator symbol "DSolve`chiG<k>". */
static const char* lie_chi_gen(int k) {
    char buf[32]; snprintf(buf, sizeof buf, "DSolve`chiG%d", k); return intern_symbol(buf);
}

/* The Pythagorean reduction rule  gc^n_ /; n>=2 :> (1 -/+ gs^2) gc^(n-2)  (minus for a
 * circular pair, plus for hyperbolic), built unevaluated for ReplaceRepeated. */
static Expr* lie_pythag_rule(const char* gc, const char* gs, bool hyperbolic) {
    const char* nn = intern_symbol("DSolve`chiN");
    Expr* patt = expr_new_function(expr_new_symbol(SYM_Pattern), (Expr*[]){
        expr_new_symbol(nn), expr_new_function(expr_new_symbol(SYM_Blank), NULL, 0) }, 2);
    Expr* lhsP = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_new_symbol(gc), patt }, 2);
    Expr* cond = expr_new_function(expr_new_symbol(SYM_GreaterEqual),
        (Expr*[]){ expr_new_symbol(nn), expr_new_integer(2) }, 2);
    Expr* lhs = expr_new_function(expr_new_symbol(SYM_Condition), (Expr*[]){ lhsP, cond }, 2);
    Expr* gs2 = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_new_symbol(gs), expr_new_integer(2) }, 2);
    Expr* base = hyperbolic ? ds_call2(SYM_Plus, expr_new_integer(1), gs2)
                            : ds_call2(SYM_Subtract, expr_new_integer(1), gs2);
    Expr* nm2 = ds_call2(SYM_Plus, expr_new_symbol(nn), expr_new_integer(-2));
    Expr* rhs = ds_call2(SYM_Times, base,
        expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ expr_new_symbol(gc), nm2 }, 2));
    return expr_new_function(expr_new_symbol(SYM_RuleDelayed), (Expr*[]){ lhs, rhs }, 2);
}

/* base^k as an evaluated Power; base borrowed. */
static Expr* lie_powk(const Expr* base, int k) {
    return powi(expr_copy((Expr*)base), k);
}

/* Enumerate the atom-monomials  prod atoms[i]^e_i  with  sum e_i <= D0  (evaluated,
 * deduplicated).  Returns a malloc'd array (caller frees each + the array). */
static Expr** lie_atom_monomials(Expr** atoms, size_t m, int D0, size_t* out_n) {
    Expr** list = NULL; size_t n = 0, cap = 0;
    Expr* one = expr_new_integer(1);
    lie_add_unique(&list, &n, &cap, one); expr_free(one);              /* seed {1} */
    /* level-by-level: multiply the previous level's monomials by each atom */
    size_t level_start = 0;
    for (int d = 1; d <= D0; d++) {
        size_t level_end = n;
        for (size_t i = level_start; i < level_end; i++)
            for (size_t a = 0; a < m; a++) {
                Expr* prod = eval_and_free(ds_call2(SYM_Times, expr_copy(list[i]), expr_copy(atoms[a])));
                lie_add_unique(&list, &n, &cap, prod);
                expr_free(prod);
            }
        level_start = level_end;
    }
    *out_n = n;
    return list;
}

/* {x, Y, gens...} as a List for CoefficientList / PolynomialQ. */
static Expr* lie_chi_varlist(const char* xv, const char* Yn, const char** gen, size_t ng) {
    Expr** a = malloc((ng + 2) * sizeof(Expr*));
    a[0] = expr_new_symbol(xv); a[1] = expr_new_symbol(Yn);
    for (size_t i = 0; i < ng; i++) a[i + 2] = expr_new_symbol(gen[i]);
    Expr* L = expr_new_function(expr_new_symbol(SYM_List), a, ng + 2);
    free(a);
    return L;
}

/* Heuristic `chi`: the eta = xi omega + chi reformulation (CPC 101, 5th algorithm). */
static Expr* lie_chi(const Expr* omega, const char* xv, const char* Yn, const char* yname) {
    if (lie_has_undefined_function(omega)) return NULL;   /* abaco2_unique_* domain */

    /* 1. collect atoms */
    Expr **circ = NULL, **hyp = NULL, **other = NULL;
    size_t nci = 0, cci = 0, nhy = 0, chy = 0, noi = 0, coi = 0;
    lie_chi_collect(omega, xv, Yn, &circ, &nci, &cci, &hyp, &nhy, &chy, &other, &noi, &coi);

    /* base atoms = {Sin,Cos}(circ) + {Sinh,Cosh}(hyp) + other; gens 1:1; pair table.
     * m == 0 -> rational omega (bivariate covers it); m > 5 -> too many atoms (cost). */
    size_t m = 2 * nci + 2 * nhy + noi;
    Expr* G = NULL;
    if (m >= 1 && m <= 5) {
    Expr** atoms = malloc(m * sizeof(Expr*));
    const char** gen = malloc(m * sizeof(char*));
    /* pair[i] = index of the partner Cos/Cosh for a Sin/Sinh (for the reduction); -1 else */
    int* pair_gc = malloc(m * sizeof(int));   /* index of Cos gen */
    int* pair_gs = malloc(m * sizeof(int));   /* index of Sin gen */
    int* pair_hyp = malloc(m * sizeof(int));
    size_t np = 0, mi = 0;
    for (size_t i = 0; i < nci; i++) {
        atoms[mi] = ds_call1(SYM_Sin, expr_copy(circ[i])); atoms[mi] = eval_and_free(atoms[mi]);
        int gs = (int)mi; gen[mi] = lie_chi_gen((int)mi + 1); mi++;
        atoms[mi] = eval_and_free(ds_call1(SYM_Cos, expr_copy(circ[i])));
        int gc = (int)mi; gen[mi] = lie_chi_gen((int)mi + 1); mi++;
        pair_gc[np] = gc; pair_gs[np] = gs; pair_hyp[np] = 0; np++;
    }
    for (size_t i = 0; i < nhy; i++) {
        atoms[mi] = eval_and_free(ds_call1(SYM_Sinh, expr_copy(hyp[i])));
        int gs = (int)mi; gen[mi] = lie_chi_gen((int)mi + 1); mi++;
        atoms[mi] = eval_and_free(ds_call1(SYM_Cosh, expr_copy(hyp[i])));
        int gc = (int)mi; gen[mi] = lie_chi_gen((int)mi + 1); mi++;
        pair_gc[np] = gc; pair_gs[np] = gs; pair_hyp[np] = 1; np++;
    }
    for (size_t i = 0; i < noi; i++) { atoms[mi] = expr_copy(other[i]); gen[mi] = lie_chi_gen((int)mi + 1); mi++; }

    /* Dtrans = product of all atoms */
    Expr* Dtrans = expr_new_integer(1);
    for (size_t i = 0; i < m; i++) Dtrans = ds_call2(SYM_Times, Dtrans, expr_copy(atoms[i]));
    Dtrans = eval_and_free(Dtrans);

    for (int D0 = 2; D0 <= 3 && !G; D0++) {
        /* 2. object monomials = atom-monomials(<=D0) x {1, x, y} */
        size_t nam = 0;
        Expr** am = lie_atom_monomials(atoms, m, D0, &nam);
        Expr* xyf[3] = { expr_new_integer(1), expr_new_symbol(xv), expr_new_symbol(Yn) };
        size_t nmon = nam * 3;
        if (nmon > 90) {                                             /* cost cap */
            for (size_t i = 0; i < nam; i++) expr_free(am[i]);
            free(am);
            for (int t = 0; t < 3; t++) expr_free(xyf[t]);
            continue;
        }
        Expr** mons = malloc(nmon * sizeof(Expr*));
        for (size_t i = 0; i < nam; i++)
            for (int t = 0; t < 3; t++)
                mons[i * 3 + t] = eval_and_free(ds_call2(SYM_Times, expr_copy(am[i]), expr_copy(xyf[t])));
        for (size_t i = 0; i < nam; i++) expr_free(am[i]);
        free(am);
        for (int t = 0; t < 3; t++) expr_free(xyf[t]);

        /* 3. chi = (sum c_k mon_k)/Dtrans, with fresh coefficient symbols c_k */
        Expr** csym = malloc(nmon * sizeof(Expr*));
        Expr* Nsum = expr_new_integer(0);
        for (size_t k = 0; k < nmon; k++) {
            char cb[40]; snprintf(cb, sizeof cb, "DSolve`chiC%zu", k + 1);
            csym[k] = expr_new_symbol(intern_symbol(cb));
            Nsum = ds_call2(SYM_Plus, Nsum, ds_call2(SYM_Times, expr_copy(csym[k]), expr_copy(mons[k])));
        }
        Nsum = eval_and_free(Nsum);
        Expr* chi = eval_and_free(ds_call2(SYM_Times, Nsum, powi(expr_copy(Dtrans), -1)));

        /* 4. E = chi_x + omega chi_y - omega_y chi */
        Expr* chix = ds_d(expr_copy(chi), expr_new_symbol(xv));
        Expr* chiy = ds_d(expr_copy(chi), expr_new_symbol(Yn));
        Expr* omy  = ds_d(expr_copy((Expr*)omega), expr_new_symbol(Yn));
        Expr* Eexpr = eval_and_free(ds_call2(SYM_Plus, chix,
                        ds_call2(SYM_Subtract,
                            ds_call2(SYM_Times, expr_copy((Expr*)omega), chiy),
                            ds_call2(SYM_Times, omy, expr_copy(chi)))));
        expr_free(chi);

        /* 5. clear denominators + fold reciprocal trig to Sin/Cos: x^(D0+3) prod atom^(D0+3) */
        Expr* clearing = powi(expr_new_symbol(xv), D0 + 3);
        for (size_t i = 0; i < m; i++)
            clearing = ds_call2(SYM_Times, clearing, lie_powk(atoms[i], D0 + 3));
        Expr* Ecl = eval_and_free(ds_call1("Expand",
                        eval_and_free(ds_call2(SYM_Times, Eexpr, eval_and_free(clearing)))));

        /* 6. substitute each atom -> its generator */
        Expr* Esub = Ecl;
        for (size_t i = 0; i < m; i++)
            Esub = ds_subst(Esub, expr_copy(atoms[i]), expr_new_symbol(gen[i]));
        /* 7. En = Numerator[Together[Esub]] and require it polynomial in x,Y,gens */
        Expr* En = eval_and_free(ds_call1(SYM_Numerator, eval_and_free(ds_call1(SYM_Together, Esub))));
        Expr* vars = lie_chi_varlist(xv, Yn, gen, m);
        /* En must be polynomial in x, y and the generators for CoefficientList to split
         * it; a leftover (Exp / special-function / un-folded reciprocal) makes it not —
         * decline this D0.  The ansatz caps (m<=5, nmon<=90, D0<=3) already bound cost,
         * so no node-budget check here (chi's determining polynomial is legitimately
         * larger than the fast helpers' inputs). */
        bool polyq = lie_ev_trueQ(eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
                        (Expr*[]){ expr_copy(En), expr_copy(vars) }, 2)));
        if (!polyq) {
            expr_free(En); expr_free(vars);
            for (size_t k=0;k<nmon;k++){expr_free(csym[k]);expr_free(mons[k]);}
            free(csym); free(mons);
            continue;
        }
        /* 8. reduce modulo the Pythagorean relations of each pair */
        for (size_t p = 0; p < np; p++) {
            Expr* rule = lie_pythag_rule(gen[pair_gc[p]], gen[pair_gs[p]], pair_hyp[p]);
            En = eval_and_free(ds_call2(SYM_ReplaceRepeated, En, rule));
        }
        /* 9. determining system: forms = Flatten[CoefficientList[En, vars]] */
        Expr* cl = eval_and_free(expr_new_function(expr_new_symbol(SYM_CoefficientList),
                        (Expr*[]){ En, vars }, 2));                    /* consumes En, vars */
        Expr* forms = eval_and_free(ds_call1(SYM_Flatten, cl));
        Expr** cargs = malloc(nmon * sizeof(Expr*));
        for (size_t k = 0; k < nmon; k++) cargs[k] = expr_copy(csym[k]);
        Expr* cvec = expr_new_function(expr_new_symbol(SYM_List), cargs, nmon);
        free(cargs);
        Expr* M = eval_and_free(expr_new_function(expr_new_symbol("Outer"),
                        (Expr*[]){ expr_new_symbol("Coefficient"), forms, cvec }, 3));
        Expr* ns = eval_and_free(ds_call1("NullSpace", M));

        /* 10. Reconstruct the nonzero chi candidates and integrate them SIMPLEST-FIRST
         * (a smaller chi has a simpler — faster, more likely elementary — Lie
         * quadrature; a messy sibling can send Integrate into a long trig search).  The
         * null vectors are exact Eq-10 solutions by construction, so lie_check is
         * skipped here (its polynomial fast-path would even mis-reject an identity-based
         * trig zero); the returned first integral is still verified by
         * dsolve_run_implicit's back-substitution. */
        if (ns && ns->type == EXPR_FUNCTION) {
            size_t nv = ns->data.function.arg_count;
            Expr** cand = malloc((nv ? nv : 1) * sizeof(Expr*));
            long*  csz  = malloc((nv ? nv : 1) * sizeof(long));
            size_t ncand = 0;
            for (size_t r = 0; r < nv; r++) {
                Expr* v = ns->data.function.args[r];
                if (!v || v->type != EXPR_FUNCTION || v->data.function.arg_count != nmon) continue;
                Expr* Nk = expr_new_integer(0);
                for (size_t k = 0; k < nmon; k++)
                    Nk = ds_call2(SYM_Plus, Nk, ds_call2(SYM_Times,
                            expr_copy(v->data.function.args[k]), expr_copy(mons[k])));
                Expr* chik = lie_ratsimp(eval_and_free(ds_call2(SYM_Times, eval_and_free(Nk),
                                 powi(expr_copy(Dtrans), -1))));
                if (lie_is_zero(expr_copy(chik))) { expr_free(chik); continue; }
                long acc = 0; lie_count_nodes(chik, &acc, 1L << 30);
                cand[ncand] = chik; csz[ncand] = acc; ncand++;
            }
            for (size_t a = 0; a < ncand; a++)          /* selection sort by size asc. */
                for (size_t b = a + 1; b < ncand; b++)
                    if (csz[b] < csz[a]) { long ts = csz[a]; csz[a] = csz[b]; csz[b] = ts;
                                           Expr* te = cand[a]; cand[a] = cand[b]; cand[b] = te; }
            Expr* zero = expr_new_integer(0);
            for (size_t a = 0; a < ncand && !G; a++)
                G = lie_first_integral(zero, cand[a], omega, xv, Yn, yname);
            expr_free(zero);
            for (size_t a = 0; a < ncand; a++) expr_free(cand[a]);
            free(cand); free(csz);
        }
        expr_free(ns);
        for (size_t k=0;k<nmon;k++){expr_free(csym[k]);expr_free(mons[k]);}
        free(csym); free(mons);
    }

    expr_free(Dtrans);
    for (size_t i = 0; i < m; i++) expr_free(atoms[i]);
    free(atoms); free(gen); free(pair_gc); free(pair_gs); free(pair_hyp);
    }  /* end if (m in [1,5]) */

    for (size_t i = 0; i < nci; i++) expr_free(circ[i]);
    for (size_t i = 0; i < nhy; i++) expr_free(hyp[i]);
    for (size_t i = 0; i < noi; i++) expr_free(other[i]);
    free(circ); free(hyp); free(other);
    return G;
}

Expr** dsolve_lie_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xv    = P->ind_names[0];

    /* omega(x, y) = the RHS of y' = omega, with y[x] -> the plain symbol Y */
    Expr* F = dsolve_solve_top_derivative(P, 1);
    if (!F) return NULL;
    const char* Yn = intern_symbol("DSolve`Y");
    Expr* omega = ds_subst(F, ds_make_funcapp(yname, 0, xv), expr_new_symbol(Yn));

    Expr* wx = ds_d(expr_copy(omega), expr_new_symbol(xv));
    Expr* wy = ds_d(expr_copy(omega), expr_new_symbol(Yn));

    /* Heuristics are tried CHEAPEST-FIRST (each yields an equally-verified first
     * integral, so the order is a pure cost choice; the first match wins):
     *   abaco1_simple   one-variable ansatze (a few derivatives + a free-of test)
     *   linear          affine ansatz, a 6-unknown determining NullSpace
     *   abaco1_product  rational-product ansatz + its inverse (two separability
     *                   free-of tests) — the quadrature ansatze chain in here
     *   function_sum    additive ansatz [F(x)+G(y), 0] + inverse (§4.2)
     *   abaco2_similar  similarity ansatz [F(x), H(x)] + inverse (§4.3) — the first
     *                   to reach irrational omega (Sqrt/root forms)
     *   abaco2_unique_unknown  [F(x),G(y)] / [G(y),F(x)] from the functions /
     *                   non-integer powers of both variables in omega (§4.4.1)
     *   bivariate       degree-2 THEN degree-3 polynomial NullSpace (up to a
     *                   20-unknown system).
     *   chi             the eta = xi omega + chi reformulation (CPC 101, 5th algorithm):
     *                   a rich-basis (transcendental atoms of omega) determining system
     *                   for chi — the richest and last-resort heuristic. */
    /* An omega carrying an undefined function (F[x^2+y^2], ...) is the domain of the
     * abaco2_unique_* heuristics, which read it off the cheap R = M_y/M_x route.  The
     * rational/algebraic classifier ansatze (abaco1_product, function_sum,
     * abaco2_similar) cannot solve such an ODE, and building their classifiers grinds
     * the general Simplify/Cancel/Integrate for tens of seconds on the transcendental
     * derivatives of an arbitrary function — so skip them here (measured: this cuts a
     * decline on y' == Tan[ArcTan[y] + F[x^2+y^2]] from ~80s to well under a second). */
    bool undef = lie_has_undefined_function(omega);
    bool skipcls = undef || lie_has_trig(omega);   /* skip the rational/algebraic classifiers */
    Expr* G = lie_abaco1_simple(omega, wx, wy, xv, Yn, yname);
    if (!G) G = lie_linear(omega, xv, Yn, yname);
    if (!G && !skipcls) G = lie_abaco1_product(omega, xv, Yn, yname);
    if (!G && !skipcls) G = lie_function_sum(omega, xv, Yn, yname);
    if (!G && !skipcls) G = lie_abaco2_similar(omega, xv, Yn, yname);
    if (!G) G = lie_abaco2_unique_unknown(omega, xv, Yn, yname);
    if (!G) G = lie_bivariate(omega, xv, Yn, yname);
    if (!G && !undef) G = lie_chi(omega, xv, Yn, yname);

    expr_free(wx); expr_free(wy); expr_free(omega);
    if (!G) return NULL;

    Expr** out = malloc(sizeof(Expr*));
    out[0] = G;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_lie(Expr* res) {
    return dsolve_method_builtin_implicit(res, dsolve_lie_try);
}

void dsolve_lie_init(void) {
    static const char* doc =
        "DSolve`LieSymmetry[eqn, y, x] integrates the first-order ODE y' == omega(x,y) "
        "by finding a Lie point symmetry (xi d/dx + eta d/dy) via a table of ansatz "
        "heuristics, then reducing to a quadrature through the integrating factor "
        "mu == 1/(eta - xi omega); returns the implicit first integral "
        "{{G(x,y[x]) == C[1]}}. The general first-order backstop of the cascade "
        "(heuristics: abaco1_simple, linear, abaco1_product, function_sum, "
        "abaco2_similar, abaco2_unique_unknown (incl. the differential-invariant-of-"
        "order-zero candidates), bivariate, and chi (the eta == xi omega + chi "
        "reformulation with a rich transcendental-atom basis, reaching genuinely "
        "transcendental chi such as Cos[y]^2/(Log[x] Sin[y]))).";
    symtab_add_builtin("DSolve`LieSymmetry", builtin_dsolve_lie);
    symtab_get_def("DSolve`LieSymmetry")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`LieSymmetry", doc);
    symtab_add_builtin("DSolve`LieGroup", builtin_dsolve_lie);   /* SymPy-style alias */
    symtab_get_def("DSolve`LieGroup")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`LieGroup",
        "DSolve`LieGroup[eqn, y, x] is an alias for DSolve`LieSymmetry (the name of "
        "SymPy's lie_group hint) — the heuristic Lie point-symmetry solver for "
        "first-order ODEs. See DSolve`LieSymmetry.");
}
