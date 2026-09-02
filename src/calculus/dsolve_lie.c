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
 * ansatz — the linear-coefficients class, via the determining-system NullSpace);
 * L3 adds `bivariate` (general degree-2/3 polynomial ansatz — the same NullSpace
 * determining system at higher degree, catching quadratic / projective symmetries);
 * the remaining heuristics chain in through the same lie_check + lie_first_integral
 * pipeline.
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

/* base^n as an evaluated Power; base consumed. */
static Expr* powi(Expr* base, int n) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                         (Expr*[]){ base, expr_new_integer(n) }, 2));
}

/* Cancel[Together[e]] — the cheap rational-function normalizer (over a common
 * denominator, then cancel the gcd).  Used in place of the far costlier general
 * Simplify on the heuristics' hot path: it is all these rational (xi, eta) forms
 * need, and it keeps the per-decline cost down to a few ms.  e consumed. */
static Expr* lie_ratsimp(Expr* e) {
    return eval_and_free(ds_call1("Cancel", eval_and_free(ds_call1(SYM_Together, e))));
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
 * the residual before the zero test so a rational S that cancels is recognized. */
static bool lie_check(const Expr* xi, const Expr* eta, const Expr* omega,
                      const char* xv, const char* Yn) {
    Expr* S  = lie_S_expr(xi, eta, omega, xv, Yn);
    Expr* St = eval_and_free(ds_call1("Together", S));   /* consumes S */
    bool z = ds_is_zero(St);
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

/* Exp[Integrate[R, var]], or NULL if that integral is non-elementary.  R consumed. */
static Expr* lie_exp_integral(Expr* R, const char* var) {
    g_integrate_quiet++;
    Expr* I = lie_integrate_signfix(R, var);                     /* consumes R */
    g_integrate_quiet--;
    if (ds_has_head(I, SYM_Integrate)) { expr_free(I); return NULL; }
    return eval_and_free(ds_call1("Exp", I));                    /* consumes I */
}

/* Fast "e is free of var" test, e consumed.  For a rational e = P/Q, e is free of
 * var iff its derivative numerator  P_var Q - P Q_var  is identically zero, which is
 * a POLYNOMIAL zero-test after Expand — far cheaper than the general zero_test on
 * the rational derivative D[e, var] (the dominant per-decline cost of the quadrature
 * heuristics: ~12ms vs ~2ms on a moderate rational).  ds_is_zero on the expanded
 * numerator still handles a transcendental e correctly (it just falls back to the
 * full zero_test there), so this only speeds up — never changes — the answer. */
static bool lie_free_of_var(Expr* e, const char* var) {
    Expr* t = eval_and_free(ds_call1(SYM_Together, e));         /* one coprime fraction */
    Expr* P = eval_and_free(ds_call1(SYM_Numerator, expr_copy(t)));
    Expr* Q = eval_and_free(ds_call1(SYM_Denominator, t));      /* consumes t */
    Expr* Pv = ds_d(expr_copy(P), expr_new_symbol(var));
    Expr* Qv = ds_d(expr_copy(Q), expr_new_symbol(var));
    Expr* num = eval_and_free(ds_call2(SYM_Subtract,            /* P_var Q - P Q_var */
                    ds_call2(SYM_Times, Pv, Q),
                    ds_call2(SYM_Times, P, Qv)));               /* consumes Pv,Q,P,Qv */
    Expr* ex = eval_and_free(ds_call1("Expand", num));          /* consumes num */
    bool z = ds_is_zero(ex);
    expr_free(ex);
    return z;
}

/* The product-separable x-factor of L: Exp[Integrate[L_x/L, x]].  Returns NULL when
 * L == 0, when L_x/L is not free of y (so L does not separate as X(x) Y(y): note
 * L_x/L free of y  <=>  d^2/dx dy log L = 0  <=>  L separable), or when the integral
 * is non-elementary.  Used by the quadrature heuristics to extract F(x) / 1/F_xx
 * from the necessary-condition expression.  L borrowed; result owned. */
static Expr* lie_sep_xfactor(const Expr* L, const char* xv, const char* Yn) {
    if (ds_is_zero(L)) return NULL;
    Expr* Lx = ds_d(expr_copy((Expr*)L), expr_new_symbol(xv));
    /* ratio = Together[L_x / L].  Together (value-preserving, so it cannot change
     * the free-of answer) collapses the L^4-denominator mess into one fraction. */
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
        Expr* ratio = ds_simplify(eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)wx),
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
        Expr* ratio = ds_simplify(eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)wy),
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
     *   bivariate       degree-2 THEN degree-3 polynomial NullSpace (up to a
     *                   20-unknown system) — the most expensive, tried last. */
    Expr* G = lie_abaco1_simple(omega, wx, wy, xv, Yn, yname);
    if (!G) G = lie_linear(omega, xv, Yn, yname);
    if (!G) G = lie_abaco1_product(omega, xv, Yn, yname);
    if (!G) G = lie_bivariate(omega, xv, Yn, yname);

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
        "(heuristics: abaco1_simple, linear, abaco1_product, bivariate).";
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
