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
    Expr* Phi = ds_integrate(Pc, expr_new_symbol(xv));           /* consumes Pc */
    g_integrate_quiet--;
    if (ds_has_head(Phi, SYM_Integrate)) { expr_free(Phi); expr_free(mu); return NULL; }

    /* corr = Q - d(Phi)/dY = mu - d(Phi)/dY  (exactness => free of x) */
    Expr* dPhiY = ds_d(expr_copy(Phi), expr_new_symbol(Yn));
    Expr* corr = ds_simplify(eval_and_free(ds_call2(SYM_Subtract, expr_copy(mu), dPhiY)));
    expr_free(mu);
    if (!ds_free_of(corr, xv)) { expr_free(corr); expr_free(Phi); return NULL; }

    g_integrate_quiet++;
    Expr* Fy = ds_integrate(corr, expr_new_symbol(Yn));          /* consumes corr */
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
    Expr* I = ds_integrate(R, expr_new_symbol(var));             /* consumes R */
    g_integrate_quiet--;
    if (ds_has_head(I, SYM_Integrate)) { expr_free(I); return NULL; }
    return eval_and_free(ds_call1("Exp", I));                    /* consumes I */
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
                if (lie_check(g, zero, omega, xv, Yn))
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

    /* L1: abaco1_simple; L2: linear (affine symmetry).  (L2/L3 rest chains here.) */
    Expr* G = lie_abaco1_simple(omega, wx, wy, xv, Yn, yname);
    if (!G) G = lie_linear(omega, xv, Yn, yname);

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
        "(heuristics: abaco1_simple, linear).";
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
