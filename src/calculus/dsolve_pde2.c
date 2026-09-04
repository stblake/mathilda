/*
 * dsolve_pde2.c — second-order constant-coefficient linear PDE via operator
 * factoring (the plan's PDEHyperbolicGeneral, realized in full generality).
 *
 * Solves the homogeneous, principal-part-only, constant-coefficient equation
 *
 *     A u_{v1 v1} + B u_{v1 v2} + C u_{v2 v2} == 0     (A, B, C constant)
 *
 * The trial u = f(v2 + lambda v1) reduces it to the characteristic quadratic
 * A lambda^2 + B lambda + C == 0, so the principal operator factors over C and
 * one method covers all three discriminant signs:
 *
 *   - distinct real roots  (hyperbolic) -> u = F(v2 + l1 v1) + G(v2 + l2 v1)
 *       e.g. the wave equation  u_tt == c^2 u_xx  ->  F(x - c t) + G(x + c t)
 *   - complex conjugate roots (elliptic) -> the complex-characteristic form,
 *       e.g. Laplace  u_xx + u_yy == 0  ->  F(y - I x) + G(y + I x)
 *   - a repeated root        (parabolic) -> u = F(w) + v1 G(w), w = v2 + l v1
 *
 * Complex/repeated roots are handled uniformly by dsolve_analyze_roots (which
 * returns distinct roots with multiplicities), so no realification is needed —
 * Mathematica returns the same complex-characteristic form for elliptic PDEs.
 *
 * LOWER-ORDER TERMS.  The method also solves the constant-coefficient equation
 * with lower-order terms
 *
 *     A u_{v1 v1} + B u_{v1 v2} + C u_{v2 v2} + D u_{v1} + E u_{v2} + F u == 0
 *
 * WHEN the full symbol A ξ² + B ξη + C η² + D ξ + E η + F factors into two
 * first-order constant-coefficient operators (∂_v1 − l_i ∂_v2 + m_i).  Then each
 * factor contributes e^{−m_i v1} C[i][v2 + l_i v1], so the general solution is
 * the exponential-damped arbitrary-function form — e.g. the distortionless
 * telegraph  u_tt − c² u_xx + a u_t + (a²/4) u == 0  ->
 * e^{−a t/2}(C[1][x + c t] + C[2][x − c t]).  When the operator does NOT factor
 * over the constants (the general telegraph b ≠ a²/4 needs Bessel functions),
 * the method declines.  D = E = F = 0 gives m_i = 0, the principal-part form.
 *
 * First cut (honest): homogeneous (zero forcing), constant coefficients, symbol
 * factoring over the constants.  Inhomogeneous forcing, the non-factorable
 * general telegraph, and the wave-IVP d'Alembert formula are future work.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

static Expr* mul(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Times, a, b)); }
static Expr* add(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Plus, a, b)); }

/* Derivative[o1,o2][u][v1,v2] */
static Expr* pde_deriv(const char* u, int o1, int o2, const char* v1, const char* v2) {
    Expr* d = expr_new_function(expr_new_symbol(SYM_Derivative),
                  (Expr*[]){ expr_new_integer(o1), expr_new_integer(o2) }, 2);
    Expr* du = expr_new_function(d, (Expr*[]){ expr_new_symbol(u) }, 1);
    return expr_new_function(du, (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2);
}

/* C[k][arg] — the generated arbitrary function of one variable (as in pde1). */
static Expr* arbfun(int k, Expr* arg) {
    return expr_new_function(ds_const(k), (Expr*[]){ arg }, 1);
}

static Expr* inv(Expr* a) { return expr_new_function(expr_new_symbol(SYM_Power),
                                (Expr*[]){ a, expr_new_integer(-1) }, 2); }

/* e^{-m dampvar} * x  (x consumed).  A zero shift m leaves x undamped. */
static Expr* damp(const Expr* m, Expr* x, const char* dampvar) {
    if (ds_is_zero((Expr*)m)) return x;
    Expr* expo = eval_and_free(ds_call1("Exp",
                     mul(mul(expr_new_integer(-1), expr_copy((Expr*)m)),
                         expr_new_symbol(dampvar))));
    return mul(expo, x);
}

/* e^{-m v1} C[k][v2 + lam v1] — the kernel of the first-order factor
 * (∂_v1 − lam ∂_v2 + m). */
static Expr* damped_term(int k, const Expr* m, const Expr* lam,
                         const char* v1, const char* v2) {
    Expr* xi = add(expr_new_symbol(v2), mul(expr_copy((Expr*)lam), expr_new_symbol(v1)));
    return damp(m, arbfun(k, xi), v1);
}

/* Build the general body from the FULL symbol A l^2 + B l + C + (lower-order
 * shifts) by factoring the constant-coefficient operator into two first-order
 * operators (∂_v1 − l_i ∂_v2 + m_i).  A != 0 assumed.  The zeroth-order shifts
 * m_i solve  m1 + m2 = D/A,  l2 m1 + l1 m2 = −E/A  (distinct l) with the
 * factorability check m1 m2 == F/A; a repeated l needs E + l D == 0, then m1,m2
 * are the roots of m^2 − (D/A) m + F/A.  The general solution is then
 * Σ e^{−m_i v1} C[i][v2 + l_i v1] (a repeated (l,m) carries the extra v1 factor).
 * Returns NULL when the operator does not factor over the constants (needs
 * special functions) or the quadratic does not fully split. */
static Expr* pde2_body(const Expr* A, const Expr* B, const Expr* C,
                       const Expr* D, const Expr* E, const Expr* F,
                       const char* v1, const char* v2) {
    const char* lam = intern_symbol("DSolve`pde2lam");
    Expr* l2 = expr_new_function(expr_new_symbol(SYM_Power),
                   (Expr*[]){ expr_new_symbol(lam), expr_new_integer(2) }, 2);
    Expr* poly = add(add(mul(expr_copy((Expr*)A), l2),
                         mul(expr_copy((Expr*)B), expr_new_symbol(lam))),
                     expr_copy((Expr*)C));
    DSolveRoots R;
    bool ok = dsolve_analyze_roots(poly, lam, 2, &R);
    expr_free(poly);
    if (!ok || R.total != 2) { dsolve_roots_free(&R); return NULL; }

    Expr* Ainv = inv(expr_copy((Expr*)A));
    Expr* body = NULL;
    if (R.ndist == 2) {
        /* distinct l: solve the linear m-system, verify m1 m2 == F/A */
        Expr* l1 = R.roots[0], *l2r = R.roots[1];
        Expr* num = mul(expr_new_integer(-1),
                        add(expr_copy((Expr*)E), mul(expr_copy(l1), expr_copy((Expr*)D))));
        Expr* den = mul(expr_copy((Expr*)A),
                        eval_and_free(ds_call2(SYM_Subtract, expr_copy(l2r), expr_copy(l1))));
        Expr* m1 = ds_simplify(mul(num, inv(den)));
        Expr* m2 = ds_simplify(eval_and_free(ds_call2(SYM_Subtract,
                       mul(expr_copy((Expr*)D), expr_copy(Ainv)), expr_copy(m1))));
        Expr* chk = eval_and_free(ds_call2(SYM_Subtract, mul(expr_copy(m1), expr_copy(m2)),
                        mul(expr_copy((Expr*)F), expr_copy(Ainv))));
        if (ds_is_zero(chk))
            body = add(damped_term(1, m1, l1, v1, v2), damped_term(2, m2, l2r, v1, v2));
        expr_free(m1); expr_free(m2); expr_free(chk);
    } else {
        /* repeated l: consistency E + l D == 0, then m from the shift quadratic */
        Expr* lam0 = R.roots[0];
        Expr* cons = eval_and_free(ds_call2(SYM_Plus, expr_copy((Expr*)E),
                         mul(expr_copy(lam0), expr_copy((Expr*)D))));
        bool consistent = ds_is_zero(cons);
        expr_free(cons);
        if (consistent) {
            const char* mv = intern_symbol("DSolve`pde2m");
            Expr* m2sq = expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ expr_new_symbol(mv), expr_new_integer(2) }, 2);
            Expr* mpoly = add(add(m2sq,
                              mul(expr_new_integer(-1),
                                  mul(mul(expr_copy((Expr*)D), expr_copy(Ainv)), expr_new_symbol(mv)))),
                              mul(expr_copy((Expr*)F), expr_copy(Ainv)));
            DSolveRoots Rm;
            bool okm = dsolve_analyze_roots(mpoly, mv, 2, &Rm);
            expr_free(mpoly);
            if (okm && Rm.total == 2) {
                if (Rm.ndist == 2) {
                    body = add(damped_term(1, Rm.roots[0], lam0, v1, v2),
                               damped_term(2, Rm.roots[1], lam0, v1, v2));
                } else {
                    /* repeated (l,m): e^{−m v1}(C[1][w] + v1 C[2][w]) */
                    Expr* w1 = add(expr_new_symbol(v2), mul(expr_copy(lam0), expr_new_symbol(v1)));
                    Expr* w2 = add(expr_new_symbol(v2), mul(expr_copy(lam0), expr_new_symbol(v1)));
                    Expr* inner = add(arbfun(1, w1), mul(expr_new_symbol(v1), arbfun(2, w2)));
                    body = damp(Rm.roots[0], inner, v1);
                }
            }
            dsolve_roots_free(&Rm);
        }
    }
    expr_free(Ainv);
    dsolve_roots_free(&R);
    return body;
}

/* Pure-mixed operator  B u_{v1 v2} + D u_{v1} + E u_{v2} + F u == 0  (A = C = 0).
 * Factors as (∂_v1 + E/B)(∂_v2 + D/B) B  iff  F == D E / B, giving
 *   u = e^{−(E/B) v1} C[1][v2] + e^{−(D/B) v2} C[2][v1].
 * Returns NULL when it does not factor. */
static Expr* pde2_body_mixed(const Expr* B, const Expr* D, const Expr* E, const Expr* F,
                             const char* v1, const char* v2) {
    Expr* Binv = inv(expr_copy((Expr*)B));
    Expr* eB = ds_simplify(mul(expr_copy((Expr*)E), expr_copy(Binv)));   /* E/B */
    Expr* dB = ds_simplify(mul(expr_copy((Expr*)D), expr_copy(Binv)));   /* D/B */
    Expr* chk = eval_and_free(ds_call2(SYM_Subtract, expr_copy((Expr*)F),
                    mul(mul(expr_copy((Expr*)D), expr_copy((Expr*)E)), expr_copy(Binv))));
    expr_free(Binv);
    Expr* body = NULL;
    if (ds_is_zero(chk)) {
        Expr* t1 = damp(eB, arbfun(1, expr_new_symbol(v2)), v1);
        Expr* t2 = damp(dB, arbfun(2, expr_new_symbol(v1)), v2);
        body = add(t1, t2);
    }
    expr_free(chk); expr_free(eB); expr_free(dB);
    return body;
}

Expr** dsolve_pde2_solve(DSolveProblem* P) {
    if (P->nfun != 1 || P->nind != 2 || P->neq != 1) return NULL;
    const char* uname = P->fun_names[0];
    const char* v1 = P->ind_names[0];
    const char* v2 = P->ind_names[1];

    /* algebraic residual: replace each derivative term with a fresh symbol */
    const char* sxx = intern_symbol("DSolve`pde2Uxx");
    const char* sxy = intern_symbol("DSolve`pde2Uxy");
    const char* syy = intern_symbol("DSolve`pde2Uyy");
    const char* sx  = intern_symbol("DSolve`pde2Ux");
    const char* sy  = intern_symbol("DSolve`pde2Uy");
    const char* su  = intern_symbol("DSolve`pde2U");
    const char* svars[6] = { sxx, sxy, syy, sx, sy, su };

    Expr* R = expr_copy(P->eq_residuals[0]);
    R = ds_subst(R, pde_deriv(uname, 2, 0, v1, v2), expr_new_symbol(sxx));
    R = ds_subst(R, pde_deriv(uname, 1, 1, v1, v2), expr_new_symbol(sxy));
    R = ds_subst(R, pde_deriv(uname, 0, 2, v1, v2), expr_new_symbol(syy));
    R = ds_subst(R, pde_deriv(uname, 1, 0, v1, v2), expr_new_symbol(sx));
    R = ds_subst(R, pde_deriv(uname, 0, 1, v1, v2), expr_new_symbol(sy));
    R = ds_subst(R, expr_new_function(expr_new_symbol(uname),
                     (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2),
                 expr_new_symbol(su));

    Expr* A  = ds_d(expr_copy(R), expr_new_symbol(sxx));
    Expr* B  = ds_d(expr_copy(R), expr_new_symbol(sxy));
    Expr* C  = ds_d(expr_copy(R), expr_new_symbol(syy));
    Expr* Dc = ds_d(expr_copy(R), expr_new_symbol(sx));
    Expr* Ec = ds_d(expr_copy(R), expr_new_symbol(sy));
    Expr* Fc = ds_d(expr_copy(R), expr_new_symbol(su));
    Expr* coefs[6] = { A, B, C, Dc, Ec, Fc };

    /* linear: every coefficient free of all s-symbols */
    bool ok = true;
    for (int c = 0; c < 6 && ok; c++)
        for (int s = 0; s < 6 && ok; s++)
            if (!ds_free_of(coefs[c], svars[s])) ok = false;

    /* constant coefficients (every coefficient free of v1,v2); lower-order terms
     * D u_{v1} + E u_{v2} + F u are now allowed (handled by operator factoring). */
    if (ok) ok = ds_free_of(A, v1) && ds_free_of(A, v2)
              && ds_free_of(B, v1) && ds_free_of(B, v2)
              && ds_free_of(C, v1) && ds_free_of(C, v2)
              && ds_free_of(Dc, v1) && ds_free_of(Dc, v2)
              && ds_free_of(Ec, v1) && ds_free_of(Ec, v2)
              && ds_free_of(Fc, v1) && ds_free_of(Fc, v2);

    /* zero forcing (homogeneous): R with all s-symbols -> 0 must vanish */
    if (ok) {
        Expr* R0 = expr_copy(R);
        for (int s = 0; s < 6; s++)
            R0 = ds_subst(R0, expr_new_symbol(svars[s]), expr_new_integer(0));
        if (!ds_is_zero(R0)) ok = false;
        expr_free(R0);
    }

    /* genuinely second order: at least one of A, B, C nonzero */
    if (ok && ds_is_zero(A) && ds_is_zero(B) && ds_is_zero(C)) ok = false;

    expr_free(R);
    if (!ok) {
        for (int i = 0; i < 6; i++) expr_free(coefs[i]);
        return NULL;
    }

    Expr* body;
    if (!ds_is_zero(A)) {
        body = pde2_body(A, B, C, Dc, Ec, Fc, v1, v2);
    } else if (!ds_is_zero(C)) {
        /* A == 0: swap v1<->v2 (and A<->C, D<->E), trial f(v1 + l v2) */
        body = pde2_body(C, B, A, Ec, Dc, Fc, v2, v1);
    } else {
        /* pure mixed  B u_{v1 v2} + D u_{v1} + E u_{v2} + F u == 0 */
        body = pde2_body_mixed(B, Dc, Ec, Fc, v1, v2);
    }
    for (int i = 0; i < 6; i++) expr_free(coefs[i]);
    if (!body) return NULL;

    Expr** bodies = malloc(sizeof(Expr*));
    bodies[0] = body;
    (void)su;
    return bodies;
}

static Expr* builtin_dsolve_pde2(Expr* res) {
    return dsolve_method_builtin_pde(res, dsolve_pde2_solve);
}

void dsolve_pde2_init(void) {
    symtab_add_builtin("DSolve`PDELinearSecondOrder", builtin_dsolve_pde2);
    symtab_get_def("DSolve`PDELinearSecondOrder")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`PDELinearSecondOrder",
        "DSolve`PDELinearSecondOrder[eqn, u, {v1, v2}] solves a homogeneous, "
        "constant-coefficient, principal-part-only second-order linear PDE "
        "A u_{v1 v1} + B u_{v1 v2} + C u_{v2 v2} == 0 by operator factoring: the "
        "trial u == f(v2 + lambda v1) gives the characteristic quadratic "
        "A lambda^2 + B lambda + C == 0.  Distinct real roots (hyperbolic) give "
        "u == C[1][v2 + l1 v1] + C[2][v2 + l2 v1] (the wave equation u_tt == c^2 u_xx "
        "-> C[1][x - c t] + C[2][x + c t], d'Alembert); complex roots (elliptic) give "
        "the complex-characteristic form (Laplace u_xx + u_yy == 0 -> "
        "C[1][y - I x] + C[2][y + I x]); a repeated root (parabolic) gives "
        "C[1][w] + v1 C[2][w], w == v2 + lambda v1.  Lower-order terms "
        "D u_{v1} + E u_{v2} + F u are handled when the full symbol factors into two "
        "first-order operators, giving the exponential-damped form "
        "e^{-m1 v1} C[1][.] + e^{-m2 v1} C[2][.] (e.g. the distortionless telegraph "
        "u_tt - c^2 u_xx + a u_t + a^2/4 u == 0 -> E^(-a t/2)(C[1][x+c t]+C[2][x-c t])).  "
        "Declines a non-factorable symbol (the general telegraph needs Bessel "
        "functions), forcing, and non-constant coefficients.");
}
