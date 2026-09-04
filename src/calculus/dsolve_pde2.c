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
 * First cut (honest): homogeneous, principal-part only (no u_{v1}, u_{v2}, u
 * terms; zero forcing), constant coefficients.  Lower-order terms (telegraph /
 * damped wave — the full symbol must factor into first-order operators),
 * inhomogeneous forcing, and the wave-IVP d'Alembert formula are future work.
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

/* Build the general body from the characteristic quadratic A l^2 + B l + C = 0
 * (trial u = f(v2 + l v1)).  A != 0 assumed.  Returns NULL if the quadratic does
 * not fully factor into two roots (total multiplicity 2). */
static Expr* pde2_body_from_quadratic(const Expr* A, const Expr* B, const Expr* C,
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

    Expr* body;
    if (R.ndist == 2) {
        /* distinct roots: C[1][v2 + l1 v1] + C[2][v2 + l2 v1] */
        Expr* a1 = add(expr_new_symbol(v2), mul(expr_copy(R.roots[0]), expr_new_symbol(v1)));
        Expr* a2 = add(expr_new_symbol(v2), mul(expr_copy(R.roots[1]), expr_new_symbol(v1)));
        body = add(arbfun(1, a1), arbfun(2, a2));
    } else {
        /* repeated root: w = v2 + l v1;  C[1][w] + v1 C[2][w] */
        Expr* w1 = add(expr_new_symbol(v2), mul(expr_copy(R.roots[0]), expr_new_symbol(v1)));
        Expr* w2 = add(expr_new_symbol(v2), mul(expr_copy(R.roots[0]), expr_new_symbol(v1)));
        body = add(arbfun(1, w1), mul(expr_new_symbol(v1), arbfun(2, w2)));
    }
    dsolve_roots_free(&R);
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

    /* constant coefficients (A,B,C free of v1,v2); homogeneous principal part
     * (no lower-order terms: Dc = Ec = Fc = 0) */
    if (ok) ok = ds_free_of(A, v1) && ds_free_of(A, v2)
              && ds_free_of(B, v1) && ds_free_of(B, v2)
              && ds_free_of(C, v1) && ds_free_of(C, v2)
              && ds_is_zero(Dc) && ds_is_zero(Ec) && ds_is_zero(Fc);

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
    expr_free(Dc); expr_free(Ec); expr_free(Fc);

    Expr* body;
    if (!ds_is_zero(A)) {
        body = pde2_body_from_quadratic(A, B, C, v1, v2);
    } else if (!ds_is_zero(C)) {
        /* A == 0: swap v1<->v2 and A<->C, trial f(v1 + l v2) */
        body = pde2_body_from_quadratic(C, B, A, v2, v1);
    } else {
        /* pure mixed  B u_{v1 v2} == 0  ->  u = F(v1) + G(v2) */
        body = add(arbfun(1, expr_new_symbol(v1)), arbfun(2, expr_new_symbol(v2)));
    }
    expr_free(A); expr_free(B); expr_free(C);
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
        "C[1][w] + v1 C[2][w], w == v2 + lambda v1.  The generated arbitrary functions "
        "are C[1][.], C[2][.].  Declines lower-order terms, forcing, and non-constant "
        "coefficients.");
}
