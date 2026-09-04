/*
 * dsolve_pdecharpit.c — DSolve`PDECharpit (Charpit's method, standard forms).
 *
 * Solves a first-order FULLY NONLINEAR PDE
 *
 *     F(v1, v2, u, p, q) == 0,     p = u_{v1},  q = u_{v2}
 *
 * for a complete integral (a two-parameter family), by Charpit's method
 * restricted to the three classic deterministic standard forms.  Runs after
 * PDEQuasilinear / PDEClairaut in the cascade, so it sees only the forms
 * nonlinear in the derivatives that those declined.
 *
 *   Type I  — F(p,q)  (free of v1, v2, u):  set p = C[1], solve F(C[1], q) == 0
 *             for q; complete integral (explicit)  u = C[1] v1 + q v2 + C[2].
 *   Type II — F(u,p,q)  (free of v1, v2):  Charpit gives q = C[1] p; solve
 *             F(u, p, C[1] p) == 0 for p = P(u, C[1]); complete integral
 *             (implicit)  Integrate[1/P, u] == v1 + C[1] v2 + C[2].
 *   Type III— separable f(v1,p) == g(v2,q)  (free of u; the (v1,p) group and the
 *             (v2,q) group decouple — all four mixed second partials vanish):
 *             solve f == C[1] for p = P, g == (F0 - C[1]) for q = Q; complete
 *             integral (explicit)  u = Integrate[P, v1] + Integrate[Q, v2] + C[2].
 *
 * Types I & III are explicit (routed through the PDEBranches wrapper — bare
 * constants C[1], C[2] survive the ordinary explicit PDE verify); Type II is the
 * implicit relation Psi == C[2] (the PDERelation wrapper, verified by the
 * implicit-function rule).  A non-elementary Integrate makes the method decline.
 *
 * First cut (honest): the three standard forms; the general Charpit
 * integrable-combination search and singular (envelope) solutions are future.
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
static Expr* sub(Expr* a, Expr* b) { return eval_and_free(ds_call2(SYM_Subtract, a, b)); }
static Expr* inv(Expr* a) { return expr_new_function(expr_new_symbol(SYM_Power),
                                (Expr*[]){ a, expr_new_integer(-1) }, 2); }
static Expr* eq0(Expr* a) { return expr_new_function(expr_new_symbol(SYM_Equal),
                                (Expr*[]){ a, expr_new_integer(0) }, 2); }

/* Derivative[o1,o2][u][v1,v2]. */
static Expr* pde_deriv(const char* u, int o1, int o2, const char* v1, const char* v2) {
    Expr* d = expr_new_function(expr_new_symbol(SYM_Derivative),
                  (Expr*[]){ expr_new_integer(o1), expr_new_integer(o2) }, 2);
    Expr* du = expr_new_function(d, (Expr*[]){ expr_new_symbol(u) }, 1);
    return expr_new_function(du, (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2);
}

/* wrap a single body/relation for the run_pde_implicit dispatcher */
static Expr** wrap1(const char* head, Expr** args, size_t n) {
    Expr* w = expr_new_function(expr_new_symbol(intern_symbol(head)), args, n);
    Expr** out = malloc(sizeof(Expr*));
    out[0] = w;
    return out;
}

/* Type I: F(p,q).  p = C[1], solve F(C[1], q) == 0 for q; branches
 * u = C[1] v1 + q v2 + C[2]. */
static Expr** charpit_type1(const Expr* F, const char* sp, const char* sq,
                            const char* v1, const char* v2) {
    Expr* Fc = ds_subst(expr_copy((Expr*)F), expr_new_symbol(sp), ds_const(1));
    Expr* sol = ds_solve(eq0(Fc), expr_new_symbol(sq));
    size_t nq = 0;
    Expr** qs = dsolve_extract_solutions(sol, sq, &nq);
    if (sol) expr_free(sol);
    if (!qs || nq == 0) { if (qs) free(qs); return NULL; }

    Expr** branches = malloc(nq * sizeof(Expr*));
    for (size_t i = 0; i < nq; i++)
        branches[i] = add(add(mul(ds_const(1), expr_new_symbol(v1)),
                              mul(qs[i], expr_new_symbol(v2))), ds_const(2));
    free(qs);
    Expr** out = wrap1("DSolve`PDEBranches", branches, nq);  /* copies the array */
    free(branches);
    return out;
}

/* Type II: F(u,p,q).  q = C[1] p, solve F(u,p,C[1] p) == 0 for p = P(u,C[1]);
 * implicit relation  Integrate[1/P, u] - v1 - C[1] v2 == C[2]. */
static Expr** charpit_type2(const Expr* F, const char* su, const char* sp,
                            const char* sq, const char* v1, const char* v2) {
    Expr* Fc = ds_subst(expr_copy((Expr*)F), expr_new_symbol(sq),
                        mul(ds_const(1), expr_new_symbol(sp)));   /* q -> C[1] p */
    Expr* sol = ds_solve(eq0(Fc), expr_new_symbol(sp));
    size_t np = 0;
    Expr** ps = dsolve_extract_solutions(sol, sp, &np);
    if (sol) expr_free(sol);
    if (!ps || np == 0) { if (ps) free(ps); return NULL; }
    Expr* P = expr_copy(ps[0]);
    for (size_t i = 0; i < np; i++) expr_free(ps[i]);
    free(ps);

    Expr* Iu = ds_integrate(inv(P), expr_new_symbol(su));         /* Integrate[1/P, u] */
    if (ds_has_head(Iu, SYM_Integrate)) { expr_free(Iu); return NULL; }
    /* Psi = Iu - v1 - C[1] v2 */
    Expr* Psi = sub(sub(Iu, expr_new_symbol(v1)),
                    mul(ds_const(1), expr_new_symbol(v2)));
    return wrap1("DSolve`PDERelation", (Expr*[]){ Psi }, 1);
}

/* Type III: separable f(v1,p) == g(v2,q).  Solve f == C[1] for p, g == (F0 - C[1])
 * for q (F0 = F at the zero reference); u = Integrate[P,v1] + Integrate[Q,v2] + C[2]. */
static Expr** charpit_type3(const Expr* F, const char* v1, const char* v2,
                            const char* sp, const char* sq) {
    /* separability: the four mixed (v1|sp)-(v2|sq) second partials vanish */
    const char* g1[2] = { v1, sp };
    const char* g2[2] = { v2, sq };
    for (int a = 0; a < 2; a++)
        for (int b = 0; b < 2; b++) {
            Expr* d2 = ds_d(ds_d(expr_copy((Expr*)F), expr_new_symbol(g1[a])),
                            expr_new_symbol(g2[b]));
            bool z = ds_is_zero(d2);
            expr_free(d2);
            if (!z) return NULL;
        }
    /* f_part = F | v2=0, sq=0 ;  g_part = F | v1=0, sp=0 ;  F0 = F | all=0 */
    Expr* fpart = ds_subst(ds_subst(expr_copy((Expr*)F), expr_new_symbol(v2), expr_new_integer(0)),
                           expr_new_symbol(sq), expr_new_integer(0));
    Expr* gpart = ds_subst(ds_subst(expr_copy((Expr*)F), expr_new_symbol(v1), expr_new_integer(0)),
                           expr_new_symbol(sp), expr_new_integer(0));
    Expr* F0 = ds_subst(ds_subst(expr_copy(gpart), expr_new_symbol(v2), expr_new_integer(0)),
                        expr_new_symbol(sq), expr_new_integer(0));

    /* Solve fpart == C[1] for p,  gpart == F0 - C[1] for q */
    Expr* Psol = ds_solve(expr_new_function(expr_new_symbol(SYM_Equal),
                     (Expr*[]){ fpart, ds_const(1) }, 2), expr_new_symbol(sp));
    Expr* Qsol = ds_solve(expr_new_function(expr_new_symbol(SYM_Equal),
                     (Expr*[]){ gpart, sub(F0, ds_const(1)) }, 2), expr_new_symbol(sq));
    size_t nP = 0, nQ = 0;
    Expr** Ps = dsolve_extract_solutions(Psol, sp, &nP);
    Expr** Qs = dsolve_extract_solutions(Qsol, sq, &nQ);
    if (Psol) expr_free(Psol);
    if (Qsol) expr_free(Qsol);
    Expr** out = NULL;
    if (Ps && nP > 0 && Qs && nQ > 0) {
        Expr* Pint = ds_integrate(expr_copy(Ps[0]), expr_new_symbol(v1));
        Expr* Qint = ds_integrate(expr_copy(Qs[0]), expr_new_symbol(v2));
        if (!ds_has_head(Pint, SYM_Integrate) && !ds_has_head(Qint, SYM_Integrate)) {
            Expr* body = add(add(Pint, Qint), ds_const(2));
            out = wrap1("DSolve`PDEBranches", (Expr*[]){ body }, 1);
        } else { expr_free(Pint); expr_free(Qint); }
    }
    if (Ps) { for (size_t i = 0; i < nP; i++) expr_free(Ps[i]); free(Ps); }
    if (Qs) { for (size_t i = 0; i < nQ; i++) expr_free(Qs[i]); free(Qs); }
    return out;
}

Expr** dsolve_pdecharpit_solve(DSolveProblem* P) {
    if (P->nfun != 1 || P->nind != 2 || P->neq != 1) return NULL;
    const char* uname = P->fun_names[0];
    const char* v1 = P->ind_names[0];
    const char* v2 = P->ind_names[1];
    const char* sp = intern_symbol("DSolve`pdeP");
    const char* sq = intern_symbol("DSolve`pdeQ");
    const char* su = intern_symbol("DSolve`pdeU");

    /* F = residual with u_{v1}->p, u_{v2}->q, u->u (bare symbols) */
    Expr* F = expr_copy(P->eq_residuals[0]);
    F = ds_subst(F, pde_deriv(uname, 1, 0, v1, v2), expr_new_symbol(sp));
    F = ds_subst(F, pde_deriv(uname, 0, 1, v1, v2), expr_new_symbol(sq));
    F = ds_subst(F, expr_new_function(expr_new_symbol(uname),
                     (Expr*[]){ expr_new_symbol(v1), expr_new_symbol(v2) }, 2),
                 expr_new_symbol(su));

    /* nonlinear in the derivatives: else quasilinear/linear owns it */
    Expr* Fp = ds_d(expr_copy(F), expr_new_symbol(sp));
    Expr* Fq = ds_d(expr_copy(F), expr_new_symbol(sq));
    bool linear = ds_free_of(Fp, sp) && ds_free_of(Fp, sq)
               && ds_free_of(Fq, sp) && ds_free_of(Fq, sq);
    expr_free(Fp); expr_free(Fq);
    if (linear) { expr_free(F); return NULL; }

    bool no_v1 = ds_free_of(F, v1), no_v2 = ds_free_of(F, v2), no_u = ds_free_of(F, su);

    Expr** out = NULL;
    if (no_v1 && no_v2 && no_u)      out = charpit_type1(F, sp, sq, v1, v2);
    else if (no_v1 && no_v2)         out = charpit_type2(F, su, sp, sq, v1, v2);
    else if (no_u)                   out = charpit_type3(F, v1, v2, sp, sq);
    expr_free(F);
    return out;
}

static Expr* builtin_dsolve_pdecharpit(Expr* res) {
    return dsolve_method_builtin_pde_implicit(res, dsolve_pdecharpit_solve);
}

void dsolve_pdecharpit_init(void) {
    symtab_add_builtin("DSolve`PDECharpit", builtin_dsolve_pdecharpit);
    symtab_get_def("DSolve`PDECharpit")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`PDECharpit",
        "DSolve`PDECharpit[eqn, u, {v1, v2}] solves a first-order fully nonlinear PDE "
        "F(v1, v2, u, u_{v1}, u_{v2}) == 0 for a complete integral by Charpit's method, "
        "restricted to the three standard forms: F(p,q) (u == C[1] v1 + q v2 + C[2], q "
        "from F(C[1],q)==0); F(u,p,q) (implicit Integrate[1/P,u] == v1 + C[1] v2 + C[2], "
        "P from q==C[1] p); and separable f(v1,p) == g(v2,q) (u == Integrate[P,v1] + "
        "Integrate[Q,v2] + C[2]).  Nonlinear in the derivatives, so the "
        "linear/quasilinear/Clairaut methods decline to it.  Declines equations "
        "depending on all of v1, v2, u (the general integrable-combination search is "
        "future) and non-elementary integrals.");
}
