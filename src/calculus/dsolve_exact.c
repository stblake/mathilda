/*
 * dsolve_exact.c — DSolve`Exact.
 *
 * Solves  M(x,y) + N(x,y) y' == 0  when it is exact (M_y == N_x), or can be made
 * exact by an integrating factor mu(x) or mu(y):
 *   - (M_y - N_x)/N free of y  ->  mu(x) = Exp[Integrate[(M_y-N_x)/N, x]]
 *   - (N_x - M_y)/M free of x  ->  mu(y) = Exp[Integrate[(N_x-M_y)/M, y]]
 * The potential F with F_x = M, F_y = N is F = Integrate[M, x] + g(y) where
 * g'(y) = N - d/dy Integrate[M, x]; the solution is F(x,y) == C[1], solved for y.
 *
 * Works on the algebraic residual R(x, Y, p) (p = y'), which must be linear in p:
 * N = R_p, M = R|p=0.  Every returned branch is still back-substitution verified
 * by the substrate, so an imperfect potential construction cannot ship a wrong
 * answer — it is dropped.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

static Expr* pw_inv(Expr* base) {
    return expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ base, expr_new_integer(-1) }, 2);
}

/* True iff `e` is a polynomial in x and Y. */
static bool poly_in_xy(const Expr* e, const char* xvar, const char* Yn) {
    Expr* vars = expr_new_function(expr_new_symbol(SYM_List),
                     (Expr*[]){ expr_new_symbol(xvar), expr_new_symbol(Yn) }, 2);
    Expr* q = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
                  (Expr*[]){ expr_copy((Expr*)e), vars }, 2));
    bool r = (q->type == EXPR_SYMBOL && q->data.symbol.name == SYM_True);
    expr_free(q);
    return r;
}

/* Propose an integrating factor mu = x^a y^b (constant exponents) and accept it
 * only after a symbolic exactness re-check.  For M + N y' == 0 the exactness
 * condition under mu = x^a y^b is
 *     b (M/y) - a (N/x) + (M_y - N_x) == 0   (identically in x, y),
 * which is linear in (a, b).  We propose (a, b) by sampling that identity at a
 * few generic non-zero points (a heuristic), then require (mu M)_y == (mu N)_x
 * symbolically before returning mu; the substrate's back-substitution verify is
 * the final backstop, so a wrong proposal is dropped, never shipped.  Solves
 * e.g. (x y - 2 x) y' == y - y^2 + 3 x^2 y^3  (mu = x^-2 y^-3).  M, N, diff
 * (= M_y - N_x) are borrowed; returns mu (owned) or NULL. */
static Expr* exact_xayb_factor(const Expr* M, const Expr* N, const Expr* diff,
                               const char* xvar, const char* Yn) {
    /* x^a y^b is a polynomial-ODE tool: gate on M, N being polynomials in x, y so
     * the speculative Solve below never runs on the heavy radical/transcendental
     * systems the AutonomousReduction recursion feeds into Exact (that Solve on
     * Sqrt-coefficient systems is what made a whole test run blow past its
     * budget). */
    if (!poly_in_xy(M, xvar, Yn) || !poly_in_xy(N, xvar, Yn)) return NULL;

    const char* aSym = intern_symbol("DSolve`exA");
    const char* bSym = intern_symbol("DSolve`exB");

    /* cond = a*(-N/x) + b*(M/y) + diff */
    Expr* Pc = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
                   expr_new_integer(-1), expr_copy((Expr*)N), pw_inv(expr_new_symbol(xvar)) }, 3));
    Expr* Qc = eval_and_free(ds_call2(SYM_Times, expr_copy((Expr*)M), pw_inv(expr_new_symbol(Yn))));
    Expr* cond = eval_and_free(ds_call2(SYM_Plus,
        ds_call2(SYM_Times, expr_new_symbol(aSym), Pc),          /* consumes Pc */
        ds_call2(SYM_Plus,
            ds_call2(SYM_Times, expr_new_symbol(bSym), Qc),      /* consumes Qc */
            expr_copy((Expr*)diff))));

    /* sample at generic non-zero points -> linear system {cond(pt_i)==0} in (a,b) */
    static const int sx[] = { 2, 3, 5 };
    static const int sy[] = { 3, 5, 7 };
    Expr* eqs[3];
    for (int i = 0; i < 3; i++) {
        Expr* ci = ds_subst(
            ds_subst(expr_copy(cond), expr_new_symbol(xvar), expr_new_integer(sx[i])),
            expr_new_symbol(Yn), expr_new_integer(sy[i]));
        eqs[i] = expr_new_function(expr_new_symbol(SYM_Equal),
                     (Expr*[]){ ci, expr_new_integer(0) }, 2);
    }
    expr_free(cond);
    Expr* eqlist  = expr_new_function(expr_new_symbol(SYM_List), eqs, 3);
    Expr* varlist = expr_new_function(expr_new_symbol(SYM_List),
                        (Expr*[]){ expr_new_symbol(aSym), expr_new_symbol(bSym) }, 2);
    Expr* solres = ds_solve(eqlist, varlist);      /* consumes eqlist, varlist */

    /* extract a, b from the first solution {{a -> ..., b -> ...}} */
    Expr* aval = NULL; Expr* bval = NULL;
    if (solres && solres->type == EXPR_FUNCTION && solres->data.function.arg_count >= 1) {
        Expr* first = solres->data.function.args[0];
        if (first && first->type == EXPR_FUNCTION) {
            for (size_t k = 0; k < first->data.function.arg_count; k++) {
                Expr* rule = first->data.function.args[k];
                if (rule && rule->type == EXPR_FUNCTION && rule->data.function.arg_count == 2
                    && rule->data.function.head->type == EXPR_SYMBOL
                    && (rule->data.function.head->data.symbol.name == SYM_Rule
                        || rule->data.function.head->data.symbol.name == SYM_RuleDelayed)) {
                    Expr* lhs = rule->data.function.args[0];
                    Expr* rhs = rule->data.function.args[1];
                    if (lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == aSym && !aval) aval = expr_copy(rhs);
                    else if (lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == bSym && !bval) bval = expr_copy(rhs);
                }
            }
        }
    }
    if (solres) expr_free(solres);
    if (!aval || !bval) { if (aval) expr_free(aval); if (bval) expr_free(bval); return NULL; }

    /* mu = x^a y^b */
    Expr* mu = eval_and_free(ds_call2(SYM_Times,
        eval_and_free(ds_call2(SYM_Power, expr_new_symbol(xvar), aval)),   /* consumes aval */
        eval_and_free(ds_call2(SYM_Power, expr_new_symbol(Yn), bval))));   /* consumes bval */

    /* symbolic exactness re-check: (mu M)_y - (mu N)_x == 0 */
    Expr* mM = eval_and_free(ds_call2(SYM_Times, expr_copy(mu), expr_copy((Expr*)M)));
    Expr* mN = eval_and_free(ds_call2(SYM_Times, expr_copy(mu), expr_copy((Expr*)N)));
    Expr* chk = eval_and_free(ds_call2(SYM_Subtract,
        ds_d(mM, expr_new_symbol(Yn)),          /* consumes mM */
        ds_d(mN, expr_new_symbol(xvar))));      /* consumes mN */
    bool exact = ds_is_zero(chk);
    expr_free(chk);
    if (!exact) { expr_free(mu); return NULL; }
    return mu;
}

Expr** dsolve_exact_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* xvar = P->ind_names[0];
    const char* Yn = intern_symbol("DSolve`Y");
    const char* Pn = intern_symbol("DSolve`p");

    Expr* R = dsolve_algebraic_residual(P, Yn, Pn);
    if (!R) return NULL;

    /* R = M + N p, N = R_p free of p, M = R|p=0 */
    Expr* N = ds_d(expr_copy(R), expr_new_symbol(Pn));
    if (!ds_free_of(N, Pn)) { expr_free(N); expr_free(R); return NULL; }
    Expr* M = ds_subst(expr_copy(R), expr_new_symbol(Pn), expr_new_integer(0));
    Expr* recon = eval_and_free(ds_call2(SYM_Plus, expr_copy(M),
                        ds_call2(SYM_Times, expr_copy(N), expr_new_symbol(Pn))));
    Expr* lin = eval_and_free(ds_call2(SYM_Subtract, expr_copy(R), recon));
    bool linear = ds_is_zero(lin);
    expr_free(lin); expr_free(R);
    if (!linear) { expr_free(M); expr_free(N); return NULL; }

    /* diff = M_y - N_x */
    Expr* My = ds_d(expr_copy(M), expr_new_symbol(Yn));
    Expr* Nx = ds_d(expr_copy(N), expr_new_symbol(xvar));
    Expr* diff = eval_and_free(ds_call2(SYM_Subtract, My, Nx));

    Expr* Mu = NULL;
    if (ds_is_zero(diff)) {
        Mu = expr_new_integer(1);
    } else {
        /* mu(x): (M_y - N_x)/N free of Y */
        Expr* r1 = eval_and_free(ds_call2(SYM_Times, expr_copy(diff), pw_inv(expr_copy(N))));
        if (ds_free_of(r1, Yn)) {
            Expr* r1int = ds_integrate(r1, expr_new_symbol(xvar));
            if (!ds_has_head(r1int, SYM_Integrate)) Mu = eval_and_free(ds_call1("Exp", r1int));
            else expr_free(r1int);
        } else {
            expr_free(r1);
            /* mu(y): (N_x - M_y)/M free of x */
            Expr* r2 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
                expr_new_integer(-1), expr_copy(diff), pw_inv(expr_copy(M)) }, 3));
            if (ds_free_of(r2, xvar)) {
                Expr* r2int = ds_integrate(r2, expr_new_symbol(Yn));
                if (!ds_has_head(r2int, SYM_Integrate)) Mu = eval_and_free(ds_call1("Exp", r2int));
                else expr_free(r2int);
            } else expr_free(r2);
        }
    }
    /* fallback: mu = x^a y^b (constant exponents) when neither mu(x) nor mu(y) fits */
    if (!Mu) Mu = exact_xayb_factor(M, N, diff, xvar, Yn);
    expr_free(diff);
    if (!Mu) { expr_free(M); expr_free(N); return NULL; }

    /* Mm = Mu M, Nn = Mu N */
    Expr* Mm = eval_and_free(ds_call2(SYM_Times, expr_copy(Mu), M));   /* consumes M */
    Expr* Nn = eval_and_free(ds_call2(SYM_Times, expr_copy(Mu), N));   /* consumes N */
    expr_free(Mu);

    /* F = Integrate[Mm, x] + g(Y),  g'(Y) = Nn - d/dY Integrate[Mm, x] */
    Expr* Fx = ds_integrate(Mm, expr_new_symbol(xvar));
    if (ds_has_head(Fx, SYM_Integrate)) { expr_free(Fx); expr_free(Nn); return NULL; }
    Expr* gp = eval_and_free(ds_call2(SYM_Subtract, Nn, ds_d(expr_copy(Fx), expr_new_symbol(Yn))));
    Expr* g = ds_integrate(gp, expr_new_symbol(Yn));
    if (ds_has_head(g, SYM_Integrate)) { expr_free(g); expr_free(Fx); return NULL; }
    Expr* Fpot = eval_and_free(ds_call2(SYM_Plus, Fx, g));

    /* solution F(x, Y) == C[1], solved for Y */
    Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ Fpot, ds_const(1) }, 2);
    Expr* solres = ds_solve(eq, expr_new_symbol(Yn));
    size_t nb = 0;
    Expr** bodies = dsolve_extract_solutions(solres, Yn, &nb);
    if (solres) expr_free(solres);
    if (!bodies) return NULL;                 /* only explicit solutions in M1 */
    *nbranch = nb;
    return bodies;
}

static Expr* builtin_dsolve_exact(Expr* res) {
    return dsolve_method_builtin(res, dsolve_exact_try);
}

void dsolve_exact_init(void) {
    symtab_add_builtin("DSolve`Exact", builtin_dsolve_exact);
    symtab_get_def("DSolve`Exact")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`Exact",
        "DSolve`Exact[eqn, y, x] solves M(x,y) + N(x,y) y' == 0 when exact "
        "(M_y == N_x), or made exact by an integrating factor mu(x) or mu(y).");
}
