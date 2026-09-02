/*
 * dsolve_autonomous.c — DSolve`AutonomousReduction.
 *
 * Solves the second-order autonomous ODE (independent variable x absent)
 *     y''[x] == f(y, y')
 * by the standard order reduction p = y' regarded as a function of y:
 *     y'' = dp/dx = (dp/dy)(dy/dx) = p p'(y),
 * turning the equation into the first-order ODE in p(y)
 *     p p'(y) == f(y, p),
 * which is solved by recursing into the scalar cascade.  With p = P(y, C[1]) the
 * remaining equation y' == P(y) is itself autonomous (hence separable) and is
 * solved by a second recursion, giving y(x, C[1], C[2]).
 *
 * The two integration constants must stay distinct across the two recursions
 * (each DSolve call independently names its constant C[1]), so the first stage's
 * constant is frozen to a private marker before the second solve and renamed to
 * C[2] afterwards.  A degenerate y = const (which trivially satisfies any
 * autonomous equation and would pass back-substitution — see the M4 lesson) is
 * rejected by requiring the final body to depend on x.
 *
 * Runs after ReductionOfOrder (missing-y) and the linear methods in the cascade,
 * so genuinely linear/constant-coefficient autonomous equations (y'' == -y, …)
 * are claimed by the cleaner method first; this one catches the nonlinear
 * remainder such as y y'' == (y')^2  →  y = C[2] E^(C[1] x).
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>

/* Solve DSolve[eqn, fname[varname], varname] and return the applied-form RHS
 * body (owned), or NULL if the sub-solve declined.  `eqn` is consumed. */
static Expr* run_dsolve_applied(Expr* eqn, const char* fname, const char* varname) {
    Expr* lhs  = ds_call1(fname, expr_new_symbol(varname));       /* fname[varname] */
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ eqn, lhs, expr_new_symbol(varname) }, 3);
    Expr* r = eval_and_free(call);
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

Expr** dsolve_autonomous_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 2) return NULL;
    const char* xvar  = P->ind_names[0];
    const char* yname = P->fun_names[0];

    Expr* F = dsolve_solve_top_derivative(P, 2);          /* y'' == F(x, y, y') */
    if (!F) return NULL;

    /* Replace the funcapps by plain markers to test for explicit x and for y. */
    const char* Ysym = intern_symbol("DSolve`arY");
    const char* Psym = intern_symbol("DSolve`arP");
    Expr* Ftmp = ds_subst(expr_copy(F), ds_make_funcapp(yname, 1, xvar), expr_new_symbol(Psym));
    Ftmp = ds_subst(Ftmp, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Ysym));
    expr_free(F);
    /* autonomous: no explicit x; and must genuinely involve y (else it is the
     * missing-y case ReductionOfOrder owns). */
    if (!ds_free_of(Ftmp, xvar) || !ds_contains(Ftmp, Ysym)) { expr_free(Ftmp); return NULL; }

    /* Stage 1: p p'(y) == f(y, p).  Independent variable is y (= Ysym), the
     * dependent function is pfun; y' (Psym) becomes pfun[Ysym]. */
    const char* pfun = intern_symbol("DSolve`arp");
    Expr* rhs1 = ds_subst(Ftmp, expr_new_symbol(Psym), ds_make_funcapp(pfun, 0, Ysym));
    Expr* lhs1 = ds_call2(SYM_Times, ds_make_funcapp(pfun, 0, Ysym), ds_make_funcapp(pfun, 1, Ysym));
    Expr* eq1  = expr_new_function(expr_new_symbol(SYM_Equal), (Expr*[]){ lhs1, rhs1 }, 2);
    Expr* pbody = run_dsolve_applied(eq1, pfun, Ysym);
    if (!pbody) return NULL;

    /* Freeze the first constant as the GENERATED constant C[2] (distinct from the
     * second solve's fresh C[1], so they cannot collide).  It must NOT be frozen to
     * a plain symbol: a plain symbolic parameter inside the stage-2 integrand sends
     * the first-order cascade down a ~27x slower path (DSolve[y'==Sqrt[a+y^4]] ~ 11s
     * vs DSolve[y'==Sqrt[C[2]+y^4]] ~ 0.4s), and stage 2 is exactly the elliptic /
     * quadrature form autonomous reduction produces.  C[2] is recognised as a
     * constant, taking the fast decline, and is already the final name (no rename). */
    pbody = ds_subst(pbody, ds_const(1), ds_const(2));

    /* Stage 2: y'[x] == P(y[x]) (autonomous, separable). */
    Expr* pOfY = ds_subst(pbody, expr_new_symbol(Ysym), ds_make_funcapp(yname, 0, xvar));
    Expr* eq2  = expr_new_function(expr_new_symbol(SYM_Equal),
                     (Expr*[]){ ds_make_funcapp(yname, 1, xvar), pOfY }, 2);
    Expr* ybody = run_dsolve_applied(eq2, yname, xvar);
    if (!ybody) return NULL;

    /* Reject a degenerate x-independent body. */
    if (ds_free_of(ybody, xvar)) { expr_free(ybody); return NULL; }

    Expr** out = malloc(sizeof(Expr*));
    out[0] = ybody;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_autonomous(Expr* res) {
    return dsolve_method_builtin(res, dsolve_autonomous_try);
}

void dsolve_autonomous_init(void) {
    symtab_add_builtin("DSolve`AutonomousReduction", builtin_dsolve_autonomous);
    symtab_get_def("DSolve`AutonomousReduction")->attributes
        |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`AutonomousReduction",
        "DSolve`AutonomousReduction[eqn, y, x] solves a second-order autonomous ODE "
        "y'' == f(y, y') (no explicit x) by the reduction p = y'(y), p p'(y) == f, "
        "then solves y' == p(y) by separation.");
}
