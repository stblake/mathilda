/*
 * dsolve_normalform.c — DSolve`NormalForm.
 *
 * Not a solver: an inspection utility that exposes the substrate helper
 * dsolve_normal_form for a homogeneous second-order linear ODE.  For
 *     y'' + P(x) y' + Q(x) y == 0
 * the substitution  y == w z  with  w == Exp[-Integrate[P/2, x]]  removes the
 * first-derivative term, leaving the reduced (Liouville normal) form
 *     z'' == r z,   r == P^2/4 + P'/2 - Q.
 * DSolve`NormalForm[eqn, y, x] returns the pair {r, w}.  This is the prerequisite
 * Kovacic and the special-function recognizers operate on, surfaced as a
 * REPL-callable builtin so the reduction can be checked directly.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../symtab.h"
#include "../attr.h"
#include <stdlib.h>

static Expr* builtin_dsolve_normalform(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count < 3) return NULL;

    DSolveProblem P;
    if (!dsolve_parse(res, &P)) { dsolve_problem_free(&P); return NULL; }

    Expr* Pc; Expr* Qc;
    if (!dsolve_second_order_PQ(&P, &Pc, &Qc)) { dsolve_problem_free(&P); return NULL; }

    const char* xvar = P.ind_names[0];
    Expr* recovery = NULL;
    Expr* r = dsolve_normal_form(Pc, Qc, xvar, &recovery);
    if (!recovery) {
        /* the integral was not elementary: return the honest unevaluated factor */
        Expr* half = eval_and_free(ds_call2(SYM_Times, expr_copy(Pc),
                         expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ expr_new_integer(2), expr_new_integer(-1) }, 2)));
        Expr* integ = ds_integrate(half, expr_new_symbol(xvar));
        recovery = eval_and_free(ds_call1("Exp",
                       eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), integ))));
    }
    expr_free(Pc); expr_free(Qc);
    dsolve_problem_free(&P);

    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ r, recovery }, 2);
}

void dsolve_normalform_init(void) {
    symtab_add_builtin("DSolve`NormalForm", builtin_dsolve_normalform);
    symtab_get_def("DSolve`NormalForm")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`NormalForm",
        "DSolve`NormalForm[eqn, y, x] returns {r, w} for the second-order linear "
        "ODE y'' + P y' + Q y == 0: the substitution y == w z with "
        "w == Exp[-Integrate[P/2, x]] removes the first-derivative term, giving the "
        "reduced form z'' == r z with r == P^2/4 + P'/2 - Q.");
}
