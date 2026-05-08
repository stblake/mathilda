/* root.c
 *
 * Mathematica-style symbolic Root and RootSum.  See root.h for the
 * representation contract.  Both heads are held forms — the C
 * builtins below intentionally return NULL for every input so the
 * evaluator leaves the call exactly as the caller wrote it.  The
 * useful work happens elsewhere:
 *
 *   - src/deriv.c   — D[RootSum[f1, f2], x] threads through the body.
 *   - src/intrat.c  — NaiveLogPart constructs RootSum nodes when the
 *                     LRT log part has no closed-form real expression.
 *
 * Splitting this module out of intrat.c keeps the held-symbolic
 * machinery available to any future caller (e.g. Solve, Reduce,
 * factorisation over algebraic extensions) that needs to name the
 * roots of a polynomial without committing to a radical form.
 */

#include "root.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"

#include <stdlib.h>

/* Held builtin: the evaluator already left the call unevaluated by
 * the time we get here, so there is nothing to compute.  Returning
 * NULL preserves the input verbatim. */
static Expr* builtin_root(Expr* res) {
    (void)res;
    return NULL;
}

static Expr* builtin_rootsum(Expr* res) {
    (void)res;
    return NULL;
}

Expr* root_make_rootsum(Expr* bvar, Expr* poly, Expr* body) {
    /* Function[bvar, poly] and Function[bvar, body].  The bvar is
     * deep-copied for the second Function so each Function node owns
     * its own subtree. */
    Expr* fn1_args[2] = { bvar, poly };
    Expr* fn1 = expr_new_function(expr_new_symbol("Function"), fn1_args, 2);

    Expr* fn2_args[2] = { expr_copy(bvar), body };
    Expr* fn2 = expr_new_function(expr_new_symbol("Function"), fn2_args, 2);

    Expr* rs_args[2] = { fn1, fn2 };
    return expr_new_function(expr_new_symbol("RootSum"), rs_args, 2);
}

void root_init(void) {
    symtab_add_builtin("Root", builtin_root);
    symtab_get_def("Root")->attributes |= ATTR_PROTECTED | ATTR_HOLDALL;
    symtab_set_docstring("Root",
        "Root[Function[t, p[t]], k]\n"
        "\tRepresents the k-th real root of the univariate polynomial p\n"
        "\tin the variable t.  Held symbolic form — the system makes no\n"
        "\tattempt to express the root as a closed-form radical.");

    symtab_add_builtin("RootSum", builtin_rootsum);
    symtab_get_def("RootSum")->attributes |= ATTR_PROTECTED | ATTR_HOLDALL;
    symtab_set_docstring("RootSum",
        "RootSum[Function[t, p[t]], Function[t, body[t]]]\n"
        "\tThe formal sum of body[\\[Alpha]] over the roots \\[Alpha] of\n"
        "\tp[\\[Alpha]] == 0.  Held symbolic form, used by the rational\n"
        "\tintegrator's NaiveLogPart fallback when the logarithmic part\n"
        "\tcannot be expressed in closed-form real elementary functions.\n"
        "\tDifferentiation threads through the body Function:\n"
        "\t  D[RootSum[f1, Function[t, body]], x]\n"
        "\t    == RootSum[f1, Function[t, D[body, x]]].");
}
