/* polarplot.c — PolarPlot[r, {theta, tmin, tmax}, opts...]
 *
 * HoldAll: r and the iterator are unevaluated on entry, exactly like
 * ParametricPlot. The implementation converts the polar body r(theta) into
 * the Cartesian pair {r*Cos[theta], r*Sin[theta]} and delegates to
 * builtin_parametricplot, so all adaptive sampling, option handling,
 * multi-curve paletting, ColorFunction, Mesh, PlotLegends, etc. come for
 * free without duplicating any logic. */

#include "polarplot.h"
#include "parametricplot.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* Build {r_expr * Cos[theta], r_expr * Sin[theta]} (caller owns result). */
static Expr* make_polar_pair(Expr* r_expr, Expr* theta) {
    Expr* cos_a[1] = { expr_copy(theta) };
    Expr* cos_t    = expr_new_function(expr_new_symbol(SYM_Cos), cos_a, 1);
    Expr* sin_a[1] = { expr_copy(theta) };
    Expr* sin_t    = expr_new_function(expr_new_symbol(SYM_Sin), sin_a, 1);
    Expr* xa[2]    = { expr_copy(r_expr), cos_t };
    Expr* ya[2]    = { expr_copy(r_expr), sin_t };
    Expr* x        = expr_new_function(expr_new_symbol(SYM_Times), xa, 2);
    Expr* y        = expr_new_function(expr_new_symbol(SYM_Times), ya, 2);
    Expr* pair[2]  = { x, y };
    return expr_new_function(expr_new_symbol(SYM_List), pair, 2);
}

Expr* builtin_polarplot(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    Expr* r_body = res->data.function.args[0];
    Expr* iter   = res->data.function.args[1];

    /* Validate: {theta, tmin, tmax} */
    if (!iter || iter->type != EXPR_FUNCTION
        || !iter->data.function.head
        || iter->data.function.head->type != EXPR_SYMBOL
        || iter->data.function.head->data.symbol.name != SYM_List
        || iter->data.function.arg_count != 3
        || iter->data.function.args[0]->type != EXPR_SYMBOL)
        return NULL;

    Expr* theta = iter->data.function.args[0]; /* borrowed */

    /* Build the Cartesian body.
     * Single curve: r -> {r*Cos[th], r*Sin[th]}.
     * Multi-curve:  {r1, r2, ...} -> {{r1*Cos[th], r1*Sin[th]}, ...}.
     * Multi is detected when arg 0 is a non-empty List — ParametricPlot's
     * own multi-curve detection then fires because each element IS a List. */
    bool is_list = (r_body->type == EXPR_FUNCTION
                    && r_body->data.function.head
                    && r_body->data.function.head->type == EXPR_SYMBOL
                    && r_body->data.function.head->data.symbol.name == SYM_List
                    && r_body->data.function.arg_count > 0);

    Expr* param_body;
    if (is_list) {
        size_t nc    = r_body->data.function.arg_count;
        Expr** pairs = malloc(sizeof(Expr*) * nc);
        for (size_t i = 0; i < nc; i++)
            pairs[i] = make_polar_pair(r_body->data.function.args[i], theta);
        param_body = expr_new_function(expr_new_symbol(SYM_List), pairs, nc);
        free(pairs);
    } else {
        param_body = make_polar_pair(r_body, theta);
    }

    /* PolarPlot's default PlotPoints is 75 (vs. ParametricPlot's 25): polar
     * curves span a full 2π rotation so they need more initial seeds to seed
     * the adaptive refiner, otherwise petals and tight curvature look angular.
     * If the user already passed PlotPoints we honour that value unchanged. */
    bool have_plot_points = false;
    for (size_t i = 2; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (a && a->type == EXPR_FUNCTION && a->data.function.arg_count == 2
            && a->data.function.head
            && a->data.function.head->type == EXPR_SYMBOL
            && (a->data.function.head->data.symbol.name == SYM_Rule
                || a->data.function.head->data.symbol.name == SYM_RuleDelayed)
            && a->data.function.args[0]->type == EXPR_SYMBOL
            && a->data.function.args[0]->data.symbol.name == SYM_PlotPoints) {
            have_plot_points = true;
            break;
        }
    }

    /* Construct ParametricPlot[param_body, iter_copy, opts..., PlotPoints->75].
     * Every sub-expression is a fresh copy so expr_free(fake_call) after the
     * delegate returns is safe — the result tree is independent of fake_call. */
    size_t extra  = have_plot_points ? 0 : 1;
    Expr** new_args = malloc(sizeof(Expr*) * (argc + extra));
    new_args[0] = param_body;
    for (size_t i = 1; i < argc; i++)
        new_args[i] = expr_copy(res->data.function.args[i]);
    if (!have_plot_points) {
        Expr* pp_a[2] = { expr_new_symbol(SYM_PlotPoints), expr_new_integer(75) };
        new_args[argc] = expr_new_function(expr_new_symbol(SYM_Rule), pp_a, 2);
    }
    Expr* fake_call = expr_new_function(
        expr_new_symbol("ParametricPlot"), new_args, argc + extra);
    free(new_args);

    Expr* result = builtin_parametricplot(fake_call);
    expr_free(fake_call);
    return result;
}
