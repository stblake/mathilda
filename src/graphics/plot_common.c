/* plot_common.c — see plot_common.h. Extracted from plot.c so Plot3D
 * (plot3d.c) reuses the exact same option/evaluation idioms as Plot. */

#include "plot_common.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif
#include <string.h>
#include <math.h>

bool expr_to_real_double(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real; return true; }
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint); return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)   { *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true; }
#endif
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational) {
        Expr* n = e->data.function.args[0];
        Expr* d = e->data.function.args[1];
        if (n->type == EXPR_INTEGER && d->type == EXPR_INTEGER && d->data.integer != 0) {
            *out = (double)n->data.integer / (double)d->data.integer;
            return true;
        }
    }
    return false;
}

bool numericize_bound(Expr* e, double* out) {
    Expr* n_arg[1] = { expr_copy(e) };
    Expr* n_call = expr_new_function(expr_new_symbol("N"), n_arg, 1);
    Expr* result = evaluate(n_call);
    expr_free(n_call);
    bool ok = expr_to_real_double(result, out) && isfinite(*out);
    expr_free(result);
    return ok;
}

bool is_rule_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    Expr* h = e->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return false;
    return (h->data.symbol == SYM_Rule || h->data.symbol == SYM_RuleDelayed)
        && e->data.function.arg_count == 2;
}

bool parse_long_value(Expr* rhs, long* out) {
    Expr* v = evaluate(expr_copy(rhs));
    bool ok = (v->type == EXPR_INTEGER);
    if (ok) *out = (long)v->data.integer;
    expr_free(v);
    return ok;
}

Expr* palette_color(size_t i) {
    static const double pal[][3] = {
        { 0.368417, 0.506779, 0.709798 },
        { 0.880722, 0.611041, 0.142051 },
        { 0.560181, 0.691569, 0.194885 },
        { 0.922526, 0.385626, 0.209179 },
        { 0.528488, 0.470624, 0.701351 },
        { 0.772079, 0.431554, 0.102387 },
        { 0.363898, 0.618501, 0.782349 },
        { 1.000000, 0.750000, 0.000000 },
        { 0.647624, 0.378160, 0.614037 },
        { 0.571589, 0.586483, 0.000000 },
    };
    size_t k = i % (sizeof(pal) / sizeof(pal[0]));
    Expr* a[3] = { expr_new_real(pal[k][0]), expr_new_real(pal[k][1]), expr_new_real(pal[k][2]) };
    return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
}

bool eval_region(Expr* region_fn, double x, double y) {
    Expr* args2[2] = { expr_new_real(x), expr_new_real(y) };
    Expr* call2 = expr_new_function(expr_copy(region_fn), args2, 2);
    Expr* r2 = evaluate(call2);
    expr_free(call2); /* evaluate() borrows its argument; the call node is ours to free */
    bool true2 = (r2->type == EXPR_SYMBOL && r2->data.symbol == SYM_True);
    bool false2 = (r2->type == EXPR_SYMBOL && r2->data.symbol == SYM_False);
    expr_free(r2);
    if (true2) return true;
    if (false2) return false;

    /* The 2-arg call didn't resolve to a boolean (likely a 1-arg function,
     * e.g. Function[x, x > 0]) -- retry with just x. */
    Expr* args1[1] = { expr_new_real(x) };
    Expr* call1 = expr_new_function(expr_copy(region_fn), args1, 1);
    Expr* r1 = evaluate(call1);
    expr_free(call1);
    bool ok = (r1->type == EXPR_SYMBOL && r1->data.symbol == SYM_True);
    expr_free(r1);
    return ok;
}

static bool is_color_head(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && (e->data.function.head->data.symbol == SYM_RGBColor
            || e->data.function.head->data.symbol == SYM_GrayLevel
            || e->data.function.head->data.symbol == SYM_Hue);
}

Expr* eval_color_function(Expr* color_fn, double x, double y,
                           double xmin, double xmax, bool scaling) {
    if (color_fn->type == EXPR_STRING && strcmp(color_fn->data.string, "Rainbow") == 0) {
        double t = (xmax > xmin) ? (x - xmin) / (xmax - xmin) : 0.0;
        Expr* a[1] = { expr_new_real(t * 0.8) }; /* stop short of wrapping back to red */
        return expr_new_function(expr_new_symbol(SYM_Hue), a, 1);
    }

    double cx = (scaling && xmax > xmin) ? (x - xmin) / (xmax - xmin) : x;

    Expr* args2[2] = { expr_new_real(cx), expr_new_real(y) };
    Expr* call2 = expr_new_function(expr_copy(color_fn), args2, 2);
    Expr* r2 = evaluate(call2);
    expr_free(call2);
    if (is_color_head(r2)) return r2;
    expr_free(r2);

    Expr* args1[1] = { expr_new_real(cx) };
    Expr* call1 = expr_new_function(expr_copy(color_fn), args1, 1);
    Expr* r1 = evaluate(call1);
    expr_free(call1);
    if (is_color_head(r1)) return r1;
    expr_free(r1);

    Expr* a[1] = { expr_new_real(0.5) };
    return expr_new_function(expr_new_symbol(SYM_GrayLevel), a, 1);
}
