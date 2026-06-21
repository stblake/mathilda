#include "list_common.h"
#include "rescale.h"

/* ====================================================================
 * Rescale[x, {min, max}]          -> (x - min)/(max - min)
 * Rescale[x, {min, max}, {y0, y1}] -> y0 + (y1 - y0)(x - min)/(max - min)
 * Rescale[list]                    -> Rescale[list, {Min[list], Max[list]}]
 *
 * Rescale threads element-wise over a List first argument and works
 * uniformly on integers, rationals, reals, arbitrary-precision numbers,
 * complex numbers, and purely symbolic quantities, because it defers all
 * arithmetic to the core evaluator rather than computing values itself.
 * ==================================================================== */


/* a - b, taking ownership of (adopting) the freshly-built a and b. */
static Expr* rescale_sub(Expr* a, Expr* b) {
    Expr* neg_b = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(-1), b }, 2);
    return expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ a, neg_b }, 2);
}

/* a / b, taking ownership of (adopting) the freshly-built a and b. */
static Expr* rescale_div(Expr* a, Expr* b) {
    Expr* inv_b = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ b, expr_new_integer(-1) }, 2);
    return expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ a, inv_b }, 2);
}

/* `Rescale::argb: Rescale called with N argument(s); between 1 and 3
 * arguments are expected.` */
static Expr* rescale_emit_argb(size_t argc) {
    fprintf(stderr,
            "Rescale::argb: Rescale called with %zu argument%s; "
            "between 1 and 3 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_rescale(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return rescale_emit_argb(argc);

    Expr* x = res->data.function.args[0];

    /* One-argument form: reduce to the two-argument form with the data
     * range {Min[x], Max[x]} and re-evaluate.  For a List this threads;
     * for a lone scalar Min[x] == Max[x] and the result is Indeterminate,
     * matching Mathematica. */
    if (argc == 1) {
        Expr* min_x = expr_new_function(expr_new_symbol(SYM_Min),
            (Expr*[]){ expr_copy(x) }, 1);
        Expr* max_x = expr_new_function(expr_new_symbol(SYM_Max),
            (Expr*[]){ expr_copy(x) }, 1);
        Expr* range = expr_new_function(expr_new_symbol(SYM_List),
            (Expr*[]){ min_x, max_x }, 2);
        Expr* call = expr_new_function(expr_new_symbol(SYM_Rescale),
            (Expr*[]){ expr_copy(x), range }, 2);
        Expr* out = evaluate(call);
        expr_free(call);
        return out;
    }

    /* The data range must be a two-element List {min, max}. */
    Expr* range = res->data.function.args[1];
    if (!is_listq(range) || range->data.function.arg_count != 2) return NULL;
    Expr* min_v = range->data.function.args[0];
    Expr* max_v = range->data.function.args[1];

    /* The optional target range must also be a two-element List. */
    Expr* yrange = NULL;
    if (argc == 3) {
        yrange = res->data.function.args[2];
        if (!is_listq(yrange) || yrange->data.function.arg_count != 2) return NULL;
    }

    /* Thread over a List first argument: Rescale[{e1, ...}, range, ...]
     * -> {Rescale[e1, range, ...], ...}. */
    if (is_listq(x)) {
        size_t n = x->data.function.arg_count;
        Expr** out_args = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            Expr* el = x->data.function.args[i];
            Expr* call;
            if (argc == 2) {
                call = expr_new_function(expr_new_symbol(SYM_Rescale),
                    (Expr*[]){ expr_copy(el), expr_copy(range) }, 2);
            } else {
                call = expr_new_function(expr_new_symbol(SYM_Rescale),
                    (Expr*[]){ expr_copy(el), expr_copy(range), expr_copy(yrange) }, 3);
            }
            out_args[i] = evaluate(call);
            expr_free(call);
        }
        Expr* out = expr_new_function(expr_new_symbol(SYM_List), out_args, n);
        free(out_args);
        return out;
    }

    /* Scalar: (x - min)/(max - min). */
    Expr* frac = rescale_div(rescale_sub(expr_copy(x), expr_copy(min_v)),
                             rescale_sub(expr_copy(max_v), expr_copy(min_v)));

    if (argc == 2) {
        Expr* out = evaluate(frac);
        expr_free(frac);
        return out;
    }

    /* Three-argument: y0 + (y1 - y0)(x - min)/(max - min). */
    Expr* ymin = yrange->data.function.args[0];
    Expr* ymax = yrange->data.function.args[1];
    Expr* scaled = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ rescale_sub(expr_copy(ymax), expr_copy(ymin)), frac }, 2);
    Expr* full = expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy(ymin), scaled }, 2);
    Expr* out = evaluate(full);
    expr_free(full);
    return out;
}
