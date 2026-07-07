/* residue.c -- Residue[expr, {z, z0}], the symbolic residue.
 *
 * The residue of f at an isolated singularity z = z0 is the coefficient of
 * (z - z0)^-1 in the Laurent expansion of f. We obtain it directly from the
 * series engine: expand f to order (z - z0)^0 (which always spans the -1 term,
 * however deep the pole), then read the coefficient at exponent -1 out of the
 * resulting SeriesData[z, z0, {coefs}, nmin, nmax, den].
 *
 * A residue is well defined only for an ordinary Laurent expansion (den == 1).
 * A fractional-power (Puiseux) expansion, den > 1, signals a branch point,
 * where the residue is undefined -- we leave the call unevaluated, matching
 * Mathematica (e.g. Residue[1/Sqrt[z], {z, 0}]).
 */

#include "residue.h"
#include "series.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/* Emit the too-few-arguments diagnostic and leave the call unevaluated. */
static Expr* residue_emit_argcount(size_t argc) {
    fprintf(stderr,
            "Residue::argm: Residue called with %zu argument%s; "
            "2 or more arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/* Build and evaluate Series[expr, {z, z0, order}]. expr_new_function copies the
 * args array (memcpy) and adopts the element pointers, so the temporary arrays
 * are freed while their contents live on in the new nodes. Returns the
 * (owned) result of evaluation, which the caller must free. */
static Expr* residue_series(Expr* expr, Expr* z, Expr* z0, int64_t order) {
    Expr** spec_args = calloc(3, sizeof(Expr*));
    if (!spec_args) return NULL;
    spec_args[0] = expr_copy(z);
    spec_args[1] = expr_copy(z0);
    spec_args[2] = expr_new_integer(order);
    Expr* series_spec = expr_new_function(expr_new_symbol(SYM_List), spec_args, 3);
    free(spec_args);

    Expr** call_args = calloc(2, sizeof(Expr*));
    if (!call_args) { expr_free(series_spec); return NULL; }
    call_args[0] = expr_copy(expr);
    call_args[1] = series_spec;
    Expr* series_call = expr_new_function(expr_new_symbol("Series"), call_args, 2);
    free(call_args);

    return eval_and_free(series_call);
}

Expr* builtin_residue(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return residue_emit_argcount(argc);
    if (argc != 2) return NULL;   /* only the two-argument form is handled */

    Expr* expr = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    /* The location spec must be List[z, z0] with z a symbol. */
    if (spec->type != EXPR_FUNCTION ||
        spec->data.function.head->type != EXPR_SYMBOL ||
        spec->data.function.head->data.symbol != SYM_List ||
        spec->data.function.arg_count != 2)
        return NULL;
    Expr* z  = spec->data.function.args[0];
    Expr* z0 = spec->data.function.args[1];
    if (z->type != EXPR_SYMBOL) return NULL;

    /* Expand to order (z-z0)^0 and read the coefficient of (z-z0)^-1. Order 0
     * spans exponent -1 for most integrands, but when the leading behaviour is
     * carried by an unknown function (e.g. f[z]/z^5) the series engine truncates
     * relative to that function's depth and the O-term can land at an exponent
     * <= -1, leaving -1 unknown. In that case we raise the requested order until
     * the -1 coefficient is inside the explicit terms (bounded, to stay safe on
     * essential singularities where it never resolves). */
    int64_t order = 0;
    for (int attempt = 0; attempt < 8; attempt++) {
        Expr* sd = residue_series(expr, z, z0, order);
        if (!sd) return NULL;

        /* Series must have produced a SeriesData; otherwise it could not expand. */
        if (!is_series_data(sd)) { expr_free(sd); return NULL; }

        Expr** a     = sd->data.function.args;
        Expr* coefs  = a[2];
        Expr* nmin_e = a[3];
        Expr* nmax_e = a[4];
        Expr* den_e  = a[5];
        if (nmin_e->type != EXPR_INTEGER || nmax_e->type != EXPR_INTEGER ||
            den_e->type != EXPR_INTEGER || coefs->type != EXPR_FUNCTION) {
            expr_free(sd);
            return NULL;
        }
        int64_t nmin = nmin_e->data.integer;
        int64_t nmax = nmax_e->data.integer;
        int64_t den  = den_e->data.integer;

        /* Fractional exponents -> branch point -> residue undefined. */
        if (den != 1) { expr_free(sd); return NULL; }

        /* Coefficient of exponent -1 lives at index (-1 - nmin) in the list;
         * it is known when nmin <= -1 <= nmax-1 (i.e. 0 <= index < len). */
        int64_t len   = (int64_t)coefs->data.function.arg_count;
        int64_t index = -1 - nmin;

        if (index < 0) {
            /* nmin > -1: f is analytic at z0 (no principal part) -> residue 0. */
            expr_free(sd);
            return expr_new_integer(0);
        }
        if (index < len) {
            Expr* result = expr_copy(coefs->data.function.args[index]);
            expr_free(sd);
            return result;
        }

        /* -1 lies at or beyond the O-term (nmax <= -1). Raise the order enough to
         * push the O-term past exponent -1, with a small margin, then retry. */
        int64_t next = order + (0 - nmax) + 2;
        if (next <= order) next = order + 1;   /* guarantee progress */
        expr_free(sd);
        if (next > 256) return NULL;           /* safety cap */
        order = next;
    }
    return NULL;
}

void residue_init(void) {
    symtab_add_builtin("Residue", builtin_residue);
    symtab_get_def("Residue")->attributes |= ATTR_PROTECTED;
}
