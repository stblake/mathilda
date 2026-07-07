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
 *
 * Algebraic pole locations.  The series engine decides whether z0 is a pole by
 * evaluating the denominator there and testing it against zero; but for a pole
 * whose location is a SUM of radicals (e.g. z0 = -2 + Sqrt[3], a root of
 * 1 + 4 z + z^2), Denominator(z0) is an expression like
 * 1 + 4 (-2 + Sqrt[3]) + (-2 + Sqrt[3])^2 that does not auto-simplify to 0, so
 * the pole is missed and the residue wrongly comes out 0.  We defeat this by
 * expanding about z0 EXPLICITLY: substitute z -> z0 + w, then Expand the
 * denominator of the result -- polynomial expansion collapses the radical
 * arithmetic (Sqrt[3]^2 -> 3, ...) so the vanishing constant term becomes a
 * literal 0 and the w-factor of the pole is exposed.  Reading the (z-z0)^-1
 * coefficient is then a plain Series-at-0 of the expanded form.
 */

#include "residue.h"
#include "series.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "internal.h"
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

/* Build head[a] and evaluate it, freeing the call.  `a` is consumed. */
static Expr* residue_eval1(const char* head, Expr* a) {
    Expr** args = calloc(1, sizeof(Expr*));
    if (!args) { expr_free(a); return NULL; }
    args[0] = a;
    Expr* call = expr_new_function(expr_new_symbol(head), args, 1);
    free(args);
    return eval_and_free(call);
}

/* Rewrite `f` so that the singularity at z = z0 is expanded about w = 0 with an
 * EXPANDED denominator: returns Numerator[Together[f/.z->z0+w]] /
 * Expand[Denominator[Together[f/.z->z0+w]]] in the fresh variable `w`.  This is
 * what makes an algebraic pole location (z0 a sum of radicals) detectable: the
 * expansion collapses the radical arithmetic in the denominator's constant term
 * to a literal 0.  Returns an owned Expr* (in the variable `w`) or NULL.  `f`,
 * `z`, `z0`, `w` are borrowed. */
static Expr* residue_shift_form(Expr* f, Expr* z, Expr* z0, Expr* w) {
    /* g = f /. z -> z0 + w */
    Expr* shift = expr_new_function(expr_new_symbol(SYM_Plus),
                      (Expr*[]){ expr_copy(z0), expr_copy(w) }, 2);
    Expr* rule  = expr_new_function(expr_new_symbol(SYM_Rule),
                      (Expr*[]){ expr_copy(z), shift }, 2);
    Expr* repl  = expr_new_function(expr_new_symbol("ReplaceAll"),
                      (Expr*[]){ expr_copy(f), rule }, 2);
    Expr* g = eval_and_free(repl);
    if (!g) return NULL;

    Expr* tog = residue_eval1("Together", g);          /* consumes g */
    if (!tog) return NULL;
    Expr* num = residue_eval1("Numerator", expr_copy(tog));
    Expr* den = residue_eval1("Denominator", tog);     /* consumes tog */
    if (!num || !den) { if (num) expr_free(num); if (den) expr_free(den); return NULL; }
    Expr* denx = residue_eval1("Expand", den);         /* consumes den */
    if (!denx) { expr_free(num); return NULL; }

    /* form = num * denx^-1 */
    Expr* denpow = expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ denx, expr_new_integer(-1) }, 2);
    Expr* form   = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ num, denpow }, 2);
    return form;
}

/* True iff PolynomialQ[p, z]. */
static bool residue_polyq(Expr* p, Expr* z) {
    Expr* call = internal_polynomialq(
        (Expr*[]){ expr_copy(p), expr_copy(z) }, 2);
    Expr* v = eval_and_free(call);
    bool ok = v && v->type == EXPR_SYMBOL && v->data.symbol == SYM_True;
    if (v) expr_free(v);
    return ok;
}

/* True iff f is a rational function of z: Together[f] has polynomial numerator
 * AND polynomial denominator in z.  The shift+Expand pole-detection preprocessing
 * is applied ONLY in this case; for transcendental / special-function integrands
 * (Cot, Zeta near its pole, unknown f[z], ...) the series engine's built-in
 * knowledge of the expansion at z0 is preferable, and shifting the argument would
 * defeat it. */
static bool residue_is_rational_in(Expr* f, Expr* z) {
    Expr* tog = residue_eval1("Together", expr_copy(f));
    if (!tog) return false;
    Expr* num = residue_eval1("Numerator", expr_copy(tog));
    Expr* den = residue_eval1("Denominator", tog);   /* consumes tog */
    if (!num || !den) { if (num) expr_free(num); if (den) expr_free(den); return false; }
    bool ok = residue_polyq(num, z) && residue_polyq(den, z);
    expr_free(num);
    expr_free(den);
    return ok;
}

/* Adaptive Laurent-coefficient extraction: expand `expr` about `spt` in `svar`
 * to order 0 (raising it when an unknown-function leading term hides exponent
 * -1) and return an owned copy of the (svar-spt)^-1 coefficient, 0 at an analytic
 * point, or NULL at a branch point / when no series can be produced. */
static Expr* residue_extract(Expr* expr, Expr* svar, Expr* spt) {
    int64_t order = 0;
    for (int attempt = 0; attempt < 8; attempt++) {
        Expr* sd = residue_series(expr, svar, spt, order);
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
            /* nmin > -1: analytic (no principal part) -> residue 0. */
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

Expr* residue_compute(Expr* f, Expr* z, Expr* z0) {
    if (!f || !z || !z0 || z->type != EXPR_SYMBOL) return NULL;

    if (residue_is_rational_in(f, z)) {
        /* Rational integrand: expand about z0 with an EXPANDED denominator so an
         * algebraic pole location (z0 a sum of radicals) is exposed as a w-factor. */
        Expr* w = expr_new_symbol("Residue`$w");
        Expr* form = residue_shift_form(f, z, z0, w);
        if (!form) { expr_free(w); return NULL; }
        Expr* zero = expr_new_integer(0);
        Expr* result = residue_extract(form, w, zero);
        expr_free(form);
        expr_free(w);
        expr_free(zero);
        return result;
    }

    /* Transcendental / special-function integrand: expand directly about z0 so
     * the series engine can use its knowledge of the function's Laurent series
     * there (e.g. Zeta at 1, Cot / 1/Sin^n at 0, unknown f[z]/z^n). */
    return residue_extract(f, z, z0);
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

    return residue_compute(expr, z, z0);
}

void residue_init(void) {
    symtab_add_builtin("Residue", builtin_residue);
    symtab_get_def("Residue")->attributes |= ATTR_PROTECTED;
}
