/*
 * integrate_interp.c
 *
 * Indefinite integration of an applied InterpolatingFunction object.
 *
 *     Integrate[InterpolatingFunction[...][x], x]
 *
 * Differentiation (D) only needs to bump the derivative-order annotation the
 * object already carries: D[ifun[x], x] reduces to a derivative-annotated
 * InterpolatingFunction whose per-window kernels then evaluate the next
 * derivative.  Integration has no such shortcut — the kernels evaluate the
 * d-th derivative for d >= 0 only and cannot produce an antiderivative — so we
 * build a genuinely new interpolant of the antiderivative:
 *
 *   1. Read the grid x-coordinates from the object's stored table.
 *   2. Sample the original function values  y_i = ifun[x_i].
 *   3. Accumulate the antiderivative node values  F_0 = 0,
 *        F_i = F_{i-1} + Integrate_{x_{i-1}}^{x_i} ifun
 *      using Gauss-Legendre quadrature of the (already piecewise-polynomial)
 *      interpolant.
 *   4. Build a Hermite InterpolatingFunction through {{x_i}, F_i, y_i}: the
 *      antiderivative's exact derivative is the original function, so supplying
 *      F'(x_i) = y_i makes D[Integrate[ifun[x], x], x] round-trip to ifun[x].
 *
 * This uniform construction also handles derivative-annotated objects (the ones
 * D produces): interp_apply evaluates whatever the object represents, and
 * integrating those samples recovers the lower-order antiderivative.
 *
 * Only the 1-D, direct case (the applied argument is the integration variable
 * itself) is reduced; Integrate[ifun[g[x]], x] is not generally expressible as
 * an InterpolatingFunction, so it is left to the caller.
 *
 * Computations use machine doubles (the dominant use of InterpolatingFunction);
 * high-precision objects are integrated at double precision.
 */
#include "integrate_interp.h"
#include "../interp.h"
#include "../eval.h"
#include "../sym_names.h"

#include <stdlib.h>
#include <gmp.h>

/* ---------------------------------------------------------------------- */
/* Small local helpers                                                    */
/* ---------------------------------------------------------------------- */

/* head symbol == name (symbols are interned, so pointer compare). */
static bool head_is(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

static bool is_list(const Expr* e) { return head_is(e, SYM_List); }

/* Coerce a real-valued atom to a double (Integer/Real/BigInt/MPFR/Rational). */
static bool real_to_double(const Expr* e, double* out) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer;       return true;
        case EXPR_REAL:    *out = e->data.real;                  return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint);     return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        default: break;
    }
    if (head_is(e, SYM_Rational) && e->data.function.arg_count == 2
        && expr_is_integer_like(e->data.function.args[0])
        && expr_is_integer_like(e->data.function.args[1])) {
        mpz_t num, den;
        expr_to_mpz(e->data.function.args[0], num);
        expr_to_mpz(e->data.function.args[1], den);
        double dd = mpz_get_d(den);
        bool ok = (dd != 0.0);
        if (ok) *out = mpz_get_d(num) / dd;
        mpz_clears(num, den, NULL);
        return ok;
    }
    return false;
}

/* The scalar x-coordinate Expr of a stored table entry, or NULL if the entry
 * is malformed / not 1-D.  Value entries are {x, y, ...} (x a scalar);
 * derivative-supplied entries are {{x}, val, der, ...} (x in a one-element
 * coordinate list).  The returned pointer is borrowed from `entry`. */
static Expr* entry_coord(Expr* entry) {
    if (!is_list(entry) || entry->data.function.arg_count < 2) return NULL;
    Expr* c = entry->data.function.args[0];
    if (is_list(c)) {
        if (c->data.function.arg_count != 1) return NULL;    /* not 1-D */
        c = c->data.function.args[0];
    }
    return c;
}

/* Evaluate ifun at the single real point `p`.  Returns false if the object
 * does not yield a real value there. */
static bool eval_ifun(Expr* obj, double p, double* out) {
    Expr* arg = expr_new_real(p);
    Expr* val = interp_apply(obj, &arg, 1);
    expr_free(arg);
    if (!val) return false;
    bool ok = real_to_double(val, out);
    expr_free(val);
    return ok;
}

/* ---------------------------------------------------------------------- */
/* Gauss-Legendre quadrature (5-point, exact through degree 9)            */
/* ---------------------------------------------------------------------- */
/* The default interpolant is piecewise-polynomial of order <= 3 and the
 * Spline/Hermite forms are cubic, so 5 nodes integrate each panel exactly.
 * Explicit high-order interpolants (InterpolationOrder -> k, k > 9) incur a
 * small quadrature error in the accumulated node values. */
static const double GL_NODE[5] = {
    -0.906179845938664, -0.538469310105683, 0.0,
     0.538469310105683,  0.906179845938664
};
static const double GL_WEIGHT[5] = {
    0.236926885056189, 0.478628670499366, 0.568888888888889,
    0.478628670499366, 0.236926885056189
};

/* Integral of ifun over [a, b], or false if a sample point cannot be
 * evaluated. */
static bool panel_integral(Expr* obj, double a, double b, double* out) {
    double half = 0.5 * (b - a), mid = 0.5 * (a + b), acc = 0.0;
    for (int k = 0; k < 5; k++) {
        double v;
        if (!eval_ifun(obj, mid + half * GL_NODE[k], &v)) return false;
        acc += GL_WEIGHT[k] * v;
    }
    *out = half * acc;
    return true;
}

/* ---------------------------------------------------------------------- */
/* Entry point                                                            */
/* ---------------------------------------------------------------------- */

Expr* integrate_interp(Expr* f, Expr* x) {
    if (!f || !x || x->type != EXPR_SYMBOL) return NULL;

    /* f must be obj[arg] with obj an InterpolatingFunction and arg == x. */
    if (f->type != EXPR_FUNCTION || f->data.function.arg_count != 1) return NULL;
    Expr* obj = f->data.function.head;
    if (!head_is(obj, SYM_InterpolatingFunction)) return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;

    if (obj->data.function.arg_count < 2) return NULL;
    Expr* domain = obj->data.function.args[0];
    Expr* table  = obj->data.function.args[1];

    /* 1-D objects only: domain = {{a, b}}. */
    if (!is_list(domain) || domain->data.function.arg_count != 1) return NULL;
    if (!is_list(table)) return NULL;
    size_t n = table->data.function.arg_count;
    if (n < 2) return NULL;

    double* xs = malloc(n * sizeof(double));
    double* ys = malloc(n * sizeof(double));
    double* Fs = malloc(n * sizeof(double));
    Expr** xe = malloc(n * sizeof(Expr*));   /* borrowed exact coordinates */
    if (!xs || !ys || !Fs || !xe) {
        free(xs); free(ys); free(Fs); free(xe); return NULL;
    }

    bool ok = true;
    for (size_t i = 0; i < n && ok; i++) {
        xe[i] = entry_coord(table->data.function.args[i]);
        if (!xe[i] || !real_to_double(xe[i], &xs[i])) ok = false;
        else if (i > 0 && xs[i] <= xs[i - 1]) ok = false;   /* must increase */
        else if (!eval_ifun(obj, xs[i], &ys[i])) ok = false;
    }

    /* Accumulate antiderivative node values, F(x_0) = 0. */
    if (ok) {
        Fs[0] = 0.0;
        for (size_t i = 1; i < n && ok; i++) {
            double seg;
            if (!panel_integral(obj, xs[i - 1], xs[i], &seg)) ok = false;
            else Fs[i] = Fs[i - 1] + seg;
        }
    }

    /* Build Hermite data {{x_i}, F_i, y_i} and let Interpolation construct the
     * antiderivative object (value + exact first derivative at each node).  The
     * coordinate is the original (exact) node Expr, so an integer grid keeps an
     * integer domain. */
    Expr** entries = ok ? malloc(n * sizeof(Expr*)) : NULL;
    if (entries) {
        for (size_t i = 0; i < n; i++) {
            Expr* xc = expr_copy(xe[i]);
            Expr* coord = expr_new_function(expr_new_symbol(SYM_List), &xc, 1);
            Expr* trip[3] = { coord, expr_new_real(Fs[i]), expr_new_real(ys[i]) };
            entries[i] = expr_new_function(expr_new_symbol(SYM_List), trip, 3);
        }
    }
    free(xs); free(ys); free(Fs); free(xe);
    if (!entries) return NULL;

    Expr* data = expr_new_function(expr_new_symbol(SYM_List), entries, n);
    free(entries);
    Expr* interp_call = expr_new_function(expr_new_symbol(SYM_Interpolation),
                                          &data, 1);
    Expr* antideriv = evaluate(interp_call);   /* evaluate() does not free arg */
    expr_free(interp_call);

    if (!head_is(antideriv, SYM_InterpolatingFunction)) {
        expr_free(antideriv);
        return NULL;
    }

    /* Return the antiderivative applied to x, mirroring D's ifun'[x] result. */
    Expr* xc = expr_copy(x);
    return expr_new_function(antideriv, &xc, 1);
}
