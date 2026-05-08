/* integrate.c
 *
 * Tiny `Integrate[f, x]` System` dispatcher.  Recognises polynomial /
 * rational integrands in `x` and forwards to the package entry
 * `Integrate`IntegrateRational[f, x]`.  Anything else returns NULL so
 * the evaluator leaves the call unevaluated (matching Mathematica's
 * unsimplified-form behaviour for non-rational integrands).
 *
 * All algorithmic content lives in intrat.c; this file is the System
 * surface area only.
 */

#include "integrate.h"
#include "intrat.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "internal.h"
#include "rationalize.h"
#include "sym_names.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* True iff `e` is the symbol `True`.  The PolynomialQ / rationalQ
 * predicates we call below return either True or False. */
static bool is_true_symbol(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == SYM_True;
}

/* Recursively scan `e` for any EXPR_REAL leaf.  picocas's rational
 * integrator cannot reason about inexact numbers; we surface them at
 * dispatch with an Integrate::inexact message rather than letting the
 * pipeline produce ComplexInfinity / 1/0 partway through. */
static bool expr_contains_real(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
    if (e->type == EXPR_FUNCTION) {
        if (expr_contains_real(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (expr_contains_real(e->data.function.args[i])) return true;
        }
    }
    return false;
}

/* Test whether `f` is a polynomial in `x`.  Calls the existing
 * PolynomialQ builtin so we get the same definition the rest of
 * picocas uses. */
static bool is_polynomial_in(Expr* f, Expr* x) {
    Expr* args[2] = { expr_copy(f), expr_copy(x) };
    Expr* call = internal_polynomialq(args, 2);
    Expr* val  = evaluate(call);
    expr_free(call);
    bool ok = is_true_symbol(val);
    expr_free(val);
    return ok;
}

/* Test whether `f` is a rational function in `x`: Together[f] must
 * have a non-trivial denominator that is itself polynomial in `x`,
 * with a polynomial numerator. */
static bool is_rational_in(Expr* f, Expr* x) {
    Expr* together = internal_together((Expr*[]){expr_copy(f)}, 1);
    Expr* combined = evaluate(together);
    expr_free(together);
    if (!combined) return false;

    Expr* num = internal_numerator((Expr*[]){expr_copy(combined)}, 1);
    Expr* den = internal_denominator((Expr*[]){expr_copy(combined)}, 1);
    expr_free(combined);
    Expr* num_v = evaluate(num);
    Expr* den_v = evaluate(den);
    expr_free(num); expr_free(den);

    bool ok = false;
    if (num_v && den_v) {
        bool num_is_poly = is_polynomial_in(num_v, x);
        bool den_is_poly = is_polynomial_in(den_v, x);
        /* Reject the trivial Denominator==1 case here: a pure polynomial
         * numerator already passed the PolynomialQ branch upstream. */
        bool den_is_unit = (den_v->type == EXPR_INTEGER && den_v->data.integer == 1);
        ok = num_is_poly && den_is_poly && !den_is_unit;
    }
    if (num_v) expr_free(num_v);
    if (den_v) expr_free(den_v);
    return ok;
}

Expr* builtin_integrate(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;

    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];

    /* The integration variable must be a (single) symbol so we can
     * compute derivatives, partial fractions, ... in it. */
    if (x->type != EXPR_SYMBOL) return NULL;

    /* Inexact integrands: try Rationalize[f] first.  Floats with an
     * exact rational equivalent within the default tolerance (4.5,
     * 3.125, ...) coerce cleanly; transcendental floats (N[Pi], ...)
     * survive Rationalize unchanged, in which case we emit
     * Integrate::inexact and bubble back unevaluated.  Hash-dedupe
     * so the eval-loop's repeated re-entry doesn't double-print. */
    Expr* coerced = NULL;
    Expr* effective_f = f;
    if (expr_contains_real(f)) {
        coerced = internal_rationalize_expr(f, 0.0, RATIONALIZE_DEFAULT);
        if (!coerced || expr_contains_real(coerced)) {
            static uint64_t last_warned_hash = 0;
            uint64_t h = expr_hash(res);
            if (h != last_warned_hash) {
                fprintf(stderr,
                    "Integrate::inexact: Integrand contains inexact numbers which "
                    "cannot be coerced to exact numbers with Rationalize[number].\n");
                last_warned_hash = h;
            }
            if (coerced) expr_free(coerced);
            return NULL;
        }
        effective_f = coerced;
    }

    bool dispatch = is_polynomial_in(effective_f, x) || is_rational_in(effective_f, x);
    if (!dispatch) {
        if (coerced) expr_free(coerced);
        return NULL;
    }

    /* Forward to the Integrate` package entry point. */
    Expr* call = expr_new_function(
        expr_new_symbol("Integrate`IntegrateRational"),
        (Expr*[]){ expr_copy(effective_f), expr_copy(x) }, 2);
    if (coerced) expr_free(coerced);
    Expr* result = evaluate(call);
    expr_free(call);
    if (!result) return NULL;

    /* Detect the unevaluated package call and bubble back as
     * Integrate[f, x] so the user-facing form stays clean.  A
     * resolved result has any other head (Plus, Log, ...). */
    if (result->type == EXPR_FUNCTION
        && result->data.function.head
        && result->data.function.head->type == EXPR_SYMBOL
        && strcmp(result->data.function.head->data.symbol,
                  "Integrate`IntegrateRational") == 0) {
        expr_free(result);
        return NULL;
    }
    return result;
}

void integrate_init(void) {
    symtab_add_builtin("Integrate", builtin_integrate);
    symtab_get_def("Integrate")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Integrate",
        "Integrate[f, x] gives the indefinite integral of f with respect to x.\n"
        "Currently dispatches to the rational-function pipeline implemented in\n"
        "the Integrate` package: polynomial term-by-term integration, Hermite\n"
        "reduction for repeated roots, derivative-recognition fast path for\n"
        "c*D'/D and c*D'/D^k, and (in later phases) the Lazard-Rioboo-Trager\n"
        "logarithmic part. Non-rational integrands return unevaluated.");

    /* Initialise the Integrate` package: HermiteReduce, IntegratePolynomial,
     * helpers, and the explicit `Integrate`IntegrateRational` entry. */
    intrat_init();
}
