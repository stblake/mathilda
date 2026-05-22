/* Normalize[v]      — v / Norm[v], leaving the zero vector unchanged.
 * Normalize[z]      — for a scalar (incl. complex) z, z / Abs[z].
 *                     Norm[z] reduces to Abs[z], so the same code path
 *                     handles scalars.
 * Normalize[expr,f] — expr / f[expr], same zero short-circuit.
 *
 * We delegate the heavy lifting to the evaluator: build
 *
 *     Times[expr, Power[f[expr], -1]]
 *
 * and evaluate it.  Times is Listable, so a single Times call distributes
 * the scalar reciprocal across every leaf of a vector / matrix / higher-
 * rank tensor input.  When the norm evaluates to an exact numeric zero
 * (Integer 0, Real 0.0, BigInt 0, or MPFR 0) we skip the division and
 * return the input untouched, matching Mathematica's documented behaviour
 * that "zero vectors are returned unchanged."
 *
 * Arity diagnostics:
 *   Normalize[]            -> `Normalize::argt`, call left unevaluated.
 *   Normalize[a, b, c, …]  -> same.
 */

#include "linalg.h"
#include "eval.h"
#include "expr.h"

#include <gmp.h>
#ifdef USE_MPFR
#  include <mpfr.h>
#endif

#include <stddef.h>
#include <stdio.h>

/* True iff `e` is an exact numeric zero.  We use *exact* equality (not a
 * symbolic-zero predicate) deliberately: short-circuiting on a symbolic
 * "Sqrt[Abs[x]^2 + Abs[y]^2]" that just *happens* to be zero would
 * silently drop the answer to a nonsense input.  Limiting to literal 0
 * keeps the contract honest -- if the norm cannot be proved zero, we
 * keep the symbolic division so the user can see what they wrote.
 */
static bool norm_is_numeric_zero(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: return e->data.integer == 0;
        case EXPR_REAL:    return e->data.real == 0.0;
        case EXPR_BIGINT:  return mpz_sgn(e->data.bigint) == 0;
#ifdef USE_MPFR
        case EXPR_MPFR:    return mpfr_zero_p(e->data.mpfr) != 0;
#endif
        default:           return false;
    }
}

Expr* builtin_normalize(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) {
        fprintf(stderr,
                "Normalize::argt: Normalize called with %zu argument%s; "
                "1 or 2 arguments are expected.\n",
                argc, argc == 1 ? "" : "s");
        return NULL;
    }

    Expr* expr = res->data.function.args[0];

    /* Build `f[expr]`.  Default `f` is Norm (which itself reduces to
     * Abs for scalars).  Owns its own copy of both head and arg so the
     * input `res` stays intact -- the evaluator frees it for us. */
    Expr* fhead;
    if (argc == 2) {
        fhead = expr_copy(res->data.function.args[1]);
    } else {
        fhead = expr_new_symbol("Norm");
    }
    Expr* norm_call = expr_new_function(
        fhead, (Expr*[]){ expr_copy(expr) }, 1);
    Expr* norm_val = eval_and_free(norm_call);

    if (norm_is_numeric_zero(norm_val)) {
        expr_free(norm_val);
        return expr_copy(expr);
    }

    /* Times[expr, Power[norm_val, -1]].  Times is Listable so this
     * threads the reciprocal across every leaf of a tensor `expr`. */
    Expr* inv = eval_and_free(expr_new_function(
        expr_new_symbol("Power"),
        (Expr*[]){ norm_val, expr_new_integer(-1) }, 2));

    return eval_and_free(expr_new_function(
        expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(expr), inv }, 2));
}
