/* jacobisymbol.c -- JacobiSymbol[].
 * Split from numbertheory.c; see numbertheory.h and
 * numbertheory_internal.h for the subsystem layout. */

#include "numbertheory.h"
#include "numbertheory_internal.h"
#include "arithmetic.h"
#include "expr.h"
#include <stdio.h>
#include <gmp.h>

/* ---- JacobiSymbol ------------------------------------------------------- */

/* Emit the wrong-argument-count diagnostic for JacobiSymbol and return NULL so
 * the call is left unevaluated. */
static Expr* jacobisymbol_emit_argrx(size_t npos) {
    if (!arith_warnings_muted()) {
        fprintf(stderr,
                "JacobiSymbol::argrx: JacobiSymbol called with %zu argument%s; "
                "2 arguments are expected.\n",
                npos, npos == 1 ? "" : "s");
    }
    return NULL;
}

/* JacobiSymbol[n, m] gives the Jacobi symbol (n/m).
 *
 * Following the Wolfram Language, this is the full Kronecker-symbol
 * generalisation of the Jacobi symbol: m need not be odd or positive, and n
 * may be negative.  For prime m it reduces to the Legendre symbol, which is
 * +-1 according to whether n is a quadratic residue modulo m (and 0 when m
 * divides n).  GMP's mpz_kronecker computes exactly this and is O((log m)^2),
 * so machine integers and arbitrary-precision bigints are handled uniformly.
 *
 * Symbolic arguments leave the call unevaluated; Listable threading over
 * lists is applied by the evaluator before this builtin is reached. */
Expr* builtin_jacobisymbol(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return jacobisymbol_emit_argrx(argc);

    Expr* n_expr = res->data.function.args[0];
    Expr* m_expr = res->data.function.args[1];
    if (!expr_is_integer_like(n_expr) || !expr_is_integer_like(m_expr))
        return NULL;

    mpz_t n, m;
    expr_to_mpz(n_expr, n);
    expr_to_mpz(m_expr, m);
    int k = mpz_kronecker(n, m);   /* full Kronecker: handles even/neg/zero m */
    mpz_clears(n, m, NULL);
    return expr_new_integer(k);    /* always -1, 0, or 1 */
}
