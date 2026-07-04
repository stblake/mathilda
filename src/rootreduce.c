/* Mathilda — RootReduce implementation.
 *
 * RootReduce[expr] rewrites an expression that lives in an algebraic tower
 *
 *     K = Q(params)(alpha_1, ..., alpha_r)
 *
 * — where each generator alpha_l is a radical Power[base, p/q] of any index
 * q >= 2 (cube roots and higher included), or a root of unity (-1)^(p/q) — into
 * a canonical, rationalised representative: a reduced polynomial in the radicals
 * over a denominator that is free of every radical / root-of-unity generator.
 *
 * The heavy lifting is the rigorous FLINT field-arithmetic engine
 * `flint_algebraic_field_canonical` (src/poly/flint_bridge.c): K is a
 * finite-dimensional Q(params)-vector space, multiplication-by-the-denominator
 * is a Q(params)-linear map (each product reduced to normal form modulo the
 * generators' minimal-polynomial ideal via fmpz_mpoly_divrem_ideal), and
 * inverting that map with a linear solve over Q(params) (gr_mat over
 * fmpz_mpoly_q) rationalises the denominator to Norm_{K/Q(params)}(D). No
 * numeric zero oracle is consulted — the reduction is exact and complete when
 * the minimal polynomials are irreducible.
 *
 * When `expr` carries no algebraic generator (or the case is out of the
 * engine's scope), RootReduce leaves it unchanged, matching the Wolfram Language
 * behaviour of returning the input for expressions it cannot canonicalise.
 */

#include "rootreduce.h"

#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "symtab.h"
#include "flint_bridge.h"

/* RootReduce[expr] — see file header. Not HoldAll: `expr` is already evaluated.
 * Ownership follows the builtin contract (return a new tree or steal from res;
 * never expr_free(res)). */
Expr* builtin_rootreduce(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

#ifdef USE_FLINT
    Expr* r = flint_algebraic_field_canonical(arg);
    if (r) return r;
#endif

    /* No canonicalisation available: RootReduce[expr] -> expr. Steal arg out of
     * res (the evaluator frees the res wrapper) to avoid a copy + double-free. */
    res->data.function.args[0] = NULL;
    return arg;
}

void rootreduce_init(void) {
    symtab_add_builtin("RootReduce", builtin_rootreduce);
    symtab_get_def("RootReduce")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_set_docstring("RootReduce",
        "RootReduce[expr] rewrites expr over the algebraic tower "
        "Q(params)(radicals) in canonical form, rationalising radical or "
        "root-of-unity denominators via rigorous FLINT field arithmetic. Leaves "
        "expr unchanged when it has no algebraic generator.");
}
