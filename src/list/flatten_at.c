#include "list_common.h"
#include "flatten_at.h"
#include "part.h"
#include "ndarray.h"

/*
 * FlattenAt[list, n]                     splice element n into its parent.
 * FlattenAt[expr, {i, j, ...}]           splice the part at one deep position.
 * FlattenAt[expr, {{i1, ...}, {i2, ...}}] splice at several positions.
 *
 * FlattenAt removes the head of the subexpression at each position and splices
 * that subexpression's arguments into the surrounding expression -- it works on
 * any head, not just List (FlattenAt[f[g[1, 2], g[3, 4]], 1] is f[1, 2, g[3, 4]]).
 *
 * Position resolution -- integers (negatives from the end), a single deep path
 * {i, j, ...}, and the list-of-paths form {{i1, ...}, {i2, ...}} -- is handled
 * by the shared walker expr_apply_at_positions (src/part.c), the same one MapAt
 * and ReplaceAt use, so the {2} (one path) vs {{2}, {4}} (two paths) distinction
 * is exactly theirs. A position that does not exist makes the walker return NULL
 * and FlattenAt stay unevaluated, matching Part[{a, b, c}, 5].
 *
 * The splice itself is deferred to the evaluator: the addressed subexpression is
 * replaced by Sequence @@ leaf, and flatten_sequences (src/eval.c) splices that
 * Sequence into the parent on the next evaluation pass -- exactly as MapAt
 * returns an unreduced f[leaf]. Because each targeted slot holds a single
 * Sequence node during the walk, arg_count never changes mid-walk, so several
 * positions resolve against the original expression with no index-shift
 * bookkeeping (no descending sort, unlike Insert / Delete).
 */

/* Mirror of part.c's is_atomic: Rational and Complex are stored as EXPR_FUNCTION
 * nodes but are atomic leaves with no arguments to splice, so FlattenAt[{1/2, x},
 * 1] must not manufacture Sequence[1, 2]. Everything that is not an EXPR_FUNCTION
 * -- symbols, numbers, strings, a packed NDArray -- is atomic too. */
static bool flatten_at_atomic(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return true;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_Complex || h == SYM_Rational) return true;
    }
    return false;
}

/* Leaf action for FlattenAt: replace the addressed subexpression g[a1, ..., ak]
 * with Sequence[a1, ..., ak], dropping its head; the evaluator splices the
 * Sequence into the parent. `leaf` is borrowed (never freed here); the returned
 * node is newly owned, per the PartLeafFn contract. */
static Expr* flatten_at_leaf(void* ctx, Expr* leaf) {
    (void)ctx;
    if (flatten_at_atomic(leaf)) return expr_copy(leaf);

    size_t n = leaf->data.function.arg_count;
    Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        args[i] = expr_copy(leaf->data.function.args[i]);
    Expr* seq = expr_new_function(expr_new_symbol(SYM_Sequence), args, n);
    free(args);   /* elements + head adopted; the array itself is ours to free */
    return seq;
}

Expr* builtin_flatten_at(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2)
        return builtin_arg_error("FlattenAt", argc, 1, 2);

    Expr* expr = res->data.function.args[0];

    /* FlattenAt[expr]: no positions given, so nothing is flattened. */
    if (argc == 1) return expr_copy(expr);

    Expr* pos = res->data.function.args[1];

    /* A visible NDArray is atomic, so the walker cannot descend into it.
     * Materialize to a nested list and flatten there. The result is ragged --
     * flattening one row of {{1, 2}, {3, 4}} gives {1, 2, {3, 4}} -- so it must
     * never be repacked (unlike MapAt, whose result can stay rectangular). */
    if (is_ndarray(expr)) {
        Expr* nested = ndarray_to_nested_list(expr);
        Expr* out = expr_apply_at_positions(nested, pos, flatten_at_leaf, NULL);
        expr_free(nested);
        return out;
    }

    return expr_apply_at_positions(expr, pos, flatten_at_leaf, NULL);
}
