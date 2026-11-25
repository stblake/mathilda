#include "list_common.h"
#include "flatten_at.h"
#include "part.h"
#include "ndarray.h"
#include "message.h"
#include "sym_intern.h"   /* intern_symbol: the splice marker */
#include "assoc.h"        /* is_association */

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
 * The addressed subexpression is replaced by a private splice marker holding
 * its arguments, and a second pass splices every marker into its parent.
 * Because each targeted slot holds a single marker node during the walk,
 * arg_count never changes mid-walk, so several positions resolve against the
 * original expression with no index-shift bookkeeping (no descending sort,
 * unlike Insert / Delete).  The splice is done here, not left to the
 * evaluator's Sequence flattening: inside HoldComplete that flattening never
 * happens (FlattenAt[HoldComplete[{a, b}], 1] is HoldComplete[a, b], not
 * HoldComplete[Sequence[a, b]]; audit FlattenAt#1), and a marker at the root
 * -- position {} -- has no parent to splice into at all (FlattenAt#2).
 */

/* The private head of a pending splice (see the note above). */
static const char* flatten_at_marker(void) {
    static const char* m = NULL;
    if (!m) m = intern_symbol("FlattenAt`Private`splice");
    return m;
}

/* Splice every marker node of `e` (owned) into its parent's arguments.  A
 * subtree with no marker is handed back as it is, so only the paths the walk
 * rebuilt are rebuilt again. */
static Expr* flatten_at_splice(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return e;
    const char* mk = flatten_at_marker();
    size_t n = e->data.function.arg_count;
    Expr** kids = malloc(sizeof(Expr*) * (n ? n : 1));
    bool changed = false;
    size_t total = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* a = e->data.function.args[i];
        kids[i] = flatten_at_splice(expr_copy(a));
        if (kids[i] != a) changed = true;
        if (head_is(kids[i], mk)) { changed = true; total += kids[i]->data.function.arg_count; }
        else total++;
    }
    if (!changed) {
        for (size_t i = 0; i < n; i++) expr_free(kids[i]);
        free(kids);
        return e;
    }
    Expr** out = malloc(sizeof(Expr*) * (total ? total : 1));
    size_t k = 0;
    for (size_t i = 0; i < n; i++) {
        if (head_is(kids[i], mk)) {
            for (size_t j = 0; j < kids[i]->data.function.arg_count; j++)
                out[k++] = expr_copy(kids[i]->data.function.args[j]);
            expr_free(kids[i]);
        } else {
            out[k++] = kids[i];
        }
    }
    free(kids);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), out, total);
    free(out);
    expr_free(e);
    return r;
}

/* The walk's result with its markers spliced; NULL (after a message) when the
 * root itself was addressed. */
static Expr* flatten_at_finish(Expr* walked, Expr* res) {
    if (!walked) return NULL;
    if (head_is(walked, flatten_at_marker())) {
        char* cs = expr_to_string(res);
        mth_message("FlattenAt", "psl", "Position {} in %s names the whole expression, "
                    "which has no parent to be flattened into.", cs ? cs : "");
        free(cs);
        expr_free(walked);
        return NULL;
    }
    return flatten_at_splice(walked);
}

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
 * with the splice marker [a1, ..., ak], dropping its head; flatten_at_splice
 * splices it into the parent. `leaf` is borrowed (never freed here); the returned
 * node is newly owned, per the PartLeafFn contract. */
static Expr* flatten_at_leaf(void* ctx, Expr* leaf) {
    (void)ctx;
    if (flatten_at_atomic(leaf)) return expr_copy(leaf);

    size_t n = leaf->data.function.arg_count;
    Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        args[i] = expr_copy(leaf->data.function.args[i]);
    Expr* seq = expr_new_function(expr_new_symbol(flatten_at_marker()), args, n);
    free(args);   /* elements + head adopted; the array itself is ours to free */
    return seq;
}

/* Does the integer path p (n indices) exist in e?  A non-integer index is
 * not judged here (true). */
static bool flatten_at_path_exists(const Expr* e, Expr** p, size_t n) {
    for (size_t d = 0; d < n; d++) {
        if (p[d]->type != EXPR_INTEGER) return true;
        if (flatten_at_atomic(e) || is_association((Expr*)e)) return true;
        int64_t len = (int64_t)e->data.function.arg_count, k = p[d]->data.integer;
        if (k < 0) k = len + k + 1;
        if (k < 1 || k > len) return false;
        e = e->data.function.args[k - 1];
    }
    return true;
}

static void flatten_at_partw_path(const Expr* e, Expr* path) {
    char* ps = expr_to_string(path);
    char* es = expr_to_string((Expr*)e);
    mth_message("FlattenAt", "partw", "Part %s of %s does not exist.", ps ? ps : "", es ? es : "");
    free(ps); free(es);
}

/* The walker declined: name the first position that does not exist,
 * FlattenAt::partw ("Part {5} of {a, b} does not exist."; audit FlattenAt#3). */
static void flatten_at_partw(const Expr* e, Expr* pos) {
    if (pos->type == EXPR_INTEGER) {
        if (flatten_at_path_exists(e, &pos, 1)) return;
        Expr* one[1] = { expr_copy(pos) };
        Expr* path = expr_new_function(expr_new_symbol(SYM_List), one, 1);
        flatten_at_partw_path(e, path);
        expr_free(path);
        return;
    }
    if (!head_is(pos, SYM_List) || pos->data.function.arg_count == 0) return;
    if (!head_is(pos->data.function.args[0], SYM_List)) {
        if (!flatten_at_path_exists(e, pos->data.function.args, pos->data.function.arg_count))
            flatten_at_partw_path(e, pos);
        return;
    }
    for (size_t i = 0; i < pos->data.function.arg_count; i++) {
        Expr* p = pos->data.function.args[i];
        if (!head_is(p, SYM_List)) return;
        if (!flatten_at_path_exists(e, p->data.function.args, p->data.function.arg_count)) {
            flatten_at_partw_path(e, p);
            return;
        }
    }
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
        return flatten_at_finish(out, res);
    }

    Expr* out = expr_apply_at_positions(expr, pos, flatten_at_leaf, NULL);
    if (!out) flatten_at_partw(expr, pos);
    return flatten_at_finish(out, res);
}
