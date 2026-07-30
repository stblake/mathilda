#include "list_common.h"
#include "ndarray.h"    /* is_ndarray */
#include "ndstruct.h"   /* ndstruct_delist_repack — packed-argument fallback */
#include "rotate.h"

/*
 * Is the rotation spec usable -- an Integer, or a List of Integers?
 *
 * rotate_rec defaults its per-level amount to 0 and only overwrites it for an
 * EXPR_INTEGER, so a SYMBOLIC spec silently rotated by zero and returned the
 * input unchanged: RotateLeft[{1,2,3}, i] answered {1,2,3} where Mathematica
 * (and Mathilda's own RotateRight, which declines) leave it unevaluated.
 *
 * That is a wrong answer on its own, and it propagates. Sum's closed-form stage
 * analyses its body for dependence on the iterator; with RotateLeft[b, {i,j}]
 * collapsing to b the body looks CONSTANT in i and j, so
 *     Sum[RotateLeft[b,{i,j}], {i,-1,1}, {j,-1,1}]
 * returned 9*b instead of the 3x3 neighbourhood sum -- which is the whole of the
 * vectorised Game of Life benchmark, computing something else entirely, slowly.
 * Mathematica gives {6,6,6} for Sum[RotateLeft[{1,2,3},i],{i,0,2}]; Mathilda
 * gave {3,6,9}.
 */
static bool rotate_spec_ok(const Expr* n_spec) {
    if (!n_spec) return true;                       /* absent: defaults to 1 */
    if (n_spec->type == EXPR_INTEGER) return true;
    if (n_spec->type == EXPR_FUNCTION
        && n_spec->data.function.head->type == EXPR_SYMBOL
        && n_spec->data.function.head->data.symbol.name == SYM_List) {
        for (size_t i = 0; i < n_spec->data.function.arg_count; i++)
            if (n_spec->data.function.args[i]->type != EXPR_INTEGER) return false;
        return true;
    }
    return false;
}

static Expr* rotate_rec(Expr* expr, Expr* n_spec, size_t level_idx) {
    if (expr->type != EXPR_FUNCTION) return expr_copy(expr);

    int64_t n = 0;
    if (n_spec->type == EXPR_INTEGER) {
        if (level_idx == 0) n = n_spec->data.integer;
    } else if (n_spec->type == EXPR_FUNCTION && n_spec->data.function.head->data.symbol.name == SYM_List) {
        if (level_idx < n_spec->data.function.arg_count) {
            Expr* sub_n = n_spec->data.function.args[level_idx];
            if (sub_n->type == EXPR_INTEGER) n = sub_n->data.integer;
        }
    }

    size_t len = expr->data.function.arg_count;
    Expr** new_args = malloc(sizeof(Expr*) * len);

    if (len > 0) {
        int64_t offset = n % (int64_t)len;
        if (offset < 0) offset += (int64_t)len;

        for (size_t i = 0; i < len; i++) {
            size_t old_idx = (i + (size_t)offset) % len;
            new_args[i] = rotate_rec(expr->data.function.args[old_idx], n_spec, level_idx + 1);
        }
    }

    Expr* result = expr_new_function(expr_copy(expr->data.function.head), new_args, len);
    free(new_args);
    return result;
}

Expr* builtin_rotateleft(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* n_spec = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    if (!rotate_spec_ok(n_spec)) return NULL;   /* symbolic amount: stay unevaluated */
    /* Native buffer rotate first: a rotation permutes contiguous blocks, so it is
     * memcpy work. ndstruct_delist_repack below is still the fallback for a spec
     * outside the fast domain, but it costs one Expr per element -- 42.6 ms on a
     * 512x512 matrix, which dominated every stencil written with RotateLeft. */
    if (is_ndarray(expr)) {
        Expr* fast = ndstruct_rotate(res, true);
        if (fast) return fast;
        return ndstruct_delist_repack(res, expr);
    }

    Expr* default_n = NULL;
    if (!n_spec) {
        default_n = expr_new_integer(1);
        n_spec = default_n;
    }

    Expr* ret = rotate_rec(expr, n_spec, 0);
    if (default_n) expr_free(default_n);
    return ret;
}

Expr* builtin_rotateright(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* n_spec = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    if (!rotate_spec_ok(n_spec)) return NULL;   /* symbolic amount: stay unevaluated */
    /* Native buffer rotate first; see the note in builtin_rotateleft. */
    if (is_ndarray(expr)) {
        Expr* fast = ndstruct_rotate(res, false);
        if (fast) return fast;
        return ndstruct_delist_repack(res, expr);
    }

    Expr* neg_n_spec = NULL;
    if (!n_spec) {
        neg_n_spec = expr_new_integer(-1);
    } else if (n_spec->type == EXPR_INTEGER) {
        neg_n_spec = expr_new_integer(-n_spec->data.integer);
    } else if (n_spec->type == EXPR_FUNCTION && n_spec->data.function.head->data.symbol.name == SYM_List) {
        Expr** neg_args = malloc(sizeof(Expr*) * n_spec->data.function.arg_count);
        for (size_t i = 0; i < n_spec->data.function.arg_count; i++) {
            if (n_spec->data.function.args[i]->type == EXPR_INTEGER) {
                neg_args[i] = expr_new_integer(-n_spec->data.function.args[i]->data.integer);
            } else {
                neg_args[i] = expr_copy(n_spec->data.function.args[i]);
            }
        }
        neg_n_spec = expr_new_function(expr_new_symbol(SYM_List), neg_args, n_spec->data.function.arg_count);
        free(neg_args);
    } else {
        return NULL;
    }

    Expr* ret = rotate_rec(expr, neg_n_spec, 0);
    expr_free(neg_n_spec);
    return ret;
}
