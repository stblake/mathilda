#include "list_common.h"
#include "ndarray.h"    /* is_ndarray */
#include "ndstruct.h"   /* ndstruct_delist_repack — packed-argument fallback */
#include "rotate.h"
#include "linkedlist_seq.h" /* LLSeq: a LinkedList rotates from its slots */
#include "message.h"        /* mth_message: ::rspec */
#include "print.h"          /* expr_to_string */

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

/* `right`: rotate the other way.  The amount is reduced modulo the length
 * BEFORE it is negated, so -2^63 rotates right correctly (audit RotateRight#1:
 * negating the int64 spec overflowed). */
static Expr* rotate_rec(Expr* expr, Expr* n_spec, size_t level_idx, bool right) {
    if (expr->type != EXPR_FUNCTION && expr->type != EXPR_LINKEDLIST) return expr_copy(expr);

    int64_t n = 0;
    if (n_spec->type == EXPR_INTEGER) {
        if (level_idx == 0) n = n_spec->data.integer;
    } else if (n_spec->type == EXPR_FUNCTION && n_spec->data.function.head->data.symbol.name == SYM_List) {
        if (level_idx < n_spec->data.function.arg_count) {
            Expr* sub_n = n_spec->data.function.args[level_idx];
            if (sub_n->type == EXPR_INTEGER) n = sub_n->data.integer;
        }
    }

    if (expr->type == EXPR_LINKEDLIST && n_spec->type == EXPR_INTEGER) {
        /* One level: the two ranges of the storage, copied span-wise (typed
         * slots move without boxing). */
        const LinkedList* src = expr->data.llist;
        int64_t len = ll_length(src);
        LinkedList* out = ll_new();
        if (len > 0) {
            int64_t offset = n % len;
            if (right) offset = -offset;
            if (offset < 0) offset += len;
            ll_append_range(out, src, offset, len);
            ll_append_range(out, src, 0, offset);
        }
        return expr_new_linkedlist(out);
    }

    LLSeq lv;
    llq_open(expr, &lv);
    size_t len = lv.n;
    Expr** new_args = malloc(sizeof(Expr*) * (len ? len : 1));

    if (len > 0) {
        int64_t offset = n % (int64_t)len;
        if (right) offset = -offset;
        if (offset < 0) offset += (int64_t)len;

        for (size_t i = 0; i < len; i++) {
            size_t old_idx = (i + (size_t)offset) % len;
            new_args[i] = rotate_rec(llq_peek(&lv, old_idx), n_spec, level_idx + 1, right);
        }
    }

    Expr* result = llq_result(&lv, new_args, len);
    free(new_args);
    llq_close(&lv);
    return result;
}

/* RotateLeft/RotateRight[expr, n]: an atom is Head::normal; an explicit
 * number that is not an integer is Head::rspec. A symbolic amount stays
 * unevaluated silently (see rotate_spec_ok). */
static bool rotate_args_ok(Expr* res, const char* head) {
    Expr* expr = res->data.function.args[0];
    Expr* n_spec = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    if (expr_is_plain_atom(expr)) { builtin_normal_error(head, res, 1); return false; }
    if (rotate_spec_ok(n_spec)) return true;
    bool numeric = expr_is_plain_atom(n_spec) && n_spec->type != EXPR_SYMBOL &&
                   n_spec->type != EXPR_STRING && n_spec->type != EXPR_BIGINT;
    if (numeric) {
        char* ns = expr_to_string(n_spec);
        mth_message(head, "rspec", "Rotation specification %s should be an integer or a "
                    "list of integers.", ns ? ns : "");
        free(ns);
    }
    return false;
}

static Expr* rotate_entry(Expr* res, bool right) {
    const char* head = right ? "RotateRight" : "RotateLeft";
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* expr = res->data.function.args[0];
    if (!rotate_args_ok(res, head)) return NULL;   /* symbolic amount: stay unevaluated */
    /* Native buffer rotate first: a rotation permutes contiguous blocks, so it is
     * memcpy work. ndstruct_delist_repack below is still the fallback for a spec
     * outside the fast domain, but it costs one Expr per element -- 42.6 ms on a
     * 512x512 matrix, which dominated every stencil written with RotateLeft. */
    if (is_ndarray(expr)) {
        Expr* fast = ndstruct_rotate(res, !right);
        if (fast) return fast;
        return ndstruct_delist_repack(res, expr);
    }
    Expr* n_spec = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    Expr* default_n = NULL;
    if (!n_spec) {
        default_n = expr_new_integer(1);
        n_spec = default_n;
    }
    Expr* ret = rotate_rec(expr, n_spec, 0, right);
    if (default_n) expr_free(default_n);
    return ret;
}

Expr* builtin_rotateleft(Expr* res) { return rotate_entry(res, false); }

Expr* builtin_rotateright(Expr* res) { return rotate_entry(res, true); }
