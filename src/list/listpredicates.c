#include "list_common.h"
#include "listpredicates.h"
#include "ndarray.h"
#include "sparse/sparse.h"
#include "linkedlist.h"

/* A sequence level: a List, or a LinkedList (which stands for one). */
static bool is_seq(const Expr* e) { return is_listq((Expr*)e) || e->type == EXPR_LINKEDLIST; }

/* test[e] is True. */
static bool lp_test(Expr* test, Expr* e) {
    Expr* call_args[1] = { expr_copy(e) };
    Expr* call = expr_new_function(expr_copy(test), call_args, 1);
    Expr* r = evaluate(call);
    bool t = (r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_True);
    expr_free(r);
    expr_free(call);
    return t;
}

Expr* builtin_listq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    if (is_listq(arg)) {
        return expr_new_symbol(SYM_True);
    }
    return expr_new_symbol(SYM_False);
}

Expr* builtin_vectorq(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count != 1 && res->data.function.arg_count != 2)) return NULL;
    Expr* arg = res->data.function.args[0];
    if (arg->type == EXPR_SPARSE) {
        bool ok = arg->data.sparse->rank == 1;
        if (ok && res->data.function.arg_count == 2) ok = sp_elements_pass(arg, res->data.function.args[1]);
        return expr_new_symbol(ok ? SYM_True : SYM_False);
    }
    if (arg->type == EXPR_NDARRAY) {
        /* NDArray elements are always machine reals, so a user-supplied
         * test is only meaningful if it accepts reals; no test means the
         * O(1) rank check alone decides it. */
        if (res->data.function.arg_count != 2) {
            return expr_new_symbol(arg->data.ndarray.rank == 1 ? SYM_True : SYM_False);
        }
        if (arg->data.ndarray.rank != 1) return expr_new_symbol(SYM_False);
        Expr* test = res->data.function.args[1];
        for (int64_t i = 0; i < arg->data.ndarray.dims[0]; i++) {
            Expr* call_args[1] = { ndarray_element_to_expr(arg, (size_t)i) };
            Expr* call = expr_new_function(expr_copy(test), call_args, 1);
            Expr* eval_res = evaluate(call);
            bool is_true = (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol.name == SYM_True);
            expr_free(eval_res);
            expr_free(call);
            if (!is_true) return expr_new_symbol(SYM_False);
        }
        return expr_new_symbol(SYM_True);
    }
    if (arg->type == EXPR_LINKEDLIST) {
        /* a vector when no element is a sequence (and each passes test) */
        Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
        LinkedList* snap = ll_copy(arg->data.llist);
        int64_t n = ll_length(snap);
        bool ok = true;
        if (test || ll_kind(snap) == LL_BOXED) {
            for (int64_t i = 0; i < n && ok; i++) {
                bool own;
                Expr* el = ll_elem_view(snap, i, &own);
                ok = test ? lp_test(test, el) : !is_seq(el);
                ll_elem_done(el, own);
            }
        }
        ll_release(snap);
        return expr_new_symbol(ok ? SYM_True : SYM_False);
    }
    if (!is_listq(arg)) return expr_new_symbol(SYM_False);

    Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;

    for (size_t i = 0; i < arg->data.function.arg_count; i++) {
        Expr* elem = arg->data.function.args[i];
        if (test == NULL) {
            if (is_seq(elem)) return expr_new_symbol(SYM_False);
        } else {
            Expr* call_args[1] = { expr_copy(elem) };
            Expr* call = expr_new_function(expr_copy(test), call_args, 1);
            Expr* eval_res = evaluate(call);
            bool is_true = (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol.name == SYM_True);
            expr_free(eval_res);
            expr_free(call);
            if (!is_true) return expr_new_symbol(SYM_False);
        }
    }
    return expr_new_symbol(SYM_True);
}

Expr* builtin_matrixq(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count != 1 && res->data.function.arg_count != 2)) return NULL;
    Expr* arg = res->data.function.args[0];
    if (arg->type == EXPR_SPARSE) {
        bool ok = arg->data.sparse->rank == 2;
        if (ok && res->data.function.arg_count == 2) ok = sp_elements_pass(arg, res->data.function.args[1]);
        return expr_new_symbol(ok ? SYM_True : SYM_False);
    }
    if (arg->type == EXPR_NDARRAY) {
        if (res->data.function.arg_count != 2) {
            return expr_new_symbol(arg->data.ndarray.rank == 2 ? SYM_True : SYM_False);
        }
        if (arg->data.ndarray.rank != 2) return expr_new_symbol(SYM_False);
        Expr* test = res->data.function.args[1];
        size_t n = ndarray_size(arg);
        for (size_t i = 0; i < n; i++) {
            Expr* call_args[1] = { ndarray_element_to_expr(arg, i) };
            Expr* call = expr_new_function(expr_copy(test), call_args, 1);
            Expr* eval_res = evaluate(call);
            bool is_true = (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol.name == SYM_True);
            expr_free(eval_res);
            expr_free(call);
            if (!is_true) return expr_new_symbol(SYM_False);
        }
        return expr_new_symbol(SYM_True);
    }
    if (arg->type == EXPR_LINKEDLIST || (is_listq(arg) && arg->data.function.arg_count > 0)) {
        /* rows that are sequences (Lists or LinkedLists) of one length, of
         * non-sequence elements; read from the slots */
        int64_t nr = ll_seq_length(arg);
        if (nr <= 0) return expr_new_symbol(SYM_False);
        bool any_ll = arg->type == EXPR_LINKEDLIST;
        for (int64_t i = 0; i < nr && !any_ll; i++)
            any_ll = arg->data.function.args[i]->type == EXPR_LINKEDLIST;
        if (any_ll) {
            Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
            Expr* src = arg->type == EXPR_LINKEDLIST ? expr_new_linkedlist(ll_copy(arg->data.llist))
                                                     : expr_copy(arg);
            int64_t cols = -1;
            bool ok = true;
            for (int64_t i = 0; i < nr && ok; i++) {
                Expr* row = ll_seq_element(src, i);
                if (row->type == EXPR_LINKEDLIST) {
                    Expr* r2 = expr_new_linkedlist(ll_copy(row->data.llist));
                    expr_free(row);
                    row = r2;
                }
                int64_t nc = is_seq(row) ? ll_seq_length(row) : -1;
                if (nc < 0 || (cols >= 0 && nc != cols)) ok = false;
                cols = nc;
                for (int64_t j = 0; j < nc && ok; j++) {
                    Expr* el = ll_seq_element(row, j);
                    ok = test ? lp_test(test, el) : !is_seq(el);
                    expr_free(el);
                }
                expr_free(row);
            }
            expr_free(src);
            return expr_new_symbol(ok ? SYM_True : SYM_False);
        }
    }
    if (!is_listq(arg)) return expr_new_symbol(SYM_False);

    if (arg->data.function.arg_count == 0) return expr_new_symbol(SYM_False);

    Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;

    size_t col_count = 0;
    bool first_row = true;

    for (size_t i = 0; i < arg->data.function.arg_count; i++) {
        Expr* row = arg->data.function.args[i];
        if (!is_listq(row)) return expr_new_symbol(SYM_False);

        if (first_row) {
            col_count = row->data.function.arg_count;
            first_row = false;
        } else {
            if (row->data.function.arg_count != col_count) return expr_new_symbol(SYM_False);
        }

        for (size_t j = 0; j < row->data.function.arg_count; j++) {
            Expr* elem = row->data.function.args[j];
            if (test == NULL) {
                if (is_listq(elem)) return expr_new_symbol(SYM_False);
            } else {
                Expr* call_args[1] = { expr_copy(elem) };
                Expr* call = expr_new_function(expr_copy(test), call_args, 1);
                Expr* eval_res = evaluate(call);
                bool is_true = (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol.name == SYM_True);
                expr_free(eval_res);
                expr_free(call);
                if (!is_true) return expr_new_symbol(SYM_False);
            }
        }
    }

    return expr_new_symbol(SYM_True);
}
