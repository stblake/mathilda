/* UpperTriangularize[m] / UpperTriangularize[m, k] and the Lower forms, for
 * dense matrices (List, packed List, NDArray).
 *
 * UpperTriangularize keeps the entries on and above the k-th diagonal
 * (column - row >= k) and LowerTriangularize those on and below it
 * (column - row <= k); every other entry becomes zero. k defaults to 0, may be
 * any machine integer, and the matrix may be rectangular.
 *
 * THE ZERO. A matrix holding a machine number anywhere gets the machine zero
 * 0. (so a Real matrix stays a Real matrix and keeps its buffer); any other
 * matrix gets the exact 0 -- including an arbitrary-precision one.
 *
 * A SparseArray argument is handled by the head's native sparse kernel
 * (src/sparse/sparse_struct.c, triangularize) before this builtin is asked.
 *
 * The buffer path writes a fresh buffer of the operand's dtype (presentation
 * inherited), one row at a time: a memcpy of the kept run and a zero fill of
 * the rest. */

#include "linalg.h"
#include "linkedlist.h"
#include "sym_names.h"
#include "common.h"
#include "pack.h"
#include "ndarray.h"
#include "message.h"
#include "print.h"
#include <stdlib.h>
#include <string.h>

static bool is_list(const Expr* e) { return head_is(e, SYM_List); }

/* Does `e` hold a machine Real anywhere (a machine Complex has Real parts)? */
static bool has_machine(const Expr* e) {
    if (e->type == EXPR_REAL) return true;
    if (e->type == EXPR_NDARRAY) {
        NDType dt = e->data.ndarray.dtype;
        return dt != NDT_INT64 && dt != NDT_BOOL;
    }
    if (e->type != EXPR_FUNCTION) return false;
    if (has_machine(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_machine(e->data.function.args[i])) return true;
    return false;
}

static void msg_with_call(const char* head, const char* tag, const char* fmt, const Expr* e) {
    char* s = expr_to_string_outputform((Expr*)e);
    mth_message(head, tag, fmt, s ? s : "");
    free(s);
}

/* The kept column range [lo, hi) of row r in a C-column matrix. */
static void kept(int64_t r, int64_t C, int64_t k, bool upper, int64_t* lo, int64_t* hi) {
    int64_t cut = r + k;                         /* |k| <= R + C: no overflow */
    if (upper) { *lo = cut < 0 ? 0 : cut > C ? C : cut; *hi = C; }
    else       { *lo = 0; *hi = cut + 1 < 0 ? 0 : cut + 1 > C ? C : cut + 1; }
}

static Expr* tri_buffer(const Expr* m, int64_t k, bool upper) {
    int64_t R = m->data.ndarray.dims[0], C = m->data.ndarray.dims[1];
    NDType dt = m->data.ndarray.dtype;
    size_t esz = ndt_elem_size(dt);
    void* buf = NULL;
    Expr* out = ndbuild_open_like(m, 2, m->data.ndarray.dims, dt, &buf);
    if (!out) return NULL;
    const char* src = m->data.ndarray.data;
    char* dst = buf;
    for (int64_t r = 0; r < R; r++) {
        int64_t lo, hi;
        kept(r, C, k, upper, &lo, &hi);
        char* row = dst + (size_t)(r * C) * esz;
        /* all-zero bits are 0, 0. and 0. + 0. I in every numeric dtype */
        memset(row, 0, (size_t)C * esz);
        if (hi > lo) memcpy(row + (size_t)lo * esz, src + (size_t)(r * C + lo) * esz, (size_t)(hi - lo) * esz);
    }
    return out;
}

static Expr* tri_list(const Expr* m, int64_t k, bool upper) {
    int64_t R = (int64_t)m->data.function.arg_count;
    int64_t C = (int64_t)m->data.function.args[0]->data.function.arg_count;
    bool mach = has_machine(m);
    Expr** rows = malloc(sizeof(Expr*) * (size_t)R);
    Expr** cells = malloc(sizeof(Expr*) * (size_t)(C ? C : 1));
    for (int64_t r = 0; r < R; r++) {
        const Expr* row = m->data.function.args[r];
        int64_t lo, hi;
        kept(r, C, k, upper, &lo, &hi);
        for (int64_t c = 0; c < C; c++)
            cells[c] = (c >= lo && c < hi) ? expr_copy(row->data.function.args[c])
                     : mach ? expr_new_real(0.0) : expr_new_integer(0);
        rows[r] = expr_new_function(expr_new_symbol(SYM_List), cells, (size_t)C);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)R);
    free(cells);
    free(rows);
    return out;
}

/* A matrix with a LinkedList level (a LinkedList of rows, or LinkedList
 * rows): tri_list's rule over the leaves read from the slots, the answer
 * shaped as the argument (a LinkedList level stays a LinkedList). NULL after
 * the message tri_list's caller would give. */
static Expr* tri_ll(Expr* res, const Expr* m, bool upper, const char* head) {
    const char* notmat = "Argument %s at position 1 is not a non-empty rectangular matrix.";
    int64_t dims[3];
    int rank = ll_tensor_dims(m, dims, 3);
    if (rank == 3) return NULL;                       /* an array of rank 3 or more */
    if (rank != 2 || dims[0] == 0) {
        if (rank == 1 && dims[0] == 0) return NULL;
        msg_with_call(head, "matrix", notmat, m);
        return NULL;
    }
    int64_t k = 0;
    if (res->data.function.arg_count == 2) {
        const Expr* ke = res->data.function.args[1];
        if (ke->type != EXPR_INTEGER) {
            msg_with_call(head, "intm", "Machine-sized integer expected at position 2 in %s.", res);
            return NULL;
        }
        k = ke->data.integer;
    }
    int64_t R = dims[0], C = dims[1];
    int64_t lim = R + C + 1;
    if (k > lim) k = lim;
    if (k < -lim) k = -lim;
    Expr** leaves = malloc(sizeof(Expr*) * (size_t)(R * C > 0 ? R * C : 1));
    size_t n = 0;
    ll_tensor_leaves(m, 2, leaves, &n);
    bool mach = false;
    for (size_t t = 0; t < n && !mach; t++) mach = has_machine(leaves[t]);
    for (int64_t r = 0; r < R; r++) {
        int64_t lo, hi;
        kept(r, C, k, upper, &lo, &hi);
        for (int64_t c = 0; c < C; c++) {
            if (c >= lo && c < hi) continue;
            expr_free(leaves[r * C + c]);
            leaves[r * C + c] = mach ? expr_new_real(0.0) : expr_new_integer(0);
        }
    }
    bool lvl[2];
    ll_tensor_levels(m, 2, lvl);
    Expr* out = ll_tensor_build(dims, 2, leaves, lvl);
    free(leaves);
    return out;
}

static Expr* triangularize(Expr* res, bool upper) {
    const char* head = upper ? "UpperTriangularize" : "LowerTriangularize";
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return builtin_arg_error(head, argc, 1, 2);
    if (argc > 2) {
        if (res->data.function.args[0]->type == EXPR_NDARRAY) return ndarray_delist_and_reeval(res);
        const Expr* x = res->data.function.args[2];
        if (head_is(x, SYM_Rule) || head_is(x, SYM_RuleDelayed)) return NULL;
        char* xs = expr_to_string_outputform((Expr*)x);
        char* cs = expr_to_string_outputform(res);
        mth_message(head, "nonopt", "Options expected (instead of %s) beyond position 2 in %s. "
                    "An option must be a rule or a list of rules.", xs ? xs : "", cs ? cs : "");
        free(xs);
        free(cs);
        return NULL;
    }
    const Expr* m = res->data.function.args[0];
    if (ll_has_ll_level1(m)) return tri_ll(res, m, upper, head);

    /* The matrix: a full array of rank 2 (a deeper array stays unevaluated). */
    const char* notmat = "Argument %s at position 1 is not a non-empty rectangular matrix.";
    bool buffer = false;
    if (m->type == EXPR_NDARRAY) {
        int rk = m->data.ndarray.rank;
        if (rk > 2) return NULL;
        /* reported on the List form, once (see ndarray_delist_and_reeval) */
        if (rk < 2) return ndarray_delist_and_reeval(res);
        buffer = m->data.ndarray.dtype != NDT_BOOL;
    } else if (is_list(m)) {
        size_t R = m->data.function.arg_count;
        if (R == 0) return NULL;
        const Expr* r0 = m->data.function.args[0];
        if (!is_list(r0)) { msg_with_call(head, "matrix", notmat, m); return NULL; }
        size_t C = r0->data.function.arg_count;
        size_t sub = 0, cells = 0;
        for (size_t i = 0; i < R; i++) {
            const Expr* row = m->data.function.args[i];
            if (!is_list(row) || row->data.function.arg_count != C) {
                msg_with_call(head, "matrix", notmat, m);
                return NULL;
            }
            for (size_t j = 0; j < C; j++, cells++)
                if (is_list(row->data.function.args[j]) || row->data.function.args[j]->type == EXPR_NDARRAY)
                    sub++;
        }
        if (sub && sub == cells) return NULL;   /* an array of rank 3 or more */
        if (sub) { msg_with_call(head, "matrix", notmat, m); return NULL; }
    } else {
        msg_with_call(head, "matrix", notmat, m);
        return NULL;
    }

    int64_t k = 0;
    if (argc == 2) {
        const Expr* ke = res->data.function.args[1];
        if (ke->type != EXPR_INTEGER) {
            if (m->type == EXPR_NDARRAY) return ndarray_delist_and_reeval(res);
            msg_with_call(head, "intm", "Machine-sized integer expected at position 2 in %s.", res);
            return NULL;
        }
        k = ke->data.integer;
    }

    Expr* plain = NULL;
    if (m->type == EXPR_NDARRAY && !buffer) m = plain = ndarray_to_nested_list(m);
    int64_t R = m->type == EXPR_NDARRAY ? m->data.ndarray.dims[0] : (int64_t)m->data.function.arg_count;
    int64_t C = m->type == EXPR_NDARRAY ? m->data.ndarray.dims[1]
                                        : (int64_t)m->data.function.args[0]->data.function.arg_count;
    /* clamp k to where it still matters, so r + k cannot overflow */
    int64_t lim = R + C + 1;
    if (k > lim) k = lim;
    if (k < -lim) k = -lim;
    Expr* out = m->type == EXPR_NDARRAY ? tri_buffer(m, k, upper) : tri_list(m, k, upper);
    if (plain) expr_free(plain);
    return out;
}

Expr* builtin_uppertriangularize(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    return triangularize(res, true);
}

Expr* builtin_lowertriangularize(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    return triangularize(res, false);
}
