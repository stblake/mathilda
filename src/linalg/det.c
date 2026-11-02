/* Determinant by Laplace expansion.
 *
 * laplace_det is also used by builtin_cross (in cross.c) to evaluate the
 * minors of the generalised cross product, so it is exported via linalg.h.
 */

#include "linalg.h"
#include "ndlinalg.h"
#include "eval.h"
#include "internal.h"   /* internal_expand: Det answers expanded */
#include "common.h"     /* head_is */
#include "print.h"
#include "sym_names.h"
#include "flint_mat_bridge.h"
#include "ludecomp_internal.h"
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

/* Minor memo for one top-level expansion.  Rows are consumed in order, so a
 * minor is fixed by the set of columns it keeps: a bitmask (original_n <= 32).
 * Cofactor expansion recomputed every minor once per path to it -- n! work --
 * and a symbolic 10x10 block-diagonal matrix (FeynCalc's FCFeynmanPrepare
 * check) took 1.6-4 s; with the memo a dense n x n costs n 2^n products, and
 * an exact-zero entry skips its minor altogether. */
typedef struct { uint32_t* keys; Expr** vals; size_t cap, n; } LapMemo;

static Expr** lap_slot(LapMemo* m, uint32_t key) {
    size_t i = ((size_t)key * 2654435761u) & (m->cap - 1);
    while (m->vals[i] && m->keys[i] != key) i = (i + 1) & (m->cap - 1);
    m->keys[i] = key;
    return &m->vals[i];
}

static void lap_grow(LapMemo* m) {
    LapMemo g = { calloc(m->cap * 2, sizeof(uint32_t)), calloc(m->cap * 2, sizeof(Expr*)), m->cap * 2, m->n };
    for (size_t i = 0; i < m->cap; i++)
        if (m->vals[i]) *lap_slot(&g, m->keys[i]) = m->vals[i];
    free(m->keys); free(m->vals);
    *m = g;
}

static bool lap_is_zero(const Expr* e) {
    return e->type == EXPR_INTEGER && e->data.integer == 0;
}

static Expr* laplace_rec(Expr** flat, int original_n, int n, int row, int* cols, LapMemo* memo) {
    if (n == 1) {
        return expr_copy(flat[row * original_n + cols[0]]);
    }
    uint32_t key = 0;
    if (memo) {
        for (int j = 0; j < n; j++) key |= (uint32_t)1 << cols[j];
        Expr** hit = lap_slot(memo, key);
        if (*hit) return expr_copy(*hit);
    }
    Expr** sum_args = malloc(sizeof(Expr*) * n);
    int* new_cols = malloc(sizeof(int) * (n - 1));
    size_t m = 0;
    for (int i = 0; i < n; i++) {
        Expr* elem = flat[row * original_n + cols[i]];
        if (lap_is_zero(elem)) continue;          /* 0 * minor: skip the minor */
        for (int j = 0, k = 0; j < n; j++) {
            if (j != i) new_cols[k++] = cols[j];
        }
        Expr* minor_det = laplace_rec(flat, original_n, n - 1, row + 1, new_cols, memo);
        if (i % 2 != 0) {
            Expr* t_args[3] = { expr_new_integer(-1), expr_copy(elem), minor_det };
            sum_args[m++] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), t_args, 3));
        } else {
            Expr* t_args[2] = { expr_copy(elem), minor_det };
            sum_args[m++] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), t_args, 2));
        }
    }
    free(new_cols);
    Expr* res = (m == 0) ? expr_new_integer(0)
              : eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), sum_args, m));
    free(sum_args);
    if (memo) {
        if (2 * (memo->n + 1) > memo->cap) lap_grow(memo);
        *lap_slot(memo, key) = expr_copy(res);
        memo->n++;
    }
    return res;
}

Expr* laplace_det(Expr** flat, int original_n, int n, int row, int* cols) {
    if (n < 4 || original_n > 32) return laplace_rec(flat, original_n, n, row, cols, NULL);
    LapMemo memo = { calloc(64, sizeof(uint32_t)), calloc(64, sizeof(Expr*)), 64, 0 };
    Expr* r = laplace_rec(flat, original_n, n, row, cols, &memo);
    for (size_t i = 0; i < memo.cap; i++)
        if (memo.vals[i]) expr_free(memo.vals[i]);
    free(memo.keys); free(memo.vals);
    return r;
}

/* Does e contain a Plus anywhere? */
static bool det_contains_plus(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Plus)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (det_contains_plus(e->data.function.args[i])) return true;
    return false;
}

/* Does e carry a negative power of a sum -- a denominator Expand would
 * distribute every numerator over? (1/k, c/b and 1/Sin[x] do not count.) */
static bool det_has_sum_denominator(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* x = e->data.function.args[1];
        bool neg = (x->type == EXPR_INTEGER && x->data.integer < 0) ||
                   (x->type == EXPR_BIGINT && mpz_sgn(x->data.bigint) < 0) ||
                   (head_is(x, SYM_Rational) && x->data.function.args[0]->type == EXPR_INTEGER &&
                    x->data.function.args[0]->data.integer < 0);
        if (neg && det_contains_plus(e->data.function.args[0])) return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (det_has_sum_denominator(e->data.function.args[i])) return true;
    return false;
}

Expr* builtin_det(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    if (linalg_call_has_ndarray(res)) return ndla_det(res);

    /* A machine-real matrix too small to have been packed reaches the same fast
     * path here. Without it a 6x6 Det runs the generic path: 1.28 s over 5000
     * repeats against 0.0012 s packed, for the identical value. */
    Expr* packed_call = linalg_pack_machine_call1(res);
    if (packed_call) {
        Expr* fast = linalg_try_fast(ndla_det, packed_call);
        expr_free(packed_call);
        if (fast) return fast;
    }

    Expr* arg = res->data.function.args[0];

    int64_t dims[64];
    int rank = get_tensor_dims(arg, dims);

    if (rank != 2 || dims[0] != dims[1] || dims[0] == 0) {
        char* arg_str = expr_to_string(arg);
        mth_message("Det", "matsq", "Argument %s at position 1 is not a non-empty square matrix.", arg_str);
        free(arg_str);
        return NULL;
    }

    int n = (int)dims[0];
    Expr** flat = malloc(sizeof(Expr*) * n * n);
    size_t idx = 0;
    flatten_tensor(arg, flat, &idx);

    /* FLINT fast path: an all-integer / all-rational matrix has an exact
     * determinant computed by fmpq_mat_det in polynomial time, versus the
     * O(n!) Laplace expansion below. Returns NULL (falling through) when any
     * entry is symbolic, or when built without FLINT. The exact value matches
     * the Laplace result. */
    Expr* fdet = flint_mat_det(flat, n);
    if (fdet) {
        for (size_t i = 0; i < idx; i++) expr_free(flat[i]);
        free(flat);
        return fdet;
    }

    /* Inexact-numeric fast path: an MPFR matrix (or a machine-real matrix that
     * did not pack) gets its determinant from an O(n^3) LU factorisation rather
     * than the O(n!) Laplace expansion below, which hangs for n >= ~12.  The
     * wide MPFR exponent also means the pivot product never spuriously
     * overflows.  Declines (NULL) for exact or symbolic input, which keep their
     * existing FLINT / Laplace paths. */
    Expr* mdet = mpfr_det_dispatch(arg, n);
    if (mdet) {
        for (size_t i = 0; i < idx; i++) expr_free(flat[i]);
        free(flat);
        return mdet;
    }

    int* cols = malloc(sizeof(int) * n);
    for (int i = 0; i < n; i++) cols[i] = i;

    Expr* det_val = laplace_det(flat, n, n, 0, cols);

    free(cols);
    for (size_t i = 0; i < idx; i++) expr_free(flat[i]);
    free(flat);

    /* The reference's Det of a symbolic matrix is the expanded polynomial
     * (Det[{{a, b, c}, ...}] is -(c e g) + b f g + ..., Det[{{a + b, c},
     * {d, e}}] is -(c d) + a e + b e); the cofactor expansion leaves it
     * nested, and FeynCalc's EpsContract squares an Eps through Det with no
     * Expand after it.  A matrix with a sum in a denominator keeps the
     * cofactor form: Expand distributes every numerator over it, and a 6x6
     * matrix of rational functions overflows (Together reduces the nested
     * form at no such cost).  Only Det itself expands: laplace_det's other callers
     * (Cross, the adjugate in Inverse) keep their forms. */
    if (det_val && det_val->type == EXPR_FUNCTION && !det_has_sum_denominator(arg)) {
        Expr* ex = eval_and_free(internal_expand((Expr*[]){ det_val }, 1));
        det_val = ex;
    }
    return det_val;
}
