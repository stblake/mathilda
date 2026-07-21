/* NDArray structural fast paths — see ndstruct.h. Each op copies from the input
 * buffer into a fresh output buffer (inputs are immutable during evaluation)
 * and returns a new EXPR_NDARRAY; anything outside the fast domain degrades to
 * ndarray_delist_and_reeval so the result matches the List path exactly. */

#include "ndstruct.h"
#include "ndarray.h"
#include "ndarray_internal.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

/* Product of dims[lo..hi). */
static size_t nd_prod(const int64_t* dims, int lo, int hi) {
    size_t p = 1;
    for (int i = lo; i < hi; i++) p *= (size_t)dims[i];
    return p;
}

/* An Integer/Real Expr as a double (Clip bounds); false for anything else. */
static bool nd_real_value(const Expr* e, double* out) {
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real; return true; }
    return false;
}

bool ndstruct_call_has_ndarray(const Expr* res) {
    return res && res->type == EXPR_FUNCTION &&
           res->data.function.arg_count >= 1 &&
           is_ndarray(res->data.function.args[0]);
}

/* ------------------------------------------------------------------- Sort */

/* Sort[a]: real rank-1 only — sort the values ascending and rebuild an NDArray
 * of the same dtype (float32 values round-trip through double exactly). A custom
 * ordering function (2-arg), a complex dtype (canonical complex order differs
 * from numeric order), or rank >= 2 (Sort orders whole rows by canonical Expr
 * order) all degrade to the List path. */
Expr* ndstruct_sort(Expr* res) {
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    NDType dt = a->data.ndarray.dtype;
    if (a->data.ndarray.rank != 1 || ndt_is_complex(dt))
        return ndarray_delist_and_reeval(res);

    size_t n = (size_t)a->data.ndarray.dims[0];
    double* s = malloc(sizeof(double) * n);
    if (!s) return ndarray_delist_and_reeval(res);
    nd_gather_real(a->data.ndarray.data, dt, 0, 1, n, s);
    nd_sort_ascending(s, n);

    void* out = malloc(ndt_elem_size(dt) * n);
    if (!out) { free(s); return ndarray_delist_and_reeval(res); }
    for (size_t i = 0; i < n; i++) ndt_set(out, i, dt, s[i], 0.0);
    free(s);
    int64_t dims[1] = { (int64_t)n };
    return expr_new_ndarray(1, dims, out, dt);
}

/* ---------------------------------------------------------------- Reverse */

/* Reverse[a]: reverse the order of the leading-axis rows (each row a contiguous
 * trailing block). A level/axis spec (2-arg) degrades. */
Expr* ndstruct_reverse(Expr* res) {
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;
    NDType dt = a->data.ndarray.dtype;
    const char* buf = (const char*)a->data.ndarray.data;

    size_t blocks = (size_t)dims[0];
    size_t rowbytes = nd_prod(dims, 1, rank) * ndt_elem_size(dt);
    char* out = malloc(rowbytes * blocks);
    if (!out) return ndarray_delist_and_reeval(res);
    for (size_t b = 0; b < blocks; b++)
        memcpy(out + b * rowbytes, buf + (blocks - 1 - b) * rowbytes, rowbytes);
    return expr_new_ndarray(rank, dims, out, dt);
}

/* -------------------------------------------------------------- Transpose */

/* Transpose[a] for a rank-2 array: out[c, r] = in[r, c]. Higher rank or an
 * explicit permutation spec degrades. */
Expr* ndstruct_transpose(Expr* res) {
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    if (a->data.ndarray.rank != 2) return ndarray_delist_and_reeval(res);
    const int64_t* dims = a->data.ndarray.dims;
    NDType dt = a->data.ndarray.dtype;
    const void* buf = a->data.ndarray.data;
    size_t R = (size_t)dims[0], C = (size_t)dims[1];

    void* out = malloc(ndt_elem_size(dt) * R * C);
    if (!out) return ndarray_delist_and_reeval(res);
    for (size_t r = 0; r < R; r++)
        for (size_t c = 0; c < C; c++) {
            double re, im;
            ndt_get(buf, r * C + c, dt, &re, &im);
            ndt_set(out, c * R + r, dt, re, im);
        }
    int64_t odims[2] = { (int64_t)C, (int64_t)R };
    return expr_new_ndarray(2, odims, out, dt);
}

/* ---------------------------------------------------------------- Flatten */

/* Flatten[a]: collapse every axis into a single rank-1 array. The buffer is
 * already contiguous row-major, so this is a straight copy + reshape. A partial
 * Flatten[a, n] (or a head arg) degrades. */
Expr* ndstruct_flatten(Expr* res) {
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    NDType dt = a->data.ndarray.dtype;
    size_t sz = ndarray_size(a);
    size_t bytes = ndt_elem_size(dt) * sz;
    void* out = malloc(bytes);
    if (!out) return ndarray_delist_and_reeval(res);
    memcpy(out, a->data.ndarray.data, bytes);
    int64_t odims[1] = { (int64_t)sz };
    return expr_new_ndarray(1, odims, out, dt);
}

/* ------------------------------------------------------------- Take / Drop */

/* Build a new NDArray of `count` leading rows starting at row `start`, copied
 * from `a` (rows are contiguous trailing blocks). */
static Expr* nd_rows(const Expr* a, size_t start, size_t count) {
    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;
    NDType dt = a->data.ndarray.dtype;
    size_t rowbytes = nd_prod(dims, 1, rank) * ndt_elem_size(dt);
    char* out = malloc(rowbytes * count);
    if (!out) return NULL;
    memcpy(out, (const char*)a->data.ndarray.data + start * rowbytes, rowbytes * count);
    int64_t odims[NDARRAY_MAX_RANK];
    for (int i = 0; i < rank; i++) odims[i] = dims[i];
    odims[0] = (int64_t)count;
    return expr_new_ndarray(rank, odims, out, dt);
}

/* Resolve a leading-axis spec into a 1-based inclusive [lo, hi] range over
 * `blocks` rows. Handles an integer n (n>0: first n; n<0: last |n|) and a list
 * {i} / {i, j} with 1-based, possibly-negative endpoints. Returns false (caller
 * degrades) for anything else or an out-of-range / empty selection. */
static bool nd_span(const Expr* spec, size_t blocks, int64_t* lo, int64_t* hi) {
    if (spec->type == EXPR_INTEGER) {
        int64_t n = spec->data.integer;
        if (n > 0)      { *lo = 1; *hi = n; }
        else if (n < 0) { *lo = (int64_t)blocks + n + 1; *hi = (int64_t)blocks; }
        else return false;                            /* n == 0 */
    } else if (spec->type == EXPR_FUNCTION &&
               spec->data.function.head->type == EXPR_SYMBOL &&
               spec->data.function.head->data.symbol.name == SYM_List &&
               (spec->data.function.arg_count == 1 || spec->data.function.arg_count == 2)) {
        Expr* e0 = spec->data.function.args[0];
        Expr* e1 = spec->data.function.args[spec->data.function.arg_count == 2 ? 1 : 0];
        if (e0->type != EXPR_INTEGER || e1->type != EXPR_INTEGER) return false;
        int64_t i = e0->data.integer, j = e1->data.integer;
        *lo = (i > 0) ? i : (int64_t)blocks + i + 1;
        *hi = (j > 0) ? j : (int64_t)blocks + j + 1;
    } else {
        return false;
    }
    if (*lo < 1 || *hi > (int64_t)blocks || *lo > *hi) return false;
    return true;
}

Expr* ndstruct_take(Expr* res) {
    if (res->data.function.arg_count != 2) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    size_t blocks = (size_t)a->data.ndarray.dims[0];
    int64_t lo, hi;
    if (!nd_span(res->data.function.args[1], blocks, &lo, &hi))
        return ndarray_delist_and_reeval(res);
    Expr* r = nd_rows(a, (size_t)(lo - 1), (size_t)(hi - lo + 1));
    return r ? r : ndarray_delist_and_reeval(res);
}

Expr* ndstruct_drop(Expr* res) {
    if (res->data.function.arg_count != 2) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;
    NDType dt = a->data.ndarray.dtype;
    size_t blocks = (size_t)dims[0];
    int64_t lo, hi;
    if (!nd_span(res->data.function.args[1], blocks, &lo, &hi))
        return ndarray_delist_and_reeval(res);

    /* Keep rows [1, lo) and (hi, blocks] (1-based); drop the middle span. */
    size_t before = (size_t)(lo - 1);
    size_t after = blocks - (size_t)hi;
    size_t kept = before + after;
    if (kept == 0) return ndarray_delist_and_reeval(res);   /* empty -> List path */
    size_t rowbytes = nd_prod(dims, 1, rank) * ndt_elem_size(dt);
    const char* buf = (const char*)a->data.ndarray.data;
    char* out = malloc(rowbytes * kept);
    if (!out) return ndarray_delist_and_reeval(res);
    memcpy(out, buf, rowbytes * before);
    memcpy(out + rowbytes * before, buf + rowbytes * (size_t)hi, rowbytes * after);
    int64_t odims[NDARRAY_MAX_RANK];
    for (int i = 0; i < rank; i++) odims[i] = dims[i];
    odims[0] = (int64_t)kept;
    return expr_new_ndarray(rank, odims, out, dt);
}

/* ------------------------------------------------------------------- Clip */

/* Clip[a] clamps to [-1, 1]; Clip[a, {min, max}] to [min, max]. Elementwise,
 * real dtype only. The 3-arg replacement form and complex dtypes degrade. */
Expr* ndstruct_clip(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    NDType dt = a->data.ndarray.dtype;
    if (ndt_is_complex(dt)) return ndarray_delist_and_reeval(res);

    double lo = -1.0, hi = 1.0;
    if (argc == 2) {
        Expr* iv = res->data.function.args[1];
        if (iv->type != EXPR_FUNCTION ||
            iv->data.function.head->type != EXPR_SYMBOL ||
            iv->data.function.head->data.symbol.name != SYM_List ||
            iv->data.function.arg_count != 2 ||
            !nd_real_value(iv->data.function.args[0], &lo) ||
            !nd_real_value(iv->data.function.args[1], &hi))
            return ndarray_delist_and_reeval(res);
    }

    const void* buf = a->data.ndarray.data;
    size_t sz = ndarray_size(a);
    void* out = malloc(ndt_elem_size(dt) * sz);
    if (!out) return ndarray_delist_and_reeval(res);
    for (size_t k = 0; k < sz; k++) {
        double r, im;
        ndt_get(buf, k, dt, &r, &im);
        if (r < lo) r = lo; else if (r > hi) r = hi;
        ndt_set(out, k, dt, r, 0.0);
    }
    return expr_new_ndarray(a->data.ndarray.rank, a->data.ndarray.dims, out, dt);
}
