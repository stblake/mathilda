/*
 * sparse_reduce.c -- native reductions and comparisons over SparseArrays:
 * Total, Mean, Variance, StandardDeviation, Tr, Max, Min, MinMax, Norm,
 * Accumulate, Equal and Unequal.
 *
 * Each is registered as the head's sparse_native kernel (see sparse.h): it
 * borrows the call and returns the answer, or NULL to decline, in which case
 * the transparency gate (sparse_gate.c) answers on the dense form exactly as
 * before. None of them densifies: a reduction reads the explicit entries once
 * and accounts for the background arithmetically -- every implicit position
 * of a reduced cell contributes the background once, so a cell holding `c`
 * explicit entries out of `size` positions adds bg * (size - c).
 *
 * EXACTNESS. A reduced value is exact exactly when every element it combines
 * is exact, as on the dense array: Total[SparseArray[{{1, 1} -> 1}, {2, 2},
 * 0.5]] is {1.5, 1.} but the column holding only explicit Integers stays an
 * Integer. One rule follows the sparse representation rather than the dense
 * one: over machine values and the exact background 0, a cell that holds no
 * explicit entry is the machine 0. (the background is summed into the same
 * machine accumulator), so Total[SparseArray[{{1, 1} -> 1.}, {2, 2}]] is
 * {1., 0.}.
 *
 * ANSWER KIND. Total and Mean give a SparseArray when the result has rank >= 2
 * and the background is the exact 0 (its explicit positions are the cells
 * that hold an explicit entry; nothing is re-minimised), and a dense List or
 * a scalar otherwise. Every other head here answers densely.
 *
 * The typed tier covers Integer and Real value buffers (SPV_INT64, SPV_F64)
 * over an Integer or Real background, with overflow-checked integer
 * accumulation, threaded by partitions of the entries for large arrays. For
 * Total and Mean everything else (symbolic, complex, exact non-integer
 * values or background) is answered on the dense form when that is small
 * enough -- so a sum mixing symbols and machine numbers rounds exactly as the
 * dense one -- and otherwise by a boxed tier that folds one level at a time.
 * The machine 2-norm of a matrix runs the Lanczos iteration on A^T A (or
 * A A^T) with CSR / CSC matrix-vector products; an exact matrix 2-norm is
 * declined.
 */

#include "equal_tol.h"
#include "sparse_internal.h"
#include "ndarray_internal.h"
#include "sym_names.h"
#include "symtab.h"
#include "eval.h"
#include "pack.h"
#include "arithmetic.h"
#include "checked_int.h"
#include "plus.h"
#include <gmp.h>
#ifdef USE_LAPACK
#include "lapack.h"     /* cblas_ddot */
#endif
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Largest dense answer (element count) built here; above it the call is
 * declined (the gate then refuses it with SparseArray::dense, as before). */
#define SR_DENSE_MAX ((int64_t)1 << 27)

/* ======================================================================== */
/* Small helpers                                                            */
/* ======================================================================== */

static bool head_is(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym;
}

static size_t nargs(const Expr* call) { return call->data.function.arg_count; }
static const Expr* arg(const Expr* call, size_t i) { return call->data.function.args[i]; }

static Expr* mk_call(const char* head, Expr** args, size_t n) {
    return expr_new_function(expr_new_symbol(head), args, n);
}
static Expr* mk_call1(const char* head, Expr* a) {
    Expr* x[1] = { a };
    return mk_call(head, x, 1);
}
static Expr* mk_call2(const char* head, Expr* a, Expr* b) {
    Expr* x[2] = { a, b };
    return mk_call(head, x, 2);
}
static Expr* eval_free(Expr* e) {
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

static bool is_infinity(const Expr* e) {
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == SYM_Infinity;
    return head_is(e, SYM_DirectedInfinity) && e->data.function.arg_count == 1 &&
           e->data.function.args[0]->type == EXPR_INTEGER &&
           e->data.function.args[0]->data.integer == 1;
}

static bool is_exact_zero(const Expr* e) {
    return e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* A machine-typed background: an Integer or a Real. */
typedef struct {
    bool    real;   /* EXPR_REAL (else EXPR_INTEGER) */
    int64_t i;
    double  d;
} Bg;

static bool bg_typed(const Expr* e, Bg* b) {
    if (e->type == EXPR_INTEGER) { b->real = false; b->i = e->data.integer; b->d = (double)b->i; return true; }
    if (e->type == EXPR_REAL) { b->real = true; b->i = 0; b->d = e->data.real; return true; }
    return false;
}

static bool typed_kind(SpValKind k) { return k == SPV_INT64 || k == SPV_F64; }

/* z = v, for any int64 (mpz_set_si takes a long, 32 bits on LLP64). */
static void z_set_i64(mpz_t z, int64_t v) {
    uint64_t m = v < 0 ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
    mpz_import(z, 1, 1, sizeof(uint64_t), 0, 0, &m);
    if (v < 0) mpz_neg(z, z);
}

/* The 0-based coordinate d of entry k (in row r). */
static inline int64_t coord_of(const SparseArray* sa, int64_t r, int64_t k, int d) {
    if (sa->rank == 1) return sp_col(sa->pat, k, 0);
    return d == 0 ? r : sp_col(sa->pat, k, d - 1);
}

/* The rows whose first entry lies in the entry range [lo, hi) -- a partition of
 * the rows that balances by entries, for nd_parallel_for over nnz. Trailing
 * empty rows (rowptr[r] == nnz) go to the chunk that ends at nnz. */
static void rows_of_chunk(const SparsePattern* p, size_t lo, size_t hi, int64_t* r0, int64_t* r1) {
    int64_t nnz = p->nnz;
    int64_t a = 0, b = p->nrows;           /* first row with rowptr >= lo */
    while (a < b) { int64_t m = a + (b - a) / 2; if (p->rowptr[m] < (int64_t)lo) a = m + 1; else b = m; }
    *r0 = (lo == 0) ? 0 : a;
    if ((int64_t)hi >= nnz) { *r1 = p->nrows; return; }
    a = 0; b = p->nrows;
    while (a < b) { int64_t m = a + (b - a) / 2; if (p->rowptr[m] < (int64_t)hi) a = m + 1; else b = m; }
    *r1 = a;
}

/* x . y over n doubles: the platform BLAS when there is one (vectorised),
 * else four independent accumulators. */
static double vdot(const double* x, const double* y, int64_t n) {
#ifdef USE_LAPACK
    if (n > 0 && n < INT32_MAX) return cblas_ddot((int)n, x, 1, y, 1);
#endif
    double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
    int64_t k = 0;
    for (; k + 4 <= n; k += 4) { s0 += x[k] * y[k]; s1 += x[k + 1] * y[k + 1]; s2 += x[k + 2] * y[k + 2]; s3 += x[k + 3] * y[k + 3]; }
    for (; k < n; k++) s0 += x[k] * y[k];
    return (s0 + s1) + (s2 + s3);
}

/* Pairwise sum of a contiguous double range (error ~ O(log n) eps). */
static double pw_sum(const double* x, size_t n) {
    if (n <= 128) {
        double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
        size_t i = 0;
        for (; i + 4 <= n; i += 4) { s0 += x[i]; s1 += x[i + 1]; s2 += x[i + 2]; s3 += x[i + 3]; }
        for (; i < n; i++) s0 += x[i];
        return (s0 + s1) + (s2 + s3);
    }
    size_t h = n / 2;
    return pw_sum(x, h) + pw_sum(x + h, n - h);
}

/* Run fn(ctx, i) for the partitions i = 0..T-1 on up to T threads. Built on
 * nd_parallel_for over a virtual range of T * NDARRAY_THREAD_THRESHOLD, so the
 * thread count follows T; a partition belongs to the chunk holding its start,
 * so each runs exactly once. True iff every partition returned true. */
typedef bool (*sr_part_fn)(void* ctx, int part);
typedef struct { sr_part_fn fn; void* ctx; int T; } SrParts;

static bool parts_chunk(void* vctx, size_t lo, size_t hi) {
    SrParts* c = vctx;
    bool ok = true;
    for (int i = 0; i < c->T; i++) {
        size_t st = (size_t)i * NDARRAY_THREAD_THRESHOLD;
        if (st >= lo && st < hi && !c->fn(c->ctx, i)) ok = false;
    }
    return ok;
}

static bool run_parts(int T, sr_part_fn fn, void* ctx) {
    SrParts c = { fn, ctx, T };
    return nd_parallel_for((size_t)T * NDARRAY_THREAD_THRESHOLD, parts_chunk, &c);
}

/* Partitions for `work` elements: one per ~2^17 elements, at most the
 * thread cap. */
static int n_parts(int64_t work) {
    int64_t t = work >> 17;
    if (t < 1) t = 1;
    if (t > NDARRAY_MAX_THREADS) t = NDARRAY_MAX_THREADS;
    return (int)t;
}

/* ======================================================================== */
/* Reduction plans and cells                                                */
/* ======================================================================== */

/* Value transform applied to every element before it is summed (Norm). */
enum { VT_ID = 0, VT_ABS, VT_SQR };

/* Summing the dimensions lo..hi (0-based, inclusive) of a rank-r array: the
 * result has the remaining dimensions, and each result cell sums `cellsize`
 * elements. */
typedef struct {
    const SparseArray* sa;
    int     lo, hi;
    int     orank;
    int64_t odims[SP_MAX_RANK];
    int64_t ncells;
    int64_t cellsize;
} Plan;

static bool plan_make(const SparseArray* sa, int lo, int hi, Plan* t) {
    if (lo < 0 || hi >= sa->rank || lo > hi) return false;
    for (int d = 0; d < sa->rank; d++) if (sa->dims[d] == 0) return false;
    t->sa = sa; t->lo = lo; t->hi = hi; t->orank = 0;
    int64_t nc = 1, cs = 1;
    for (int d = 0; d < sa->rank; d++) {
        int64_t n = sa->dims[d];
        if (d >= lo && d <= hi) {
            if (n != 0 && cs > INT64_MAX / n) return false;
            cs *= n;
        } else {
            t->odims[t->orank++] = n;
            if (n != 0 && nc > INT64_MAX / n) return false;
            nc *= n;
        }
    }
    t->ncells = nc;
    t->cellsize = cs;
    return true;
}

static inline int64_t cell_of(const Plan* t, int64_t r, int64_t k) {
    const SparseArray* sa = t->sa;
    int64_t idx = 0;
    for (int d = 0; d < sa->rank; d++) {
        if (d >= t->lo && d <= t->hi) continue;
        idx = idx * sa->dims[d] + coord_of(sa, r, k, d);
    }
    return idx;
}

/* Per-cell typed accumulators. Exactly one of isum / fsum is used (by the
 * value kind); cnt holds the explicit-entry count per cell. */
typedef struct {
    SpValKind kind;
    int vt;
    int64_t* isum;
    double*  fsum;
    int64_t* cnt;
    bool     borrowed;   /* isum / fsum is the caller's (zeroed) buffer */
} Acc;

/* Idempotent: the pointers are cleared. */
static void acc_free(Acc* a) {
    if (!a->borrowed) { free(a->isum); free(a->fsum); }
    free(a->cnt);
    a->isum = NULL; a->fsum = NULL; a->cnt = NULL;
}

/* The transformed value of entry k: exact (INT64, false on overflow) or machine. */
static inline bool tv_i(int vt, int64_t v, int64_t* out) {
    if (vt == VT_ID) { *out = v; return true; }
    if (vt == VT_ABS) { if (v == INT64_MIN) return false; *out = v < 0 ? -v : v; return true; }
    return !ci_mul_i64(v, v, out);
}
static inline double tv_d(int vt, double v) {
    return vt == VT_ID ? v : vt == VT_ABS ? fabs(v) : v * v;
}

/* --- parallel helpers for the two common matrix shapes ------------------- */

typedef struct {
    const SparseArray* sa;
    int vt;
    int64_t* isum;
    double*  fsum;
    int64_t* cnt;
} RowCtx;

/* Row sums (cell = row) over the rows of an entry chunk. */
static bool row_sum_chunk(void* vctx, size_t lo, size_t hi) {
    RowCtx* c = vctx;
    const SparsePattern* p = c->sa->pat;
    int64_t r0, r1;
    rows_of_chunk(p, lo, hi, &r0, &r1);
    for (int64_t r = r0; r < r1; r++) {
        int64_t a = p->rowptr[r], b = p->rowptr[r + 1];
        if (c->cnt) c->cnt[r] = b - a;
        if (c->sa->kind == SPV_F64) {
            const double* v = (const double*)c->sa->vals + a;
            int64_t len = b - a;
            if (c->vt == VT_ID && len <= 8) {
                double s0 = 0;
                for (int64_t k = 0; k < len; k++) s0 += v[k];
                c->fsum[r] = s0;
            } else if (c->vt == VT_ID) c->fsum[r] = pw_sum(v, (size_t)len);
            else {
                double s = 0;
                for (int64_t k = 0; k < b - a; k++) s += tv_d(c->vt, v[k]);
                c->fsum[r] = s;
            }
        } else {
            const int64_t* v = (const int64_t*)c->sa->vals;
            int64_t s = 0, x;
            for (int64_t k = a; k < b; k++) {
                if (!tv_i(c->vt, v[k], &x) || ci_add_i64(s, x, &s)) return false;
            }
            c->isum[r] = s;
        }
    }
    return true;
}

typedef struct {
    const SparseArray* sa;
    int vt, T;
    double part[NDARRAY_MAX_THREADS];
} AllCtx;

/* Partition i of the whole-buffer sum (pairwise within a partition). */
static bool all_sum_part(void* vctx, int i) {
    AllCtx* c = vctx;
    int64_t nnz = c->sa->pat->nnz;
    size_t lo = (size_t)(nnz * i / c->T), hi = (size_t)(nnz * (i + 1) / c->T);
    const double* v = (const double*)c->sa->vals;
    if (c->vt == VT_ID) { c->part[i] = pw_sum(v + lo, hi - lo); return true; }
    double s = 0;
    for (size_t k = lo; k < hi; k++) s += tv_d(c->vt, v[k]);
    c->part[i] = s;
    return true;
}

/* Scatter entries [k0, k1) of a matrix into column sums f / is (and counts). */
static bool col_scatter_range(const SparseArray* sa, int vt, int64_t k0, int64_t k1,
                              double* f, int64_t* is, int64_t* cnt) {
    const SparsePattern* p = sa->pat;
    if (sa->kind == SPV_F64) {
        const double* v = sa->vals;
        if (p->idx32) {
            const int32_t* ci = p->colidx;
            if (vt == VT_ID) for (int64_t k = k0; k < k1; k++) f[ci[k]] += v[k];
            else for (int64_t k = k0; k < k1; k++) f[ci[k]] += tv_d(vt, v[k]);
            if (cnt) for (int64_t k = k0; k < k1; k++) cnt[ci[k]]++;
        } else {
            const int64_t* ci = p->colidx;
            for (int64_t k = k0; k < k1; k++) f[ci[k]] += tv_d(vt, v[k]);
            if (cnt) for (int64_t k = k0; k < k1; k++) cnt[ci[k]]++;
        }
        return true;
    }
    const int64_t* v = sa->vals;
    int64_t x;
    for (int64_t k = k0; k < k1; k++) {
        int64_t c = sp_col(p, k, 0);
        if (!tv_i(vt, v[k], &x) || ci_add_i64(is[c], x, &is[c])) return false;
        if (cnt) cnt[c]++;
    }
    return true;
}

typedef struct {
    const SparseArray* sa;
    int vt, T;
    int64_t ncols;
    void* a;                         /* Acc* (declared below) */
    double*  fb[NDARRAY_MAX_THREADS];
    int64_t* ib[NDARRAY_MAX_THREADS];
    int64_t* cb[NDARRAY_MAX_THREADS];
} ColCtx;

static bool col_scatter_part(void* vctx, int i) {
    ColCtx* c = vctx;
    int64_t nnz = c->sa->pat->nnz;
    int64_t k0 = nnz * i / c->T, k1 = nnz * (i + 1) / c->T;
    return col_scatter_range(c->sa, c->vt, k0, k1, c->fb[i], c->ib[i], c->cb[i]);
}

/* Merge the partitions' columns [j0, j1) into partition 0's buffers. */
static bool col_merge_part(void* vctx, int i) {
    ColCtx* c = vctx;
    int64_t j0 = c->ncols * i / c->T, j1 = c->ncols * (i + 1) / c->T;
    for (int q = 1; q < c->T; q++) {
        if (c->fb[0]) { double* d = c->fb[0]; const double* s = c->fb[q]; for (int64_t j = j0; j < j1; j++) d[j] += s[j]; }
        else {
            int64_t* d = c->ib[0];
            const int64_t* s = c->ib[q];
            for (int64_t j = j0; j < j1; j++) if (ci_add_i64(d[j], s[j], &d[j])) return false;
        }
        if (c->cb[0]) { int64_t* d = c->cb[0]; const int64_t* s = c->cb[q]; for (int64_t j = j0; j < j1; j++) d[j] += s[j]; }
    }
    return true;
}

/* Accumulate the explicit entries of `t->sa` into per-cell sums. `need_cnt`
 * asks for per-cell explicit counts. False on integer overflow or allocation
 * failure. */
static bool accumulate_into(const Plan* t, int vt, bool need_cnt, Acc* a, void* zeroed) {
    const SparseArray* sa = t->sa;
    const SparsePattern* p = sa->pat;
    size_t nc = (size_t)(t->ncells ? t->ncells : 1);
    memset(a, 0, sizeof *a);
    a->kind = sa->kind;
    a->vt = vt;
    a->borrowed = zeroed != NULL;
    if (sa->kind == SPV_INT64) a->isum = zeroed ? zeroed : calloc(nc, sizeof(int64_t));
    else a->fsum = zeroed ? zeroed : calloc(nc, sizeof(double));
    if (need_cnt) a->cnt = calloc(nc, sizeof(int64_t));
    if ((!a->isum && !a->fsum) || (need_cnt && !a->cnt)) { acc_free(a); return false; }
    int64_t nnz = p->nnz;

    /* Every element into one cell. */
    if (t->ncells == 1) {
        if (a->cnt) a->cnt[0] = nnz;
        if (nnz == 0) {
            /* calloc'd zero sums */
        } else if (sa->kind == SPV_F64) {
            AllCtx c = { sa, vt, n_parts(nnz), { 0 } };
            run_parts(c.T, all_sum_part, &c);
            double s = 0;
            for (int i = 0; i < c.T; i++) s += c.part[i];
            a->fsum[0] = s;
        } else {
            const int64_t* v = sa->vals;
            int64_t s = 0, x;
            for (int64_t k = 0; k < nnz; k++)
                if (!tv_i(vt, v[k], &x) || ci_add_i64(s, x, &s)) { acc_free(a); return false; }
            a->isum[0] = s;
        }
        return true;
    }
    /* Matrix row sums: one cell per row, in parallel by entries. */
    if (sa->rank == 2 && t->lo == 1) {
        RowCtx c = { sa, vt, a->isum, a->fsum, a->cnt };
        if (!nd_parallel_for((size_t)nnz, row_sum_chunk, &c)) { acc_free(a); return false; }
        return true;
    }
    /* Matrix column sums: cell = column. Large arrays scatter into private
     * per-partition columns and merge in partition order. */
    if (sa->rank == 2 && t->lo == 0) {
        int T = n_parts(nnz);
        if (T > 1 && nnz >= 4 * t->ncells) {
            ColCtx c;
            memset(&c, 0, sizeof c);
            c.sa = sa; c.vt = vt; c.T = T; c.ncols = t->ncells; c.a = a;
            bool ok = true;
            for (int i = 1; i < T && ok; i++) {
                if (sa->kind == SPV_F64) ok = (c.fb[i] = calloc(nc, sizeof(double))) != NULL;
                else ok = (c.ib[i] = calloc(nc, sizeof(int64_t))) != NULL;
                if (ok && a->cnt) ok = (c.cb[i] = calloc(nc, sizeof(int64_t))) != NULL;
            }
            c.fb[0] = a->fsum; c.ib[0] = a->isum; c.cb[0] = a->cnt;
            ok = ok && run_parts(T, col_scatter_part, &c) && run_parts(T, col_merge_part, &c);
            for (int i = 1; i < T; i++) { free(c.fb[i]); free(c.ib[i]); free(c.cb[i]); }
            if (!ok) { acc_free(a); return false; }
            return true;
        }
        if (!col_scatter_range(sa, vt, 0, nnz, a->fsum, a->isum, a->cnt)) { acc_free(a); return false; }
        return true;
    }
    /* General: any rank, any contiguous level range. */
    for (int64_t r = 0; r < p->nrows; r++)
        for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1]; k++) {
            int64_t c = cell_of(t, r, k);
            if (sa->kind == SPV_F64) a->fsum[c] += tv_d(vt, ((const double*)sa->vals)[k]);
            else {
                int64_t x;
                if (!tv_i(vt, ((const int64_t*)sa->vals)[k], &x) || ci_add_i64(a->isum[c], x, &a->isum[c])) {
                    acc_free(a);
                    return false;
                }
            }
            if (a->cnt) a->cnt[c]++;
        }
    return true;
}

static bool accumulate(const Plan* t, int vt, bool need_cnt, Acc* a) {
    return accumulate_into(t, vt, need_cnt, a, NULL);
}

/* One finished cell: exact Integer or machine Real. */
typedef struct {
    bool    m;
    int64_t i;
    double  d;
} Cell;

/* The background, transformed. False on overflow. */
static bool bg_transform(int vt, const Bg* b, Bg* out) {
    *out = *b;
    if (b->real) { out->d = tv_d(vt, b->d); return true; }
    if (!tv_i(vt, b->i, &out->i)) return false;
    out->d = (double)out->i;
    return true;
}

/* Combine cell c's explicit sum with the background's share. False on
 * integer overflow. See the EXACTNESS note at the top of the file. */
static bool cell_finish(const Acc* a, int64_t c, int64_t cellsize, const Bg* b, Cell* out) {
    int64_t cnt = a->cnt ? a->cnt[c] : cellsize;   /* no counts: bg is exact 0 */
    int64_t k = cellsize - cnt;
    if (a->kind == SPV_F64) {
        if (cnt == 0 && !b->real && b->i != 0) {
            int64_t v;
            if (ci_mul_i64(b->i, cellsize, &v)) return false;
            out->m = false; out->i = v; return true;
        }
        out->m = true;
        out->d = a->fsum[c] + (k ? b->d * (double)k : 0.0);
        return true;
    }
    if (b->real && k > 0) {
        out->m = true;
        out->d = (double)a->isum[c] + b->d * (double)k;
        return true;
    }
    int64_t bk, s;
    if (ci_mul_i64(b->i, k, &bk) || ci_add_i64(a->isum[c], bk, &s)) return false;
    out->m = false; out->i = s;
    return true;
}

static Expr* cell_expr(const Cell* c) {
    return c->m ? expr_new_real(c->d) : expr_new_integer(c->i);
}

/* Fold a flat row-major array of Exprs into nested Lists (moves the items). */
static Expr* fold_items(Expr** items, int rank, const int64_t* dims, int level, int64_t* off) {
    int64_t n = dims[level];
    Expr** it = malloc(sizeof(Expr*) * (size_t)(n ? n : 1));
    for (int64_t i = 0; i < n; i++)
        it[i] = (level + 1 == rank) ? items[(*off)++] : fold_items(items, rank, dims, level + 1, off);
    Expr* out = mk_call(SYM_List, it, (size_t)n);
    free(it);
    return out;
}

static Expr* dense_from_items(Expr** items, int rank, const int64_t* dims) {
    int64_t off = 0;
    return fold_items(items, rank, dims, 0, &off);
}

/* A dense array (rank >= 1) of cells: packed when every cell has one type
 * (and packing accepts it), nested Lists otherwise. Divides each cell by
 * `div` (> 1: Mean) on the way. */
static Expr* dense_from_cells(const Cell* cells, int64_t n, int rank, const int64_t* dims, int64_t div) {
    bool any_m = false, any_x = false;
    for (int64_t i = 0; i < n; i++) { if (cells[i].m) any_m = true; else any_x = true; }
    if (any_m && !any_x) {
        void* buf = NULL;
        Expr* node = ndbuild_open(rank, dims, NDT_FLOAT64, &buf);
        if (node) {
            double* b = buf;
            if (div > 1) { double dd = (double)div; for (int64_t i = 0; i < n; i++) b[i] = cells[i].d / dd; }
            else for (int64_t i = 0; i < n; i++) b[i] = cells[i].d;
            return node;
        }
    } else if (any_x && !any_m && div <= 1) {
        void* buf = NULL;
        Expr* node = ndbuild_open(rank, dims, NDT_INT64, &buf);
        if (node) {
            int64_t* b = buf;
            for (int64_t i = 0; i < n; i++) b[i] = cells[i].i;
            return node;
        }
    }
    Expr** items = malloc(sizeof(Expr*) * (size_t)(n ? n : 1));
    if (!items) return NULL;
    for (int64_t i = 0; i < n; i++) {
        if (div > 1)
            items[i] = cells[i].m ? expr_new_real(cells[i].d / (double)div) : make_rational(cells[i].i, div);
        else
            items[i] = cell_expr(&cells[i]);
    }
    Expr* out = dense_from_items(items, rank, dims);
    free(items);
    return out;
}

/* ======================================================================== */
/* Grouping explicit entries by result cell                                 */
/* ======================================================================== */

/* The explicit entries ordered by result cell (stably), with the distinct
 * cells: groups g = 0..ng-1 cover order[start[g] .. start[g+1]) and all fall
 * in cell gcell[g] (increasing). */
typedef struct {
    int64_t  ng;
    int64_t* gcell;
    int64_t* start;
    int64_t* order;
} Groups;

static void groups_free(Groups* g) { free(g->gcell); free(g->start); free(g->order); }

typedef struct { int64_t cell, k; } CellKey;

static int cellkey_cmp(const void* x, const void* y) {
    const CellKey* a = x;
    const CellKey* b = y;
    if (a->cell != b->cell) return a->cell < b->cell ? -1 : 1;
    return (a->k > b->k) - (a->k < b->k);
}

static bool group_entries(const Plan* t, Groups* g) {
    const SparsePattern* p = t->sa->pat;
    int64_t nnz = p->nnz;
    memset(g, 0, sizeof *g);
    g->order = malloc(sizeof(int64_t) * (size_t)(nnz ? nnz : 1));
    g->gcell = malloc(sizeof(int64_t) * (size_t)(nnz ? nnz : 1));
    g->start = malloc(sizeof(int64_t) * (size_t)(nnz + 1));
    if (!g->order || !g->gcell || !g->start) { groups_free(g); return false; }
    CellKey* ck = malloc(sizeof(CellKey) * (size_t)(nnz ? nnz : 1));
    if (!ck) { groups_free(g); return false; }
    for (int64_t r = 0; r < p->nrows; r++)
        for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1]; k++) {
            ck[k].cell = cell_of(t, r, k);
            ck[k].k = k;
        }
    qsort(ck, (size_t)nnz, sizeof(CellKey), cellkey_cmp);
    int64_t ng = 0;
    for (int64_t i = 0; i < nnz; i++) {
        if (i == 0 || ck[i].cell != ck[i - 1].cell) {
            g->gcell[ng] = ck[i].cell;
            g->start[ng] = i;
            ng++;
        }
        g->order[i] = ck[i].k;
    }
    g->start[ng] = nnz;
    g->ng = ng;
    free(ck);
    return true;
}

/* Does `e` contain a machine or arbitrary-precision number? */
static bool has_inexact(const Expr* e) {
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_inexact(e->data.function.args[i])) return true;
    return false;
}

/* Any inexact value among the explicit entries? */
static bool any_inexact_value(const SparseArray* sa) {
    if (sa->kind == SPV_F64 || sa->kind == SPV_C128) return true;
    if (sa->kind != SPV_EXPR) return false;
    Expr* const* v = sa->vals;
    for (int64_t k = 0; k < sa->pat->nnz; k++) if (has_inexact(v[k])) return true;
    return false;
}

/* Plus over args (adopted) in the given order, exactly as Total's own level
 * fold adds them: builtin_plus first (no Orderless re-sort, so machine
 * numbers are added in element order), then a finishing evaluation. */
static Expr* plus_ordered(Expr** args, size_t n) {
    if (n == 0) return expr_new_integer(0);
    Expr* pe = mk_call(SYM_Plus, args, n);
    Expr* c = builtin_plus(pe);
    if (c) { expr_free(pe); return eval_free(c); }
    return eval_free(pe);
}

/* One fold of the dense Total adds the elements of a cell in element order
 * whichever level it folds: the innermost level hands them to Plus
 * directly, an outer one adds whole sub-arrays left to right. */
static Expr* plus_fold(Expr** args, size_t n, bool innermost) {
    (void)innermost;
    return plus_ordered(args, n);
}

/* The sum of group gi's explicit values (exact 0 background: the implicit
 * positions add nothing). */
static Expr* group_sum(const SparseArray* sa, const Groups* g, int64_t gi, bool innermost) {
    int64_t a = g->start[gi], b = g->start[gi + 1];
    Expr** it = malloc(sizeof(Expr*) * (size_t)(b - a));
    for (int64_t i = a; i < b; i++) it[i - a] = sp_value_expr(sa, g->order[i]);
    Expr* r = plus_fold(it, (size_t)(b - a), innermost);
    free(it);
    return r;
}

static Expr* scale_by(Expr* x, int64_t div) {
    if (div <= 1) return x;
    return eval_free(mk_call2(SYM_Times, x, make_rational(1, div)));
}

/* The row holding entry k. */
static int64_t entry_row(const SparsePattern* p, int64_t k) {
    int64_t a = 0, b = p->nrows;
    while (b - a > 1) { int64_t m = a + (b - a) / 2; if (p->rowptr[m] <= k) a = m; else b = m; }
    return a;
}

/* ======================================================================== */
/* Total / Mean                                                             */
/* ======================================================================== */

/* Sparse answer (rank >= 2, exact 0 background): one explicit entry per cell
 * that holds an explicit entry. NULL to decline. */
static Expr* reduce_sparse(const Plan* t, int64_t div) {
    const SparseArray* sa = t->sa;
    Groups g;
    if (!group_entries(t, &g)) return NULL;
    int64_t ng = g.ng;
    int orank = t->orank;
    int width = orank - 1;
    int64_t inner = t->ncells / (t->odims[0] ? t->odims[0] : 1);
    SparsePattern* pat = sp_pattern_new(t->odims[0], ng, width, sp_dims_idx32(orank, t->odims));
    if (!pat) { groups_free(&g); return NULL; }
    /* rowptr and tuples from the cell indices */
    for (int64_t r = 0; r <= t->odims[0]; r++) pat->rowptr[r] = 0;
    for (int64_t i = 0; i < ng; i++) {
        int64_t row = g.gcell[i] / (inner ? inner : 1);
        pat->rowptr[row + 1]++;
        int64_t rem = g.gcell[i] - row * inner;
        for (int j = width - 1; j >= 0; j--) {
            int64_t dj = t->odims[j + 1];
            int64_t cj = dj ? rem % dj : 0;
            rem = dj ? rem / dj : 0;
            if (pat->idx32) ((int32_t*)pat->colidx)[i * width + j] = (int32_t)cj;
            else ((int64_t*)pat->colidx)[i * width + j] = cj;
        }
    }
    for (int64_t r = 0; r < t->odims[0]; r++) pat->rowptr[r + 1] += pat->rowptr[r];

    SpValKind kind = SPV_EXPR;
    void* vals = NULL;
    bool ok = true;
    if (sa->kind == SPV_F64) {
        double* v = malloc(sizeof(double) * (size_t)(ng ? ng : 1));
        const double* sv = sa->vals;
        for (int64_t i = 0; i < ng; i++) {
            double s = 0;
            for (int64_t j = g.start[i]; j < g.start[i + 1]; j++) s += sv[g.order[j]];
            v[i] = div > 1 ? s / (double)div : s;
        }
        kind = SPV_F64; vals = v;
    } else if (sa->kind == SPV_INT64) {
        const int64_t* sv = sa->vals;
        int64_t* v = malloc(sizeof(int64_t) * (size_t)(ng ? ng : 1));
        for (int64_t i = 0; i < ng && ok; i++) {
            int64_t s = 0;
            for (int64_t j = g.start[i]; j < g.start[i + 1]; j++)
                if (ci_add_i64(s, sv[g.order[j]], &s)) { ok = false; break; }
            v[i] = s;
        }
        if (ok && div > 1) {
            Expr** e = malloc(sizeof(Expr*) * (size_t)(ng ? ng : 1));
            for (int64_t i = 0; i < ng; i++) e[i] = make_rational(v[i], div);
            free(v);
            kind = SPV_EXPR; vals = e;
        } else if (ok) {
            kind = SPV_INT64; vals = v;
        } else {
            free(v);
            ok = true;   /* overflow: the boxed sums below take the GMP path */
        }
    }
    if (!vals) {
        /* Boxed: element order within a cell is the dense fold's order only
         * for one level; over several, inexact sums would round differently. */
        if (t->lo < t->hi && any_inexact_value(sa)) {
            groups_free(&g);
            sp_pattern_release(pat);
            return NULL;
        }
        Expr** e = malloc(sizeof(Expr*) * (size_t)(ng ? ng : 1));
        for (int64_t i = 0; i < ng; i++) e[i] = scale_by(group_sum(sa, &g, i, t->hi == sa->rank - 1), div);
        kind = SPV_EXPR; vals = e;
    }
    groups_free(&g);
    Expr* bg = any_inexact_value(sa) ? expr_new_real(0.0) : expr_new_integer(0);
    SparseArray* out = sp_new(orank, t->odims, bg, pat, kind, vals);
    return out ? sp_to_expr(out) : NULL;
}

/* Dense / scalar answer through the boxed tier (any value kind, any
 * background), folding one level at a time from the deepest, each fold a Plus
 * over the elements in order -- the dense Total's own evaluation order, so
 * machine rounding agrees with it too. NULL to decline. */
static Expr* reduce_boxed(const Plan* t, int64_t div) {
    const SparseArray* sa = t->sa;
    const SparsePattern* p = sa->pat;
    /* first fold: dimension hi, straight from the sparse entries */
    Plan p1;
    if (!plan_make(sa, t->hi, t->hi, &p1) || p1.ncells > SR_DENSE_MAX) return NULL;
    int64_t L = sa->dims[t->hi];
    bool zero_bg = is_exact_zero(sa->bg);
    bool interleave = !zero_bg && (has_inexact(sa->bg) || any_inexact_value(sa));
    if (interleave && (p1.ncells > SR_DENSE_MAX / (L ? L : 1))) return NULL;
    Groups g;
    if (!group_entries(&p1, &g)) return NULL;
    int64_t n1 = p1.ncells;
    Expr** cur = malloc(sizeof(Expr*) * (size_t)(n1 ? n1 : 1));
    Expr** args = malloc(sizeof(Expr*) * (size_t)(L + 1));
    int64_t gi = 0;
    for (int64_t c = 0; c < n1; c++) {
        int64_t a = 0, b = 0;
        if (gi < g.ng && g.gcell[gi] == c) { a = g.start[gi]; b = g.start[gi + 1]; gi++; }
        size_t q = 0;
        if (interleave) {
            int64_t j = a;
            for (int64_t pos = 0; pos < L; pos++) {
                int64_t k = j < b ? g.order[j] : -1;
                if (k >= 0 && coord_of(sa, entry_row(p, k), k, t->hi) == pos) {
                    args[q++] = sp_value_expr(sa, k);
                    j++;
                } else {
                    args[q++] = expr_copy(sa->bg);
                }
            }
        } else {
            for (int64_t j = a; j < b; j++) args[q++] = sp_value_expr(sa, g.order[j]);
            int64_t k = L - (b - a);
            if (k > 0 && !zero_bg)
                args[q++] = mk_call2(SYM_Times, expr_new_integer(k), expr_copy(sa->bg));
        }
        cur[c] = plus_fold(args, q, t->hi == sa->rank - 1);
    }
    free(args);
    groups_free(&g);
    /* further folds over the dense intermediate, dimensions hi-1 .. lo */
    int crank = p1.orank;
    int64_t cdims[SP_MAX_RANK];
    for (int i = 0; i < crank; i++) cdims[i] = p1.odims[i];
    for (int d = t->hi - 1; d >= t->lo; d--) {
        int64_t outer = 1, inner = 1, Ld = cdims[d];
        for (int i = 0; i < d; i++) outer *= cdims[i];
        for (int i = d + 1; i < crank; i++) inner *= cdims[i];
        int64_t nn = outer * inner;
        Expr** nxt = malloc(sizeof(Expr*) * (size_t)(nn ? nn : 1));
        Expr** a2 = malloc(sizeof(Expr*) * (size_t)(Ld ? Ld : 1));
        for (int64_t o = 0; o < outer; o++)
            for (int64_t i = 0; i < inner; i++) {
                for (int64_t j = 0; j < Ld; j++) a2[j] = cur[(o * Ld + j) * inner + i];
                nxt[o * inner + i] = plus_fold(a2, (size_t)Ld, d == crank - 1);
            }
        free(a2);
        free(cur);
        cur = nxt;
        for (int i = d; i + 1 < crank; i++) cdims[i] = cdims[i + 1];
        crank--;
    }
    int64_t nout = t->ncells;
    for (int64_t c = 0; c < nout; c++) cur[c] = scale_by(cur[c], div);
    Expr* out = crank == 0 ? cur[0] : dense_from_items(cur, crank, cdims);
    free(cur);
    return out;
}

/* Total (div = 1) or Mean (div = the summed length) over the plan. NULL to
 * decline. */
/* The boxed tier for an array small enough to densify: Total / Mean on the
 * dense form itself, re-sparsified when the answer is to be sparse. Symbolic
 * cells mixing machine numbers are then rounded exactly as the dense fold
 * rounds them (its order of additions is not a simple function of the
 * element order). NULL to decline. */
static Expr* reduce_dense_eval(const Plan* t, int64_t div, bool resparse) {
    Expr* d = sp_normal(t->sa);
    if (!d) return NULL;
    Expr* call;
    if (div > 1) call = mk_call1(SYM_Mean, d);
    else {
        Expr* lv = mk_call2(SYM_List, expr_new_integer(t->lo + 1), expr_new_integer(t->hi + 1));
        call = mk_call2(SYM_Total, d, lv);
    }
    Expr* r = eval_free(call);
    if (head_is(r, SYM_Total) || head_is(r, SYM_Mean)) { expr_free(r); return NULL; }
    if (!resparse) return r;
    Expr* sp = eval_free(mk_call1(SYM_SparseArray, r));
    if (is_sparse(sp)) return sp;
    expr_free(sp);
    return NULL;
}

static Expr* reduce_plan(const Plan* t, int64_t div) {
    const SparseArray* sa = t->sa;
    int64_t total = sp_total(sa->rank, sa->dims);
    bool small = total >= 0 && total <= SR_DENSE_MAX;
    if (t->orank >= 2 && is_exact_zero(sa->bg)) {
        if (!typed_kind(sa->kind) && small) return reduce_dense_eval(t, div, true);
        return reduce_sparse(t, div);
    }
    if (t->ncells > SR_DENSE_MAX) return NULL;
    Bg b;
    if (typed_kind(sa->kind) && bg_typed(sa->bg, &b) && t->orank >= 1) {
        /* One cell type throughout -- Reals over a Real or 0 background,
         * Integers over an Integer one: sum straight into the answer's packed
         * buffer. */
        bool f64 = sa->kind == SPV_F64 && (b.real || b.i == 0);
        bool i64 = sa->kind == SPV_INT64 && !b.real && div <= 1;
        void* buf = NULL;
        Expr* node = (f64 || i64) ? ndbuild_open_zeroed(t->orank, t->odims, f64 ? NDT_FLOAT64 : NDT_INT64, &buf)
                                  : NULL;
        if (node) {
            Acc a;
            bool need_cnt = b.real || b.i != 0;
            bool ok = accumulate_into(t, VT_ID, need_cnt, &a, buf);
            int64_t n = t->ncells;
            if (ok && f64) {
                double* f = buf;
                if (need_cnt) for (int64_t c = 0; c < n; c++) f[c] += b.d * (double)(t->cellsize - a.cnt[c]);
                if (div > 1) { double dd = (double)div; for (int64_t c = 0; c < n; c++) f[c] /= dd; }
            } else if (ok) {
                int64_t* is = buf;
                for (int64_t c = 0; ok && need_cnt && c < n; c++) {
                    int64_t bk;
                    ok = !ci_mul_i64(b.i, t->cellsize - a.cnt[c], &bk) && !ci_add_i64(is[c], bk, &is[c]);
                }
            }
            acc_free(&a);
            if (ok) return node;
            expr_free(node);
        }
    }
    if (typed_kind(sa->kind) && bg_typed(sa->bg, &b)) {
        bool need_cnt = b.real || b.i != 0;
        Acc a;
        if (accumulate(t, VT_ID, need_cnt, &a)) {
            int64_t n = t->ncells;
            Cell* cells = malloc(sizeof(Cell) * (size_t)(n ? n : 1));
            bool ok = cells != NULL;
            for (int64_t c = 0; ok && c < n; c++) ok = cell_finish(&a, c, t->cellsize, &b, &cells[c]);
            acc_free(&a);
            if (ok) {
                Expr* out;
                if (t->orank == 0) {
                    if (div > 1) out = cells[0].m ? expr_new_real(cells[0].d / (double)div)
                                                  : make_rational(cells[0].i, div);
                    else out = cell_expr(&cells[0]);
                } else {
                    out = dense_from_cells(cells, n, t->orank, t->odims, div);
                }
                free(cells);
                return out;
            }
            free(cells);
        }
    }
    if (small) return reduce_dense_eval(t, div, false);
    return reduce_boxed(t, div);
}

/* A level index (1-based, negative from the end) as a 0-based dimension. */
static bool level_dim(const Expr* e, int rank, int* out) {
    if (e->type != EXPR_INTEGER) return false;
    int64_t k = e->data.integer;
    if (k < 0) k = rank + 1 + k;
    if (k < 1 || k > rank) return false;
    *out = (int)(k - 1);
    return true;
}

/* Total[s], Total[s, n], Total[s, Infinity], Total[s, {k}], Total[s, {k1, k2}]. */
static Expr* native_total(const Expr* call) {
    size_t n = nargs(call);
    if (n < 1 || n > 2 || !is_sparse(arg(call, 0))) return NULL;
    const SparseArray* sa = sp_of(arg(call, 0));
    int lo = 0, hi = 0;
    if (n == 2) {
        const Expr* spec = arg(call, 1);
        if (is_infinity(spec)) {
            hi = sa->rank - 1;
        } else if (spec->type == EXPR_INTEGER) {
            int64_t k = spec->data.integer;
            if (k == 0) return expr_copy((Expr*)arg(call, 0));
            if (k < 0) k = sa->rank + 1 + k;
            if (k < 1) return NULL;
            hi = (int)(k > sa->rank ? sa->rank : k) - 1;
        } else if (head_is(spec, SYM_List) && spec->data.function.arg_count >= 1 &&
                   spec->data.function.arg_count <= 2) {
            if (!level_dim(spec->data.function.args[0], sa->rank, &lo)) return NULL;
            hi = lo;
            if (spec->data.function.arg_count == 2 &&
                !level_dim(spec->data.function.args[1], sa->rank, &hi)) return NULL;
            if (hi < lo) return NULL;
        } else {
            return NULL;
        }
    }
    Plan t;
    if (!plan_make(sa, lo, hi, &t)) return NULL;
    return reduce_plan(&t, 1);
}

/* Mean[s]: the column means (level 1). */
static Expr* native_mean(const Expr* call) {
    if (nargs(call) != 1 || !is_sparse(arg(call, 0))) return NULL;
    const SparseArray* sa = sp_of(arg(call, 0));
    if (sa->dims[0] < 1) return NULL;
    Plan t;
    if (!plan_make(sa, 0, 0, &t)) return NULL;
    if (sa->dims[0] == 1) return reduce_plan(&t, 1);
    return reduce_plan(&t, sa->dims[0]);
}

/* ======================================================================== */
/* Variance / StandardDeviation                                             */
/* ======================================================================== */

/* Column-wise (level 1) variance of a typed vector or matrix: per cell the
 * exact (n Sx2 - Sx^2) / (n (n - 1)), or the two-pass machine value. */
static Expr* variance_core(const Expr* call, bool sd) {
    if (nargs(call) != 1 || !is_sparse(arg(call, 0))) return NULL;
    const SparseArray* sa = sp_of(arg(call, 0));
    Bg b;
    if (sa->rank > 2 || !typed_kind(sa->kind) || !bg_typed(sa->bg, &b)) return NULL;
    int64_t n = sa->dims[0];
    if (n < 2) return NULL;
    Plan t;
    if (!plan_make(sa, 0, 0, &t) || t.ncells > SR_DENSE_MAX) return NULL;
    Acc s1, s2;
    if (!accumulate(&t, VT_ID, true, &s1)) return NULL;
    bool have_s2 = sa->kind == SPV_INT64 && accumulate(&t, VT_SQR, false, &s2);
    int64_t nc = t.ncells;
    Expr** items = malloc(sizeof(Expr*) * (size_t)(nc ? nc : 1));
    double* mean = malloc(sizeof(double) * (size_t)(nc ? nc : 1));
    double* ss = calloc((size_t)(nc ? nc : 1), sizeof(double));
    bool* mach = malloc(sizeof(bool) * (size_t)(nc ? nc : 1));
    bool ok = items && mean && ss && mach;
    /* which cells are machine, and their means */
    for (int64_t c = 0; ok && c < nc; c++) {
        int64_t cnt = s1.cnt[c], k = n - cnt;
        mach[c] = sa->kind == SPV_F64 ? (cnt > 0 || b.real) : (b.real && k > 0);
        double tot = sa->kind == SPV_F64 ? s1.fsum[c] : (double)s1.isum[c];
        mean[c] = (tot + b.d * (double)k) / (double)n;
        if (mach[c] && k > 0) { double dv = b.d - mean[c]; ss[c] = dv * dv * (double)k; }
    }
    if (ok) {
        /* second pass over the explicit entries of machine cells */
        const SparsePattern* p = sa->pat;
        for (int64_t r = 0; r < p->nrows; r++)
            for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1]; k++) {
                int64_t c = cell_of(&t, r, k);
                if (!mach[c]) continue;
                double x = sa->kind == SPV_F64 ? ((const double*)sa->vals)[k]
                                               : (double)((const int64_t*)sa->vals)[k];
                double dv = x - mean[c];
                ss[c] += dv * dv;
            }
    }
    for (int64_t c = 0; ok && c < nc; c++) {
        if (mach[c]) {
            double v = ss[c] / (double)(n - 1);
            items[c] = expr_new_real(sd ? sqrt(v) : v);
            continue;
        }
        /* exact: values are Integers (INT64 cell) or only the background */
        int64_t cnt = s1.cnt[c], k = n - cnt;
        mpz_t S1, S2, tmp, num, den;
        mpz_inits(S1, S2, tmp, num, den, NULL);
        if (sa->kind == SPV_INT64) {
            if (!have_s2) ok = false;
            else { z_set_i64(S1, s1.isum[c]); z_set_i64(S2, s2.isum[c]); }
        }
        if (ok) {
            mpz_t bb, bk;
            mpz_inits(bb, bk, NULL);
            z_set_i64(bb, b.i);
            z_set_i64(bk, k);
            mpz_mul(bk, bk, bb);             /* b k */
            mpz_add(S1, S1, bk);
            mpz_mul(bk, bk, bb);             /* b^2 k */
            mpz_add(S2, S2, bk);
            z_set_i64(den, n);
            mpz_mul(num, den, S2);
            mpz_mul(tmp, S1, S1);
            mpz_sub(num, num, tmp);          /* n S2 - S1^2 */
            mpz_sub_ui(tmp, den, 1);
            mpz_mul(den, den, tmp);          /* n (n - 1) */
            Expr* v = make_rational_mpz(num, den);
            if (sd) v = eval_free(mk_call1(SYM_Sqrt, v));
            items[c] = v;
            mpz_clears(bb, bk, NULL);
        }
        mpz_clears(S1, S2, tmp, num, den, NULL);
        if (!ok) { for (int64_t j = 0; j < c; j++) expr_free(items[j]); }
    }
    acc_free(&s1);
    if (have_s2) acc_free(&s2);
    free(mean); free(ss); free(mach);
    if (!ok) { free(items); return NULL; }
    Expr* out;
    if (sa->rank == 1) out = items[0];
    else {
        out = dense_from_items(items, 1, &t.odims[0]);
        out = pack_offer(out);
    }
    free(items);
    return out;
}

static Expr* native_variance(const Expr* call) { return variance_core(call, false); }
static Expr* native_sd(const Expr* call) { return variance_core(call, true); }

/* ======================================================================== */
/* Max / Min / MinMax                                                       */
/* ======================================================================== */

/* The extreme of a typed value buffer, in parallel partitions. */
typedef struct {
    const SparseArray* sa;
    bool is_max;
    int T;
    double  fm[NDARRAY_MAX_THREADS];
    int64_t im[NDARRAY_MAX_THREADS];
    bool    nan[NDARRAY_MAX_THREADS];
} ExtCtx;

static bool ext_part(void* vctx, int i) {
    ExtCtx* c = vctx;
    int64_t nnz = c->sa->pat->nnz;
    int64_t k0 = nnz * i / c->T, k1 = nnz * (i + 1) / c->T;
    if (c->sa->kind == SPV_F64) {
        /* four independent lanes: one compare chain is latency-bound */
        const double* v = c->sa->vals;
        double m0 = v[k0], m1 = m0, m2 = m0, m3 = m0;
        int nan = 0;
        int64_t k = k0;
        if (c->is_max) {
            for (; k + 4 <= k1; k += 4) {
                nan |= (v[k] != v[k]) | (v[k + 1] != v[k + 1]) | (v[k + 2] != v[k + 2]) | (v[k + 3] != v[k + 3]);
                m0 = v[k] > m0 ? v[k] : m0;         m1 = v[k + 1] > m1 ? v[k + 1] : m1;
                m2 = v[k + 2] > m2 ? v[k + 2] : m2; m3 = v[k + 3] > m3 ? v[k + 3] : m3;
            }
            for (; k < k1; k++) { nan |= v[k] != v[k]; m0 = v[k] > m0 ? v[k] : m0; }
            m0 = m1 > m0 ? m1 : m0; m2 = m3 > m2 ? m3 : m2; m0 = m2 > m0 ? m2 : m0;
        } else {
            for (; k + 4 <= k1; k += 4) {
                nan |= (v[k] != v[k]) | (v[k + 1] != v[k + 1]) | (v[k + 2] != v[k + 2]) | (v[k + 3] != v[k + 3]);
                m0 = v[k] < m0 ? v[k] : m0;         m1 = v[k + 1] < m1 ? v[k + 1] : m1;
                m2 = v[k + 2] < m2 ? v[k + 2] : m2; m3 = v[k + 3] < m3 ? v[k + 3] : m3;
            }
            for (; k < k1; k++) { nan |= v[k] != v[k]; m0 = v[k] < m0 ? v[k] : m0; }
            m0 = m1 < m0 ? m1 : m0; m2 = m3 < m2 ? m3 : m2; m0 = m2 < m0 ? m2 : m0;
        }
        c->fm[i] = m0; c->nan[i] = nan || m0 != m0;
    } else {
        const int64_t* v = c->sa->vals;
        int64_t m = v[k0];
        if (c->is_max) for (int64_t k = k0; k < k1; k++) m = v[k] > m ? v[k] : m;
        else for (int64_t k = k0; k < k1; k++) m = v[k] < m ? v[k] : m;
        c->im[i] = m; c->nan[i] = false;
    }
    return true;
}

/* Append the candidates an array contributes to Max (is_max) or Min: its
 * extreme explicit value (all of them for a boxed array) and, when some
 * position is implicit, the background. False to decline. */
static bool extreme_parts(const SparseArray* sa, bool is_max, Expr*** items, size_t* n, size_t* cap) {
    int64_t total = sp_total(sa->rank, sa->dims);
    if (total == 0) return false;
    if (head_is(sa->bg, SYM_Complex) || sa->bg->type == EXPR_SPARSE || head_is(sa->bg, SYM_List)) return false;
    int64_t nnz = sa->pat->nnz;
    size_t need = *n + 2 + (sa->kind == SPV_EXPR ? (size_t)nnz : 0);
    if (need > *cap) {
        *cap = need * 2;
        *items = realloc(*items, sizeof(Expr*) * *cap);
    }
    if (nnz > 0) {
        if (sa->kind == SPV_F64 || sa->kind == SPV_INT64) {
            ExtCtx c;
            c.sa = sa; c.is_max = is_max; c.T = n_parts(nnz);
            if (c.T > nnz) c.T = 1;
            run_parts(c.T, ext_part, &c);
            if (sa->kind == SPV_F64) {
                double m = c.fm[0];
                for (int i = 0; i < c.T; i++) {
                    if (c.nan[i]) return false;   /* NaN: let the dense path decide */
                    if (is_max ? c.fm[i] > m : c.fm[i] < m) m = c.fm[i];
                }
                (*items)[(*n)++] = expr_new_real(m);
            } else {
                int64_t m = c.im[0];
                for (int i = 0; i < c.T; i++) if (is_max ? c.im[i] > m : c.im[i] < m) m = c.im[i];
                (*items)[(*n)++] = expr_new_integer(m);
            }
        } else if (sa->kind == SPV_EXPR) {
            Expr* const* v = sa->vals;
            /* a symbolic value stays a candidate (Max[x, 0] is Max[0, x], as
             * on the dense array); a complex number or a nested array does not */
            for (int64_t k = 0; k < nnz; k++)
                if (head_is(v[k], SYM_Complex) || head_is(v[k], SYM_List) ||
                    v[k]->type == EXPR_NDARRAY || v[k]->type == EXPR_SPARSE) return false;
            for (int64_t k = 0; k < nnz; k++) (*items)[(*n)++] = expr_copy(v[k]);
        } else {
            return false;
        }
    }
    if (total < 0 || total > nnz) (*items)[(*n)++] = expr_copy(sa->bg);
    return true;
}

/* Max[args] / Min[args] with every SparseArray argument replaced by its
 * candidates; the result is returned unevaluated for Max / Min to finish. */
static Expr* extreme_call(const Expr* call, bool is_max) {
    size_t na = nargs(call), n = 0, cap = na + 8;
    Expr** items = malloc(sizeof(Expr*) * cap);
    bool ok = true;
    for (size_t i = 0; i < na && ok; i++) {
        const Expr* a = arg(call, i);
        if (is_sparse(a)) ok = extreme_parts(sp_of(a), is_max, &items, &n, &cap);
        else {
            if (n + 1 > cap) { cap = 2 * cap + 1; items = realloc(items, sizeof(Expr*) * cap); }
            items[n++] = expr_copy((Expr*)a);
        }
    }
    if (!ok) {
        for (size_t i = 0; i < n; i++) expr_free(items[i]);
        free(items);
        return NULL;
    }
    Expr* out = mk_call(is_max ? SYM_Max : SYM_Min, items, n);
    free(items);
    return out;
}

static Expr* native_max(const Expr* call) { return extreme_call(call, true); }
static Expr* native_min(const Expr* call) { return extreme_call(call, false); }

static Expr* native_minmax(const Expr* call) {
    if (nargs(call) != 1 || !is_sparse(arg(call, 0))) return NULL;
    Expr* mn = extreme_call(call, false);
    if (!mn) return NULL;
    Expr* mx = extreme_call(call, true);
    if (!mx) { expr_free(mn); return NULL; }
    return mk_call2(SYM_List, mn, mx);
}

/* ======================================================================== */
/* Tr                                                                       */
/* ======================================================================== */

/* The entry index of (i, i, ..., i) in row i, or -1. */
static int64_t diag_entry(const SparseArray* sa, int64_t i) {
    const SparsePattern* p = sa->pat;
    int64_t a = p->rowptr[i], b = p->rowptr[i + 1];
    int w = p->width;
    while (a < b) {
        int64_t m = a + (b - a) / 2;
        int cmp = 0;
        for (int j = 0; j < w && !cmp; j++) {
            int64_t c = sp_col(p, m, j);
            cmp = c < i ? -1 : c > i ? 1 : 0;
        }
        if (cmp == 0) return m;
        if (cmp < 0) a = m + 1; else b = m;
    }
    return -1;
}

static Expr* native_tr(const Expr* call) {
    size_t na = nargs(call);
    if (na < 1 || na > 2 || !is_sparse(arg(call, 0))) return NULL;
    const SparseArray* sa = sp_of(arg(call, 0));
    if (sa->rank == 1) {
        if (na == 2) return NULL;
        Plan t;
        if (!plan_make(sa, 0, 0, &t)) return NULL;
        return reduce_plan(&t, 1);
    }
    int64_t m = sa->dims[0];
    for (int d = 1; d < sa->rank; d++) if (sa->dims[d] < m) m = sa->dims[d];
    Bg b;
    if (na == 1 && typed_kind(sa->kind) && bg_typed(sa->bg, &b)) {
        /* one cell of m diagonal elements */
        Acc a;
        memset(&a, 0, sizeof a);
        a.kind = sa->kind;
        int64_t isum = 0, cnt = 0;
        double fsum = 0;
        bool ok = true;
        for (int64_t i = 0; i < m && ok; i++) {
            int64_t k = diag_entry(sa, i);
            if (k < 0) continue;
            cnt++;
            if (sa->kind == SPV_F64) fsum += ((const double*)sa->vals)[k];
            else ok = !ci_add_i64(isum, ((const int64_t*)sa->vals)[k], &isum);
        }
        if (ok) {
            a.isum = &isum; a.fsum = &fsum; a.cnt = &cnt;
            Cell c;
            if (cell_finish(&a, 0, m, &b, &c)) return cell_expr(&c);
        }
    }
    if (m > SR_DENSE_MAX) return NULL;
    /* f[d1, ..., dm] (Plus by default), for the evaluator to finish */
    Expr** it = malloc(sizeof(Expr*) * (size_t)(m ? m : 1));
    for (int64_t i = 0; i < m; i++) {
        int64_t k = diag_entry(sa, i);
        it[i] = k >= 0 ? sp_value_expr(sa, k) : expr_copy(sa->bg);
    }
    Expr* f = na == 2 ? expr_copy((Expr*)arg(call, 1)) : expr_new_symbol(SYM_Plus);
    Expr* out = expr_new_function(f, it, (size_t)m);
    free(it);
    return out;
}

/* ======================================================================== */
/* Norm                                                                     */
/* ======================================================================== */

/* Machine 2-norm of a vector given as explicit values plus k copies of bg,
 * scaled against overflow / underflow. */
static void sumsq_chunk(void* vctx, size_t lo, size_t hi, double* slot) {
    const SparseArray* sa = vctx;
    double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
    size_t k = lo;
    if (sa->kind == SPV_F64) {
        const double* v = (const double*)sa->vals + lo;
        s0 = vdot(v, v, (int64_t)(hi - lo));
    } else {
        const int64_t* v = sa->vals;
        for (; k < hi; k++) { double x = (double)v[k]; s0 += x * x; }
    }
    slot[0] = (s0 + s1) + (s2 + s3);
}

typedef struct { const SparseArray* sa; int T; double part[NDARRAY_MAX_THREADS]; } SqCtx;

static bool sumsq_part(void* vctx, int i) {
    SqCtx* c = vctx;
    int64_t nnz = c->sa->pat->nnz;
    sumsq_chunk((void*)c->sa, (size_t)(nnz * i / c->T), (size_t)(nnz * (i + 1) / c->T), &c->part[i]);
    return true;
}

static double norm2_machine(const SparseArray* sa, double bgv, int64_t k) {
    int64_t nnz = sa->pat->nnz;
    double s = 0, scale = fabs(bgv) * (k > 0);
    if (nnz > 0) {
        /* memory-bound: a thread per 2^18 values */
        SqCtx c;
        c.sa = sa;
        c.T = (int)(nnz >> 18);
        if (c.T < 1) c.T = 1;
        if (c.T > NDARRAY_MAX_THREADS) c.T = NDARRAY_MAX_THREADS;
        run_parts(c.T, sumsq_part, &c);
        for (int i = 0; i < c.T; i++) s += c.part[i];
    }
    s += bgv * bgv * (double)k;
    if (isfinite(s) && s > 1e-290) return sqrt(s);
    for (int64_t i = 0; i < nnz; i++) {
        double x = sa->kind == SPV_F64 ? ((const double*)sa->vals)[i] : (double)((const int64_t*)sa->vals)[i];
        if (fabs(x) > scale) scale = fabs(x);
    }
    if (scale == 0 || !isfinite(scale)) return sqrt(s);
    s = 0;
    for (int64_t i = 0; i < nnz; i++) {
        double x = (sa->kind == SPV_F64 ? ((const double*)sa->vals)[i] : (double)((const int64_t*)sa->vals)[i]) / scale;
        s += x * x;
    }
    double xb = bgv / scale;
    s += xb * xb * (double)k;
    return scale * sqrt(s);
}

/* Integer p >= 3 vector norm over exact Integers: (Sum |x|^p)^(1/p). */
static Expr* normp_exact(const SparseArray* sa, const Bg* b, int64_t k, unsigned long p) {
    mpz_t s, t, x;
    mpz_inits(s, t, x, NULL);
    const int64_t* v = sa->vals;
    for (int64_t i = 0; i < sa->pat->nnz; i++) {
        z_set_i64(x, v[i]);
        mpz_abs(x, x);
        mpz_pow_ui(t, x, p);
        mpz_add(s, s, t);
    }
    if (k > 0) {
        z_set_i64(x, b->i);
        mpz_abs(x, x);
        mpz_pow_ui(t, x, p);
        z_set_i64(x, k);
        mpz_mul(t, t, x);
        mpz_add(s, s, t);
    }
    mpz_set_ui(t, 1);
    Expr* sum = make_rational_mpz(s, t);
    mpz_clears(s, t, x, NULL);
    return eval_free(mk_call2(SYM_Power, sum, make_rational(1, (int64_t)p)));
}

/* Single-cell transformed sum (VT_ABS / VT_SQR) of every element, or NULL. */
static bool total_cell(const SparseArray* sa, int vt, const Bg* b, Cell* out) {
    Plan t;
    if (!plan_make(sa, 0, sa->rank - 1, &t)) return false;
    Bg tb;
    if (!bg_transform(vt, b, &tb)) return false;
    Acc a;
    if (!accumulate(&t, vt, tb.real || tb.i != 0, &a)) return false;
    bool ok = cell_finish(&a, 0, t.cellsize, &tb, out);
    acc_free(&a);
    return ok;
}

static Expr* sqrt_cell(const Cell* c) {
    if (c->m) return expr_new_real(sqrt(c->d));
    return eval_free(mk_call1(SYM_Sqrt, expr_new_integer(c->i)));
}

static double lanczos_norm2(const SparseArray* sa, bool* ok);

static Expr* native_norm(const Expr* call) {
    size_t na = nargs(call);
    if (na < 1 || na > 2 || !is_sparse(arg(call, 0))) return NULL;
    const SparseArray* sa = sp_of(arg(call, 0));
    Bg b;
    if (!typed_kind(sa->kind) || !bg_typed(sa->bg, &b)) return NULL;
    const Expr* p = na == 2 ? arg(call, 1) : NULL;
    bool frob = p && p->type == EXPR_STRING && strcmp(p->data.string, "Frobenius") == 0;
    bool inf = p && is_infinity(p);
    int64_t total = sp_total(sa->rank, sa->dims);
    if (total < 0) return NULL;
    int64_t k = total - sa->pat->nnz;

    if (sa->rank == 1) {
        if (!p || frob || (p->type == EXPR_INTEGER && p->data.integer == 2)) {
            bool exact = sa->kind == SPV_INT64 && (!b.real || k == 0);
            if (!exact) return expr_new_real(norm2_machine(sa, b.d, k));
            Cell c;
            if (!total_cell(sa, VT_SQR, &b, &c)) return NULL;
            return sqrt_cell(&c);
        }
        if (inf) {
            /* Max of the absolute values */
            Bg tb;
            if (!bg_transform(VT_ABS, &b, &tb)) return NULL;
            int64_t nnz = sa->pat->nnz;
            Expr* parts[2];
            size_t q = 0;
            if (nnz > 0) {
                if (sa->kind == SPV_F64) {
                    const double* v = sa->vals;
                    double m = 0;
                    for (int64_t i = 0; i < nnz; i++) { if (v[i] != v[i]) return NULL; if (fabs(v[i]) > m) m = fabs(v[i]); }
                    parts[q++] = expr_new_real(m);
                } else {
                    const int64_t* v = sa->vals;
                    int64_t m = 0;
                    for (int64_t i = 0; i < nnz; i++) { if (v[i] == INT64_MIN) return NULL; int64_t a = v[i] < 0 ? -v[i] : v[i]; if (a > m) m = a; }
                    parts[q++] = expr_new_integer(m);
                }
            }
            if (k > 0) parts[q++] = tb.real ? expr_new_real(tb.d) : expr_new_integer(tb.i);
            if (q == 0) return NULL;
            return mk_call(SYM_Max, parts, q);
        }
        if (p->type == EXPR_INTEGER && p->data.integer == 1) {
            Cell c;
            if (!total_cell(sa, VT_ABS, &b, &c)) return NULL;
            return cell_expr(&c);
        }
        if (p->type == EXPR_INTEGER && p->data.integer > 2 && p->data.integer < 4096) {
            bool exact = sa->kind == SPV_INT64 && (!b.real || k == 0);
            if (exact) return normp_exact(sa, &b, k, (unsigned long)p->data.integer);
        }
        if ((p->type == EXPR_INTEGER && p->data.integer >= 1) || (p->type == EXPR_REAL && p->data.real >= 1)) {
            double pp = p->type == EXPR_INTEGER ? (double)p->data.integer : p->data.real;
            if (!isfinite(pp)) return NULL;
            int64_t nnz = sa->pat->nnz;
            double s = 0;
            for (int64_t i = 0; i < nnz; i++) {
                double x = sa->kind == SPV_F64 ? ((const double*)sa->vals)[i] : (double)((const int64_t*)sa->vals)[i];
                s += pow(fabs(x), pp);
            }
            if (k > 0) s += pow(fabs(b.d), pp) * (double)k;
            if (!isfinite(s)) return NULL;
            return expr_new_real(pow(s, 1.0 / pp));
        }
        return NULL;
    }
    if (sa->rank != 2) return NULL;
    if (frob) {
        bool exact = sa->kind == SPV_INT64 && (!b.real || k == 0);
        if (!exact) return expr_new_real(norm2_machine(sa, b.d, k));
        Cell c;
        if (!total_cell(sa, VT_SQR, &b, &c)) return NULL;
        return sqrt_cell(&c);
    }
    if (p && (inf || (p->type == EXPR_INTEGER && p->data.integer == 1))) {
        /* Max over the absolute column sums (1) or row sums (Infinity) */
        Plan t;
        if (!plan_make(sa, inf ? 1 : 0, inf ? 1 : 0, &t) || t.ncells > SR_DENSE_MAX || t.ncells == 0) return NULL;
        Bg tb;
        if (!bg_transform(VT_ABS, &b, &tb)) return NULL;
        Acc a;
        if (!accumulate(&t, VT_ABS, tb.real || tb.i != 0, &a)) return NULL;
        bool have_m = false, have_x = false, ok = true;
        double bm = 0;
        int64_t bx = 0;
        for (int64_t c = 0; c < t.ncells && ok; c++) {
            Cell cl;
            ok = cell_finish(&a, c, t.cellsize, &tb, &cl);
            if (!ok) break;
            if (cl.m) { if (!have_m || cl.d > bm) bm = cl.d; have_m = true; }
            else { if (!have_x || cl.i > bx) bx = cl.i; have_x = true; }
        }
        acc_free(&a);
        if (!ok) return NULL;
        if (have_m && have_x) return mk_call2(SYM_Max, expr_new_real(bm), expr_new_integer(bx));
        return have_m ? expr_new_real(bm) : expr_new_integer(bx);
    }
    if (!p || (p->type == EXPR_INTEGER && p->data.integer == 2)) {
        /* the largest singular value: machine values over a zero background */
        if (sa->kind != SPV_F64) return NULL;
        if (!(is_exact_zero(sa->bg) || (sa->bg->type == EXPR_REAL && sa->bg->data.real == 0.0))) return NULL;
        bool ok = false;
        double v = lanczos_norm2(sa, &ok);
        return ok ? expr_new_real(v) : NULL;
    }
    return NULL;
}

/* ---- the matrix 2-norm: Lanczos on A^T A --------------------------------- */

/* A CSR matrix of doubles with 64-bit or 32-bit column indices. */
typedef struct {
    int64_t nr, nc;
    const int64_t* rp;
    const void* ci;
    bool idx32;
    const double* v;
} Csr;

typedef struct {
    const Csr* A;
    const double* x;
    double* y;
} SpmvCtx;

static bool spmv_chunk(void* vctx, size_t lo, size_t hi) {
    SpmvCtx* c = vctx;
    const Csr* A = c->A;
    int64_t nnz = A->rp[A->nr];
    /* rows whose first entry lies in [lo, hi) */
    int64_t a = 0, b = A->nr;
    while (a < b) { int64_t m = a + (b - a) / 2; if (A->rp[m] < (int64_t)lo) a = m + 1; else b = m; }
    int64_t r0 = lo == 0 ? 0 : a, r1;
    if ((int64_t)hi >= nnz) r1 = A->nr;
    else {
        a = 0; b = A->nr;
        while (a < b) { int64_t m = a + (b - a) / 2; if (A->rp[m] < (int64_t)hi) a = m + 1; else b = m; }
        r1 = a;
    }
    const double* x = c->x;
    for (int64_t r = r0; r < r1; r++) {
        double s = 0;
        if (A->idx32) {
            const int32_t* ci = A->ci;
            for (int64_t k = A->rp[r]; k < A->rp[r + 1]; k++) s += A->v[k] * x[ci[k]];
        } else {
            const int64_t* ci = A->ci;
            for (int64_t k = A->rp[r]; k < A->rp[r + 1]; k++) s += A->v[k] * x[ci[k]];
        }
        c->y[r] = s;
    }
    return true;
}

static void spmv(const Csr* A, const double* x, double* y) {
    SpmvCtx c = { A, x, y };
    int64_t nnz = A->rp[A->nr];
    if (nnz == 0) { memset(y, 0, sizeof(double) * (size_t)A->nr); return; }
    nd_parallel_for((size_t)nnz, spmv_chunk, &c);
}

/* The transpose of A as CSR (that is, A in CSC), by a counting sort. */
static bool csr_transpose(const Csr* A, Csr* T, int64_t** rp, int64_t** ci, double** v) {
    int64_t nnz = A->rp[A->nr];
    *rp = calloc((size_t)A->nc + 1, sizeof(int64_t));
    *ci = malloc(sizeof(int64_t) * (size_t)(nnz ? nnz : 1));
    *v = malloc(sizeof(double) * (size_t)(nnz ? nnz : 1));
    if (!*rp || !*ci || !*v) { free(*rp); free(*ci); free(*v); return false; }
    int64_t* p = *rp;
    for (int64_t k = 0; k < nnz; k++) {
        int64_t c = A->idx32 ? ((const int32_t*)A->ci)[k] : ((const int64_t*)A->ci)[k];
        p[c + 1]++;
    }
    for (int64_t c = 0; c < A->nc; c++) p[c + 1] += p[c];
    int64_t* pos = malloc(sizeof(int64_t) * (size_t)(A->nc ? A->nc : 1));
    if (!pos) { free(*rp); free(*ci); free(*v); return false; }
    memcpy(pos, p, sizeof(int64_t) * (size_t)A->nc);
    for (int64_t r = 0; r < A->nr; r++)
        for (int64_t k = A->rp[r]; k < A->rp[r + 1]; k++) {
            int64_t c = A->idx32 ? ((const int32_t*)A->ci)[k] : ((const int64_t*)A->ci)[k];
            int64_t d = pos[c]++;
            (*ci)[d] = r;
            (*v)[d] = A->v[k];
        }
    free(pos);
    T->nr = A->nc; T->nc = A->nr; T->rp = p; T->ci = *ci; T->idx32 = false; T->v = *v;
    return true;
}

/* Number of eigenvalues of the symmetric tridiagonal (a[0..n), b[0..n-1))
 * strictly less than x (Sturm count by the LDL^T pivots). */
static int64_t sturm_count(const double* a, const double* b, int64_t n, double x) {
    int64_t cnt = 0;
    double d = 1.0;
    for (int64_t i = 0; i < n; i++) {
        double bb = i ? b[i - 1] * b[i - 1] : 0.0;
        d = (a[i] - x) - (i ? bb / d : 0.0);
        if (d == 0.0) d = -DBL_EPSILON * (fabs(x) + DBL_MIN);
        if (d < 0) cnt++;
    }
    return cnt;
}

/* The j-th largest eigenvalue (j = 1, 2) of the tridiagonal, by bisection to
 * full precision within [lo, hi]. */
static double tri_eig_from_top(const double* a, const double* b, int64_t n, int64_t j, double lo, double hi) {
    for (int it = 0; it < 200; it++) {
        double mid = 0.5 * (lo + hi);
        if (mid <= lo || mid >= hi) break;
        if (sturm_count(a, b, n, mid) <= n - j) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

/* Last component of the unit eigenvector of the tridiagonal for its largest
 * eigenvalue th (inverse iteration with a shift just above th, where T - sI
 * is negative definite, so the LDL^T solve needs no pivoting). */
static double tri_last_component(const double* a, const double* b, int64_t n, double th, double* work) {
    double s = th + fmax(fabs(th), 1e-300) * 64 * DBL_EPSILON;
    double* d = work;        /* n pivots */
    double* x = work + n;    /* n */
    double* l = work + 2 * n;
    for (int64_t i = 0; i < n; i++) {
        d[i] = (a[i] - s) - (i ? l[i - 1] * b[i - 1] : 0.0);
        if (d[i] == 0.0) d[i] = -1e-300;
        if (i + 1 < n) l[i] = b[i] / d[i];
    }
    for (int64_t i = 0; i < n; i++) x[i] = 1.0;
    for (int pass = 0; pass < 3; pass++) {
        for (int64_t i = 1; i < n; i++) x[i] -= l[i - 1] * x[i - 1];      /* L y = x */
        for (int64_t i = 0; i < n; i++) x[i] /= d[i];                     /* D z = y */
        for (int64_t i = n - 2; i >= 0; i--) x[i] -= l[i] * x[i + 1];     /* L^T w = z */
        double nn = 0;
        for (int64_t i = 0; i < n; i++) nn += x[i] * x[i];
        nn = sqrt(nn);
        if (!(nn > 0) || !isfinite(nn)) return 1.0;
        for (int64_t i = 0; i < n; i++) x[i] /= nn;
    }
    return x[n - 1];
}

/* sigma_max(A) by the Lanczos iteration on M = A^T A (or A A^T, whichever is
 * smaller), without reorthogonalisation: the extreme Ritz value converges to
 * the extreme eigenvalue regardless (loss of orthogonality only duplicates
 * converged values). Stops on a residual bound at machine precision. */
static double lanczos_norm2(const SparseArray* sa, bool* ok) {
    *ok = false;
    const SparsePattern* p = sa->pat;
    Csr A = { sa->dims[0], sa->dims[1], p->rowptr, p->colidx, p->idx32, sa->vals };
    if (p->nnz == 0) { *ok = true; return 0.0; }
    Csr T;
    int64_t *trp, *tci;
    double* tv;
    if (!csr_transpose(&A, &T, &trp, &tci, &tv)) return 0.0;
    /* M v = B^T (B v) with B the factor whose columns are the smaller side */
    const Csr* B  = A.nc <= A.nr ? &A : &T;   /* B: m x n, n = smaller */
    const Csr* Bt = A.nc <= A.nr ? &T : &A;
    int64_t n = B->nc, m = B->nr;
    int64_t maxit = n + 2 < 4000 ? n + 2 : 4000;
    double* v  = malloc(sizeof(double) * (size_t)n);
    double* vp = calloc((size_t)n, sizeof(double));
    double* w  = malloc(sizeof(double) * (size_t)n);
    double* u  = malloc(sizeof(double) * (size_t)m);
    double* al = malloc(sizeof(double) * (size_t)maxit);
    double* be = malloc(sizeof(double) * (size_t)maxit);
    double* work = malloc(sizeof(double) * 3 * (size_t)maxit);
    double result = 0.0;
    if (!v || !vp || !w || !u || !al || !be || !work) goto done;
    {
        /* a fixed, generic start vector (deterministic answers) */
        uint64_t st = UINT64_C(0x9E3779B97F4A7C15);
        double nn = 0;
        for (int64_t i = 0; i < n; i++) {
            st = st * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
            v[i] = 0.5 + (double)(st >> 11) * (1.0 / 9007199254740992.0);
            nn += v[i] * v[i];
        }
        nn = sqrt(nn);
        for (int64_t i = 0; i < n; i++) v[i] /= nn;
        double beta = 0, theta = 0, prev = -1;
        int stable = 0;
        int64_t j, next_check = 65;
        for (j = 0; j < maxit; j++) {
            spmv(B, v, u);
            spmv(Bt, u, w);
            double alpha = 0;
            alpha = vdot(v, w, n);
            for (int64_t i = 0; i < n; i++) w[i] -= alpha * v[i] + beta * vp[i];
            /* one local reorthogonalisation against v */
            double c = 0;
            c = vdot(v, w, n);
            for (int64_t i = 0; i < n; i++) w[i] -= c * v[i];
            alpha += c;
            double bn = 0;
            bn = vdot(w, w, n);
            bn = sqrt(bn);
            al[j] = alpha;
            be[j] = bn;
            int64_t sz = j + 1;
            /* Ritz checks cost O(sz) each: every step early, then spaced
             * geometrically so they stay a small fraction of the work */
            bool check = sz <= 64 || sz >= next_check || bn == 0;
            if (check && sz > 64) next_check = sz + sz / 8;
            if (check || bn == 0) {
                /* Gershgorin bounds of T_sz */
                double lo = al[0], hi = al[0];
                for (int64_t i = 0; i < sz; i++) {
                    double r = (i ? fabs(be[i - 1]) : 0) + (i + 1 < sz ? fabs(be[i]) : 0);
                    if (al[i] - r < lo) lo = al[i] - r;
                    if (al[i] + r > hi) hi = al[i] + r;
                }
                theta = tri_eig_from_top(al, be, sz, 1, lo, hi);
                double s = sz > 1 ? tri_last_component(al, be, sz, theta, work) : 1.0;
                double res = bn * fabs(s);
                double tol = 4 * DBL_EPSILON * fabs(theta);
                bool conv = res <= tol || bn <= tol;
                if (!conv && sz > 1) {
                    double th2 = tri_eig_from_top(al, be, sz, 2, lo, theta);
                    double gap = theta - th2;
                    if (gap > 1e-10 * fabs(theta) && res * res / gap <= DBL_EPSILON * fabs(theta)) conv = true;
                }
                if (prev >= 0 && fabs(theta - prev) <= DBL_EPSILON * fabs(theta)) stable++;
                else stable = 0;
                prev = theta;
                if (conv || stable >= 3) break;
            }
            if (bn == 0) break;
            for (int64_t i = 0; i < n; i++) { vp[i] = v[i]; v[i] = w[i] / bn; }
            beta = bn;
        }
        if (j < maxit && theta >= 0 && isfinite(theta)) { result = sqrt(theta); *ok = true; }
    }
done:
    free(v); free(vp); free(w); free(u); free(al); free(be); free(work);
    free(trp); free(tci); free(tv);
    return result;
}

/* ======================================================================== */
/* Equal / Unequal                                                          */
/* ======================================================================== */

enum { EQ_FALSE = 0, EQ_TRUE = 1, EQ_UNKNOWN = 2 };

/* Equal on two machine numbers, as the scalar Equal decides it: exact on two
 * Integers, otherwise within Internal`$EqualTolerance (equal_tol.h). */
static inline int eq_num(bool am, int64_t ai, double ad, bool bm, int64_t bi, double bd) {
    if (!am && !bm) return ai == bi ? EQ_TRUE : EQ_FALSE;
    double x = am ? ad : (double)ai, y = bm ? bd : (double)bi;
    if (x != x || y != y) return EQ_FALSE;
    if (x == y) return EQ_TRUE;
    return mth_machine_tolerant_cmp(x, y, mth_machine_equal_rel()) == 0 ? EQ_TRUE : EQ_FALSE;
}

/* Equal on two arbitrary scalars, by the Equal builtin. */
static int eq_expr(const Expr* a, const Expr* b) {
    if (expr_eq((Expr*)a, (Expr*)b)) {
        if (a->type == EXPR_REAL && a->data.real != a->data.real) return EQ_FALSE;
        if (a->type != EXPR_SYMBOL || a->data.symbol.name != SYM_Indeterminate) return EQ_TRUE;
    }
    Expr* r = eval_free(mk_call2(SYM_Equal, expr_copy((Expr*)a), expr_copy((Expr*)b)));
    int v = EQ_UNKNOWN;
    if (r->type == EXPR_SYMBOL) {
        if (r->data.symbol.name == SYM_True) v = EQ_TRUE;
        else if (r->data.symbol.name == SYM_False) v = EQ_FALSE;
    }
    expr_free(r);
    return v;
}

/* One element source: explicit value k of a sparse array, a background, or
 * an element of a dense buffer, read as a machine number when it is one. */
typedef struct {
    bool num;       /* machine Integer / Real */
    bool m;         /* Real */
    int64_t i;
    double d;
} Num;

static bool expr_num(const Expr* e, Num* n) {
    if (e->type == EXPR_INTEGER) { n->num = true; n->m = false; n->i = e->data.integer; n->d = 0; return true; }
    if (e->type == EXPR_REAL) { n->num = true; n->m = true; n->d = e->data.real; n->i = 0; return true; }
    n->num = false;
    return false;
}

static inline bool sa_num(const SparseArray* sa, int64_t k, Num* n) {
    if (sa->kind == SPV_INT64) { n->num = true; n->m = false; n->i = ((const int64_t*)sa->vals)[k]; return true; }
    if (sa->kind == SPV_F64) { n->num = true; n->m = true; n->d = ((const double*)sa->vals)[k]; return true; }
    if (sa->kind == SPV_EXPR) return expr_num(((Expr* const*)sa->vals)[k], n);
    n->num = false;
    return false;
}

static Expr* sa_elem_expr(const SparseArray* sa, int64_t k) {
    return k >= 0 ? sp_value_expr(sa, k) : expr_copy(sa->bg);
}

/* Compare two elements, each an explicit entry (k >= 0) or the background. */
static int eq_elems(const SparseArray* a, int64_t ka, const SparseArray* b, int64_t kb) {
    Num x = { false, false, 0, 0.0 }, y = x;
    bool nx = ka >= 0 ? sa_num(a, ka, &x) : expr_num(a->bg, &x);
    bool ny = kb >= 0 ? sa_num(b, kb, &y) : expr_num(b->bg, &y);
    if (nx && ny) return eq_num(x.m, x.i, x.d, y.m, y.i, y.d);
    Expr* ea = sa_elem_expr(a, ka);
    Expr* eb = sa_elem_expr(b, kb);
    int v = eq_expr(ea, eb);
    expr_free(ea);
    expr_free(eb);
    return v;
}

/* Lexicographic comparison of the tuples of entry ka of a and kb of b. */
static inline int tuple_cmp2(const SparsePattern* pa, int64_t ka, const SparsePattern* pb, int64_t kb) {
    for (int j = 0; j < pa->width; j++) {
        int64_t x = sp_col(pa, ka, j), y = sp_col(pb, kb, j);
        if (x != y) return x < y ? -1 : 1;
    }
    return 0;
}

/* Two SparseArrays of equal dims: merge the explicit positions row by row. */
static int eq_sparse_sparse(const SparseArray* a, const SparseArray* b) {
    for (int d = 0; d < a->rank; d++) if (a->dims[d] != b->dims[d]) return EQ_FALSE;
    int64_t total = sp_total(a->rank, a->dims);
    const SparsePattern* pa = a->pat;
    const SparsePattern* pb = b->pat;
    int verdict = EQ_TRUE;
    /* Same pattern and both machine: a straight buffer comparison. */
    if (sp_pattern_same(pa, pb) && typed_kind(a->kind) && typed_kind(b->kind)) {
        for (int64_t k = 0; k < pa->nnz; k++) {
            Num x = { false, false, 0, 0.0 }, y = x;
            sa_num(a, k, &x); sa_num(b, k, &y);
            if (eq_num(x.m, x.i, x.d, y.m, y.i, y.d) == EQ_FALSE) return EQ_FALSE;
        }
        if (total < 0 || total > pa->nnz) return eq_elems(a, -1, b, -1);
        return EQ_TRUE;
    }
    int64_t uni = 0;
    for (int64_t r = 0; r < pa->nrows; r++) {
        int64_t i = pa->rowptr[r], ie = pa->rowptr[r + 1];
        int64_t j = pb->rowptr[r], je = pb->rowptr[r + 1];
        while (i < ie || j < je) {
            int c = i >= ie ? 1 : j >= je ? -1 : tuple_cmp2(pa, i, pb, j);
            int v;
            if (c == 0) v = eq_elems(a, i++, b, j++);
            else if (c < 0) v = eq_elems(a, i++, b, -1);
            else v = eq_elems(a, -1, b, j++);
            uni++;
            if (v == EQ_FALSE) return EQ_FALSE;
            if (v == EQ_UNKNOWN) verdict = EQ_UNKNOWN;
        }
    }
    if (total < 0 || total > uni) {
        int v = eq_elems(a, -1, b, -1);
        if (v == EQ_FALSE) return EQ_FALSE;
        if (v == EQ_UNKNOWN) verdict = EQ_UNKNOWN;
    }
    return verdict;
}

/* A dense operand: a packed buffer (Integer / Real) or a full nested List,
 * flattened into borrowed element pointers. */
typedef struct {
    int rank;
    int64_t dims[SP_MAX_RANK];
    bool packed;
    NDType dt;
    const void* buf;
    const Expr** items;     /* boxed: borrowed, row-major */
    int64_t n;
} Dense;

/* Is `e` a full nested List of the dims in d (its leaves not themselves
 * arrays)? */
static bool list_check(const Expr* e, int depth, const Dense* d) {
    if (depth == d->rank)
        return !head_is(e, SYM_List) && e->type != EXPR_NDARRAY && e->type != EXPR_SPARSE;
    if (!head_is(e, SYM_List) || (int64_t)e->data.function.arg_count != d->dims[depth]) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!list_check(e->data.function.args[i], depth + 1, d)) return false;
    return true;
}

static void list_flatten(const Expr* e, int depth, Dense* d) {
    if (depth == d->rank) { d->items[d->n++] = e; return; }
    for (size_t i = 0; i < e->data.function.arg_count; i++) list_flatten(e->data.function.args[i], depth + 1, d);
}

/* False when `e` is not a full array (or not one handled here). */
static bool dense_open(const Expr* e, Dense* d) {
    memset(d, 0, sizeof *d);
    if (e->type == EXPR_NDARRAY) {
        const NDArrayData* nd = &e->data.ndarray;
        if (nd->dtype != NDT_FLOAT64 && nd->dtype != NDT_INT64) return false;
        if (nd->rank > SP_MAX_RANK) return false;
        d->rank = nd->rank;
        for (int i = 0; i < nd->rank; i++) d->dims[i] = nd->dims[i];
        d->packed = true; d->dt = nd->dtype; d->buf = nd->data;
        return true;
    }
    if (!head_is(e, SYM_List)) return false;
    /* dims along the first-element spine, then a full check */
    const Expr* x = e;
    while (head_is(x, SYM_List)) {
        if (d->rank >= SP_MAX_RANK) return false;
        d->dims[d->rank++] = (int64_t)x->data.function.arg_count;
        if (x->data.function.arg_count == 0) break;
        x = x->data.function.args[0];
    }
    if (!list_check(e, 0, d)) return false;
    int64_t total = sp_total(d->rank, d->dims);
    if (total < 0 || total > SR_DENSE_MAX) return false;
    d->items = malloc(sizeof(Expr*) * (size_t)(total ? total : 1));
    if (!d->items) return false;
    list_flatten(e, 0, d);
    return true;
}

static int eq_dense_elem(const Dense* d, int64_t off, const SparseArray* sa, int64_t k) {
    Num x = { false, false, 0, 0.0 }, y = x;
    bool nx;
    if (d->packed) {
        if (d->dt == NDT_INT64) { x.m = false; x.i = ((const int64_t*)d->buf)[off]; }
        else { x.m = true; x.d = ((const double*)d->buf)[off]; }
        nx = true;
    } else {
        nx = expr_num(d->items[off], &x);
    }
    bool ny = k >= 0 ? sa_num(sa, k, &y) : expr_num(sa->bg, &y);
    if (nx && ny) return eq_num(x.m, x.i, x.d, y.m, y.i, y.d);
    Expr* ea = d->packed ? (x.m ? expr_new_real(x.d) : expr_new_integer(x.i)) : expr_copy((Expr*)d->items[off]);
    Expr* eb = sa_elem_expr(sa, k);
    int v = eq_expr(ea, eb);
    expr_free(ea);
    expr_free(eb);
    return v;
}

/* Row-major offset of entry k (row r). */
static inline int64_t entry_off(const SparseArray* sa, int64_t r, int64_t k) {
    const SparsePattern* p = sa->pat;
    int64_t off = sa->rank >= 2 ? r : 0;
    int first = sa->rank >= 2 ? 1 : 0;
    for (int j = 0; j < p->width; j++) off = off * sa->dims[first + j] + sp_col(p, k, j);
    return off;
}

static int eq_sparse_dense(const SparseArray* sa, const Dense* d) {
    if (d->rank != sa->rank) return EQ_FALSE;
    for (int i = 0; i < sa->rank; i++) if (d->dims[i] != sa->dims[i]) return EQ_FALSE;
    int64_t total = sp_total(sa->rank, sa->dims);
    if (total < 0) return EQ_UNKNOWN;
    const SparsePattern* p = sa->pat;
    int verdict = EQ_TRUE;
    int64_t off = 0;
    for (int64_t r = 0; r < p->nrows; r++)
        for (int64_t k = p->rowptr[r]; k <= p->rowptr[r + 1]; k++) {
            /* the implicit run before entry k (or the end of the row's span) */
            int64_t stop;
            if (k < p->rowptr[r + 1]) stop = entry_off(sa, r, k);
            else stop = (sa->rank >= 2) ? (r + 1) * (total / (sa->dims[0] ? sa->dims[0] : 1)) : total;
            for (; off < stop; off++) {
                int v = eq_dense_elem(d, off, sa, -1);
                if (v == EQ_FALSE) return EQ_FALSE;
                if (v == EQ_UNKNOWN) verdict = EQ_UNKNOWN;
            }
            if (k < p->rowptr[r + 1]) {
                int v = eq_dense_elem(d, off++, sa, k);
                if (v == EQ_FALSE) return EQ_FALSE;
                if (v == EQ_UNKNOWN) verdict = EQ_UNKNOWN;
            }
        }
    return verdict;
}

static bool is_array_operand(const Expr* e) {
    return is_sparse(e) || e->type == EXPR_NDARRAY || head_is(e, SYM_List);
}

/* Equal on one pair, at least one of them sparse. */
static int eq_pair(const Expr* a, const Expr* b) {
    if (is_sparse(a) && is_sparse(b)) {
        const SparseArray* x = sp_of(a);
        const SparseArray* y = sp_of(b);
        /* Unequal lengths are unequal; equal lengths but a different rank
         * compare a row with an element -- {1, 2} == {{1, 2}, {1, 2}} stays
         * unevaluated, as on the dense forms. */
        if (x->dims[0] != y->dims[0]) return EQ_FALSE;
        if (x->rank != y->rank) return EQ_UNKNOWN;
        return eq_sparse_sparse(x, y);
    }
    const Expr* s = is_sparse(a) ? a : b;
    const Expr* o = is_sparse(a) ? b : a;
    Dense d;
    if (!dense_open(o, &d)) return EQ_UNKNOWN;
    int v = eq_sparse_dense(sp_of(s), &d);
    free(d.items);
    return v;
}

static Expr* native_equal(const Expr* call) {
    size_t n = nargs(call);
    if (n < 2) return NULL;
    for (size_t i = 0; i < n; i++) if (!is_array_operand(arg(call, i))) return NULL;
    int verdict = EQ_TRUE;
    for (size_t i = 0; i + 1 < n; i++) {
        const Expr* a = arg(call, i);
        const Expr* b = arg(call, i + 1);
        int v;
        if (!is_sparse(a) && !is_sparse(b)) {
            v = eq_expr(a, b);
        } else {
            v = eq_pair(a, b);
        }
        if (v == EQ_FALSE) return expr_new_symbol(SYM_False);
        if (v == EQ_UNKNOWN) verdict = EQ_UNKNOWN;
    }
    return verdict == EQ_TRUE ? expr_new_symbol(SYM_True) : NULL;
}

static Expr* native_unequal(const Expr* call) {
    size_t n = nargs(call);
    if (n < 2) return NULL;
    for (size_t i = 0; i < n; i++) if (!is_array_operand(arg(call, i))) return NULL;
    bool all_unequal = true;
    for (size_t i = 0; i < n; i++)
        for (size_t j = i + 1; j < n; j++) {
            const Expr* a = arg(call, i);
            const Expr* b = arg(call, j);
            int v = (!is_sparse(a) && !is_sparse(b)) ? eq_expr(a, b) : eq_pair(a, b);
            if (v == EQ_TRUE) return expr_new_symbol(SYM_False);
            if (v == EQ_UNKNOWN) all_unequal = false;
        }
    return all_unequal ? expr_new_symbol(SYM_True) : NULL;
}

/* ======================================================================== */
/* Accumulate                                                               */
/* ======================================================================== */

/* Accumulate[s]: running sums along level 1, streamed from the rows into a
 * dense answer. Typed: Integer values over an Integer background (exact,
 * overflow-checked), or Real values over a Real or exact 0 background. */
static Expr* native_accumulate(const Expr* call) {
    if (nargs(call) != 1 || !is_sparse(arg(call, 0))) return NULL;
    const SparseArray* sa = sp_of(arg(call, 0));
    Bg b;
    if (!typed_kind(sa->kind) || !bg_typed(sa->bg, &b)) return NULL;
    bool f64 = sa->kind == SPV_F64;
    if (f64 ? !(b.real || b.i == 0) : b.real) return NULL;
    int64_t total = sp_total(sa->rank, sa->dims);
    if (total < 0 || total > SR_DENSE_MAX || total == 0) return NULL;
    void* buf = NULL;
    NDType dt = f64 ? NDT_FLOAT64 : NDT_INT64;
    Expr* node = ndbuild_open(sa->rank, sa->dims, dt, &buf);
    void* own = NULL;
    if (!node) {
        own = malloc((size_t)total * 8);
        if (!own) return NULL;
        buf = own;
    }
    const SparsePattern* p = sa->pat;
    bool ok = true;
    if (sa->rank == 1) {
        int64_t n = sa->dims[0], k = 0, kend = p->rowptr[1];
        if (f64) {
            double* o = buf;
            const double* v = sa->vals;
            double s = 0, bd = b.d;
            for (int64_t i = 0; i < n; i++) {
                if (k < kend && sp_col(p, k, 0) == i) s += v[k++];
                else s += bd;
                o[i] = s;
            }
        } else {
            int64_t* o = buf;
            const int64_t* v = sa->vals;
            int64_t s = 0, bi = b.i;
            for (int64_t i = 0; i < n && ok; i++) {
                int64_t x = (k < kend && sp_col(p, k, 0) == i) ? v[k++] : bi;
                if (ci_add_i64(s, x, &s)) ok = false;
                o[i] = s;
            }
        }
    } else {
        /* row i of the answer = row i-1 + row i of s (rows of `inner` elements) */
        int64_t nr = sa->dims[0], inner = total / (nr ? nr : 1);
        for (int64_t r = 0; r < nr && ok; r++) {
            int64_t base = r * inner;
            if (f64) {
                double* o = buf;
                double bd = b.d;
                if (r == 0) for (int64_t i = 0; i < inner; i++) o[i] = bd;
                else for (int64_t i = 0; i < inner; i++) o[base + i] = o[base - inner + i] + bd;
                const double* v = sa->vals;
                for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1]; k++) {
                    int64_t off = entry_off(sa, r, k);
                    o[off] = r > 0 ? o[off - inner] + v[k] : v[k];
                }
            } else {
                int64_t* o = buf;
                int64_t bi = b.i;
                if (r == 0) for (int64_t i = 0; i < inner; i++) o[i] = bi;
                else for (int64_t i = 0; i < inner && ok; i++)
                    if (ci_add_i64(o[base - inner + i], bi, &o[base + i])) ok = false;
                const int64_t* v = sa->vals;
                for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1] && ok; k++) {
                    int64_t off = entry_off(sa, r, k);
                    if (r > 0) { if (ci_add_i64(o[off - inner], v[k], &o[off])) ok = false; }
                    else o[off] = v[k];
                }
            }
        }
    }
    if (!ok) {
        if (node) expr_free(node);
        free(own);
        return NULL;
    }
    if (node) return node;
    Expr** items = malloc(sizeof(Expr*) * (size_t)total);
    for (int64_t i = 0; i < total; i++)
        items[i] = f64 ? expr_new_real(((double*)own)[i]) : expr_new_integer(((int64_t*)own)[i]);
    free(own);
    Expr* out = dense_from_items(items, sa->rank, sa->dims);
    free(items);
    return out;
}

/* ======================================================================== */
/* Registration                                                             */
/* ======================================================================== */

void sparse_reduce_init(void) {
    symtab_set_sparse_native("Total", (const void*)native_total);
    symtab_set_sparse_native("Mean", (const void*)native_mean);
    symtab_set_sparse_native("Variance", (const void*)native_variance);
    symtab_set_sparse_native("StandardDeviation", (const void*)native_sd);
    symtab_set_sparse_native("Max", (const void*)native_max);
    symtab_set_sparse_native("Min", (const void*)native_min);
    symtab_set_sparse_native("MinMax", (const void*)native_minmax);
    symtab_set_sparse_native("Tr", (const void*)native_tr);
    symtab_set_sparse_native("Norm", (const void*)native_norm);
    symtab_set_sparse_native("Equal", (const void*)native_equal);
    symtab_set_sparse_native("Unequal", (const void*)native_unequal);
    symtab_set_sparse_native("Accumulate", (const void*)native_accumulate);
}
