/*
 * linkedlist_ops.c -- native numeric algorithms over LinkedList storage.
 *
 * Every routine here walks the chunks directly through the ll_span_at cursor:
 * no List is built, no Expr is created per element, nothing is copied to a
 * contiguous buffer. A full scan costs O(n) slot reads plus O(log n) cursor
 * calls (chunk capacities double), so it runs at the speed of the packed
 * kernels it mirrors.
 *
 * Results agree with the same operation on Normal[ll] to the bit. A machine
 * float sum reproduces ndreduce.c's summation tree exactly -- the same
 * nd_parallel_reduce split into thread ranges, the same pairwise halving down
 * to ND_PAIRWISE_BLOCK, and the same nd_sum_block_f64 leaf -- only over
 * logical indices instead of buffer offsets. A leaf block that straddles two
 * chunks (or lies in the downward-growing front run) is gathered into a
 * 128-element stack buffer first; every other leaf sums the chunk in place.
 */
#include "linkedlist.h"
#include "ndreduce.h"          /* ND_PAIRWISE_BLOCK, nd_sum_block_f64, ndred_result_degenerate */
#include "ndarray_internal.h"  /* nd_parallel_reduce */
#include "checked_int.h"       /* ci_add_i64 */
#include "sym_names.h"

#include <gmp.h>
#include <string.h>

/* ------------------------------------------------------------------------
 *  Machine float sum (LL_REAL, LL_CPLX)
 * ---------------------------------------------------------------------- */

/* Elements [base, base + count) of a REAL (w = 1) or CPLX (w = 2) list as a
 * contiguous run of doubles in logical order: the slot itself when one chunk
 * holds them in ascending order, else gathered into `tmp`. */
static const double* ll_run_f64(const LinkedList* ll, int64_t base, int64_t count,
                                size_t w, double* tmp) {
    const void* p;
    int dir;
    int64_t k = ll_span_at(ll, base, &p, &dir);
    if (dir > 0 && k >= count) return (const double*)p;
    int64_t got = 0;
    while (got < count) {
        k = ll_span_at(ll, base + got, &p, &dir);
        if (k > count - got) k = count - got;
        const double* d = (const double*)p;
        if (dir > 0) {
            memcpy(tmp + (size_t)got * w, d, (size_t)k * w * sizeof(double));
        } else {
            for (int64_t t = 0; t < k; t++)
                memcpy(tmp + (size_t)(got + t) * w, d - (ptrdiff_t)((size_t)t * w),
                       w * sizeof(double));
        }
        got += k;
    }
    return tmp;
}

/* ndreduce.c's nd_sum_strided (stride 1), over logical indices. */
static void ll_sum_range(const LinkedList* ll, bool cplx, int64_t base, int64_t count,
                         double* re, double* im) {
    if ((size_t)count <= ND_PAIRWISE_BLOCK) {
        double tmp[2 * ND_PAIRWISE_BLOCK];
        size_t w = cplx ? 2 : 1;
        const double* d = ll_run_f64(ll, base, count, w, tmp);
        if (!cplx) {
            *re = nd_sum_block_f64(d, (size_t)count);
            *im = 0.0;
        } else {
            double sr = 0.0, si = 0.0;
            for (int64_t i = 0; i < count; i++) { sr += d[2 * i]; si += d[2 * i + 1]; }
            *re = sr; *im = si;
        }
        return;
    }
    int64_t half = count / 2;
    double r1, i1, r2, i2;
    ll_sum_range(ll, cplx, base, half, &r1, &i1);
    ll_sum_range(ll, cplx, base + half, count - half, &r2, &i2);
    *re = r1 + r2; *im = i1 + i2;
}

typedef struct { const LinkedList* ll; bool cplx; } LLSumCtx;

static void ll_sum_reduce(void* c, size_t lo, size_t hi, double* slot) {
    const LLSumCtx* x = (const LLSumCtx*)c;
    ll_sum_range(x->ll, x->cplx, (int64_t)lo, (int64_t)(hi - lo), &slot[0], &slot[1]);
}

/* ------------------------------------------------------------------------
 *  Exact int64 sum (LL_I64), promoting to GMP on overflow
 * ---------------------------------------------------------------------- */

static Expr* ll_sum_i64(const LinkedList* ll) {
    int64_t n = ll_length(ll), i = 0, acc = 0;
    while (i < n) {
        const void* p;
        int dir;
        int64_t k = ll_span_at(ll, i, &p, &dir);
        const int64_t* v = (const int64_t*)p;
        for (int64_t t = 0; t < k; t++) {
            int64_t x = v[dir > 0 ? t : -t], nx;
            if (!ci_add_i64(acc, x, &nx)) { acc = nx; continue; }
            {
                /* Overflow: carry on exactly in GMP from this element. */
                mpz_t z, y;
                mth_mpz_init_set_i64(z, acc);
                mpz_init(y);
                for (int64_t j = i + t; j < n; ) {
                    const void* q;
                    int dj;
                    int64_t kj = ll_span_at(ll, j, &q, &dj);
                    const int64_t* w = (const int64_t*)q;
                    for (int64_t s = 0; s < kj; s++) {
                        mth_mpz_set_i64(y, w[dj > 0 ? s : -s]);
                        mpz_add(z, z, y);
                    }
                    j += kj;
                }
                Expr* r = expr_bigint_normalize(expr_new_bigint_from_mpz(z));
                mpz_clear(z);
                mpz_clear(y);
                return r;
            }
        }
        i += k;
    }
    return expr_new_integer(acc);
}

Expr* ll_sum_machine(const LinkedList* ll) {
    switch (ll_kind(ll)) {
        case LL_I64:
            return ll_sum_i64(ll);
        case LL_REAL:
        case LL_CPLX: {
            LLSumCtx c = { ll, ll_kind(ll) == LL_CPLX };
            double slots[NDARRAY_MAX_THREADS * 2];
            int k = nd_parallel_reduce((size_t)ll_length(ll), ll_sum_reduce, &c, 2, slots);
            double re = 0.0, im = 0.0;
            for (int t = 0; t < k; t++) { re += slots[2 * t]; im += slots[2 * t + 1]; }
            Expr* r;
            if (c.cplx) {
                Expr* a[2] = { expr_new_real(re), expr_new_real(im) };
                r = expr_new_function(expr_new_symbol(SYM_Complex), a, 2);
            } else {
                r = expr_new_real(re);
            }
            /* Overflow / subnormal: decline, as the packed Total does, so the
             * caller's Plus keeps the value in range at machine precision. */
            if (ndred_result_degenerate(r)) { expr_free(r); return NULL; }
            return r;
        }
        default:
            return NULL;
    }
}
