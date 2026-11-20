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
#include "ndarray.h"           /* NDUnaryKernel, ndk_unary_lane_mode, nd_lane_degenerate1 */
#include "symtab.h"            /* symtab_lookup */
#include "ndreduce.h"          /* ND_PAIRWISE_BLOCK, nd_sum_block_f64, ndred_result_degenerate */
#include "ndarray_internal.h"  /* nd_parallel_reduce */
#include "checked_int.h"       /* ci_add_i64 */
#include "sym_names.h"

#include <float.h>
#include <gmp.h>
#include <math.h>
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

/* ------------------------------------------------------------------------
 *  Listable threading straight from the chunks into the chunks
 *
 *  Sin[ll], ll + 1., 2. ll, ll1 ll2 over machine reals. The streamed path
 *  (ll_thread_listable) gathers the list into a packed block, runs the packed
 *  kernel and copies the answer back: two extra passes over memory, which for a
 *  cheap kernel cost more than the kernel. Here the kernel reads the source
 *  slots and writes the destination slots directly, split across threads by
 *  nd_parallel_for exactly as the packed map is.
 *
 *  Bit-identical to the packed List by construction: the unary arm is
 *  ndarray.c's ndu_hot_chunk loop (the same k->real per element and the same
 *  lane check), and the binary arm is one IEEE add or multiply per element,
 *  which is all the packed Plus/Times do on float64 buffers. Any element that
 *  path would not answer at machine precision -- a kernel decline, a result
 *  outside the normal range (the numeric tower then wants MPFR) -- declines
 *  the whole call, and the streamed path answers exactly as the List would.
 * ---------------------------------------------------------------------- */


typedef struct {
    const LinkedList* a;     /* source (machine reals)                     */
    const LinkedList* b;     /* second operand list, or NULL               */
    LinkedList*       out;   /* reserved, same length                      */
    const NDUnaryKernel* k;  /* unary arm, or NULL for the binary arm      */
    int    mode;             /* the kernel's lane mode                     */
    bool   times;            /* binary: multiply, else add                 */
    double s;                /* binary: the scalar operand (when !b)       */
} LLMapCtx;

static inline bool ll_real_degenerate(double y) {
    double m = fabs(y);
    return !(m >= DBL_MIN && m <= DBL_MAX) && y != 0.0;
}

static bool ll_map_chunk(void* c, size_t lo, size_t hi) {
    const LLMapCtx* x = (const LLMapCtx*)c;
    for (int64_t i = (int64_t)lo; i < (int64_t)hi; ) {
        const void *pa, *pb = NULL, *po;
        int da, db = 1, dout;
        int64_t k = ll_span_at(x->a, i, &pa, &da);
        int64_t ko = ll_span_at(x->out, i, &po, &dout);      /* back run: dout = +1 */
        if (ko < k) k = ko;
        if (x->b) {
            int64_t kb = ll_span_at(x->b, i, &pb, &db);
            if (kb < k) k = kb;
        }
        if (k > (int64_t)hi - i) k = (int64_t)hi - i;
        const double* sa = (const double*)pa;
        const double* sb = (const double*)pb;
        double* o = (double*)po;
        if (x->k) {
            for (int64_t t = 0; t < k; t++) {
                double in = sa[da > 0 ? t : -t];
                if (!x->k->real(in, &o[t])) return false;
                if (x->mode != ND_LANE_NONE) {
                    double y = o[t], m = fabs(y);
                    if (!(m >= DBL_MIN && m <= DBL_MAX) &&
                        nd_lane_degenerate1(x->mode, in, 0.0, y, 0.0)) return false;
                }
            }
        } else {
            for (int64_t t = 0; t < k; t++) {
                double u = sa[da > 0 ? t : -t];
                double v = sb ? sb[db > 0 ? t : -t] : x->s;
                double y = x->times ? u * v : u + v;
                if (ll_real_degenerate(y)) return false;
                o[t] = y;
            }
        }
        i += k;
    }
    return true;
}

static const LinkedList* ll_real_operand(const Expr* e) {
    const LinkedList* l = expr_linkedlist(e);
    return (l && ll_kind(l) == LL_REAL) ? l : NULL;
}

/* A machine scalar for the binary arm: a machine Real, or an Integer the
 * conversion carries exactly. Exact 0 declines: Plus[x, 0] is x itself (it
 * keeps -0.) and Times[0, x] the exact 0, which the general path supplies. */
static bool ll_real_scalar(const Expr* e, double* v) {
    if (e->type == EXPR_REAL) { *v = e->data.real; return true; }
    if (e->type == EXPR_INTEGER && e->data.integer != 0 &&
        e->data.integer >= -(INT64_C(1) << 53) && e->data.integer <= (INT64_C(1) << 53)) {
        *v = (double)e->data.integer;
        return true;
    }
    return false;
}

Expr* ll_listable_direct(Expr* res) {
    const Expr* h = res->data.function.head;
    if (h->type != EXPR_SYMBOL) return NULL;
    const char* name = h->data.symbol.name;
    size_t argc = res->data.function.arg_count;
    LLMapCtx c;
    memset(&c, 0, sizeof c);
    if (argc == 1) {
        const SymbolDef* d = symtab_lookup(name);
        c.a = ll_real_operand(res->data.function.args[0]);
        if (!c.a || !d || d->down_values || !d->ndarray_unary_kernel) return NULL;
        c.k = (const NDUnaryKernel*)d->ndarray_unary_kernel;
        if (!c.k->real || c.k->to_int || c.k->to_real) return NULL;
        c.mode = ndk_unary_lane_mode(c.k);
    } else if (argc == 2 && (name == SYM_Plus || name == SYM_Times)) {
        const Expr* x = res->data.function.args[0];
        const Expr* y = res->data.function.args[1];
        c.times = name == SYM_Times;
        if ((c.a = ll_real_operand(x)) != NULL) {
            if (!(c.b = ll_real_operand(y)) && !ll_real_scalar(y, &c.s)) return NULL;
        } else if ((c.a = ll_real_operand(y)) != NULL) {
            if (!ll_real_scalar(x, &c.s)) return NULL;
        } else {
            return NULL;
        }
        if (c.b && ll_length(c.b) != ll_length(c.a)) return NULL;   /* tdlen: general path */
    } else {
        return NULL;
    }
    int64_t n = ll_length(c.a);
    LinkedList* out = ll_new();
    if (n > 0) {
        ll_reserve_typed(out, LL_REAL, n);
        c.out = out;
        if (!nd_parallel_for((size_t)n, ll_map_chunk, &c)) { ll_release(out); return NULL; }
    }
    return expr_new_linkedlist(out);
}
