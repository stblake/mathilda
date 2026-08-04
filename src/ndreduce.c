/* NDArray reduction & order-statistic fast paths — see ndreduce.h.
 *
 * Every summation routes through a pairwise reducer (nd_sum_strided /
 * nd_sumsq_strided), not a naive left fold: over ~1e7 float64 a left fold's
 * rounding error grows ~O(n)·eps and drifts from numpy's pairwise sum, which
 * would make the "matches the List path" parity tests flaky. Pairwise keeps the
 * error ~O(log n)·eps, matching numpy. All loops are serial: numpy's
 * sum/mean/std are single-threaded and memory-bound, so a cache-friendly serial
 * pass over the contiguous buffer already meets or beats them. */

#include "ndreduce.h"
#include "ndarray.h"
#include "ndarray_internal.h"
#include "sym_names.h"
#include "arithmetic.h"       /* make_rational — an exact integer Mean/Median */
#include "checked_int.h"      /* ci_add_i64 — exact accumulate, bail on overflow */
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>           /* Commonest::dstlms — the same message the List path prints */
#include <inttypes.h>        /* PRId64, for that message's argument */
#include <math.h>

/* Below this many summands a strided range is summed by a plain loop; above it
 * the range is split in two and each half summed recursively (pairwise). */
#define ND_PAIRWISE_BLOCK ((size_t)128)

/* Below this, an introsort's lower constant factor wins; above it an 8-pass LSD
 * radix sort (memory-bound, no comparisons) is far faster. Used by both the
 * double and the int64 sorts. */
#define ND_RADIX_MIN ((size_t)2048)

/* ------------------------------------------------------------------ helpers */

bool ndred_call_has_ndarray(const Expr* res) {
    return res && res->type == EXPR_FUNCTION &&
           res->data.function.arg_count >= 1 &&
           is_ndarray(res->data.function.args[0]);
}

/* A scalar result leaf: a bare Real, or Complex[re, im] for a complex value. */
static Expr* nd_scalar(double re, double im, bool cplx) {
    if (cplx) {
        Expr* a[2] = { expr_new_real(re), expr_new_real(im) };
        return expr_new_function(expr_new_symbol(SYM_Complex), a, 2);
    }
    return expr_new_real(re);
}

/* The real-dtype counterpart of `dt`, preserving component width (a complex
 * reduction that produces a real magnitude — Variance/Std/RMS — writes here). */
static NDType nd_real_of(NDType dt) {
    return (ndt_comp_size(dt) == sizeof(double)) ? NDT_FLOAT64 : NDT_FLOAT32;
}

/* ------------------------------------------------- exact int64 reductions
 *
 * Automatic packing infers NDT_INT64 for an exact list, so Total[Range[10^6]]
 * now arrives here as a buffer where it used to arrive as a List. Everything
 * above accumulates through `double`, which is exact only to 2^53 and answers
 * with a Real -- both wrong for an integer list, whose Total is an Integer that
 * promotes past int64 when it has to.
 *
 * So the int64 arms below are not an optimisation: without them the gate in
 * evaluate_step has to materialise every integer buffer before any reduction,
 * and Total[Range[10^6]] measured 1.55x SLOWER packed than plain. With them the
 * buffer is summed in place and the answer is the interpreter's.
 *
 * The rule everywhere: produce the EXACT answer, or return false and let
 * ndarray_delist_and_reeval re-run on the materialised List where GMP is
 * waiting. Never wrap, never round. Same contract as src/compile/. */
static bool nd_sum_i64(const void* buf, NDType dt, size_t base, size_t stride,
                       size_t count, int64_t* out) {
    const int64_t* p = (const int64_t*)buf;
    int64_t acc = 0;
    if (dt != NDT_INT64) return false;
    for (size_t i = 0; i < count; i++)
        if (ci_add_i64(acc, p[base + i * stride], &acc)) return false;
    *out = acc;
    return true;
}

static int nd_i64_cmp(const void* a, const void* b) {
    int64_t x = *(const int64_t*)a, y = *(const int64_t*)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

/* Ascending sort of an int64 array, in place.
 *
 * Radix, like the double path, and for the same reason: qsort's indirect
 * comparator call per comparison made Sort on a 10^6-element integer buffer ~9x
 * slower than the float64 radix sort on the same shape, which showed up as
 * packing being slower than a plain List for Sort[Reverse[Range[10^6]]].
 *
 * The key transform is the signed analogue of nd_d2key's: flipping the sign bit
 * maps int64 order onto uint64 order exactly (INT64_MIN -> 0, -1 -> 2^63-1,
 * 0 -> 2^63, INT64_MAX -> 2^64-1), so eight unsigned byte passes sort signed
 * values with no special-casing and no loss. */
static uint64_t nd_i2key(int64_t v) { return (uint64_t)v ^ 0x8000000000000000ULL; }
static int64_t  nd_key2i(uint64_t k) { return (int64_t)(k ^ 0x8000000000000000ULL); }

void nd_sort_i64_asc(int64_t* v, size_t n) {
    if (n < 2) return;
    if (n < ND_RADIX_MIN) { qsort(v, n, sizeof(int64_t), nd_i64_cmp); return; }
    uint64_t* a = malloc(n * sizeof(uint64_t));
    uint64_t* b = malloc(n * sizeof(uint64_t));
    if (!a || !b) { free(a); free(b); qsort(v, n, sizeof(int64_t), nd_i64_cmp); return; }
    for (size_t i = 0; i < n; i++) a[i] = nd_i2key(v[i]);
    uint64_t* src = a; uint64_t* dst = b;
    for (int pass = 0; pass < 8; pass++) {          /* one byte per pass */
        size_t count[256] = {0};
        int sh = pass * 8;
        for (size_t i = 0; i < n; i++) count[(src[i] >> sh) & 0xff]++;
        size_t off = 0;
        for (int d = 0; d < 256; d++) { size_t c = count[d]; count[d] = off; off += c; }
        for (size_t i = 0; i < n; i++) dst[count[(src[i] >> sh) & 0xff]++] = src[i];
        uint64_t* t = src; src = dst; dst = t;       /* 8 passes: ends back in `a` */
    }
    for (size_t i = 0; i < n; i++) v[i] = nd_key2i(src[i]);
    free(a); free(b);
}
#define nd_sort_i64 nd_sort_i64_asc

/* Stable argsort via bottom-up merge sort of an index array.
 *
 * The value radix sorts above cannot serve Ordering: the permutation has to be
 * carried alongside the keys, and Ordering must break ties by original position.
 * A bottom-up merge is stable BY CONSTRUCTION -- each merge takes from the left
 * run on a tie, and the left run always holds the lower original indices -- so it
 * delivers the ties-by-index order for free, with no separate tie-break pass. It
 * is O(n log n) with one O(n) scratch buffer, and its keys are read indirectly
 * through `idx`, which is exactly what an argsort is; the radix sorts' advantage
 * (no comparisons, contiguous keys) does not transfer to a permutation.
 *
 * KEYCMP(x, y) must be a strict "x sorts before y" test. Written as a macro so
 * the int64 and real variants share the merge and differ only in the comparison,
 * with no function-pointer call per comparison. */
#define ND_ARGSORT_BODY(KEYCMP)                                               \
    do {                                                                      \
        for (size_t i = 0; i < n; i++) idx[i] = (int64_t)i;                   \
        if (n < 2) return true;                                              \
        int64_t* scratch = malloc(n * sizeof(int64_t));                       \
        if (!scratch) return false;                                          \
        int64_t* src = idx;                                                   \
        int64_t* dst = scratch;                                               \
        for (size_t width = 1; width < n; width *= 2) {                       \
            for (size_t lo = 0; lo < n; lo += 2 * width) {                    \
                size_t mid = lo + width;   if (mid > n) mid = n;              \
                size_t hi  = lo + 2*width; if (hi  > n) hi  = n;              \
                size_t i = lo, j = mid, k = lo;                              \
                while (i < mid && j < hi) {                                   \
                    /* strict "right before left" -> ties keep left (lower i) */ \
                    if (KEYCMP(src[j], src[i])) dst[k++] = src[j++];          \
                    else                         dst[k++] = src[i++];         \
                }                                                             \
                while (i < mid) dst[k++] = src[i++];                          \
                while (j < hi)  dst[k++] = src[j++];                          \
            }                                                                 \
            int64_t* t = src; src = dst; dst = t;                             \
        }                                                                     \
        if (src != idx) memcpy(idx, src, n * sizeof(int64_t));               \
        free(scratch);                                                        \
        return true;                                                          \
    } while (0)

bool nd_argsort_i64(const int64_t* v, size_t n, int64_t* idx) {
#define ND_I64_LT(pa, pb) (v[(pa)] < v[(pb)])
    ND_ARGSORT_BODY(ND_I64_LT);
#undef ND_I64_LT
}

bool nd_argsort_real(const double* v, size_t n, int64_t* idx) {
#define ND_REAL_LT(pa, pb) (v[(pa)] < v[(pb)])
    ND_ARGSORT_BODY(ND_REAL_LT);
#undef ND_REAL_LT
}

#undef ND_ARGSORT_BODY

/* Max/Min over an int64 range: always exact, no overflow to worry about. */
static int64_t nd_extreme_i64(const int64_t* p, size_t n, bool want_max) {
    int64_t b = p[0];
    for (size_t i = 1; i < n; i++)
        b = want_max ? (p[i] > b ? p[i] : b) : (p[i] < b ? p[i] : b);
    return b;
}

/* Pairwise sum of the `count` elements buf[base], buf[base+stride], ... into
 * (*re, *im). Strided so the same routine sums a contiguous run (stride == 1,
 * full/flatten reductions) or one column of a leading-axis reduction
 * (stride == trailing size). float64 gets a tight raw-pointer inner loop. */
static void nd_sum_strided(const void* buf, NDType dt, size_t base,
                           size_t stride, size_t count, double* re, double* im) {
    if (count <= ND_PAIRWISE_BLOCK) {
        double sr = 0.0, si = 0.0;
        if (dt == NDT_FLOAT64) {
            const double* d = (const double*)buf;
            if (stride == 1) {
                /* Contiguous: 8 independent accumulators break the loop-carried
                 * add dependency so the compiler can vectorize/pipeline the
                 * reduction without needing -ffast-math (FP add isn't
                 * associative, so a single accumulator would stay scalar). */
                const double* p = d + base;
                double a0=0,a1=0,a2=0,a3=0,a4=0,a5=0,a6=0,a7=0;
                size_t i = 0;
                for (; i + 8 <= count; i += 8) {
                    a0+=p[i];   a1+=p[i+1]; a2+=p[i+2]; a3+=p[i+3];
                    a4+=p[i+4]; a5+=p[i+5]; a6+=p[i+6]; a7+=p[i+7];
                }
                for (; i < count; i++) a0 += p[i];
                sr = ((a0+a1)+(a2+a3)) + ((a4+a5)+(a6+a7));
            } else {
                for (size_t i = 0; i < count; i++) sr += d[base + i * stride];
            }
        } else {
            for (size_t i = 0; i < count; i++) {
                double r, m;
                ndt_get(buf, base + i * stride, dt, &r, &m);
                sr += r; si += m;
            }
        }
        *re = sr; *im = si;
        return;
    }
    size_t half = count / 2;
    double r1, i1, r2, i2;
    nd_sum_strided(buf, dt, base, stride, half, &r1, &i1);
    nd_sum_strided(buf, dt, base + half * stride, stride, count - half, &r2, &i2);
    *re = r1 + r2; *im = i1 + i2;
}

/* Pairwise sum of squared deviations |buf[k] - (mr,mi)|^2 over the same strided
 * range (mr = mi = 0 gives the plain sum of |buf[k]|^2, used by RMS). Returns a
 * real accumulation. */
static double nd_sumsq_strided(const void* buf, NDType dt, size_t base,
                               size_t stride, size_t count, double mr, double mi) {
    if (count <= ND_PAIRWISE_BLOCK) {
        double s = 0.0;
        if (dt == NDT_FLOAT64 && stride == 1) {       /* contiguous real */
            const double* p = (const double*)buf + base;
            double a0=0,a1=0,a2=0,a3=0;               /* independent accumulators */
            size_t i = 0;
            for (; i + 4 <= count; i += 4) {
                double d0=p[i]-mr, d1=p[i+1]-mr, d2=p[i+2]-mr, d3=p[i+3]-mr;
                a0+=d0*d0; a1+=d1*d1; a2+=d2*d2; a3+=d3*d3;
            }
            for (; i < count; i++) { double dr = p[i] - mr; a0 += dr * dr; }
            return (a0 + a1) + (a2 + a3);
        }
        for (size_t i = 0; i < count; i++) {
            double r, m;
            ndt_get(buf, base + i * stride, dt, &r, &m);
            double dr = r - mr, di = m - mi;
            s += dr * dr + di * di;
        }
        return s;
    }
    size_t half = count / 2;
    return nd_sumsq_strided(buf, dt, base, stride, half, mr, mi) +
           nd_sumsq_strided(buf, dt, base + half * stride, stride, count - half, mr, mi);
}

/* Product of dims[lo..hi). */
static size_t nd_dim_prod(const int64_t* dims, int lo, int hi) {
    size_t p = 1;
    for (int i = lo; i < hi; i++) p *= (size_t)dims[i];
    return p;
}

/* ------------------------------- parallel reduction over a contiguous run */
/* Large flat reductions (a vector Total/Mean/Variance/…, or a full flatten) are
 * memory-bound: a single core tops out well below the machine's memory
 * bandwidth, so we split the range across threads (each folds its own private
 * partial) and combine — reaching bandwidth on 3-4 cores. Small arrays and
 * thread-less builds fall back to one serial chunk inside nd_parallel_reduce. */

typedef struct { const void* buf; NDType dt; } nd_sum_ctx;
static void nd_sum_reduce(void* c, size_t lo, size_t hi, double* slot) {
    const nd_sum_ctx* x = (const nd_sum_ctx*)c;
    nd_sum_strided(x->buf, x->dt, lo, 1, hi - lo, &slot[0], &slot[1]);
}
/* Sum of `n` contiguous elements, threaded, into (*re, *im). */
static void nd_full_sum(const void* buf, NDType dt, size_t n, double* re, double* im) {
    nd_sum_ctx c = { buf, dt };
    double slots[NDARRAY_MAX_THREADS * 2];
    int k = nd_parallel_reduce(n, nd_sum_reduce, &c, 2, slots);
    double r = 0.0, m = 0.0;
    for (int t = 0; t < k; t++) { r += slots[2 * t]; m += slots[2 * t + 1]; }
    *re = r; *im = m;
}

typedef struct { const void* buf; NDType dt; double mr, mi; } nd_sq_ctx;
static void nd_sq_reduce(void* c, size_t lo, size_t hi, double* slot) {
    const nd_sq_ctx* x = (const nd_sq_ctx*)c;
    slot[0] = nd_sumsq_strided(x->buf, x->dt, lo, 1, hi - lo, x->mr, x->mi);
}
/* Sum of |x-(mr,mi)|^2 over `n` contiguous elements, threaded. */
static double nd_full_sumsq(const void* buf, NDType dt, size_t n, double mr, double mi) {
    nd_sq_ctx c = { buf, dt, mr, mi };
    double slots[NDARRAY_MAX_THREADS];
    int k = nd_parallel_reduce(n, nd_sq_reduce, &c, 1, slots);
    double s = 0.0;
    for (int t = 0; t < k; t++) s += slots[t];
    return s;
}

/* ---------------------------------------- shared buffer helpers (internal) */

void nd_gather_real(const void* buf, NDType dt, size_t base, size_t stride,
                    size_t count, double* out) {
    if (dt == NDT_FLOAT64 && stride == 1) {
        const double* d = (const double*)buf;
        for (size_t i = 0; i < count; i++) out[i] = d[base + i];
        return;
    }
    for (size_t i = 0; i < count; i++) {
        double r, m;
        ndt_get(buf, base + i * stride, dt, &r, &m);
        out[i] = r;
    }
}

double nd_select_kth(double* s, size_t n, size_t k) {
    size_t lo = 0, hi = n - 1;
    while (lo < hi) {
        double pivot = s[lo + (hi - lo) / 2];
        size_t i = lo, j = hi;
        while (i <= j) {                      /* Hoare partition */
            while (s[i] < pivot) i++;
            while (s[j] > pivot) j--;
            if (i <= j) {
                double t = s[i]; s[i] = s[j]; s[j] = t;
                i++;
                if (j == 0) break;            /* size_t underflow guard */
                j--;
            }
        }
        if (k <= j) hi = j;
        else if (k >= i) lo = i;
        else break;                           /* j < k < i: s[k] is settled */
    }
    return s[k];
}

/* Inlined quicksort on doubles — a median-of-three introsort with an
 * insertion-sort cutoff. Avoids qsort's per-comparison indirect call (the
 * dominant cost when sorting millions of machine doubles). */
#define ND_SORT_INSERTION 24
static void nd_qsort(double* s, size_t lo, size_t hi) {
    while (hi - lo > ND_SORT_INSERTION) {
        size_t mid = lo + (hi - lo) / 2;          /* median-of-3 -> s[lo] pivot */
        if (s[mid] < s[lo]) { double t = s[mid]; s[mid] = s[lo]; s[lo] = t; }
        if (s[hi] < s[lo])  { double t = s[hi];  s[hi]  = s[lo]; s[lo] = t; }
        if (s[hi] < s[mid]) { double t = s[hi];  s[hi]  = s[mid]; s[mid] = t; }
        double pivot = s[mid];
        size_t i = lo, j = hi;
        for (;;) {
            do { i++; } while (s[i] < pivot);
            do { j--; } while (s[j] > pivot);
            if (i >= j) break;
            double t = s[i]; s[i] = s[j]; s[j] = t;
        }
        /* Recurse into the smaller partition, loop on the larger (bounded stack). */
        if (j - lo < hi - j) { nd_qsort(s, lo, j); lo = j + 1; }
        else                 { nd_qsort(s, j + 1, hi); hi = j; }
    }
    for (size_t i = lo + 1; i <= hi; i++) {        /* insertion sort the tail */
        double v = s[i];
        size_t k = i;
        while (k > lo && s[k - 1] > v) { s[k] = s[k - 1]; k--; }
        s[k] = v;
    }
}

/* Order-preserving double <-> uint64 key: flip the sign bit for positives, and
 * all bits for negatives, so unsigned ascending order == double ascending order
 * (handles -0.0/+0.0; NaNs cluster at an end). */
static uint64_t nd_d2key(double d) {
    uint64_t u; memcpy(&u, &d, sizeof u);
    return (u >> 63) ? ~u : (u | 0x8000000000000000ULL);
}
static double nd_key2d(uint64_t k) {
    uint64_t u = (k & 0x8000000000000000ULL) ? (k ^ 0x8000000000000000ULL) : ~k;
    double d; memcpy(&d, &u, sizeof d);
    return d;
}

/* LSD radix sort of the doubles in `s` (ascending). Falls back to introsort if
 * the scratch allocation fails. */
static void nd_radix_sort(double* s, size_t n) {
    uint64_t* a = malloc(n * sizeof(uint64_t));
    uint64_t* b = malloc(n * sizeof(uint64_t));
    if (!a || !b) { free(a); free(b); nd_qsort(s, 0, n - 1); return; }
    for (size_t i = 0; i < n; i++) a[i] = nd_d2key(s[i]);
    uint64_t* src = a; uint64_t* dst = b;
    for (int pass = 0; pass < 8; pass++) {          /* one byte per pass */
        size_t count[256] = {0};
        int sh = pass * 8;
        for (size_t i = 0; i < n; i++) count[(src[i] >> sh) & 0xff]++;
        size_t off = 0;
        for (int d = 0; d < 256; d++) { size_t c = count[d]; count[d] = off; off += c; }
        for (size_t i = 0; i < n; i++) dst[count[(src[i] >> sh) & 0xff]++] = src[i];
        uint64_t* t = src; src = dst; dst = t;      /* 8 passes: ends back in `a` */
    }
    for (size_t i = 0; i < n; i++) s[i] = nd_key2d(src[i]);
    free(a); free(b);
}

void nd_sort_ascending(double* s, size_t n) {
    if (n < 2) return;
    if (n < ND_RADIX_MIN) nd_qsort(s, 0, n - 1);
    else nd_radix_sort(s, n);
}

/* ----------------------------------------------------------------- Median */

/* Median of `count` reals gathered at (base, stride). Uses quickselect (O(n))
 * rather than a full sort. Even length averages the two central order stats. */
static double nd_median_of(const void* buf, NDType dt, size_t base,
                           size_t stride, size_t count) {
    double* s = malloc(sizeof(double) * count);
    nd_gather_real(buf, dt, base, stride, count, s);
    double med;
    if (count % 2 == 1) {
        med = nd_select_kth(s, count, count / 2);
    } else {
        double hi = nd_select_kth(s, count, count / 2);      /* upper middle */
        double lo = s[0];                                    /* lower half max */
        for (size_t i = 1; i < count / 2; i++)
            if (s[i] > lo) lo = s[i];
        med = 0.5 * (lo + hi);
    }
    free(s);
    return med;
}

Expr* ndred_median(Expr* res) {
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    int rank = a->data.ndarray.rank;
    NDType dt = a->data.ndarray.dtype;
    if (ndt_is_complex(dt) || rank > 2) return ndarray_delist_and_reeval(res);
    const int64_t* dims = a->data.ndarray.dims;
    const void* buf = a->data.ndarray.data;

    if (rank == 1) {
        if (dt == NDT_INT64) {
            /* Median[Range[300]] is 301/2. Quickselect over doubles would be
             * both inexact past 2^53 and the wrong head, so gather exactly. */
            size_t n = (size_t)dims[0];
            if (n == 0) return ndarray_delist_and_reeval(res);
            const int64_t* p = (const int64_t*)buf;
            int64_t* t = malloc(sizeof(int64_t) * n);
            if (!t) return ndarray_delist_and_reeval(res);
            memcpy(t, p, sizeof(int64_t) * n);
            nd_sort_i64(t, n);
            Expr* out;
            if (n % 2 == 1) out = expr_new_integer(t[n / 2]);
            else {
                int64_t s;
                out = ci_add_i64(t[n / 2 - 1], t[n / 2], &s)
                        ? NULL : make_rational(s, 2);
            }
            free(t);
            return out ? out : ndarray_delist_and_reeval(res);
        }
        return expr_new_real(nd_median_of(buf, dt, 0, 1, (size_t)dims[0]));
    }
    if (dt == NDT_INT64) return ndarray_delist_and_reeval(res);   /* columnwise */

    /* Matrix: columnwise median -> rank-1 vector of length ncols. */
    size_t nrows = (size_t)dims[0], ncols = (size_t)dims[1];
    NDType odt = nd_real_of(dt);
    void* out = malloc(ndt_elem_size(odt) * ncols);
    if (!out) return ndarray_delist_and_reeval(res);
    for (size_t j = 0; j < ncols; j++)
        ndt_set(out, j, odt, nd_median_of(buf, dt, j, ncols, nrows), 0.0);
    int64_t odims[1] = { (int64_t)ncols };
    return expr_new_ndarray_like(a, 1, odims, out, odt);
}

/* ---------------------------------------------------- RankedMin / RankedMax
 *
 * RankedMin[v, n] selects the n-th SMALLEST element of the rank-1 machine
 * vector v (n<0 counts from the largest); RankedMax the n-th LARGEST. The two
 * differ only by the sign of n (RankedMax[v, n] == RankedMin[v, -n]), so both
 * reduce to one ASCENDING rank r in [1, m] and the r-th order statistic taken
 * straight off the buffer: an int64 vector selects EXACTLY (nd_sort_i64_asc, so
 * the answer is an Integer, never a rounded Real, past 2^53), a real vector via
 * O(m) quickselect (nd_select_kth). A complex dtype, rank > 1, non-integer or
 * out-of-range n degrades to ndarray_delist_and_reeval — the List path then
 * answers, identically. These two entry points are also the Compile ND_REDS
 * delegates for the heads. */
static Expr* ndred_ranked(Expr* res, bool is_max) {
    if (res->data.function.arg_count != 2) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    Expr* nexpr = res->data.function.args[1];
    if (nexpr->type != EXPR_INTEGER || nexpr->data.integer == 0)
        return ndarray_delist_and_reeval(res);
    NDType dt = a->data.ndarray.dtype;
    if (ndt_is_complex(dt) || a->data.ndarray.rank != 1)
        return ndarray_delist_and_reeval(res);
    size_t m = (size_t)a->data.ndarray.dims[0];
    if (m == 0) return ndarray_delist_and_reeval(res);

    /* Ascending rank r (1-based). RankedMax negates n first, so both heads
     * share one formula: n>0 -> r=n, n<0 -> r=m+n+1. */
    int64_t n = is_max ? -nexpr->data.integer : nexpr->data.integer;
    int64_t r = (n > 0) ? n : ((int64_t)m + n + 1);
    if (r < 1 || r > (int64_t)m) return ndarray_delist_and_reeval(res);
    size_t k = (size_t)(r - 1);   /* 0-based order statistic */

    const void* buf = a->data.ndarray.data;
    if (dt == NDT_INT64) {
        const int64_t* p = (const int64_t*)buf;
        int64_t* t = malloc(sizeof(int64_t) * m);
        if (!t) return ndarray_delist_and_reeval(res);
        memcpy(t, p, sizeof(int64_t) * m);
        nd_sort_i64_asc(t, m);
        Expr* out = expr_new_integer(t[k]);
        free(t);
        return out;
    }
    double* s = malloc(sizeof(double) * m);
    if (!s) return ndarray_delist_and_reeval(res);
    nd_gather_real(buf, dt, 0, 1, m, s);
    double val = nd_select_kth(s, m, k);
    free(s);
    return expr_new_real(val);
}

Expr* ndred_ranked_min(Expr* res) { return ndred_ranked(res, false); }
Expr* ndred_ranked_max(Expr* res) { return ndred_ranked(res, true); }

/* --------------------------------------------------------------- Quartiles */

/* The reductions whose exact integer answer is a Rational or a root -- so no
 * int64 buffer can hold it and there is nothing to gain from trying. The gate in
 * evaluate_step already materialises an integer argument for every head not on
 * its INT64_OK list, so this guard is for the surface the user opts into
 * directly: NDArray[..., DataType -> "int64"], and a buffer arriving from
 * Compile[]. Without it Variance[intArray] would answer with a Real. */
static bool nd_int64_degrade(const Expr* a) {
    return a && a->type == EXPR_NDARRAY && a->data.ndarray.dtype == NDT_INT64;
}

/* Mathematica's default Quartiles parameters {{1/2, 0}, {0, 1}}: for the k-th
 * quantile q, h = 1/2 + n*q, clamped, then linear interpolation between the
 * bracketing order statistics of the ascending-sorted data. */
Expr* ndred_quartiles(Expr* res) {
    /* Only the default (1-arg) rank-1 real case is fast-pathed; a custom
     * parameter matrix or higher rank degrades to the exact List method. */
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    if (a->data.ndarray.rank != 1 || ndt_is_complex(a->data.ndarray.dtype) ||
        nd_int64_degrade(a))
        return ndarray_delist_and_reeval(res);

    size_t n = (size_t)a->data.ndarray.dims[0];
    NDType dt = a->data.ndarray.dtype;
    double* s = malloc(sizeof(double) * n);
    nd_gather_real(a->data.ndarray.data, dt, 0, 1, n, s);
    nd_sort_ascending(s, n);

    const double q[3] = { 0.25, 0.5, 0.75 };
    NDType odt = nd_real_of(dt);
    void* out = malloc(ndt_elem_size(odt) * 3);
    for (int k = 0; k < 3; k++) {
        double h = 0.5 + (double)n * q[k];
        double val;
        if (h <= 1.0) val = s[0];
        else if (h >= (double)n) val = s[n - 1];
        else {
            int64_t j = (int64_t)floor(h);
            if (j < 1) j = 1;
            if (j >= (int64_t)n) j = (int64_t)n - 1;
            double g = h - (double)j;                    /* fractional part */
            val = s[j - 1] + g * (s[j] - s[j - 1]);       /* 1-based order stats */
        }
        ndt_set(out, (size_t)k, odt, val, 0.0);
    }
    free(s);
    int64_t odims[1] = { 3 };
    return expr_new_ndarray_like(a, 1, odims, out, odt);
}

/* ---------------------------------------------------- Moving statistics */

/* Shared preamble for MovingAverage/MovingMedian: require a 2-arg call on a
 * real rank-1 array with a positive integer window r in [1, n]. Returns r via
 * *r_out and the array via *a_out, or false (caller degrades). */
static bool nd_moving_window(Expr* res, Expr** a_out, size_t* r_out) {
    if (res->data.function.arg_count != 2) return false;
    Expr* a = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];
    if (a->data.ndarray.rank != 1 || ndt_is_complex(a->data.ndarray.dtype) ||
        nd_int64_degrade(a))
        return false;
    if (spec->type != EXPR_INTEGER) return false;   /* weight list / bignum -> List path */
    int64_t r = spec->data.integer;
    size_t n = (size_t)a->data.ndarray.dims[0];
    if (r < 1 || (size_t)r > n) return false;
    *a_out = a; *r_out = (size_t)r;
    return true;
}

/* --- Tally ------------------------------------------------------------
 *
 * The generic Tally already hashes, so this is not an algorithmic fix -- it is
 * a boxing fix. Tally[list] reads list->args[i], so a packed buffer had to be
 * materialised into one Expr per element before the builtin could start, and
 * then every probe hashed and compared Expr NODES. Measured over 10^6 packed
 * int64: 147 ms, of which ~51 ms was the materialisation and ~96 ms the boxed
 * hash walk -- and 96 ms was ALSO what the same data cost as a plain List, so
 * being packed made Tally slower than not being packed.
 *
 * Here the key is the machine word itself. Open addressing, linear probing,
 * power-of-two capacity; the unique values are appended in first-appearance
 * order and the table stores indices into that array, which is what makes the
 * output order match the List path by construction rather than by a sort.
 */

/* 64-bit finalizer from splitmix64: cheap, and it avalanches the low bits,
 * which matters because the mask below keeps only those. Small consecutive
 * integers -- exactly what Tally is usually given -- would otherwise collide in
 * long runs under linear probing. */
static uint64_t tally_mix(uint64_t z) {
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

/* Integer keys almost always come from a bounded range -- RandomInteger bounds,
 * category codes, image samples, digits -- and for those the hash is pure
 * overhead: the value IS an index. counts[v - min] costs one random write into a
 * table sized by the RANGE rather than by the distinct count, with no mixing, no
 * probe loop and no key comparison, so the pass runs at memory bandwidth.
 *
 * Measured over 10^7 int64 (n = 10^7), against the hash path below:
 *
 *      distinct     hash      direct
 *         10^2    23.5 ms    17.9 ms
 *         10^4   105.1 ms    22.3 ms      4.7x
 *         10^6  1239.2 ms    66.4 ms     18.7x
 *
 * The 10^6 row is the one that matters: there the hash table (16 MB) and its key
 * array (8 MB) both fall out of L3, so every element paid TWO dependent DRAM
 * misses -- the probe, then the key comparison it depends on. Direct indexing
 * has one table and no dependent second load.
 *
 * The cost is a min/max pass before anything can be sized, so this reads the
 * buffer twice. That is why the 10^2 row gains so little: two bandwidth-bound
 * passes (~19 ms for 80 MB) is already the floor there, and the hash was nearly
 * at it. It is still the right default -- the floor does not degrade with the
 * distinct count, which is the whole failure mode being fixed.
 *
 * Declines (returns NULL, caller falls through to the hash) when the range is
 * wider than the input it would replace, so the table is never bigger than the
 * buffer already in hand. Returns the {value, count} List on success. */
#define TALLY_DIRECT_MAX_RANGE ((uint64_t)1 << 25)   /* 32M bins; see guard */

/* Both halves of an int64 tally -- the value and its count -- are int64, so the
 * answer is a rank-2 machine matrix, not nuniq boxed pairs. That is not a
 * micro-optimisation at this scale: building 10^6 {value, count} pairs as Expr
 * nodes costs 0.65 s on its own, which was the ENTIRE remaining cost of the
 * 10^6-distinct tally once the counting had been fixed. Interleaving into one
 * buffer is a memcpy by comparison, and it keeps the result packed for whatever
 * consumes it next. expr_new_ndarray_like inherits `src`'s presentation, so a
 * packed List in gives a packed List out (never a visible NDArray[...]).
 *
 * Float64 tallies deliberately keep the boxed-pair List: there the value is a
 * Real and the count an Integer, and one buffer cannot hold two heads without
 * turning the counts into Reals. */
static Expr* tally_pack_i64(const Expr* src, int64_t* buf, size_t nuniq) {
    int64_t dims[2] = { (int64_t)nuniq, 2 };
    return expr_new_ndarray_like(src, 2, dims, buf, NDT_INT64);
}

/* The counted distinct values of a rank-1 machine buffer, in FIRST-APPEARANCE
 * order. This is the substrate under both Tally and Commonest: they need the
 * same three arrays and differ only in what they build from them. One
 * implementation is not merely tidier -- Commonest breaks a count tie by first
 * appearance, so a second counting routine that inserted in a different order
 * would answer differently on ties without either being "wrong", and nothing
 * downstream would catch it. */
typedef struct {
    uint64_t* keys;    /* the distinct machine words, first-appearance order */
    int64_t*  cnts;    /* keys[u]'s multiplicity */
    size_t    nuniq;
} NDWordTally;

static void ndwt_free(NDWordTally* t) {
    free(t->keys); free(t->cnts);
    t->keys = NULL; t->cnts = NULL; t->nuniq = 0;
}

/* Rank 1 only: a tally of a matrix tallies its ROWS, which are not machine
 * words. Complex dtypes go the generic way for the same reason, and the narrow
 * dtypes because widening them to a key would merge values the List path keeps
 * apart. */
static bool nd_word_keyable(const Expr* a) {
    NDType dt = a->data.ndarray.dtype;
    return a->data.ndarray.rank == 1 && !ndt_is_complex(dt) &&
           (dt == NDT_INT64 || dt == NDT_FLOAT64);
}

static bool nd_tally_direct_i64(const int64_t* iv, size_t n, NDWordTally* out) {
    int64_t mn = iv[0], mx = iv[0];
    for (size_t i = 1; i < n; i++) {
        if (iv[i] < mn) mn = iv[i];
        if (iv[i] > mx) mx = iv[i];
    }
    /* Unsigned subtraction: mx - mn as int64 overflows for a genuinely full-width
     * spread (mn near INT64_MIN, mx near INT64_MAX), and that is UB, not a large
     * number. In uint64 the difference of two int64 is exact whenever mx >= mn. */
    uint64_t span = (uint64_t)mx - (uint64_t)mn;
    if (span >= TALLY_DIRECT_MAX_RANGE || span + 1 > (uint64_t)n) return false;
    uint64_t range = span + 1;

    int64_t*  cnt = calloc((size_t)range, sizeof(int64_t));
    size_t    ocap = 1024, nuniq = 0;
    uint32_t* ord = malloc(sizeof(uint32_t) * ocap);   /* offsets, in first-seen order */
    if (!cnt || !ord) { free(cnt); free(ord); return false; }

    for (size_t i = 0; i < n; i++) {
        uint64_t idx = (uint64_t)iv[i] - (uint64_t)mn;
        if (cnt[idx]++ == 0) {                         /* first appearance */
            if (nuniq == ocap) {
                size_t nc = ocap * 2;
                uint32_t* o2 = realloc(ord, sizeof(uint32_t) * nc);
                if (!o2) { free(cnt); free(ord); return false; }
                ord = o2; ocap = nc;
            }
            ord[nuniq++] = (uint32_t)idx;              /* idx < 2^25, fits */
        }
    }

    uint64_t* keys = malloc(sizeof(uint64_t) * (nuniq ? nuniq : 1));
    int64_t*  cs   = malloc(sizeof(int64_t)  * (nuniq ? nuniq : 1));
    if (!keys || !cs) { free(cnt); free(ord); free(keys); free(cs); return false; }
    for (size_t u = 0; u < nuniq; u++) {
        keys[u] = (uint64_t)(mn + (int64_t)ord[u]);
        cs[u]   = cnt[ord[u]];
    }
    free(cnt); free(ord);
    out->keys = keys; out->cnts = cs; out->nuniq = nuniq;
    return true;
}

/* One table entry. Keeping the key beside its count is the point: the probe
 * loads the key and increments the count in the SAME cache line, where the
 * previous index-into-a-separate-key-array made every element pay a second
 * dependent miss. Worth 2.8x on its own at 10^6 distinct (1239 -> 442 ms).
 * cnt == 0 marks an empty slot -- every live entry has cnt >= 1. */
typedef struct { uint64_t key; int64_t cnt; } TallyEnt;

/* Returns false when the buffer cannot be keyed faithfully or an allocation
 * fails; the caller then hands the whole call back to the List path. */
static bool nd_tally_hash(const void* data, NDType dt, size_t n, NDWordTally* out) {
    const int64_t* iv = (const int64_t*)data;
    const double*  dv = (const double*)data;

    /* Sized by the DISTINCT count, which is discovered, not by n. Sizing for n
     * up front is the obvious thing and it is badly wrong at the scale this
     * exists for: a 10^7-element tally over 10^4 distinct values wanted a
     * 2^25-slot table -- 268 MB to zero and then random-probe, for 10^4 live
     * entries -- and that allocation, not the hashing, was most of the time.
     * Growing costs one rehash of the LIVE keys per doubling, so the whole
     * sequence of rehashes is O(nuniq), against O(n) probes either way. */
    size_t cap = 1024, nuniq = 0, ocap = 1024;
    size_t    mask = cap - 1;
    TallyEnt* tab = calloc(cap, sizeof(TallyEnt));
    uint32_t* ord = malloc(sizeof(uint32_t) * ocap);   /* table slots, first-seen order */
    if (!tab || !ord) { free(tab); free(ord); return false; }
    for (size_t i = 0; i < n; i++) {
        uint64_t k;
        if (dt == NDT_INT64) {
            k = (uint64_t)iv[i];
        } else {
            double x = dv[i];
            /* A non-finite element cannot be keyed on its bits and stay
             * faithful: NaN != NaN, so the Expr path gives every NaN its own
             * bucket while identical bit patterns would share one. Rare enough
             * that handing the whole call back is the right trade. */
            if (!(x == x) || x > 1.7976931348623157e308 || x < -1.7976931348623157e308) {
                free(tab); free(ord);
                return false;
            }
            /* The key is the bit pattern, deliberately UNnormalised, so -0.0
             * and 0.0 stay two tallies. Folding them together is the obvious
             * thing to do and it is wrong here: Mathilda's List path compares
             * them with expr_eq, which distinguishes them, so normalising made
             * Tally[{0., -0., 1.}] answer {{0.,2},{1.,1}} over a buffer and
             * {{0.,1},{-0.,1},{1.,1}} over the identical plain list. Whether
             * that List behaviour is itself right is a separate question -- what
             * cannot happen is the two disagreeing. The differential test in
             * tests/test_packed_list.c is what caught it. */
            memcpy(&k, &x, sizeof(k));
        }
        size_t h = (size_t)tally_mix(k) & mask;
        for (;;) {
            TallyEnt* e = &tab[h];        /* key and count in ONE cache line */
            if (e->cnt == 0) {            /* empty slot -> first sighting */
                e->key = k; e->cnt = 1;
                if (nuniq == ocap) {
                    size_t nc = ocap * 2;
                    uint32_t* o2 = realloc(ord, sizeof(uint32_t) * nc);
                    if (!o2) { free(tab); free(ord); return false; }
                    ord = o2; ocap = nc;
                }
                ord[nuniq++] = (uint32_t)h;
                break;
            }
            if (e->key == k) { e->cnt++; break; }
            h = (h + 1) & mask;
        }

        if (nuniq * 10 >= cap * 7) {               /* grow the table, rehash */
            size_t nc = cap * 2, nmask = nc - 1;
            /* `ord` holds slot indices as uint32. Passing 2^32 slots needs ~3e9
             * DISTINCT keys, so this is unreachable for any array that fits in
             * memory -- but truncating silently would corrupt the output order,
             * so hand the call back rather than assume. */
            if (nc > 0xFFFFFFFFu) { free(tab); free(ord); return false; }
            TallyEnt* t2 = calloc(nc, sizeof(TallyEnt));
            if (!t2) { free(tab); free(ord); return false; }
            /* Rehash through `ord` so first-appearance order is carried across
             * the resize -- the slot indices it holds are all about to move. */
            for (size_t u = 0; u < nuniq; u++) {
                TallyEnt src = tab[ord[u]];
                size_t g = (size_t)tally_mix(src.key) & nmask;
                while (t2[g].cnt) g = (g + 1) & nmask;
                t2[g] = src;
                ord[u] = (uint32_t)g;
            }
            free(tab);
            tab = t2; cap = nc; mask = nmask;
        }
    }

    uint64_t* keys = malloc(sizeof(uint64_t) * (nuniq ? nuniq : 1));
    int64_t*  cs   = malloc(sizeof(int64_t)  * (nuniq ? nuniq : 1));
    if (!keys || !cs) { free(tab); free(ord); free(keys); free(cs); return false; }
    for (size_t u = 0; u < nuniq; u++) {
        keys[u] = tab[ord[u]].key;
        cs[u]   = tab[ord[u]].cnt;
    }
    free(tab); free(ord);
    out->keys = keys; out->cnts = cs; out->nuniq = nuniq;
    return true;
}

/* Integer keys try direct indexing first; it declines on a range too wide to
 * index, and only then is the hash the right structure. Borrows `a`; on true
 * the caller owns *out and must ndwt_free it. */
static bool nd_tally_words(const Expr* a, NDWordTally* out) {
    NDType      dt   = a->data.ndarray.dtype;
    size_t      n    = (size_t)a->data.ndarray.dims[0];
    const void* data = a->data.ndarray.data;
    out->keys = NULL; out->cnts = NULL; out->nuniq = 0;
    if (dt == NDT_INT64 && nd_tally_direct_i64((const int64_t*)data, n, out))
        return true;
    return nd_tally_hash(data, dt, n, out);
}

Expr* ndred_tally(Expr* res) {
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    if (!nd_word_keyable(a)) return ndarray_delist_and_reeval(res);

    size_t n = (size_t)a->data.ndarray.dims[0];
    if (n == 0) return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);

    NDWordTally t;
    if (!nd_tally_words(a, &t)) return ndarray_delist_and_reeval(res);
    size_t nuniq = t.nuniq;

    if (a->data.ndarray.dtype == NDT_INT64) {   /* uniform int64 -> packed matrix */
        int64_t* buf = malloc(sizeof(int64_t) * 2 * (nuniq ? nuniq : 1));
        if (!buf) { ndwt_free(&t); return ndarray_delist_and_reeval(res); }
        for (size_t u = 0; u < nuniq; u++) {
            buf[2 * u]     = (int64_t)t.keys[u];
            buf[2 * u + 1] = t.cnts[u];
        }
        ndwt_free(&t);
        return tally_pack_i64(a, buf, nuniq);
    }

    Expr** out = malloc(sizeof(Expr*) * (nuniq ? nuniq : 1));
    if (!out) { ndwt_free(&t); return ndarray_delist_and_reeval(res); }
    for (size_t u = 0; u < nuniq; u++) {
        Expr* pair[2];
        double x;
        memcpy(&x, &t.keys[u], sizeof(x));
        pair[0] = expr_new_real(x);
        pair[1] = expr_new_integer(t.cnts[u]);
        out[u] = expr_new_function(expr_new_symbol(SYM_List), pair, 2);
    }
    ndwt_free(&t);
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), out, nuniq);
    free(out);
    return list;
}

/* ------------------------------------------------------------ Commonest ---
 *
 * Commonest IS a tally plus a selection, and before this it paid for the tally
 * the expensive way. On the workload that prompted it -- 10^7 int64, four
 * distinct values -- the transparency gate materialised 20 million boxed Expr
 * across the two calls, and Commonest cost 1.19 s against 22.6 ms for Tally of
 * the identical buffer. Sharing nd_tally_words closes the whole gap: what
 * remains is one O(u log u) sort over the DISTINCT values, u = 4 here, which is
 * unmeasurable beside the count pass.
 *
 * The selection is the List path's, transcribed: order the distinct values by
 * count descending and first appearance ascending, take the first target_n, then
 * restore first-appearance order among those -- which is why the two agree on
 * ties rather than merely on the multiset. */
typedef struct { uint64_t key; int64_t cnt; size_t first; } CommonWord;

static int commonest_by_count(const void* pa, const void* pb) {
    const CommonWord* a = (const CommonWord*)pa;
    const CommonWord* b = (const CommonWord*)pb;
    if (a->cnt != b->cnt) return (b->cnt > a->cnt) ? 1 : -1;
    return (a->first > b->first) ? 1 : -1;
}

static int commonest_by_first(const void* pa, const void* pb) {
    const CommonWord* a = (const CommonWord*)pa;
    const CommonWord* b = (const CommonWord*)pb;
    if (a->first == b->first) return 0;
    return (a->first > b->first) ? 1 : -1;
}

/* Reads `n` out of Commonest's second argument. Sets *n to the count and
 * *upto to whether the UpTo[] wrapper suppressed the too-few message; returns
 * false for anything that is not an exact Integer or UpTo[Integer], which is
 * not ours to interpret and goes back to the List path. */
static bool commonest_count_arg(const Expr* na, int64_t* n, bool* upto) {
    if (na->type == EXPR_INTEGER) { *n = na->data.integer; *upto = false; return true; }
    if (na->type == EXPR_FUNCTION && na->data.function.head->type == EXPR_SYMBOL &&
        na->data.function.head->data.symbol.name == SYM_UpTo &&
        na->data.function.arg_count == 1 &&
        na->data.function.args[0]->type == EXPR_INTEGER) {
        *n = na->data.function.args[0]->data.integer;
        *upto = true;
        return true;
    }
    return false;
}

Expr* ndred_commonest(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    if (!nd_word_keyable(a)) return ndarray_delist_and_reeval(res);

    int64_t n = 0;
    bool have_n = (argc == 2), upto = false;
    if (have_n && !commonest_count_arg(res->data.function.args[1], &n, &upto))
        return ndarray_delist_and_reeval(res);

    size_t len = (size_t)a->data.ndarray.dims[0];
    if (len == 0) return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);

    NDWordTally t;
    if (!nd_tally_words(a, &t)) return ndarray_delist_and_reeval(res);
    size_t nuniq = t.nuniq;              /* >= 1: len > 0 and every element counts */

    /* The chosen keys, in first-appearance order -- which is the order the answer
     * comes back in, for both spellings. */
    uint64_t* sel = malloc(sizeof(uint64_t) * nuniq);
    size_t nsel = 0;
    if (!sel) { ndwt_free(&t); return ndarray_delist_and_reeval(res); }

    if (!have_n) {
        /* Everything attaining the top count. t.keys is ALREADY in
         * first-appearance order, so this is a max scan and a filter -- no sort.
         * Worth separating from the branch below: the distinct count is the
         * thing that can be large (10^6 distinct values in a 10^7 buffer is an
         * ordinary case), and sorting all of it to read off a leading run is
         * O(u log u) for an answer that costs O(u). */
        int64_t top = t.cnts[0];
        for (size_t u = 1; u < nuniq; u++) if (t.cnts[u] > top) top = t.cnts[u];
        for (size_t u = 0; u < nuniq; u++) if (t.cnts[u] == top) sel[nsel++] = t.keys[u];
    } else {
        size_t target;
        if (n < 0) {
            target = 0;
        } else if ((uint64_t)n > (uint64_t)nuniq) {
            if (!upto)                   /* UpTo[n] asks for "at most", and is silent */
                printf("Commonest::dstlms: The requested number of elements %" PRId64
                       " is greater than the number of distinct elements %zu."
                       " Only %zu elements will be returned.\n", n, nuniq, nuniq);
            target = nuniq;
        } else {
            target = (size_t)n;
        }
        if (target > 0) {
            /* A genuine ranking is needed here, so: order by count descending
             * and first appearance ascending, keep the leading `target`, then
             * restore first-appearance order among those. Both sorts are the
             * List path's, transcribed, which is why the two agree on ties. */
            CommonWord* w = malloc(sizeof(CommonWord) * nuniq);
            if (!w) { free(sel); ndwt_free(&t); return ndarray_delist_and_reeval(res); }
            for (size_t u = 0; u < nuniq; u++) {
                w[u].key = t.keys[u]; w[u].cnt = t.cnts[u]; w[u].first = u;
            }
            qsort(w, nuniq, sizeof(CommonWord), commonest_by_count);
            qsort(w, target, sizeof(CommonWord), commonest_by_first);
            for (size_t i = 0; i < target; i++) sel[nsel++] = w[i].key;
            free(w);
        }
    }
    ndwt_free(&t);

    if (nsel == 0) { free(sel); return expr_new_function(expr_new_symbol(SYM_List), NULL, 0); }

    /* Every selected element came OUT of the buffer, so they share its dtype and
     * the answer is one head -- unlike Tally, whose {value, count} pairs pack
     * only when the value is itself an integer. Keeping it packed matters
     * because Commonest is usually read by something else. */
    NDType   dt = a->data.ndarray.dtype;
    int64_t  dims[1] = { (int64_t)nsel };
    void*    buf = malloc(ndt_elem_size(dt) * nsel);
    if (!buf) { free(sel); return ndarray_delist_and_reeval(res); }
    if (dt == NDT_INT64) {
        int64_t* iv = (int64_t*)buf;
        for (size_t i = 0; i < nsel; i++) iv[i] = (int64_t)sel[i];
    } else {
        memcpy(buf, sel, sizeof(double) * nsel);   /* the keys ARE the bit patterns */
    }
    free(sel);
    return expr_new_ndarray_like(a, 1, dims, buf, dt);
}

Expr* ndred_moving_average(Expr* res) {
    Expr* a; size_t r;
    if (!nd_moving_window(res, &a, &r)) return ndarray_delist_and_reeval(res);
    NDType dt = a->data.ndarray.dtype;
    const void* buf = a->data.ndarray.data;
    size_t n = (size_t)a->data.ndarray.dims[0];
    size_t L = n - r + 1;

    void* out = malloc(ndt_elem_size(dt) * L);
    if (!out) return ndarray_delist_and_reeval(res);
    double inv = 1.0 / (double)r;
    for (size_t i = 0; i < L; i++) {
        double sr, si;
        nd_sum_strided(buf, dt, i, 1, r, &sr, &si);   /* window [i, i+r) */
        ndt_set(out, i, dt, sr * inv, si * inv);
    }
    int64_t odims[1] = { (int64_t)L };
    return expr_new_ndarray_like(a, 1, odims, out, dt);
}

Expr* ndred_moving_median(Expr* res) {
    Expr* a; size_t r;
    if (!nd_moving_window(res, &a, &r)) return ndarray_delist_and_reeval(res);
    NDType dt = a->data.ndarray.dtype;
    const void* buf = a->data.ndarray.data;
    size_t n = (size_t)a->data.ndarray.dims[0];
    size_t L = n - r + 1;

    NDType odt = nd_real_of(dt);
    void* out = malloc(ndt_elem_size(odt) * L);
    if (!out) return ndarray_delist_and_reeval(res);
    for (size_t i = 0; i < L; i++)
        ndt_set(out, i, odt, nd_median_of(buf, dt, i, 1, r), 0.0);
    int64_t odims[1] = { (int64_t)L };
    return expr_new_ndarray_like(a, 1, odims, out, odt);
}

/* ExponentialMovingAverage[a, alpha]: r[0] = x[0], r[i] = alpha*x[i] +
 * (1-alpha)*r[i-1]. Same length/dtype. Real rank-1 with a real numeric alpha;
 * a complex array or symbolic/complex alpha degrades. */
Expr* ndred_ema(Expr* res) {
    if (res->data.function.arg_count != 2) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    Expr* alpha = res->data.function.args[1];
    NDType dt = a->data.ndarray.dtype;
    if (a->data.ndarray.rank != 1 || ndt_is_complex(dt) || nd_int64_degrade(a))
        return ndarray_delist_and_reeval(res);
    double al;
    if (alpha->type == EXPR_REAL) al = alpha->data.real;
    else if (alpha->type == EXPR_INTEGER) al = (double)alpha->data.integer;
    else return ndarray_delist_and_reeval(res);

    const void* buf = a->data.ndarray.data;
    size_t n = (size_t)a->data.ndarray.dims[0];
    void* out = malloc(ndt_elem_size(dt) * n);
    if (!out) return ndarray_delist_and_reeval(res);
    double prev, im0;
    ndt_get(buf, 0, dt, &prev, &im0);
    ndt_set(out, 0, dt, prev, 0.0);
    for (size_t i = 1; i < n; i++) {
        double x, im;
        ndt_get(buf, i, dt, &x, &im);
        prev = al * x + (1.0 - al) * prev;
        ndt_set(out, i, dt, prev, 0.0);
    }
    return expr_new_ndarray_like(a, 1, a->data.ndarray.dims, out, dt);
}

/* ------------------------------------------------------------------- Total */

/* Sum the leading `m` axes of `a` (1 <= m <= rank). m == rank collapses to a
 * scalar; otherwise the result is a rank-(rank-m) NDArray of the trailing dims.
 * dtype is preserved (a real sum stays real, a complex sum stays complex). */
/* Columnwise (strided) sum: each output column j = sum over its `blocks`
 * summands at stride T. Parallelized over disjoint output columns. */
typedef struct { const void* buf; NDType dt; void* out; size_t T, blocks; double inv; } nd_col_ctx;
static bool nd_total_cols(void* c, size_t lo, size_t hi) {
    const nd_col_ctx* x = (const nd_col_ctx*)c;
    for (size_t j = lo; j < hi; j++) {
        double re, im;
        nd_sum_strided(x->buf, x->dt, j, x->T, x->blocks, &re, &im);
        ndt_set(x->out, j, x->dt, re * x->inv, im * x->inv);
    }
    return true;
}

static Expr* nd_total_leading(const Expr* a, int m) {
    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;
    NDType dt = a->data.ndarray.dtype;
    const void* buf = a->data.ndarray.data;
    bool cplx = ndt_is_complex(dt);

    size_t T = nd_dim_prod(dims, m, rank);       /* trailing (output) size + stride */
    size_t blocks = nd_dim_prod(dims, 0, m);     /* number of summands per output */

    if (m == rank) {                             /* full reduction -> scalar */
        if (dt == NDT_INT64) {
            int64_t s;
            if (!nd_sum_i64(buf, dt, 0, 1, blocks, &s)) return NULL;  /* overflow */
            return expr_new_integer(s);
        }
        double re, im;
        nd_full_sum(buf, dt, blocks, &re, &im);
        return nd_scalar(re, im, cplx);
    }

    void* out = malloc(ndt_elem_size(dt) * T);
    if (!out) return NULL;
    if (dt == NDT_INT64) {
        /* Serial and un-pairwise on purpose: integer addition is exact and
         * associative, so there is nothing to compensate for -- and one
         * overflow anywhere has to abandon the whole result. */
        for (size_t j = 0; j < T; j++) {
            int64_t s;
            if (!nd_sum_i64(buf, dt, j, T, blocks, &s)) { free(out); return NULL; }
            ndt_set_i(out, j, dt, s);
        }
        return expr_new_ndarray_like(a, rank - m, dims + m, out, dt);
    }
    nd_col_ctx c = { buf, dt, out, T, blocks, 1.0 };
    nd_parallel_for(T, nd_total_cols, &c);
    return expr_new_ndarray_like(a, rank - m, dims + m, out, dt); /* adopts out */
}

/* Sum a CONTIGUOUS RANGE of axes [p, q) with p > 0, i.e. Total[a, {n1, n2}]
 * where the summed levels do not start at the top. The surviving axes are
 * dims[0..p) followed by dims[q..rank).
 *
 * With outer = prod(dims[0..p)), mid = prod(dims[p..q)) and
 * inner = prod(dims[q..rank)), a row-major buffer gives
 *
 *     out[o*inner + i] = sum over k < mid of buf[(o*mid + k)*inner + i]
 *
 * so every output element is again a strided sum -- the same shape of work
 * nd_total_leading does, just with a base that walks two indices instead of
 * one. p == 0 is left to nd_total_leading, which is the hot path (plain
 * Total[a]) and does not need the div/mod. */
typedef struct { const void* buf; NDType dt; void* out; size_t inner, mid; } nd_axes_ctx;
static bool nd_total_axes_cols(void* c, size_t lo, size_t hi) {
    const nd_axes_ctx* x = (const nd_axes_ctx*)c;
    for (size_t j = lo; j < hi; j++) {
        size_t o = j / x->inner, i = j % x->inner;
        double re, im;
        nd_sum_strided(x->buf, x->dt, o * x->mid * x->inner + i,
                       x->inner, x->mid, &re, &im);
        ndt_set(x->out, j, x->dt, re, im);
    }
    return true;
}

static Expr* nd_total_axes(const Expr* a, int p, int q) {
    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;
    NDType dt = a->data.ndarray.dtype;
    const void* buf = a->data.ndarray.data;

    size_t outer = nd_dim_prod(dims, 0, p);
    size_t mid   = nd_dim_prod(dims, p, q);
    size_t inner = nd_dim_prod(dims, q, rank);
    size_t T = outer * inner;                    /* output element count */

    int64_t odims[64];
    int orank = 0;
    for (int i = 0; i < p; i++)        odims[orank++] = dims[i];
    for (int i = q; i < rank; i++)     odims[orank++] = dims[i];

    void* out = malloc(ndt_elem_size(dt) * (T ? T : 1));
    if (!out) return NULL;
    if (dt == NDT_INT64) {
        /* Serial and un-pairwise, as in nd_total_leading: integer addition is
         * exact and associative, and one overflow anywhere abandons the whole
         * result rather than silently widening. */
        for (size_t j = 0; j < T; j++) {
            int64_t s;
            if (!nd_sum_i64(buf, dt, (j / inner) * mid * inner + (j % inner),
                            inner, mid, &s)) { free(out); return NULL; }
            ndt_set_i(out, j, dt, s);
        }
        return expr_new_ndarray_like(a, orank, odims, out, dt);
    }
    nd_axes_ctx c = { buf, dt, out, inner, mid };
    nd_parallel_for(T, nd_total_axes_cols, &c);
    return expr_new_ndarray_like(a, orank, odims, out, dt);   /* adopts out */
}

Expr* ndred_total(Expr* res) {
    Expr* a = res->data.function.args[0];
    int rank = a->data.ndarray.rank;
    int m = 1;                                   /* Total[a] sums level 1 */

    if (res->data.function.arg_count == 2) {
        Expr* spec = res->data.function.args[1];
        if (spec->type == EXPR_INTEGER) {
            int64_t n = spec->data.integer;
            if (n < 1 || n > rank) return ndarray_delist_and_reeval(res);
            m = (int)n;
        } else if (spec->type == EXPR_SYMBOL && spec->data.symbol.name == SYM_Infinity) {
            m = rank;                            /* Total[a, Infinity] flattens fully */
        } else if (spec->type == EXPR_FUNCTION && spec->data.function.head &&
                   spec->data.function.head->type == EXPR_SYMBOL &&
                   spec->data.function.head->data.symbol.name == SYM_List) {
            /* Total[a, {n}] sums level n only; Total[a, {n1, n2}] sums levels
             * n1..n2. Both are a contiguous axis range, which the buffer can do
             * directly -- and had to, because these were the ONLY spellings
             * falling through to the List path. Total[m] and Total[m, {1}] are
             * the same value by definition and read 1.36 ms against 258 ms, a
             * 190x penalty for writing the level as a list. Negative and
             * out-of-range levels still take the List path, which is always
             * right. */
            size_t nspec = spec->data.function.arg_count;
            if (nspec < 1 || nspec > 2) return ndarray_delist_and_reeval(res);
            Expr* lo = spec->data.function.args[0];
            if (lo->type != EXPR_INTEGER) return ndarray_delist_and_reeval(res);
            int64_t n1 = lo->data.integer, n2 = n1;
            if (nspec == 2) {
                Expr* hi = spec->data.function.args[1];
                if (hi->type == EXPR_SYMBOL && hi->data.symbol.name == SYM_Infinity)
                    n2 = rank;
                else if (hi->type == EXPR_INTEGER)
                    n2 = hi->data.integer;
                else
                    return ndarray_delist_and_reeval(res);
            }
            if (n1 < 1 || n2 < n1 || n2 > rank) return ndarray_delist_and_reeval(res);
            Expr* r = (n1 == 1) ? nd_total_leading(a, (int)n2)
                                : nd_total_axes(a, (int)n1 - 1, (int)n2);
            return r ? r : ndarray_delist_and_reeval(res);
        } else {
            return ndarray_delist_and_reeval(res);
        }
    } else if (res->data.function.arg_count != 1) {
        return ndarray_delist_and_reeval(res);
    }

    Expr* r = nd_total_leading(a, m);
    return r ? r : ndarray_delist_and_reeval(res);
}

Expr* ndred_total_all(const Expr* a) {
    if (!a || a->type != EXPR_NDARRAY) return NULL;
    return nd_total_leading(a, a->data.ndarray.rank);
}

/* -------------------------------------------------------------------- Mean */

Expr* ndred_mean(Expr* res) {
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;
    NDType dt = a->data.ndarray.dtype;
    const void* buf = a->data.ndarray.data;
    bool cplx = ndt_is_complex(dt);

    size_t blocks = (size_t)dims[0];             /* leading-axis length = divisor */
    if (blocks == 0) return ndarray_delist_and_reeval(res);
    size_t T = nd_dim_prod(dims, 1, rank);
    double inv = 1.0 / (double)blocks;

    if (rank == 1) {                             /* vector -> scalar */
        if (dt == NDT_INT64) {
            /* Mean[Range[10]] is 11/2, not 5.5 -- an exact Rational, reduced. */
            int64_t s;
            if (!nd_sum_i64(buf, dt, 0, 1, blocks, &s))
                return ndarray_delist_and_reeval(res);
            Expr* q = make_rational(s, (int64_t)blocks);
            return q ? q : ndarray_delist_and_reeval(res);
        }
        double re, im;
        nd_full_sum(buf, dt, blocks, &re, &im);
        return nd_scalar(re * inv, im * inv, cplx);
    }
    /* A columnwise integer mean is a vector of Rationals, which no buffer can
     * hold; the List path builds it exactly. */
    if (dt == NDT_INT64) return ndarray_delist_and_reeval(res);

    void* out = malloc(ndt_elem_size(dt) * T);
    if (!out) return ndarray_delist_and_reeval(res);
    nd_col_ctx c = { buf, dt, out, T, blocks, inv };  /* inv scales sum -> mean */
    nd_parallel_for(T, nd_total_cols, &c);
    return expr_new_ndarray_like(a, rank - 1, dims + 1, out, dt);
}

/* ------------------------------------------------- Variance / Std / RMS */

/* One column's second moment at (base, stride) over `blocks` summands, serial
 * (used per-column in the columnwise path). mode: 0=Variance, 1=Std, 2=RMS. */
static double nd_moment_at(const void* buf, NDType dt, size_t base, size_t stride,
                           size_t blocks, int mode) {
    if (mode == 2) {
        double ss = nd_sumsq_strided(buf, dt, base, stride, blocks, 0.0, 0.0);
        return sqrt(ss / (double)blocks);
    }
    double sr, si;
    nd_sum_strided(buf, dt, base, stride, blocks, &sr, &si);
    double mr = sr / (double)blocks, mi = si / (double)blocks;
    double ss = nd_sumsq_strided(buf, dt, base, stride, blocks, mr, mi);
    double var_ = ss / (double)(blocks - 1);
    return (mode == 1) ? sqrt(var_) : var_;
}

typedef struct {
    const void* buf; NDType dt; void* out; NDType odt; size_t T, blocks; int mode;
} nd_mom_ctx;
static bool nd_moment_cols(void* c, size_t lo, size_t hi) {
    const nd_mom_ctx* x = (const nd_mom_ctx*)c;
    for (size_t j = lo; j < hi; j++)
        ndt_set(x->out, j, x->odt, nd_moment_at(x->buf, x->dt, j, x->T, x->blocks, x->mode), 0.0);
    return true;
}

/* Columnwise (leading-axis) second-moment reduction. `mode`:
 *   0 = Variance (Sum|x-mean|^2 / (n-1)),
 *   1 = StandardDeviation (sqrt of Variance),
 *   2 = RootMeanSquare (sqrt(Sum|x|^2 / n)).
 * Always produces a REAL result (a complex column yields a real spread), so the
 * output dtype is the real dtype of matching component width. The rank-1 case
 * runs the two passes (mean, then Σ|x−μ|²) through the threaded full reducers. */
static Expr* nd_moment_leading(Expr* res, int mode) {
    Expr* a = res->data.function.args[0];
    /* Variance[Range[10]] is 55/6 and RootMeanSquare[{1,2}] is Sqrt[5/2]; both
     * exact answers are outside a machine buffer. */
    if (nd_int64_degrade(a)) return ndarray_delist_and_reeval(res);
    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;
    NDType dt = a->data.ndarray.dtype;
    const void* buf = a->data.ndarray.data;

    size_t blocks = (size_t)dims[0];
    /* Variance/Std need n >= 2 (n-1 divisor); RMS only n >= 1. */
    if (blocks < (mode == 2 ? 1u : 2u)) return ndarray_delist_and_reeval(res);
    size_t T = nd_dim_prod(dims, 1, rank);
    NDType odt = nd_real_of(dt);

    if (rank == 1) {                             /* vector -> scalar (threaded) */
        double v;
        if (mode == 2) {
            double ss = nd_full_sumsq(buf, dt, blocks, 0.0, 0.0);
            v = sqrt(ss / (double)blocks);
        } else {
            double sr, si;
            nd_full_sum(buf, dt, blocks, &sr, &si);
            double mr = sr / (double)blocks, mi = si / (double)blocks;
            double ss = nd_full_sumsq(buf, dt, blocks, mr, mi);
            double var_ = ss / (double)(blocks - 1);
            v = (mode == 1) ? sqrt(var_) : var_;
        }
        return expr_new_real(v);
    }

    void* out = malloc(ndt_elem_size(odt) * T);
    if (!out) return ndarray_delist_and_reeval(res);
    nd_mom_ctx c = { buf, dt, out, odt, T, blocks, mode };
    nd_parallel_for(T, nd_moment_cols, &c);
    return expr_new_ndarray_like(a, rank - 1, dims + 1, out, odt);
}

Expr* ndred_variance(Expr* res) {
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    return nd_moment_leading(res, 0);
}
Expr* ndred_std(Expr* res) {
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    return nd_moment_leading(res, 1);
}
Expr* ndred_rms(Expr* res) {
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    return nd_moment_leading(res, 2);
}

/* ------------------------------------------------------------- CentralMoment */

/* Integer power by squaring on a double (exponent >= 0). */
static double nd_ipow(double b, int64_t e) {
    double r = 1.0;
    while (e > 0) { if (e & 1) r *= b; b *= b; e >>= 1; }
    return r;
}

/* Pairwise Sum[(x - mean)^r] over the strided range (real component only; the
 * leading kernel degrades int64/complex before ever calling this). The pairwise
 * split mirrors nd_sumsq_strided so the accumulation order — hence rounding —
 * matches the other reductions across the packed / visible / List surfaces. */
static double nd_summoment_strided(const void* buf, NDType dt, size_t base,
                                   size_t stride, size_t count, double mean, int64_t r) {
    if (count <= ND_PAIRWISE_BLOCK) {
        double s = 0.0;
        if (dt == NDT_FLOAT64 && stride == 1) {          /* contiguous real */
            const double* p = (const double*)buf + base;
            for (size_t i = 0; i < count; i++) s += nd_ipow(p[i] - mean, r);
            return s;
        }
        for (size_t i = 0; i < count; i++) {
            double re, im;
            ndt_get(buf, base + i * stride, dt, &re, &im);
            s += nd_ipow(re - mean, r);
        }
        return s;
    }
    size_t half = count / 2;
    return nd_summoment_strided(buf, dt, base, stride, half, mean, r) +
           nd_summoment_strided(buf, dt, base + half * stride, stride, count - half, mean, r);
}

typedef struct { const void* buf; NDType dt; double mean; int64_t r; } nd_cmfull_ctx;
static void nd_cm_reduce(void* c, size_t lo, size_t hi, double* slot) {
    const nd_cmfull_ctx* x = (const nd_cmfull_ctx*)c;
    slot[0] = nd_summoment_strided(x->buf, x->dt, lo, 1, hi - lo, x->mean, x->r);
}
/* Sum[(x - mean)^r] over `n` contiguous elements, threaded (rank-1 path). */
static double nd_full_moment(const void* buf, NDType dt, size_t n, double mean, int64_t r) {
    nd_cmfull_ctx c = { buf, dt, mean, r };
    double slots[NDARRAY_MAX_THREADS];
    int k = nd_parallel_reduce(n, nd_cm_reduce, &c, 1, slots);
    double s = 0.0;
    for (int t = 0; t < k; t++) s += slots[t];
    return s;
}

/* One column's r-th central moment at (base, stride) over `blocks` summands. */
static double nd_cm_at(const void* buf, NDType dt, size_t base, size_t stride,
                       size_t blocks, int64_t r) {
    double sr, si;
    nd_sum_strided(buf, dt, base, stride, blocks, &sr, &si);
    double mean = sr / (double)blocks;
    return nd_summoment_strided(buf, dt, base, stride, blocks, mean, r) / (double)blocks;
}

typedef struct {
    const void* buf; NDType dt; void* out; NDType odt; size_t T, blocks; int64_t r;
} nd_cm_ctx;
static bool nd_cm_cols(void* c, size_t lo, size_t hi) {
    const nd_cm_ctx* x = (const nd_cm_ctx*)c;
    for (size_t j = lo; j < hi; j++)
        ndt_set(x->out, j, x->odt, nd_cm_at(x->buf, x->dt, j, x->T, x->blocks, x->r), 0.0);
    return true;
}

/* CentralMoment[a, r]: the r-th moment about the mean. Two passes (mean, then
 * Sum[(x-mean)^r]) over the leading axis, dividing by n — like Variance but /n
 * (not n-1), power r (not a square), no Conjugate, n >= 1. Always a REAL result.
 * Real dtypes and a non-negative integer r only: an int64 buffer (the exact
 * answer is a Rational), a complex buffer ((x-mu)^r is complex, no real slot for
 * it), and a list-valued / non-integer / negative r all degrade to the List path,
 * which answers exactly or symbolically. */
static Expr* nd_central_moment_leading(Expr* res) {
    if (res->data.function.arg_count != 2) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    Expr* rexpr = res->data.function.args[1];
    if (rexpr->type != EXPR_INTEGER || rexpr->data.integer < 0)
        return ndarray_delist_and_reeval(res);
    int64_t r = rexpr->data.integer;
    if (nd_int64_degrade(a) || ndt_is_complex(a->data.ndarray.dtype))
        return ndarray_delist_and_reeval(res);

    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;
    NDType dt = a->data.ndarray.dtype;
    const void* buf = a->data.ndarray.data;

    size_t blocks = (size_t)dims[0];
    if (blocks < 1u) return ndarray_delist_and_reeval(res);   /* Mean of empty is undefined */
    size_t T = nd_dim_prod(dims, 1, rank);
    NDType odt = nd_real_of(dt);

    if (rank == 1) {                             /* vector -> scalar (threaded) */
        double sr, si;
        nd_full_sum(buf, dt, blocks, &sr, &si);
        double mean = sr / (double)blocks;
        double ss = nd_full_moment(buf, dt, blocks, mean, r);
        return expr_new_real(ss / (double)blocks);
    }

    void* out = malloc(ndt_elem_size(odt) * T);
    if (!out) return ndarray_delist_and_reeval(res);
    nd_cm_ctx c = { buf, dt, out, odt, T, blocks, r };
    nd_parallel_for(T, nd_cm_cols, &c);
    return expr_new_ndarray_like(a, rank - 1, dims + 1, out, odt);
}

Expr* ndred_central_moment(Expr* res) {
    return nd_central_moment_leading(res);
}

/* ------------------------------------------------- Skewness / Kurtosis */

/* A standardized central moment m_p / m_2^(p/2): Skewness is p=3
 * (m_3 / m_2^(3/2)), Kurtosis is p=4 (m_4 / m_2^2). Three passes over the
 * leading axis (mean, then m_2 and m_p about it), reusing the CentralMoment
 * summation helpers. Real dtypes only; int64 / complex degrade to the exact
 * List path (an integer sample's standardized moment is a radical, e.g.
 * Skewness[{1,2,3,10}] = 18 Sqrt[2]/25 — no machine slot holds it). */
static double nd_stdmoment_at(const void* buf, NDType dt, size_t base, size_t stride,
                              size_t blocks, int p, double half_p) {
    double sr, si;
    nd_sum_strided(buf, dt, base, stride, blocks, &sr, &si);
    double mean = sr / (double)blocks;
    double m2 = nd_summoment_strided(buf, dt, base, stride, blocks, mean, 2) / (double)blocks;
    double mp = nd_summoment_strided(buf, dt, base, stride, blocks, mean, p) / (double)blocks;
    return mp / pow(m2, half_p);
}

typedef struct {
    const void* buf; NDType dt; void* out; NDType odt; size_t T, blocks; int p; double half_p;
} nd_sm_ctx;
static bool nd_sm_cols(void* c, size_t lo, size_t hi) {
    const nd_sm_ctx* x = (const nd_sm_ctx*)c;
    for (size_t j = lo; j < hi; j++)
        ndt_set(x->out, j, x->odt,
                nd_stdmoment_at(x->buf, x->dt, j, x->T, x->blocks, x->p, x->half_p), 0.0);
    return true;
}

static Expr* nd_stdmoment_leading(Expr* res, int p) {
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    if (nd_int64_degrade(a) || ndt_is_complex(a->data.ndarray.dtype))
        return ndarray_delist_and_reeval(res);

    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;
    NDType dt = a->data.ndarray.dtype;
    const void* buf = a->data.ndarray.data;

    size_t blocks = (size_t)dims[0];
    if (blocks < 1u) return ndarray_delist_and_reeval(res);
    size_t T = nd_dim_prod(dims, 1, rank);
    NDType odt = nd_real_of(dt);
    double half_p = (double)p / 2.0;

    if (rank == 1) {
        double sr, si;
        nd_full_sum(buf, dt, blocks, &sr, &si);
        double mean = sr / (double)blocks;
        double m2 = nd_full_moment(buf, dt, blocks, mean, 2) / (double)blocks;
        double mp = nd_full_moment(buf, dt, blocks, mean, p) / (double)blocks;
        return expr_new_real(mp / pow(m2, half_p));
    }

    void* out = malloc(ndt_elem_size(odt) * T);
    if (!out) return ndarray_delist_and_reeval(res);
    nd_sm_ctx c = { buf, dt, out, odt, T, blocks, p, half_p };
    nd_parallel_for(T, nd_sm_cols, &c);
    return expr_new_ndarray_like(a, rank - 1, dims + 1, out, odt);
}

Expr* ndred_skewness(Expr* res) { return nd_stdmoment_leading(res, 3); }
Expr* ndred_kurtosis(Expr* res) { return nd_stdmoment_leading(res, 4); }

/* --------------------------------------------------------------- Max / Min */

/* Extreme of the flat range [lo, hi): *best gets the max/min, *sawnan is set if
 * any element was NaN. float64 uses four independent running extrema so the
 * ternary compiles to a vector max/min reduction (NaN flagged branchlessly via
 * v != v); other dtypes go through ndt_get. */
static void nd_range_extreme(const void* buf, NDType dt, size_t lo, size_t hi,
                             bool want_max, double* best, int* sawnan) {
    double b = want_max ? -INFINITY : INFINITY;
    int nan = 0;
    if (dt == NDT_FLOAT64) {
        const double* p = (const double*)buf;
        double b0=b, b1=b, b2=b, b3=b;
        size_t k = lo;
        for (; k + 4 <= hi; k += 4) {
            double v0=p[k], v1=p[k+1], v2=p[k+2], v3=p[k+3];
            nan |= (v0!=v0)|(v1!=v1)|(v2!=v2)|(v3!=v3);
            if (want_max) { b0=v0>b0?v0:b0; b1=v1>b1?v1:b1; b2=v2>b2?v2:b2; b3=v3>b3?v3:b3; }
            else          { b0=v0<b0?v0:b0; b1=v1<b1?v1:b1; b2=v2<b2?v2:b2; b3=v3<b3?v3:b3; }
        }
        for (; k < hi; k++) {
            double v = p[k]; nan |= (v != v);
            b = want_max ? (v > b ? v : b) : (v < b ? v : b);
        }
        double u = want_max ? (b0>b1?b0:b1) : (b0<b1?b0:b1);
        double w = want_max ? (b2>b3?b2:b3) : (b2<b3?b2:b3);
        double m = want_max ? (u>w?u:w) : (u<w?u:w);
        b = want_max ? (m>b?m:b) : (m<b?m:b);
    } else {
        for (size_t k = lo; k < hi; k++) {
            double r, im; ndt_get(buf, k, dt, &r, &im);
            nan |= (r != r);
            b = want_max ? (r > b ? r : b) : (r < b ? r : b);
        }
    }
    *best = b; *sawnan = nan;
}

typedef struct { const void* buf; NDType dt; bool want_max; } nd_ext_ctx;
static void nd_ext_reduce(void* c, size_t lo, size_t hi, double* slot) {
    const nd_ext_ctx* x = (const nd_ext_ctx*)c;
    double b; int nan;
    nd_range_extreme(x->buf, x->dt, lo, hi, x->want_max, &b, &nan);
    slot[0] = b; slot[1] = (double)nan;         /* slot: {extreme, nan-flag} */
}

/* Max (want_max=true) / Min over every element (full flatten) -> real scalar.
 * Complex has no order and a NaN element would make the List result symbolic,
 * so both degrade. Threaded over the flat range. */
static Expr* nd_extreme(Expr* res, bool want_max) {
    Expr* a = res->data.function.args[0];
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    NDType dt = a->data.ndarray.dtype;
    if (ndt_is_complex(dt)) return ndarray_delist_and_reeval(res);
    size_t sz = ndarray_size(a);

    if (dt == NDT_INT64) {
        /* Exact, and the head matters: Max[Range[300]] is the Integer 300. */
        if (sz == 0) return ndarray_delist_and_reeval(res);
        return expr_new_integer(
            nd_extreme_i64((const int64_t*)a->data.ndarray.data, sz, want_max));
    }

    nd_ext_ctx c = { a->data.ndarray.data, dt, want_max };
    double slots[NDARRAY_MAX_THREADS * 2];
    int k = nd_parallel_reduce(sz, nd_ext_reduce, &c, 2, slots);
    double best = want_max ? -INFINITY : INFINITY;
    int sawnan = 0;
    for (int t = 0; t < k; t++) {
        double b = slots[2 * t];
        best = want_max ? (b > best ? b : best) : (b < best ? b : best);
        sawnan |= (slots[2 * t + 1] != 0.0);
    }
    if (sawnan) return ndarray_delist_and_reeval(res);
    return expr_new_real(best);
}

Expr* ndred_max(Expr* res) { return nd_extreme(res, true); }
Expr* ndred_min(Expr* res) { return nd_extreme(res, false); }

/* ------------------------------------------------------------- Accumulate */

/* Prefix sum along the leading axis, same shape/dtype:
 *   out[b, j] = Sum_{b' <= b} in[b', j],  j over the trailing block. */
/* ------------------------------------------------------------------- scans */

/* Recognise the scan operator — see ndreduce.h. */
bool ndred_scan_op_for(const Expr* f, NDScanOp* op) {
    const char* nm = NULL;
    if (f->type == EXPR_SYMBOL) {
        nm = f->data.symbol.name;
    } else if (f->type == EXPR_FUNCTION &&
               f->data.function.head->type == EXPR_SYMBOL &&
               f->data.function.head->data.symbol.name == SYM_Function &&
               f->data.function.arg_count == 1) {
        /* Function[h[Slot[1], Slot[2]]] — the `h[#1, #2] &` spelling. */
        const Expr* b = f->data.function.args[0];
        if (b->type != EXPR_FUNCTION || b->data.function.arg_count != 2 ||
            b->data.function.head->type != EXPR_SYMBOL) return false;
        for (int i = 0; i < 2; i++) {
            const Expr* s = b->data.function.args[i];
            if (s->type != EXPR_FUNCTION || s->data.function.arg_count != 1 ||
                s->data.function.head->type != EXPR_SYMBOL ||
                s->data.function.head->data.symbol.name != SYM_Slot ||
                s->data.function.args[0]->type != EXPR_INTEGER ||
                s->data.function.args[0]->data.integer != i + 1) return false;
        }
        nm = b->data.function.head->data.symbol.name;
    } else {
        return false;
    }
    if (nm == SYM_Plus)  { *op = ND_SCAN_PLUS;  return true; }
    if (nm == SYM_Times) { *op = ND_SCAN_TIMES; return true; }
    if (nm == SYM_Max)   { *op = ND_SCAN_MAX;   return true; }
    if (nm == SYM_Min)   { *op = ND_SCAN_MIN;   return true; }
    return false;
}

Expr* ndred_scan(const Expr* a, NDScanOp op, const Expr* seed, bool as_list) {
    if (!a || !is_ndarray(a) || a->data.ndarray.rank != 1) return NULL;
    NDType dt = a->data.ndarray.dtype;
    if (ndt_is_complex(dt)) return NULL;
    size_t n = (size_t)a->data.ndarray.dims[0];

    /* Uniformity of the ANSWER decides whether a buffer can hold it. Max and Min
     * hand back one of their arguments unchanged, and Plus/Times over exact
     * inputs stay exact, so seed and elements must agree in exactness — a Real
     * seed over an int64 buffer would give Reals only where the seed won. */
    bool want_exact = (dt == NDT_INT64);
    if (seed) {
        if (want_exact) { if (seed->type != EXPR_INTEGER) return NULL; }
        else            { if (seed->type != EXPR_REAL)    return NULL; }
    } else {
        if (n == 0) return NULL;      /* FoldList[f, {}] -> {}; List path */
    }
    /* Seedless (2-arg) form scans elements 1.. from element 0. */
    size_t start = seed ? 0 : 1;
    size_t outn = seed ? n + 1 : n;
    if (outn == 0) return NULL;

    if (want_exact) {
        const int64_t* in = (const int64_t*)a->data.ndarray.data;
        int64_t acc = seed ? seed->data.integer : in[0];
        int64_t* out = as_list ? malloc(sizeof(int64_t) * outn) : NULL;
        if (as_list && !out) return NULL;
        if (as_list) out[0] = acc;
        for (size_t i = start, k = 1; i < n; i++, k++) {
            int64_t x = in[i];
            switch (op) {
            /* Overflow abandons the whole scan: the List answer promotes to a
             * bigint and no int64 buffer can hold it. */
            case ND_SCAN_PLUS:  if (ci_add_i64(acc, x, &acc)) { free(out); return NULL; } break;
            case ND_SCAN_TIMES: if (ci_mul_i64(acc, x, &acc)) { free(out); return NULL; } break;
            case ND_SCAN_MAX:   if (x > acc) acc = x; break;
            case ND_SCAN_MIN:   if (x < acc) acc = x; break;
            }
            if (as_list) out[k] = acc;
        }
        if (!as_list) return expr_new_integer(acc);
        int64_t dims[1] = { (int64_t)outn };
        return expr_new_ndarray_like(a, 1, dims, out, dt);
    }

    /* Real dtypes. float64 reads the buffer directly (the common case and the
     * one that vectorises); float32 goes through the ndt_get/ndt_set choke
     * point so the stored precision is respected. */
    double* out = as_list ? malloc(sizeof(double) * outn) : NULL;
    if (as_list && !out) return NULL;
    double acc;
    if (seed) acc = seed->data.real;
    else { double im; ndt_get(a->data.ndarray.data, 0, dt, &acc, &im); }
    if (as_list) out[0] = acc;

    if (dt == NDT_FLOAT64) {
        const double* in = (const double*)a->data.ndarray.data;
        for (size_t i = start, k = 1; i < n; i++, k++) {
            double x = in[i];
            switch (op) {
            case ND_SCAN_PLUS:  acc += x; break;
            case ND_SCAN_TIMES: acc *= x; break;
            case ND_SCAN_MAX:   if (x > acc) acc = x; break;
            case ND_SCAN_MIN:   if (x < acc) acc = x; break;
            }
            if (as_list) out[k] = acc;
        }
    } else {
        for (size_t i = start, k = 1; i < n; i++, k++) {
            double x, im;
            ndt_get(a->data.ndarray.data, i, dt, &x, &im);
            switch (op) {
            case ND_SCAN_PLUS:  acc += x; break;
            case ND_SCAN_TIMES: acc *= x; break;
            case ND_SCAN_MAX:   if (x > acc) acc = x; break;
            case ND_SCAN_MIN:   if (x < acc) acc = x; break;
            }
            if (as_list) out[k] = acc;
        }
    }
    if (!as_list) { free(out); return expr_new_real(acc); }
    /* A float32 source keeps its dtype, so the partials are stored at the
     * source's precision — same as every other structural path here. */
    int64_t dims[1] = { (int64_t)outn };
    if (dt == NDT_FLOAT64)
        return expr_new_ndarray_like(a, 1, dims, out, dt);
    void* nb = malloc(ndt_elem_size(dt) * outn);
    if (!nb) { free(out); return NULL; }
    for (size_t i = 0; i < outn; i++) ndt_set(nb, i, dt, out[i], 0.0);
    free(out);
    return expr_new_ndarray_like(a, 1, dims, nb, dt);
}

Expr* ndred_accumulate(Expr* res) {
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;
    NDType dt = a->data.ndarray.dtype;
    const void* buf = a->data.ndarray.data;

    size_t blocks = (size_t)dims[0];
    size_t T = nd_dim_prod(dims, 1, rank);
    size_t sz = blocks * T;

    void* out = malloc(ndt_elem_size(dt) * sz);
    if (!out) return ndarray_delist_and_reeval(res);
    if (dt == NDT_INT64) {
        const int64_t* in = (const int64_t*)buf;
        int64_t* o = (int64_t*)out;
        for (size_t j = 0; j < T; j++) o[j] = in[j];
        for (size_t b = 1; b < blocks; b++)
            for (size_t j = 0; j < T; j++)
                if (ci_add_i64(o[(b - 1) * T + j], in[b * T + j], &o[b * T + j])) {
                    free(out);
                    return ndarray_delist_and_reeval(res);   /* bigint answer */
                }
        return expr_new_ndarray_like(a, rank, dims, out, dt);
    }
    if (dt == NDT_FLOAT64) {
        /* Direct on the buffer: the generic arm below pays two indirect calls
         * per element through the ndt_get/ndt_set choke point, which is what
         * kept Accumulate at 4x np.cumsum on data that is one subtraction away
         * from a plain double loop. Same values, same order of summation. */
        const double* in = (const double*)buf;
        double* o = (double*)out;
        if (T == 1) {
            /* Rank 1 -- the overwhelmingly common Accumulate. Written flat
             * rather than as the T-inner loop below, because T is a RUNTIME
             * value: at T == 1 the compiler still emits the whole loop prologue,
             * the index multiplies and the trip test around a single add, which
             * is ~10 operations per element instead of one. That, not the serial
             * dependency, was the 4x against np.cumsum. */
            double acc = in[0];
            o[0] = acc;
            for (size_t b = 1; b < blocks; b++) { acc += in[b]; o[b] = acc; }
            return expr_new_ndarray_like(a, rank, dims, out, dt);
        }
        for (size_t j = 0; j < T; j++) o[j] = in[j];
        for (size_t b = 1; b < blocks; b++)
            for (size_t j = 0; j < T; j++)
                o[b * T + j] = o[(b - 1) * T + j] + in[b * T + j];
        return expr_new_ndarray_like(a, rank, dims, out, dt);
    }
    /* First block copies through; each later block adds the running total. */
    for (size_t j = 0; j < T; j++) {
        double r, im;
        ndt_get(buf, j, dt, &r, &im);
        ndt_set(out, j, dt, r, im);
    }
    for (size_t b = 1; b < blocks; b++) {
        for (size_t j = 0; j < T; j++) {
            double pr, pi, cr, ci;
            ndt_get(out, (b - 1) * T + j, dt, &pr, &pi);
            ndt_get(buf, b * T + j, dt, &cr, &ci);
            ndt_set(out, b * T + j, dt, pr + cr, pi + ci);
        }
    }
    return expr_new_ndarray_like(a, rank, dims, out, dt);
}
