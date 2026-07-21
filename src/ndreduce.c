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
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

/* Below this many summands a strided range is summed by a plain loop; above it
 * the range is split in two and each half summed recursively (pairwise). */
#define ND_PAIRWISE_BLOCK ((size_t)128)

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

/* Below this, the introsort's lower constant factor wins; above it an 8-pass
 * LSD radix sort (memory-bound, no comparisons) is far faster on machine
 * doubles. */
#define ND_RADIX_MIN ((size_t)2048)

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

    if (rank == 1)
        return expr_new_real(nd_median_of(buf, dt, 0, 1, (size_t)dims[0]));

    /* Matrix: columnwise median -> rank-1 vector of length ncols. */
    size_t nrows = (size_t)dims[0], ncols = (size_t)dims[1];
    NDType odt = nd_real_of(dt);
    void* out = malloc(ndt_elem_size(odt) * ncols);
    if (!out) return ndarray_delist_and_reeval(res);
    for (size_t j = 0; j < ncols; j++)
        ndt_set(out, j, odt, nd_median_of(buf, dt, j, ncols, nrows), 0.0);
    int64_t odims[1] = { (int64_t)ncols };
    return expr_new_ndarray(1, odims, out, odt);
}

/* --------------------------------------------------------------- Quartiles */

/* Mathematica's default Quartiles parameters {{1/2, 0}, {0, 1}}: for the k-th
 * quantile q, h = 1/2 + n*q, clamped, then linear interpolation between the
 * bracketing order statistics of the ascending-sorted data. */
Expr* ndred_quartiles(Expr* res) {
    /* Only the default (1-arg) rank-1 real case is fast-pathed; a custom
     * parameter matrix or higher rank degrades to the exact List method. */
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    if (a->data.ndarray.rank != 1 || ndt_is_complex(a->data.ndarray.dtype))
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
    return expr_new_ndarray(1, odims, out, odt);
}

/* ---------------------------------------------------- Moving statistics */

/* Shared preamble for MovingAverage/MovingMedian: require a 2-arg call on a
 * real rank-1 array with a positive integer window r in [1, n]. Returns r via
 * *r_out and the array via *a_out, or false (caller degrades). */
static bool nd_moving_window(Expr* res, Expr** a_out, size_t* r_out) {
    if (res->data.function.arg_count != 2) return false;
    Expr* a = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];
    if (a->data.ndarray.rank != 1 || ndt_is_complex(a->data.ndarray.dtype))
        return false;
    if (spec->type != EXPR_INTEGER) return false;   /* weight list / bignum -> List path */
    int64_t r = spec->data.integer;
    size_t n = (size_t)a->data.ndarray.dims[0];
    if (r < 1 || (size_t)r > n) return false;
    *a_out = a; *r_out = (size_t)r;
    return true;
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
    return expr_new_ndarray(1, odims, out, dt);
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
    return expr_new_ndarray(1, odims, out, odt);
}

/* ExponentialMovingAverage[a, alpha]: r[0] = x[0], r[i] = alpha*x[i] +
 * (1-alpha)*r[i-1]. Same length/dtype. Real rank-1 with a real numeric alpha;
 * a complex array or symbolic/complex alpha degrades. */
Expr* ndred_ema(Expr* res) {
    if (res->data.function.arg_count != 2) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    Expr* alpha = res->data.function.args[1];
    NDType dt = a->data.ndarray.dtype;
    if (a->data.ndarray.rank != 1 || ndt_is_complex(dt))
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
    return expr_new_ndarray(1, a->data.ndarray.dims, out, dt);
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
        double re, im;
        nd_full_sum(buf, dt, blocks, &re, &im);
        return nd_scalar(re, im, cplx);
    }

    void* out = malloc(ndt_elem_size(dt) * T);
    if (!out) return NULL;
    nd_col_ctx c = { buf, dt, out, T, blocks, 1.0 };
    nd_parallel_for(T, nd_total_cols, &c);
    return expr_new_ndarray(rank - m, dims + m, out, dt); /* adopts out */
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
        } else {
            return ndarray_delist_and_reeval(res); /* {k}, {n1,n2}, ... -> List path */
        }
    } else if (res->data.function.arg_count != 1) {
        return ndarray_delist_and_reeval(res);
    }

    Expr* r = nd_total_leading(a, m);
    return r ? r : ndarray_delist_and_reeval(res);
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
        double re, im;
        nd_full_sum(buf, dt, blocks, &re, &im);
        return nd_scalar(re * inv, im * inv, cplx);
    }

    void* out = malloc(ndt_elem_size(dt) * T);
    if (!out) return ndarray_delist_and_reeval(res);
    nd_col_ctx c = { buf, dt, out, T, blocks, inv };  /* inv scales sum -> mean */
    nd_parallel_for(T, nd_total_cols, &c);
    return expr_new_ndarray(rank - 1, dims + 1, out, dt);
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
    return expr_new_ndarray(rank - 1, dims + 1, out, odt);
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
    return expr_new_ndarray(rank, dims, out, dt);
}
