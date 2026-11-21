/* stats_ll.c -- the descriptive statistics over LinkedList data.
 *
 * A LinkedList (src/linkedlist.h) is computed on directly, never converted to
 * a List: H[ll] answers exactly what H[Normal[ll]] answers. Every statistics
 * builtin calls stats_ll_arm() before its List code; it takes three shapes:
 *
 *   1. A machine-real (or machine-complex) LinkedList vector. Normal[ll] of a
 *      typed list is a PACKED List, so the answer must be the packed kernel's
 *      to the bit. The reductions (Mean, Variance, StandardDeviation,
 *      RootMeanSquare, Moment, CentralMoment, Skewness, Kurtosis) walk the
 *      chunks through ll_pairwise_reduce with ndreduce.c's own pairwise leaves,
 *      so the summation tree is the packed one over logical indices; the order
 *      statistics and windows (Median, Quartiles, MovingAverage, MovingMedian)
 *      need random access and run ndreduce.c's raw kernels over an ll_gather
 *      scratch buffer; ExponentialMovingAverage streams the cursor. Whatever
 *      the packed kernel declines (an int64 list, a complex Moment, a result
 *      outside the machine range, ...) is declined here too, and the builtin's
 *      element path answers, as the List path answers the packed decline.
 *
 *   2. A matrix with a LinkedList in it: a LinkedList of rows, or a List whose
 *      rows include a LinkedList, all of one length. The statistic is
 *      columnwise. When every row is a machine-real LinkedList the rows are
 *      gathered row-major into a raw buffer and ndreduce.c's strided column
 *      kernels run on it -- the code the packed rank-2 array runs. Otherwise
 *      each column is built as a fresh LinkedList (slots copied, typed stays
 *      typed) and the head is applied to it.
 *
 *   3. Anything else: not handled; a LinkedList vector continues into the
 *      builtin's element path, which reads it through LLElems (linkedlist.h).
 *
 * Result container: a vector derived element for element from the argument
 * (MovingAverage, ...) is a LinkedList; a columnwise statistic is in the image
 * of the rows (LinkedList rows, LinkedList answer); a scalar stays a scalar and
 * Quartiles stays the fixed-shape List of three.
 */

#include "stats.h"
#include "stats_common.h"
#include "linkedlist.h"
#include "ndreduce.h"
#include "ndarray.h"
#include "eval.h"
#include "sym_names.h"
#include "arithmetic.h"   /* make_rational */
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------
 *  Pairwise leaves (ndreduce.c's, over one contiguous run)
 * ---------------------------------------------------------------------- */

static void leaf_sum(const double* p, size_t c, void* ctx, double* out) {
    if (!*(const bool*)ctx) { out[0] = nd_sum_block_f64(p, c); out[1] = 0.0; }
    else nd_sum_strided(p, NDT_COMPLEX64, 0, 1, c, &out[0], &out[1]);
}

typedef struct { NDType dt; double mr, mi; } SqCtx;
static void leaf_sumsq(const double* p, size_t c, void* ctx, double* out) {
    const SqCtx* x = (const SqCtx*)ctx;
    out[0] = nd_sumsq_strided(p, x->dt, 0, 1, c, x->mr, x->mi);
}

typedef struct { double mean; int64_t r; } MomCtx;
static void leaf_moment(const double* p, size_t c, void* ctx, double* out) {
    const MomCtx* x = (const MomCtx*)ctx;
    out[0] = nd_summoment_strided(p, NDT_FLOAT64, 0, 1, c, x->mean, x->r);
}

static void ll_full_sum(const LinkedList* ll, double* re, double* im) {
    bool cplx = ll_kind(ll) == LL_CPLX;
    double o[2];
    ll_pairwise_reduce(ll, leaf_sum, &cplx, 2, o);
    *re = o[0]; *im = o[1];
}

static double ll_full_sumsq(const LinkedList* ll, double mr, double mi) {
    SqCtx c = { ll_kind(ll) == LL_CPLX ? NDT_COMPLEX64 : NDT_FLOAT64, mr, mi };
    double o;
    ll_pairwise_reduce(ll, leaf_sumsq, &c, 1, &o);
    return o;
}

static double ll_full_moment(const LinkedList* ll, double mean, int64_t r) {
    MomCtx c = { mean, r };
    double o;
    ll_pairwise_reduce(ll, leaf_moment, &c, 1, &o);
    return o;
}

/* ------------------------------------------------------------------------
 *  Shape 1: a machine vector
 * ---------------------------------------------------------------------- */

static Expr* scalar_of(double re, double im, bool cplx) {
    if (!cplx) return expr_new_real(re);
    Expr* a[2] = { expr_new_real(re), expr_new_real(im) };
    return expr_new_function(expr_new_symbol(SYM_Complex), a, 2);
}

/* A non-negative machine-integer order argument. */
static bool order_arg(const Expr* e, int64_t* r) {
    if (e->type != EXPR_INTEGER || e->data.integer < 0) return false;
    *r = e->data.integer;
    return true;
}

/* The elements of a REAL list as a malloc'd double buffer (the random-access
 * scratch of the order statistics). */
static double* gather_real(const LinkedList* ll) {
    int64_t n = ll_length(ll);
    double* s = malloc(sizeof(double) * (size_t)(n ? n : 1));
    if (s) ll_gather(ll, 0, n, s);
    return s;
}

/* A REAL LinkedList of the n doubles in buf (copied). NULL when `guard` is set
 * and one is outside the machine range: the packed kernel's nd_guard declines
 * there too. */
static Expr* real_ll(const double* buf, size_t n, bool guard) {
    if (guard)
        for (size_t i = 0; i < n; i++)
            if (!stats_machine_ok(buf[i])) return NULL;
    return expr_new_linkedlist(ll_from_buffer(LL_REAL, buf, (int64_t)n));
}

static Expr* machine_vector(Expr* res, StatsOp op, const LinkedList* ll) {
    LLKind k = ll_kind(ll);
    if (k != LL_REAL && k != LL_CPLX) return NULL;
    /* A list too short to pack has a plain Normal: the element path is its
     * reference (the boxed algorithms answer a plain List). */
    if (!ll_normal_packs(ll)) return NULL;
    bool cplx = k == LL_CPLX;
    size_t argc = res->data.function.arg_count;
    Expr* a1 = argc >= 2 ? res->data.function.args[1] : NULL;
    size_t n = (size_t)ll_length(ll);
    if (n == 0) return NULL;
    double dn = (double)n;
    Expr* r = NULL;

    switch (op) {
    case STATS_MEAN: {
        if (argc != 1) return NULL;
        double re, im, inv = 1.0 / dn;
        ll_full_sum(ll, &re, &im);
        r = scalar_of(re * inv, im * inv, cplx);
        break;
    }
    case STATS_VARIANCE: case STATS_STD: case STATS_RMS: {
        if (argc != 1 || n < (op == STATS_RMS ? 1u : 2u)) return NULL;
        double v;
        if (op == STATS_RMS) {
            v = sqrt(ll_full_sumsq(ll, 0.0, 0.0) / dn);
        } else {
            double sr, si;
            ll_full_sum(ll, &sr, &si);
            double var_ = ll_full_sumsq(ll, sr / dn, si / dn) / (double)(n - 1);
            v = (op == STATS_STD) ? sqrt(var_) : var_;
        }
        r = expr_new_real(v);
        break;
    }
    case STATS_MOMENT: case STATS_CENTRAL_MOMENT: {
        int64_t ord;
        if (argc != 2 || cplx || !order_arg(a1, &ord)) return NULL;
        double mean = 0.0;
        if (op == STATS_CENTRAL_MOMENT) {
            double sr, si;
            ll_full_sum(ll, &sr, &si);
            mean = sr / dn;
        }
        r = expr_new_real(ll_full_moment(ll, mean, ord) / dn);
        break;
    }
    case STATS_SKEWNESS: case STATS_KURTOSIS: {
        if (argc != 1 || cplx) return NULL;
        int p = op == STATS_SKEWNESS ? 3 : 4;
        double sr, si;
        ll_full_sum(ll, &sr, &si);
        double mean = sr / dn;
        double m2 = ll_full_moment(ll, mean, 2) / dn;
        double mp = ll_full_moment(ll, mean, p) / dn;
        r = expr_new_real(mp / pow(m2, (double)p / 2.0));
        break;
    }
    case STATS_MEDIAN: {
        if (argc != 1 || cplx) return NULL;
        double* s = gather_real(ll);
        if (!s) return NULL;
        r = expr_new_real(nd_median_of(s, NDT_FLOAT64, 0, 1, n));
        free(s);
        return r;                             /* unguarded, as ndred_median */
    }
    case STATS_QUARTILES: {
        if (argc != 1 || cplx) return NULL;
        double* s = gather_real(ll);
        if (!s) return NULL;
        double q[3];
        nd_quartiles_f64(s, n, q);
        free(s);
        Expr* e[3] = { expr_new_real(q[0]), expr_new_real(q[1]), expr_new_real(q[2]) };
        return expr_new_function(expr_new_symbol(SYM_List), e, 3);
    }
    case STATS_MOVING_AVERAGE: case STATS_MOVING_MEDIAN: {
        if (argc != 2 || cplx || a1->type != EXPR_INTEGER) return NULL;
        int64_t w = a1->data.integer;
        if (w < 1 || (size_t)w > n) return NULL;
        size_t L = n - (size_t)w + 1;
        double* s = gather_real(ll);
        double* out = malloc(sizeof(double) * L);
        if (s && out) {
            if (op == STATS_MOVING_AVERAGE) {
                nd_moving_average_f64(s, n, (size_t)w, out);
                r = real_ll(out, L, true);
            } else {
                for (size_t i = 0; i < L; i++)
                    out[i] = nd_median_of(s, NDT_FLOAT64, i, 1, (size_t)w);
                r = real_ll(out, L, false);
            }
        }
        free(s);
        free(out);
        return r;
    }
    case STATS_EMA: {
        if (argc != 2 || cplx) return NULL;
        double al;
        if (a1->type == EXPR_REAL) al = a1->data.real;
        else if (a1->type == EXPR_INTEGER) al = (double)a1->data.integer;
        else return NULL;
        double* out = malloc(sizeof(double) * n);
        if (!out) return NULL;
        LLCursor c;
        const void* p;
        int dir;
        int64_t base, kk;
        double prev = 0.0;
        ll_cursor(&c, ll, 0, (int64_t)n);
        while ((kk = ll_cursor_next(&c, &p, &dir, &base)) > 0)
            for (int64_t t = 0; t < kk; t++) {
                double x = LL_SLOT(double, p, dir, t);
                prev = (base + t == 0) ? x : al * x + (1.0 - al) * prev;
                out[base + t] = prev;
            }
        r = real_ll(out, n, true);
        free(out);
        return r;
    }
    default:
        return NULL;
    }
    /* The packed reductions' nd_guard: a result outside the machine range is
     * recomputed on the element path. */
    if (ndred_result_degenerate(r)) { expr_free(r); return NULL; }
    return r;
}

/* ------------------------------------------------------------------------
 *  Shape 2: a matrix with a LinkedList in it
 * ---------------------------------------------------------------------- */

static bool is_seq_elem(const Expr* e) { return ll_is_seq(e); }

/* No element of the row is itself a sequence (the rows are vectors). */
static bool row_is_vector(const Expr* row) {
    const LinkedList* ll = expr_linkedlist(row);
    if (ll) {
        if (ll_kind(ll) != LL_BOXED) return true;
        for (int64_t i = 0, n = ll_length(ll); i < n; i++)
            if (is_seq_elem(ll_peek_boxed(ll, i))) return false;
        return true;
    }
    if (is_ndarray(row)) return row->data.ndarray.rank == 1;
    for (size_t i = 0; i < row->data.function.arg_count; i++)
        if (is_seq_elem(row->data.function.args[i])) return false;
    return true;
}

typedef struct {
    Expr** rows;      /* borrowed (List) or owned (LinkedList data) */
    LLElems e;
    size_t nrows, ncols;
    bool packed_ref;  /* the reference matrix is a PACKED machine-real array */
} Matrix;

/* A row whose Normal is a packed machine-real vector: a machine-real
 * LinkedList long enough to pack, or a packed real List. The List of the rows'
 * Normals -- the reference -- is then itself a packed rank-2 array, and the
 * packed kernels are its algorithms; any other row makes it a plain List of
 * rows, answered column by column on plain columns. */
static bool row_packs_real(const Expr* r) {
    const LinkedList* rl = expr_linkedlist(r);
    if (rl) return ll_kind(rl) == LL_REAL && ll_normal_packs(rl);
    return is_ndarray(r) && r->data.ndarray.rank == 1 && r->data.ndarray.dtype == NDT_FLOAT64;
}

/* `data` a List or LinkedList of >= 1 rows, each a List / LinkedList / rank-1
 * array of one common length >= 1, with a LinkedList somewhere. */
static bool matrix_open(const Expr* data, Matrix* m) {
    if (!ll_is_seq(data)) return false;
    bool any_ll = data->type == EXPR_LINKEDLIST;
    if (!ll_elems_open(data, &m->e)) return false;
    m->rows = m->e.v;
    m->nrows = m->e.n;
    m->packed_ref = true;
    if (m->nrows == 0) goto no;
    int64_t len = -1;
    for (size_t i = 0; i < m->nrows; i++) {
        const Expr* r = m->rows[i];
        int64_t k = ll_seq_length(r);
        if (k < 1 || (len >= 0 && k != len)) goto no;
        len = k;
        if (r->type == EXPR_LINKEDLIST) any_ll = true;
        if (!row_packs_real(r)) m->packed_ref = false;
        if (!row_is_vector(r)) goto no;
    }
    if (!any_ll) goto no;
    m->ncols = (size_t)len;
    return true;
no:
    ll_elems_close(&m->e);
    return false;
}

static void matrix_close(Matrix* m) { ll_elems_close(&m->e); }

/* Column j as a fresh BOXED LinkedList: the counterpart of the plain column
 * the List code's Transpose hands each per-column call. */
static Expr* matrix_column(const Matrix* m, size_t j) {
    Expr** v = malloc(sizeof(Expr*) * m->nrows);
    for (size_t i = 0; i < m->nrows; i++) v[i] = ll_seq_element(m->rows[i], (int64_t)j);
    LinkedList* c = ll_new_boxed(v, m->nrows);
    for (size_t i = 0; i < m->nrows; i++) expr_free(v[i]);
    free(v);
    return expr_new_linkedlist(c);
}

/* The column answers out[0..n) (consumed) in the image of the rows. */
static Expr* matrix_answer(const Matrix* m, Expr** out, size_t n) {
    return ll_seq_like(m->rows[0], out, n);
}

/* The packed rank-2 kernels: ndreduce.c's strided column reductions over the
 * rows gathered row-major. NULL declines (the per-column path answers). */
static Expr* machine_columns(Expr* res, StatsOp op, const Matrix* m) {
    size_t argc = res->data.function.arg_count;
    size_t R = m->nrows, C = m->ncols;
    int64_t ord = 0;
    switch (op) {
    case STATS_MEAN: case STATS_SKEWNESS: case STATS_KURTOSIS: case STATS_MEDIAN:
        if (argc != 1) return NULL;
        break;
    case STATS_VARIANCE: case STATS_STD:
        if (argc != 1 || R < 2) return NULL;
        break;
    case STATS_RMS:
        if (argc != 1) return NULL;
        break;
    case STATS_MOMENT: case STATS_CENTRAL_MOMENT:
        if (argc != 2 || !order_arg(res->data.function.args[1], &ord)) return NULL;
        break;
    default:
        return NULL;
    }
    double* buf = malloc(sizeof(double) * R * C);
    double* out = malloc(sizeof(double) * C);
    Expr* r = NULL;
    if (buf && out) {
        for (size_t i = 0; i < R; i++) {
            const LinkedList* rl = expr_linkedlist(m->rows[i]);
            if (rl) ll_gather(rl, 0, (int64_t)C, buf + i * C);
            else memcpy(buf + i * C, m->rows[i]->data.ndarray.data, sizeof(double) * C);
        }
        for (size_t j = 0; j < C; j++) {
            double re, im;
            switch (op) {
            case STATS_MEAN:
                nd_sum_strided(buf, NDT_FLOAT64, j, C, R, &re, &im);
                out[j] = re * (1.0 / (double)R);
                break;
            case STATS_VARIANCE: out[j] = nd_moment_at(buf, NDT_FLOAT64, j, C, R, 0); break;
            case STATS_STD:      out[j] = nd_moment_at(buf, NDT_FLOAT64, j, C, R, 1); break;
            case STATS_RMS:      out[j] = nd_moment_at(buf, NDT_FLOAT64, j, C, R, 2); break;
            case STATS_MOMENT:   out[j] = nd_rawmoment_at(buf, NDT_FLOAT64, j, C, R, ord); break;
            case STATS_CENTRAL_MOMENT: out[j] = nd_cm_at(buf, NDT_FLOAT64, j, C, R, ord); break;
            case STATS_SKEWNESS: out[j] = nd_stdmoment_at(buf, NDT_FLOAT64, j, C, R, 3, 1.5); break;
            case STATS_KURTOSIS: out[j] = nd_stdmoment_at(buf, NDT_FLOAT64, j, C, R, 4, 2.0); break;
            case STATS_MEDIAN:   out[j] = nd_median_of(buf, NDT_FLOAT64, j, C, R); break;
            default: break;
            }
        }
        bool ok = true;
        for (size_t j = 0; ok && op != STATS_MEDIAN && j < C; j++) ok = stats_machine_ok(out[j]);
        if (ok) {
            Expr** cells = malloc(sizeof(Expr*) * C);
            for (size_t j = 0; j < C; j++) cells[j] = expr_new_real(out[j]);
            r = matrix_answer(m, cells, C);
            free(cells);
        }
    }
    free(buf);
    free(out);
    return r;
}

/* The columnwise heads (the List code's MatrixQ arm). */
static bool columnwise_op(StatsOp op) {
    switch (op) {
    case STATS_MEAN: case STATS_VARIANCE: case STATS_STD: case STATS_RMS:
    case STATS_MOMENT: case STATS_CENTRAL_MOMENT: case STATS_SKEWNESS:
    case STATS_KURTOSIS: case STATS_MEDIAN: case STATS_QUARTILES:
    case STATS_IQR: case STATS_MEAN_DEVIATION: case STATS_MEDIAN_DEVIATION:
        return true;
    default:
        return false;
    }
}

/* Column j of a moment statistic as the List code reduces a matrix: through
 * the exact / symbolic trees of the whole matrix (stats_cm_elements,
 * stats_moment_elements) -- not the vector fast loops a column would take --
 * and Skewness / Kurtosis as m_p / m_2^(p/2) of those, as
 * stats_standardized_moment threads Power and Divide over the columns. */
static Expr* moment_column(Expr* res, StatsOp op, const Matrix* m, size_t j) {
    Expr** v = malloc(sizeof(Expr*) * m->nrows);
    for (size_t i = 0; i < m->nrows; i++) v[i] = ll_seq_element(m->rows[i], (int64_t)j);
    Expr* r;
    if (op == STATS_MOMENT) {
        r = stats_moment_elements(v, m->nrows, res->data.function.args[1]);
    } else if (op == STATS_CENTRAL_MOMENT) {
        r = stats_cm_elements(v, m->nrows, res->data.function.args[1]);
    } else {
        int p = op == STATS_SKEWNESS ? 3 : 4;
        Expr* pe = expr_new_integer(p);
        Expr* two = expr_new_integer(2);
        Expr* mp = stats_cm_elements(v, m->nrows, pe);
        Expr* m2 = stats_cm_elements(v, m->nrows, two);
        expr_free(pe);
        expr_free(two);
        if (!mp || !m2) {
            expr_free(mp); expr_free(m2);
            r = NULL;
        } else {
            Expr* denom = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                              (Expr*[]){ m2, make_rational(p, 2) }, 2));
            r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Divide),
                    (Expr*[]){ mp, denom }, 2));
        }
    }
    for (size_t i = 0; i < m->nrows; i++) expr_free(v[i]);
    free(v);
    return r;
}

/* head[column_j, rest...] for every column, in the image of the rows. */
static Expr* per_column(Expr* res, StatsOp op, const Matrix* m) {
    size_t argc = res->data.function.arg_count;
    Expr** out = malloc(sizeof(Expr*) * m->ncols);
    if (!out) return NULL;
    bool moment = op == STATS_MOMENT || op == STATS_CENTRAL_MOMENT ||
                  op == STATS_SKEWNESS || op == STATS_KURTOSIS;
    for (size_t j = 0; j < m->ncols; j++) {
        if (moment) {
            out[j] = moment_column(res, op, m, j);
            if (!out[j]) {                           /* irreducible: unevaluated */
                for (size_t k = 0; k < j; k++) expr_free(out[k]);
                free(out);
                return NULL;
            }
            continue;
        }
        Expr** a = malloc(sizeof(Expr*) * argc);
        a[0] = matrix_column(m, j);
        for (size_t t = 1; t < argc; t++) a[t] = expr_copy(res->data.function.args[t]);
        Expr* call = expr_new_function(expr_copy(res->data.function.head), a, argc);
        free(a);
        out[j] = eval_and_free(call);
    }
    Expr* r = matrix_answer(m, out, m->ncols);
    free(out);
    return r;
}

/* ------------------------------------------------------------------------
 *  Entry point
 * ---------------------------------------------------------------------- */

bool stats_ll_arm(Expr* res, StatsOp op, Expr** out) {
    *out = NULL;
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return false;
    Expr* data = res->data.function.args[0];
    if (!ll_is_seq(data)) return false;

    if (columnwise_op(op)) {
        Matrix m;
        if (matrix_open(data, &m)) {
            Expr* r = m.packed_ref ? machine_columns(res, op, &m) : NULL;
            if (!r) r = per_column(res, op, &m);
            matrix_close(&m);
            *out = r;
            return true;
        }
    }
    const LinkedList* ll = expr_linkedlist(data);
    if (!ll) return false;
    Expr* r = machine_vector(res, op, ll);
    if (r) { *out = r; return true; }
    return false;
}

Expr* stats_apply_head(const char* head, Expr* seq) {
    const LinkedList* ll = expr_linkedlist(seq);
    if (ll) {
        Expr* r = ll_call_with_elements(expr_new_symbol(head), ll);
        expr_free(seq);
        return r;
    }
    return expr_new_function(expr_new_symbol(SYM_Apply),
                             (Expr*[]){ expr_new_symbol(head), seq }, 2);
}
