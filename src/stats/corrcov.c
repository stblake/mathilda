/* corrcov.c -- Covariance[] and Correlation[].
 *
 * Covariance[v, w]   covariance between two length-n vectors (a scalar)
 * Covariance[a, b]   p x q cross-covariance of the columns of two n-row matrices
 * Covariance[a]      p x p auto-covariance of a matrix, i.e. Covariance[a, a]
 * Correlation[...]   the same three shapes, normalized by the standard deviations
 *
 * For length-n vectors the covariance is
 *   (1/(n-1)) Sum_i (v_i - Mean[v]) Conjugate[w_i - Mean[w]]
 * (the conjugate is on the SECOND argument), and the correlation divides that by
 * StandardDeviation[v] StandardDeviation[w] (the (n-1) factors cancel). The
 * matrix forms apply the vector definition to each pair of columns.
 *
 * Following variance.c, the exact/complex/symbolic work is built as
 * sub-expressions and evaluated, so exact input yields exact output, complex
 * yields complex, and symbolic yields symbolic — with no int64-overflow risk. A
 * fast machine-double path covers real numeric vectors; an NDArray / packed-array
 * argument takes the buffer fast path in src/linalg/ndcorrcov.c.
 *
 * See stats.h and stats_common.h for the subsystem layout. */

#include "stats.h"
#include "stats_common.h"
#include <math.h>
#include "eval.h"
#include "arithmetic.h"   /* make_rational */
#include "sym_names.h"
#include "common.h"       /* builtin_arg_error */
#include "ndreduce.h"     /* nd_covariance / nd_correlation */
#include "ndarray.h"      /* is_ndarray */
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "linkedlist.h" /* LLReader */

/* ------------------------------------------------------------ small helpers */

/* Evaluate pred[x] and test it for the symbol True. */
static bool eval_pred(const char* pred, Expr* x) {
    Expr* call = expr_new_function(expr_new_symbol(pred), (Expr*[]){ expr_copy(x) }, 1);
    Expr* r = evaluate(call);
    expr_free(call);
    bool t = (r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_True);
    expr_free(r);
    return t;
}

/* StandardDeviation[v], evaluated. */
static Expr* g_std(Expr* v) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_StandardDeviation),
                                           (Expr*[]){ expr_copy(v) }, 1));
}

static bool corrcov_has_ndarray(const Expr* res) {
    for (size_t i = 0; i < res->data.function.arg_count; i++)
        if (is_ndarray(res->data.function.args[i])) return true;
    return false;
}

/* ------------------------------------------------------------ vector core */

/* Scalar covariance of two vectors v, w (both EXPR_FUNCTION List, equal length
 * n >= 2), or NULL to leave the call unevaluated. Mirrors the numeric split of
 * builtin_variance: a real-double fast path, otherwise build & evaluate the
 * exact/complex/symbolic expression. */
static Expr* g_cov_core(LLReader* ve, LLReader* we, size_t n);

/* Covariance of two vectors, each a List or a LinkedList (read through
 * LLReader, so both containers run g_cov_core). */
static Expr* g_cov_vectors(Expr* v, Expr* w) {
    LLReader ev, ew;
    if (!ll_rd_open(v, &ev)) return NULL;
    if (!ll_rd_open(w, &ew)) { ll_rd_close(&ev); return NULL; }
    Expr* r = (ev.n == ew.n && ev.n >= 2) ? g_cov_core(&ev, &ew, ev.n) : NULL;
    ll_rd_close(&ev);
    ll_rd_close(&ew);
    return r;
}

static Expr* g_cov_core(LLReader* ve, LLReader* we, size_t n) {

    bool all_numeric = true, has_real = false, has_complex = false;
    for (size_t i = 0; i < n; i++) {
        bool cplx = false;
        if (!stats_is_numeric(ll_rd_at(ve, i), NULL, &cplx)) { all_numeric = false; break; }
        if (cplx) has_complex = true;
        else if (ll_rd_at(ve, i)->type == EXPR_REAL) has_real = true;
        cplx = false;
        if (!stats_is_numeric(ll_rd_at(we, i), NULL, &cplx)) { all_numeric = false; break; }
        if (cplx) has_complex = true;
        else if (ll_rd_at(we, i)->type == EXPR_REAL) has_real = true;
    }

    if (all_numeric && has_real && !has_complex) {
        double mv = 0.0, mw = 0.0;
        for (size_t i = 0; i < n; i++) {
            double xv = 0.0, xw = 0.0;
            stats_is_numeric(ll_rd_at(ve, i), &xv, NULL);
            stats_is_numeric(ll_rd_at(we, i), &xw, NULL);
            mv += xv; mw += xw;
        }
        mv /= (double)n; mw /= (double)n;
        double s = 0.0;
        for (size_t i = 0; i < n; i++) {
            double xv = 0.0, xw = 0.0;
            stats_is_numeric(ll_rd_at(ve, i), &xv, NULL);
            stats_is_numeric(ll_rd_at(we, i), &xw, NULL);
            s += (xv - mv) * (xw - mw);
        }
        double c = s / (double)(n - 1);
        if (stats_machine_ok(c) && isfinite(mv) && isfinite(mw))
            return expr_new_real(c);   /* else decline to the exact path (C3) */
    }

    /* Exact / complex / symbolic: (1/(n-1)) Sum (v_i - mu_v) Conjugate[w_i - mu_w]. */
    Expr* mu_v = stats_mean_vec(ve);
    Expr* mu_w = stats_mean_vec(we);
    Expr** terms = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* dv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ expr_copy(ll_rd_at(ve, i)),
                       expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ expr_new_integer(-1), expr_copy(mu_v) }, 2) }, 2));
        Expr* dw = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ expr_copy(ll_rd_at(we, i)),
                       expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ expr_new_integer(-1), expr_copy(mu_w) }, 2) }, 2));
        Expr* cw = eval_and_free(expr_new_function(expr_new_symbol(SYM_Conjugate),
                                                   (Expr*[]){ dw }, 1));
        terms[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                                                   (Expr*[]){ dv, cw }, 2));
    }
    expr_free(mu_v);
    expr_free(mu_w);

    Expr* sum = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms, n));
    free(terms);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ make_rational(1, (int64_t)n - 1), sum }, 2));
}

/* Scalar correlation of two vectors: Cov[v,w] / (StandardDeviation[v] StandardDeviation[w]). */
static Expr* g_corr_vectors(Expr* v, Expr* w) {
    Expr* cov = g_cov_vectors(v, w);
    if (!cov) return NULL;
    Expr* sv = g_std(v);
    Expr* sw = g_std(w);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ cov,
                   expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ sv, expr_new_integer(-1) }, 2),
                   expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ sw, expr_new_integer(-1) }, 2) }, 3));
}

/* ------------------------------------------------------------ matrix core */

/* p x q (co)variance / correlation matrix over the columns of a and b (each an
 * n-row matrix). `correlate` selects Correlation; `auto_form` marks the
 * one-argument Covariance[a] / Correlation[a] case (b == a) — its correlation
 * diagonal is set to exactly 1. Returns NULL (unevaluated) if the transpose is
 * not a matrix or there are fewer than two observations. */
static Expr* g_matrix(Expr* a, Expr* b, int correlate, int auto_form) {
    Expr* ta = eval_and_free(expr_new_function(expr_new_symbol(SYM_Transpose),
                                               (Expr*[]){ expr_copy(a) }, 1));
    if (ta->type != EXPR_FUNCTION || ta->data.function.arg_count == 0
        || ta->data.function.args[0]->type != EXPR_FUNCTION
        || ta->data.function.args[0]->data.function.arg_count < 2) {
        expr_free(ta);
        return NULL;                 /* not a matrix, or fewer than 2 observations */
    }
    Expr* tb;
    if (auto_form) {
        tb = ta;
    } else {
        tb = eval_and_free(expr_new_function(expr_new_symbol(SYM_Transpose),
                                             (Expr*[]){ expr_copy(b) }, 1));
        if (tb->type != EXPR_FUNCTION) { expr_free(ta); expr_free(tb); return NULL; }
    }
    size_t p = ta->data.function.arg_count;
    size_t q = tb->data.function.arg_count;

    /* Auto-correlation of real data yields an all-real matrix, so its unit
     * diagonal must be the Real 1. (not the exact Integer 1) to keep the matrix
     * one uniform type — a mixed exact/real matrix leaves SymmetricMatrixQ,
     * PositiveSemidefiniteMatrixQ and == unevaluated. Symbolic/exact data keeps
     * the exact 1, matching Mathematica's Correlation[{{a,b},{c,d}}]. */
    bool real_diag = common_has_machine_real(a);

    Expr** sds_a = NULL;
    Expr** sds_b = NULL;
    if (correlate) {
        sds_a = malloc(sizeof(Expr*) * p);
        for (size_t i = 0; i < p; i++) sds_a[i] = g_std(ta->data.function.args[i]);
        if (auto_form) {
            sds_b = sds_a;
        } else {
            sds_b = malloc(sizeof(Expr*) * q);
            for (size_t j = 0; j < q; j++) sds_b[j] = g_std(tb->data.function.args[j]);
        }
    }

    Expr** rows = malloc(sizeof(Expr*) * p);
    for (size_t i = 0; i < p; i++) {
        Expr** cells = malloc(sizeof(Expr*) * q);
        for (size_t j = 0; j < q; j++) {
            if (correlate && auto_form && i == j) {
                cells[j] = real_diag ? expr_new_real(1.0) : expr_new_integer(1);
                continue;
            }
            Expr* cov = g_cov_vectors(ta->data.function.args[i], tb->data.function.args[j]);
            if (!cov) cov = expr_new_integer(0);      /* degenerate column; keep shape */
            if (correlate) {
                cells[j] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ cov,
                               expr_new_function(expr_new_symbol(SYM_Power),
                                   (Expr*[]){ expr_copy(sds_a[i]), expr_new_integer(-1) }, 2),
                               expr_new_function(expr_new_symbol(SYM_Power),
                                   (Expr*[]){ expr_copy(sds_b[j]), expr_new_integer(-1) }, 2) }, 3));
            } else {
                cells[j] = cov;
            }
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), cells, q);
        free(cells);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), rows, p);
    free(rows);

    if (correlate) {
        for (size_t i = 0; i < p; i++) expr_free(sds_a[i]);
        free(sds_a);
        if (!auto_form) {
            for (size_t j = 0; j < q; j++) expr_free(sds_b[j]);
            free(sds_b);
        }
    }
    expr_free(ta);
    if (!auto_form) expr_free(tb);
    return result;
}

/* ------------------------------------------------------------ dispatch */

/* ------------------------------------------------------------------------
 * LinkedList arguments. A vector is a List or LinkedList none of whose
 * elements is a sequence; a matrix a List or LinkedList of such vectors, all
 * of one length. Two machine-real LinkedList vectors take the packed kernel
 * (nd_corrcov_vec_f64) on raw gathers of their slots, a matrix of machine-real
 * LinkedList rows the packed matrix kernel (nd_corrcov_matrix_f64) on the rows
 * gathered row-major; everything else reads the elements (g_cov_vectors) and
 * the columns as fresh LinkedLists. A matrix answer is in the image of the
 * rows (LinkedList rows, LinkedList of LinkedLists).
 * ---------------------------------------------------------------------- */

static bool cc_is_vector(const Expr* e) {
    const LinkedList* ll = expr_linkedlist(e);
    if (ll) {
        if (ll_kind(ll) != LL_BOXED) return true;
        for (int64_t i = 0, n = ll_length(ll); i < n; i++)
            if (ll_is_seq(ll_peek_boxed(ll, i))) return false;
        return true;
    }
    if (!ll_is_seq(e)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (ll_is_seq(e->data.function.args[i])) return false;
    return true;
}

/* Rows of a matrix (owned through *rows); *cols their common length. */
static bool cc_matrix(const Expr* e, LLReader* rows, size_t* cols) {
    if (!ll_rd_open(e, rows)) return false;
    bool ok = rows->n >= 1;
    int64_t len = -1;
    for (size_t i = 0; ok && i < rows->n; i++) {
        int64_t k = ll_seq_length(ll_rd_at(rows, i));
        if (k < 1 || (len >= 0 && k != len) || !cc_is_vector(ll_rd_at(rows, i))) ok = false;
        len = k;
    }
    if (!ok) { ll_rd_close(rows); return false; }
    *cols = (size_t)len;
    return true;
}

static bool cc_involves_ll(const Expr* e) {
    if (e->type == EXPR_LINKEDLIST) return true;
    if (!ll_is_seq(e)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (e->data.function.args[i]->type == EXPR_LINKEDLIST) return true;
    return false;
}

/* A machine-real vector whose Normal is PACKED (a typed LinkedList long enough
 * to pack, or a packed real List): the packed kernels are then the reference
 * algorithms; a plain Normal is answered on the elements. */
static bool cc_real_packed(const Expr* e) {
    const LinkedList* ll = expr_linkedlist(e);
    if (ll) return ll_kind(ll) == LL_REAL && ll_normal_packs(ll);
    return is_ndarray(e) && e->data.ndarray.rank == 1 && e->data.ndarray.dtype == NDT_FLOAT64;
}

/* Rows gathered row-major (n x cols) when every one packs as machine reals --
 * the reference matrix is then a packed array -- else NULL. */
static double* cc_gather_rows(LLReader* rows, size_t cols) {
    for (size_t i = 0; i < rows->n; i++)
        if (!cc_real_packed(ll_rd_at(rows, i))) return NULL;
    double* A = malloc(sizeof(double) * rows->n * cols);
    if (!A) return NULL;
    for (size_t i = 0; i < rows->n; i++) {
        const LinkedList* rl = expr_linkedlist(ll_rd_at(rows, i));
        if (rl) ll_gather(rl, 0, (int64_t)cols, A + i * cols);
        else memcpy(A + i * cols, ll_rd_at(rows, i)->data.ndarray.data, sizeof(double) * cols);
    }
    return A;
}

/* Column j of a matrix as a fresh BOXED LinkedList: the counterpart of the
 * plain column g_matrix's Transpose makes. */
static Expr* cc_column(LLReader* rows, size_t j) {
    LinkedList* c = ll_boxed_begin();
    for (size_t i = 0; i < rows->n; i++)
        ll_push_back(c, ll_seq_element(ll_rd_at(rows, i), (int64_t)j));
    ll_boxed_end(c);
    return expr_new_linkedlist(c);
}

/* A p x q matrix in the image of the rows (`proto_row` is one of them):
 * cells[i*q + j] consumed. */
static Expr* cc_build_matrix(const Expr* proto_row, Expr** cells, size_t p, size_t q) {
    Expr** rows = malloc(sizeof(Expr*) * p);
    for (size_t i = 0; i < p; i++) rows[i] = ll_seq_like(proto_row, cells + i * q, q);
    Expr* r = ll_seq_like(proto_row, rows, p);
    free(rows);
    return r;
}

static Expr* cc_matrix_ll(Expr* a, Expr* b, int correlate, int auto_form) {
    LLReader ra, rb;
    size_t p, q;
    if (!cc_matrix(a, &ra, &p)) return NULL;
    if (auto_form) { rb = ra; q = p; }
    else if (!cc_matrix(b, &rb, &q)) { ll_rd_close(&ra); return NULL; }
    Expr* result = NULL;
    size_t n = ra.n;
    if (n != rb.n || n < 2) goto done;

    /* The packed matrix kernel, when every row is machine-real. */
    double* A = cc_gather_rows(&ra, p);
    double* B = auto_form ? A : (A ? cc_gather_rows(&rb, q) : NULL);
    if (A && B) {
        double* C = malloc(sizeof(double) * p * q);
        if (C && nd_corrcov_matrix_f64(A, B, (int)n, (int)p, (int)q, correlate, C)) {
            Expr** cells = malloc(sizeof(Expr*) * p * q);
            for (size_t k = 0; k < p * q; k++) cells[k] = expr_new_real(C[k]);
            result = cc_build_matrix(ll_rd_at(&ra, 0), cells, p, q);
            free(cells);
        }
        free(C);
    }
    if (A) free(A);
    if (B && B != A) free(B);
    if (result) goto done;

    {   /* Generic: the columns as LinkedLists, cell by cell as g_matrix. */
        Expr** ca = malloc(sizeof(Expr*) * p);
        Expr** cb = auto_form ? ca : malloc(sizeof(Expr*) * q);
        for (size_t i = 0; i < p; i++) ca[i] = cc_column(&ra, i);
        if (!auto_form) for (size_t j = 0; j < q; j++) cb[j] = cc_column(&rb, j);
        bool real_diag = false;
        for (size_t i = 0; i < n && !real_diag; i++) {
            const LinkedList* rl = expr_linkedlist(ll_rd_at(&ra, i));
            if (rl) {
                if (ll_kind(rl) == LL_REAL || ll_kind(rl) == LL_CPLX) real_diag = true;
                else if (ll_kind(rl) == LL_BOXED)
                    for (int64_t j = 0; j < ll_length(rl) && !real_diag; j++)
                        real_diag = common_has_machine_real(ll_peek_boxed(rl, j));
            } else {
                real_diag = common_has_machine_real(ll_rd_at(&ra, i));
            }
        }
        Expr** cells = malloc(sizeof(Expr*) * p * q);
        for (size_t i = 0; i < p; i++)
            for (size_t j = 0; j < q; j++) {
                Expr** cell = &cells[i * q + j];
                if (correlate && auto_form && i == j) {
                    *cell = real_diag ? expr_new_real(1.0) : expr_new_integer(1);
                    continue;
                }
                Expr* cov = g_cov_vectors(ca[i], cb[j]);
                if (!cov) cov = expr_new_integer(0);      /* degenerate column; keep shape */
                if (correlate) {
                    *cell = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ cov,
                                   expr_new_function(expr_new_symbol(SYM_Power),
                                       (Expr*[]){ g_std(ca[i]), expr_new_integer(-1) }, 2),
                                   expr_new_function(expr_new_symbol(SYM_Power),
                                       (Expr*[]){ g_std(cb[j]), expr_new_integer(-1) }, 2) }, 3));
                } else {
                    *cell = cov;
                }
            }
        result = cc_build_matrix(ll_rd_at(&ra, 0), cells, p, q);
        free(cells);
        for (size_t i = 0; i < p; i++) expr_free(ca[i]);
        free(ca);
        if (!auto_form) { for (size_t j = 0; j < q; j++) expr_free(cb[j]); free(cb); }
    }
done:
    ll_rd_close(&ra);
    if (!auto_form) ll_rd_close(&rb);
    return result;
}

/* The n doubles of a cc_real_packed vector. */
static void cc_gather_vec(const Expr* e, double* out, int64_t n) {
    const LinkedList* ll = expr_linkedlist(e);
    if (ll) ll_gather(ll, 0, n, out);
    else memcpy(out, e->data.ndarray.data, sizeof(double) * (size_t)n);
}

/* Two vectors with a LinkedList among them. */
static Expr* cc_vectors_ll(Expr* a, Expr* b, int correlate) {
    if (cc_real_packed(a) && cc_real_packed(b)) {
        int64_t n = ll_seq_length(a);
        if (n == ll_seq_length(b)) {
            double* v = malloc(sizeof(double) * (size_t)(n ? n : 1));
            double* w = malloc(sizeof(double) * (size_t)(n ? n : 1));
            double r;
            bool ok = v && w;
            if (ok) { cc_gather_vec(a, v, n); cc_gather_vec(b, w, n); }
            ok = ok && nd_corrcov_vec_f64(v, w, (size_t)n, correlate, &r);
            free(v);
            free(w);
            if (ok) return expr_new_real(r);
        }
    }
    return correlate ? g_corr_vectors(a, b) : g_cov_vectors(a, b);
}

static bool corrcov_ll(Expr* res, int correlate, Expr** out) {
    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;
    if (argc == 2 && (cc_involves_ll(a[0]) || cc_involves_ll(a[1]))) {
        if (cc_is_vector(a[0]) && cc_is_vector(a[1])) *out = cc_vectors_ll(a[0], a[1], correlate);
        else *out = cc_matrix_ll(a[0], a[1], correlate, 0);
        return true;
    }
    if (argc == 1 && cc_involves_ll(a[0])) {
        *out = cc_is_vector(a[0]) ? NULL : cc_matrix_ll(a[0], a[0], correlate, 1);
        return true;
    }
    return false;
}

static Expr* corrcov_dispatch(Expr* res, int correlate, const char* head) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return builtin_arg_error(head, argc, 1, 3);

    /* A LinkedList argument: computed on its slots (above). */
    { Expr* r; if (corrcov_ll(res, correlate, &r)) return r; }

    /* NDArray / packed-array buffer fast path. */
    if (corrcov_has_ndarray(res))
        return correlate ? nd_correlation(res) : nd_covariance(res);

    if (argc == 2) {
        Expr* a = res->data.function.args[0];
        Expr* b = res->data.function.args[1];
        if (eval_pred(SYM_MatrixQ, a) && eval_pred(SYM_MatrixQ, b))
            return g_matrix(a, b, correlate, 0);
        if (eval_pred(SYM_VectorQ, a) && eval_pred(SYM_VectorQ, b))
            return correlate ? g_corr_vectors(a, b) : g_cov_vectors(a, b);
        return NULL;
    }
    if (argc == 1) {
        Expr* a = res->data.function.args[0];
        if (eval_pred(SYM_MatrixQ, a)) return g_matrix(a, a, correlate, 1);
        return NULL;                 /* a single vector / scalar: unevaluated */
    }
    return NULL;                     /* argc == 3: unhandled, leave unevaluated */
}

Expr* builtin_covariance(Expr* res)  { return corrcov_dispatch(res, 0, "Covariance"); }
Expr* builtin_correlation(Expr* res) { return corrcov_dispatch(res, 1, "Correlation"); }
