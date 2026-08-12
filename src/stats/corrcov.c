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
#include "eval.h"
#include "arithmetic.h"   /* make_rational */
#include "sym_names.h"
#include "common.h"       /* builtin_arg_error */
#include "ndreduce.h"     /* nd_covariance / nd_correlation */
#include "ndarray.h"      /* is_ndarray */
#include <stdlib.h>
#include <stdbool.h>

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
static Expr* g_cov_vectors(Expr* v, Expr* w) {
    if (v->type != EXPR_FUNCTION || w->type != EXPR_FUNCTION) return NULL;
    size_t n = v->data.function.arg_count;
    if (n != w->data.function.arg_count || n < 2) return NULL;

    bool all_numeric = true, has_real = false, has_complex = false;
    for (size_t i = 0; i < n; i++) {
        bool cplx = false;
        if (!stats_is_numeric(v->data.function.args[i], NULL, &cplx)) { all_numeric = false; break; }
        if (cplx) has_complex = true;
        else if (v->data.function.args[i]->type == EXPR_REAL) has_real = true;
        cplx = false;
        if (!stats_is_numeric(w->data.function.args[i], NULL, &cplx)) { all_numeric = false; break; }
        if (cplx) has_complex = true;
        else if (w->data.function.args[i]->type == EXPR_REAL) has_real = true;
    }

    if (all_numeric && has_real && !has_complex) {
        double mv = 0.0, mw = 0.0;
        for (size_t i = 0; i < n; i++) {
            double xv = 0.0, xw = 0.0;
            stats_is_numeric(v->data.function.args[i], &xv, NULL);
            stats_is_numeric(w->data.function.args[i], &xw, NULL);
            mv += xv; mw += xw;
        }
        mv /= (double)n; mw /= (double)n;
        double s = 0.0;
        for (size_t i = 0; i < n; i++) {
            double xv = 0.0, xw = 0.0;
            stats_is_numeric(v->data.function.args[i], &xv, NULL);
            stats_is_numeric(w->data.function.args[i], &xw, NULL);
            s += (xv - mv) * (xw - mw);
        }
        return expr_new_real(s / (double)(n - 1));
    }

    /* Exact / complex / symbolic: (1/(n-1)) Sum (v_i - mu_v) Conjugate[w_i - mu_w]. */
    Expr* mu_v = eval_and_free(expr_new_function(expr_new_symbol(SYM_Mean),
                                                 (Expr*[]){ expr_copy(v) }, 1));
    Expr* mu_w = eval_and_free(expr_new_function(expr_new_symbol(SYM_Mean),
                                                 (Expr*[]){ expr_copy(w) }, 1));
    Expr** terms = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* dv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ expr_copy(v->data.function.args[i]),
                       expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ expr_new_integer(-1), expr_copy(mu_v) }, 2) }, 2));
        Expr* dw = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ expr_copy(w->data.function.args[i]),
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

static Expr* corrcov_dispatch(Expr* res, int correlate, const char* head) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return builtin_arg_error(head, argc, 1, 3);

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
