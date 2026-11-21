/* central_moment.c -- CentralMoment[].
 * Split from stats.c; see stats.h and stats_common.h for the subsystem layout.
 *
 * CentralMoment[data, r]              — the r-th moment about the mean,
 *                                       mu~_r = (1/n) Sum[(x_i - mu_1)^r], where
 *                                       mu_1 = Mean[data]. For a matrix / array the
 *                                       reduction is columnwise over the first axis
 *                                       (equivalently ArrayReduce[CentralMoment[#,r]&, x, 1]).
 * CentralMoment[data, {r1, ..., rm}]  — the multivariate mixed central moment,
 *                                       (1/n) Sum_i Product_j (x[[i,j]] - mu_1[[j]])^r_j,
 *                                       summing the first axis and taking a product over
 *                                       the second (its length must equal Length[{r1,...}]).
 *
 * The design mirrors Variance (a central moment is Variance without the n/(n-1)
 * bias correction): divide by n (not n-1), raise to the power r (not square),
 * n >= 1 suffices, and there is no Conjugate — a central moment is (x-mu)^r, not
 * |x-mu|^2. Numeric real vectors take a tight C loop (and packed / NDArray inputs
 * a machine-buffer fast path via ndred_central_moment); every other case — exact,
 * symbolic, matrix/array, multivariate — is built as an expression and handed to
 * the evaluator, which already knows how to be exact or symbolic. */

#include "stats.h"
#include "stats_common.h"
#include "eval.h"
#include "arithmetic.h"   /* is_rational */
#include "sym_names.h"
#include "assoc.h"
#include "ndreduce.h"
#include "linkedlist.h" /* LLReader */
#include <stdlib.h>

/* Integer power by squaring on a double (exponent >= 0). */
static double cm_ipow(double b, int64_t e) {
    double r = 1.0;
    while (e > 0) {
        if (e & 1) r *= b;
        b *= b;
        e >>= 1;
    }
    return r;
}

/* True when e is a List[...] expression (a vector, matrix row, or subarray). */
static bool cm_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/* Exact / symbolic vector base case: Mean[(data - Mean[data])^r]. Plus and Power
 * are Listable, so (data - mu) and (...)^r thread over the vector elementwise and
 * the outer Mean divides by n. Reuses Mean's exact-rational and symbolic paths. */
Expr* stats_cm_elements(LLReader* rd, Expr* order) {
    size_t n = rd->n;
    Expr* mu = stats_mean_vec(rd);
    if (!mu) return NULL;
    Expr* neg_mu = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_integer(-1), mu }, 2));
    Expr** p = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* c = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                      (Expr*[]){ expr_copy(ll_rd_at(rd, i)), expr_copy(neg_mu) }, 2));
        p[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                   (Expr*[]){ c, expr_copy(order) }, 2));
    }
    expr_free(neg_mu);
    Expr* r = stats_mean_of(p, n);
    for (size_t i = 0; i < n; i++) expr_free(p[i]);
    free(p);
    return r;
}

static Expr* cm_vector_symbolic(Expr* data, Expr* order) {
    if (data->type == EXPR_LINKEDLIST) {
        /* The same tree over a LinkedList's elements (stats_cm_elements). */
        LLReader e;
        ll_rd_open(data, &e);
        Expr* r = stats_cm_elements(&e, order);
        ll_rd_close(&e);
        return r;
    }
    Expr* mu = eval_and_free(expr_new_function(expr_new_symbol(SYM_Mean),
                   (Expr*[]){ expr_copy(data) }, 1));
    Expr* neg_mu = expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_integer(-1), mu }, 2);  /* adopts mu */
    Expr* centered = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                         (Expr*[]){ expr_copy(data), neg_mu }, 2));
    Expr* powered = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ centered, expr_copy(order) }, 2));
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Mean),
               (Expr*[]){ powered }, 1));
}

/* Scalar order r over an array of depth >= 2: reduce the first axis columnwise.
 * mu = Mean[data] is rank k-1; for each first-axis slice sub, term = (sub - mu)^r
 * (elementwise, rank k-1); the result is Mean over the slices — Mean[{term_i}]
 * divides the first axis by n. For a matrix this is the columnwise vector of
 * moments; for higher rank the columnwise array. Correct for every rank because mu
 * is subtracted per slice (data - Mean[data] would thread row-wise, not column). */
static Expr* cm_columnwise(Expr* data, Expr* order) {
    size_t n = data->data.function.arg_count;
    Expr* mu = eval_and_free(expr_new_function(expr_new_symbol(SYM_Mean),
                   (Expr*[]){ expr_copy(data) }, 1));
    Expr** terms = malloc(sizeof(Expr*) * n);
    if (!terms) { expr_free(mu); return NULL; }
    for (size_t i = 0; i < n; i++) {
        Expr* sub = data->data.function.args[i];              /* borrowed, rank k-1 */
        Expr* neg_mu = expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ expr_new_integer(-1), expr_copy(mu) }, 2);
        Expr* centered = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                             (Expr*[]){ expr_copy(sub), neg_mu }, 2));
        terms[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ centered, expr_copy(order) }, 2));
    }
    expr_free(mu);
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), terms, n);
    free(terms);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Mean),
               (Expr*[]){ list }, 1));
}

/* Multivariate order {r1, ..., rm} over an array of depth >= 2. mu = Mean[data]
 * (rank k-1). For each first-axis slice sub, term = Times @@ (sub - mu)^rvec: the
 * elementwise Power threads rvec over the slice's second axis (raising block j to
 * r_j) and Times @@ takes the product over that axis, leaving a rank k-2 block.
 * The result is Mean over the slices (dividing the first axis by n). Requires the
 * slice's second-axis length to equal Length[rvec]. */
static Expr* cm_multivariate(Expr* data, Expr* rvec) {
    LLReader rows;
    if (!ll_rd_open(data, &rows)) return NULL;
    Expr* first = ll_rd_at(&rows, 0);
    size_t m = rvec->data.function.arg_count;
    if (!ll_is_seq(first) || ll_seq_length(first) != (int64_t)m) {   /* depth >= 2, axis 2 = order */
        ll_rd_close(&rows);
        return NULL;
    }

    size_t n = rows.n;
    Expr* mu = eval_and_free(expr_new_function(expr_new_symbol(SYM_Mean),
                   (Expr*[]){ expr_copy(data) }, 1));
    Expr** terms = malloc(sizeof(Expr*) * n);
    if (!terms) { expr_free(mu); ll_rd_close(&rows); return NULL; }
    for (size_t i = 0; i < n; i++) {
        Expr* sub = ll_rd_at(&rows, i);                               /* borrowed, rank k-1 */
        Expr* neg_mu = expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ expr_new_integer(-1), expr_copy(mu) }, 2);
        Expr* centered = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                             (Expr*[]){ expr_copy(sub), neg_mu }, 2));
        Expr* powered = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                            (Expr*[]){ centered, expr_copy(rvec) }, 2));  /* threads r over axis 2 */
        terms[i] = eval_and_free(stats_apply_head(SYM_Times, powered));   /* Times @@ powered */
    }
    expr_free(mu);
    ll_rd_close(&rows);
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), terms, n);
    free(terms);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Mean),
               (Expr*[]){ list }, 1));
}

Expr* builtin_central_moment(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* data = res->data.function.args[0];
    Expr* order = res->data.function.args[1];

    /* CentralMoment[assoc, r] works on the association's values. */
    if (is_association(data)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    /* NDArray / packed fast path (scalar integer order); see ndreduce.c. A
     * list-valued or non-integer order, an int64 buffer, or a complex buffer
     * degrades back to a nested List and re-enters here. */
    if (ndred_call_has_ndarray(res)) return ndred_central_moment(res);

    /* LinkedList: the packed kernel on the chunks, or columnwise (stats_ll.c). */
    { Expr* r; if (!cm_is_list(order) && stats_ll_arm(res, STATS_CENTRAL_MOMENT, &r)) return r; }

    if (!ll_is_seq(data)) return NULL;
    if (ll_seq_length(data) == 0) return NULL;

    /* Multivariate order {r1, ..., rm}. */
    if (cm_is_list(order)) return cm_multivariate(data, order);

    LLReader e;
    ll_rd_open(data, &e);
    size_t n = e.n;

    /* Scalar order over an array (depth >= 2): columnwise over the first axis. */
    if (cm_is_list(ll_rd_at(&e, 0))) { ll_rd_close(&e); return cm_columnwise(data, order); }

    /* Scalar order over a vector (depth 1). Fast real path: every element is a
     * real number, at least one is inexact (approximate input -> approximate
     * output), and the order is a non-negative machine integer. */
    if (order->type == EXPR_INTEGER && order->data.integer >= 0) {
        bool all_numeric = true, has_real = false;
        for (size_t i = 0; i < n; i++) {
            Expr* x = ll_rd_at(&e, i);
            if (x->type == EXPR_REAL) has_real = true;
            else if (x->type == EXPR_INTEGER || is_rational(x, NULL, NULL)) { /* numeric */ }
            else { all_numeric = false; break; }
        }
        if (all_numeric && has_real) {
            int64_t r = order->data.integer;
            double mean = 0.0;
            for (size_t i = 0; i < n; i++) {
                double v = 0.0;
                stats_is_numeric(ll_rd_at(&e, i), &v, NULL);
                mean += v;
            }
            mean /= (double)n;
            double acc = 0.0;
            for (size_t i = 0; i < n; i++) {
                double v = 0.0;
                stats_is_numeric(ll_rd_at(&e, i), &v, NULL);
                acc += cm_ipow(v - mean, r);
            }
            double m = acc / (double)n;
            if (stats_machine_ok(m)) { ll_rd_close(&e); return expr_new_real(m); }   /* else decline (C3) */
        }
    }
    ll_rd_close(&e);

    /* Exact / symbolic vector, or symbolic / non-integer order. */
    return cm_vector_symbolic(data, order);
}
