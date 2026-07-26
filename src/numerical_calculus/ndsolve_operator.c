/* Mathilda — build a compiled linear operator from a symbolic MoL RHS. */
#include "ndsolve_operator.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../numeric.h"
#include <stdlib.h>
#include <string.h>
#ifdef USE_LAPACK
#include "../linalg/lapack.h"    /* cblas_dgbmv / cblas_dgemv */
#endif

/* out = A·Y for the compiled operator: banded BLAS dgbmv when a packed band is
 * available, dense BLAS dgemv otherwise, and a scalar loop with no BLAS. */
void nd_operator_matvec(const NdOperator* op, const double* Y, double* out) {
    size_t n = op->n;
#ifdef USE_LAPACK
    if (op->banded && op->AB) {
        /* dgbmv: y = alpha*A*x + beta*y, A in band storage (kl+ku+1 x n). */
        cblas_dgbmv(CblasColMajor, CblasNoTrans, (int)n, (int)n, op->kl, op->ku,
                    1.0, op->AB, op->kl + op->ku + 1, Y, 1, 0.0, out, 1);
        return;
    }
    /* A is row-major n x n; dgemv on it as column-major A^T with Trans reads the
     * same memory as row-major A, giving A·Y. */
    cblas_dgemv(CblasColMajor, CblasTrans, (int)n, (int)n, 1.0, op->A, (int)n,
                Y, 1, 0.0, out, 1);
    return;
#else
    for (size_t i = 0; i < n; i++) {
        const double* Ai = &op->A[i * n];
        double acc = 0.0;
        for (size_t j = 0; j < n; j++) acc += Ai[j] * Y[j];
        out[i] = acc;
    }
#endif
}

/* Collect the reduced-state indices (symbols named "NDSolve`w<k>") that appear
 * in `e`, marking seen[k] = true.  This bounds the coupling of a row to just the
 * states it actually references, so only those columns are differentiated. */
static void nd_collect_state_indices(const Expr* e, size_t d, bool* seen) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        const char* tick = nm ? strrchr(nm, '`') : NULL;
        if (tick && tick[1] == 'w') {
            char* end = NULL;
            long k = strtol(tick + 2, &end, 10);
            if (end && *end == '\0' && k >= 0 && (size_t)k < d) seen[(size_t)k] = true;
        }
        return;
    }
    if (e->type == EXPR_FUNCTION) {
        nd_collect_state_indices(e->data.function.head, d, seen);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            nd_collect_state_indices(e->data.function.args[i], d, seen);
    }
}

void nd_operator_free(NdOperator* op) {
    if (!op) return;
    free(op->A);
    free(op->AB);
    free(op->s0);
    if (op->st) { for (size_t i = 0; i < op->n; i++) expr_free(op->st[i]); free(op->st); }
    free(op);
}

NdOperator* nd_operator_try_build(NdProblem* P) {
    size_t d = P->d;
    if (d == 0) return NULL;
    NumericSpec spec = P->spec;

    double* A = calloc(d * d, sizeof(double));
    Expr** st = malloc(sizeof(Expr*) * d);        /* per-node forcing s_i(t)   */
    bool* seen = malloc(sizeof(bool) * d);
    if (!A || !st || !seen) { free(A); free(st); free(seen); return NULL; }
    for (size_t i = 0; i < d; i++) st[i] = NULL;

    bool ok = true;
    for (size_t i = 0; i < d && ok; i++) {
        memset(seen, 0, sizeof(bool) * d);
        nd_collect_state_indices(P->f[i], d, seen);
        /* A[i][j] = d f_i / d y_j for each referenced j; must be a constant. */
        Expr** lit0 = malloc(sizeof(Expr*) * d);   /* state -> 0 (for forcing)  */
        Expr** sub0 = malloc(sizeof(Expr*) * d);
        size_t nz = 0;
        for (size_t j = 0; j < d && ok; j++) {
            if (!seen[j]) continue;
            Expr* dargs[2] = { expr_copy(P->f[i]), expr_copy(P->ysym[j]) };
            Expr* dc = eval_and_free(expr_new_function(expr_new_symbol(SYM_D), dargs, 2));
            double v;
            if (dc && nd_eval_to_double(dc, spec, &v)) {
                A[i * d + j] = v;
            } else {
                ok = false;    /* nonlinear or t-dependent coefficient */
            }
            expr_free(dc);
            lit0[nz] = expr_copy(P->ysym[j]);
            sub0[nz] = expr_new_integer(0);
            nz++;
        }
        if (ok) {
            /* forcing s_i(t) = f_i with every referenced state set to 0 */
            Expr* si = nz ? nd_replace_all(expr_copy(P->f[i]), lit0, sub0, nz)
                          : expr_copy(P->f[i]);
            st[i] = si;
        }
        for (size_t k = 0; k < nz; k++) { expr_free(lit0[k]); expr_free(sub0[k]); }
        free(lit0); free(sub0);
    }
    free(seen);

    if (!ok) {
        free(A);
        for (size_t i = 0; i < d; i++) expr_free(st[i]);
        free(st);
        return NULL;
    }

    /* bandwidth of A */
    int kl = 0, ku = 0;
    for (size_t i = 0; i < d; i++)
        for (size_t j = 0; j < d; j++)
            if (A[i * d + j] != 0.0) {
                if ((long)i - (long)j > kl) kl = (int)((long)i - (long)j);
                if ((long)j - (long)i > ku) ku = (int)((long)j - (long)i);
            }

    /* try to fold the forcing to a constant vector s0 (autonomous / const BCs) */
    double* s0 = malloc(sizeof(double) * d);
    bool const_forcing = true;
    for (size_t i = 0; i < d; i++) {
        double v;
        if (st[i] && nd_eval_to_double(st[i], spec, &v)) s0[i] = v;
        else { const_forcing = false; break; }
    }

    NdOperator* op = malloc(sizeof(NdOperator));
    op->n = d; op->A = A; op->kl = kl; op->ku = ku;
    op->banded = ((size_t)(kl + ku + 1) <= d / 2 + 2);
    op->AB = NULL;
#ifdef USE_LAPACK
    /* Pack A into BLAS band storage (kl+ku+1 x d, col-major) for dgbmv; the
     * implicit solve reuses it to seed the dgbtrf factor band. */
    if (op->banded) {
        int ldab = kl + ku + 1;
        double* AB = calloc((size_t)ldab * d, sizeof(double));
        if (AB) {
            for (size_t j = 0; j < d; j++) {
                size_t i0 = (j > (size_t)ku) ? j - (size_t)ku : 0;
                size_t i1 = (j + (size_t)kl < d - 1) ? j + (size_t)kl : d - 1;
                for (size_t i = i0; i <= i1; i++)
                    AB[(size_t)(ku + (long)i - (long)j) + j * (size_t)ldab] = A[i * d + j];
            }
            op->AB = AB;
        }
    }
#endif
    if (const_forcing) {
        op->s0 = s0; op->st = NULL; op->time_forcing = false;
        for (size_t i = 0; i < d; i++) expr_free(st[i]);
        free(st);
    } else {
        free(s0);
        op->s0 = NULL; op->st = st; op->time_forcing = true;
    }
    return op;
}
