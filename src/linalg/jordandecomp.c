/*
 * jordandecomp.c -- JordanDecomposition[m].
 *
 * JordanDecomposition[m] gives {s, j} where j is the Jordan canonical form of
 * the square matrix m and s is the similarity matrix, so that
 *
 *     m == s . j . Inverse[s].
 *
 * j is block-diagonal in Jordan blocks J_k(lambda) (lambda on the diagonal,
 * 1 on the superdiagonal); the columns of s are the generalized eigenvectors
 * grouped into Jordan chains and ordered to match j's blocks.
 *
 * ------------------------------------------------------------------------
 * Two code paths, one chain engine.
 *
 *   Exact / symbolic (jd_exact_core):  the characteristic polynomial
 *   (eigen_char_poly_faddeev) is solved for the distinct eigenvalues with
 *   their algebraic multiplicities.  For each lambda, N = m - lambda*I, and
 *   the nullity sequence nu_i = dim ker(N^i) drives a top-down chain-top
 *   selection / bottom-up chain construction.  Every structural step runs in
 *   the field of the matrix entries (exact rational, or free-symbolic) using
 *   eigen_null_space + a MatrixRank-based "extend a spanning set to a basis"
 *   primitive.  A generalized eigenspace that cannot be spanned (an
 *   irrational *defective* eigenvalue, where is_zero_poly cannot decide the
 *   RowReduce pivots) leaves the whole call unevaluated rather than returning
 *   a wrong answer.
 *
 *   Numeric (machine / MPFR, real / complex):  a generic numeric matrix has
 *   distinct eigenvalues and is diagonalizable, so jd_numeric_fast reads the
 *   eigenvectors straight off the numeric eigensolver (the Eigenvectors head,
 *   LAPACK-style Direct kernel / MPFR twin) as the columns of s and puts the
 *   per-column Rayleigh eigenvalue on the diagonal of j.  This is what every
 *   "diagonal j" example produces, and it scales (a 100x100 random matrix
 *   never touches exact arithmetic).  A numerically *defective* matrix (its
 *   eigenvectors are rank-deficient) falls back to rationalize -> exact core
 *   -> numericalize, which recovers the genuine block structure.
 *
 * This file reuses the internal eigen helpers (eigen_char_poly_faddeev,
 * eigen_solve_poly, eigen_null_space, ...) declared in eigen_internal.h.
 * JordanDecomposition is part of the eigen/linalg cluster, so pulling in that
 * header here is deliberate.
 *
 * Memory contract: standard builtin ownership (SPEC.md §4).  This file does
 * NOT free `res`; it returns a fresh Expr* the evaluator owns, or NULL to
 * leave the call unevaluated.
 */

#include "jordandecomp.h"
#include "eigen.h"
#include "eigen_internal.h"
#include "linalg.h"
#include "ndlinalg.h"
#include "common.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "print.h"
#include "pack.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ------------------------------------------------------------------------ *
 *  Small matrix / vector primitives (evaluate through the ordinary heads so *
 *  every arithmetic field -- exact, rational, symbolic -- is handled).      *
 * ------------------------------------------------------------------------ */

/* A . B  (matrix.matrix or matrix.vector), returned as a nested List.  dot2
 * of a machine operand answers with a packed NDArray; unpack it so the
 * caller can index rows/entries as data.function.args.  Caller owns. */
static Expr* jd_dot(Expr* a, Expr* b) {
    bool err = false;
    Expr* p = dot2(a, b, &err);
    if (!p) return NULL;
    Expr* r = eval_and_free(p);
    Expr* u = pack_unpack(r);
    if (u) { expr_free(r); return u; }
    return r;
}

/* m - lambda*I, entrywise, evaluated.  Returns a fresh n×n nested List. */
static Expr* jd_shift(Expr* m, Expr* lambda, int n) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) {
        Expr* mrow = m->data.function.args[i];
        Expr** cells = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int j = 0; j < n; j++) {
            Expr* mij = mrow->data.function.args[j];
            if (i == j) {
                Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(lambda) }, 2);
                cells[j] = eval_and_free(expr_new_function(
                    expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(mij), neg }, 2));
            } else {
                cells[j] = expr_copy(mij);
            }
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), cells, (size_t)n);
        free(cells);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    free(rows);
    return out;
}

/* Rank of the matrix whose rows are rows[0..k-1] (each a length-n List),
 * via MatrixRank.  Returns -1 if the answer is not a machine Integer. */
static int jd_rank_of_rows(Expr** rows, size_t k) {
    if (k == 0) return 0;
    Expr** rc = (Expr**)malloc(sizeof(Expr*) * k);
    for (size_t i = 0; i < k; i++) rc[i] = expr_copy(rows[i]);
    Expr* mat = expr_new_function(expr_new_symbol(SYM_List), rc, k);
    free(rc);
    Expr* call = expr_new_function(expr_new_symbol("MatrixRank"),
                                   (Expr*[]){ mat }, 1);
    Expr* r = eval_and_free(call);
    int rank = -1;
    if (r->type == EXPR_INTEGER) rank = (int)r->data.integer;
    expr_free(r);
    return rank;
}

/* Select from cand[0..nc-1] the vectors that are linearly independent modulo
 * the span of base[0..nb-1] (greedy, by rank increase).  Returns freshly
 * copied chosen vectors; *out holds the count.  A nullspace basis is already
 * independent, so when base is empty every candidate is taken directly. */
static Expr** jd_extend_basis(Expr** base, size_t nb,
                              Expr** cand, size_t nc, size_t* out) {
    *out = 0;
    if (nc == 0) return NULL;
    if (nb == 0) {
        Expr** chosen = (Expr**)malloc(sizeof(Expr*) * nc);
        for (size_t i = 0; i < nc; i++) chosen[i] = expr_copy(cand[i]);
        *out = nc;
        return chosen;
    }
    Expr** running = (Expr**)malloc(sizeof(Expr*) * (nb + nc));
    size_t nr = 0;
    for (size_t i = 0; i < nb; i++) running[nr++] = expr_copy(base[i]);
    int rrank = jd_rank_of_rows(running, nr);
    Expr** chosen = (Expr**)malloc(sizeof(Expr*) * nc);
    size_t nchosen = 0;
    for (size_t i = 0; i < nc; i++) {
        running[nr] = expr_copy(cand[i]);
        int nrank = jd_rank_of_rows(running, nr + 1);
        if (nrank > rrank) { nr++; rrank = nrank; chosen[nchosen++] = expr_copy(cand[i]); }
        else { expr_free(running[nr]); }
    }
    for (size_t i = 0; i < nr; i++) expr_free(running[i]);
    free(running);
    *out = nchosen;
    return chosen;
}

/* ------------------------------------------------------------------------ *
 *  Result assembly.                                                         *
 * ------------------------------------------------------------------------ */

/* Build the n×n similarity matrix s from its columns cols[0..n-1] (each a
 * length-n List).  s[[i, c]] = cols[c][[i]]. */
static Expr* jd_build_S(Expr** cols, int n) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) {
        Expr** cells = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int c = 0; c < n; c++)
            cells[c] = expr_copy(cols[c]->data.function.args[i]);
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), cells, (size_t)n);
        free(cells);
    }
    Expr* S = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    free(rows);
    return S;
}

/* Build the n×n Jordan matrix from the block list: block b has eigenvalue
 * bl[b] on the diagonal and size bs[b], with 1's on the block superdiagonal. */
static Expr* jd_build_J(const int* bs, Expr** bl, int nblk, int n) {
    Expr** flat = (Expr**)malloc(sizeof(Expr*) * (size_t)n * (size_t)n);
    for (int i = 0; i < n * n; i++) flat[i] = expr_new_integer(0);
    int off = 0;
    for (int b = 0; b < nblk; b++) {
        int s = bs[b];
        for (int k = 0; k < s; k++) {
            int d = (off + k) * n + (off + k);
            expr_free(flat[d]);
            flat[d] = expr_copy(bl[b]);
            if (k < s - 1) {
                int sup = (off + k) * n + (off + k + 1);
                expr_free(flat[sup]);
                flat[sup] = expr_new_integer(1);
            }
        }
        off += s;
    }
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) {
        Expr** cells = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int j = 0; j < n; j++) cells[j] = flat[i * n + j];
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), cells, (size_t)n);
        free(cells);
    }
    free(flat);
    Expr* J = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    free(rows);
    return J;
}

/* ------------------------------------------------------------------------ *
 *  Exact / symbolic eigenvalue collection and ordering.                     *
 * ------------------------------------------------------------------------ */

/* N[lambda] as a double (real part collapse) for ordering.  Returns false
 * when lambda does not reduce to a concrete real number. */
static bool jd_num_key(Expr* lambda, double* out) {
    Expr* nn = eval_and_free(expr_new_function(expr_new_symbol("N"),
                             (Expr*[]){ expr_copy(lambda) }, 1));
    bool ok = true; double d = 0;
    if (nn->type == EXPR_REAL) d = nn->data.real;
    else if (nn->type == EXPR_INTEGER) d = (double)nn->data.integer;
    else if (nn->type == EXPR_BIGINT) d = mpz_get_d(nn->data.bigint);
#ifdef USE_MPFR
    else if (nn->type == EXPR_MPFR) d = mpfr_get_d(nn->data.mpfr, MPFR_RNDN);
#endif
    else ok = false;
    expr_free(nn);
    if (ok) *out = d;
    return ok;
}

/* Distinct eigenvalues (owned) with algebraic multiplicities, sorted
 * ascending by value when all are concrete-real (matches the spec examples).
 * Returns false and leaves outputs untouched on failure. */
static bool jd_collect_eigenvalues(Expr* m, int n,
                                   Expr*** out_lams, int** out_mult, int* out_nd) {
    const char* lam = eigen_lambda_name();
    Expr* poly = eigen_char_poly_faddeev(m, lam, n);
    if (!poly) return false;
    Expr* sols = eigen_solve_poly(poly, lam, false, false);
    expr_free(poly);
    if (!sols) return false;
    size_t vc = 0;
    Expr** vals = eigen_extract_values(sols, &vc);
    expr_free(sols);
    if (!vals || vc == 0) { free(vals); return false; }

    Expr** lams = (Expr**)malloc(sizeof(Expr*) * vc);
    int* mult = (int*)malloc(sizeof(int) * vc);
    int nd = 0;
    for (size_t i = 0; i < vc; i++) {
        int found = -1;
        for (int t = 0; t < nd; t++)
            if (expr_eq(lams[t], vals[i])) { found = t; break; }
        if (found >= 0) { mult[found]++; expr_free(vals[i]); }
        else { lams[nd] = vals[i]; mult[nd] = 1; nd++; }
    }
    free(vals);

    /* Sort ascending by numeric key when every eigenvalue is concrete-real. */
    double* keys = (double*)malloc(sizeof(double) * (size_t)nd);
    bool all_real = true;
    for (int i = 0; i < nd; i++)
        if (!jd_num_key(lams[i], &keys[i])) { all_real = false; break; }
    if (all_real) {
        for (int i = 1; i < nd; i++) {
            Expr* lv = lams[i]; int mv = mult[i]; double kv = keys[i];
            int j = i - 1;
            while (j >= 0 && keys[j] > kv) {
                lams[j + 1] = lams[j]; mult[j + 1] = mult[j]; keys[j + 1] = keys[j];
                j--;
            }
            lams[j + 1] = lv; mult[j + 1] = mv; keys[j + 1] = kv;
        }
    }
    free(keys);

    *out_lams = lams; *out_mult = mult; *out_nd = nd;
    return true;
}

/* ------------------------------------------------------------------------ *
 *  Per-eigenvalue Jordan-chain construction (the core of the exact path).   *
 * ------------------------------------------------------------------------ */

/* Append the Jordan chains for eigenvalue `lambda` (algebraic multiplicity
 * `mult`) to the running column / block accumulators.  Returns false if the
 * generalized eigenspace cannot be spanned exactly (the call must then be
 * left unevaluated). */
static bool jd_process_eigenvalue(Expr* m, Expr* lambda, int mult, int n,
                                  Expr** cols, int* pncols,
                                  int* blk_size, Expr** blk_lam, int* pnblk) {
    Expr* N = jd_shift(m, lambda, n);
    int cap = n + 2;
    Expr*** Bk = (Expr***)calloc((size_t)cap, sizeof(Expr**));  /* ker(N^k) bases */
    size_t*  dk = (size_t*)calloc((size_t)cap, sizeof(size_t)); /* nullities      */
    Expr*** tops  = (Expr***)calloc((size_t)cap, sizeof(Expr**));
    size_t*  ntops = (size_t*)calloc((size_t)cap, sizeof(size_t));

    bool fail = false;
    int p = 0;

    /* Nullity sequence: grow k until dim ker(N^k) reaches the algebraic
     * multiplicity (the generalized eigenspace is then complete). */
    Expr* Nk = expr_copy(N);
    int k = 1;
    size_t prev_d = 0;
    while (1) {
        size_t dc = 0;
        Expr** basis = eigen_null_space(Nk, n, &dc);
        if (dc == 0 && k == 1) {          /* irrational simple root recovery */
            free(basis);
            basis = eigen_null_space_algebraic(m, NULL, lambda, n, &dc);
        }
        Bk[k] = basis; dk[k] = dc;
        if ((int)dc >= mult) { p = k; expr_free(Nk); break; }
        if (dc == prev_d) { p = k; expr_free(Nk); fail = true; break; }  /* stalled */
        prev_d = dc;
        if (k >= n) { p = k; expr_free(Nk); fail = true; break; }
        Expr* Nk1 = jd_dot(Nk, N);
        expr_free(Nk);
        if (!Nk1) { p = k; fail = true; break; }
        Nk = Nk1; k++;
    }
    if (!fail && (int)dk[p] != mult) fail = true;

    /* Chain-top selection: length-p tops extend ker(N^{p-1}) to ker(N^p). */
    if (!fail) {
        Expr** basePm1 = (p >= 2) ? Bk[p - 1] : NULL;
        size_t nbPm1   = (p >= 2) ? dk[p - 1] : 0;
        tops[p] = jd_extend_basis(basePm1, nbPm1, Bk[p], dk[p], &ntops[p]);

        for (int i = p - 1; i >= 1 && !fail; i--) {
            size_t base_n = (i >= 2) ? dk[i - 1] : 0;
            size_t img_n = 0;
            for (int j = i + 1; j <= p; j++) img_n += ntops[j];
            Expr** spanned = (Expr**)malloc(sizeof(Expr*) * (base_n + img_n + 1));
            size_t sc = 0;
            if (i >= 2)
                for (size_t t = 0; t < dk[i - 1]; t++)
                    spanned[sc++] = expr_copy(Bk[i - 1][t]);
            for (int j = i + 1; j <= p && !fail; j++) {
                for (size_t t = 0; t < ntops[j]; t++) {
                    Expr* v = expr_copy(tops[j][t]);
                    for (int s = 0; s < (j - i) && v; s++) {
                        Expr* nv = jd_dot(N, v); expr_free(v); v = nv;
                    }
                    if (!v) { fail = true; break; }
                    spanned[sc++] = v;
                }
            }
            if (!fail)
                tops[i] = jd_extend_basis(spanned, sc, Bk[i], dk[i], &ntops[i]);
            for (size_t t = 0; t < sc; t++) expr_free(spanned[t]);
            free(spanned);
        }
    }

    /* Generate chains: each length-i top w yields columns
     * (N^{i-1} w, ..., N w, w) and one Jordan block J_i(lambda). */
    if (!fail) {
        for (int i = 1; i <= p && !fail; i++) {
            for (size_t t = 0; t < ntops[i]; t++) {
                if (*pncols + i > n) { fail = true; break; }
                Expr** c = (Expr**)malloc(sizeof(Expr*) * (size_t)i);
                c[0] = expr_copy(tops[i][t]);
                bool cok = true;
                for (int s = 1; s < i; s++) {
                    c[s] = jd_dot(N, c[s - 1]);
                    if (!c[s]) { cok = false; break; }
                }
                if (!cok) {
                    for (int s = 0; s < i; s++) if (c[s]) expr_free(c[s]);
                    free(c); fail = true; break;
                }
                for (int s = i - 1; s >= 0; s--) cols[(*pncols)++] = c[s];
                free(c);
                blk_size[*pnblk] = i;
                blk_lam[*pnblk] = expr_copy(lambda);
                (*pnblk)++;
            }
        }
    }

    for (int i = 1; i < cap; i++) {
        if (Bk[i]) { for (size_t t = 0; t < dk[i]; t++) expr_free(Bk[i][t]); free(Bk[i]); }
        if (tops[i]) { for (size_t t = 0; t < ntops[i]; t++) expr_free(tops[i][t]); free(tops[i]); }
    }
    free(Bk); free(dk); free(tops); free(ntops);
    expr_free(N);
    return !fail;
}

/* Exact / symbolic Jordan decomposition of the (exact-entry) matrix m.
 * Returns {s, j}, or NULL if the structure cannot be recovered exactly. */
static Expr* jd_exact_core(Expr* m, int n) {
    Expr** lams; int* mult; int nd;
    if (!jd_collect_eigenvalues(m, n, &lams, &mult, &nd)) return NULL;

    Expr** cols = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    int ncols = 0;
    int* blk_size = (int*)malloc(sizeof(int) * (size_t)n);
    Expr** blk_lam = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    int nblk = 0;
    bool ok = true;
    for (int e = 0; e < nd && ok; e++)
        ok = jd_process_eigenvalue(m, lams[e], mult[e], n,
                                   cols, &ncols, blk_size, blk_lam, &nblk);

    for (int e = 0; e < nd; e++) expr_free(lams[e]);
    free(lams); free(mult);

    if (!ok || ncols != n) {
        for (int i = 0; i < ncols; i++) expr_free(cols[i]);
        for (int b = 0; b < nblk; b++) expr_free(blk_lam[b]);
        free(cols); free(blk_size); free(blk_lam);
        return NULL;
    }

    Expr* S = jd_build_S(cols, n);
    Expr* J = jd_build_J(blk_size, blk_lam, nblk, n);
    for (int i = 0; i < ncols; i++) expr_free(cols[i]);
    for (int b = 0; b < nblk; b++) expr_free(blk_lam[b]);
    free(cols); free(blk_size); free(blk_lam);

    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ S, J }, 2);
}

/* ------------------------------------------------------------------------ *
 *  Numeric diagonalizable fast path.                                        *
 * ------------------------------------------------------------------------ */

/* Owned nested-List view of a possibly-packed matrix. */
static Expr* jd_as_boxed(Expr* e) {
    Expr* u = pack_unpack(e);
    return u ? u : expr_copy(e);
}

/* Concrete real value of a numeric leaf (Real / Integer / BigInt / MPFR /
 * Rational).  Returns false when e is not one of those. */
static bool jd_real_double(const Expr* e, double* out) {
    switch (e->type) {
        case EXPR_REAL:    *out = e->data.real;              return true;
        case EXPR_INTEGER: *out = (double)e->data.integer;   return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        default: break;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rational
        && e->data.function.arg_count == 2) {
        double p, q;
        if (jd_real_double(e->data.function.args[0], &p)
            && jd_real_double(e->data.function.args[1], &q) && q != 0.0) {
            *out = p / q; return true;
        }
    }
    return false;
}

/* (re, im) of a numeric leaf, including Complex[a, b]; false if not numeric. */
static bool jd_entry_complex(const Expr* e, double* re, double* im) {
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Complex
        && e->data.function.arg_count == 2) {
        return jd_real_double(e->data.function.args[0], re)
            && jd_real_double(e->data.function.args[1], im);
    }
    *im = 0.0;
    return jd_real_double(e, re);
}

/* True iff the n eigenvalues in `evals` (a List) are pairwise distinct to a
 * relative tolerance.  Distinct eigenvalues  <=>  the matrix is diagonalizable,
 * which is what makes the fast path valid without any inverse check. */
static bool jd_all_distinct(Expr* evals, int n) {
    double* re = (double*)malloc(sizeof(double) * (size_t)n);
    double* im = (double*)malloc(sizeof(double) * (size_t)n);
    bool ok = true;
    for (int i = 0; i < n && ok; i++)
        if (!jd_entry_complex(evals->data.function.args[i], &re[i], &im[i])) ok = false;
    if (ok) {
        for (int i = 0; i < n && ok; i++)
            for (int j = i + 1; j < n && ok; j++) {
                double dr = re[i] - re[j], di = im[i] - im[j];
                double gap = sqrt(dr * dr + di * di);
                double sc = 1.0 + sqrt(re[i]*re[i] + im[i]*im[i]);
                if (gap <= 1e-10 * sc) ok = false;   /* a repeated eigenvalue */
            }
    }
    free(re); free(im);
    return ok;
}

/* Numeric fast path: when the spectrum is distinct the matrix is
 * diagonalizable, so s = the numeric eigenvectors as columns and
 * j = DiagonalMatrix[eigenvalues] (paired positionally -- Eigenvalues and
 * Eigenvectors share the numeric kernel's deterministic order).  A repeated
 * numeric eigenvalue (defective, or diagonalizable-with-multiplicity) is
 * ambiguous at machine precision, so we decline (return NULL) and let the
 * caller rationalize and take the exact core, which resolves it correctly.
 * No inverse / product is needed here, so this stays fast even for a matrix
 * with complex eigenvalues (which Mathilda has no packed inverse for). */
static Expr* jd_numeric_fast(Expr* m, int n) {
    Expr* evals_e = eval_and_free(expr_new_function(expr_new_symbol("Eigenvalues"),
                     (Expr*[]){ expr_copy(m) }, 1));
    Expr* evals = jd_as_boxed(evals_e);
    expr_free(evals_e);
    if (!evals || evals->type != EXPR_FUNCTION
        || evals->data.function.head->type != EXPR_SYMBOL
        || evals->data.function.head->data.symbol.name != SYM_List
        || (int)evals->data.function.arg_count != n
        || !jd_all_distinct(evals, n)) {
        if (evals) expr_free(evals);
        return NULL;
    }

    Expr* evecs_e = eval_and_free(expr_new_function(expr_new_symbol("Eigenvectors"),
                     (Expr*[]){ expr_copy(m) }, 1));
    Expr* evecs = jd_as_boxed(evecs_e);
    expr_free(evecs_e);
    if (!evecs || evecs->type != EXPR_FUNCTION
        || evecs->data.function.head->type != EXPR_SYMBOL
        || evecs->data.function.head->data.symbol.name != SYM_List
        || (int)evecs->data.function.arg_count != n) {
        if (evecs) expr_free(evecs);
        expr_free(evals);
        return NULL;
    }

    /* s: eigenvectors as columns; j: diagonal of the matching eigenvalues. */
    Expr** cols = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    Expr** lam  = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    int* bs = (int*)calloc((size_t)n, sizeof(int));
    for (int k = 0; k < n; k++) {
        cols[k] = expr_copy(evecs->data.function.args[k]);
        lam[k]  = expr_copy(evals->data.function.args[k]);
        bs[k] = 1;
    }
    expr_free(evecs); expr_free(evals);
    Expr* S = jd_build_S(cols, n);
    Expr* J = jd_build_J(bs, lam, n, n);
    for (int k = 0; k < n; k++) { expr_free(cols[k]); expr_free(lam[k]); }
    free(cols); free(lam); free(bs);

    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ S, J }, 2);
}

/* ------------------------------------------------------------------------ *
 *  Builtin entry point.                                                     *
 * ------------------------------------------------------------------------ */

Expr* builtin_jordandecomposition(Expr* res) {
    if (linalg_call_has_ndarray(res)) return linalg_delist_and_reeval(res);
    if (res->type != EXPR_FUNCTION) return NULL;

    size_t argc = res->data.function.arg_count;
    if (argc != 1) {
        fprintf(stderr,
                "JordanDecomposition::argx: JordanDecomposition called with "
                "%zu argument%s; 1 argument is expected.\n",
                argc, argc == 1 ? "" : "s");
        return NULL;
    }

    Expr* mexpr = res->data.function.args[0];
    int64_t dims[64];
    int rank = get_tensor_dims(mexpr, dims);
    if (rank != 2 || dims[0] == 0 || dims[1] == 0 || dims[0] != dims[1]) {
        char* s = expr_to_string_fullform(mexpr);
        fprintf(stderr,
                "JordanDecomposition::matsq: Argument %s at position 1 is not "
                "a non-empty square matrix.\n", s);
        free(s);
        return NULL;
    }
    int n = (int)dims[0];

    if (!eigen_matrix_is_inexact(mexpr))
        return jd_exact_core(mexpr, n);

    /* Inexact: diagonalizable fast path, else rationalize -> exact -> numericalize. */
    Expr* fast = jd_numeric_fast(mexpr, n);
    if (fast) return fast;

    CommonInexactInfo info = common_scan_inexact(mexpr);
    long bits = info.min_bits ? info.min_bits : 53;
    Expr* m_rat = common_rationalize_input(mexpr, bits);
    Expr* exact = jd_exact_core(m_rat, n);
    expr_free(m_rat);
    if (!exact) return NULL;
    Expr* num = common_numericalize_result(exact, bits);
    expr_free(exact);
    return num;
}

void jordandecomp_init(void) {
    symtab_add_builtin("JordanDecomposition", builtin_jordandecomposition);
    symtab_get_def("JordanDecomposition")->attributes |= ATTR_PROTECTED;
}
