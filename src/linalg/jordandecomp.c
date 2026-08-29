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
#include "ndarray.h"
#include "numarray.h"
#ifdef USE_LAPACK
#include "lapack.h"
#endif

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

/* (re, im) of a numeric leaf, including Complex[a, b] and the bare imaginary
 * unit `I`; false if not a concrete number. */
static bool jd_entry_complex(const Expr* e, double* re, double* im) {
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Complex
        && e->data.function.arg_count == 2) {
        return jd_real_double(e->data.function.args[0], re)
            && jd_real_double(e->data.function.args[1], im);
    }
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_I) {
        *re = 0.0; *im = 1.0; return true;
    }
    *im = 0.0;
    return jd_real_double(e, re);
}

/* True iff the eigenvalue (re, im) arrays are pairwise distinct to a relative
 * tolerance.  Distinct spectrum <=> diagonalizable, which validates the fast
 * path with no inverse check. */
static bool jd_distinct(const double* re, const double* im, int n) {
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            double dr = re[i] - re[j], di = im[i] - im[j];
            double gap = sqrt(dr * dr + di * di);
            double sc = 1.0 + sqrt(re[i]*re[i] + im[i]*im[i]);
            if (gap <= 1e-10 * sc) return false;   /* a repeated eigenvalue */
        }
    return true;
}

/* Read an n×n numeric matrix into row-major (re, im) double buffers; false if
 * any entry is not a concrete number. */
static bool jd_read_matrix(Expr* m, int n, double* re, double* im) {
    for (int i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        if (row->type != EXPR_FUNCTION || (int)row->data.function.arg_count != n)
            return false;
        for (int j = 0; j < n; j++)
            if (!jd_entry_complex(row->data.function.args[j],
                                  &re[i*n + j], &im[i*n + j])) return false;
    }
    return true;
}

/* Eigenvalue of eigenvector `v` recovered from `m` (given as re/im buffers) as
 * the component ratio (m.v)_p / v_p at v's largest-magnitude component p.  This
 * is O(n) — it replaces a whole second eigensolve, since v already came from
 * one.  Writes (*lr, *li) and returns the value as an owned Expr (chopped so a
 * real eigenvalue stays Real), or NULL on failure. */
static Expr* jd_recover_lambda(const double* mre, const double* mim, int n,
                               Expr* v, double* lr, double* li) {
    if (v->type != EXPR_FUNCTION || (int)v->data.function.arg_count != n) return NULL;
    double* vre = (double*)malloc(sizeof(double) * (size_t)n);
    double* vim = (double*)malloc(sizeof(double) * (size_t)n);
    bool ok = true;
    for (int j = 0; j < n; j++)
        if (!jd_entry_complex(v->data.function.args[j], &vre[j], &vim[j])) { ok = false; break; }
    int p = 0; double best = -1.0;
    if (ok) for (int j = 0; j < n; j++) {
        double a = vre[j]*vre[j] + vim[j]*vim[j];
        if (a > best) { best = a; p = j; }
    }
    if (!ok || best <= 0.0) { free(vre); free(vim); return NULL; }
    double wr = 0.0, wi = 0.0;                       /* (m.v)_p */
    for (int j = 0; j < n; j++) {
        double ar = mre[p*n + j], ai = mim[p*n + j], br = vre[j], bi = vim[j];
        wr += ar*br - ai*bi;
        wi += ar*bi + ai*br;
    }
    double dr = vre[p], di = vim[p], den = dr*dr + di*di;
    double xr = (wr*dr + wi*di) / den, xi = (wi*dr - wr*di) / den;   /* w / v_p */
    free(vre); free(vim);
    double th = 1e-12 * (fabs(xr) + fabs(xi)) + 1e-14;
    if (fabs(xi) < th) xi = 0.0;
    if (fabs(xr) < th) xr = 0.0;
    *lr = xr; *li = xi;
    if (xi == 0.0) return expr_new_real(xr);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_new_real(xr),
                   expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_real(xi), expr_new_symbol(SYM_I) }, 2) }, 2));
}

/* Numeric fast path: when the spectrum is distinct the matrix is
 * diagonalizable, so s = the numeric eigenvectors as columns and
 * j = DiagonalMatrix[eigenvalues].  A repeated numeric eigenvalue (defective,
 * or ambiguous at machine precision) makes us decline (return NULL) so the
 * caller rationalizes and takes the exact core.  No inverse / product is
 * formed, so this stays fast even with complex eigenvalues (Mathilda has no
 * packed inverse for those).
 *
 * Eigensolve fusion: at machine precision the one `Eigenvectors` solve is
 * enough — each eigenvalue is recovered from its eigenvector by an O(n)
 * component ratio, so the second (`Eigenvalues`) solve is dropped (~24% of the
 * numeric cost at n=200).  Arbitrary-precision input keeps the MPFR
 * `Eigenvalues` solve so the eigenvalues carry full precision (those matrices
 * are small, so the second solve is cheap). */
static Expr* jd_numeric_fast(Expr* m, int n) {
    Expr* evecs_e = eval_and_free(expr_new_function(expr_new_symbol("Eigenvectors"),
                     (Expr*[]){ expr_copy(m) }, 1));
    Expr* evecs = jd_as_boxed(evecs_e);
    expr_free(evecs_e);
    if (!evecs || evecs->type != EXPR_FUNCTION
        || evecs->data.function.head->type != EXPR_SYMBOL
        || evecs->data.function.head->data.symbol.name != SYM_List
        || (int)evecs->data.function.arg_count != n) {
        if (evecs) expr_free(evecs);
        return NULL;
    }

    Expr** lam = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    double* lre = (double*)malloc(sizeof(double) * (size_t)n);
    double* lim = (double*)malloc(sizeof(double) * (size_t)n);
    for (int k = 0; k < n; k++) lam[k] = NULL;
    bool ok = true;

    CommonInexactInfo info = common_scan_inexact(m);
    if (info.min_bits <= 53) {
        double* mre = (double*)malloc(sizeof(double) * (size_t)n * (size_t)n);
        double* mim = (double*)malloc(sizeof(double) * (size_t)n * (size_t)n);
        Expr* m_num = NULL;
        /* Direct read for a concrete-numeric matrix; a symbolic-but-numeric
         * entry (Pi, E, ...) needs one N[] pass, matching what the eigensolver
         * did internally to produce the eigenvectors. */
        if (!jd_read_matrix(m, n, mre, mim)) {
            Expr* nn = eval_and_free(expr_new_function(expr_new_symbol("N"),
                        (Expr*[]){ expr_copy(m) }, 1));
            m_num = jd_as_boxed(nn);
            expr_free(nn);
            if (!m_num || m_num->type != EXPR_FUNCTION
                || (int)m_num->data.function.arg_count != n
                || !jd_read_matrix(m_num, n, mre, mim)) ok = false;
        }
        for (int k = 0; k < n && ok; k++) {
            lam[k] = jd_recover_lambda(mre, mim, n, evecs->data.function.args[k],
                                       &lre[k], &lim[k]);
            if (!lam[k]) ok = false;
        }
        if (m_num) expr_free(m_num);
        free(mre); free(mim);
    } else {
        Expr* evals_e = eval_and_free(expr_new_function(expr_new_symbol("Eigenvalues"),
                         (Expr*[]){ expr_copy(m) }, 1));
        Expr* evals = jd_as_boxed(evals_e);
        expr_free(evals_e);
        if (!evals || evals->type != EXPR_FUNCTION
            || evals->data.function.head->type != EXPR_SYMBOL
            || evals->data.function.head->data.symbol.name != SYM_List
            || (int)evals->data.function.arg_count != n) {
            ok = false;
        } else {
            for (int k = 0; k < n && ok; k++) {
                lam[k] = expr_copy(evals->data.function.args[k]);
                if (!jd_entry_complex(evals->data.function.args[k], &lre[k], &lim[k]))
                    ok = false;
            }
        }
        if (evals) expr_free(evals);
    }

    if (ok) ok = jd_distinct(lre, lim, n);
    free(lre); free(lim);

    if (!ok) {
        for (int k = 0; k < n; k++) if (lam[k]) expr_free(lam[k]);
        free(lam); expr_free(evecs);
        return NULL;
    }

    Expr** cols = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    int* bs = (int*)calloc((size_t)n, sizeof(int));
    for (int k = 0; k < n; k++) { cols[k] = expr_copy(evecs->data.function.args[k]); bs[k] = 1; }
    expr_free(evecs);
    Expr* S = jd_build_S(cols, n);
    Expr* J = jd_build_J(bs, lam, n, n);
    for (int k = 0; k < n; k++) { expr_free(cols[k]); expr_free(lam[k]); }
    free(cols); free(lam); free(bs);

    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ S, J }, 2);
}

/* ------------------------------------------------------------------------ *
 *  LAPACK buffer fast path (the packed / NDArray numeric kernel).           *
 *                                                                           *
 *  Reads the matrix into a column-major double buffer -- straight off the   *
 *  NDArray's float64 payload when the argument is packed, with no delist to  *
 *  boxed Exprs -- calls LAPACK dgeev ONCE (eigenvalues + right eigenvectors),*
 *  and builds s (eigenvectors as columns) and j directly from the raw       *
 *  buffers.  This skips the whole chain the boxed path pays: the delist, the *
 *  intermediate boxed eigenvector list that Eigenvectors materialises, and a *
 *  second (Eigenvalues) solve.  Real matrices only (dgeev); a complex-entry  *
 *  matrix, a repeated spectrum, or a build without LAPACK returns NULL and   *
 *  the caller falls back to the boxed numeric / exact paths.                 *
 * ------------------------------------------------------------------------ */

/* Load an n×n REAL numeric matrix (packed NDArray or boxed List) into a
 * column-major double buffer Ac (Ac[i + j*n] = m[i][j]).  Returns false for a
 * complex, non-concrete, or wrong-shaped matrix. */
static bool jd_load_real_colmajor(Expr* m, int n, double* Ac) {
    if (is_ndarray(m)) {
        const NDArrayData* nd = &m->data.ndarray;
        if (nd->rank != 2 || nd->dims[0] != n || nd->dims[1] != n) return false;
        if (ndt_is_complex(nd->dtype)) return false;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double re, im;
                ndt_get(nd->data, (size_t)i * (size_t)n + (size_t)j, nd->dtype, &re, &im);
                Ac[i + j * n] = re;
            }
        return true;
    }
    if (m->type != EXPR_FUNCTION || (int)m->data.function.arg_count != n) return false;
    for (int i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        if (row->type != EXPR_FUNCTION || (int)row->data.function.arg_count != n) return false;
        for (int j = 0; j < n; j++) {
            double re, im;
            if (!jd_entry_complex(row->data.function.args[j], &re, &im) || im != 0.0)
                return false;
            Ac[i + j * n] = re;
        }
    }
    return true;
}

/* Build s (n×n) from column buffers: column c is the eigenvector
 * (Vr[c*n+i], Vi[c*n+i]); entries with zero imaginary part stay Real. */
static Expr* jd_matrix_from_columns(const double* Vr, const double* Vi, int n) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) {
        Expr** cells = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int c = 0; c < n; c++) {
            double r = Vr[c * n + i], im = Vi[c * n + i];
            if (im == 0.0) cells[c] = expr_new_real(r);
            else {
                Expr* a[2] = { expr_new_real(r), expr_new_real(im) };
                cells[c] = expr_new_function(expr_new_symbol(SYM_Complex), a, 2);
            }
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), cells, (size_t)n);
        free(cells);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    free(rows);
    return out;
}

/* Build the diagonal j (n×n) from eigenvalue buffers (er, ei), chopping tiny
 * noise so a real eigenvalue stays Real; off-diagonal is exact 0. */
static Expr* jd_diag_from_buffers(const double* er, const double* ei, int n) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) {
        Expr** cells = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int j = 0; j < n; j++) {
            if (i != j) { cells[j] = expr_new_integer(0); continue; }
            double r = er[i], im = ei[i];
            double th = 1e-12 * (fabs(r) + fabs(im)) + 1e-14;
            if (fabs(im) < th) im = 0.0;
            if (fabs(r) < th) r = 0.0;
            if (im == 0.0) cells[j] = expr_new_real(r);
            else {
                Expr* a[2] = { expr_new_real(r), expr_new_real(im) };
                cells[j] = expr_new_function(expr_new_symbol(SYM_Complex), a, 2);
            }
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), cells, (size_t)n);
        free(cells);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    free(rows);
    return out;
}

static Expr* jd_numeric_lapack(Expr* m, int n) {
#ifdef USE_LAPACK
    if (!mathilda_lapack_probe()) return NULL;
    double* Ac = (double*)malloc(sizeof(double) * (size_t)n * (size_t)n);
    if (!Ac) return NULL;
    if (!jd_load_real_colmajor(m, n, Ac)) { free(Ac); return NULL; }

    double* wr = (double*)malloc(sizeof(double) * (size_t)n);
    double* wi = (double*)malloc(sizeof(double) * (size_t)n);
    double* VR = (double*)malloc(sizeof(double) * (size_t)n * (size_t)n);
    Expr* result = NULL;
    if (wr && wi && VR
        && mat_lapack_dgeev(n, Ac, n, wr, wi, VR, n) == 0) {
        size_t* perm = (size_t*)malloc(sizeof(size_t) * (size_t)n);
        double* er = (double*)malloc(sizeof(double) * (size_t)n);
        double* ei = (double*)malloc(sizeof(double) * (size_t)n);
        direct_sort_perm_desc_abs_complex(wr, wi, n, perm);
        for (int s = 0; s < n; s++) { er[s] = wr[perm[s]]; ei[s] = wi[perm[s]]; }
        /* Distinct spectrum <=> diagonalizable; else hand back to the exact core. */
        if (jd_distinct(er, ei, n)) {
            double* Vr = (double*)malloc(sizeof(double) * (size_t)n * (size_t)n);
            double* Vi = (double*)malloc(sizeof(double) * (size_t)n * (size_t)n);
            for (int s = 0; s < n; s++) {
                size_t jj = perm[s];
                double norm2 = 0.0;
                for (int i = 0; i < n; i++) {
                    double re, im;               /* dgeev packs conjugate pairs */
                    if (wi[jj] == 0.0)     { re = VR[i + jj * n];       im = 0.0; }
                    else if (wi[jj] > 0.0) { re = VR[i + jj * n];       im = VR[i + (jj + 1) * n]; }
                    else                   { re = VR[i + (jj - 1) * n]; im = -VR[i + jj * n]; }
                    Vr[s * n + i] = re; Vi[s * n + i] = im;
                    norm2 += re * re + im * im;
                }
                double inv = (norm2 > 0.0) ? 1.0 / sqrt(norm2) : 1.0;
                for (int i = 0; i < n; i++) { Vr[s * n + i] *= inv; Vi[s * n + i] *= inv; }
            }
            bool real_spec = true;
            for (int s = 0; s < n; s++) if (ei[s] != 0.0) { real_spec = false; break; }
            Expr* S = NULL; Expr* J = NULL;
            if (real_spec && is_ndarray(m)) {
                /* Real spectrum + packed input -> packed s, j inheriting the
                 * input presentation (na_build_matrix; complex stays boxed). */
                double* sb = (double*)malloc((size_t)n * (size_t)n * sizeof(double));
                double* jb = (double*)calloc((size_t)n * (size_t)n, sizeof(double));
                if (sb && jb) {
                    for (int i = 0; i < n; i++)
                        for (int c = 0; c < n; c++) sb[(size_t)i*n + c] = Vr[(size_t)c*n + i];
                    for (int i = 0; i < n; i++) jb[(size_t)i*n + i] = er[i];
                    S = na_build_matrix(sb, n, n, false, /*colmajor=*/false);
                    J = na_build_matrix(jb, n, n, false, /*colmajor=*/false);
                    NDPresentation pres = m->data.ndarray.present_as;
                    if (S && is_ndarray(S)) S->data.ndarray.present_as = pres;
                    if (J && is_ndarray(J)) J->data.ndarray.present_as = pres;
                }
                free(sb); free(jb);
            }
            if (!S) S = jd_matrix_from_columns(Vr, Vi, n);
            if (!J) J = jd_diag_from_buffers(er, ei, n);
            result = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ S, J }, 2);
            free(Vr); free(Vi);
        }
        free(perm); free(er); free(ei);
    }
    free(Ac); free(wr); free(wi); free(VR);
    return result;
#else
    (void)m; (void)n;
    return NULL;
#endif
}

/* Square-matrix order of `arg` (packed NDArray or boxed List); -1 otherwise. */
static int jd_matrix_order(Expr* arg) {
    if (is_ndarray(arg)) {
        const NDArrayData* nd = &arg->data.ndarray;
        if (nd->rank == 2 && nd->dims[0] > 0 && nd->dims[0] == nd->dims[1])
            return (int)nd->dims[0];
        return -1;
    }
    int64_t dims[64];
    int rank = get_tensor_dims(arg, dims);
    if (rank == 2 && dims[0] > 0 && dims[1] > 0 && dims[0] == dims[1])
        return (int)dims[0];
    return -1;
}

/* Inexact iff a machine/arbitrary-precision matrix (float / complex NDArray, or
 * a boxed matrix carrying a Real / MPFR leaf). */
static bool jd_matrix_is_inexact(Expr* arg) {
    if (is_ndarray(arg)) {
        NDType dt = arg->data.ndarray.dtype;
        return dt == NDT_FLOAT64 || dt == NDT_FLOAT32
            || dt == NDT_COMPLEX64 || dt == NDT_COMPLEX32;
    }
    return eigen_matrix_is_inexact(arg);
}

/* ------------------------------------------------------------------------ *
 *  Builtin entry point.                                                     *
 * ------------------------------------------------------------------------ */

Expr* builtin_jordandecomposition(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;

    size_t argc = res->data.function.arg_count;
    if (argc != 1) {
        fprintf(stderr,
                "JordanDecomposition::argx: JordanDecomposition called with "
                "%zu argument%s; 1 argument is expected.\n",
                argc, argc == 1 ? "" : "s");
        return NULL;
    }

    Expr* arg = res->data.function.args[0];
    int n = jd_matrix_order(arg);
    if (n < 0) {
        char* s = expr_to_string_fullform(arg);
        fprintf(stderr,
                "JordanDecomposition::matsq: Argument %s at position 1 is not "
                "a non-empty square matrix.\n", s);
        free(s);
        return NULL;
    }

    bool inexact = jd_matrix_is_inexact(arg);

    /* Packed / NDArray numeric kernel: read the buffer directly, one dgeev.
     * Machine precision only -- dgeev is double, so an arbitrary-precision
     * (MPFR) matrix must keep the full-precision path below. A float / complex
     * NDArray is machine by dtype; a boxed matrix is machine iff every inexact
     * leaf is <= 53 bits. */
    if (inexact) {
        bool machine = is_ndarray(arg) || common_scan_inexact(arg).min_bits <= 53;
        if (machine) {
            Expr* fast = jd_numeric_lapack(arg, n);
            if (fast) return fast;
        }
    }

    /* Everything else needs a boxed matrix: materialise a packed NDArray once. */
    Expr* mb = is_ndarray(arg) ? pack_unpack(arg) : NULL;
    Expr* m = mb ? mb : arg;
    Expr* result = NULL;

    if (!inexact) {
        result = jd_exact_core(m, n);
    } else {
        result = jd_numeric_fast(m, n);
        if (!result) {
            CommonInexactInfo info = common_scan_inexact(m);
            long bits = info.min_bits ? info.min_bits : 53;
            Expr* m_rat = common_rationalize_input(m, bits);
            Expr* exact = jd_exact_core(m_rat, n);
            expr_free(m_rat);
            if (exact) {
                result = common_numericalize_result(exact, bits);
                expr_free(exact);
            }
        }
    }

    if (mb) expr_free(mb);
    return result;
}

void jordandecomp_init(void) {
    symtab_add_builtin("JordanDecomposition", builtin_jordandecomposition);
    symtab_get_def("JordanDecomposition")->attributes |= ATTR_PROTECTED;
}
