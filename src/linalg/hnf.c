/*
 * hnf.c -- Hermite Normal Form over Z, and the HermiteDecomposition builtin.
 *
 * `linalg_hnf` computes a unimodular P and row-HNF R with P*A == R, tracking
 * every integer row operation on P.  The elimination in each column uses the
 * extended-gcd 2x2 unimodular transform
 *
 *     [ s   t ] [row_r]      [ g*... ]           s*a_r + t*a_i = g
 *     [-b   a ] [row_i]  ->  [   0   ]  in col c, a = a_r/g, b = a_i/g,
 *
 * whose determinant is s*a + t*b = (s*a_r + t*a_i)/g = 1, so P stays
 * unimodular.  Pivots are then made positive and entries above each pivot are
 * reduced into [0, pivot).  This is the reusable integer primitive behind
 * exact linear Diophantine system solving (src/solve/solveint_linear.c).
 */
#include "hnf.h"
#include "linalg.h"
#include "ndlinalg.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"
#include "sym_names.h"
#include "eval.h"
#include "message.h"   /* mth_message: Quiet/Check funnel */

#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef USE_FLINT
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#endif

/* ------------------------------------------------------------------ *
 *  Row HNF core.                                                     *
 * ------------------------------------------------------------------ */

static mpz_t* hnf_mat_alloc(int count) {
    mpz_t* v = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)count);
    if (!v) return NULL;
    for (int i = 0; i < count; i++) mpz_init(v[i]);
    return v;
}

void linalg_hnf_free(mpz_t* mat, int count) {
    if (!mat) return;
    for (int i = 0; i < count; i++) mpz_clear(mat[i]);
    free(mat);
}

/* [Row_r; Row_i] <- [[s, t], [-b, a]] [Row_r; Row_i] over `width` columns.
 * `nr`, `ni`, `tmp` are caller-owned scratch. */
static void hnf_row_combine(mpz_t* Row_r, mpz_t* Row_i, int width,
                            const mpz_t s, const mpz_t t,
                            const mpz_t a, const mpz_t b,
                            mpz_t nr, mpz_t ni, mpz_t tmp) {
    for (int j = 0; j < width; j++) {
        mpz_mul(nr, s, Row_r[j]);
        mpz_mul(tmp, t, Row_i[j]);
        mpz_add(nr, nr, tmp);            /* nr = s*Row_r + t*Row_i         */
        mpz_mul(ni, a, Row_i[j]);
        mpz_mul(tmp, b, Row_r[j]);
        mpz_sub(ni, ni, tmp);            /* ni = a*Row_i - b*Row_r         */
        mpz_set(Row_r[j], nr);
        mpz_set(Row_i[j], ni);
    }
}

#ifdef USE_FLINT
/* linalg_hnf through fmpz_mat_hnf_transform (H = U A). Returns the rank, or
 * -1 on allocation failure (the caller then runs the classical path). */
static int linalg_hnf_flint(const mpz_t* A, int m, int n, mpz_t** R_out, mpz_t** P_out) {
    mpz_t* R = hnf_mat_alloc(m * n);
    mpz_t* P = hnf_mat_alloc(m * m);
    if (!R || !P) { linalg_hnf_free(R, m * n); linalg_hnf_free(P, m * m); return -1; }
    fmpz_mat_t FA, FH, FU;
    fmpz_mat_init(FA, m, n);
    fmpz_mat_init(FH, m, n);
    fmpz_mat_init(FU, m, m);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) fmpz_set_mpz(fmpz_mat_entry(FA, i, j), A[i * n + j]);
    fmpz_mat_hnf_transform(FH, FU, FA);
    int rank = 0;
    for (int i = 0; i < m; i++) {
        int nz = 0;
        for (int j = 0; j < n; j++) {
            fmpz_get_mpz(R[i * n + j], fmpz_mat_entry(FH, i, j));
            if (mpz_sgn(R[i * n + j]) != 0) nz = 1;
        }
        rank += nz;
        for (int j = 0; j < m; j++) fmpz_get_mpz(P[i * m + j], fmpz_mat_entry(FU, i, j));
    }
    fmpz_mat_clear(FA); fmpz_mat_clear(FH); fmpz_mat_clear(FU);
    *R_out = R;
    *P_out = P;
    return rank;
}
#endif

int linalg_hnf(const mpz_t* A, int m, int n, mpz_t** R_out, mpz_t** P_out) {
    if (R_out) *R_out = NULL;
    if (P_out) *P_out = NULL;
    if (m <= 0 || n <= 0) return -1;

#ifdef USE_FLINT
    /* Past a small size, FLINT's fmpz_mat_hnf_transform: the gcd elimination
     * below never reduces as it goes, so its intermediate entries explode (a
     * 40 x 40 matrix of 3-digit entries took seconds). The HNF is unique, so R
     * is the same matrix; P is one valid unimodular transform. Small inputs keep
     * the classical path and its transform. */
    if ((int64_t)m * n >= 64) {
        int r = linalg_hnf_flint(A, m, n, R_out, P_out);
        if (r >= 0) return r;
    }
#endif

    mpz_t* R = hnf_mat_alloc(m * n);
    mpz_t* P = hnf_mat_alloc(m * m);
    if (!R || !P) { linalg_hnf_free(R, m * n); linalg_hnf_free(P, m * m); return -1; }

    for (int i = 0; i < m * n; i++) mpz_set(R[i], A[i]);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++) mpz_set_ui(P[i * m + j], (i == j) ? 1u : 0u);

    mpz_t g, s, t, a, b, q, nr, ni, tmp;
    mpz_inits(g, s, t, a, b, q, nr, ni, tmp, NULL);

    int r = 0;                                   /* next pivot row */
    for (int c = 0; c < n && r < m; c++) {
        /* Locate a pivot at or below row r in column c. */
        int p = -1;
        for (int i = r; i < m; i++)
            if (mpz_sgn(R[i * n + c]) != 0) { p = i; break; }
        if (p < 0) continue;                     /* free column */

        if (p != r) {                            /* bring pivot to row r */
            for (int j = 0; j < n; j++) mpz_swap(R[r * n + j], R[p * n + j]);
            for (int j = 0; j < m; j++) mpz_swap(P[r * m + j], P[p * m + j]);
        }

        /* Eliminate every entry below the pivot in column c via xgcd. */
        for (int i = r + 1; i < m; i++) {
            if (mpz_sgn(R[i * n + c]) == 0) continue;
            mpz_gcdext(g, s, t, R[r * n + c], R[i * n + c]);   /* s*Rr + t*Ri = g */
            mpz_divexact(a, R[r * n + c], g);
            mpz_divexact(b, R[i * n + c], g);
            hnf_row_combine(R + (size_t)r * n, R + (size_t)i * n, n, s, t, a, b, nr, ni, tmp);
            hnf_row_combine(P + (size_t)r * m, P + (size_t)i * m, m, s, t, a, b, nr, ni, tmp);
        }

        /* Normalise the pivot to be positive. */
        if (mpz_sgn(R[r * n + c]) < 0) {
            for (int j = 0; j < n; j++) mpz_neg(R[r * n + j], R[r * n + j]);
            for (int j = 0; j < m; j++) mpz_neg(P[r * m + j], P[r * m + j]);
        }

        /* Reduce entries above the pivot into [0, pivot). */
        for (int i = 0; i < r; i++) {
            if (mpz_sgn(R[i * n + c]) == 0) continue;
            mpz_fdiv_q(q, R[i * n + c], R[r * n + c]);        /* floor -> remainder >= 0 */
            if (mpz_sgn(q) == 0) continue;
            for (int j = 0; j < n; j++) mpz_submul(R[i * n + j], q, R[r * n + j]);
            for (int j = 0; j < m; j++) mpz_submul(P[i * m + j], q, P[r * m + j]);
        }

        r++;
    }

    mpz_clears(g, s, t, a, b, q, nr, ni, tmp, NULL);
    *R_out = R;
    *P_out = P;
    return r;                                    /* rank */
}

/* ------------------------------------------------------------------ *
 *  HermiteDecomposition[A] builtin.                                  *
 * ------------------------------------------------------------------ */

/* Build a matrix Expr (List of Lists) from a flat row-major mpz array. */
static Expr* hnf_mat_to_expr(const mpz_t* M, int rows, int cols) {
    Expr** rowv = (Expr**)malloc(sizeof(Expr*) * (size_t)rows);
    for (int i = 0; i < rows; i++) {
        Expr** cells = (Expr**)malloc(sizeof(Expr*) * (size_t)cols);
        for (int j = 0; j < cols; j++)
            cells[j] = expr_bigint_normalize(expr_new_bigint_from_mpz(M[i * cols + j]));
        rowv[i] = expr_new_function(expr_new_symbol(SYM_List), cells, (size_t)cols);
        free(cells);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rowv, (size_t)rows);
    free(rowv);
    return out;
}

Expr* builtin_hermite_decomposition(Expr* res) {
    if (linalg_call_has_ndarray(res)) return linalg_delist_and_reeval(res);
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) {
        mth_message("HermiteDecomposition", "argx", "HermiteDecomposition called with %zu arguments; 1 argument is expected.", argc);
        return NULL;
    }

    Expr* mexpr = res->data.function.args[0];
    int64_t dims[64];
    int rank = get_tensor_dims(mexpr, dims);
    if (rank != 2 || dims[0] == 0 || dims[1] == 0) {
        char* s = expr_to_string(mexpr);
        mth_message("HermiteDecomposition", "matrix", "Argument %s at position 1 is not a non-empty rectangular matrix.", s);
        free(s);
        return NULL;
    }
    int m = (int)dims[0];
    int n = (int)dims[1];

    Expr** flat = (Expr**)malloc(sizeof(Expr*) * (size_t)m * (size_t)n);
    { size_t idx = 0; flatten_tensor(mexpr, flat, &idx); }

    mpz_t* A = hnf_mat_alloc(m * n);
    bool ok = (A != NULL);
    for (int i = 0; ok && i < m * n; i++) {
        if (!expr_is_integer_like(flat[i])) ok = false;
        else expr_to_mpz(flat[i], A[i]);
    }
    for (int i = 0; i < m * n; i++) expr_free(flat[i]);
    free(flat);

    if (!ok) {
        mth_message("HermiteDecomposition", "intm", "The matrix must have integer entries.");
        linalg_hnf_free(A, m * n);
        return NULL;
    }

    mpz_t* R = NULL; mpz_t* P = NULL;
    int rr = linalg_hnf(A, m, n, &R, &P);
    linalg_hnf_free(A, m * n);
    if (rr < 0) { linalg_hnf_free(R, m * n); linalg_hnf_free(P, m * m); return NULL; }

    Expr* U = hnf_mat_to_expr(P, m, m);
    Expr* Rm = hnf_mat_to_expr(R, m, n);
    linalg_hnf_free(R, m * n);
    linalg_hnf_free(P, m * m);

    Expr* out = expr_new_function(expr_new_symbol(SYM_List),
                                  (Expr*[]){ U, Rm }, 2);
    return out;
}
