/* Mathilda -- sparse linear algebra: matrix predicates in O(nnz).
 *
 * SquareMatrixQ, SymmetricMatrixQ, HermitianMatrixQ, DiagonalMatrixQ,
 * UpperTriangularMatrixQ, LowerTriangularMatrixQ and PositiveDefiniteMatrixQ
 * of a SparseArray, answered from the pattern and the stored values without
 * densifying (the gate densified them: O(n^2), and refused past the 2^27
 * densify cap).
 *
 * The answers are those of the dense builtins (src/list/matrixq.c,
 * src/linalg/posdef_q.c) on Normal[s]: structural equality of mirrored
 * entries (SymmetricMatrixQ), equality with the conjugate (HermitianMatrixQ),
 * exact zeros off the diagonal / outside the triangle. A stored entry of a
 * SparseArray is never equal to its background (construction drops those), so
 * with a zero background "zero" means "not stored". Anything with an option, a
 * diagonal offset k, a nonzero background, or boxed values where the dense
 * test evaluates (HermitianMatrixQ's Conjugate) is declined to the gate.
 *
 * PositiveDefiniteMatrixQ of a machine SparseArray above order 128 tests the
 * Hermitian part (s + s^H)/2 -- the quadratic form of a non-Hermitian matrix --
 * with a Cholesky attempt (CHOLMOD when linked, else the dense potrf within
 * the densify cap); smaller ones keep the dense path's eigenvalue test.
 */
#include "sparse_la.h"
#include "sparse_internal.h"
#include "symtab.h"
#include "eval.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

static Expr* bool_expr(bool b) { return expr_new_symbol(b ? SYM_True : SYM_False); }

/* The single sparse matrix argument with a zero background and no other
 * arguments; NULL to decline. */
static const SparseArray* only_matrix(const Expr* call) {
    if (call->data.function.arg_count != 1) return NULL;
    const SparseArray* sa = sp_of(call->data.function.args[0]);
    if (!sa || sa->rank != 2 || !sp_is_zero(sa->bg)) return NULL;
    return sa;
}

static int64_t col_of(const SparsePattern* p, int64_t k) { return sp_col(p, k, 0); }

/* Index of the stored entry (i, j), or -1. */
static int64_t find_ij(const SparsePattern* p, int64_t i, int64_t j) {
    int64_t lo = p->rowptr[i], hi = p->rowptr[i + 1];
    while (lo < hi) {
        int64_t mid = lo + (hi - lo) / 2;
        if (col_of(p, mid) < j) lo = mid + 1;
        else hi = mid;
    }
    return (lo < p->rowptr[i + 1] && col_of(p, lo) == j) ? lo : -1;
}

/* mode 0: v(k) equal to v(t); mode 1: v(k) equal to conj(v(t)). False for a
 * comparison this file does not make (boxed value under mode 1). */
static bool vals_match(const SparseArray* sa, int64_t k, int64_t t, int mode, bool* known) {
    *known = true;
    switch (sa->kind) {
        case SPV_INT64: return ((const int64_t*)sa->vals)[k] == ((const int64_t*)sa->vals)[t];
        case SPV_F64:   return ((const double*)sa->vals)[k] == ((const double*)sa->vals)[t];
        case SPV_C128: {
            const double* v = sa->vals;
            return v[2 * k] == v[2 * t] && v[2 * k + 1] == (mode ? -v[2 * t + 1] : v[2 * t + 1]);
        }
        default:
            if (mode) { *known = false; return false; }
            return expr_eq(((Expr* const*)sa->vals)[k], ((Expr* const*)sa->vals)[t]);
    }
}

/* SymmetricMatrixQ (mode 0) / HermitianMatrixQ (mode 1); -1 to decline. */
static int mirror_test(const SparseArray* sa, int mode) {
    if (sa->dims[0] != sa->dims[1]) return 0;
    const SparsePattern* p = sa->pat;
    for (int64_t i = 0; i < sa->dims[0]; i++)
        for (int64_t k = p->rowptr[i]; k < p->rowptr[i + 1]; k++) {
            int64_t j = col_of(p, k);
            int64_t t = find_ij(p, j, i);
            if (t < 0) return 0;
            bool known;
            bool ok = vals_match(sa, k, t, mode, &known);
            if (!known) return -1;
            if (!ok) return 0;
        }
    return 1;
}

static Expr* nat_symmetric(const Expr* call) {
    const SparseArray* sa = only_matrix(call);
    if (!sa) return NULL;
    int r = mirror_test(sa, 0);
    return r < 0 ? NULL : bool_expr(r == 1);
}

static Expr* nat_hermitian(const Expr* call) {
    const SparseArray* sa = only_matrix(call);
    if (!sa) return NULL;
    int r = mirror_test(sa, 1);
    return r < 0 ? NULL : bool_expr(r == 1);
}

static Expr* nat_square(const Expr* call) {
    if (call->data.function.arg_count != 1) return NULL;
    const SparseArray* sa = sp_of(call->data.function.args[0]);
    if (!sa) return NULL;
    return bool_expr(sa->rank == 2 && sa->dims[0] == sa->dims[1]);
}

/* which: 0 diagonal, 1 upper, 2 lower -- no stored entry outside the band. */
static Expr* band_test(const Expr* call, int which) {
    const SparseArray* sa = only_matrix(call);
    if (!sa) return NULL;
    const SparsePattern* p = sa->pat;
    for (int64_t i = 0; i < sa->dims[0]; i++)
        for (int64_t k = p->rowptr[i]; k < p->rowptr[i + 1]; k++) {
            int64_t j = col_of(p, k);
            if ((which == 0 && j != i) || (which == 1 && j < i) || (which == 2 && j > i))
                return bool_expr(false);
        }
    return bool_expr(true);
}
static Expr* nat_diagonal(const Expr* call) { return band_test(call, 0); }
static Expr* nat_upper(const Expr* call) { return band_test(call, 1); }
static Expr* nat_lower(const Expr* call) { return band_test(call, 2); }

static Expr* nat_posdef(const Expr* call) {
    const SparseArray* sa = only_matrix(call);
    if (!sa || sa->dims[0] != sa->dims[1] || sa->dims[0] <= 128) return NULL;
    if (sa->kind != SPV_F64 && sa->kind != SPV_C128) return NULL;
    /* the Hermitian part, by the native sparse arithmetic */
    Expr* s = expr_copy((Expr*)call->data.function.args[0]);
    Expr* ct_args[1] = { expr_copy(s) };
    Expr* ct = expr_new_function(expr_new_symbol("ConjugateTranspose"), ct_args, 1);
    Expr* pl_args[2] = { s, ct };
    Expr* sum = expr_new_function(expr_new_symbol("Plus"), pl_args, 2);
    Expr* h_args[2] = { expr_new_real(0.5), sum };
    Expr* hc = expr_new_function(expr_new_symbol("Times"), h_args, 2);
    Expr* h = evaluate(hc);
    expr_free(hc);
    SpMat A;
    bool ok = h && h->type == EXPR_SPARSE && spla_view(sp_of(h), true, &A);
    Expr* r = NULL;
    if (ok) {
        SpfStatus st;
        SpFactor* f = spf_factor(&A, SPF_WANT_CHOL, &st);
        if (f) r = bool_expr(true);
        else if (st == SPF_NOTPD) r = bool_expr(false);
        spf_free(f);
        spla_view_free(&A);
    }
    expr_free(h);
    return r;
}

void sparse_la_pred_init(void) {
    symtab_set_sparse_native("SymmetricMatrixQ", (const void*)nat_symmetric);
    symtab_set_sparse_native("HermitianMatrixQ", (const void*)nat_hermitian);
    symtab_set_sparse_native("SquareMatrixQ", (const void*)nat_square);
    symtab_set_sparse_native("DiagonalMatrixQ", (const void*)nat_diagonal);
    symtab_set_sparse_native("UpperTriangularMatrixQ", (const void*)nat_upper);
    symtab_set_sparse_native("LowerTriangularMatrixQ", (const void*)nat_lower);
    symtab_set_sparse_native("PositiveDefiniteMatrixQ", (const void*)nat_posdef);
}
