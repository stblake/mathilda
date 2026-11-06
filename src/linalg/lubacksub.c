/* src/linalg/lubacksub.c
 *
 * LUBackSubstitution[{lu, p, c}, b] -- solve m . x == b from the triple
 * LUDecomposition[m] returns, for a vector b or a matrix b (one solution
 * column per column of b).
 *
 * lu packs the unit lower-triangular L (strictly below the diagonal) and the
 * upper-triangular U (on and above it); p is the 1-indexed row permutation
 * with m[[p]] == L . U.  So m . x == b is L . U . x == b[[p]]: one forward
 * substitution through L, one back substitution through U -- O(n^2) per
 * right-hand side, against O(n^3) for a fresh elimination.
 *
 * Three kernels, chosen by the precision of lu and b together (as for any
 * arithmetic, a machine number on either side makes the answer machine):
 *
 *   machine   doubles (real, or interleaved complex) read straight off a
 *             packed / NDArray buffer or a boxed list; the answer is packed
 *             and keeps b's presentation;
 *   MPFR      an MpfrLU (ludecomp_mpfr.c) at the minimum precision of the
 *             operands -- the same solver a LinearSolveFunction of an
 *             arbitrary-precision matrix uses;
 *   exact     the substitutions through the evaluator, so exact and symbolic
 *             entries stay exact (no simplification beyond evaluation).
 *
 * Messages: ::matsq when lu is not a non-empty square matrix, ::lslc when b
 * does not have n rows, ::singlu when U has an exact (or machine / MPFR)
 * zero on its diagonal; the call then stays unevaluated.
 *
 * Compile[]: not lowered.  The first argument is the heterogeneous triple
 * {matrix, permutation, condition number}, which no compiled register
 * carries (LUDecomposition, which produces it, is exempt for the same
 * reason).
 *
 * Ownership: standard builtin contract -- `res` is borrowed; every buffer
 * allocated here is freed on every path.
 */

#include "lubacksub.h"
#include "ludecomp_internal.h"
#include "la_dense.h"
#include "numarray.h"
#include "linalg.h"
#include "ndarray.h"
#include "common.h"
#include "eval.h"
#include "expr.h"
#include "print.h"
#include "symtab.h"
#include "sym_names.h"
#include "attr.h"
#include "message.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ---- argument reading --------------------------------------------------- */

static bool lub_is_list(const Expr* e)
{
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/* The permutation p as 1-indexed ints; false unless it is a permutation of
 * 1..n (a List of Integers or an int64 buffer). */
static bool lub_read_perm(const Expr* p, int n, int* out)
{
    if (p->type == EXPR_NDARRAY) {
        if (p->data.ndarray.rank != 1 || p->data.ndarray.dims[0] != n ||
            p->data.ndarray.dtype != NDT_INT64)
            return false;
        for (int i = 0; i < n; i++) {
            Expr* e = ndarray_element_to_expr(p, (size_t)i);
            bool ok = e && e->type == EXPR_INTEGER;
            if (ok) out[i] = (int)e->data.integer;
            expr_free(e);
            if (!ok) return false;
        }
    } else {
        if (!lub_is_list(p) || (int)p->data.function.arg_count != n) return false;
        for (int i = 0; i < n; i++) {
            const Expr* e = p->data.function.args[i];
            if (e->type != EXPR_INTEGER) return false;
            out[i] = (int)e->data.integer;
        }
    }
    char* seen = (char*)calloc((size_t)n, 1);
    bool ok = seen != NULL;
    for (int i = 0; i < n && ok; i++) {
        if (out[i] < 1 || out[i] > n || seen[out[i] - 1]) ok = false;
        else seen[out[i] - 1] = 1;
    }
    free(seen);
    return ok;
}

/* Rows of b, and its column count (0 for a vector); false when b is neither
 * a vector nor a matrix. */
static bool lub_rhs_shape(const Expr* b, int* rows, int* k)
{
    if (b->type == EXPR_NDARRAY) {
        if (b->data.ndarray.rank == 1) { *rows = (int)b->data.ndarray.dims[0]; *k = 0; return true; }
        if (b->data.ndarray.rank == 2) {
            *rows = (int)b->data.ndarray.dims[0];
            *k = (int)b->data.ndarray.dims[1];
            return *k > 0;
        }
        return false;
    }
    if (!lub_is_list(b) || b->data.function.arg_count == 0) return false;
    *rows = (int)b->data.function.arg_count;
    if (la_looks_like_matrix(b)) {
        int r, c;
        if (!la_matrix_dims(b, &r, &c) || c == 0) return false;
        *k = c;
        return true;
    }
    for (int i = 0; i < *rows; i++)
        if (lub_is_list(b->data.function.args[i])) return false;
    *k = 0;
    return true;
}

static void lub_msg_expr(const char* tag, const char* fmt, const Expr* e)
{
    char* s = expr_to_string((Expr*)e);
    mth_message("LUBackSubstitution", tag, fmt, s ? s : "");
    free(s);
}

static void lub_singular(const Expr* lu)
{
    lub_msg_expr("singlu", "The LU decomposition matrix %s is singular.", lu);
}

/* ---- machine kernel ----------------------------------------------------- */

/* Solve in place: x (n x w row-major, real or interleaved complex) holds
 * b[[p]] on entry and the solution on exit.  LU is n x n row-major.  False
 * on a zero pivot. */
static bool lub_solve_double(const double* LU, int n, bool cplx, double* x, int w)
{
    if (!cplx) {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < i; j++) {
                double l = LU[(size_t)i * n + j];
                if (l == 0.0) continue;
                double* xi = x + (size_t)i * w;
                const double* xj = x + (size_t)j * w;
                for (int c = 0; c < w; c++) xi[c] -= l * xj[c];
            }
        for (int i = n - 1; i >= 0; i--) {
            double* xi = x + (size_t)i * w;
            for (int j = i + 1; j < n; j++) {
                double u = LU[(size_t)i * n + j];
                if (u == 0.0) continue;
                const double* xj = x + (size_t)j * w;
                for (int c = 0; c < w; c++) xi[c] -= u * xj[c];
            }
            double d = LU[(size_t)i * n + i];
            if (d == 0.0) return false;
            for (int c = 0; c < w; c++) xi[c] /= d;
        }
        return true;
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < i; j++) {
            double lr = LU[2 * ((size_t)i * n + j)], li = LU[2 * ((size_t)i * n + j) + 1];
            if (lr == 0.0 && li == 0.0) continue;
            double* xi = x + 2 * (size_t)i * w;
            const double* xj = x + 2 * (size_t)j * w;
            for (int c = 0; c < w; c++) {
                double yr = xj[2 * c], yi = xj[2 * c + 1];
                xi[2 * c]     -= lr * yr - li * yi;
                xi[2 * c + 1] -= lr * yi + li * yr;
            }
        }
    for (int i = n - 1; i >= 0; i--) {
        double* xi = x + 2 * (size_t)i * w;
        for (int j = i + 1; j < n; j++) {
            double ur = LU[2 * ((size_t)i * n + j)], ui = LU[2 * ((size_t)i * n + j) + 1];
            if (ur == 0.0 && ui == 0.0) continue;
            const double* xj = x + 2 * (size_t)j * w;
            for (int c = 0; c < w; c++) {
                double yr = xj[2 * c], yi = xj[2 * c + 1];
                xi[2 * c]     -= ur * yr - ui * yi;
                xi[2 * c + 1] -= ur * yi + ui * yr;
            }
        }
        double dr = LU[2 * ((size_t)i * n + i)], di = LU[2 * ((size_t)i * n + i) + 1];
        double den = dr * dr + di * di;
        if (den == 0.0) return false;
        for (int c = 0; c < w; c++) {
            double ar = xi[2 * c], ai = xi[2 * c + 1];
            xi[2 * c]     = (ar * dr + ai * di) / den;
            xi[2 * c + 1] = (ai * dr - ar * di) / den;
        }
    }
    return true;
}

/* NULL to decline (a leaf that is not a number); *singular on a zero pivot. */
static Expr* lub_machine(const Expr* lu, const int* perm, int n, const Expr* b, int k,
                         bool* singular)
{
    *singular = false;
    int r = 0, c = 0;
    double* LU = NULL;
    bool cplx = false;
    if (!na_load_matrix(lu, false, false, &r, &c, &LU)) {
        if (!na_load_matrix(lu, true, false, &r, &c, &LU)) return NULL;
        cplx = true;
    }
    int w = k ? k : 1;
    double* B = NULL;
    int bn = 0, br = 0, bc = 0;
    bool ok = k ? na_load_matrix(b, cplx, false, &br, &bc, &B)
                : na_load_vector(b, cplx, &bn, &B);
    if (!ok && !cplx) {
        /* a complex right-hand side on a real factor */
        free(LU);
        LU = NULL;
        if (!na_load_matrix(lu, true, false, &r, &c, &LU)) return NULL;
        cplx = true;
        ok = k ? na_load_matrix(b, true, false, &br, &bc, &B) : na_load_vector(b, true, &bn, &B);
    }
    if (!ok) { free(LU); return NULL; }
    size_t sz = (cplx ? 2 : 1) * (size_t)n * (size_t)w;
    double* X = (double*)malloc(sizeof(double) * sz);
    if (!X) { free(LU); free(B); return NULL; }
    size_t rowlen = (cplx ? 2 : 1) * (size_t)w;
    for (int i = 0; i < n; i++)
        memcpy(X + (size_t)i * rowlen, B + (size_t)(perm[i] - 1) * rowlen, sizeof(double) * rowlen);
    free(B);
    Expr* out = NULL;
    if (lub_solve_double(LU, n, cplx, X, w)) {
        NDPresentation pres = na_result_presentation(b);
        out = k ? na_build_matrix_as(X, n, k, cplx, false, pres)
                : na_build_vector_as(X, n, cplx, pres);
    } else {
        *singular = true;
    }
    free(X);
    free(LU);
    return out;
}

/* ---- exact / symbolic kernel -------------------------------------------- */

static Expr* lub_ev(Expr* e)
{
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

/* x_i := s - sum_j a_j y_j, evaluated (the terms are borrowed). */
static Expr* lub_axpy_sum(Expr* s, Expr** a, Expr** y, int cnt)
{
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * (size_t)(cnt + 1));
    int t = 0;
    terms[t++] = s;
    for (int j = 0; j < cnt; j++) {
        Expr* f[3] = { expr_new_integer(-1), expr_copy(a[j]), expr_copy(y[j]) };
        terms[t++] = expr_new_function(expr_new_symbol(SYM_Times), f, 3);
    }
    Expr* sum = t == 1 ? terms[0] : lub_ev(expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)t));
    free(terms);
    return sum;
}

static bool lub_is_exact_zero(const Expr* e)
{
    return e->type == EXPR_INTEGER && e->data.integer == 0;
}

static Expr* lub_exact(const Expr* lu, const int* perm, int n, const Expr* b0, int k,
                       bool* singular)
{
    *singular = false;
    for (int i = 0; i < n; i++)
        if (lub_is_exact_zero(la_at(lu, i, i))) { *singular = true; return NULL; }
    Expr* b = la_boxed_copy(b0);
    int w = k ? k : 1;
    Expr** x = (Expr**)malloc(sizeof(Expr*) * (size_t)n * (size_t)w);
    Expr** a = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    Expr** y = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) {
        const Expr* row = b->data.function.args[perm[i] - 1];
        for (int c = 0; c < w; c++)
            x[(size_t)i * w + c] = expr_copy(k ? row->data.function.args[c] : (Expr*)row);
    }
    for (int c = 0; c < w; c++) {
        /* forward: unit L */
        for (int i = 1; i < n; i++) {
            int cnt = 0;
            for (int j = 0; j < i; j++) {
                Expr* l = la_at(lu, i, j);
                if (lub_is_exact_zero(l)) continue;
                a[cnt] = l;
                y[cnt++] = x[(size_t)j * w + c];
            }
            if (cnt) x[(size_t)i * w + c] = lub_axpy_sum(x[(size_t)i * w + c], a, y, cnt);
        }
        /* back: U */
        for (int i = n - 1; i >= 0; i--) {
            int cnt = 0;
            for (int j = i + 1; j < n; j++) {
                Expr* u = la_at(lu, i, j);
                if (lub_is_exact_zero(u)) continue;
                a[cnt] = u;
                y[cnt++] = x[(size_t)j * w + c];
            }
            Expr* s = cnt ? lub_axpy_sum(x[(size_t)i * w + c], a, y, cnt) : x[(size_t)i * w + c];
            Expr* pw[2] = { expr_copy(la_at(lu, i, i)), expr_new_integer(-1) };
            Expr* f[2] = { s, expr_new_function(expr_new_symbol(SYM_Power), pw, 2) };
            x[(size_t)i * w + c] = lub_ev(expr_new_function(expr_new_symbol(SYM_Times), f, 2));
        }
    }
    free(a);
    free(y);
    expr_free(b);
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++)
        rows[i] = k ? expr_new_function(expr_new_symbol(SYM_List), x + (size_t)i * w, (size_t)k)
                    : x[i];
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    free(rows);
    free(x);
    return out;
}

/* ---- entry -------------------------------------------------------------- */

Expr* builtin_lubacksubstitution(Expr* res)
{
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    const Expr* data = res->data.function.args[0];
    const Expr* b = res->data.function.args[1];
    if (!lub_is_list(data) || data->data.function.arg_count != 3) return NULL;
    const Expr* lu = data->data.function.args[0];
    const Expr* p = data->data.function.args[1];

    int n = 0, nc = 0;
    if (!la_matrix_dims(lu, &n, &nc) || n == 0 || n != nc) {
        lub_msg_expr("matsq", "Argument %s at position 1 is not a non-empty square matrix.", lu);
        return NULL;
    }
    int* perm = (int*)malloc(sizeof(int) * (size_t)n);
    if (!perm) return NULL;
    if (!lub_read_perm(p, n, perm)) { free(perm); return NULL; }

    int brows = 0, k = 0;
    if (!lub_rhs_shape(b, &brows, &k)) { free(perm); return NULL; }
    if (brows != n) {
        mth_message("LUBackSubstitution", "lslc",
                    "Coefficient matrix and target vector(s) or matrix do not have the same dimensions.");
        free(perm);
        return NULL;
    }

    long bits_lu = 53, bits_b = 53;
    LaPrecClass cl = la_prec_class(lu, &bits_lu);
    LaPrecClass cb = la_prec_class(b, &bits_b);
    Expr* out = NULL;
    bool singular = false;

    if (cl == LA_MACHINE || cb == LA_MACHINE) {
        out = lub_machine(lu, perm, n, b, k, &singular);
    } else if (cl == LA_MPFR || cb == LA_MPFR) {
        long bits = cl == LA_MPFR ? bits_lu : bits_b;
        if (cl == LA_MPFR && cb == LA_MPFR && bits_b < bits) bits = bits_b;
        Expr* lub = la_boxed_copy(lu);
        Expr* bb = la_boxed_copy(b);
        MpfrLU* f = mpfr_lu_from_parts(lub, perm, n, bits);
        if (f) {
            bool zero = false;
            for (int i = 0; i < n && !zero; i++) {
                const Expr* d = la_at(lub, i, i);
                if (d->type == EXPR_INTEGER ? d->data.integer == 0
                                            : (d->type == EXPR_REAL && d->data.real == 0.0))
                    zero = true;
            }
#ifdef USE_MPFR
            if (!zero) {
                /* an MPFR zero pivot: Det[lu] == 0 to working precision */
                for (int i = 0; i < n && !zero; i++) {
                    const Expr* d = la_at(lub, i, i);
                    if (d->type == EXPR_MPFR && mpfr_zero_p(d->data.mpfr)) zero = true;
                }
            }
#endif
            if (zero) singular = true;
            else out = mpfr_lu_solve(f, bb);
            mpfr_lu_free(f);
        }
        expr_free(lub);
        expr_free(bb);
    }
    if (!out && !singular) {
        Expr* lub = la_boxed_copy(lu);
        out = lub_exact(lub, perm, n, b, k, &singular);
        expr_free(lub);
    }
    free(perm);
    if (singular) {
        lub_singular(lu);
        expr_free(out);
        return NULL;
    }
    return out;
}

void lubacksub_init(void)
{
    symtab_add_builtin("LUBackSubstitution", builtin_lubacksubstitution);
    symtab_get_def("LUBackSubstitution")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("LUBackSubstitution",
        "LUBackSubstitution[{lu, p, c}, b]\n"
        "\tSolves m . x == b from the LU decomposition {lu, p, c} that\n"
        "\tLUDecomposition[m] returns: forward substitution through the unit\n"
        "\tlower triangle of lu, back substitution through its upper triangle,\n"
        "\tapplied to b[[p]]. b may be a vector or a matrix (one solution per\n"
        "\tcolumn). Machine, arbitrary-precision and exact or symbolic entries\n"
        "\tare supported; the answer has the precision of lu and b together.");
}
