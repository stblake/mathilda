/*
 * flint_mat_bridge.c
 * ------------------
 * FLINT-backed exact dense linear algebra over integer / rational matrices.
 * See flint_mat_bridge.h. Mirrors the Expr<->fmpq conversion pattern of
 * src/poly/flint_bridge.c. Everything is behind USE_FLINT; the no-FLINT build
 * provides stubs so the module links either way.
 */

#include "flint_mat_bridge.h"
#include "linalg.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

#ifdef USE_FLINT

#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>

/* ------------------------------------------------------------------ */
/*  Expr <-> fmpq (integer / rational scalars only)                    */
/* ------------------------------------------------------------------ */

/* Fill an fmpz from an integer-like Expr (EXPR_INTEGER or EXPR_BIGINT). */
static int mb_fmpz_from_int_expr(fmpz_t out, const Expr* e) {
    if (e->type == EXPR_INTEGER) { fmpz_set_si(out, e->data.integer); return 1; }
    if (e->type == EXPR_BIGINT)  { fmpz_set_mpz(out, e->data.bigint); return 1; }
    return 0;
}

/* Fill `out` from a scalar Expr that is Integer, BigInt, or Rational[n, d].
 * Returns 0 (leaving the matrix to the classical path) for anything else. */
static int mb_expr_to_fmpq(const Expr* e, fmpq_t out) {
    if (e->type == EXPR_INTEGER) { fmpq_set_si(out, e->data.integer, 1); return 1; }
    if (e->type == EXPR_BIGINT) {
        fmpz_set_mpz(fmpq_numref(out), e->data.bigint);
        fmpz_one(fmpq_denref(out));
        return 1;
    }
    if (e->type == EXPR_FUNCTION) {
        const Expr* h = e->data.function.head;
        if (h && h->type == EXPR_SYMBOL
            && strcmp(h->data.symbol, "Rational") == 0
            && e->data.function.arg_count == 2) {
            if (!mb_fmpz_from_int_expr(fmpq_numref(out), e->data.function.args[0]) ||
                !mb_fmpz_from_int_expr(fmpq_denref(out), e->data.function.args[1]))
                return 0;
            fmpq_canonicalise(out);
            return 1;
        }
    }
    return 0;
}

static Expr* mb_expr_from_mpz(const mpz_t z) {
    return expr_bigint_normalize(expr_new_bigint_from_mpz(z));
}

static Expr* mb_expr_from_fmpq(const fmpq_t q) {
    mpz_t nz, dz;
    mpz_init(nz);
    mpz_init(dz);
    fmpz_get_mpz(nz, fmpq_numref(q));
    fmpz_get_mpz(dz, fmpq_denref(q));
    Expr* r;
    if (mpz_cmp_si(dz, 1) == 0) {
        r = mb_expr_from_mpz(nz);
    } else {
        Expr* args[2] = { mb_expr_from_mpz(nz), mb_expr_from_mpz(dz) };
        r = expr_new_function(expr_new_symbol("Rational"), args, 2);
    }
    mpz_clear(nz);
    mpz_clear(dz);
    return r;
}

/* Fill an r×c fmpq_mat from a row-major flat[] Expr array. Returns 1 on
 * success, 0 as soon as any entry is not integer/rational (matrix left to the
 * classical path). Does not take ownership of `flat`. */
static int mb_fill_mat(fmpq_mat_t M, Expr** flat, int r, int c) {
    for (int i = 0; i < r * c; i++)
        if (!mb_expr_to_fmpq(flat[i], fmpq_mat_entry(M, i / c, i % c))) return 0;
    return 1;
}

/* r×c fmpq_mat -> List-of-Lists Expr (caller owns). */
static Expr* mb_mat_to_expr(const fmpq_mat_t M, int r, int c) {
    Expr** rows = malloc(sizeof(Expr*) * (size_t)r);
    for (int i = 0; i < r; i++) {
        Expr** el = malloc(sizeof(Expr*) * (size_t)c);
        for (int j = 0; j < c; j++)
            el[j] = mb_expr_from_fmpq(fmpq_mat_entry(M, i, j));
        rows[i] = expr_new_function(expr_new_symbol("List"), el, (size_t)c);
        free(el);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), rows, (size_t)r);
    free(rows);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Determinant                                                        */
/* ------------------------------------------------------------------ */

Expr* flint_mat_det(Expr** flat, int n) {
    if (n <= 0) return NULL;

    fmpq_mat_t M;
    fmpq_mat_init(M, n, n);

    int ok = 1;
    for (int i = 0; i < n * n && ok; i++) {
        ok = mb_expr_to_fmpq(flat[i], fmpq_mat_entry(M, i / n, i % n));
    }

    Expr* out = NULL;
    if (ok) {
        fmpq_t d; fmpq_init(d);
        fmpq_mat_det(d, M);
        out = mb_expr_from_fmpq(d);
        fmpq_clear(d);
    }

    fmpq_mat_clear(M);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Inverse / RowReduce / LinearSolve / MatrixRank                     */
/* ------------------------------------------------------------------ */

Expr* flint_mat_inverse(Expr** flat, int n) {
    if (n <= 0) return NULL;
    fmpq_mat_t A, B;
    fmpq_mat_init(A, n, n);
    fmpq_mat_init(B, n, n);
    Expr* out = NULL;
    if (mb_fill_mat(A, flat, n, n) && fmpq_mat_inv(B, A))
        out = mb_mat_to_expr(B, n, n);
    fmpq_mat_clear(A);
    fmpq_mat_clear(B);
    return out;  /* NULL: non-rational entry OR singular -> classical path */
}

Expr* flint_mat_rref(Expr** flat, int r, int c) {
    if (r <= 0 || c <= 0) return NULL;
    fmpq_mat_t A, B;
    fmpq_mat_init(A, r, c);
    fmpq_mat_init(B, r, c);
    Expr* out = NULL;
    if (mb_fill_mat(A, flat, r, c)) {
        fmpq_mat_rref(B, A);
        out = mb_mat_to_expr(B, r, c);
    }
    fmpq_mat_clear(A);
    fmpq_mat_clear(B);
    return out;
}

Expr** flint_mat_solve(Expr** flat_m, int r, int c, Expr** flat_b, int k) {
    if (r <= 0 || c <= 0 || k <= 0 || r != c) return NULL;
    fmpq_mat_t A, B, X;
    fmpq_mat_init(A, r, c);
    fmpq_mat_init(B, r, k);
    fmpq_mat_init(X, c, k);
    Expr** out = NULL;
    if (mb_fill_mat(A, flat_m, r, c) && mb_fill_mat(B, flat_b, r, k)
        && fmpq_mat_solve(X, A, B)) {
        out = malloc(sizeof(Expr*) * (size_t)c * (size_t)k);
        for (int i = 0; i < c * k; i++)
            out[i] = mb_expr_from_fmpq(fmpq_mat_entry(X, i / k, i % k));
    }
    fmpq_mat_clear(A);
    fmpq_mat_clear(B);
    fmpq_mat_clear(X);
    return out;  /* NULL: non-rational / non-square / singular -> classical */
}

int flint_mat_rank(Expr** flat, int r, int c) {
    if (r <= 0 || c <= 0) return -1;
    fmpq_mat_t A, B;
    fmpq_mat_init(A, r, c);
    fmpq_mat_init(B, r, c);
    int rank = -1;
    if (mb_fill_mat(A, flat, r, c))
        rank = (int)fmpq_mat_rref(B, A);
    fmpq_mat_clear(A);
    fmpq_mat_clear(B);
    return rank;
}

/* ------------------------------------------------------------------ */
/*  FLINT` context builtins                                            */
/* ------------------------------------------------------------------ */

/* Shared: validate a rank-2 tensor argument and flatten it row-major.
 * On success sets *out_flat (malloc'd, caller frees each + array), *out_r,
 * *out_c and returns 1; returns 0 (nothing allocated) otherwise. */
static int mb_take_matrix(Expr* arg, Expr*** out_flat, int* out_r, int* out_c) {
    int64_t dims[64];
    if (get_tensor_dims(arg, dims) != 2 || dims[0] == 0 || dims[1] == 0) return 0;
    int r = (int)dims[0], c = (int)dims[1];
    Expr** flat = malloc(sizeof(Expr*) * (size_t)r * (size_t)c);
    size_t idx = 0;
    flatten_tensor(arg, flat, &idx);
    *out_flat = flat; *out_r = r; *out_c = c;
    return 1;
}

static void mb_free_flat(Expr** flat, int count) {
    for (int i = 0; i < count; i++) expr_free(flat[i]);
    free(flat);
}

static Expr* builtin_flint_det(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    int64_t dims[64];
    int rank = get_tensor_dims(arg, dims);
    if (rank != 2 || dims[0] != dims[1] || dims[0] == 0) return NULL;

    int n = (int)dims[0];
    Expr** flat = malloc(sizeof(Expr*) * (size_t)n * (size_t)n);
    size_t idx = 0;
    flatten_tensor(arg, flat, &idx);

    Expr* out = flint_mat_det(flat, n);

    for (size_t i = 0; i < idx; i++) expr_free(flat[i]);
    free(flat);
    return out;  /* NULL (non-rational entry) -> stays unevaluated */
}

static Expr* builtin_flint_inverse(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr** flat; int r, c;
    if (!mb_take_matrix(res->data.function.args[0], &flat, &r, &c)) return NULL;
    Expr* out = (r == c) ? flint_mat_inverse(flat, r) : NULL;
    mb_free_flat(flat, r * c);
    return out;
}

static Expr* builtin_flint_rowreduce(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr** flat; int r, c;
    if (!mb_take_matrix(res->data.function.args[0], &flat, &r, &c)) return NULL;
    Expr* out = flint_mat_rref(flat, r, c);
    mb_free_flat(flat, r * c);
    return out;
}

static Expr* builtin_flint_matrixrank(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr** flat; int r, c;
    if (!mb_take_matrix(res->data.function.args[0], &flat, &r, &c)) return NULL;
    int rank = flint_mat_rank(flat, r, c);
    mb_free_flat(flat, r * c);
    return (rank < 0) ? NULL : expr_new_integer(rank);
}

/* FLINT`LinearSolve[m, b]: m an r×c matrix, b a length-r vector (result is a
 * length-c vector) or an r×k matrix (result is c×k). Square unique-solution
 * case only; anything else stays unevaluated. */
static Expr* builtin_flint_linearsolve(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* m = res->data.function.args[0];
    Expr* b = res->data.function.args[1];

    int64_t md[64], bd[64];
    if (get_tensor_dims(m, md) != 2 || md[0] == 0 || md[1] == 0) return NULL;
    int r = (int)md[0], c = (int)md[1];
    int brank = get_tensor_dims(b, bd);
    int k, vector_rhs;
    if (brank == 1 && bd[0] == r)                    { k = 1; vector_rhs = 1; }
    else if (brank == 2 && bd[0] == r && bd[1] != 0) { k = (int)bd[1]; vector_rhs = 0; }
    else return NULL;

    Expr** flat_m = malloc(sizeof(Expr*) * (size_t)r * (size_t)c);
    Expr** flat_b = malloc(sizeof(Expr*) * (size_t)r * (size_t)k);
    { size_t i = 0; flatten_tensor(m, flat_m, &i); }
    { size_t i = 0; flatten_tensor(b, flat_b, &i); }

    Expr** sol = flint_mat_solve(flat_m, r, c, flat_b, k);
    mb_free_flat(flat_m, r * c);
    mb_free_flat(flat_b, r * k);
    if (!sol) return NULL;

    Expr* out;
    if (vector_rhs) {
        Expr** el = malloc(sizeof(Expr*) * (size_t)c);
        for (int j = 0; j < c; j++) el[j] = sol[j];   /* k == 1 */
        out = expr_new_function(expr_new_symbol("List"), el, (size_t)c);
        free(el);
    } else {
        Expr** rows = malloc(sizeof(Expr*) * (size_t)c);
        for (int i = 0; i < c; i++) {
            Expr** el = malloc(sizeof(Expr*) * (size_t)k);
            for (int j = 0; j < k; j++) el[j] = sol[i * k + j];
            rows[i] = expr_new_function(expr_new_symbol("List"), el, (size_t)k);
            free(el);
        }
        out = expr_new_function(expr_new_symbol("List"), rows, (size_t)c);
        free(rows);
    }
    free(sol);
    return out;
}

void flint_mat_bridge_init(void) {
    symtab_add_builtin(SYM_FLINT_Det, builtin_flint_det);
    symtab_get_def(SYM_FLINT_Det)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_Det,
        "FLINT`Det[m] gives the determinant of the square matrix m when every "
        "entry is an integer or rational, computed exactly and directly via "
        "FLINT (fmpq_mat_det). Returns unevaluated for a matrix with any "
        "non-rational entry.");

    symtab_add_builtin(SYM_FLINT_Inverse, builtin_flint_inverse);
    symtab_get_def(SYM_FLINT_Inverse)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_Inverse,
        "FLINT`Inverse[m] gives the inverse of the square matrix m when every "
        "entry is an integer or rational, computed exactly via FLINT "
        "(fmpq_mat_inv). Returns unevaluated for a singular or non-rational "
        "matrix.");

    symtab_add_builtin(SYM_FLINT_RowReduce, builtin_flint_rowreduce);
    symtab_get_def(SYM_FLINT_RowReduce)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_RowReduce,
        "FLINT`RowReduce[m] gives the reduced row echelon form of the matrix m "
        "when every entry is an integer or rational, computed exactly via FLINT "
        "(fmpq_mat_rref). Returns unevaluated for a matrix with any non-rational "
        "entry.");

    symtab_add_builtin(SYM_FLINT_LinearSolve, builtin_flint_linearsolve);
    symtab_get_def(SYM_FLINT_LinearSolve)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_LinearSolve,
        "FLINT`LinearSolve[m, b] solves the square system m.x == b exactly via "
        "FLINT (fmpq_mat_solve) when m is a nonsingular rational matrix and b a "
        "rational vector or matrix. Returns unevaluated for a non-square, "
        "singular, or non-rational system.");

    symtab_add_builtin(SYM_FLINT_MatrixRank, builtin_flint_matrixrank);
    symtab_get_def(SYM_FLINT_MatrixRank)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_MatrixRank,
        "FLINT`MatrixRank[m] gives the rank of the matrix m when every entry is "
        "an integer or rational, computed exactly via FLINT (fmpq_mat_rref). "
        "Returns unevaluated for a matrix with any non-rational entry.");
}

#else /* !USE_FLINT */

Expr* flint_mat_det(Expr** flat, int n) { (void)flat; (void)n; return NULL; }
Expr* flint_mat_inverse(Expr** flat, int n) { (void)flat; (void)n; return NULL; }
Expr* flint_mat_rref(Expr** flat, int r, int c) { (void)flat; (void)r; (void)c; return NULL; }
Expr** flint_mat_solve(Expr** fm, int r, int c, Expr** fb, int k) {
    (void)fm; (void)r; (void)c; (void)fb; (void)k; return NULL;
}
int flint_mat_rank(Expr** flat, int r, int c) { (void)flat; (void)r; (void)c; return -1; }
void  flint_mat_bridge_init(void) { /* no FLINT: nothing to register */ }

#endif /* USE_FLINT */
