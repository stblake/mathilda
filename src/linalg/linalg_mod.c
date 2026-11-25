/* linalg_mod.c - exact linear algebra over Z/p: the Modulus -> p option.
 *
 *   Det[m, Modulus -> p]        Det[m] reduced mod p (any p >= 2)
 *   Inverse[m, Modulus -> p]    the inverse over Z/p; Inverse::sing if none
 *   NullSpace[m, Modulus -> p]  a basis of the kernel over Z/p
 *   RowReduce[m, Modulus -> p]  the reduced row echelon form over Z/p
 *
 * Entries must be integers or rationals whose denominator is invertible mod p;
 * results are residues in [0, p). Modulus -> 0 is the ordinary computation and
 * never reaches this file. The elimination is Gauss-Jordan on int64 residues
 * (p <= 2^31 - 1, so a product of two residues fits). A composite p is a ring,
 * not a field: a pivot must be a unit, and a column whose only nonzero entries
 * are zero divisors leaves the call unevaluated rather than guessing.
 *
 * Det does not eliminate at all: the exact determinant of an integer (or
 * rational) matrix reduced mod p is the determinant over Z/p for every p, prime
 * or not, and the exact Det already has a polynomial-time FLINT path.
 *
 * Memory (SPEC section 4): every argument is borrowed; results are fresh.
 */

#include "linalg_mod.h"
#include "linalg.h"
#include "eval.h"
#include "message.h"
#include "print.h"
#include "sym_names.h"
#include "common.h"
#include <gmp.h>
#include <stdbool.h>
#include <stdlib.h>

int linalg_mod_parse_option(const Expr* opt, int64_t* p) {
    if (!opt || opt->type != EXPR_FUNCTION || opt->data.function.arg_count != 2) return 0;
    if (!head_is(opt, SYM_Rule) && !head_is(opt, SYM_RuleDelayed)) return 0;
    const Expr* k = opt->data.function.args[0];
    if (k->type != EXPR_SYMBOL || k->data.symbol.name != SYM_Modulus) return 0;
    const Expr* v = opt->data.function.args[1];
    if (v->type != EXPR_INTEGER || v->data.integer < 0 || v->data.integer == 1 ||
        v->data.integer > LINALG_MOD_MAX)
        return -1;
    *p = v->data.integer;
    return 1;
}

/* a mod p in [0, p). */
static int64_t mod_norm(int64_t a, int64_t p) {
    int64_t r = a % p;
    return r < 0 ? r + p : r;
}

/* The inverse of a mod p, or -1 when gcd(a, p) != 1. */
static int64_t mod_inv(int64_t a, int64_t p) {
    int64_t t = 0, nt = 1, r = p, nr = mod_norm(a, p);
    while (nr != 0) {
        int64_t q = r / nr, tmp;
        tmp = t - q * nt; t = nt; nt = tmp;
        tmp = r - q * nr; r = nr; nr = tmp;
    }
    if (r != 1) return -1;
    return mod_norm(t, p);
}

/* An integer / BigInt / Rational entry as a residue; false otherwise. */
static bool mod_entry(const Expr* e, int64_t p, int64_t* out) {
    if (e->type == EXPR_INTEGER) { *out = mod_norm(e->data.integer, p); return true; }
    if (e->type == EXPR_BIGINT) {
        *out = (int64_t)mpz_fdiv_ui(e->data.bigint, (unsigned long)p);
        return true;
    }
    if (head_is(e, SYM_Rational) && e->data.function.arg_count == 2) {
        int64_t a, b;
        if (!mod_entry(e->data.function.args[0], p, &a) ||
            !mod_entry(e->data.function.args[1], p, &b)) return false;
        int64_t bi = mod_inv(b, p);
        if (bi < 0) return false;
        *out = (a * bi) % p;
        return true;
    }
    return false;
}

/* m as an r x c row-major residue array (malloc'd), or NULL. */
static int64_t* mod_read(Expr* m, int64_t p, int* r, int* c) {
    int64_t dims[64];
    if (get_tensor_dims(m, dims) != 2 || dims[0] == 0 || dims[1] == 0) return NULL;
    *r = (int)dims[0];
    *c = (int)dims[1];
    int64_t* a = malloc(sizeof(int64_t) * (size_t)(*r) * (size_t)(*c));
    if (!a) return NULL;
    for (int i = 0; i < *r; i++) {
        const Expr* row = m->data.function.args[i];
        for (int j = 0; j < *c; j++)
            if (!mod_entry(row->data.function.args[j], p, &a[(size_t)i * (size_t)(*c) + (size_t)j])) {
                free(a);
                return NULL;
            }
    }
    return a;
}

/* In-place Gauss-Jordan RREF of a[r][c] mod p over the first `ncols`
 * columns. Fills piv[k] = pivot column of row k and returns the rank, or -1
 * when a column has nonzero entries but no unit among them (composite p). */
static int mod_rref(int64_t* a, int r, int c, int ncols, int64_t p, int* piv) {
    int rank = 0;
    for (int col = 0; col < ncols && rank < r; col++) {
        int pr = -1;
        bool nonzero = false;
        int64_t inv = -1;
        for (int i = rank; i < r; i++) {
            int64_t v = a[(size_t)i * c + col];
            if (v == 0) continue;
            nonzero = true;
            inv = mod_inv(v, p);
            if (inv >= 0) { pr = i; break; }
        }
        if (pr < 0) {
            if (nonzero) return -1;
            continue;
        }
        if (pr != rank)
            for (int j = 0; j < c; j++) {
                int64_t t = a[(size_t)pr * c + j];
                a[(size_t)pr * c + j] = a[(size_t)rank * c + j];
                a[(size_t)rank * c + j] = t;
            }
        for (int j = 0; j < c; j++)
            a[(size_t)rank * c + j] = (a[(size_t)rank * c + j] * inv) % p;
        for (int i = 0; i < r; i++) {
            if (i == rank) continue;
            int64_t f = a[(size_t)i * c + col];
            if (f == 0) continue;
            for (int j = 0; j < c; j++)
                a[(size_t)i * c + j] =
                    mod_norm(a[(size_t)i * c + j] - (f * a[(size_t)rank * c + j]) % p, p);
        }
        piv[rank++] = col;
    }
    return rank;
}

/* Rows [r0, r0+nr) x columns [c0, c0+nc) of a (row stride c) as a List matrix. */
static Expr* mod_to_matrix(const int64_t* a, int c, int r0, int nr, int c0, int nc) {
    Expr** rows = malloc(sizeof(Expr*) * (size_t)(nr > 0 ? nr : 1));
    for (int i = 0; i < nr; i++) {
        Expr** el = malloc(sizeof(Expr*) * (size_t)(nc > 0 ? nc : 1));
        for (int j = 0; j < nc; j++)
            el[j] = expr_new_integer(a[(size_t)(r0 + i) * c + (c0 + j)]);
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), el, (size_t)nc);
        free(el);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)nr);
    free(rows);
    return out;
}

Expr* linalg_mod_det(Expr* m, int64_t p) {
    int64_t dims[64];
    if (get_tensor_dims(m, dims) != 2 || dims[0] != dims[1] || dims[0] == 0) return NULL;
    int r, c;
    int64_t* a = mod_read(m, p, &r, &c);       /* entry check only */
    if (!a) return NULL;
    free(a);
    Expr* arg = expr_copy(m);
    Expr* d = eval_and_free(expr_new_function(expr_new_symbol(SYM_Det), &arg, 1));
    int64_t v;
    bool ok = mod_entry(d, p, &v);
    expr_free(d);
    return ok ? expr_new_integer(v) : NULL;
}

Expr* linalg_mod_rowreduce(Expr* m, int64_t p, const char* head) {
    (void)head;
    int r, c;
    int64_t* a = mod_read(m, p, &r, &c);
    if (!a) return NULL;
    int* piv = malloc(sizeof(int) * (size_t)r);
    int rank = piv ? mod_rref(a, r, c, c, p, piv) : -1;
    Expr* out = (rank >= 0) ? mod_to_matrix(a, c, 0, r, 0, c) : NULL;
    free(piv);
    free(a);
    return out;
}

Expr* linalg_mod_inverse(Expr* m, int64_t p) {
    int r, c;
    int64_t* a0 = mod_read(m, p, &r, &c);
    if (!a0) return NULL;
    if (r != c) { free(a0); return NULL; }
    int n = r, w = 2 * n;
    int64_t* a = calloc((size_t)n * (size_t)w, sizeof(int64_t));
    int* piv = malloc(sizeof(int) * (size_t)n);
    if (!a || !piv) { free(a0); free(a); free(piv); return NULL; }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) a[(size_t)i * w + j] = a0[(size_t)i * n + j];
        a[(size_t)i * w + n + i] = 1;
    }
    free(a0);
    int rank = mod_rref(a, n, w, n, p, piv);
    Expr* out = NULL;
    if (rank == n) {
        out = mod_to_matrix(a, w, 0, n, n, n);
    } else if (rank >= 0) {
        char* s = expr_to_string(m);
        mth_message("Inverse", "sing", "Matrix %s is singular.", s);
        free(s);
    }
    free(a); free(piv);
    return out;
}

Expr* linalg_mod_nullspace(Expr* m, int64_t p) {
    int r, c;
    int64_t* a = mod_read(m, p, &r, &c);
    if (!a) return NULL;
    int* piv = malloc(sizeof(int) * (size_t)(r > 0 ? r : 1));
    char* is_piv = calloc((size_t)c, 1);
    int rank = (piv && is_piv) ? mod_rref(a, r, c, c, p, piv) : -1;
    Expr* out = NULL;
    if (rank >= 0) {
        for (int k = 0; k < rank; k++) is_piv[piv[k]] = 1;
        int nfree = c - rank;
        Expr** vecs = malloc(sizeof(Expr*) * (size_t)(nfree > 0 ? nfree : 1));
        int nv = 0;
        /* One basis vector per free column, the last free column first (the
         * order of the exact NullSpace). */
        for (int f = c - 1; f >= 0; f--) {
            if (is_piv[f]) continue;
            Expr** el = malloc(sizeof(Expr*) * (size_t)c);
            for (int j = 0; j < c; j++) el[j] = NULL;
            el[f] = expr_new_integer(1);
            for (int k = 0; k < rank; k++)
                el[piv[k]] = expr_new_integer(mod_norm(-a[(size_t)k * c + f], p));
            for (int j = 0; j < c; j++) if (!el[j]) el[j] = expr_new_integer(0);
            vecs[nv++] = expr_new_function(expr_new_symbol(SYM_List), el, (size_t)c);
            free(el);
        }
        out = expr_new_function(expr_new_symbol(SYM_List), vecs, (size_t)nv);
        free(vecs);
    }
    free(piv); free(is_piv); free(a);
    return out;
}
