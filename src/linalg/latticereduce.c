/* latticereduce.c
 *
 * LatticeReduce[{v1, v2, ...}] -- an LLL-reduced basis for the lattice
 * spanned by the row vectors v_i.
 *
 *   LatticeReduce[m]   m an n x d matrix (List of n equal-length Lists).
 *                      Returns an n x d matrix whose rows form a reduced
 *                      basis of the same lattice (same Z- / Z[i]-module).
 *
 * Entries may be:
 *   - integers (machine int64 or GMP bigint),
 *   - rationals (Rational[p, q]),
 *   - Gaussian integers / Gaussian rationals (Complex[a, b] with a, b
 *     integer or rational).
 *
 * Algorithm: the classical Lenstra-Lenstra-Lovas (LLL) reduction with
 * Lovas parameter delta = 3/4, run in EXACT arithmetic.  Every quantity
 * is an exact Gaussian rational stored as a pair of GMP `mpq_t`
 * (`GRat`), so the routine is correct for both machine-size and
 * arbitrary-precision (bignum) inputs -- floating point is never used,
 * which matters because LatticeReduce is most often used to discover
 * integer relations where a rounding error would yield a wrong relation.
 *
 * The Gram-Schmidt orthogonalisation is generalised to the complex
 * (Hermitian) inner product  <x, y> = sum_k x_k conj(y_k),  so the same
 * code path handles real and Gaussian lattices.  Size reduction rounds
 * the mu coefficients to the nearest Gaussian integer; the Gram-Schmidt
 * data (mu, |b*|^2) is maintained incrementally -- computed once up
 * front, updated in place on size reduction, and updated on a Lovas swap
 * via the conjugate-aware Cohen swap formulas (no full recomputation).
 *
 * Because every basis transformation is an integer (Z, or Z[i] in the
 * Gaussian case) row operation or a row swap, the lattice -- and hence
 * Abs[Det] and every linear relation in the right null space -- is
 * preserved exactly.
 *
 * Linearly independent rows are required (every documented use of
 * LatticeReduce -- integer-relation finding, basis reduction -- supplies
 * a full-rank generating set).  A dependent generating set is detected
 * during Gram-Schmidt and the call is left unevaluated with a
 * diagnostic.
 *
 * Memory ownership: standard builtin contract.  This file does NOT free
 * `res` on success or failure -- the evaluator owns it (see MEMORY.md /
 * SPEC.md S4.1).
 */

#include "linalg.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"
#include "sym_names.h"
#include "expr.h"

#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ *
 *  Exact Gaussian-rational scalar:  re + im * I,  re, im in Q.         *
 * ------------------------------------------------------------------ */
typedef struct { mpq_t re, im; } GRat;

static void gr_init(GRat* z)  { mpq_init(z->re); mpq_init(z->im); }
static void gr_clear(GRat* z) { mpq_clear(z->re); mpq_clear(z->im); }
static void gr_set(GRat* d, const GRat* s) { mpq_set(d->re, s->re); mpq_set(d->im, s->im); }
static bool gr_is_zero(const GRat* z) { return mpq_sgn(z->re) == 0 && mpq_sgn(z->im) == 0; }

/* d = a + b  (d may alias a and/or b). */
static void gr_add(GRat* d, const GRat* a, const GRat* b) {
    mpq_t r, i; mpq_init(r); mpq_init(i);
    mpq_add(r, a->re, b->re);
    mpq_add(i, a->im, b->im);
    mpq_set(d->re, r); mpq_set(d->im, i);
    mpq_clear(r); mpq_clear(i);
}

/* d = a - b  (d may alias a and/or b). */
static void gr_sub(GRat* d, const GRat* a, const GRat* b) {
    mpq_t r, i; mpq_init(r); mpq_init(i);
    mpq_sub(r, a->re, b->re);
    mpq_sub(i, a->im, b->im);
    mpq_set(d->re, r); mpq_set(d->im, i);
    mpq_clear(r); mpq_clear(i);
}

/* d = a * b  (complex multiply; d may alias a and/or b). */
static void gr_mul(GRat* d, const GRat* a, const GRat* b) {
    mpq_t t1, t2, rr, ii;
    mpq_init(t1); mpq_init(t2); mpq_init(rr); mpq_init(ii);
    mpq_mul(t1, a->re, b->re);
    mpq_mul(t2, a->im, b->im);
    mpq_sub(rr, t1, t2);                 /* re = a.re b.re - a.im b.im */
    mpq_mul(t1, a->re, b->im);
    mpq_mul(t2, a->im, b->re);
    mpq_add(ii, t1, t2);                 /* im = a.re b.im + a.im b.re */
    mpq_set(d->re, rr); mpq_set(d->im, ii);
    mpq_clear(t1); mpq_clear(t2); mpq_clear(rr); mpq_clear(ii);
}

/* d = conj(a)  (d may alias a). */
static void gr_conj(GRat* d, const GRat* a) {
    mpq_set(d->re, a->re);
    mpq_neg(d->im, a->im);
}

/* d = a * s,  s a real rational  (d may alias a). */
static void gr_scale(GRat* d, const GRat* a, const mpq_t s) {
    mpq_mul(d->re, a->re, s);
    mpq_mul(d->im, a->im, s);
}

/* out = |z|^2 = z.re^2 + z.im^2  (real, >= 0). out must be initialised. */
static void gr_norm2(mpq_t out, const GRat* z) {
    mpq_t t; mpq_init(t);
    mpq_mul(out, z->re, z->re);
    mpq_mul(t, z->im, z->im);
    mpq_add(out, out, t);
    mpq_clear(t);
}

/* out = nearest integer to x (ties rounded up), as an integer-valued
 * rational.  out must be initialised. */
static void mpq_round(mpq_t out, const mpq_t x) {
    /* floor(x + 1/2) = floor((2*num + den) / (2*den)),  den > 0. */
    mpz_t n2, d2, q;
    mpz_init(n2); mpz_init(d2); mpz_init(q);
    mpz_mul_ui(n2, mpq_numref(x), 2);
    mpz_add(n2, n2, mpq_denref(x));
    mpz_mul_ui(d2, mpq_denref(x), 2);
    mpz_fdiv_q(q, n2, d2);
    mpq_set_z(out, q);
    mpz_clear(n2); mpz_clear(d2); mpz_clear(q);
}

/* d = nearest Gaussian integer to a  (round re and im independently). */
static void gr_round_gauss(GRat* d, const GRat* a) {
    mpq_round(d->re, a->re);
    mpq_round(d->im, a->im);
}

/* ------------------------------------------------------------------ *
 *  Vectors (row-major GRat blocks) and the Hermitian inner product.   *
 * ------------------------------------------------------------------ */

/* out = <x, y> = sum_k x_k conj(y_k).  out must be initialised. */
static void vec_inner(GRat* out, const GRat* x, const GRat* y, int d) {
    GRat acc, cy, prod;
    gr_init(&acc); gr_init(&cy); gr_init(&prod);
    for (int k = 0; k < d; k++) {
        gr_conj(&cy, &y[k]);
        gr_mul(&prod, &x[k], &cy);
        gr_add(&acc, &acc, &prod);
    }
    gr_set(out, &acc);
    gr_clear(&acc); gr_clear(&cy); gr_clear(&prod);
}

/* ------------------------------------------------------------------ *
 *  Expr <-> exact-number conversions.                                 *
 * ------------------------------------------------------------------ */

/* Parse an exact real rational Expr (Integer / BigInt / Rational of
 * integer-likes) into the pre-initialised `out`.  Returns false for
 * anything else.  On a false return `out` may have been written but is
 * left in a valid (clearable) state. */
static bool expr_to_mpq(Expr* e, mpq_t out) {
    if (expr_is_integer_like(e)) {
        mpz_t z; expr_to_mpz(e, z);
        mpq_set_z(out, z);
        mpz_clear(z);
        return true;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational
        && e->data.function.arg_count == 2) {
        Expr* n = e->data.function.args[0];
        Expr* dn = e->data.function.args[1];
        if (!expr_is_integer_like(n) || !expr_is_integer_like(dn)) return false;
        mpz_t zn, zd; expr_to_mpz(n, zn); expr_to_mpz(dn, zd);
        bool ok = (mpz_sgn(zd) != 0);
        if (ok) {
            mpq_set_num(out, zn);
            mpq_set_den(out, zd);
            mpq_canonicalize(out);
        }
        mpz_clear(zn); mpz_clear(zd);
        return ok;
    }
    return false;
}

/* Convert a matrix entry to a Gaussian rational `z` (already gr_init'd).
 * Accepts integers, rationals, the symbol I, and Complex[re, im] with
 * exact rational parts.  Returns false otherwise. */
static bool entry_to_grat(Expr* e, GRat* z) {
    if (expr_to_mpq(e, z->re)) {           /* real integer / rational */
        mpq_set_ui(z->im, 0, 1);
        return true;
    }
    if (e->type == EXPR_SYMBOL && strcmp(e->data.symbol, "I") == 0) {
        mpq_set_ui(z->re, 0, 1);
        mpq_set_si(z->im, 1, 1);
        return true;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Complex
        && e->data.function.arg_count == 2) {
        if (!expr_to_mpq(e->data.function.args[0], z->re)) return false;
        if (!expr_to_mpq(e->data.function.args[1], z->im)) return false;
        return true;
    }
    return false;
}

/* Build an exact Integer / BigInt / Rational Expr from a rational q. */
static Expr* mpq_to_expr(const mpq_t q) {
    Expr* num = expr_bigint_normalize(expr_new_bigint_from_mpz(mpq_numref(q)));
    if (mpz_cmp_ui(mpq_denref(q), 1) == 0) return num;
    Expr* den = expr_bigint_normalize(expr_new_bigint_from_mpz(mpq_denref(q)));
    Expr* args[2] = { num, den };
    return expr_new_function(expr_new_symbol("Rational"), args, 2);
}

/* Build an exact Integer / Rational / Complex Expr from a GRat. */
static Expr* grat_to_expr(const GRat* z) {
    if (mpq_sgn(z->im) == 0) return mpq_to_expr(z->re);
    Expr* re = mpq_to_expr(z->re);
    Expr* im = mpq_to_expr(z->im);
    Expr* args[2] = { re, im };
    return expr_new_function(expr_new_symbol("Complex"), args, 2);
}

/* ------------------------------------------------------------------ *
 *  LLL core.                                                           *
 *                                                                     *
 *  B      : n x d basis, row-major GRat block.   ROW(B, r) = first    *
 *           entry of row r.                                            *
 *  mu     : n x n GRat block; only mu[i*n + j], j < i, is meaningful.  *
 *  Bnorm  : length-n array of mpq_t, Bnorm[i] = |b*_i|^2 (real, > 0).  *
 * ------------------------------------------------------------------ */
#define ROW(B, r, d) (&(B)[(size_t)(r) * (size_t)(d)])

/* Initial Gram-Schmidt.  Fills mu and Bnorm from B, using a scratch
 * orthogonal basis Bstar.  Returns 0 on success, 1 if the rows are
 * linearly dependent (some |b*_i|^2 == 0). */
static int gso_init(const GRat* B, GRat* Bstar, GRat* mu, mpq_t* Bnorm,
                    int n, int d) {
    GRat cij, prod;
    gr_init(&cij); gr_init(&prod);
    mpq_t invn, t;
    mpq_init(invn); mpq_init(t);

    int dependent = 0;
    for (int i = 0; i < n && !dependent; i++) {
        for (int c = 0; c < d; c++) gr_set(&ROW(Bstar, i, d)[c], &ROW(B, i, d)[c]);
        for (int j = 0; j < i; j++) {
            vec_inner(&cij, ROW(B, i, d), ROW(Bstar, j, d), d);
            mpq_inv(invn, Bnorm[j]);                 /* Bnorm[j] > 0 here */
            gr_scale(&mu[i * n + j], &cij, invn);    /* mu_ij = <b_i,b*_j>/|b*_j|^2 */
            for (int c = 0; c < d; c++) {
                gr_mul(&prod, &mu[i * n + j], &ROW(Bstar, j, d)[c]);
                gr_sub(&ROW(Bstar, i, d)[c], &ROW(Bstar, i, d)[c], &prod);
            }
        }
        mpq_set_ui(Bnorm[i], 0, 1);
        for (int c = 0; c < d; c++) {
            gr_norm2(t, &ROW(Bstar, i, d)[c]);
            mpq_add(Bnorm[i], Bnorm[i], t);
        }
        if (mpq_sgn(Bnorm[i]) == 0) dependent = 1;
    }

    gr_clear(&cij); gr_clear(&prod);
    mpq_clear(invn); mpq_clear(t);
    return dependent;
}

/* Size reduction RED(k, l):  subtract round(mu[k][l]) * b_l from b_k and
 * propagate the change to the mu table. */
static void red(GRat* B, GRat* mu, int k, int l, int n, int d) {
    GRat* mkl = &mu[k * n + l];
    GRat q; gr_init(&q);
    gr_round_gauss(&q, mkl);
    if (!gr_is_zero(&q)) {
        GRat prod; gr_init(&prod);
        for (int c = 0; c < d; c++) {                /* b_k -= q b_l */
            gr_mul(&prod, &q, &ROW(B, l, d)[c]);
            gr_sub(&ROW(B, k, d)[c], &ROW(B, k, d)[c], &prod);
        }
        for (int j = 0; j < l; j++) {                /* mu_kj -= q mu_lj */
            gr_mul(&prod, &q, &mu[l * n + j]);
            gr_sub(&mu[k * n + j], &mu[k * n + j], &prod);
        }
        gr_sub(mkl, mkl, &q);                         /* mu_kl -= q */
        gr_clear(&prod);
    }
    gr_clear(&q);
}

/* Lovas swap of rows k and k-1, with the conjugate-aware incremental
 * Gram-Schmidt update (Cohen, "A Course in Computational Algebraic
 * Number Theory", Alg. 2.6.3, extended to the Hermitian inner product). */
static void swap_rows(GRat* B, GRat* mu, mpq_t* Bnorm, int k, int n, int d) {
    GRat mu_old, conj_old, tnum, tprod;
    gr_init(&mu_old); gr_init(&conj_old); gr_init(&tnum); gr_init(&tprod);
    mpq_t mu2, Bk, Bk1, Bnew, ratio;
    mpq_init(mu2); mpq_init(Bk); mpq_init(Bk1); mpq_init(Bnew); mpq_init(ratio);

    gr_set(&mu_old, &mu[k * n + (k - 1)]);
    gr_norm2(mu2, &mu_old);                          /* |mu|^2 */
    mpq_set(Bk, Bnorm[k]);
    mpq_set(Bk1, Bnorm[k - 1]);
    mpq_mul(Bnew, mu2, Bk1);
    mpq_add(Bnew, Bnew, Bk);                          /* Bnew = Bk + |mu|^2 Bk1 */

    /* new mu[k][k-1] = conj(mu) * Bk1 / Bnew */
    mpq_div(ratio, Bk1, Bnew);
    gr_conj(&conj_old, &mu_old);
    gr_scale(&mu[k * n + (k - 1)], &conj_old, ratio);

    mpq_mul(Bnorm[k], Bk1, Bk);
    mpq_div(Bnorm[k], Bnorm[k], Bnew);               /* Bnorm[k]   = Bk1 Bk / Bnew */
    mpq_set(Bnorm[k - 1], Bnew);                     /* Bnorm[k-1] = Bnew */

    /* swap the (already-computed) mu rows for columns j < k-1 */
    for (int j = 0; j < k - 1; j++) {
        gr_set(&tnum, &mu[(k - 1) * n + j]);
        gr_set(&mu[(k - 1) * n + j], &mu[k * n + j]);
        gr_set(&mu[k * n + j], &tnum);
    }
    /* swap the basis rows themselves */
    for (int c = 0; c < d; c++) {
        gr_set(&tnum, &ROW(B, k, d)[c]);
        gr_set(&ROW(B, k, d)[c], &ROW(B, k - 1, d)[c]);
        gr_set(&ROW(B, k - 1, d)[c], &tnum);
    }
    /* update mu[i][k-1], mu[i][k] for the rows i > k */
    for (int i = k + 1; i < n; i++) {
        gr_set(&tnum, &mu[i * n + k]);                       /* t = old mu[i][k]   */
        gr_mul(&tprod, &mu_old, &tnum);
        gr_sub(&mu[i * n + k], &mu[i * n + (k - 1)], &tprod);/* mu[i][k]  = a - mu t */
        gr_mul(&tprod, &mu[k * n + (k - 1)], &mu[i * n + k]);
        gr_add(&mu[i * n + (k - 1)], &tnum, &tprod);         /* mu[i][k-1]= t + mu_new mu[i][k] */
    }

    gr_clear(&mu_old); gr_clear(&conj_old); gr_clear(&tnum); gr_clear(&tprod);
    mpq_clear(mu2); mpq_clear(Bk); mpq_clear(Bk1); mpq_clear(Bnew); mpq_clear(ratio);
}

/* In-place LLL reduction of the n x d basis B.  Returns 0 on success,
 * 1 if the rows are linearly dependent. */
static int lll_reduce(GRat* B, int n, int d) {
    if (n <= 1) return 0;                             /* already reduced */

    GRat* mu = malloc(sizeof(GRat) * (size_t)n * (size_t)n);
    GRat* Bstar = malloc(sizeof(GRat) * (size_t)n * (size_t)d);
    mpq_t* Bnorm = malloc(sizeof(mpq_t) * (size_t)n);
    for (int i = 0; i < n * n; i++) gr_init(&mu[i]);
    for (int i = 0; i < n * d; i++) gr_init(&Bstar[i]);
    for (int i = 0; i < n; i++) mpq_init(Bnorm[i]);

    int dependent = gso_init(B, Bstar, mu, Bnorm, n, d);

    for (int i = 0; i < n * d; i++) gr_clear(&Bstar[i]);
    free(Bstar);

    if (!dependent) {
        mpq_t delta, mu2, rhs;
        mpq_init(delta); mpq_init(mu2); mpq_init(rhs);
        mpq_set_si(delta, 3, 4);                     /* Lovas parameter */

        int k = 1;
        while (k < n) {
            red(B, mu, k, k - 1, n, d);
            gr_norm2(mu2, &mu[k * n + (k - 1)]);
            mpq_sub(rhs, delta, mu2);
            mpq_mul(rhs, rhs, Bnorm[k - 1]);         /* (delta - |mu|^2) |b*_{k-1}|^2 */
            if (mpq_cmp(Bnorm[k], rhs) >= 0) {       /* Lovas condition holds */
                for (int l = k - 2; l >= 0; l--) red(B, mu, k, l, n, d);
                k++;
            } else {
                swap_rows(B, mu, Bnorm, k, n, d);
                k = (k - 1 >= 1) ? (k - 1) : 1;
            }
        }
        mpq_clear(delta); mpq_clear(mu2); mpq_clear(rhs);
    }

    for (int i = 0; i < n * n; i++) gr_clear(&mu[i]);
    for (int i = 0; i < n; i++) mpq_clear(Bnorm[i]);
    free(mu);
    free(Bnorm);
    return dependent;
}

/* ------------------------------------------------------------------ *
 *  Public builtin.                                                    *
 * ------------------------------------------------------------------ */
Expr* builtin_latticereduce(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) {
        fprintf(stderr,
                "LatticeReduce::argx: LatticeReduce called with %zu "
                "argument%s; 1 argument is expected.\n",
                argc, argc == 1 ? "" : "s");
        return NULL;
    }

    Expr* m = res->data.function.args[0];
    int64_t dims[64];
    int rank = get_tensor_dims(m, dims);
    if (rank != 2 || dims[0] == 0 || dims[1] == 0) {
        char* s = expr_to_string_fullform(m);
        fprintf(stderr,
                "LatticeReduce::matrix: Argument %s at position 1 is not a "
                "non-empty rectangular matrix.\n", s);
        free(s);
        return NULL;
    }

    int n = (int)dims[0];
    int d = (int)dims[1];

    Expr** flat = malloc(sizeof(Expr*) * (size_t)n * (size_t)d);
    {
        size_t idx = 0;
        flatten_tensor(m, flat, &idx);
    }

    GRat* B = malloc(sizeof(GRat) * (size_t)n * (size_t)d);
    for (int i = 0; i < n * d; i++) gr_init(&B[i]);

    bool ok = true;
    for (int i = 0; i < n * d; i++) {
        if (!entry_to_grat(flat[i], &B[i])) { ok = false; break; }
    }
    for (int i = 0; i < n * d; i++) expr_free(flat[i]);
    free(flat);

    if (!ok) {
        fprintf(stderr,
                "LatticeReduce::latm: Matrix contains an entry that is not "
                "rational.\n");
        for (int i = 0; i < n * d; i++) gr_clear(&B[i]);
        free(B);
        return NULL;
    }

    int dependent = lll_reduce(B, n, d);
    if (dependent) {
        fprintf(stderr,
                "LatticeReduce::dep: The rows of the argument are linearly "
                "dependent; a basis (linearly independent generating set) is "
                "expected.\n");
        for (int i = 0; i < n * d; i++) gr_clear(&B[i]);
        free(B);
        return NULL;
    }

    Expr** rows = malloc(sizeof(Expr*) * (size_t)n);
    for (int r = 0; r < n; r++) {
        Expr** cells = malloc(sizeof(Expr*) * (size_t)d);
        for (int c = 0; c < d; c++) cells[c] = grat_to_expr(&ROW(B, r, d)[c]);
        rows[r] = expr_new_function(expr_new_symbol("List"), cells, (size_t)d);
        free(cells);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), rows, (size_t)n);
    free(rows);

    for (int i = 0; i < n * d; i++) gr_clear(&B[i]);
    free(B);
    return out;
}
