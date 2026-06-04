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
#include "eval.h"
#include "numeric.h"
#include "zero_test.h"

#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
 * 1 if the rows are linearly dependent.
 *
 * If `min_gso2` is non-NULL (and must be a pre-initialised mpq_t), it is
 * set to min_i |b*_i|^2 of the final reduced basis -- the rigorous lower
 * bound lambda_1(L)^2 >= min_i |b*_i|^2 on the shortest lattice vector,
 * which FindIntegerNullVector uses to certify that no short integer
 * relation exists.  It is left untouched when the rows are dependent. */
static int lll_reduce(GRat* B, int n, int d, mpq_t* min_gso2) {
    if (n <= 1) {                                     /* already reduced */
        if (min_gso2 && n == 1) {
            mpq_t t; mpq_init(t);
            mpq_set_ui(*min_gso2, 0, 1);
            for (int c = 0; c < d; c++) { gr_norm2(t, &ROW(B, 0, d)[c]); mpq_add(*min_gso2, *min_gso2, t); }
            mpq_clear(t);
        }
        return 0;
    }

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

    if (min_gso2 && !dependent) {
        mpq_set(*min_gso2, Bnorm[0]);
        for (int i = 1; i < n; i++)
            if (mpq_cmp(Bnorm[i], *min_gso2) < 0) mpq_set(*min_gso2, Bnorm[i]);
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

    int dependent = lll_reduce(B, n, d, NULL);
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

/* ================================================================== *
 *  FindIntegerNullVector  (PSLQ / integer-relation detection)         *
 *                                                                     *
 *  Given x = {x_1, ..., x_n}, find integers a = {a_1, ..., a_n}, not  *
 *  all zero, with a_1 x_1 + ... + a_n x_n == 0.  For complex x_i the  *
 *  a_i are Gaussian integers.                                         *
 *                                                                     *
 *  Method: build the integer-relation lattice whose rows are          *
 *      r_i = ( e_i | round(2^b x_i) )         (Gaussian round if      *
 *  complex), LLL-reduce it exactly (reusing the machinery above), and *
 *  read the relation off the shortest reduced row.  The certified     *
 *  lower bound on the norm of *any* integer relation is               *
 *      B = sqrt( M2 / (1 + (sqrt(n)/2)^2) ),   M2 = min_i |b*_i|^2,    *
 *  the rigorous LLL Gram-Schmidt bound lambda_1(L)^2 >= M2 combined    *
 *  with the worst-case |a . round(2^b x)| <= (sqrt(n)/2)||a|| coming   *
 *  from the rounding of the scaled reals.  B drives the no-relation    *
 *  diagnostics (norel / lgrelb / rnfb / rnfu).                        *
 * ================================================================== */

/* True if `e` contains any inexact (Real / MPFR) leaf. */
static bool finv_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type == EXPR_FUNCTION) {
        if (finv_is_inexact(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (finv_is_inexact(e->data.function.args[i])) return true;
    }
    return false;
}

#ifdef USE_MPFR
/* Largest MPFR precision (bits) carried by any leaf of `e`; 0 if none. */
static long finv_max_prec_bits(const Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_MPFR) return (long)mpfr_get_prec(e->data.mpfr);
    long m = 0;
    if (e->type == EXPR_FUNCTION) {
        long h = finv_max_prec_bits(e->data.function.head);
        if (h > m) m = h;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            long c = finv_max_prec_bits(e->data.function.args[i]);
            if (c > m) m = c;
        }
    }
    return m;
}
#endif

/* Real double value of a numeric Expr (Integer / BigInt / Rational /
 * Real / MPFR).  Returns false for anything else. */
static bool finv_num_to_double(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real;            return true; }
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint); return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    { *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true; }
#endif
    {
        mpq_t q; mpq_init(q);
        if (expr_to_mpq((Expr*)e, q)) { *out = mpq_get_d(q); mpq_clear(q); return true; }
        mpq_clear(q);
    }
    return false;
}

#ifndef USE_MPFR
/* Real + imaginary double parts of a numericalised entry (machine path). */
static bool finv_extract_double_pair(const Expr* e, double* re, double* im) {
    if (finv_num_to_double(e, re)) { *im = 0.0; return true; }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Complex
        && e->data.function.arg_count == 2) {
        if (finv_num_to_double(e->data.function.args[0], re)
         && finv_num_to_double(e->data.function.args[1], im))
            return true;
    }
    return false;
}
#endif

/* Set (Xre, Xim) = round(2^work_bits * value(entry)).  Returns false if
 * `entry` does not numericalise to a number.  Sets *any_complex when the
 * imaginary part is nonzero. */
static bool finv_scaled_value(const Expr* entry, long work_bits,
                              mpz_t Xre, mpz_t Xim, bool* any_complex) {
#ifdef USE_MPFR
    long guard = 64;
    NumericSpec spec; spec.mode = NUMERIC_MODE_MPFR; spec.bits = work_bits + guard;
    Expr* num = numericalize(entry, spec);
    mpfr_t re, im;
    mpfr_init2(re, work_bits + guard);
    mpfr_init2(im, work_bits + guard);
    bool inexact;
    bool ok = get_approx_mpfr(num, re, im, &inexact);
    if (ok) {
        mpfr_mul_2si(re, re, work_bits, MPFR_RNDN);
        mpfr_mul_2si(im, im, work_bits, MPFR_RNDN);
        mpfr_get_z(Xre, re, MPFR_RNDN);
        mpfr_get_z(Xim, im, MPFR_RNDN);
        if (mpz_sgn(Xim) != 0) *any_complex = true;
    }
    mpfr_clear(re); mpfr_clear(im);
    expr_free(num);
    return ok;
#else
    NumericSpec spec = numeric_machine_spec();
    Expr* num = numericalize(entry, spec);
    double re, im;
    bool ok = finv_extract_double_pair(num, &re, &im);
    if (ok) {
        double sre = round(ldexp(re, (int)work_bits));
        double sim = round(ldexp(im, (int)work_bits));
        mpz_set_d(Xre, sre);
        mpz_set_d(Xim, sim);
        if (im != 0.0) *any_complex = true;
    }
    expr_free(num);
    return ok;
#endif
}

/* Read a real-valued system variable ($MachinePrecision, $MaxExtraPrecision)
 * as a double, falling back to `fallback` if unset / non-real. */
static double finv_sysreal(const char* name, double fallback) {
    Expr* v = eval_and_free(expr_new_symbol(name));
    double out = fallback;
    if (v) { (void)finv_num_to_double(v, &out); expr_free(v); }
    return out;
}

/* Build the exact dot product Plus[Times[a_i, x_i], ...] for validation. */
static Expr* finv_build_dot(Expr** a_list, Expr** entries, size_t n) {
    Expr** terms = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* factors[2] = { expr_copy(a_list[i]), expr_copy(entries[i]) };
        terms[i] = expr_new_function(expr_new_symbol("Times"), factors, 2);
    }
    Expr* dot = expr_new_function(expr_new_symbol("Plus"), terms, n);
    free(terms);
    return dot;
}

/* Apply a user ZeroTest function f to the residual: returns 1 (zero),
 * 0 (nonzero), or -1 (undecided). */
static int finv_user_zerotest(Expr* fn, Expr* residual) {
    Expr** av = malloc(sizeof(Expr*));
    av[0] = expr_copy(residual);
    Expr* call = expr_new_function(expr_copy(fn), av, 1);
    free(av);
    Expr* r = eval_and_free(call);
    int verdict = -1;
    if (r && r->type == EXPR_SYMBOL) {
        if (strcmp(r->data.symbol, "True") == 0) verdict = 1;
        else if (strcmp(r->data.symbol, "False") == 0) verdict = 0;
    }
    expr_free(r);
    return verdict;
}

/* Emit a FindIntegerNullVector::tag diagnostic naming the input list. */
static void finv_msg(const char* tag, const char* body, Expr* list) {
    char* s = expr_to_string_fullform(list);
    fprintf(stderr, "FindIntegerNullVector::%s: %s %s.\n", tag, body, s);
    free(s);
}

Expr* builtin_findintegernullvector(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** av = res->data.function.args;

    /* --- peel trailing options (WorkingPrecision, ZeroTest) --- */
    size_t pos_end = argc;
    Expr* zerotest_fn = NULL;       /* borrowed; NULL = Automatic       */
    bool wp_explicit = false;
    long wp_bits = 0;               /* explicit working precision (bits)*/
    while (pos_end > 0) {
        Expr* o = av[pos_end - 1];
        if (o->type != EXPR_FUNCTION
            || o->data.function.head->type != EXPR_SYMBOL
            || (o->data.function.head->data.symbol != SYM_Rule
             && o->data.function.head->data.symbol != SYM_RuleDelayed)
            || o->data.function.arg_count != 2
            || o->data.function.args[0]->type != EXPR_SYMBOL)
            break;
        const char* opt = o->data.function.args[0]->data.symbol;
        Expr* val = o->data.function.args[1];
        if (opt == SYM_WorkingPrecision) {
            if (!(val->type == EXPR_SYMBOL && val->data.symbol == SYM_Automatic)) {
                double digits;
                if (!finv_num_to_double(val, &digits) || digits <= 0.0) return NULL;
                wp_explicit = true;
                wp_bits = numeric_digits_to_bits(digits);
                if (wp_bits < 8) wp_bits = 8;
            }
        } else if (strcmp(opt, "ZeroTest") == 0) {
            if (!(val->type == EXPR_SYMBOL && val->data.symbol == SYM_Automatic))
                zerotest_fn = val;
        } else {
            break;                  /* unknown trailing rule: not an option */
        }
        pos_end--;
    }

    /* --- positional arguments: list, optional norm bound d --- */
    if (pos_end < 1 || pos_end > 2) return NULL;
    Expr* list = av[0];
    if (list->type != EXPR_FUNCTION
        || list->data.function.head->type != EXPR_SYMBOL
        || list->data.function.head->data.symbol != SYM_List)
        return NULL;
    size_t n = list->data.function.arg_count;
    if (n < 2) return NULL;
    Expr** entries = list->data.function.args;

    bool has_bound = (pos_end == 2);
    double d_val = 0.0;
    if (has_bound && !finv_num_to_double(av[1], &d_val)) return NULL;

    bool inexact = false;
    for (size_t i = 0; i < n; i++)
        if (finv_is_inexact(entries[i])) { inexact = true; break; }

    /* --- precision schedule --- */
    double mach_digits = finv_sysreal("$MachinePrecision",
                                      NUMERIC_MACHINE_PRECISION_DIGITS);
    long start_bits, cap_bits;
    bool escalate;
    if (wp_explicit) {
        start_bits = cap_bits = wp_bits;
        escalate = false;
    } else if (inexact) {
        long ib = 0;
#ifdef USE_MPFR
        for (size_t i = 0; i < n; i++) {
            long b = finv_max_prec_bits(entries[i]);
            if (b > ib) ib = b;
        }
#endif
        if (ib < 8) ib = numeric_digits_to_bits(mach_digits);
        start_bits = cap_bits = ib;
        escalate = false;
    } else {
        double extra = finv_sysreal("$MaxExtraPrecision", 50.0);
        start_bits = numeric_digits_to_bits(mach_digits);
        cap_bits   = numeric_digits_to_bits(mach_digits + (extra > 0 ? extra : 0));
        if (cap_bits < start_bits) cap_bits = start_bits;
        escalate = true;
    }

    /* --- search loop --- */
    mpz_t* Xre = malloc(sizeof(mpz_t) * n);
    mpz_t* Xim = malloc(sizeof(mpz_t) * n);
    for (size_t i = 0; i < n; i++) { mpz_init(Xre[i]); mpz_init(Xim[i]); }

    Expr* result = NULL;            /* the relation vector, if accepted */
    double a_norm = 0.0, B_bound = 0.0;
    bool accepted = false, ztest1 = false, bail = false;
    long bits = start_bits;

    for (;;) {
        bool any_complex = false;
        bool numeric_ok = true;
        for (size_t i = 0; i < n; i++) {
            if (!finv_scaled_value(entries[i], bits, Xre[i], Xim[i], &any_complex)) {
                numeric_ok = false; break;
            }
        }
        if (!numeric_ok) { bail = true; break; }

        int nc = (int)n, dcol = (int)n + 1;
        GRat* B = malloc(sizeof(GRat) * (size_t)nc * (size_t)dcol);
        for (int i = 0; i < nc * dcol; i++) gr_init(&B[i]);
        for (int i = 0; i < nc; i++) {
            for (int c = 0; c < nc; c++)
                if (c == i) mpq_set_ui(ROW(B, i, dcol)[c].re, 1, 1);
            mpq_set_z(ROW(B, i, dcol)[nc].re, Xre[i]);
            mpq_set_z(ROW(B, i, dcol)[nc].im, Xim[i]);
        }

        mpq_t min_gso2; mpq_init(min_gso2);
        int dep = lll_reduce(B, nc, dcol, &min_gso2);

        /* pick the shortest reduced row */
        int best = 0;
        mpq_t bestnorm, rownorm, t;
        mpq_init(bestnorm); mpq_init(rownorm); mpq_init(t);
        for (int r = 0; r < nc; r++) {
            mpq_set_ui(rownorm, 0, 1);
            for (int c = 0; c < dcol; c++) { gr_norm2(t, &ROW(B, r, dcol)[c]); mpq_add(rownorm, rownorm, t); }
            if (r == 0 || mpq_cmp(rownorm, bestnorm) < 0) { mpq_set(bestnorm, rownorm); best = r; }
        }

        /* a-part norm^2 (integer) and relation vector */
        mpq_t anorm2; mpq_init(anorm2);
        for (int c = 0; c < nc; c++) { gr_norm2(t, &ROW(B, best, dcol)[c]); mpq_add(anorm2, anorm2, t); }
        Expr** a_cells = malloc(sizeof(Expr*) * n);
        for (int c = 0; c < nc; c++) a_cells[c] = grat_to_expr(&ROW(B, best, dcol)[c]);
        Expr* A = expr_new_function(expr_new_symbol("List"), a_cells, n);

        /* certified lower bound B = sqrt(M2/(1+n/4)) */
        double M2 = dep ? 0.0 : mpq_get_d(min_gso2);
        double kE = 0.5 * sqrt((double)n);
        B_bound = (M2 > 0.0 && isfinite(M2)) ? sqrt(M2 / (1.0 + kE * kE)) : INFINITY;
        a_norm = sqrt(mpq_get_d(anorm2));

        /* validate */
        bool ok;
        if (inexact) {
            /* Inexact input: the relation holds only to the precision of the
             * input, so the shortest vector is returned as-is (Wolfram does
             * the same -- it may be a spurious relation at low precision). */
            ok = true;
        } else {
            /* Confidence guard.  A spurious shortest vector of a no-relation
             * lattice has norm ~ 2^(bits/n); only a genuine relation is
             * anomalously short.  Require n*log2(|a|^2) < bits so that the
             * precision comfortably exceeds the relation's information
             * content -- this rejects large-coefficient artifacts (which can
             * fool an exact zero test by holding to tens of digits) and forces
             * escalation, matching Wolfram's conservative behaviour.
             *
             * A spurious shortest vector sits near |a| ~ 2^(bits/n), i.e.
             * n*log2(|a|^2) ~ 2*bits; a genuine relation is far below that.
             * The 1.35 threshold separates the two with margin on both sides
             * (verified against {E,Pi} -> rnfu and the degree-30 relation
             * found only at WorkingPrecision -> 300). */
            double confidence = (double)n * log2(mpq_get_d(anorm2));
            bool confident = confidence < 1.35 * (double)bits;
            if (!confident) {
                ok = false;
            } else {
                Expr* dot = eval_and_free(finv_build_dot(a_cells, entries, n));
                int verdict;
                if (zerotest_fn) verdict = finv_user_zerotest(zerotest_fn, dot);
                else {
                    ZeroTestResult zt = zero_test_decide(dot);
                    verdict = (zt == ZERO_TEST_TRUE) ? 1 : (zt == ZERO_TEST_FALSE) ? 0 : -1;
                }
                expr_free(dot);
                if (verdict == 0)      { ok = false; }
                else if (verdict == 1) { ok = true;  ztest1 = false; }
                else                   { ok = true;  ztest1 = true;  } /* assume zero */
            }
        }
        /* a_cells' element pointers were transferred into A by
         * expr_new_function (which copies only the array); free the array. */
        free(a_cells);

        if (ok) { result = A; accepted = true; }
        else { expr_free(A); ztest1 = false; }

        mpq_clear(bestnorm); mpq_clear(rownorm); mpq_clear(t);
        mpq_clear(anorm2); mpq_clear(min_gso2);
        for (int i = 0; i < nc * dcol; i++) gr_clear(&B[i]);
        free(B);

        if (accepted) break;
        if (!escalate || bits >= cap_bits) break;
        bits = (long)((double)bits * 1.6) + 16;
        if (bits >= cap_bits) bits = cap_bits;     /* one final try at the cap */
    }

    for (size_t i = 0; i < n; i++) { mpz_clear(Xre[i]); mpz_clear(Xim[i]); }
    free(Xre); free(Xim);

    if (bail) { if (result) expr_free(result); return NULL; }

    /* --- decision tree --- */
    double d2 = d_val * d_val, a2 = a_norm * a_norm;
    if (accepted) {
        if (!has_bound || a2 <= d2 * (1.0 + 1e-12)) {
            if (ztest1)
                fprintf(stderr, "FindIntegerNullVector::ztest1: Unable to decide "
                        "whether the residual is zero; assuming it is.\n");
            return result;                                  /* success */
        }
        /* a relation was found but it exceeds the norm bound d */
        if (d_val < B_bound) {
            finv_msg("norel", "There is no integer null vector for", list);
            expr_free(result);
            return NULL;
        }
        if (inexact) {
            char buf[160];
            snprintf(buf, sizeof buf,
                     "has not found an integer null vector (proven none with norm "
                     "<= %.6g) for", B_bound);
            finv_msg("rnfb", buf, list);
            expr_free(result);
            return NULL;
        }
        {
            char buf[200];
            snprintf(buf, sizeof buf,
                     "has proven no integer null vector with norm <= %.6g exists; "
                     "returning a null vector with larger norm for", B_bound);
            finv_msg("lgrelb", buf, list);
        }
        return result;                                       /* larger vector */
    }

    /* no relation found (exact input, precision cap reached) */
    if (!has_bound) {
        finv_msg("rnfu", "has not found an integer null vector for", list);
        return NULL;
    }
    if (d_val <= B_bound) {
        finv_msg("norel", "There is no integer null vector for", list);
    } else {
        char buf[160];
        snprintf(buf, sizeof buf,
                 "has not found an integer null vector (proven none with norm "
                 "<= %.6g) for", B_bound);
        finv_msg("rnfb", buf, list);
    }
    return NULL;
}
