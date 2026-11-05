/* rref_numeric.c -- tolerance-aware Gauss-Jordan on inexact matrices.
 * See rref_numeric.h for the contract and the zero tolerance. */

#include "rref_numeric.h"
#include "pack.h"
#include "la_dense.h"
#include "common.h"
#include "expr.h"
#include "sym_names.h"
#include "eval.h"
#include "ndarray.h"
#include "print.h"     /* expr_to_string: the ::luc text */
#include "message.h"   /* mth_message */

#include <stdlib.h>
#include <string.h>

#ifdef USE_MPFR
#include <mpfr.h>
#include "numeric.h"          /* get_approx_mpfr */
#include "numeric_complex.h"  /* numeric_mpfr_make_complex */

/* Entry state: what the reference prints for an entry the elimination has
 * not turned into a computed number. */
enum { NM_NUM = 0,      /* a computed (inexact) number                    */
       NM_ZERO = 1,     /* exact 0: an exact input zero, below a pivot, a */
                        /* zero row                                       */
       NM_ONE = 2,      /* exact 1: a pivot                               */
       NM_ELIM = 3 };   /* zeroed above a pivot: inexact 0 for a machine  */
                        /* matrix, exact 0 for an arbitrary-precision one */

typedef struct {
    int r, c;
    bool cplx;
    mpfr_prec_t bits;
    mpfr_t* re;           /* r * c, row-major */
    mpfr_t* im;           /* NULL when real   */
    unsigned char* st;    /* NM_* per entry   */
} NMat;

#define NM_AT(A, i, j) ((size_t)(i) * (size_t)(A)->c + (size_t)(j))

static void nm_free(NMat* A)
{
    size_t n = (size_t)A->r * (size_t)A->c;
    for (size_t k = 0; k < n; k++) {
        mpfr_clear(A->re[k]);
        if (A->im) mpfr_clear(A->im[k]);
    }
    free(A->re); free(A->im); free(A->st);
    A->re = A->im = NULL; A->st = NULL;
}

static bool nm_alloc(NMat* A, int r, int c, bool cplx, mpfr_prec_t bits)
{
    size_t n = (size_t)r * (size_t)c;
    A->r = r; A->c = c; A->cplx = cplx; A->bits = bits;
    A->re = (mpfr_t*)malloc(sizeof(mpfr_t) * (n ? n : 1));
    A->im = cplx ? (mpfr_t*)malloc(sizeof(mpfr_t) * (n ? n : 1)) : NULL;
    A->st = (unsigned char*)calloc(n ? n : 1, 1);
    if (!A->re || (cplx && !A->im) || !A->st) {
        free(A->re); free(A->im); free(A->st);
        return false;
    }
    for (size_t k = 0; k < n; k++) {
        mpfr_init2(A->re[k], bits);
        mpfr_set_zero(A->re[k], 1);
        if (cplx) { mpfr_init2(A->im[k], bits); mpfr_set_zero(A->im[k], 1); }
        A->st[k] = NM_ZERO;
    }
    return true;
}

/* Exact integer zero? (an input 0 the reference keeps exact) */
static bool nm_is_exact_zero(const Expr* e)
{
    return e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* Rows of a boxed matrix (or a vector, as one column) as Expr rows. */
static bool nm_shape(const Expr* m, int* r, int* c, bool* is_vec)
{
    *is_vec = false;
    if (!la_looks_like_matrix(m)) {
        if (!head_is(m, SYM_List) || m->data.function.arg_count == 0) return false;
        for (size_t i = 0; i < m->data.function.arg_count; i++)
            if (head_is(m->data.function.args[i], SYM_List)) return false;
        *r = (int)m->data.function.arg_count; *c = 1; *is_vec = true;
        return true;
    }
    return la_matrix_dims(m, r, c) && *r > 0 && *c > 0;
}

/* Load the boxed matrix (or vector) `m` into columns [c0, c0 + width) of A. */
static bool nm_load_block(NMat* A, const Expr* m, int c0)
{
    int r, c; bool is_vec;
    if (!nm_shape(m, &r, &c, &is_vec) || r != A->r) return false;
    mpfr_t im; mpfr_init2(im, A->bits);
    bool ok = true;
    for (int i = 0; i < r && ok; i++) {
        for (int j = 0; j < c && ok; j++) {
            const Expr* e = is_vec ? m->data.function.args[i] : la_at(m, i, j);
            size_t k = NM_AT(A, i, c0 + j);
            bool inexact = false;
            mpfr_set_zero(im, 1);
            if (!get_approx_mpfr(e, A->re[k], im, &inexact)) { ok = false; break; }
            if (A->cplx) mpfr_set(A->im[k], im, MPFR_RNDN);
            else if (!mpfr_zero_p(im)) { ok = false; break; }
            A->st[k] = nm_is_exact_zero(e) ? NM_ZERO : NM_NUM;
        }
    }
    mpfr_clear(im);
    return ok;
}

/* |A[k]| into out. */
static void nm_abs(const NMat* A, size_t k, mpfr_t out)
{
    if (A->cplx) mpfr_hypot(out, A->re[k], A->im[k], MPFR_RNDN);
    else         mpfr_abs(out, A->re[k], MPFR_RNDN);
}

/* maxabs = max |A[:, 0:pc]|; tol = f max(r, c) eps maxabs with f = 1 for
 * a machine matrix -- the reference's own rule, probed on diagonal matrices of
 * order 2, 3 and 6 (diag(1, 1, e) has rank 3 exactly when e > 3 eps) -- and
 * f = 64 at arbitrary precision, where the reference tracks the error
 * instead and a fixed few-ulp margin would call rounding residue a pivot. */
static void nm_tolerance(const NMat* A, int pc, mpfr_t tol, mpfr_t maxabs)
{
    mpfr_t t; mpfr_init2(t, A->bits);
    mpfr_set_zero(maxabs, 1);
    for (int i = 0; i < A->r; i++)
        for (int j = 0; j < pc; j++) {
            size_t k = NM_AT(A, i, j);
            if (A->st[k] == NM_ZERO) continue;
            nm_abs(A, k, t);
            if (mpfr_cmp(t, maxabs) > 0) mpfr_set(maxabs, t, MPFR_RNDN);
        }
    int dim = pc > A->r ? pc : A->r;
    mpfr_mul_ui(tol, maxabs, (unsigned long)((A->bits <= 53 ? 1 : 64) * dim), MPFR_RNDN);
    mpfr_mul_2si(tol, tol, 1 - (long)A->bits, MPFR_RNDN);   /* eps = 2^(1 - bits) */
    mpfr_clear(t);
}

static void nm_swap_rows(NMat* A, int a, int b)
{
    if (a == b) return;
    for (int j = 0; j < A->c; j++) {
        size_t x = NM_AT(A, a, j), y = NM_AT(A, b, j);
        mpfr_swap(A->re[x], A->re[y]);
        if (A->cplx) mpfr_swap(A->im[x], A->im[y]);
        unsigned char s = A->st[x]; A->st[x] = A->st[y]; A->st[y] = s;
    }
}

/* Gauss-Jordan to reduced row echelon form, pivots searched in the first
 * `pc` columns only (the coefficient block of an augmented system).  Fills
 * pivcol[0 .. rank) and returns the rank.  *ill (when non-NULL) is set when a
 * pivot was used whose magnitude is below 10^6 eps of the largest entry: the
 * reference's RowReduce::luc (probed on diag(1, e): e = 1e-10 warns, 3e-10
 * does not, at every scale). */
static int nm_rref(NMat* A, int pc, int* pivcol, bool* ill)
{
    mpfr_prec_t bits = A->bits;
    mpfr_t tol, best, t, ir, ii, fr, fi, pr, pi, tr, ti, maxabs, illtol;
    mpfr_inits2(bits, tol, best, t, ir, ii, fr, fi, pr, pi, tr, ti, maxabs, illtol,
                (mpfr_ptr)0);
    nm_tolerance(A, pc, tol, maxabs);
    mpfr_mul_d(illtol, maxabs, 1.0e6 * 2.220446049250313e-16, MPFR_RNDN);
    if (ill) *ill = false;
    int rank = 0;
    for (int col = 0; col < pc && rank < A->r; col++) {
        int p = -1;
        mpfr_set_zero(best, 1);
        for (int i = rank; i < A->r; i++) {
            size_t k = NM_AT(A, i, col);
            if (A->st[k] == NM_ZERO) continue;
            nm_abs(A, k, t);
            if (mpfr_cmp(t, best) > 0) { mpfr_set(best, t, MPFR_RNDN); p = i; }
        }
        if (p < 0 || mpfr_cmp(best, tol) <= 0) {
            /* a free column: what is left below the pivots is rounding */
            for (int i = rank; i < A->r; i++) {
                size_t k = NM_AT(A, i, col);
                mpfr_set_zero(A->re[k], 1);
                if (A->cplx) mpfr_set_zero(A->im[k], 1);
                A->st[k] = NM_ZERO;
            }
            continue;
        }
        if (ill && mpfr_cmp(best, illtol) < 0) *ill = true;
        nm_swap_rows(A, p, rank);
        /* inverse of the pivot */
        size_t kp = NM_AT(A, rank, col);
        if (A->cplx) {
            /* 1/(a + bi) = (a - bi)/(a^2 + b^2) */
            mpfr_sqr(t, A->re[kp], MPFR_RNDN);
            mpfr_sqr(tr, A->im[kp], MPFR_RNDN);
            mpfr_add(t, t, tr, MPFR_RNDN);
            mpfr_div(ir, A->re[kp], t, MPFR_RNDN);
            mpfr_div(ii, A->im[kp], t, MPFR_RNDN);
            mpfr_neg(ii, ii, MPFR_RNDN);
        } else {
            mpfr_ui_div(ir, 1, A->re[kp], MPFR_RNDN);
        }
        /* scale the pivot row */
        for (int j = col + 1; j < A->c; j++) {
            size_t k = NM_AT(A, rank, j);
            if (A->st[k] == NM_ZERO) continue;
            if (A->cplx) {
                mpfr_mul(tr, A->re[k], ir, MPFR_RNDN);
                mpfr_mul(t, A->im[k], ii, MPFR_RNDN);
                mpfr_sub(tr, tr, t, MPFR_RNDN);
                mpfr_mul(ti, A->re[k], ii, MPFR_RNDN);
                mpfr_mul(t, A->im[k], ir, MPFR_RNDN);
                mpfr_add(A->im[k], ti, t, MPFR_RNDN);
                mpfr_set(A->re[k], tr, MPFR_RNDN);
            } else {
                mpfr_mul(A->re[k], A->re[k], ir, MPFR_RNDN);
            }
            A->st[k] = NM_NUM;
        }
        mpfr_set_ui(A->re[kp], 1, MPFR_RNDN);
        if (A->cplx) mpfr_set_zero(A->im[kp], 1);
        A->st[kp] = NM_ONE;
        /* eliminate the column from every other row */
        for (int i = 0; i < A->r; i++) {
            if (i == rank) continue;
            size_t kf = NM_AT(A, i, col);
            if (A->st[kf] == NM_ZERO) continue;
            mpfr_set(fr, A->re[kf], MPFR_RNDN);
            if (A->cplx) mpfr_set(fi, A->im[kf], MPFR_RNDN);
            for (int j = col + 1; j < A->c; j++) {
                size_t ks = NM_AT(A, rank, j);
                if (A->st[ks] == NM_ZERO) continue;
                size_t k = NM_AT(A, i, j);
                if (A->cplx) {
                    /* a -= f * s */
                    mpfr_mul(pr, fr, A->re[ks], MPFR_RNDN);
                    mpfr_mul(t, fi, A->im[ks], MPFR_RNDN);
                    mpfr_sub(pr, pr, t, MPFR_RNDN);
                    mpfr_mul(pi, fr, A->im[ks], MPFR_RNDN);
                    mpfr_mul(t, fi, A->re[ks], MPFR_RNDN);
                    mpfr_add(pi, pi, t, MPFR_RNDN);
                    mpfr_sub(A->re[k], A->re[k], pr, MPFR_RNDN);
                    mpfr_sub(A->im[k], A->im[k], pi, MPFR_RNDN);
                } else {
                    mpfr_mul(pr, fr, A->re[ks], MPFR_RNDN);
                    mpfr_sub(A->re[k], A->re[k], pr, MPFR_RNDN);
                }
                A->st[k] = NM_NUM;
            }
            mpfr_set_zero(A->re[kf], 1);
            if (A->cplx) mpfr_set_zero(A->im[kf], 1);
            A->st[kf] = i < rank ? NM_ELIM : NM_ZERO;
        }
        pivcol[rank++] = col;
    }
    mpfr_clears(tol, best, t, ir, ii, fr, fi, pr, pi, tr, ti, maxabs, illtol, (mpfr_ptr)0);
    return rank;
}

/* One output entry. A MACHINE matrix answers machine numbers throughout --
 * the project rule (a machine matrix in, a machine matrix out, one element
 * head: check-array-exactness), a documented deviation from the reference,
 * which writes exact 1 / 0 for the structural entries
 * (RowReduce[{{2., 0.}, {0., 4.}}] is {{1, 0.}, {0, 1}} there). An
 * arbitrary-precision matrix keeps the reference's exact pivots and zeros. */
static Expr* nm_entry(const NMat* A, size_t k)
{
    bool machine = A->bits <= 53;
    switch (A->st[k]) {
        case NM_ONE:  return machine ? expr_new_real(1.0) : expr_new_integer(1);
        case NM_ZERO: return machine ? expr_new_real(0.0) : expr_new_integer(0);
        case NM_ELIM: return machine ? expr_new_real(0.0) : expr_new_integer(0);
        default: break;
    }
    if (machine) {
        double re = mpfr_get_d(A->re[k], MPFR_RNDN);
        double im = A->cplx ? mpfr_get_d(A->im[k], MPFR_RNDN) : 0.0;
        if (im == 0.0) return expr_new_real(re);
        return expr_new_function(expr_new_symbol(SYM_Complex),
                                 (Expr*[]){ expr_new_real(re), expr_new_real(im) }, 2);
    }
    if (A->cplx && !mpfr_zero_p(A->im[k]))
        return numeric_mpfr_make_complex(A->re[k], A->im[k]);
    return expr_new_mpfr_copy(A->re[k]);
}

/* An inexact zero at the working precision (a free variable). */
static Expr* nm_zero(mpfr_prec_t bits)
{
    if (bits <= 53) return expr_new_real(0.0);
    return expr_new_mpfr_from_si(0, bits);
}

/* Load a numeric boxed copy of `m` (an NDArray is unpacked). */
static Expr* nm_boxed(const Expr* m)
{
    return is_ndarray(m) ? la_boxed_copy(m) : expr_copy((Expr*)m);
}

long rref_numeric_bits(const Expr* m, const Expr* b)
{
    if (!m || !la_all_numeric(m) || (b && !la_all_numeric(b))) return 0;
    CommonInexactInfo im = common_scan_inexact(m);
    CommonInexactInfo ib = b ? common_scan_inexact(b) : (CommonInexactInfo){ false, 0 };
    if (!im.has_inexact && !ib.has_inexact) return 0;
    long bits = 0;
    if (im.has_inexact) bits = im.min_bits;
    if (ib.has_inexact && (bits == 0 || ib.min_bits < bits)) bits = ib.min_bits;
    return bits <= 53 ? 53 : bits;
}

Expr* rref_numeric_rowreduce(const Expr* m0, long bits)
{
    Expr* m = nm_boxed(m0);
    int r, c; bool is_vec;
    if (!nm_shape(m, &r, &c, &is_vec) || is_vec) { expr_free(m); return NULL; }
    NMat A;
    if (!nm_alloc(&A, r, c, la_has_complex(m), (mpfr_prec_t)bits)) { expr_free(m); return NULL; }
    if (!nm_load_block(&A, m, 0)) { nm_free(&A); expr_free(m); return NULL; }
    expr_free(m);
    int* piv = (int*)malloc(sizeof(int) * (size_t)(r < c ? r : c) + sizeof(int));
    bool ill = false;
    int rank = nm_rref(&A, c, piv, &ill);
    if (ill && bits <= 53) {
        char* s = expr_to_string((Expr*)m0);
        mth_message("RowReduce", "luc", "Result for RowReduce of badly conditioned "
                    "matrix %s may contain significant numerical errors.", s ? s : "");
        free(s);
    }
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)r);
    for (int i = 0; i < r; i++) {
        Expr** el = (Expr**)malloc(sizeof(Expr*) * (size_t)c);
        for (int j = 0; j < c; j++)
            el[j] = i >= rank ? (A.bits <= 53 ? expr_new_real(0.0) : expr_new_integer(0))
                              : nm_entry(&A, NM_AT(&A, i, j));
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), el, (size_t)c);
        free(el);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)r);
    /* a machine answer has one machine head throughout, so it packs (as the
     * LAPACK path's did) */
    if (A.bits <= 53) out = pack_offer(out);
    free(rows); free(piv);
    nm_free(&A);
    return out;
}

Expr* rref_numeric_linearsolve(const Expr* m0, const Expr* b0, long bits, bool* nosol)
{
    *nosol = false;
    Expr* m = nm_boxed(m0);
    Expr* b = nm_boxed(b0);
    int r, c, br, bc; bool mv, bv;
    if (!nm_shape(m, &r, &c, &mv) || mv || !nm_shape(b, &br, &bc, &bv) || br != r) {
        expr_free(m); expr_free(b); return NULL;
    }
    NMat A;
    bool cplx = la_has_complex(m) || la_has_complex(b);
    if (!nm_alloc(&A, r, c + bc, cplx, (mpfr_prec_t)bits)) { expr_free(m); expr_free(b); return NULL; }
    if (!nm_load_block(&A, m, 0) || !nm_load_block(&A, b, c)) {
        nm_free(&A); expr_free(m); expr_free(b); return NULL;
    }
    /* the right-hand side's own scale, for the consistency test */
    mpfr_t bmax, amax, t, tol;
    mpfr_inits2((mpfr_prec_t)bits, bmax, amax, t, tol, (mpfr_ptr)0);
    mpfr_set_zero(bmax, 1); mpfr_set_zero(amax, 1);
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c + bc; j++) {
            nm_abs(&A, NM_AT(&A, i, j), t);
            mpfr_ptr dst = j < c ? amax : bmax;
            if (mpfr_cmp(t, dst) > 0) mpfr_set(dst, t, MPFR_RNDN);
        }
    }
    expr_free(m); expr_free(b);
    int* piv = (int*)malloc(sizeof(int) * (size_t)(r < c ? r : c) + sizeof(int));
    int rank = nm_rref(&A, c, piv, NULL);
    /* consistency: the right-hand side left in the zero rows must be rounding,
     * measured against max(|b|, |A| |x|) */
    mpfr_t xmax; mpfr_init2(xmax, (mpfr_prec_t)bits);
    mpfr_set_zero(xmax, 1);
    for (int i = 0; i < rank; i++)
        for (int j = c; j < c + bc; j++) {
            nm_abs(&A, NM_AT(&A, i, j), t);
            if (mpfr_cmp(t, xmax) > 0) mpfr_set(xmax, t, MPFR_RNDN);
        }
    mpfr_mul(tol, amax, xmax, MPFR_RNDN);
    if (mpfr_cmp(bmax, tol) > 0) mpfr_set(tol, bmax, MPFR_RNDN);
    int dim = r > c ? r : c;
    mpfr_mul_ui(tol, tol, (unsigned long)(64 * dim), MPFR_RNDN);
    mpfr_mul_2si(tol, tol, -(long)bits, MPFR_RNDN);
    bool consistent = true;
    for (int i = rank; i < r && consistent; i++)
        for (int j = c; j < c + bc; j++) {
            size_t k = NM_AT(&A, i, j);
            if (A.st[k] == NM_ZERO) continue;
            nm_abs(&A, k, t);
            if (mpfr_cmp(t, tol) > 0) { consistent = false; break; }
        }
    mpfr_clears(bmax, amax, t, tol, xmax, (mpfr_ptr)0);
    if (!consistent) {
        *nosol = true;
        free(piv); nm_free(&A);
        return NULL;
    }
    /* x[pivcol[i]] = R[i, c + j]; every free variable 0 */
    int* row_of = (int*)malloc(sizeof(int) * (size_t)c);
    for (int j = 0; j < c; j++) row_of[j] = -1;
    for (int i = 0; i < rank; i++) row_of[piv[i]] = i;
    Expr** xs = (Expr**)malloc(sizeof(Expr*) * (size_t)c);
    for (int v = 0; v < c; v++) {
        Expr** el = (Expr**)malloc(sizeof(Expr*) * (size_t)bc);
        for (int j = 0; j < bc; j++) {
            if (row_of[v] < 0) { el[j] = nm_zero((mpfr_prec_t)bits); continue; }
            size_t k = NM_AT(&A, row_of[v], c + j);
            el[j] = A.st[k] == NM_NUM ? nm_entry(&A, k) : nm_zero((mpfr_prec_t)bits);
        }
        if (bv) { xs[v] = el[0]; free(el); }
        else {
            xs[v] = expr_new_function(expr_new_symbol(SYM_List), el, (size_t)bc);
            free(el);
        }
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), xs, (size_t)c);
    free(xs); free(row_of); free(piv);
    nm_free(&A);
    return out;
}

Expr* rref_numeric_inverse(const Expr* m0, long bits, bool* singular)
{
    *singular = false;
    Expr* m = nm_boxed(m0);
    int r, c; bool is_vec;
    if (!nm_shape(m, &r, &c, &is_vec) || is_vec || r != c) { expr_free(m); return NULL; }
    NMat A;
    if (!nm_alloc(&A, r, 2 * c, la_has_complex(m), (mpfr_prec_t)bits)) { expr_free(m); return NULL; }
    if (!nm_load_block(&A, m, 0)) { nm_free(&A); expr_free(m); return NULL; }
    expr_free(m);
    for (int i = 0; i < r; i++) {
        size_t k = NM_AT(&A, i, c + i);
        mpfr_set_ui(A.re[k], 1, MPFR_RNDN);
        A.st[k] = NM_NUM;
    }
    int* piv = (int*)malloc(sizeof(int) * (size_t)r);
    int rank = nm_rref(&A, c, piv, NULL);
    free(piv);
    if (rank < r) { *singular = true; nm_free(&A); return NULL; }
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)r);
    for (int i = 0; i < r; i++) {
        Expr** el = (Expr**)malloc(sizeof(Expr*) * (size_t)c);
        for (int j = 0; j < c; j++) {
            size_t k = NM_AT(&A, i, c + j);
            /* an entry the elimination never touched stays exact 0, as in the
             * reference (Inverse[N[{{2, 0}, {0, 4}}, 30]] has exact zeros) */
            el[j] = A.st[k] == NM_NUM ? nm_entry(&A, k)
                  : (bits <= 53 ? expr_new_real(0.0) : expr_new_integer(0));
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), el, (size_t)c);
        free(el);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)r);
    /* a machine answer has one machine head throughout, so it packs (as the
     * LAPACK path's did) */
    if (A.bits <= 53) out = pack_offer(out);
    free(rows);
    nm_free(&A);
    return out;
}

Expr* rref_numeric_nullspace(const Expr* m0, long bits)
{
    Expr* m = nm_boxed(m0);
    int r, c; bool is_vec;
    if (!nm_shape(m, &r, &c, &is_vec) || is_vec) { expr_free(m); return NULL; }
    bool cplx = la_has_complex(m);
    NMat A;
    if (!nm_alloc(&A, r, c, cplx, (mpfr_prec_t)bits)) { expr_free(m); return NULL; }
    if (!nm_load_block(&A, m, 0)) { nm_free(&A); expr_free(m); return NULL; }
    expr_free(m);
    int* piv = (int*)malloc(sizeof(int) * (size_t)(r < c ? r : c) + sizeof(int));
    int rank = nm_rref(&A, c, piv, NULL);
    int nfree = c - rank;
    if (nfree == 0) { free(piv); nm_free(&A); return la_empty_list(); }
    bool* is_piv = (bool*)calloc((size_t)c, sizeof(bool));
    for (int i = 0; i < rank; i++) is_piv[piv[i]] = true;
    /* basis vector per free column f: v[f] = 1, v[piv[i]] = -R[i, f] */
    NMat V;
    nm_alloc(&V, nfree, c, cplx, (mpfr_prec_t)bits);
    int q = 0;
    for (int f = 0; f < c; f++) {
        if (is_piv[f]) continue;
        size_t kf = NM_AT(&V, q, f);
        mpfr_set_ui(V.re[kf], 1, MPFR_RNDN);
        V.st[kf] = NM_NUM;
        for (int i = 0; i < rank; i++) {
            size_t ks = NM_AT(&A, i, f), kd = NM_AT(&V, q, piv[i]);
            if (A.st[ks] != NM_NUM) continue;
            mpfr_neg(V.re[kd], A.re[ks], MPFR_RNDN);
            if (cplx) mpfr_neg(V.im[kd], A.im[ks], MPFR_RNDN);
            V.st[kd] = NM_NUM;
        }
        q++;
    }
    free(is_piv); free(piv);
    nm_free(&A);
    /* orthonormalise (modified Gram-Schmidt, Hermitian inner product) */
    mpfr_t sr, si, t, nrm;
    mpfr_inits2((mpfr_prec_t)bits, sr, si, t, nrm, (mpfr_ptr)0);
    for (int a = 0; a < nfree; a++) {
        for (int p = 0; p < a; p++) {
            /* s = <v_p, v_a> = sum conj(v_p) v_a */
            mpfr_set_zero(sr, 1); mpfr_set_zero(si, 1);
            for (int j = 0; j < c; j++) {
                size_t kp = NM_AT(&V, p, j), ka = NM_AT(&V, a, j);
                mpfr_mul(t, V.re[kp], V.re[ka], MPFR_RNDN); mpfr_add(sr, sr, t, MPFR_RNDN);
                if (cplx) {
                    mpfr_mul(t, V.im[kp], V.im[ka], MPFR_RNDN); mpfr_add(sr, sr, t, MPFR_RNDN);
                    mpfr_mul(t, V.re[kp], V.im[ka], MPFR_RNDN); mpfr_add(si, si, t, MPFR_RNDN);
                    mpfr_mul(t, V.im[kp], V.re[ka], MPFR_RNDN); mpfr_sub(si, si, t, MPFR_RNDN);
                }
            }
            /* v_a -= s v_p */
            for (int j = 0; j < c; j++) {
                size_t kp = NM_AT(&V, p, j), ka = NM_AT(&V, a, j);
                mpfr_mul(t, sr, V.re[kp], MPFR_RNDN); mpfr_sub(V.re[ka], V.re[ka], t, MPFR_RNDN);
                if (cplx) {
                    mpfr_mul(t, si, V.im[kp], MPFR_RNDN); mpfr_add(V.re[ka], V.re[ka], t, MPFR_RNDN);
                    mpfr_mul(t, sr, V.im[kp], MPFR_RNDN); mpfr_sub(V.im[ka], V.im[ka], t, MPFR_RNDN);
                    mpfr_mul(t, si, V.re[kp], MPFR_RNDN); mpfr_sub(V.im[ka], V.im[ka], t, MPFR_RNDN);
                }
                V.st[ka] = NM_NUM;
            }
        }
        mpfr_set_zero(nrm, 1);
        for (int j = 0; j < c; j++) {
            size_t k = NM_AT(&V, a, j);
            mpfr_sqr(t, V.re[k], MPFR_RNDN); mpfr_add(nrm, nrm, t, MPFR_RNDN);
            if (cplx) { mpfr_sqr(t, V.im[k], MPFR_RNDN); mpfr_add(nrm, nrm, t, MPFR_RNDN); }
        }
        mpfr_sqrt(nrm, nrm, MPFR_RNDN);
        for (int j = 0; j < c; j++) {
            size_t k = NM_AT(&V, a, j);
            if (V.st[k] != NM_NUM) continue;
            mpfr_div(V.re[k], V.re[k], nrm, MPFR_RNDN);
            if (cplx) mpfr_div(V.im[k], V.im[k], nrm, MPFR_RNDN);
        }
    }
    mpfr_clears(sr, si, t, nrm, (mpfr_ptr)0);
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)nfree);
    for (int a = 0; a < nfree; a++) {
        Expr** el = (Expr**)malloc(sizeof(Expr*) * (size_t)c);
        for (int j = 0; j < c; j++) el[j] = nm_entry(&V, NM_AT(&V, a, j));
        rows[a] = expr_new_function(expr_new_symbol(SYM_List), el, (size_t)c);
        free(el);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)nfree);
    free(rows);
    nm_free(&V);
    return out;
}

#else  /* !USE_MPFR */

long rref_numeric_bits(const Expr* m, const Expr* b) { (void)m; (void)b; return 0; }
Expr* rref_numeric_rowreduce(const Expr* m, long bits) { (void)m; (void)bits; return NULL; }
Expr* rref_numeric_nullspace(const Expr* m, long bits) { (void)m; (void)bits; return NULL; }
Expr* rref_numeric_linearsolve(const Expr* m, const Expr* b, long bits, bool* nosol)
{ (void)m; (void)b; (void)bits; *nosol = false; return NULL; }
Expr* rref_numeric_inverse(const Expr* m, long bits, bool* singular)
{ (void)m; (void)bits; *singular = false; return NULL; }

#endif /* USE_MPFR */
