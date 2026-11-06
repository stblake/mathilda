/* src/linalg/dot_mpfr.c -- Dot of arbitrary-precision (MPFR) vectors and
 * matrices as fused dot products.
 *
 * The generic Dot builds Times[a_ik, b_kj] and Plus[...] for every output
 * element and evaluates them: a 256 x 256 matrix times a vector at 800 digits
 * spent 0.35 s, ~5 us per multiply-add, nearly all of it in Expr allocation,
 * evaluation and per-operation MPFR temporaries.  Here every output element is
 * one fused dot product: with FLINT, arb_approx_dot over copies of the leaves'
 * mantissas (~1.6x faster than mpfr_dot at 500-800 digits: the products are
 * accumulated in fixed point and rounded once), else mpfr_dot over pointers
 * straight into the leaves (exact products, one correct rounding).  Either is
 * at least as accurate as the pairwise Times / Plus it replaces.  Large
 * products are split by rows over worker threads.
 *
 * Domain: both operands of rank 1 or 2 (rectangular Lists), at least one of
 * them entirely MPFR reals, the other MPFR reals or exact Integers /
 * Rationals.  So every product has an MPFR factor and every output element is
 * an MPFR number, exactly as the generic path would give.  An exact 0 is
 * admitted: Times[0, x] of an MPFR x is an MPFR zero at x's precision, so it
 * contributes nothing to the sum and imposes x's precision, which is what the
 * fused sum with that leaf read as an MPFR zero does (a DiagonalMatrix or any
 * other sparse exact operand: m . DiagonalMatrix[w] declined to n^3 symbolic
 * Times / Plus before, 1.4 s at 64 x 64).  Anything else -- a machine Real
 * (contagion to machine precision), a Complex, a symbol -- declines to the
 * generic path.
 *
 * Precision follows the scalar arithmetic: Times and Plus of MPFR numbers
 * answer at the smaller operand precision (an exact factor imposes none), so
 * output (r, s) gets min(min_k prec a_rk, min_k prec b_ks) over the MPFR
 * entries.  Cancellation in each sum is
 * reported to numeric_note_cancellation_bits (bits lost and bits kept), as a
 * cancelling Plus would, so it takes part in N[expr, p]'s certification.
 *
 * Threads: mp_parallel_for (mp_parallel.c); the workers touch nothing but
 * MPFR / arb numbers (no Expr, no evaluator).
 */
#include "dot_mpfr.h"
#include "linalg.h"
#include "expr.h"
#include "numeric.h"
#include "sym_names.h"
#include "mp_parallel.h"

#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>


#ifndef USE_MPFR

Expr* dot_mpfr(const Expr* a, const Expr* b) { (void)a; (void)b; return NULL; }

#else

#include <gmp.h>
#include <mpfr.h>
#ifdef USE_FLINT
#include <flint/flint.h>
#include <flint/arb.h>
#endif


/* Leaf classes. */
enum { DL_BAD = 0, DL_MPFR = 1, DL_EXACT = 2 };

static int dotm_leaf(const Expr* e)
{
    if (e->type == EXPR_MPFR) return DL_MPFR;
    if (e->type == EXPR_INTEGER) return DL_EXACT;
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint) != 0 ? DL_EXACT : DL_BAD;
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational)
        return DL_EXACT;                     /* a Rational is never 0 */
    return DL_BAD;
}

static bool dotm_is_list(const Expr* e)
{
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/* An operand as a rows x cols grid of leaves (a vector is 1 x n for the left
 * operand, n x 1 for the right).  Fills `leaf` row-major; returns false for
 * anything that is not a rectangular rank-1/2 List of admissible leaves.
 * *all_mpfr reports whether every leaf is MPFR. */
typedef struct {
    int rank;
    int64_t rows, cols;
    const Expr** leaf;
    bool all_mpfr;
} DotmOp;

static bool dotm_read(const Expr* e, bool left, DotmOp* op)
{
    op->leaf = NULL;
    if (!dotm_is_list(e) || e->data.function.arg_count == 0) return false;
    size_t n = e->data.function.arg_count;
    const Expr* r0 = e->data.function.args[0];
    op->all_mpfr = true;
    if (!dotm_is_list(r0)) {
        op->rank = 1;
        op->rows = left ? 1 : (int64_t)n;
        op->cols = left ? (int64_t)n : 1;
        op->leaf = (const Expr**)malloc(sizeof(Expr*) * n);
        if (!op->leaf) return false;
        for (size_t i = 0; i < n; i++) {
            int c = dotm_leaf(e->data.function.args[i]);
            if (c == DL_BAD) { free(op->leaf); op->leaf = NULL; return false; }
            if (c != DL_MPFR) op->all_mpfr = false;
            op->leaf[i] = e->data.function.args[i];
        }
        return true;
    }
    size_t m = r0->data.function.arg_count;
    if (m == 0) return false;
    op->rank = 2;
    op->rows = (int64_t)n;
    op->cols = (int64_t)m;
    op->leaf = (const Expr**)malloc(sizeof(Expr*) * n * m);
    if (!op->leaf) return false;
    for (size_t i = 0; i < n; i++) {
        const Expr* ri = e->data.function.args[i];
        if (!dotm_is_list(ri) || ri->data.function.arg_count != m) {
            free(op->leaf); op->leaf = NULL; return false;
        }
        for (size_t j = 0; j < m; j++) {
            const Expr* x = ri->data.function.args[j];
            int c = dotm_leaf(x);
            if (c == DL_BAD) { free(op->leaf); op->leaf = NULL; return false; }
            if (c != DL_MPFR) op->all_mpfr = false;
            op->leaf[i * m + j] = x;
        }
    }
    return true;
}

/* The work: C[r, s] = sum_k A[r, k] B[k, s]. */
typedef struct {
    mpfr_ptr* A;        /* R x K row-major pointers */
    mpfr_ptr* Bt;       /* S x K: column s of B, contiguous */
#ifdef USE_FLINT
    arb_srcptr fBt;     /* Bt as arb midpoints (A is converted a row at a
                         * time by the worker that owns the row) */
#endif
    int64_t K, S;
    mpfr_t* C;          /* R x S, initialised by the caller at the right bits */
    long lost[MP_PAR_MAX];  /* worst cancellation per chunk */
    long left[MP_PAR_MAX];  /* fewest bits a cancelling sum kept, per chunk */
} DotmJob;

/* Rows [lo, hi) of C. */
static void dotm_rows(void* p, int64_t lo, int64_t hi, int chunk)
{
    DotmJob* j = (DotmJob*)p;
    long worst = 0, fewest = LONG_MAX;
#ifdef USE_FLINT
    arb_ptr row = _arb_vec_init((slong)j->K);
    arb_t acc;
    arb_init(acc);
#endif
    for (int64_t r = lo; r < hi; r++) {
        mpfr_ptr* ar = j->A + r * j->K;
#ifdef USE_FLINT
        for (int64_t k = 0; k < j->K; k++) arf_set_mpfr(arb_midref(row + k), ar[k]);
#endif
        for (int64_t s = 0; s < j->S; s++) {
            mpfr_ptr* bs = j->Bt + s * j->K;
            mpfr_ptr out = j->C[r * j->S + s];
#ifdef USE_FLINT
            arb_approx_dot(acc, NULL, 0, row, 1, j->fBt + s * j->K, 1,
                           (slong)j->K, (slong)mpfr_get_prec(out));
            arf_get_mpfr(out, arb_midref(acc), MPFR_RNDN);
#else
            mpfr_dot(out, ar, bs, (unsigned long)j->K, MPFR_RNDN);
#endif
            /* leading bits lost: the largest product's exponent against the
             * sum's (all of them when the sum is exactly zero) */
            long top = LONG_MIN;
            for (int64_t k = 0; k < j->K; k++) {
                if (!mpfr_regular_p(ar[k]) || !mpfr_regular_p(bs[k])) continue;
                long ex = (long)mpfr_get_exp(ar[k]) + (long)mpfr_get_exp(bs[k]);
                if (ex > top) top = ex;
            }
            if (top == LONG_MIN) continue;
            long l = mpfr_zero_p(out) ? (long)mpfr_get_prec(out)
                   : mpfr_regular_p(out) ? top - (long)mpfr_get_exp(out) : 0;
            if (l > worst) worst = l;
            if (l > 0) {
                long kept = mpfr_zero_p(out) ? 0 : (long)mpfr_get_prec(out) - l;
                if (kept < fewest) fewest = kept;
            }
        }
    }
#ifdef USE_FLINT
    arb_clear(acc);
    _arb_vec_clear(row, (slong)j->K);
#endif
    j->lost[chunk] = worst;
    j->left[chunk] = fewest;
}

Expr* dot_mpfr(const Expr* a, const Expr* b)
{
    DotmOp A, B;
    if (!dotm_read(a, true, &A)) return NULL;
    if (!dotm_read(b, false, &B)) { free(A.leaf); return NULL; }
    Expr* result = NULL;
    int64_t R = A.rows, K = A.cols, S = B.cols;
    if (K != B.rows || !(A.all_mpfr || B.all_mpfr)) goto done;
    {
        /* precisions: the smaller MPFR precision, as Times / Plus give */
        long* rbits = (long*)malloc(sizeof(long) * (size_t)R);
        long* cbits = (long*)malloc(sizeof(long) * (size_t)S);
        long maxbits = 0;
        for (int64_t r = 0; r < R; r++) rbits[r] = LONG_MAX;
        for (int64_t s = 0; s < S; s++) cbits[s] = LONG_MAX;
        for (int64_t r = 0; r < R; r++)
            for (int64_t k = 0; k < K; k++) {
                const Expr* x = A.leaf[r * K + k];
                if (x->type != EXPR_MPFR) continue;
                long p = (long)mpfr_get_prec(x->data.mpfr);
                if (p < rbits[r]) rbits[r] = p;
                if (p > maxbits) maxbits = p;
            }
        for (int64_t k = 0; k < K; k++)
            for (int64_t s = 0; s < S; s++) {
                const Expr* x = B.leaf[k * S + s];
                if (x->type != EXPR_MPFR) continue;
                long p = (long)mpfr_get_prec(x->data.mpfr);
                if (p < cbits[s]) cbits[s] = p;
                if (p > maxbits) maxbits = p;
            }

        /* exact leaves become MPFR numbers at the widest precision in play */
        size_t nexact = 0;
        for (int64_t i = 0; i < R * K; i++) if (A.leaf[i]->type != EXPR_MPFR) nexact++;
        for (int64_t i = 0; i < K * S; i++) if (B.leaf[i]->type != EXPR_MPFR) nexact++;
        mpfr_t* ex = nexact ? (mpfr_t*)malloc(sizeof(mpfr_t) * nexact) : NULL;
        mpfr_ptr* Ap = (mpfr_ptr*)malloc(sizeof(mpfr_ptr) * (size_t)(R * K));
        mpfr_ptr* Bt = (mpfr_ptr*)malloc(sizeof(mpfr_ptr) * (size_t)(K * S));
        mpfr_t* C = (mpfr_t*)malloc(sizeof(mpfr_t) * (size_t)(R * S));
        size_t ei = 0;
        mpfr_t im;
        mpfr_init2(im, 2);
        for (int64_t i = 0; i < R * K; i++) {
            const Expr* x = A.leaf[i];
            if (x->type == EXPR_MPFR) { Ap[i] = (mpfr_ptr)x->data.mpfr; continue; }
            mpfr_init2(ex[ei], (mpfr_prec_t)maxbits);
            get_approx_mpfr(x, ex[ei], im, NULL);
            Ap[i] = ex[ei++];
        }
        for (int64_t k = 0; k < K; k++)
            for (int64_t s = 0; s < S; s++) {
                const Expr* x = B.leaf[k * S + s];
                if (x->type == EXPR_MPFR) { Bt[s * K + k] = (mpfr_ptr)x->data.mpfr; continue; }
                mpfr_init2(ex[ei], (mpfr_prec_t)maxbits);
                get_approx_mpfr(x, ex[ei], im, NULL);
                Bt[s * K + k] = ex[ei++];
            }
        mpfr_clear(im);
        for (int64_t r = 0; r < R; r++)
            for (int64_t s = 0; s < S; s++) {
                /* one operand is all MPFR, so at least one side is finite */
                long p = rbits[r] < cbits[s] ? rbits[r] : cbits[s];
                mpfr_init2(C[r * S + s], (mpfr_prec_t)p);
            }

        bool swap = R == 1 && S > 1;
#ifdef USE_FLINT
        /* the right-hand grid, shared by every worker */
        int64_t nb = swap ? R * K : K * S;
        arb_ptr fB = _arb_vec_init((slong)nb);
        for (int64_t i = 0; i < nb; i++) arf_set_mpfr(arb_midref(fB + i), swap ? Ap[i] : Bt[i]);
#endif
        DotmJob job;
        if (swap) {
            /* vector . matrix: split over the columns instead, as the rows of
             * B^T (Bt is exactly B^T row-major) times the vector */
            job.A = Bt; job.Bt = Ap;
            job.K = K; job.S = 1;
        } else {
            job.A = Ap; job.Bt = Bt;
            job.K = K; job.S = S;
        }
#ifdef USE_FLINT
        job.fBt = fB;
#endif
        int64_t nrows = swap ? S : R;
        job.C = C;
        double limbs = (double)maxbits / 64.0;
        int nc = mp_parallel_for(dotm_rows, &job, nrows,
                                 (double)nrows * (double)job.S * (double)K * limbs * limbs);
        long worst = 0, fewest = LONG_MAX;
        for (int i = 0; i < nc; i++) {
            if (job.lost[i] > worst) worst = job.lost[i];
            if (job.left[i] < fewest) fewest = job.left[i];
        }
        if (worst > 0) numeric_note_cancellation_bits(worst, fewest);
#ifdef USE_FLINT
        _arb_vec_clear(fB, (slong)nb);
#endif

        /* shape: vector.vector -> scalar; matrix.vector / vector.matrix ->
         * vector; matrix.matrix -> matrix */
        if (A.rank == 1 && B.rank == 1) {
            result = expr_new_mpfr_move(C[0]);
        } else if (A.rank == 1 || B.rank == 1) {
            int64_t len = A.rank == 1 ? S : R;
            Expr** el = (Expr**)malloc(sizeof(Expr*) * (size_t)len);
            for (int64_t i = 0; i < len; i++) el[i] = expr_new_mpfr_move(C[i]);
            result = expr_new_function(expr_new_symbol(SYM_List), el, (size_t)len);
            free(el);
        } else {
            Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)R);
            for (int64_t r = 0; r < R; r++) {
                Expr** el = (Expr**)malloc(sizeof(Expr*) * (size_t)S);
                for (int64_t s = 0; s < S; s++) el[s] = expr_new_mpfr_move(C[r * S + s]);
                rows[r] = expr_new_function(expr_new_symbol(SYM_List), el, (size_t)S);
                free(el);
            }
            result = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)R);
            free(rows);
        }
        for (size_t i = 0; i < ei; i++) mpfr_clear(ex[i]);
        free(ex);
        free(Ap);
        free(Bt);
        free(C);
        free(rbits);
        free(cbits);
    }
done:
    free(A.leaf);
    free(B.leaf);
    return result;
}

#endif /* USE_MPFR */
