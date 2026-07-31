/* Tensor contraction (Dot) and the shared dot2 helper.
 *
 * dot2 contracts the trailing axis of a with the leading axis of b and is
 * called from MatrixPower and the eigenvalue solver as well as from
 * builtin_dot itself.
 */

#include "linalg.h"
#include "eval.h"
#include "print.h"
#include "sym_names.h"
#include "ndarray.h"
#include "pack.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef USE_LAPACK
#include "lapack.h"   /* pulls in Accelerate.h / cblas.h (CBLAS declarations) */
#include "ndarray_internal.h"   /* nd_parallel_reduce — threaded inner product */

/* Threaded inner product: each chunk is its own cblas_ddot, the partials are
 * summed afterwards.
 *
 * Accelerate's cblas_ddot appears to be single-threaded. A 10^7 inner product
 * touches 160 MB, and one core reached 20 GB/s (7.9 ms) against Mathematica's
 * 5.19 ms == 31 GB/s -- the shape of a memory-bound reduction that has not been
 * split. This is the same treatment nd_full_sum already gives Total, for the
 * same reason, and reuses the same machinery.
 *
 * Summation order therefore differs from a single ddot, as it already does for
 * Total. nd_parallel_reduce runs ONE serial chunk below its threading
 * threshold, so small vectors are bit-identical to the serial path and only
 * large ones -- where floating-point summation order was never pinned down
 * anyway -- reassociate. */
typedef struct { const double* x; const double* y; } nd_ddot_ctx;

static void nd_ddot_reduce(void* c, size_t lo, size_t hi, double* slot) {
    const nd_ddot_ctx* d = (const nd_ddot_ctx*)c;
    slot[0] = cblas_ddot((int)(hi - lo), d->x + lo, 1, d->y + lo, 1);
}

static double nd_ddot_threaded(const double* x, const double* y, size_t n) {
    nd_ddot_ctx c = { x, y };
    double slots[NDARRAY_MAX_THREADS];
    int k = nd_parallel_reduce(n, nd_ddot_reduce, &c, 1, slots);
    double acc = 0.0;
    for (int t = 0; t < k; t++) acc += slots[t];
    return acc;
}

/* Real contraction via the platform BLAS (Apple Accelerate or system CBLAS —
 * multithreaded + vectorized), for every rank combination BLAS covers:
 *
 *     rank 2 . rank 2   ->  dgemm   matrix * matrix
 *     rank 2 . rank 1   ->  dgemv   matrix * vector
 *     rank 1 . rank 2   ->  dgemv   vector * matrix (transposed)
 *     rank 1 . rank 1   ->  ddot    inner product, a scalar
 *
 * Both operands must be float64 NDArrays with compatible inner dims; anything
 * else returns NULL so the caller falls back to ndarray_dot2 (complex dtypes,
 * higher ranks, shape errors, which ndarray_dot2 reports).
 *
 * The three non-dgemm shapes were added 2026-07-30. cblas_ddot has been bound
 * since the LAPACK bridge landed but was reached only by the runtime
 * availability probe: an ordinary `x . y` on two rank-1 buffers fell through to
 * ndarray_dot2's scalar loop, measured at 12.3 ms for 10^7 elements against
 * Mathematica's 5.19 ms. Nothing needs marshalling here -- an NDArray is already
 * a contiguous row-major float64 buffer, which is exactly what BLAS wants. */
static Expr* nd_blas_matmul(const Expr* a, const Expr* b) {
    if (a->data.ndarray.dtype != NDT_FLOAT64 ||
        b->data.ndarray.dtype != NDT_FLOAT64) return NULL;
    int ra = a->data.ndarray.rank, rb = b->data.ndarray.rank;
    const double* A = (const double*)a->data.ndarray.data;
    const double* B = (const double*)b->data.ndarray.data;

    if (ra == 1 && rb == 1) {                          /* inner product */
        int64_t k = a->data.ndarray.dims[0];
        if (k != b->data.ndarray.dims[0]) return NULL;
        /* A scalar, so there is no presentation to inherit. */
        return expr_new_real(nd_ddot_threaded(A, B, (size_t)k));
    }

    if (ra == 2 && rb == 1) {                          /* matrix * vector */
        int64_t m = a->data.ndarray.dims[0], k = a->data.ndarray.dims[1];
        if (k != b->data.ndarray.dims[0]) return NULL;
        double* out = malloc(sizeof(double) * (size_t)m);
        if (!out) return NULL;
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (int)m, (int)k,
                    1.0, A, (int)k, B, 1, 0.0, out, 1);
        int64_t dims[1] = { m };
        return expr_new_ndarray_like(nd_present_src2(a, b), 1, dims, out, NDT_FLOAT64);
    }

    if (ra == 1 && rb == 2) {                          /* vector * matrix */
        int64_t k = a->data.ndarray.dims[0];
        int64_t k2 = b->data.ndarray.dims[0], n = b->data.ndarray.dims[1];
        if (k != k2) return NULL;
        double* out = malloc(sizeof(double) * (size_t)n);
        if (!out) return NULL;
        /* x^T B == (B^T x), so gemv over B with the transpose flag. */
        cblas_dgemv(CblasRowMajor, CblasTrans, (int)k2, (int)n,
                    1.0, B, (int)n, A, 1, 0.0, out, 1);
        int64_t dims[1] = { n };
        return expr_new_ndarray_like(nd_present_src2(a, b), 1, dims, out, NDT_FLOAT64);
    }

    if (ra != 2 || rb != 2) return NULL;
    int64_t m = a->data.ndarray.dims[0], k = a->data.ndarray.dims[1];
    int64_t k2 = b->data.ndarray.dims[0], n = b->data.ndarray.dims[1];
    if (k != k2) return NULL;                 /* shape error: let ndarray_dot2 report */
    double* out = malloc(sizeof(double) * (size_t)m * (size_t)n);
    if (!out) return NULL;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                (int)m, (int)n, (int)k, 1.0, A, (int)k, B, (int)n, 0.0, out, (int)n);
    int64_t dims[2] = { m, n };
    return expr_new_ndarray_like(nd_present_src2(a, b), 2, dims, out, NDT_FLOAT64); /* adopts out */
}
#endif

/* Build a nested-List tensor of shape dims[0..rank-1] by consuming leaves
 * from `flat` in row-major order.  Each leaf is deep-copied. */
static Expr* build_tensor(int64_t* dims, int rank, Expr** flat, size_t* idx) {
    if (rank == 0) {
        return expr_copy(flat[(*idx)++]);
    }
    int64_t len = dims[0];
    Expr** args = NULL;
    if (len > 0) args = malloc(sizeof(Expr*) * len);
    for (int64_t i = 0; i < len; i++) {
        args[i] = build_tensor(dims + 1, rank - 1, flat, idx);
    }
    Expr* res = expr_new_function(expr_new_symbol(SYM_List), args, len);
    if (args) free(args);
    return res;
}

/* Pack a plain numeric List operand so a mixed pair can take the buffer path.
 * Returns a NEW packed Expr the caller owns, or NULL when `list` is not
 * packable at all (a symbolic list, a Rational, a jagged nesting).
 *
 * pack_force rather than pack_offer, i.e. ignoring PACK_MIN_ELEMENTS: the
 * threshold exists to bound blast radius on values too small for the
 * representation to pay for itself, and it decides that by looking at the list
 * IN ISOLATION. Dot is the case where that reasoning does not hold, because the
 * work is O(R*S*K) in BOTH operands -- a 40-element vector is far below the
 * threshold and entirely right to be, but contracting it against a 20000x40
 * matrix is 800,000 multiplies either way. Packing 40 elements to buy those is
 * never a bad trade. */
static Expr* dot_pack_operand(Expr* list) {
    if (!list || list->type != EXPR_FUNCTION) return NULL;
    Expr* cand = pack_force(expr_copy(list), false, NDT_FLOAT64);
    if (cand && cand->type == EXPR_NDARRAY) return cand;
    if (cand) expr_free(cand);
    return NULL;
}

Expr* dot2(Expr* a, Expr* b, bool* error_printed) {
    /* One operand packed and the other an ordinary List: pack the plain one and
     * let the pair take the buffer path.
     *
     * Without this a single unpacked operand disables the fast path for the
     * whole product, and the cost is charged in proportion to the OTHER operand.
     * The measured case is gradient descent -- w = ConstantArray[0., 32] is 32
     * elements, correctly below the packing threshold, and `X . w` for a
     * 20000x40 X then ran the generic Times/Plus-per-element loop at 320 ms
     * against 0.40 ms for the same product with both sides packed. 798x, from a
     * representation decision made about the small operand and paid for by the
     * large one.
     *
     * Guarded by pack_any_created() so a session that never packs anything does
     * not pay the probe, and declining is O(1) -- pack_force's sniff bails at
     * the first element that is not a machine number, so a symbolic Dot costs
     * one pointer test. */
    Expr* tmp_a = NULL;
    Expr* tmp_b = NULL;
    if (pack_any_created()) {
        if (a->type == EXPR_NDARRAY && b->type != EXPR_NDARRAY)
            tmp_b = dot_pack_operand(b);
        else if (b->type == EXPR_NDARRAY && a->type != EXPR_NDARRAY)
            tmp_a = dot_pack_operand(a);
    }
    Expr* fa = tmp_a ? tmp_a : a;
    Expr* fb = tmp_b ? tmp_b : b;

    /* Fast path: both operands are dense NDArray values of rank <= 2 —
     * contract directly over the flat double buffers, no symbolic
     * Times/Plus per element. Higher-rank NDArray operands fall through to
     * the generic tensor path below via a nested-List conversion. */
    if (fa->type == EXPR_NDARRAY && fb->type == EXPR_NDARRAY) {
#ifdef USE_LAPACK
        Expr* bl = nd_blas_matmul(fa, fb);   /* rank-2 f64 matmul -> BLAS dgemm */
        if (bl) { expr_free(tmp_a); expr_free(tmp_b); return bl; }
#endif
        bool shape_error = false;
        Expr* fast = ndarray_dot2(fa, fb, &shape_error);
        if (fast) { expr_free(tmp_a); expr_free(tmp_b); return fast; }
        if (shape_error) {
            if (!*error_printed) {
                char* a_str = expr_to_string_fullform(fa);
                char* b_str = expr_to_string_fullform(fb);
                fprintf(stderr, "Dot::dotsh: Tensors %s and %s have incompatible shapes.\n", a_str, b_str);
                free(a_str);
                free(b_str);
                *error_printed = true;
            }
            expr_free(tmp_a);
            expr_free(tmp_b);
            return NULL;
        }
    }
    /* Fell through: the generic path below reads the ORIGINAL operands, so the
     * packed temporaries are dead here. */
    expr_free(tmp_a);
    expr_free(tmp_b);

    Expr* conv_a = NULL;
    Expr* conv_b = NULL;
    if (a->type == EXPR_NDARRAY) { conv_a = ndarray_to_nested_list(a); a = conv_a; }
    if (b->type == EXPR_NDARRAY) { conv_b = ndarray_to_nested_list(b); b = conv_b; }

    int64_t dimsA[64];
    int rankA = get_tensor_dims(a, dimsA);
    if (rankA <= 0) { // Not a tensor, or jagged
        if (conv_a) expr_free(conv_a);
        if (conv_b) expr_free(conv_b);
        return NULL;
    }

    int64_t dimsB[64];
    int rankB = get_tensor_dims(b, dimsB);
    if (rankB <= 0) {
        if (conv_a) expr_free(conv_a);
        if (conv_b) expr_free(conv_b);
        return NULL;
    }

    int64_t K = dimsA[rankA - 1];
    if (K != dimsB[0]) {
        if (!*error_printed) {
            char* a_str = expr_to_string_fullform(a);
            char* b_str = expr_to_string_fullform(b);
            fprintf(stderr, "Dot::dotsh: Tensors %s and %s have incompatible shapes.\n", a_str, b_str);
            free(a_str);
            free(b_str);
            *error_printed = true;
        }
        if (conv_a) expr_free(conv_a);
        if (conv_b) expr_free(conv_b);
        return NULL;
    }

    int64_t N_A = 1; for(int i=0; i<rankA; i++) N_A *= dimsA[i];
    int64_t N_B = 1; for(int i=0; i<rankB; i++) N_B *= dimsB[i];

    Expr** flatA = NULL; if (N_A > 0) flatA = malloc(sizeof(Expr*) * N_A);
    Expr** flatB = NULL; if (N_B > 0) flatB = malloc(sizeof(Expr*) * N_B);
    size_t idxA = 0; if (N_A > 0) flatten_tensor(a, flatA, &idxA);
    size_t idxB = 0; if (N_B > 0) flatten_tensor(b, flatB, &idxB);

    int64_t R = K == 0 ? N_A : N_A / K;
    int64_t S = K == 0 ? N_B : N_B / K;

    Expr** flatC = NULL;
    if (R * S > 0) flatC = malloc(sizeof(Expr*) * (R * S));

    for (int64_t r = 0; r < R; r++) {
        for (int64_t s = 0; s < S; s++) {
            Expr** sum_args = NULL;
            if (K > 0) sum_args = malloc(sizeof(Expr*) * K);
            for (int64_t k = 0; k < K; k++) {
                Expr* a_elem = flatA[r * K + k];
                Expr* b_elem = flatB[k * S + s];
                Expr* t_args[2] = { expr_copy(a_elem), expr_copy(b_elem) };
                sum_args[k] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), t_args, 2));
            }
            if (K == 0) {
                flatC[r * S + s] = expr_new_integer(0);
            } else if (K == 1) {
                flatC[r * S + s] = sum_args[0];
            } else {
                flatC[r * S + s] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), sum_args, K));
            }
            if (sum_args) free(sum_args);
        }
    }

    if (flatA) { for(size_t i=0; i<idxA; i++) expr_free(flatA[i]); free(flatA); }
    if (flatB) { for(size_t i=0; i<idxB; i++) expr_free(flatB[i]); free(flatB); }

    int64_t dimsC[64];
    int rankC = rankA + rankB - 2;
    for (int i = 0; i < rankA - 1; i++) dimsC[i] = dimsA[i];
    for (int i = 0; i < rankB - 1; i++) dimsC[rankA - 1 + i] = dimsB[i + 1];

    size_t c_idx = 0;
    Expr* result = build_tensor(dimsC, rankC, flatC, &c_idx);
    if (flatC) {
        for (int64_t i = 0; i < R * S; i++) expr_free(flatC[i]);
        free(flatC);
    }

    if (conv_a) expr_free(conv_a);
    if (conv_b) expr_free(conv_b);
    return result;
}

Expr* builtin_dot(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t count = res->data.function.arg_count;
    if (count == 0) return NULL;
    if (count == 1) return expr_copy(res->data.function.args[0]);

    Expr** new_args = malloc(sizeof(Expr*) * count);
    for (size_t i=0; i<count; i++) new_args[i] = expr_copy(res->data.function.args[i]);
    size_t new_count = count;

    bool changed = false;
    bool error_printed = false;
    for (size_t i = 0; i < new_count - 1; i++) {
        Expr* a = new_args[i];
        Expr* b = new_args[i+1];

        int64_t dA[64], dB[64];
        bool a_is_tensor = (a->type == EXPR_NDARRAY) || get_tensor_dims(a, dA) > 0;
        bool b_is_tensor = (b->type == EXPR_NDARRAY) || get_tensor_dims(b, dB) > 0;
        if (a_is_tensor && b_is_tensor) {
            Expr* d = dot2(a, b, &error_printed);
            if (d) {
                expr_free(a);
                expr_free(b);
                new_args[i] = d;
                for (size_t j = i + 2; j < new_count; j++) {
                    new_args[j - 1] = new_args[j];
                }
                new_count--;
                changed = true;
                i--; // re-check the new element at index i with i+1 (which shifted)
                if (i == (size_t)-1) i = 0; // In case we backed up past 0 (wait, i-- makes i=2^64-1, then i++ in loop makes it 0, correct)
            } else if (error_printed) {
                for (size_t j=0; j<new_count; j++) expr_free(new_args[j]);
                free(new_args);
                return NULL;
            }
        }
    }

    if (!changed) {
        for (size_t j=0; j<new_count; j++) expr_free(new_args[j]);
        free(new_args);
        return NULL;
    }

    Expr* final_res;
    if (new_count == 1) {
        final_res = new_args[0];
    } else {
        final_res = expr_new_function(expr_new_symbol(SYM_Dot), new_args, new_count);
    }

    free(new_args);
    return final_res;
}
