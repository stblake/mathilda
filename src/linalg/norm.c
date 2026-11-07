/* Vector and matrix norms.
 *
 *   Norm[x]       scalar  -> Abs[x]
 *   Norm[v]       vector  -> 2-norm
 *   Norm[v, p]    vector  -> p-norm, also Infinity
 *   Norm[m, "Frobenius"]  -> Frobenius norm (works on any tensor rank)
 *
 *   Norm[m] / Norm[m, 2]  -> the spectral norm, the largest singular value:
 *       machine   LAPACK (ndla_matrix_norm_direct),
 *       MPFR      svd_mpfr_values at the input precision,
 *       exact / symbolic   Sqrt of the largest eigenvalue of the smaller Gram
 *                 matrix m.m^H or m^H.m, as the reference answers it:
 *                 Norm[{{4, 1}, {0, 5}}] is Sqrt[21 + Sqrt[41]], a 3 x 3
 *                 gives Sqrt[Root[...]], and a symbolic matrix a Max of the
 *                 Sqrt of every eigenvalue.
 */

#include "linalg.h"
#include "ndlinalg.h"
#include "linsolve.h"  /* matsol_is_inexact — the exact/inexact gate below */
#include "ndarray.h"   /* ndarray_from_nested_list, NDT_* */
#include "common.h"    /* common_scan_inexact — machine vs MPFR gate below */
#include "eval.h"
#include "sym_names.h"
#include "la_dense.h"         /* la_all_numeric */
#include "svdecomp_internal.h" /* svd_mpfr_values: the MPFR spectral norm */
#ifdef USE_MPFR
#include <mpfr.h>
#endif
#include <stdlib.h>
#include <string.h>

/* True iff p is a norm order the ND vector fast path (ndla_norm) handles:
 * absent (default 2), a positive Integer/Real, or Infinity.  The "Frobenius"
 * string and everything else keep the symbolic path. */
static bool norm_p_is_nd_compatible(const Expr* p) {
    if (!p) return true;
    if (p->type == EXPR_SYMBOL && p->data.symbol.name == SYM_Infinity) return true;
    if (p->type == EXPR_INTEGER && p->data.integer > 0) return true;
    if (p->type == EXPR_REAL && p->data.real > 0.0) return true;
    return false;
}

Expr* builtin_norm(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count != 1 && res->data.function.arg_count != 2)) return NULL;
    if (linalg_call_has_ndarray(res)) return ndla_norm(res);

    Expr* expr = res->data.function.args[0];
    Expr* p = NULL;
    if (res->data.function.arg_count == 2) {
        p = res->data.function.args[1];
    }

    int64_t dims[64];
    int rank = get_tensor_dims(expr, dims);
    if (rank < 0) return NULL; // jagged array

    if (rank == 0) {
        // Scalar: Norm[x] -> Abs[x]
        if (!p) {
            Expr* args[1] = { expr_copy(expr) };
            return eval_and_free(expr_new_function(expr_new_symbol(SYM_Abs), args, 1));
        }
        return NULL;
    }

    /* Machine-precision vector fast path.  The element-wise Power[Abs[x], p]
     * pipeline below overflows (Abs[1e200]^2 -> inf) for well-scaled inputs
     * whose true norm is representable, so Norm[{1e200, 1e200}] returned inf.
     * Route a machine-numeric vector through ndla_norm, whose scaled 2-norm /
     * p-norm are overflow-safe.  Restricted to machine precision (min_bits <=
     * 53): an EXACT vector must keep its exact symbolic answer (Norm[{1,2}] ==
     * Sqrt[5]), and an MPFR vector is already safe on the symbolic path (its
     * wide exponent keeps the sum finite) and must not be narrowed to float64. */
    if (rank == 1 && norm_p_is_nd_compatible(p)) {
        CommonInexactInfo ii = common_scan_inexact(expr);
        if (ii.has_inexact && ii.min_bits <= 53) {
            Expr* packed = ndarray_from_nested_list(expr, NDT_FLOAT64);
            if (!packed) packed = ndarray_from_nested_list(expr, NDT_COMPLEX64);
            if (packed) {
                Expr* call = p
                    ? expr_new_function(expr_new_symbol(SYM_Norm),
                          (Expr*[]){ packed, expr_copy(p) }, 2)
                    : expr_new_function(expr_new_symbol(SYM_Norm),
                          (Expr*[]){ packed }, 1);
                Expr* out = ndla_norm(call);
                expr_free(call);
                if (out) return out;
            }
        }
    }

    if (rank == 1 || (rank >= 2 && p && p->type == EXPR_STRING && strcmp(p->data.string, "Frobenius") == 0)) {
        int64_t N = 1;
        for (int i = 0; i < rank; i++) N *= dims[i];

        Expr** flat = NULL;
        if (N > 0) {
            flat = malloc(sizeof(Expr*) * N);
            size_t idx = 0;
            flatten_tensor(expr, flat, &idx);
        }

        Expr* result = NULL;

        if (p && p->type == EXPR_SYMBOL && p->data.symbol.name == SYM_Infinity) {
            if (N == 0) {
                result = expr_new_integer(0);
            } else {
                Expr** max_args = malloc(sizeof(Expr*) * N);
                for (int64_t i = 0; i < N; i++) {
                    max_args[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Abs), (Expr*[]){expr_copy(flat[i])}, 1));
                }
                result = eval_and_free(expr_new_function(expr_new_symbol(SYM_Max), max_args, N));
                free(max_args);
            }
        } else {
            Expr* norm_p = NULL;
            if (!p || (p->type == EXPR_STRING && strcmp(p->data.string, "Frobenius") == 0)) {
                norm_p = expr_new_integer(2);
            } else {
                norm_p = expr_copy(p);
            }

            if (N == 0) {
                result = expr_new_integer(0);
            } else {
                Expr** plus_args = malloc(sizeof(Expr*) * N);
                for (int64_t i = 0; i < N; i++) {
                    Expr* abs_val = eval_and_free(expr_new_function(expr_new_symbol(SYM_Abs), (Expr*[]){expr_copy(flat[i])}, 1));
                    plus_args[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){abs_val, expr_copy(norm_p)}, 2));
                }
                Expr* sum = NULL;
                if (N == 1) {
                    sum = plus_args[0];
                } else {
                    sum = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), plus_args, N));
                }
                free(plus_args);

                Expr* inv_p = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){expr_copy(norm_p), expr_new_integer(-1)}, 2));
                result = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){sum, inv_p}, 2));
            }
            expr_free(norm_p);
        }

        if (flat) {
            for (int64_t i = 0; i < N; i++) expr_free(flat[i]);
            free(flat);
        }
        return result;
    }

    /*
     * Matrix norms other than Frobenius -- the induced 1-, 2- and Infinity-norms
     * -- have no symbolic implementation here, but ndla_norm does have one
     * (LAPACK dlange/zlange, and an SVD for the 2-norm). So rather than stay
     * unevaluated, pack a machine-numeric matrix and use it.
     *
     * This closed an observable inconsistency rather than adding a feature:
     * once ordinary Lists pack automatically, Norm[m] on a 300x300 matrix took
     * the ND path and answered, while the SAME value written small enough to
     * stay unpacked answered Norm[{{...}}] unevaluated. A differential sweep
     * (every head, packed against MATHILDA_NO_PACK=1) is what surfaced it; the
     * values agree with Mathematica to the last digit on the five test matrices.
     *
     * A symbolic or exact-but-non-machine matrix still returns NULL and stays
     * symbolic, which is the documented behaviour for those.
     */
    /* INEXACT matrices only. Norm of an EXACT matrix is an exact algebraic
     * number -- Mathematica answers Norm[{{1,2},{4,5}}] with
     * Sqrt[23 + 2 Sqrt[130]] -- and Mathilda has no symbolic SVD to produce it,
     * so the honest answer for exact input is to stay unevaluated, as before.
     * Packing it to float64 and returning 6.76783 would be a machine-real answer
     * to an exact question; normalize_tests caught exactly that, via
     * Normalize[{{1,2},{4,5}}, Norm]. */
#ifdef USE_MPFR
    /* Arbitrary precision: the largest singular value at the input precision
     * (packing to float64 below would silently drop it to machine). */
    /* Norm[m] and Norm[m, 2] alike: the spectral norm (an explicit 2 used to
     * fall through to the float64 packing below, which refuses MPFR entries,
     * and the call stayed unevaluated) */
    bool spectral = !p || (p->type == EXPR_INTEGER && p->data.integer == 2);
    if (rank == 2 && spectral && dims[0] > 0 && dims[1] > 0 && la_all_numeric(expr)) {
        CommonInexactInfo info = common_scan_inexact(expr);
        if (info.has_inexact && info.min_bits > 53) {
            int n = (int)dims[0], q = (int)dims[1];
            Expr* work = expr_copy(expr);
            if (q > n) {          /* the singular values of m^H are m's */
                work = eval_and_free(expr_new_function(
                    expr_new_symbol(SYM_ConjugateTranspose), (Expr*[]){ work }, 1));
                int t = n; n = q; q = t;
            }
            mpfr_t* sig = svd_mpfr_values(work, n, q, (mpfr_prec_t)info.min_bits);
            expr_free(work);
            if (sig) {
                Expr* out = expr_new_mpfr_copy(sig[0]);
                svd_mpfr_values_free(sig, q);
                return out;
            }
        }
    }
#endif
    if (rank == 2 && matsol_is_inexact(expr)) {
        Expr* packed = ndarray_from_nested_list(expr, NDT_FLOAT64);
        if (!packed) packed = ndarray_from_nested_list(expr, NDT_COMPLEX64);
        if (packed) {
            /* ndla_matrix_norm_direct, NOT ndla_norm: ndla_norm's fallback is
             * linalg_delist_and_reeval, which would rebuild Norm[plainMatrix]
             * and evaluate it -- straight back into this branch, forever. The
             * direct entry returns NULL instead of deferring, so the recursion
             * cannot form. */
            Expr* out = ndla_matrix_norm_direct(packed, p);
            expr_free(packed);
            if (out) return out;
        }
    }

    /* The induced 1- and Infinity-norms need no SVD: the largest absolute
     * column sum and row sum. Exact and symbolic matrices get them here
     * (Norm[{{1, 2}, {-3, 4}}, 1] is 6, Norm[{{a, 0}, {0, b}}, 1] is
     * Max[Abs[a], Abs[b]]). */
    if (rank == 2 && dims[0] > 0 && dims[1] > 0 && p &&
        ((p->type == EXPR_INTEGER && p->data.integer == 1) ||
         (p->type == EXPR_SYMBOL && p->data.symbol.name == SYM_Infinity))) {
        bool rows = p->type == EXPR_SYMBOL;
        Expr* absm = eval_and_free(expr_new_function(expr_new_symbol(SYM_Abs),
                                                      (Expr*[]){ expr_copy(expr) }, 1));
        Expr* tot = rows
            ? expr_new_function(expr_new_symbol(SYM_Total),
                  (Expr*[]){ absm, expr_new_function(expr_new_symbol(SYM_List),
                                                      (Expr*[]){ expr_new_integer(2) }, 1) }, 2)
            : expr_new_function(expr_new_symbol(SYM_Total), (Expr*[]){ absm }, 1);
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Max), (Expr*[]){ tot }, 1));
    }

    /* The spectral norm of an exact or symbolic matrix: Sqrt of the largest
     * eigenvalue of the smaller Gram matrix (m.m^H when m is wide or square,
     * else m^H.m).  The Gram matrix is positive semidefinite, so for a numeric
     * matrix the eigenvalue of largest magnitude (Eigenvalues[g, 1]) is the
     * largest; a symbolic one has no order, so the answer is the Max over
     * all of them. */
    if (rank == 2 && dims[0] > 0 && dims[1] > 0 &&
        (!p || (p->type == EXPR_INTEGER && p->data.integer == 2)) &&
        !matsol_is_inexact(expr)) {
        Expr* ct = eval_and_free(expr_new_function(expr_new_symbol(SYM_ConjugateTranspose),
                                                   (Expr*[]){ expr_copy(expr) }, 1));
        Expr* g = dims[0] <= dims[1]
            ? expr_new_function(expr_new_symbol(SYM_Dot), (Expr*[]){ expr_copy(expr), ct }, 2)
            : expr_new_function(expr_new_symbol(SYM_Dot), (Expr*[]){ ct, expr_copy(expr) }, 2);
        g = eval_and_free(g);
        bool numeric = la_all_numeric(expr);
        Expr* ev = numeric
            ? expr_new_function(expr_new_symbol(SYM_Eigenvalues), (Expr*[]){ g, expr_new_integer(1) }, 2)
            : expr_new_function(expr_new_symbol(SYM_Eigenvalues), (Expr*[]){ g }, 1);
        ev = eval_and_free(ev);
        if (!head_is(ev, SYM_List) || ev->data.function.arg_count == 0) {
            expr_free(ev);
            return NULL;
        }
        Expr* out;
        if (numeric) {
            Expr* top = expr_copy(ev->data.function.args[0]);
            expr_free(ev);
            /* the reference's form of a quadratic surd: Sqrt[21 + Sqrt[41]],
             * Sqrt[(91 + Sqrt[8065])/2] -- one fraction, content cancelled */
            top = eval_and_free(expr_new_function(expr_new_symbol(SYM_Together),
                                                  (Expr*[]){ top }, 1));
            out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Sqrt), (Expr*[]){ top }, 1));
        } else {
            Expr* sq = eval_and_free(expr_new_function(expr_new_symbol(SYM_Sqrt), (Expr*[]){ ev }, 1));
            out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Apply),
                                    (Expr*[]){ expr_new_symbol(SYM_Max), sq }, 2));
        }
        return out;
    }

    return NULL;
}
