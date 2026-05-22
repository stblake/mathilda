/* ludecomp.c
 *
 * LUDecomposition[m] -- {lu, p, c} Doolittle factorisation with row
 * pivoting.
 *
 * Strategy.  One algorithmic core - Doolittle's algorithm with partial
 * pivoting driven through the Mathilda evaluator - serves every input
 * family:
 *
 *   - exact integer / rational / complex / free-symbolic matrices
 *     run the pipeline as-is.  Pivoting is zero-only (advance to the
 *     next non-zero pivot if the natural choice is provably zero).
 *
 *   - inexact matrices (Real or MPFR leaves) follow the
 *     rationalise -> exact pipeline -> numericalise round-trip used by
 *     PseudoInverse / Eigenvalues / QRDecomposition.  Output precision
 *     matches the minimum precision present in the input.
 *
 * Conventions.  We work internally on a flat n x n row-major buffer:
 *
 *     LU[i * n + k] = (i, k) entry of the combined factorisation
 *
 * with the strict lower triangle storing L (unit diagonal implicit)
 * and the upper triangle storing U.  perm[k] is the 1-indexed original
 * row that was placed at row k by the pivoting steps, so the public
 * identity m[[perm]] == l . u holds:
 *
 *     l = LowerTriangularize[lu, -1] + IdentityMatrix[n]
 *     u = UpperTriangularize[lu]
 *
 * For exact / symbolic inputs the condition-number slot is the exact
 * Integer 0 (matching the Mathematica example).  The machine and MPFR
 * kernels emit a Real / MPFR L-infinity estimate; in this top-level
 * file we only need the exact-zero default.
 *
 * Memory contract.  Standard builtin contract.  This file does NOT
 * call expr_free(res) - the evaluator owns `res` and frees it on a
 * non-NULL return (MEMORY.md / SPEC.md §4.1).  Every intermediate
 * allocation is matched by a free along every exit path.
 */

#include "ludecomp.h"
#include "ludecomp_internal.h"
#include "linalg.h"
#include "linsolve.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"
#include "poly.h"
#include "expand.h"
#include "sym_names.h"
#include "common.h"

#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ *
 *  Small evaluator helpers (lifted from qrdecomp.c -- kept local so   *
 *  the two modules stay independent translation units).               *
 * ------------------------------------------------------------------ */
static Expr* eval_plus(Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(
        expr_new_symbol("Plus"), (Expr*[]){a, b}, 2));
}
static Expr* eval_times(Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(
        expr_new_symbol("Times"), (Expr*[]){a, b}, 2));
}
static Expr* eval_power(Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(
        expr_new_symbol("Power"), (Expr*[]){a, b}, 2));
}
static Expr* eval_together(Expr* a) {
    return eval_and_free(expr_new_function(
        expr_new_symbol("Together"), (Expr*[]){a}, 1));
}

/* True iff e canonicalises to zero.  Mirrors qrdecomp.c's
 * is_definitely_zero: Together first, then the polynomial zero test. */
static bool is_definitely_zero(Expr* e) {
    Expr* simp = eval_together(expr_copy(e));
    bool z = is_zero_poly(simp);
    expr_free(simp);
    return z;
}

/* ------------------------------------------------------------------ *
 *  Doolittle core (symbolic).                                          *
 *                                                                       *
 *  Pivoting rule: at step k, scan rows k..n-1 of column k for the      *
 *  first row whose entry is NOT provably zero.  If none exists, leave  *
 *  the zero in place, mark singular, and continue (the matrix is        *
 *  singular but the factorisation still completes -- matching          *
 *  Mathematica's LUDecomposition::sing behaviour).                      *
 *                                                                       *
 *  Update step (after pivot is in place at LU[k, k]):                  *
 *     for i in [k+1, n):                                                *
 *         LU[i, k] = LU[i, k] / pivot               (L entry)           *
 *         for j in [k+1, n):                                            *
 *             LU[i, j] = LU[i, j] - LU[i, k] * LU[k, j]                *
 *                                                                       *
 *  Each arithmetic primitive goes through the Mathilda evaluator so    *
 *  symbolic / rational / Sqrt entries all just work.                    *
 * ------------------------------------------------------------------ */
bool lu_symbolic_core(Expr** A_flat, int n,
                      Expr*** out_LU_flat, int** out_perm,
                      bool* out_singular)
{
    /* Working LU: deep copy of A so we never alias the caller's
     * entries; we mutate freely in place. */
    Expr** LU = (Expr**)malloc(sizeof(Expr*) * (size_t)n * (size_t)n);
    for (int t = 0; t < n * n; t++) LU[t] = expr_copy(A_flat[t]);

    int* perm = (int*)malloc(sizeof(int) * (size_t)n);
    for (int k = 0; k < n; k++) perm[k] = k + 1;   /* 1-indexed identity */

    bool singular = false;

    for (int k = 0; k < n; k++) {
        /* Pivot selection: first non-zero entry in column k from row k
         * down.  For symbolic / exact inputs this is the conservative
         * choice -- a magnitude-based rule isn't meaningful when
         * entries contain free variables.  Matches the spec example
         *     LUDecomposition[{{a, b}, {c, d}}] -> perm = {1, 2}. */
        int pivot_row = -1;
        for (int i = k; i < n; i++) {
            if (!is_definitely_zero(LU[i * n + k])) { pivot_row = i; break; }
        }
        if (pivot_row < 0) {
            /* Whole column from row k is zero -- singular at this stage.
             * Leave LU[k, k] = 0 and continue. */
            singular = true;
            continue;
        }
        if (pivot_row != k) {
            /* Swap rows pivot_row and k in LU. */
            for (int j = 0; j < n; j++) {
                Expr* tmp = LU[k * n + j];
                LU[k * n + j] = LU[pivot_row * n + j];
                LU[pivot_row * n + j] = tmp;
            }
            int tp = perm[k]; perm[k] = perm[pivot_row]; perm[pivot_row] = tp;
        }

        /* Eliminate.  pivot is the upper-triangle diagonal entry; the
         * L entries below it are stored as ratios. */
        Expr* pivot = LU[k * n + k];        /* borrowed */
        for (int i = k + 1; i < n; i++) {
            Expr* l_ik = eval_times(expr_copy(LU[i * n + k]),
                                    eval_power(expr_copy(pivot),
                                               expr_new_integer(-1)));
            expr_free(LU[i * n + k]);
            LU[i * n + k] = l_ik;          /* now owns the L entry */

            for (int j = k + 1; j < n; j++) {
                Expr* prod = eval_times(expr_copy(LU[i * n + k]),
                                        expr_copy(LU[k * n + j]));
                Expr* neg  = eval_times(expr_new_integer(-1), prod);
                Expr* new_v = eval_plus(LU[i * n + j], neg); /* consumes LHS */
                LU[i * n + j] = new_v;
            }
        }
    }

    *out_LU_flat  = LU;
    *out_perm     = perm;
    *out_singular = singular;
    return true;
}

/* ------------------------------------------------------------------ *
 *  Wrap a flat n*n row-major buffer into a List[List[...]].  Steals    *
 *  each entry: the caller must not free buf's entries (only the buf    *
 *  array itself).                                                       *
 * ------------------------------------------------------------------ */
static Expr* wrap_matrix(Expr** buf, int n) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) {
        Expr** elems = NULL;
        if (n > 0) elems = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int j = 0; j < n; j++) elems[j] = buf[i * n + j];   /* steal */
        rows[i] = expr_new_function(
            expr_new_symbol("List"), elems, (size_t)n);
        if (elems) free(elems);
    }
    Expr* out = expr_new_function(
        expr_new_symbol("List"), rows, (size_t)n);
    free(rows);
    return out;
}

/* Wrap a 1-indexed length-n int array into a Mathilda List of
 * Integers. */
static Expr* wrap_perm(const int* perm, int n) {
    Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int k = 0; k < n; k++) elems[k] = expr_new_integer(perm[k]);
    Expr* out = expr_new_function(
        expr_new_symbol("List"), elems, (size_t)n);
    free(elems);
    return out;
}

/* Element-wise Together to canonicalise small cancellations in the LU
 * matrix.  Together is Listable so threading is automatic. */
static Expr* tidy_matrix(Expr* m) {
    return eval_and_free(expr_new_function(
        expr_new_symbol("Together"), (Expr*[]){expr_copy(m)}, 1));
}

/* One-shot singular warning. */
static void lu_warn_singular_once(uint64_t* counter, Expr* m)
{
    if (*counter) return;
    *counter = 1;
    char* s = expr_to_string(m);
    fprintf(stderr,
        "LUDecomposition::sing: Matrix %s is singular.\n", s);
    free(s);
}

/* ------------------------------------------------------------------ *
 *  Symbolic kernel dispatcher.                                          *
 * ------------------------------------------------------------------ */
Expr* lu_symbolic_dispatch(Expr* m, int n)
{
    static uint64_t sing_warn_counter = 0;

    /* Inexact-preprocessing pipeline (same shape as PseudoInverse /
     * QRDecomposition). */
    CommonInexactInfo info = common_scan_inexact(m);
    Expr* m_rat = NULL;
    long prec_bits = 53;
    Expr* matrix_to_use = m;
    if (info.has_inexact) {
        prec_bits = info.min_bits ? info.min_bits : 53;
        m_rat = common_rationalize_input(m, prec_bits);
        matrix_to_use = m_rat;
    }

    /* Flatten m (row-major) into Expr** for the Doolittle worker. */
    Expr** A_flat = (Expr**)malloc(sizeof(Expr*) * (size_t)n * (size_t)n);
    {
        size_t idx = 0;
        flatten_tensor(matrix_to_use, A_flat, &idx);
    }

    Expr** LU_flat = NULL;
    int*   perm    = NULL;
    bool   singular = false;
    bool ok = lu_symbolic_core(A_flat, n, &LU_flat, &perm, &singular);

    for (int t = 0; t < n * n; t++) expr_free(A_flat[t]);
    free(A_flat);
    if (m_rat) expr_free(m_rat);

    if (!ok) {
        if (LU_flat) {
            for (int t = 0; t < n * n; t++) expr_free(LU_flat[t]);
            free(LU_flat);
        }
        if (perm) free(perm);
        return NULL;
    }

    if (singular) lu_warn_singular_once(&sing_warn_counter, m);

    Expr* lu = wrap_matrix(LU_flat, n);
    free(LU_flat);
    Expr* p  = wrap_perm(perm, n);
    free(perm);

    Expr* lu_t = tidy_matrix(lu); expr_free(lu); lu = lu_t;

    /* Numericalise lu back to the input precision when we rationalised.
     * For symbolic / exact inputs we never rationalised, so this is a
     * no-op.  The condition-number slot is the exact Integer 0 for
     * symbolic input -- matching Mathematica's example -- and for the
     * fall-through-from-fast-kernel inexact case the symbolic kernel
     * also returns 0 (the upper-level dispatcher recovers a real
     * condition number when the LAPACK / MPFR kernels succeed). */
    Expr* c = expr_new_integer(0);
    if (info.has_inexact) {
        Expr* lu_n = common_numericalize_result(lu, prec_bits);
        expr_free(lu); lu = lu_n;
        /* Leave c as exact Integer 0 to flag "no estimate available"
         * from the symbolic fallback. */
    }

    Expr** items = (Expr**)malloc(sizeof(Expr*) * 3);
    items[0] = lu; items[1] = p; items[2] = c;
    Expr* result = expr_new_function(expr_new_symbol("List"), items, 3);
    free(items);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Top-level kernel router.                                             *
 *                                                                       *
 *  Inexact, min_bits <= 53  -> machine LAPACK kernel.                   *
 *  Inexact, min_bits > 53   -> MPFR kernel.                              *
 *  Anything else            -> symbolic dispatcher.                     *
 *                                                                       *
 *  On any soft failure (USE_LAPACK=0, USE_MPFR=0, non-numeric leaf,     *
 *  fatal LAPACK info, etc.) the chosen kernel returns NULL and we fall  *
 *  through to the symbolic dispatcher -- which understands the inexact  *
 *  input via the rationalise / numericalise round-trip.                 *
 * ------------------------------------------------------------------ */
Expr* lu_dispatch(Expr* m, int n)
{
    CommonInexactInfo info = common_scan_inexact(m);
    if (info.has_inexact) {
        if (info.min_bits <= 53) {
            Expr* fast = lu_machine_dispatch(m, n);
            if (fast) return fast;
        } else {
            Expr* fast = lu_mpfr_dispatch(m, n);
            if (fast) return fast;
        }
    }
    return lu_symbolic_dispatch(m, n);
}

/* ------------------------------------------------------------------ *
 *  Public entry.                                                        *
 * ------------------------------------------------------------------ */
Expr* builtin_ludecomposition(Expr* res)
{
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return NULL;

    Expr* m = res->data.function.args[0];

    int64_t dims[64];
    int trank = get_tensor_dims(m, dims);
    if (trank != 2 || dims[0] == 0 || dims[1] == 0 || dims[0] != dims[1]) {
        char* s = expr_to_string(m);
        fprintf(stderr,
                "LUDecomposition::matsq: Argument %s at position 1 is "
                "not a non-empty square matrix.\n", s);
        free(s);
        return NULL;
    }
    int n = (int)dims[0];

    return lu_dispatch(m, n);
}

void ludecomp_init(void)
{
    symtab_add_builtin("LUDecomposition", builtin_ludecomposition);
    symtab_get_def("LUDecomposition")->attributes |= ATTR_PROTECTED;
}
