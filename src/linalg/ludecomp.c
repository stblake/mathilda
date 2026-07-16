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
#include "ndlinalg.h"
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
        expr_new_symbol(SYM_Plus), (Expr*[]){a, b}, 2));
}
static Expr* eval_times(Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Times), (Expr*[]){a, b}, 2));
}
static Expr* eval_power(Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Power), (Expr*[]){a, b}, 2));
}
static Expr* eval_together(Expr* a) {
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Together), (Expr*[]){a}, 1));
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
 *  Magnitude-squared as an exact non-negative rational.               *
 *                                                                       *
 *  Recognised numeric forms:                                            *
 *    Integer / BigInt n        -> |n|^2 = n^2                          *
 *    Rational[p, q]            -> p^2 / q^2                            *
 *    Complex[re, im]           -> |re|^2 + |im|^2 if re, im numeric    *
 *                                                                       *
 *  Returns true and writes |e|^2 into out_sq (which the caller must    *
 *  have mpq_init'd already); false for anything we can't handle        *
 *  exactly (symbols, Sqrt, transcendentals, etc.).                     *
 * ------------------------------------------------------------------ */
static bool numeric_abs_sq_as_mpq(Expr* e, mpq_t out_sq) {
    if (!e) return false;

    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) {
        mpz_t n;
        mpz_init(n);
        expr_to_mpz(e, n);
        mpz_t n_sq;
        mpz_init(n_sq);
        mpz_mul(n_sq, n, n);
        mpq_set_z(out_sq, n_sq);
        mpz_clear(n);
        mpz_clear(n_sq);
        return true;
    }

    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;

        if (h == SYM_Rational && e->data.function.arg_count == 2) {
            Expr* num = e->data.function.args[0];
            Expr* den = e->data.function.args[1];
            if (!expr_is_integer_like(num) || !expr_is_integer_like(den)) {
                return false;
            }
            mpz_t p, q;
            mpz_init(p); mpz_init(q);
            expr_to_mpz(num, p);
            expr_to_mpz(den, q);
            mpz_t p_sq, q_sq;
            mpz_init(p_sq); mpz_init(q_sq);
            mpz_mul(p_sq, p, p);
            mpz_mul(q_sq, q, q);
            mpq_set_num(out_sq, p_sq);
            mpq_set_den(out_sq, q_sq);
            mpq_canonicalize(out_sq);
            mpz_clear(p); mpz_clear(q);
            mpz_clear(p_sq); mpz_clear(q_sq);
            return true;
        }

        if (h == SYM_Complex && e->data.function.arg_count == 2) {
            mpq_t re_sq, im_sq;
            mpq_init(re_sq); mpq_init(im_sq);
            bool ok = numeric_abs_sq_as_mpq(e->data.function.args[0], re_sq)
                   && numeric_abs_sq_as_mpq(e->data.function.args[1], im_sq);
            if (ok) {
                mpq_add(out_sq, re_sq, im_sq);
            }
            mpq_clear(re_sq); mpq_clear(im_sq);
            return ok;
        }
    }

    return false;
}

/* True iff every entry in column k of LU, rows [k, rows), is either
 * provably zero or recognised as an exact numeric (Integer / BigInt /
 * Rational / Complex of those).  When true, the pivot rule below
 * upgrades from "first non-zero" to "smallest absolute value" --
 * matching Mathematica's behaviour for exact numeric input. */
static bool lu_column_all_numeric(Expr** LU, int rows, int cols, int k) {
    mpq_t scratch;
    mpq_init(scratch);
    for (int i = k; i < rows; i++) {
        Expr* e = LU[i * cols + k];
        if (is_definitely_zero(e)) continue;
        if (!numeric_abs_sq_as_mpq(e, scratch)) {
            mpq_clear(scratch);
            return false;
        }
    }
    mpq_clear(scratch);
    return true;
}

/* ------------------------------------------------------------------ *
 *  Doolittle core (symbolic).                                          *
 *                                                                       *
 *  Pivoting rule: at step k, scan rows k..rows-1 of column k for the   *
 *  first row whose entry is NOT provably zero.  If none exists, leave  *
 *  the zero in place, mark singular, and continue (the matrix is        *
 *  singular but the factorisation still completes -- matching          *
 *  Mathematica's LUDecomposition::sing behaviour).                      *
 *                                                                       *
 *  Update step (after pivot is in place at LU[k, k]):                  *
 *     for i in [k+1, rows):                                             *
 *         LU[i, k] = LU[i, k] / pivot               (L entry)           *
 *         for j in [k+1, cols):                                         *
 *             LU[i, j] = LU[i, j] - LU[i, k] * LU[k, j]                *
 *                                                                       *
 *  Each arithmetic primitive goes through the Mathilda evaluator so    *
 *  symbolic / rational / Sqrt entries all just work.                    *
 *                                                                       *
 *  Rectangular shape note: the elimination terminates at step          *
 *  min(rows, cols) - 1, after which any remaining rows or columns are  *
 *  already in their final form (extra rows of U are zero by            *
 *  construction; extra columns of U are filled by the Schur update     *
 *  during earlier steps and left alone afterwards).                    *
 * ------------------------------------------------------------------ */
bool lu_symbolic_core(Expr** A_flat, int rows, int cols,
                      Expr*** out_LU_flat, int** out_perm,
                      bool* out_singular)
{
    int steps = (rows < cols) ? rows : cols;

    /* Working LU: deep copy of A so we never alias the caller's
     * entries; we mutate freely in place. */
    size_t total = (size_t)rows * (size_t)cols;
    Expr** LU = (Expr**)malloc(sizeof(Expr*) * total);
    for (size_t t = 0; t < total; t++) LU[t] = expr_copy(A_flat[t]);

    /* perm is the full row permutation -- length `rows`.  Rows past
     * `steps - 1` are never touched by elimination and keep their
     * identity values, matching Mathematica's contract for tall input
     * (e.g. LUDecomposition[{{1,2},{3,4},{5,6}}] returns p = {1,2,3},
     * not {1,2}). */
    int* perm = (int*)malloc(sizeof(int) * (size_t)rows);
    for (int i = 0; i < rows; i++) perm[i] = i + 1;

    bool singular = false;

    for (int k = 0; k < steps; k++) {
        /* Pivot selection.  Two regimes:
         *
         *   (a) Column k of the working LU is entirely exact-numeric
         *       (Integer / BigInt / Rational / Complex of those).
         *       Pick the row with the smallest |entry|^2 among the
         *       non-zeros -- matching Mathematica's exact-numeric
         *       behaviour (probed empirically, e.g.
         *       LUDecomposition[{{1/2, 1/3}, {1/5, 1/7}}] picks the
         *       1/5 pivot, not the 1/2).  Smallest |pivot| tends to
         *       keep intermediate L entries integer for integer
         *       input and avoids unnecessary fraction expansion for
         *       rational input.
         *
         *   (b) Otherwise (any free symbol, Sqrt, transcendental,
         *       or other non-exact-numeric entry in the column):
         *       fall back to "first non-zero" -- matching the spec
         *       example LUDecomposition[{{a, b}, {c, d}}] -> p =
         *       {1, 2}.  A magnitude rule isn't meaningful when
         *       entries can contain free variables.
         *
         *   In both regimes, columns whose entries are all provably
         *   zero from row k down leave LU[k, k] = 0 and flag the
         *   matrix singular. */
        int pivot_row = -1;
        if (lu_column_all_numeric(LU, rows, cols, k)) {
            mpq_t best_sq, mag_sq;
            mpq_init(best_sq);
            mpq_init(mag_sq);
            bool have_best = false;
            for (int i = k; i < rows; i++) {
                if (is_definitely_zero(LU[i * cols + k])) continue;
                bool ok = numeric_abs_sq_as_mpq(LU[i * cols + k], mag_sq);
                (void)ok; /* guaranteed by lu_column_all_numeric */
                if (!have_best || mpq_cmp(mag_sq, best_sq) < 0) {
                    pivot_row = i;
                    mpq_set(best_sq, mag_sq);
                    have_best = true;
                }
            }
            mpq_clear(best_sq);
            mpq_clear(mag_sq);
        } else {
            for (int i = k; i < rows; i++) {
                if (!is_definitely_zero(LU[i * cols + k])) {
                    pivot_row = i;
                    break;
                }
            }
        }
        if (pivot_row < 0) {
            /* Whole column from row k is zero -- singular at this stage.
             * Leave LU[k, k] = 0 and continue. */
            singular = true;
            continue;
        }
        if (pivot_row != k) {
            /* Swap rows pivot_row and k in LU. */
            for (int j = 0; j < cols; j++) {
                Expr* tmp = LU[k * cols + j];
                LU[k * cols + j] = LU[pivot_row * cols + j];
                LU[pivot_row * cols + j] = tmp;
            }
            int tp = perm[k]; perm[k] = perm[pivot_row]; perm[pivot_row] = tp;
        }

        /* Eliminate.  pivot is the upper-triangle diagonal entry; the
         * L entries below it are stored as ratios. */
        Expr* pivot = LU[k * cols + k];        /* borrowed */
        for (int i = k + 1; i < rows; i++) {
            Expr* l_ik = eval_times(expr_copy(LU[i * cols + k]),
                                    eval_power(expr_copy(pivot),
                                               expr_new_integer(-1)));
            expr_free(LU[i * cols + k]);
            LU[i * cols + k] = l_ik;          /* now owns the L entry */

            for (int j = k + 1; j < cols; j++) {
                Expr* prod = eval_times(expr_copy(LU[i * cols + k]),
                                        expr_copy(LU[k * cols + j]));
                Expr* neg  = eval_times(expr_new_integer(-1), prod);
                Expr* new_v = eval_plus(LU[i * cols + j], neg); /* consumes LHS */
                LU[i * cols + j] = new_v;
            }
        }
    }

    *out_LU_flat  = LU;
    *out_perm     = perm;
    *out_singular = singular;
    return true;
}

/* ------------------------------------------------------------------ *
 *  Wrap a flat rows*cols row-major buffer into a List[List[...]].     *
 *  Steals each entry: the caller must not free buf's entries (only    *
 *  the buf array itself).                                              *
 * ------------------------------------------------------------------ */
static Expr* wrap_matrix(Expr** buf, int rows, int cols) {
    Expr** row_exprs = (Expr**)malloc(sizeof(Expr*) * (size_t)rows);
    for (int i = 0; i < rows; i++) {
        Expr** elems = NULL;
        if (cols > 0) elems = (Expr**)malloc(sizeof(Expr*) * (size_t)cols);
        for (int j = 0; j < cols; j++) elems[j] = buf[i * cols + j]; /* steal */
        row_exprs[i] = expr_new_function(
            expr_new_symbol(SYM_List), elems, (size_t)cols);
        if (elems) free(elems);
    }
    Expr* out = expr_new_function(
        expr_new_symbol(SYM_List), row_exprs, (size_t)rows);
    free(row_exprs);
    return out;
}

/* Wrap a 1-indexed length-n int array into a Mathilda List of
 * Integers. */
static Expr* wrap_perm(const int* perm, int n) {
    Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int k = 0; k < n; k++) elems[k] = expr_new_integer(perm[k]);
    Expr* out = expr_new_function(
        expr_new_symbol(SYM_List), elems, (size_t)n);
    free(elems);
    return out;
}

/* Element-wise Together to canonicalise small cancellations in the LU
 * matrix.  Together is Listable so threading is automatic. */
static Expr* tidy_matrix(Expr* m) {
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Together), (Expr*[]){expr_copy(m)}, 1));
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
Expr* lu_symbolic_dispatch(Expr* m, int rows, int cols)
{
    static uint64_t sing_warn_counter = 0;

    size_t total = (size_t)rows * (size_t)cols;

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
    Expr** A_flat = (Expr**)malloc(sizeof(Expr*) * total);
    {
        size_t idx = 0;
        flatten_tensor(matrix_to_use, A_flat, &idx);
    }

    Expr** LU_flat = NULL;
    int*   perm    = NULL;
    bool   singular = false;
    bool ok = lu_symbolic_core(A_flat, rows, cols,
                               &LU_flat, &perm, &singular);

    for (size_t t = 0; t < total; t++) expr_free(A_flat[t]);
    free(A_flat);
    if (m_rat) expr_free(m_rat);

    if (!ok) {
        if (LU_flat) {
            for (size_t t = 0; t < total; t++) expr_free(LU_flat[t]);
            free(LU_flat);
        }
        if (perm) free(perm);
        return NULL;
    }

    if (singular) lu_warn_singular_once(&sing_warn_counter, m);

    Expr* lu = wrap_matrix(LU_flat, rows, cols);
    free(LU_flat);
    Expr* p  = wrap_perm(perm, rows);
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
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), items, 3);
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
Expr* lu_dispatch(Expr* m, int rows, int cols)
{
    CommonInexactInfo info = common_scan_inexact(m);
    if (info.has_inexact) {
        if (info.min_bits <= 53) {
            Expr* fast = lu_machine_dispatch(m, rows, cols);
            if (fast) return fast;
        } else {
            Expr* fast = lu_mpfr_dispatch(m, rows, cols);
            if (fast) return fast;
        }
    }
    return lu_symbolic_dispatch(m, rows, cols);
}

/* ------------------------------------------------------------------ *
 *  Public entry.                                                        *
 *                                                                       *
 *  Accepts any non-empty rectangular matrix (rows >= 1, cols >= 1).    *
 *  Emits LUDecomposition::matsq and returns NULL only for non-list,    *
 *  empty, or higher-rank input -- matching Mathematica's behaviour     *
 *  on rectangular m x n input where it returns a partial Doolittle     *
 *  factorisation with perm of length min(m, n).                         *
 * ------------------------------------------------------------------ */
Expr* builtin_ludecomposition(Expr* res)
{
    if (linalg_call_has_ndarray(res)) return linalg_delist_and_reeval(res);
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return NULL;

    Expr* m = res->data.function.args[0];

    int64_t dims[64];
    int trank = get_tensor_dims(m, dims);
    if (trank != 2 || dims[0] == 0 || dims[1] == 0) {
        char* s = expr_to_string(m);
        fprintf(stderr,
                "LUDecomposition::matsq: Argument %s at position 1 is "
                "not a non-empty rectangular matrix.\n", s);
        free(s);
        return NULL;
    }
    int rows = (int)dims[0];
    int cols = (int)dims[1];

    return lu_dispatch(m, rows, cols);
}

void ludecomp_init(void)
{
    symtab_add_builtin("LUDecomposition", builtin_ludecomposition);
    symtab_get_def("LUDecomposition")->attributes |= ATTR_PROTECTED;
}
