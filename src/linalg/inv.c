/* matinv.c
 *
 * Inverse and PseudoInverse.
 *
 *   Inverse[m]           -- exact / fraction-free Gauss-Jordan inversion
 *                           of a non-empty square matrix.  Lifted verbatim
 *                           from the previous src/linalg.c implementation.
 *
 *   PseudoInverse[m]
 *   PseudoInverse[m,
 *       Tolerance -> t]  -- Moore-Penrose pseudoinverse of a rectangular
 *                           (or rank-deficient square) matrix.
 *
 * Algorithm (PseudoInverse):
 *
 *   For an m x n matrix A with rank r > 0, compute the reduced row-echelon
 *   form R of A.  The first r non-zero rows of R, taken together as a
 *   r x n matrix C, span the row space of A; the columns of A at the pivot
 *   positions, taken together as a m x r matrix B, span the column space.
 *   This gives a full-rank decomposition A = B . C with rank(B) = rank(C) = r.
 *   The Moore-Penrose pseudoinverse is then
 *
 *       A^+ = C^H . (C . C^H)^-1 . (B^H . B)^-1 . B^H
 *
 *   where ^H is the conjugate transpose.  When A is invertible (m == n, r == n)
 *   the formula collapses to the standard inverse.  For the zero matrix
 *   (r == 0) the pseudoinverse is the n x m zero matrix.
 *
 *   For inexact (Real / MPFR) matrices we rationalise the input at the
 *   minimum precision present (the common_rationalize_input pipeline used
 *   throughout the system), do every step in exact rational arithmetic so
 *   the rank is well-defined, then numericalise the final result back to
 *   that precision.  Tolerance -> Automatic uses the input precision.
 *
 * Memory ownership follows the standard builtin contract: this file owns
 * the `res` argument on success and frees it; on failure (returning NULL)
 * the caller (evaluator) retains ownership.  Every intermediate matrix
 * allocated by Dot/Inverse/Transpose/Conjugate via eval_and_free is
 * explicitly released.
 */

#include "inv.h"
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

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ *
 *  Inverse workers.                                                   *
 *                                                                     *
 *  Each worker assumes its caller has already validated the shape of  *
 *  `arg` (rank-2, non-empty, square; `n` == dims[0] == dims[1]).      *
 *  Returns the inverse on success, or NULL on a singular matrix       *
 *  (after printing the Inverse::sing message to stderr).              *
 *                                                                     *
 *  Three methods are exposed via the public dispatcher below:         *
 *                                                                     *
 *    inverse_divfree    Bareiss-like fraction-free Gauss-Jordan on    *
 *                        the augmented matrix [A | I].  The original  *
 *                        implementation, used when                    *
 *                        Method -> Automatic or                       *
 *                        Method -> "DivisionFreeRowReduction".        *
 *                                                                     *
 *    inverse_onestep    Classical Gauss-Jordan on [A | I] with one    *
 *                        division per pivot per row entry.  Each      *
 *                        entry is canonicalised via Together so that  *
 *                        symbolic cancellations are still detected.   *
 *                                                                     *
 *    inverse_cofactor   Adjugate / determinant formula:               *
 *                        A^{-1}[i,j] = (-1)^(i+j) det(M_{j,i})/det(A) *
 *                        with det computed via Laplace cofactor       *
 *                        expansion (the same primitive used by Det).  *
 *                        Time O(n!*n^2); use only for small n.        *
 * ------------------------------------------------------------------ */
static Expr* inverse_divfree(Expr* arg, int n) {
    int cols = 2 * n;

    /* Build augmented matrix [A | I] */
    Expr** matrix = malloc(sizeof(Expr*) * n * cols);
    size_t idx = 0;
    Expr** flat_a = malloc(sizeof(Expr*) * n * n);
    flatten_tensor(arg, flat_a, &idx);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i * cols + j] = flat_a[i * n + j]; /* take ownership */
        }
        for (int j = 0; j < n; j++) {
            matrix[i * cols + n + j] = expr_new_integer(i == j ? 1 : 0);
        }
    }
    free(flat_a); /* elements transferred to matrix */

    /* Fraction-free Gauss-Jordan elimination (Bareiss-like, same as RowReduce) */
    Expr* P = expr_new_integer(1);
    int r = 0;

    for (int c = 0; c < n && r < n; c++) {
        /* Find pivot */
        int pivot_row = -1;
        for (int i = r; i < n; i++) {
            if (!is_zero_poly(matrix[i * cols + c])) {
                pivot_row = i;
                break;
            }
        }
        if (pivot_row == -1) {
            /* Singular matrix */
            char* arg_str = expr_to_string(arg);
            fprintf(stderr, "Inverse::sing: Matrix %s is singular.\n", arg_str);
            free(arg_str);
            expr_free(P);
            for (int i = 0; i < n * cols; i++) expr_free(matrix[i]);
            free(matrix);
            return NULL;
        }

        /* Swap rows if needed */
        if (pivot_row != r) {
            for (int j = 0; j < cols; j++) {
                Expr* tmp = matrix[r * cols + j];
                matrix[r * cols + j] = matrix[pivot_row * cols + j];
                matrix[pivot_row * cols + j] = tmp;
            }
        }

        Expr* pivot = matrix[r * cols + c];

        /* Eliminate column c in all other rows */
        for (int i = 0; i < n; i++) {
            if (i == r) continue;
            Expr* M_ic = matrix[i * cols + c];
            if (is_zero_poly(M_ic)) {
                for (int j = 0; j < cols; j++) {
                    if (j == c) continue;
                    Expr* num_eval = eval_and_free(expr_new_function(
                        expr_new_symbol("Times"),
                        (Expr*[]){expr_copy(pivot), expr_copy(matrix[i * cols + j])}, 2));
                    Expr* new_val = exact_div_wrapper(num_eval, P);
                    expr_free(num_eval);
                    expr_free(matrix[i * cols + j]);
                    matrix[i * cols + j] = new_val;
                }
            } else {
                for (int j = 0; j < cols; j++) {
                    if (j == c) continue;
                    Expr* t1 = eval_and_free(expr_new_function(
                        expr_new_symbol("Times"),
                        (Expr*[]){expr_copy(pivot), expr_copy(matrix[i * cols + j])}, 2));
                    Expr* t2 = eval_and_free(expr_new_function(
                        expr_new_symbol("Times"),
                        (Expr*[]){expr_copy(M_ic), expr_copy(matrix[r * cols + j])}, 2));
                    Expr* t2_neg = eval_and_free(expr_new_function(
                        expr_new_symbol("Times"),
                        (Expr*[]){expr_new_integer(-1), t2}, 2));
                    Expr* num_eval = eval_and_free(expr_new_function(
                        expr_new_symbol("Plus"),
                        (Expr*[]){t1, t2_neg}, 2));
                    Expr* new_val = exact_div_wrapper(num_eval, P);
                    expr_free(num_eval);
                    expr_free(matrix[i * cols + j]);
                    matrix[i * cols + j] = new_val;
                }
            }
            expr_free(matrix[i * cols + c]);
            matrix[i * cols + c] = expr_new_integer(0);
        }

        expr_free(P);
        P = expr_copy(pivot);
        r++;
    }
    expr_free(P);

    /* If we didn't get n pivots, matrix is singular */
    if (r < n) {
        char* arg_str = expr_to_string(arg);
        fprintf(stderr, "Inverse::sing: Matrix %s is singular.\n", arg_str);
        free(arg_str);
        for (int i = 0; i < n * cols; i++) expr_free(matrix[i]);
        free(matrix);
        return NULL;
    }

    /* Normalize: divide each row of the right half by its diagonal element */
    for (int i = 0; i < n; i++) {
        Expr* diag = matrix[i * cols + i];
        for (int j = n; j < cols; j++) {
            if (is_zero_poly(matrix[i * cols + j])) continue;

            Expr* num_val = matrix[i * cols + j];
            Expr* den_val = expr_copy(diag);

            /* Cancel common factors via PolynomialGCD */
            Expr* g = eval_and_free(expr_new_function(
                expr_new_symbol("PolynomialGCD"),
                (Expr*[]){expr_copy(num_val), expr_copy(den_val)}, 2));
            Expr* new_num = exact_div_wrapper(num_val, g);
            Expr* new_den = exact_div_wrapper(den_val, g);
            expr_free(g);

            /* Normalize sign so denominator is positive */
            if (new_den->type == EXPR_INTEGER && new_den->data.integer < 0) {
                Expr* t = eval_and_free(expr_new_function(
                    expr_new_symbol("Times"),
                    (Expr*[]){new_num, expr_new_integer(-1)}, 2));
                new_num = expr_expand(t);
                expr_free(t);
                /* Replace, don't mutate: the integer atom may be shared. */
                int64_t v = -new_den->data.integer;
                expr_free(new_den);
                new_den = expr_new_integer(v);
            }

            if (new_den->type == EXPR_INTEGER && new_den->data.integer == 1) {
                expr_free(new_den);
                matrix[i * cols + j] = new_num;
            } else {
                Expr* inv_den = eval_and_free(expr_new_function(
                    expr_new_symbol("Power"),
                    (Expr*[]){new_den, expr_new_integer(-1)}, 2));
                Expr* final_val = eval_and_free(expr_new_function(
                    expr_new_symbol("Times"),
                    (Expr*[]){new_num, inv_den}, 2));
                matrix[i * cols + j] = expr_expand(final_val);
                expr_free(final_val);
            }
            expr_free(num_val);
            expr_free(den_val);
        }
    }

    /* Extract the right half [I | A^{-1}] -> A^{-1} */
    Expr** rows = malloc(sizeof(Expr*) * n);
    for (int i = 0; i < n; i++) {
        Expr** row_elems = malloc(sizeof(Expr*) * n);
        for (int j = 0; j < n; j++) {
            row_elems[j] = matrix[i * cols + n + j]; /* take ownership */
        }
        rows[i] = expr_new_function(expr_new_symbol("List"), row_elems, n);
        free(row_elems);
    }
    Expr* result = expr_new_function(expr_new_symbol("List"), rows, n);
    free(rows);

    /* Free the left half and the flat matrix array */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            expr_free(matrix[i * cols + j]);
        }
        /* right half elements were transferred to result */
    }
    free(matrix);

    return result;
}

/* ------------------------------------------------------------------ *
 *  Inverse -- OneStepRowReduction.                                    *
 *                                                                     *
 *  Classical Gauss-Jordan on the augmented matrix [A | I]: per pivot, *
 *  divide the pivot row through by the pivot value (so the pivot      *
 *  becomes 1) and then subtract a scaled copy of the pivot row from   *
 *  every other row to zero column c.  Reads the right half off as     *
 *  the inverse with no post-normalisation pass.                       *
 *                                                                     *
 *  Per-entry canonicalisation through matsol_canon_entry (Together)   *
 *  ensures that symbolic cancellations are still detected by          *
 *  is_zero_poly, so a singular symbolic matrix produces a definite    *
 *  Inverse::sing rather than a silently wrong answer.                 *
 * ------------------------------------------------------------------ */
static Expr* inverse_onestep(Expr* arg, int n) {
    int cols = 2 * n;

    Expr** matrix = malloc(sizeof(Expr*) * (size_t)n * (size_t)cols);
    Expr** flat_a = malloc(sizeof(Expr*) * (size_t)n * (size_t)n);
    {
        size_t idx = 0;
        flatten_tensor(arg, flat_a, &idx);
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i * cols + j] = flat_a[i * n + j]; /* take ownership */
        }
        for (int j = 0; j < n; j++) {
            matrix[i * cols + n + j] = expr_new_integer(i == j ? 1 : 0);
        }
    }
    free(flat_a);

    for (int c = 0; c < n; c++) {
        int pivot_row = -1;
        for (int i = c; i < n; i++) {
            if (!is_zero_poly(matrix[i * cols + c])) {
                pivot_row = i;
                break;
            }
        }
        if (pivot_row == -1) {
            char* arg_str = expr_to_string(arg);
            fprintf(stderr, "Inverse::sing: Matrix %s is singular.\n", arg_str);
            free(arg_str);
            for (int i = 0; i < n * cols; i++) expr_free(matrix[i]);
            free(matrix);
            return NULL;
        }

        if (pivot_row != c) {
            for (int j = 0; j < cols; j++) {
                Expr* tmp = matrix[c * cols + j];
                matrix[c * cols + j] = matrix[pivot_row * cols + j];
                matrix[pivot_row * cols + j] = tmp;
            }
        }

        /* Normalise the pivot row: pivot cell -> 1, others -> entry/pivot. */
        Expr* pivot = expr_copy(matrix[c * cols + c]);
        for (int j = 0; j < cols; j++) {
            if (j == c) {
                expr_free(matrix[c * cols + j]);
                matrix[c * cols + j] = expr_new_integer(1);
            } else if (!is_zero_poly(matrix[c * cols + j])) {
                Expr* old = matrix[c * cols + j];
                matrix[c * cols + j] = matsol_div_entry(old, pivot);
                expr_free(old);
            }
        }
        expr_free(pivot);

        /* Eliminate column c in every other row. */
        for (int i = 0; i < n; i++) {
            if (i == c) continue;
            Expr* factor = matrix[i * cols + c];
            if (is_zero_poly(factor)) continue;

            for (int j = 0; j < cols; j++) {
                if (j == c) continue;
                Expr* term = eval_and_free(expr_new_function(
                    expr_new_symbol("Times"),
                    (Expr*[]){expr_copy(factor), expr_copy(matrix[c * cols + j])}, 2));
                Expr* neg_term = eval_and_free(expr_new_function(
                    expr_new_symbol("Times"),
                    (Expr*[]){expr_new_integer(-1), term}, 2));
                Expr* updated = eval_and_free(expr_new_function(
                    expr_new_symbol("Plus"),
                    (Expr*[]){expr_copy(matrix[i * cols + j]), neg_term}, 2));
                Expr* updated_c = matsol_canon_entry(updated);
                expr_free(updated);
                expr_free(matrix[i * cols + j]);
                matrix[i * cols + j] = updated_c;
            }
            expr_free(matrix[i * cols + c]);
            matrix[i * cols + c] = expr_new_integer(0);
        }
    }

    /* Extract the right half [I | A^{-1}] -> A^{-1} */
    Expr** rows = malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) {
        Expr** row_elems = malloc(sizeof(Expr*) * (size_t)n);
        for (int j = 0; j < n; j++) {
            row_elems[j] = matrix[i * cols + n + j]; /* take ownership */
        }
        rows[i] = expr_new_function(expr_new_symbol("List"), row_elems, (size_t)n);
        free(row_elems);
    }
    Expr* result = expr_new_function(expr_new_symbol("List"), rows, (size_t)n);
    free(rows);

    /* Free the left (identity) half; right-half ownership transferred. */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            expr_free(matrix[i * cols + j]);
        }
    }
    free(matrix);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Inverse -- CofactorExpansion.                                      *
 *                                                                     *
 *  Computes the inverse via the classical adjugate / determinant      *
 *  formula                                                            *
 *                                                                     *
 *      A^{-1}[i, j] = (-1)^(i + j) det(M_{j, i}) / det(A)             *
 *                                                                     *
 *  where M_{j, i} is the (n-1) x (n-1) minor obtained by deleting     *
 *  row j and column i from A.  det is computed via Laplace cofactor   *
 *  expansion (linalg.c:laplace_det), the same primitive used by Det   *
 *  and by the LinearSolve / RowReduce "CofactorExpansion" workers.    *
 *                                                                     *
 *  Time complexity is O(n! * n^2) -- intended for small n and for     *
 *  the closed-form symbolic inverse of small matrices.                *
 * ------------------------------------------------------------------ */
static Expr* inverse_cofactor(Expr* arg, int n) {
    Expr** flat = malloc(sizeof(Expr*) * (size_t)n * (size_t)n);
    {
        size_t idx = 0;
        flatten_tensor(arg, flat, &idx);
    }

    int* cols_idx = malloc(sizeof(int) * (size_t)n);
    for (int i = 0; i < n; i++) cols_idx[i] = i;
    Expr* det_a = laplace_det(flat, n, n, 0, cols_idx);
    free(cols_idx);

    if (is_zero_poly(det_a)) {
        char* arg_str = expr_to_string(arg);
        fprintf(stderr, "Inverse::sing: Matrix %s is singular.\n", arg_str);
        free(arg_str);
        expr_free(det_a);
        for (int i = 0; i < n * n; i++) expr_free(flat[i]);
        free(flat);
        return NULL;
    }

    /* 1x1 fast path: A^{-1} = {{1/A[0,0]}}.  The minor is 0 x 0 and
     * laplace_det's recursion base is n == 1, so we can't go through
     * the generic loop below. */
    if (n == 1) {
        Expr* one = expr_new_integer(1);
        Expr* inv = matsol_div_entry(one, det_a);
        expr_free(one);
        Expr** row_elems = malloc(sizeof(Expr*));
        row_elems[0] = inv;
        Expr** rows = malloc(sizeof(Expr*));
        rows[0] = expr_new_function(expr_new_symbol("List"), row_elems, 1);
        free(row_elems);
        Expr* result = expr_new_function(expr_new_symbol("List"), rows, 1);
        free(rows);
        expr_free(det_a);
        for (int i = 0; i < n * n; i++) expr_free(flat[i]);
        free(flat);
        return result;
    }

    /* Build A^{-1} row by row.  result[i][j] = cof(j, i) / det(A). */
    int m = n - 1;
    Expr** minor_flat = malloc(sizeof(Expr*) * (size_t)m * (size_t)m);
    int* minor_cols = malloc(sizeof(int) * (size_t)m);
    for (int j = 0; j < m; j++) minor_cols[j] = j;

    Expr** rows = malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) {
        Expr** row_elems = malloc(sizeof(Expr*) * (size_t)n);
        for (int j = 0; j < n; j++) {
            /* Build minor M_{j, i}: delete row j and column i from A. */
            int dst = 0;
            for (int r = 0; r < n; r++) {
                if (r == j) continue;
                int dst_col = 0;
                for (int c = 0; c < n; c++) {
                    if (c == i) continue;
                    minor_flat[dst * m + dst_col] = expr_copy(flat[r * n + c]);
                    dst_col++;
                }
                dst++;
            }

            Expr* minor_det = laplace_det(minor_flat, m, m, 0, minor_cols);

            /* Sign: (-1)^(i + j). */
            if (((i + j) & 1) != 0) {
                Expr* neg = eval_and_free(expr_new_function(
                    expr_new_symbol("Times"),
                    (Expr*[]){expr_new_integer(-1), minor_det}, 2));
                minor_det = neg;
            }

            Expr* entry = matsol_div_entry(minor_det, det_a);
            expr_free(minor_det);

            /* Free minor_flat entries for next iteration. */
            for (int k = 0; k < m * m; k++) expr_free(minor_flat[k]);

            row_elems[j] = entry;
        }
        rows[i] = expr_new_function(expr_new_symbol("List"), row_elems, (size_t)n);
        free(row_elems);
    }

    Expr* result = expr_new_function(expr_new_symbol("List"), rows, (size_t)n);
    free(rows);
    free(minor_cols);
    free(minor_flat);
    expr_free(det_a);
    for (int i = 0; i < n * n; i++) expr_free(flat[i]);
    free(flat);

    return result;
}

/* ------------------------------------------------------------------ *
 *  Inverse -- public dispatcher.                                      *
 *                                                                     *
 *      Inverse[m]                                                     *
 *      Inverse[m, Method -> "<name>"]                                 *
 *                                                                     *
 *  Method names match the RowReduce / LinearSolve set so that any     *
 *  matrix routine can be driven through the same option.  Default is  *
 *  Method -> Automatic, which currently resolves to                   *
 *  DivisionFreeRowReduction (the historical behaviour).  An unknown   *
 *  Method value emits Inverse::method and leaves the call             *
 *  unevaluated.  Non-square or empty input emits Inverse::matsq.      *
 * ------------------------------------------------------------------ */
Expr* builtin_inverse(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    Expr* arg = res->data.function.args[0];

    MatsolMethod method = MATSOL_AUTOMATIC;
    if (argc == 2) {
        method = matsol_parse_method_option(res->data.function.args[1]);
        if (method == MATSOL_INVALID) {
            static uint64_t last_warned = 0;
            matsol_warn_once(&last_warned, res,
                "Inverse::method: Method option value is not one of "
                "\"Automatic\", \"DivisionFreeRowReduction\", "
                "\"OneStepRowReduction\", \"CofactorExpansion\".\n");
            return NULL;
        }
    }

    int64_t dims[64];
    int trank = get_tensor_dims(arg, dims);
    if (trank != 2 || dims[0] != dims[1] || dims[0] == 0) {
        char* arg_str = expr_to_string(arg);
        fprintf(stderr,
                "Inverse::matsq: Argument %s at position 1 is not a "
                "non-empty square matrix.\n",
                arg_str);
        free(arg_str);
        return NULL;
    }

    int n = (int)dims[0];
    switch (method) {
        case MATSOL_AUTOMATIC:
        case MATSOL_DIVFREE:    return inverse_divfree(arg, n);
        case MATSOL_ONESTEP:    return inverse_onestep(arg, n);
        case MATSOL_COFACTOR:   return inverse_cofactor(arg, n);
        case MATSOL_INVALID:    return NULL;  /* unreachable */
    }
    return NULL;
}

/* ------------------------------------------------------------------ *
 *  PseudoInverse helpers                                              *
 * ------------------------------------------------------------------ */

/* Return Conjugate[Transpose[m]] as a fresh expression.  Caller owns. */
static Expr* hermitian_transpose(Expr* m) {
    Expr* trans = eval_and_free(expr_new_function(
        expr_new_symbol("Transpose"),
        (Expr*[]){expr_copy(m)}, 1));
    Expr* conj = eval_and_free(expr_new_function(
        expr_new_symbol("Conjugate"),
        (Expr*[]){trans}, 1));
    return conj;
}

/* Matrix product a . b via Dot, fully evaluated.  Caller owns. */
static Expr* mat_mult(Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(
        expr_new_symbol("Dot"),
        (Expr*[]){expr_copy(a), expr_copy(b)}, 2));
}

/* Inverse[m] via the registered builtin, fully evaluated.  Caller owns.
 * Returns NULL (with an Inverse::sing error already on stderr) when m is
 * singular; otherwise a fresh n x n matrix. */
static Expr* mat_invert(Expr* m) {
    Expr* call = expr_new_function(expr_new_symbol("Inverse"),
        (Expr*[]){expr_copy(m)}, 1);
    Expr* result = evaluate(call);
    /* If Inverse failed, the evaluator returns the call unchanged. */
    if (result && result->type == EXPR_FUNCTION
        && result->data.function.head->type == EXPR_SYMBOL
        && result->data.function.head->data.symbol == SYM_Inverse) {
        expr_free(call);
        expr_free(result);
        return NULL;
    }
    expr_free(call);
    return result;
}

/* RowReduce[m] via the registered builtin, fully evaluated.  Caller owns. */
static Expr* mat_rref(Expr* m) {
    return eval_and_free(expr_new_function(
        expr_new_symbol("RowReduce"),
        (Expr*[]){expr_copy(m)}, 1));
}

/* Locate pivot columns of a row-reduced matrix.  Writes the column index
 * of each pivot into pivot_cols (must be sized >= min(m, n)) and stores
 * the rank.  A row is considered a pivot row when its first non-zero entry
 * is encountered; subsequent rows whose entries up to that point are all
 * zero are recognised by scanning from the left. */
static void find_pivots(Expr* rref, int m, int n, int* rank_out, int* pivot_cols) {
    int r = 0;
    int next_search_col = 0;
    for (int i = 0; i < m; i++) {
        Expr* row = rref->data.function.args[i];
        int found = -1;
        for (int j = next_search_col; j < n; j++) {
            if (!is_zero_poly(row->data.function.args[j])) {
                found = j;
                break;
            }
        }
        if (found < 0) continue;
        pivot_cols[r++] = found;
        next_search_col = found + 1;
    }
    *rank_out = r;
}

/* Build a fresh m x r matrix by picking the columns of `a` (m x n) named
 * in col_indices.  Owns the returned expression; `a` is untouched. */
static Expr* extract_columns(Expr* a, int m, int r, const int* col_indices) {
    Expr** rows = malloc(sizeof(Expr*) * (size_t)m);
    for (int i = 0; i < m; i++) {
        Expr* a_row = a->data.function.args[i];
        Expr** elems = NULL;
        if (r > 0) elems = malloc(sizeof(Expr*) * (size_t)r);
        for (int k = 0; k < r; k++) {
            elems[k] = expr_copy(a_row->data.function.args[col_indices[k]]);
        }
        rows[i] = expr_new_function(expr_new_symbol("List"), elems, (size_t)r);
        if (elems) free(elems);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), rows, (size_t)m);
    free(rows);
    return out;
}

/* Build a fresh r x n matrix containing the first r rows of `rref`.
 * Caller owns the result. */
static Expr* extract_rows(Expr* rref, int r, int n) {
    Expr** rows = malloc(sizeof(Expr*) * (size_t)r);
    for (int i = 0; i < r; i++) {
        Expr* src_row = rref->data.function.args[i];
        Expr** elems = NULL;
        if (n > 0) elems = malloc(sizeof(Expr*) * (size_t)n);
        for (int j = 0; j < n; j++) {
            elems[j] = expr_copy(src_row->data.function.args[j]);
        }
        rows[i] = expr_new_function(expr_new_symbol("List"), elems, (size_t)n);
        if (elems) free(elems);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), rows, (size_t)r);
    free(rows);
    return out;
}

/* n x m zero matrix.  Caller owns. */
static Expr* zero_matrix(int rows, int cols) {
    Expr** row_arr = malloc(sizeof(Expr*) * (size_t)rows);
    for (int i = 0; i < rows; i++) {
        Expr** elems = NULL;
        if (cols > 0) elems = malloc(sizeof(Expr*) * (size_t)cols);
        for (int j = 0; j < cols; j++) elems[j] = expr_new_integer(0);
        row_arr[i] = expr_new_function(expr_new_symbol("List"), elems, (size_t)cols);
        if (elems) free(elems);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), row_arr, (size_t)rows);
    free(row_arr);
    return out;
}

/* Apply Expand element-wise (Listable threads it) to keep results in a
 * tidy normal form. */
static Expr* expand_matrix(Expr* m) {
    return eval_and_free(expr_new_function(
        expr_new_symbol("Expand"),
        (Expr*[]){expr_copy(m)}, 1));
}

/* ------------------------------------------------------------------ *
 *  Tolerance-option parsing                                           *
 *  Accepted forms:                                                    *
 *      Tolerance -> Automatic                                         *
 *      Tolerance -> <non-negative number>                             *
 *  Returns true on a recognised option, false otherwise.              *
 * ------------------------------------------------------------------ */
static bool parse_tolerance_option(Expr* opt, bool* is_automatic_out,
                                   double* tol_out) {
    *is_automatic_out = true;
    *tol_out = 0.0;
    if (opt->type != EXPR_FUNCTION) return false;
    if (opt->data.function.head->type != EXPR_SYMBOL) return false;
    const char* hd = opt->data.function.head->data.symbol;
    if ((hd != SYM_Rule && hd != SYM_RuleDelayed) ||
        opt->data.function.arg_count != 2) return false;
    Expr* lhs = opt->data.function.args[0];
    Expr* rhs = opt->data.function.args[1];
    if (lhs->type != EXPR_SYMBOL) return false;
    if (lhs->data.symbol != SYM_Tolerance) return false;

    if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic) {
        *is_automatic_out = true;
        return true;
    }
    if (rhs->type == EXPR_INTEGER) {
        if (rhs->data.integer < 0) return false;
        *is_automatic_out = false;
        *tol_out = (double)rhs->data.integer;
        return true;
    }
    if (rhs->type == EXPR_REAL) {
        if (rhs->data.real < 0.0) return false;
        *is_automatic_out = false;
        *tol_out = rhs->data.real;
        return true;
    }
    if (rhs->type == EXPR_FUNCTION
        && rhs->data.function.head->type == EXPR_SYMBOL
        && rhs->data.function.head->data.symbol == SYM_Rational
        && rhs->data.function.arg_count == 2) {
        Expr* num = rhs->data.function.args[0];
        Expr* den = rhs->data.function.args[1];
        if (num->type != EXPR_INTEGER || den->type != EXPR_INTEGER) return false;
        if (num->data.integer < 0) return false;
        if (den->data.integer <= 0) return false;
        *is_automatic_out = false;
        *tol_out = (double)num->data.integer / (double)den->data.integer;
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  PseudoInverse -- exact pipeline                                    *
 * ------------------------------------------------------------------ */
static Expr* pseudoinverse_exact(Expr* a, int m, int n) {
    /* Row reduce to identify rank and pivot columns. */
    Expr* rref = mat_rref(a);
    if (!rref) return NULL;

    int* pivot_cols = malloc(sizeof(int) * (size_t)(m < n ? m : n));
    int rank = 0;
    find_pivots(rref, m, n, &rank, pivot_cols);

    /* Zero matrix -> n x m zero matrix. */
    if (rank == 0) {
        Expr* z = zero_matrix(n, m);
        expr_free(rref);
        free(pivot_cols);
        return z;
    }

    /* Full-rank decomposition A = B . C with B m x r and C r x n. */
    Expr* B = extract_columns(a, m, rank, pivot_cols);
    Expr* C = extract_rows(rref, rank, n);
    expr_free(rref);
    free(pivot_cols);

    /* A^+ = C^H . (C C^H)^-1 . (B^H B)^-1 . B^H */
    Expr* BH = hermitian_transpose(B);
    Expr* CH = hermitian_transpose(C);
    Expr* BHB = mat_mult(BH, B);  /* r x r */
    Expr* CCH = mat_mult(C, CH);  /* r x r */
    expr_free(B);
    expr_free(C);

    Expr* BHB_inv = mat_invert(BHB);
    Expr* CCH_inv = mat_invert(CCH);
    expr_free(BHB);
    expr_free(CCH);
    if (!BHB_inv || !CCH_inv) {
        if (BHB_inv) expr_free(BHB_inv);
        if (CCH_inv) expr_free(CCH_inv);
        expr_free(BH);
        expr_free(CH);
        return NULL;
    }

    Expr* mid = mat_mult(CCH_inv, BHB_inv);  /* r x r */
    expr_free(BHB_inv);
    expr_free(CCH_inv);
    Expr* left = mat_mult(CH, mid);          /* n x r */
    expr_free(CH);
    expr_free(mid);
    Expr* result = mat_mult(left, BH);       /* n x m */
    expr_free(left);
    expr_free(BH);

    Expr* simplified = expand_matrix(result);
    expr_free(result);
    return simplified;
}

/* ------------------------------------------------------------------ *
 *  PseudoInverse builtin                                              *
 * ------------------------------------------------------------------ */
Expr* builtin_pseudoinverse(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    Expr* arg = res->data.function.args[0];

    /* Parse Tolerance option(s).  Unknown options leave the call
     * unevaluated so the user sees a clean PseudoInverse[...] back. */
    bool tol_automatic = true;
    double tol_value = 0.0;
    for (size_t i = 1; i < argc; i++) {
        bool is_auto = true;
        double tv = 0.0;
        if (!parse_tolerance_option(res->data.function.args[i],
                                    &is_auto, &tv)) {
            return NULL;
        }
        if (!is_auto) {
            tol_automatic = false;
            tol_value = tv;
        }
    }
    (void)tol_automatic;
    (void)tol_value;

    /* Validate matrix shape (rank-2 tensor, non-empty). */
    int64_t dims[64];
    int trank = get_tensor_dims(arg, dims);
    if (trank != 2 || dims[0] == 0 || dims[1] == 0) {
        char* arg_str = expr_to_string(arg);
        fprintf(stderr,
                "PseudoInverse::matrix: Argument %s at position 1 is not a "
                "non-empty rectangular matrix.\n",
                arg_str);
        free(arg_str);
        return NULL;
    }
    int m = (int)dims[0];
    int n = (int)dims[1];

    /* Inexact preprocessing: rationalise -> exact pipeline -> numericalise.
     * This mirrors the Eigenvalues / Eigenvectors / Solve flow so that
     * approximate matrices yield approximate answers at the input precision
     * while still benefiting from a well-defined rank in the row reduction. */
    CommonInexactInfo info = common_scan_inexact(arg);
    Expr* a_rat = NULL;
    long prec_bits = 53;
    Expr* matrix_to_use = arg;
    if (info.has_inexact) {
        prec_bits = info.min_bits ? info.min_bits : 53;
        a_rat = common_rationalize_input(arg, prec_bits);
        matrix_to_use = a_rat;
    }

    Expr* result = pseudoinverse_exact(matrix_to_use, m, n);

    if (a_rat) expr_free(a_rat);

    if (!result) {
        /* The exact pipeline must always succeed for a well-formed matrix
         * because the rank-r mid-products are constructed to be invertible
         * by definition.  If we got NULL anyway something genuinely odd
         * happened; let the call stay unevaluated. */
        return NULL;
    }

    if (info.has_inexact) {
        Expr* numeric = common_numericalize_result(result, prec_bits);
        expr_free(result);
        result = numeric;
    }

    return result;
}

void matinv_init(void) {
    symtab_add_builtin("Inverse", builtin_inverse);
    symtab_get_def("Inverse")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("PseudoInverse", builtin_pseudoinverse);
    symtab_get_def("PseudoInverse")->attributes |= ATTR_PROTECTED;
}
