/* matsol.c
 *
 * Method-aware dispatcher for `RowReduce[m]` and `LinearSolve[m, b]`.
 *
 * Both builtins accept an optional `Method -> "<name>"` argument
 * (RowReduce as arg 2, LinearSolve as arg 3) and route to one of
 * three explicit algorithms:
 *
 *   "Automatic"                  -- alias for "DivisionFreeRowReduction"
 *   "DivisionFreeRowReduction"   -- Bareiss-like fraction-free Gauss-Jordan
 *   "OneStepRowReduction"        -- classical Gauss-Jordan, one division per
 *                                   pivot per element of the pivot row
 *   "CofactorExpansion"          -- LinearSolve: Cramer's rule via Laplace
 *                                   cofactor expansion;
 *                                   RowReduce: identity-if-invertible (for a
 *                                   non-singular square matrix) with fallback
 *                                   to DivisionFreeRowReduction on
 *                                   singular / rectangular / empty input.
 *
 * The DivisionFreeRowReduction workers are direct lifts of the
 * previous bodies of `builtin_rowreduce` and `builtin_linearsolve` in
 * src/linalg.c, so existing behaviour is preserved bit-for-bit when
 * the user does not supply a Method option.
 *
 * Algorithm choice for the four supported matrix families:
 *   - Machine precision (Real)            -> OneStep is fastest; DivFree fine.
 *   - Bignum integer                      -> DivFree avoids GCD blow-up.
 *   - MPFR / arbitrary precision          -> Same as machine precision.
 *   - Symbolic                            -> DivFree avoids algebraic growth;
 *                                            CofactorExpansion gives the
 *                                            textbook Cramer form for small n.
 */

#include "linsolve.h"
#include "linalg.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"
#include "poly.h"
#include "expand.h"
#include "sym_names.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------
 * Method-option parsing.  Mirrors the canonical SYM_Method / SYM_Rule
 * idiom (see src/integrate.c:199 and src/list.c:1480-1491).  Exposed
 * via matsol.h so src/matinv.c (Inverse) can reuse the same grammar.
 * ------------------------------------------------------------------ */
MatsolMethod matsol_parse_method_option(Expr* opt) {
    if (opt->type != EXPR_FUNCTION) return MATSOL_INVALID;
    if (opt->data.function.head->type != EXPR_SYMBOL) return MATSOL_INVALID;
    const char* hd = opt->data.function.head->data.symbol;
    if ((hd != SYM_Rule && hd != SYM_RuleDelayed) ||
        opt->data.function.arg_count != 2) return MATSOL_INVALID;
    Expr* lhs = opt->data.function.args[0];
    Expr* rhs = opt->data.function.args[1];
    if (lhs->type != EXPR_SYMBOL || lhs->data.symbol != SYM_Method)
        return MATSOL_INVALID;

    /* Accept Method -> Automatic (the symbol). */
    if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic)
        return MATSOL_AUTOMATIC;
    if (rhs->type != EXPR_STRING) return MATSOL_INVALID;
    if (strcmp(rhs->data.string, "Automatic")                == 0) return MATSOL_AUTOMATIC;
    if (strcmp(rhs->data.string, "DivisionFreeRowReduction") == 0) return MATSOL_DIVFREE;
    if (strcmp(rhs->data.string, "OneStepRowReduction")      == 0) return MATSOL_ONESTEP;
    if (strcmp(rhs->data.string, "CofactorExpansion")        == 0) return MATSOL_COFACTOR;
    return MATSOL_INVALID;
}

/* Rate-limit a per-call warning so test loops don't spew. */
void matsol_warn_once(uint64_t* last_hash, Expr* key, const char* msg) {
    uint64_t h = expr_hash(key);
    if (h == *last_hash) return;
    *last_hash = h;
    fputs(msg, stderr);
}

/* ------------------------------------------------------------------
 * RowReduce -- DivisionFreeRowReduction (Bareiss-like)
 *
 * Direct lift of the previous body of `builtin_rowreduce` from
 * src/linalg.c.  The running pivot P is freed before being replaced
 * by a fresh expr_copy of the current pivot cell; do not change the
 * order of those two operations -- the pivot cell is *not* zeroed
 * until later in the same loop iteration.
 * ------------------------------------------------------------------ */
static Expr* rowreduce_divfree(Expr* arg) {
    int64_t dims[64];
    int rank = get_tensor_dims(arg, dims);
    if (rank != 2 || dims[0] == 0 || dims[1] == 0) {
        return expr_copy(arg);
    }

    int m = (int)dims[0];
    int n = (int)dims[1];

    Expr** matrix = malloc(sizeof(Expr*) * (size_t)m * (size_t)n);
    size_t idx = 0;
    flatten_tensor(arg, matrix, &idx);

    Expr* P = expr_new_integer(1);
    int r = 0;

    for (int c = 0; c < n && r < m; c++) {
        int pivot_row = -1;
        for (int i = r; i < m; i++) {
            if (!is_zero_poly(matrix[i * n + c])) {
                pivot_row = i;
                break;
            }
        }
        if (pivot_row == -1) continue;

        if (pivot_row != r) {
            for (int j = 0; j < n; j++) {
                Expr* tmp = matrix[r * n + j];
                matrix[r * n + j] = matrix[pivot_row * n + j];
                matrix[pivot_row * n + j] = tmp;
            }
        }

        Expr* pivot = matrix[r * n + c];

        for (int i = 0; i < m; i++) {
            if (i == r) continue;
            Expr* M_ic = matrix[i * n + c];
            if (is_zero_poly(M_ic)) {
                for (int j = 0; j < n; j++) {
                    if (j == c) continue;
                    Expr* num_eval = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_copy(pivot), expr_copy(matrix[i * n + j])}, 2));
                    Expr* new_val = exact_div_wrapper(num_eval, P);
                    expr_free(num_eval);
                    expr_free(matrix[i * n + j]);
                    matrix[i * n + j] = new_val;
                }
            } else {
                for (int j = 0; j < n; j++) {
                    if (j == c) continue;
                    Expr* t1 = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_copy(pivot), expr_copy(matrix[i * n + j])}, 2));
                    Expr* t2 = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_copy(M_ic), expr_copy(matrix[r * n + j])}, 2));
                    Expr* t2_neg = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), t2}, 2));
                    Expr* num_eval = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){t1, t2_neg}, 2));
                    Expr* new_val = exact_div_wrapper(num_eval, P);
                    expr_free(num_eval);
                    expr_free(matrix[i * n + j]);
                    matrix[i * n + j] = new_val;
                }
            }
            expr_free(matrix[i * n + c]);
            matrix[i * n + c] = expr_new_integer(0);
        }

        expr_free(P);
        P = expr_copy(pivot);
        r++;
    }
    expr_free(P);

    /* Final pass: scale each pivot row so its leading entry is 1, in
     * the same minimally-canonical way the previous code used. */
    for (int i = 0; i < m; i++) {
        Expr* leading = NULL;
        int lead_j = -1;
        for (int j = 0; j < n; j++) {
            if (!is_zero_poly(matrix[i * n + j])) {
                leading = expr_copy(matrix[i * n + j]);
                lead_j = j;
                break;
            }
        }
        if (leading) {
            for (int j = lead_j; j < n; j++) {
                if (j == lead_j) {
                    expr_free(matrix[i * n + j]);
                    matrix[i * n + j] = expr_new_integer(1);
                } else if (!is_zero_poly(matrix[i * n + j])) {
                    Expr* num_val = matrix[i * n + j];
                    Expr* den_val = expr_copy(leading);

                    Expr* g = eval_and_free(expr_new_function(expr_new_symbol("PolynomialGCD"), (Expr*[]){expr_copy(num_val), expr_copy(den_val)}, 2));
                    Expr* new_num = exact_div_wrapper(num_val, g);
                    Expr* new_den = exact_div_wrapper(den_val, g);
                    expr_free(g);

                    if (new_den->type == EXPR_INTEGER && new_den->data.integer < 0) {
                        Expr* t = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){new_num, expr_new_integer(-1)}, 2));
                        new_num = expr_expand(t);
                        expr_free(t);
                        /* Replace, don't mutate: the integer atom may be
                         * shared (M3 atom-sharing). */
                        int64_t v = -new_den->data.integer;
                        expr_free(new_den);
                        new_den = expr_new_integer(v);
                    }

                    if (new_den->type == EXPR_INTEGER && new_den->data.integer == 1) {
                        expr_free(new_den);
                        matrix[i * n + j] = new_num;
                    } else {
                        Expr* inv_den = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){new_den, expr_new_integer(-1)}, 2));
                        Expr* final_val = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){new_num, inv_den}, 2));
                        matrix[i * n + j] = expr_expand(final_val);
                        expr_free(final_val);
                    }
                    expr_free(num_val);
                    expr_free(den_val);
                }
            }
            expr_free(leading);
        }
    }

    Expr** rows = malloc(sizeof(Expr*) * (size_t)m);
    for (int i = 0; i < m; i++) {
        Expr** row_elems = malloc(sizeof(Expr*) * (size_t)n);
        for (int j = 0; j < n; j++) {
            row_elems[j] = matrix[i * n + j];
        }
        rows[i] = expr_new_function(expr_new_symbol("List"), row_elems, (size_t)n);
        free(row_elems);
    }

    Expr* result = expr_new_function(expr_new_symbol("List"), rows, (size_t)m);
    free(rows);
    free(matrix);

    return result;
}

/* ------------------------------------------------------------------
 * Shared helpers used by the OneStep workers.
 *
 * `matsol_canon_entry(e)` -- canonicalise a matrix entry so subsequent
 * is_zero_poly checks behave correctly even on symbolic rationals.
 * For numeric input (Integer / Real / Bignum / Rational) Together is
 * a no-op cheaply; for symbolic input it combines everything over a
 * common denominator so e.g. `c + f - c - d c/a - f + d c/a` reduces
 * to 0 rather than persisting as a non-trivial sum that the
 * structural zero-check would miss.  Without this, OneStep on a
 * singular *symbolic* matrix mis-pivots and yields a wrong RREF.
 *
 * `matsol_div_entry(num, den)` -- num / den, canonicalised.  Prefers exact
 * polynomial division via exact_div_wrapper when applicable.
 * ------------------------------------------------------------------ */
Expr* matsol_canon_entry(Expr* e) {
    return eval_and_free(expr_new_function(
        expr_new_symbol("Together"), (Expr*[]){expr_copy(e)}, 1));
}

Expr* matsol_div_entry(Expr* num, Expr* den) {
    Expr* q = exact_div_wrapper(num, den);
    if (!q) {
        Expr* inv = eval_and_free(expr_new_function(
            expr_new_symbol("Power"),
            (Expr*[]){expr_copy(den), expr_new_integer(-1)}, 2));
        q = eval_and_free(expr_new_function(
            expr_new_symbol("Times"),
            (Expr*[]){expr_copy(num), inv}, 2));
    }
    Expr* qc = matsol_canon_entry(q);
    expr_free(q);
    return qc;
}

/* ------------------------------------------------------------------
 * RowReduce -- OneStepRowReduction (textbook Gauss-Jordan with division)
 *
 * Per pivot:
 *   1. Find a non-zero pivot in column c at row r' >= r; swap.
 *   2. Divide pivot row by pivot, in place (one division per entry).
 *      The pivot cell becomes 1 (so we skip dividing it -- store 1
 *      directly).
 *   3. For every other row i, subtract M[i,c] * pivot_row from row i,
 *      writing the result back.  After elimination M[i,c] is 0.
 *
 * Result is already RREF; no final scaling pass.
 *
 * Note: pivot detection still uses `is_zero_poly`, which is purely
 * structural.  For Real-input matrices a floating-cancellation pivot
 * of magnitude 1e-18 will not be detected as zero -- same property
 * as the DivisionFreeRowReduction path.
 * ------------------------------------------------------------------ */
static Expr* rowreduce_onestep(Expr* arg) {
    int64_t dims[64];
    int rank = get_tensor_dims(arg, dims);
    if (rank != 2 || dims[0] == 0 || dims[1] == 0) {
        return expr_copy(arg);
    }

    int m = (int)dims[0];
    int n = (int)dims[1];

    Expr** matrix = malloc(sizeof(Expr*) * (size_t)m * (size_t)n);
    size_t idx = 0;
    flatten_tensor(arg, matrix, &idx);

    int r = 0;
    for (int c = 0; c < n && r < m; c++) {
        int pivot_row = -1;
        for (int i = r; i < m; i++) {
            if (!is_zero_poly(matrix[i * n + c])) {
                pivot_row = i;
                break;
            }
        }
        if (pivot_row == -1) continue;

        if (pivot_row != r) {
            for (int j = 0; j < n; j++) {
                Expr* tmp = matrix[r * n + j];
                matrix[r * n + j] = matrix[pivot_row * n + j];
                matrix[pivot_row * n + j] = tmp;
            }
        }

        /* Capture pivot value, then normalise the pivot row in place.
         * Pivot cell itself becomes 1. */
        Expr* pivot = expr_copy(matrix[r * n + c]);
        for (int j = 0; j < n; j++) {
            if (j == c) {
                expr_free(matrix[r * n + j]);
                matrix[r * n + j] = expr_new_integer(1);
            } else if (!is_zero_poly(matrix[r * n + j])) {
                Expr* old = matrix[r * n + j];
                matrix[r * n + j] = matsol_div_entry(old, pivot);
                expr_free(old);
            }
        }

        /* Eliminate column c in every other row. */
        for (int i = 0; i < m; i++) {
            if (i == r) continue;
            Expr* factor = matrix[i * n + c];
            if (is_zero_poly(factor)) continue;

            for (int j = 0; j < n; j++) {
                if (j == c) continue;
                /* M[i,j] <- M[i,j] - factor * M[r,j] */
                Expr* term = eval_and_free(expr_new_function(
                    expr_new_symbol("Times"),
                    (Expr*[]){expr_copy(factor), expr_copy(matrix[r * n + j])}, 2));
                Expr* neg_term = eval_and_free(expr_new_function(
                    expr_new_symbol("Times"),
                    (Expr*[]){expr_new_integer(-1), term}, 2));
                Expr* updated = eval_and_free(expr_new_function(
                    expr_new_symbol("Plus"),
                    (Expr*[]){expr_copy(matrix[i * n + j]), neg_term}, 2));
                Expr* updated_c = matsol_canon_entry(updated);
                expr_free(updated);
                expr_free(matrix[i * n + j]);
                matrix[i * n + j] = updated_c;
            }
            expr_free(matrix[i * n + c]);
            matrix[i * n + c] = expr_new_integer(0);
        }

        expr_free(pivot);
        r++;
    }

    Expr** rows = malloc(sizeof(Expr*) * (size_t)m);
    for (int i = 0; i < m; i++) {
        Expr** row_elems = malloc(sizeof(Expr*) * (size_t)n);
        for (int j = 0; j < n; j++) {
            row_elems[j] = matrix[i * n + j];
        }
        rows[i] = expr_new_function(expr_new_symbol("List"), row_elems, (size_t)n);
        free(row_elems);
    }
    Expr* result = expr_new_function(expr_new_symbol("List"), rows, (size_t)m);
    free(rows);
    free(matrix);
    return result;
}

/* ------------------------------------------------------------------
 * RowReduce -- CofactorExpansion
 *
 * For a square non-singular matrix the row-reduced echelon form is
 * the identity matrix.  We verify non-singularity by computing
 * Det[m] via Laplace cofactor expansion and checking it is not
 * structurally zero.
 *
 * For singular / rectangular / empty input the cofactor approach
 * does not give a closed-form RREF -- we emit a warning once per
 * distinct call and fall back to DivisionFreeRowReduction so the
 * user still gets a sensible answer.
 * ------------------------------------------------------------------ */
static Expr* rowreduce_cofactor(Expr* arg) {
    int64_t dims[64];
    int rank = get_tensor_dims(arg, dims);
    if (rank != 2 || dims[0] == 0 || dims[1] == 0) {
        return rowreduce_divfree(arg);
    }
    if (dims[0] != dims[1]) {
        static uint64_t last_hash = 0;
        matsol_warn_once(&last_hash, arg,
            "RowReduce::cofnsq: Method -> \"CofactorExpansion\" requires "
            "a non-singular square matrix; falling back to "
            "\"DivisionFreeRowReduction\".\n");
        return rowreduce_divfree(arg);
    }

    int n = (int)dims[0];
    Expr** flat = malloc(sizeof(Expr*) * (size_t)n * (size_t)n);
    size_t idx = 0;
    flatten_tensor(arg, flat, &idx);

    int* cols = malloc(sizeof(int) * (size_t)n);
    for (int i = 0; i < n; i++) cols[i] = i;
    Expr* det = laplace_det(flat, n, n, 0, cols);
    free(cols);

    bool singular = is_zero_poly(det);
    expr_free(det);
    for (size_t i = 0; i < idx; i++) expr_free(flat[i]);
    free(flat);

    if (singular) {
        static uint64_t last_hash = 0;
        matsol_warn_once(&last_hash, arg,
            "RowReduce::cofnsq: Method -> \"CofactorExpansion\" requires "
            "a non-singular square matrix; falling back to "
            "\"DivisionFreeRowReduction\".\n");
        return rowreduce_divfree(arg);
    }

    /* Build the n x n identity matrix in place. */
    Expr** rows = malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) {
        Expr** row_elems = malloc(sizeof(Expr*) * (size_t)n);
        for (int j = 0; j < n; j++) {
            row_elems[j] = expr_new_integer(i == j ? 1 : 0);
        }
        rows[i] = expr_new_function(expr_new_symbol("List"), row_elems, (size_t)n);
        free(row_elems);
    }
    Expr* result = expr_new_function(expr_new_symbol("List"), rows, (size_t)n);
    free(rows);
    return result;
}

/* ------------------------------------------------------------------
 * Build a nested List of shape {dims[0], ..., dims[rank-1]} from a
 * row-major flat array.  Consumes one element of `flat` per leaf,
 * advancing `*idx`; ownership of every consumed Expr* transfers into
 * the returned tree.
 * ------------------------------------------------------------------ */
static Expr* build_nested_list(Expr** flat, int64_t* dims, int rank,
                                size_t* idx) {
    if (rank == 0) {
        return flat[(*idx)++];
    }
    size_t n = (size_t)dims[0];
    Expr** items = NULL;
    if (n > 0) items = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        items[i] = build_nested_list(flat, dims + 1, rank - 1, idx);
    }
    Expr* result = expr_new_function(expr_new_symbol("List"), items, n);
    if (items) free(items);
    return result;
}

/* Assemble the final LinearSolve result from a flat answer array.
 *
 * `flat` has length `c * k` (k = product of trail_dims, or 1 when
 * trail_rank == 0) and is indexed `flat[j * k + kk]` for unknown
 * index j (0..c-1) and trailing flat index kk (0..k-1).  The
 * resulting List has shape `{c, trail_dims[0], ..., trail_dims[trail_rank-1]}`
 * when `trail_rank >= 1`, and is a flat List of length c when
 * `trail_rank == 0` (legacy vector-RHS shape).  Takes ownership of
 * every entry in `flat`.
 * ------------------------------------------------------------------ */
static Expr* build_solution_tensor(Expr** flat, int c,
                                    int64_t* trail_dims, int trail_rank) {
    int64_t full_dims[64];
    full_dims[0] = (int64_t)c;
    for (int i = 0; i < trail_rank; i++) full_dims[i + 1] = trail_dims[i];
    size_t idx = 0;
    return build_nested_list(flat, full_dims, trail_rank + 1, &idx);
}

/* ------------------------------------------------------------------
 * LinearSolve -- DivisionFreeRowReduction (Bareiss-like)
 *
 * Direct lift of the previous body of `builtin_linearsolve` from
 * src/linalg.c (see that file's commit history and comments for the
 * algorithmic rationale and the under-determined / inconsistent
 * handling).
 *
 * `trail_dims` / `trail_rank` describe the trailing (non-leading)
 * dimensions of `b` after the leading dims that match `m`'s shape.
 * For a legacy vector RHS (b is rank-1 against a rank-2 m) the trail
 * is empty (trail_rank == 0, k == 1) and the answer is returned as a
 * flat List of c entries.  For a matrix RHS (trail_rank == 1) the
 * answer has shape {c, k}.  For higher-rank b the answer has shape
 * {c, trail_dims[0], ..., trail_dims[trail_rank-1]}.
 * ------------------------------------------------------------------ */
static Expr* linearsolve_divfree(Expr* m, Expr* b,
                                  int64_t* trail_dims, int trail_rank,
                                  int r, int c, int k) {
    int cols = c + k;

    Expr** matrix = malloc(sizeof(Expr*) * (size_t)r * (size_t)cols);
    Expr** flat_m = malloc(sizeof(Expr*) * (size_t)r * (size_t)c);
    {
        size_t idx = 0;
        flatten_tensor(m, flat_m, &idx);
    }
    Expr** flat_b = malloc(sizeof(Expr*) * (size_t)r * (size_t)k);
    {
        size_t idx = 0;
        flatten_tensor(b, flat_b, &idx);
    }

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            matrix[i * cols + j] = flat_m[i * c + j];
        }
        for (int j = 0; j < k; j++) {
            matrix[i * cols + c + j] = flat_b[i * k + j];
        }
    }
    free(flat_m);
    free(flat_b);

    int* pivot_col_for_row = malloc(sizeof(int) * (size_t)r);
    for (int i = 0; i < r; i++) pivot_col_for_row[i] = -1;

    Expr* P = expr_new_integer(1);
    int row = 0;
    for (int col = 0; col < c && row < r; col++) {
        int pivot_row = -1;
        for (int i = row; i < r; i++) {
            if (!is_zero_poly(matrix[i * cols + col])) {
                pivot_row = i;
                break;
            }
        }
        if (pivot_row == -1) continue;

        if (pivot_row != row) {
            for (int j = 0; j < cols; j++) {
                Expr* tmp = matrix[row * cols + j];
                matrix[row * cols + j] = matrix[pivot_row * cols + j];
                matrix[pivot_row * cols + j] = tmp;
            }
        }

        Expr* pivot = matrix[row * cols + col];

        for (int i = 0; i < r; i++) {
            if (i == row) continue;
            Expr* M_ic = matrix[i * cols + col];
            if (is_zero_poly(M_ic)) {
                for (int j = 0; j < cols; j++) {
                    if (j == col) continue;
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
                    if (j == col) continue;
                    Expr* t1 = eval_and_free(expr_new_function(
                        expr_new_symbol("Times"),
                        (Expr*[]){expr_copy(pivot), expr_copy(matrix[i * cols + j])}, 2));
                    Expr* t2 = eval_and_free(expr_new_function(
                        expr_new_symbol("Times"),
                        (Expr*[]){expr_copy(M_ic), expr_copy(matrix[row * cols + j])}, 2));
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
            expr_free(matrix[i * cols + col]);
            matrix[i * cols + col] = expr_new_integer(0);
        }

        expr_free(P);
        P = expr_copy(pivot);
        pivot_col_for_row[row] = col;
        row++;
    }
    expr_free(P);

    int pivots = row;

    bool inconsistent = false;
    for (int i = pivots; i < r && !inconsistent; i++) {
        for (int j = 0; j < k; j++) {
            if (!is_zero_poly(matrix[i * cols + c + j])) {
                inconsistent = true;
                break;
            }
        }
    }
    if (inconsistent) {
        fprintf(stderr,
                "LinearSolve::nosol: Linear equation encountered "
                "that has no solution.\n");
        for (int i = 0; i < r * cols; i++) expr_free(matrix[i]);
        free(matrix);
        free(pivot_col_for_row);
        return NULL;
    }

    /* Normalise pivot rows' RHS by pivot value, cancelling a common
     * polynomial factor as RowReduce / Inverse do. */
    for (int i = 0; i < pivots; i++) {
        int pc = pivot_col_for_row[i];
        Expr* pivot = matrix[i * cols + pc];
        for (int j = c; j < cols; j++) {
            if (is_zero_poly(matrix[i * cols + j])) continue;

            Expr* num_val = matrix[i * cols + j];
            Expr* den_val = expr_copy(pivot);

            Expr* g = eval_and_free(expr_new_function(
                expr_new_symbol("PolynomialGCD"),
                (Expr*[]){expr_copy(num_val), expr_copy(den_val)}, 2));
            Expr* new_num = exact_div_wrapper(num_val, g);
            Expr* new_den = exact_div_wrapper(den_val, g);
            expr_free(g);

            if (new_den->type == EXPR_INTEGER && new_den->data.integer < 0) {
                Expr* t = eval_and_free(expr_new_function(
                    expr_new_symbol("Times"),
                    (Expr*[]){new_num, expr_new_integer(-1)}, 2));
                new_num = expr_expand(t);
                expr_free(t);
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

    int* pivot_for_col = malloc(sizeof(int) * (size_t)c);
    for (int j = 0; j < c; j++) pivot_for_col[j] = -1;
    for (int i = 0; i < pivots; i++) {
        pivot_for_col[pivot_col_for_row[i]] = i;
    }

    Expr** flat_answers = malloc(sizeof(Expr*) * (size_t)c * (size_t)k);
    for (int j = 0; j < c; j++) {
        for (int kk = 0; kk < k; kk++) {
            if (pivot_for_col[j] == -1) {
                flat_answers[j * k + kk] = expr_new_integer(0);
            } else {
                int pr = pivot_for_col[j];
                flat_answers[j * k + kk] = matrix[pr * cols + c + kk];
                matrix[pr * cols + c + kk] = NULL;
            }
        }
    }
    Expr* result = build_solution_tensor(flat_answers, c, trail_dims, trail_rank);
    free(flat_answers);

    for (int i = 0; i < r * cols; i++) {
        if (matrix[i]) expr_free(matrix[i]);
    }
    free(matrix);
    free(pivot_col_for_row);
    free(pivot_for_col);

    return result;
}

/* ------------------------------------------------------------------
 * LinearSolve -- OneStepRowReduction
 *
 * Standard Gauss-Jordan with one division per pivot on the augmented
 * matrix [m | b].  After elimination each pivot row's pivot entry is
 * 1 and the rest of its pivot column is 0; reading off the answer is
 * direct -- the i-th pivot column's variable is matrix[pivot_row, c+kk]
 * for rhs column kk.
 * ------------------------------------------------------------------ */
static Expr* linearsolve_onestep(Expr* m, Expr* b,
                                  int64_t* trail_dims, int trail_rank,
                                  int r, int c, int k) {
    int cols = c + k;

    Expr** matrix = malloc(sizeof(Expr*) * (size_t)r * (size_t)cols);
    Expr** flat_m = malloc(sizeof(Expr*) * (size_t)r * (size_t)c);
    {
        size_t idx = 0;
        flatten_tensor(m, flat_m, &idx);
    }
    Expr** flat_b = malloc(sizeof(Expr*) * (size_t)r * (size_t)k);
    {
        size_t idx = 0;
        flatten_tensor(b, flat_b, &idx);
    }
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++)
            matrix[i * cols + j] = flat_m[i * c + j];
        for (int j = 0; j < k; j++)
            matrix[i * cols + c + j] = flat_b[i * k + j];
    }
    free(flat_m);
    free(flat_b);

    int* pivot_col_for_row = malloc(sizeof(int) * (size_t)r);
    for (int i = 0; i < r; i++) pivot_col_for_row[i] = -1;

    int row = 0;
    for (int col = 0; col < c && row < r; col++) {
        int pivot_row = -1;
        for (int i = row; i < r; i++) {
            if (!is_zero_poly(matrix[i * cols + col])) {
                pivot_row = i;
                break;
            }
        }
        if (pivot_row == -1) continue;

        if (pivot_row != row) {
            for (int j = 0; j < cols; j++) {
                Expr* tmp = matrix[row * cols + j];
                matrix[row * cols + j] = matrix[pivot_row * cols + j];
                matrix[pivot_row * cols + j] = tmp;
            }
        }

        Expr* pivot = expr_copy(matrix[row * cols + col]);
        for (int j = 0; j < cols; j++) {
            if (j == col) {
                expr_free(matrix[row * cols + j]);
                matrix[row * cols + j] = expr_new_integer(1);
            } else if (!is_zero_poly(matrix[row * cols + j])) {
                Expr* old = matrix[row * cols + j];
                matrix[row * cols + j] = matsol_div_entry(old, pivot);
                expr_free(old);
            }
        }

        for (int i = 0; i < r; i++) {
            if (i == row) continue;
            Expr* factor = matrix[i * cols + col];
            if (is_zero_poly(factor)) continue;

            for (int j = 0; j < cols; j++) {
                if (j == col) continue;
                Expr* term = eval_and_free(expr_new_function(
                    expr_new_symbol("Times"),
                    (Expr*[]){expr_copy(factor), expr_copy(matrix[row * cols + j])}, 2));
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
            expr_free(matrix[i * cols + col]);
            matrix[i * cols + col] = expr_new_integer(0);
        }

        expr_free(pivot);
        pivot_col_for_row[row] = col;
        row++;
    }

    int pivots = row;

    /* Consistency check: rows past the last pivot must have all-zero RHS. */
    bool inconsistent = false;
    for (int i = pivots; i < r && !inconsistent; i++) {
        for (int j = 0; j < k; j++) {
            if (!is_zero_poly(matrix[i * cols + c + j])) {
                inconsistent = true;
                break;
            }
        }
    }
    if (inconsistent) {
        fprintf(stderr,
                "LinearSolve::nosol: Linear equation encountered "
                "that has no solution.\n");
        for (int i = 0; i < r * cols; i++) expr_free(matrix[i]);
        free(matrix);
        free(pivot_col_for_row);
        return NULL;
    }

    int* pivot_for_col = malloc(sizeof(int) * (size_t)c);
    for (int j = 0; j < c; j++) pivot_for_col[j] = -1;
    for (int i = 0; i < pivots; i++) {
        pivot_for_col[pivot_col_for_row[i]] = i;
    }

    Expr** flat_answers = malloc(sizeof(Expr*) * (size_t)c * (size_t)k);
    for (int j = 0; j < c; j++) {
        for (int kk = 0; kk < k; kk++) {
            if (pivot_for_col[j] == -1) {
                flat_answers[j * k + kk] = expr_new_integer(0);
            } else {
                int pr = pivot_for_col[j];
                flat_answers[j * k + kk] = matrix[pr * cols + c + kk];
                matrix[pr * cols + c + kk] = NULL;
            }
        }
    }
    Expr* result = build_solution_tensor(flat_answers, c, trail_dims, trail_rank);
    free(flat_answers);

    for (int i = 0; i < r * cols; i++) {
        if (matrix[i]) expr_free(matrix[i]);
    }
    free(matrix);
    free(pivot_col_for_row);
    free(pivot_for_col);

    return result;
}

/* ------------------------------------------------------------------
 * LinearSolve -- CofactorExpansion (Cramer's rule)
 *
 * Restricted to square m.  Computes det(m) via Laplace cofactor
 * expansion (the same primitive used by `Det`).  For each rhs column
 * kk and each variable j, x[j, kk] = det(m with column j replaced by
 * b[*, kk]) / det(m), computed via exact_div_wrapper.
 *
 * Errors:
 *   - non-square m  -> LinearSolve::cofnsq, return NULL.
 *   - det(m) == 0   -> LinearSolve::cofsng, return NULL.
 *
 * Time complexity: O(n! * (n columns of b + 1)).  Use only for small n.
 * ------------------------------------------------------------------ */
static Expr* linearsolve_cofactor(Expr* m, Expr* b,
                                   int64_t* trail_dims, int trail_rank,
                                   int r, int c, int k, Expr* call_key) {
    if (r != c) {
        static uint64_t last_hash = 0;
        matsol_warn_once(&last_hash, call_key,
            "LinearSolve::cofnsq: Method -> \"CofactorExpansion\" requires "
            "a non-singular square matrix.\n");
        return NULL;
    }
    int n = c;

    Expr** flat_m = malloc(sizeof(Expr*) * (size_t)n * (size_t)n);
    {
        size_t idx = 0;
        flatten_tensor(m, flat_m, &idx);
    }
    Expr** flat_b = malloc(sizeof(Expr*) * (size_t)r * (size_t)k);
    {
        size_t idx = 0;
        flatten_tensor(b, flat_b, &idx);
    }

    int* cols = malloc(sizeof(int) * (size_t)n);
    for (int i = 0; i < n; i++) cols[i] = i;
    Expr* det_m = laplace_det(flat_m, n, n, 0, cols);

    bool singular = is_zero_poly(det_m);
    if (singular) {
        static uint64_t last_hash = 0;
        matsol_warn_once(&last_hash, call_key,
            "LinearSolve::cofsng: Method -> \"CofactorExpansion\" requires "
            "a non-singular matrix; det is structurally zero.\n");
        expr_free(det_m);
        free(cols);
        for (int i = 0; i < n * n; i++) expr_free(flat_m[i]);
        free(flat_m);
        for (int i = 0; i < r * k; i++) expr_free(flat_b[i]);
        free(flat_b);
        return NULL;
    }

    /* Build a scratch n x n matrix that initially mirrors m; we will
     * swap column j with b's kk-th column, compute det, then swap
     * back so the structure is reusable for the next (j, kk). */
    Expr** scratch = malloc(sizeof(Expr*) * (size_t)n * (size_t)n);
    for (int i = 0; i < n * n; i++) scratch[i] = expr_copy(flat_m[i]);

    /* Compute x[j, kk] = det(scratch_j) / det_m. */
    Expr** results = malloc(sizeof(Expr*) * (size_t)n * (size_t)k);

    for (int kk = 0; kk < k; kk++) {
        for (int j = 0; j < n; j++) {
            /* Replace column j with b[:, kk]. */
            Expr** saved = malloc(sizeof(Expr*) * (size_t)n);
            for (int i = 0; i < n; i++) {
                saved[i] = scratch[i * n + j];
                scratch[i * n + j] = expr_copy(flat_b[i * k + kk]);
            }

            Expr* dj = laplace_det(scratch, n, n, 0, cols);
            Expr* xj = matsol_div_entry(dj, det_m);
            expr_free(dj);
            results[j * k + kk] = xj;

            /* Restore column j. */
            for (int i = 0; i < n; i++) {
                expr_free(scratch[i * n + j]);
                scratch[i * n + j] = saved[i];
            }
            free(saved);
        }
    }

    /* Free scratch + flats + det. */
    for (int i = 0; i < n * n; i++) expr_free(scratch[i]);
    free(scratch);
    for (int i = 0; i < n * n; i++) expr_free(flat_m[i]);
    free(flat_m);
    for (int i = 0; i < r * k; i++) expr_free(flat_b[i]);
    free(flat_b);
    expr_free(det_m);
    free(cols);

    /* Assemble result.  `results[j * k + kk]` is the flat answer array
     * already in the layout expected by build_solution_tensor. */
    Expr* result = build_solution_tensor(results, n, trail_dims, trail_rank);
    free(results);
    return result;
}

/* ------------------------------------------------------------------
 * Dispatch shells
 * ------------------------------------------------------------------ */
Expr* builtin_rowreduce(Expr* res) {
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
                "RowReduce::method: Method option value is not one of "
                "\"Automatic\", \"DivisionFreeRowReduction\", "
                "\"OneStepRowReduction\", \"CofactorExpansion\".\n");
            return NULL;
        }
    }

    switch (method) {
        case MATSOL_AUTOMATIC:
        case MATSOL_DIVFREE:    return rowreduce_divfree(arg);
        case MATSOL_ONESTEP:    return rowreduce_onestep(arg);
        case MATSOL_COFACTOR:   return rowreduce_cofactor(arg);
        case MATSOL_INVALID:    return NULL;  /* unreachable */
    }
    return NULL;
}

Expr* builtin_linearsolve(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;

    Expr* m = res->data.function.args[0];
    Expr* b = res->data.function.args[1];

    MatsolMethod method = MATSOL_AUTOMATIC;
    if (argc == 3) {
        method = matsol_parse_method_option(res->data.function.args[2]);
        if (method == MATSOL_INVALID) {
            static uint64_t last_warned = 0;
            matsol_warn_once(&last_warned, res,
                "LinearSolve::method: Method option value is not one of "
                "\"Automatic\", \"DivisionFreeRowReduction\", "
                "\"OneStepRowReduction\", \"CofactorExpansion\".\n");
            return NULL;
        }
    }

    /* Shape validation -- shared by every method.
     *
     * Mathematica's LinearSolve interprets a rank-N matrix `m`
     * (N >= 2) with dimensions {d_1, ..., d_{N-1}, n} as a
     * (d_1*...*d_{N-1}) x n linear system: the last dim is the
     * number of unknowns, the leading dims combine into rows.  The
     * right-hand side `b` then has dimensions
     * {d_1, ..., d_{N-1}, e_1, ..., e_p}: leading dims must match
     * `m`'s leading dims; trailing dims become the trailing dims of
     * the result, which has shape {n, e_1, ..., e_p}.  When `p == 0`
     * (b is rank N-1) the answer is a flat vector of length n. */
    int64_t m_dims[64];
    int m_rank = get_tensor_dims(m, m_dims);
    if (m_rank < 2) {
        char* m_str = expr_to_string(m);
        fprintf(stderr,
                "LinearSolve::matrix: Argument %s at position 1 is "
                "not a non-empty rectangular matrix.\n",
                m_str);
        free(m_str);
        return NULL;
    }
    for (int i = 0; i < m_rank; i++) {
        if (m_dims[i] == 0) {
            char* m_str = expr_to_string(m);
            fprintf(stderr,
                    "LinearSolve::matrix: Argument %s at position 1 is "
                    "not a non-empty rectangular matrix.\n",
                    m_str);
            free(m_str);
            return NULL;
        }
    }

    int lead = m_rank - 1;                /* leading dims of m  */
    int r = 1;
    for (int i = 0; i < lead; i++) r *= (int)m_dims[i];
    int c = (int)m_dims[lead];

    int64_t b_dims[64];
    int b_rank = get_tensor_dims(b, b_dims);
    if (b_rank < lead) {
        char* b_str = expr_to_string(b);
        fprintf(stderr,
                "LinearSolve::lvec: %s is neither a vector nor a matrix.\n",
                b_str);
        free(b_str);
        return NULL;
    }
    for (int i = 0; i < lead; i++) {
        if (b_dims[i] != m_dims[i]) {
            char* m_str = expr_to_string(m);
            char* b_str = expr_to_string(b);
            fprintf(stderr,
                    "LinearSolve::lvec1: Coefficient matrix and target "
                    "vector %s . x == %s do not have the same dimensions.\n",
                    m_str, b_str);
            free(m_str);
            free(b_str);
            return NULL;
        }
    }

    int trail_rank = b_rank - lead;
    int64_t* trail_dims = b_dims + lead;
    int k = 1;
    for (int i = 0; i < trail_rank; i++) k *= (int)trail_dims[i];

    switch (method) {
        case MATSOL_AUTOMATIC:
        case MATSOL_DIVFREE:
            return linearsolve_divfree(m, b, trail_dims, trail_rank, r, c, k);
        case MATSOL_ONESTEP:
            return linearsolve_onestep(m, b, trail_dims, trail_rank, r, c, k);
        case MATSOL_COFACTOR:
            return linearsolve_cofactor(m, b, trail_dims, trail_rank, r, c, k, res);
        case MATSOL_INVALID:
            return NULL;  /* unreachable */
    }
    return NULL;
}

/* ------------------------------------------------------------------
 * Registration
 * ------------------------------------------------------------------ */
void matsol_init(void) {
    symtab_add_builtin("RowReduce", builtin_rowreduce);
    symtab_get_def("RowReduce")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("LinearSolve", builtin_linearsolve);
    symtab_get_def("LinearSolve")->attributes |= ATTR_PROTECTED;
}
