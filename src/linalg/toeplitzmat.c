/* ToeplitzMatrix — a matrix that is constant along its diagonals.
 *
 *   ToeplitzMatrix[n]            n x n Toeplitz matrix whose first row and
 *                                first column are the successive integers
 *                                1..n (symmetric: entry (i, j) is |i - j| + 1).
 *   ToeplitzMatrix[{c1,...,cn}]  n x n symmetric Toeplitz matrix whose first
 *                                column (and first row) is the given list.
 *   ToeplitzMatrix[{c1,...,cm},  m x n Toeplitz matrix with the first list
 *                  {r1,...,rn}]  down the first column and the second list
 *                                across the first row.
 *
 * The entry (i, j) is c_{i-j+1} when i >= j, and r_{j-i+1} otherwise (1-based
 * indices).  The shared corner c_1 and r_1 should be equal; if they differ,
 * the column element is used and a ToeplitzMatrix::crs warning is emitted (the
 * formula always reads c_1 on the diagonal, never r_1).  Entries are copied
 * verbatim, so symbolic, complex, exact and inexact entries all flow through
 * unchanged; arbitrary precision comes from the entries themselves (e.g.
 * `1`20`) or from wrapping in N.
 *
 * Diagnostics mirror Wolfram's surface text:
 *   - zero arguments  ->  ToeplitzMatrix::argb
 *   - mismatched corner element  ->  ToeplitzMatrix::crs (warning; still builds)
 */

#include "linalg.h"
#include "ndlinalg.h"
#include "sym_names.h"
#include "print.h"
#include <stdio.h>
#include <stdlib.h>

/* Recognise a `head[args...]` whose head is the symbol `sym`. */
static bool tz_is_call(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym;
}

/* Parse a single positive machine integer (Integer > 0). */
static bool tz_positive_int(const Expr* e, int64_t* out) {
    if (e->type == EXPR_INTEGER && e->data.integer > 0) {
        *out = e->data.integer;
        return true;
    }
    return false;
}

/* Build the m x n Toeplitz matrix as a List of Lists.
 *
 * Entry (i, j) (0-based here) is:
 *   - the integer |i - j| + 1 when integer_form;
 *   - cvals[i - j] when i >= j (the (i-j+1)-th column entry);
 *   - rvals[j - i] otherwise (the (j-i+1)-th row entry).
 * Source entries are deep-copied; the inputs keep their ownership.  The c
 * index reaches m-1 (needs |c| = m) and the r index reaches n-1 (needs
 * |r| = n), which the callers guarantee. */
static Expr* tz_build(int64_t m, int64_t n,
                      Expr* const* cvals, Expr* const* rvals,
                      bool integer_form) {
    Expr** rows = malloc(sizeof(Expr*) * (size_t)m);
    for (int64_t i = 0; i < m; i++) {
        Expr** cells = malloc(sizeof(Expr*) * (size_t)n);
        for (int64_t j = 0; j < n; j++) {
            if (integer_form) {
                int64_t d = i >= j ? i - j : j - i;
                cells[j] = expr_new_integer(d + 1);
            } else if (i >= j) {
                cells[j] = expr_copy(cvals[i - j]);
            } else {
                cells[j] = expr_copy(rvals[j - i]);
            }
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), cells, (size_t)n);
        free(cells);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)m);
    free(rows);
    return result;
}

Expr* builtin_toeplitzmatrix(Expr* res) {
    if (linalg_call_has_ndarray(res)) return linalg_delist_and_reeval(res);
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 0) {
        fprintf(stderr,
                "ToeplitzMatrix::argb: ToeplitzMatrix called with 0 arguments; "
                "between 1 and 3 arguments are expected.\n");
        return NULL;
    }

    Expr* a0 = res->data.function.args[0];

    /* Form 1: ToeplitzMatrix[n] — successive integers, symmetric. */
    int64_t nint;
    if (argc == 1 && tz_positive_int(a0, &nint)) {
        return tz_build(nint, nint, NULL, NULL, true);
    }

    /* Form 2: ToeplitzMatrix[{c1,...,cn}] — symmetric, first column given
     * (the first row is the same list). */
    if (argc == 1 && tz_is_call(a0, SYM_List)) {
        size_t nn = a0->data.function.arg_count;
        if (nn == 0) return NULL;  /* empty list: leave unevaluated */
        return tz_build((int64_t)nn, (int64_t)nn,
                        a0->data.function.args, a0->data.function.args, false);
    }

    /* Form 3: ToeplitzMatrix[{c1,...,cm}, {r1,...,rn}] — first column and
     * first row given. */
    if (argc == 2 && tz_is_call(a0, SYM_List) &&
        tz_is_call(res->data.function.args[1], SYM_List)) {
        Expr* a1 = res->data.function.args[1];
        size_t mm = a0->data.function.arg_count;
        size_t nn = a1->data.function.arg_count;
        if (mm == 0 || nn == 0) return NULL;  /* leave unevaluated */

        /* The shared corner: c_1 must match r_1.  If not, warn and use the
         * column element (the formula always reads c_1, never r_1). */
        Expr* cfirst = a0->data.function.args[0];
        Expr* rfirst = a1->data.function.args[0];
        if (!expr_eq(cfirst, rfirst)) {
            char* cs = expr_to_string(cfirst);
            char* rs = expr_to_string(rfirst);
            fprintf(stderr,
                    "ToeplitzMatrix::crs: Warning: the column element %s and row "
                    "element %s at positions 1 and 1 are not the same. "
                    "Using column element.\n",
                    cs ? cs : "?", rs ? rs : "?");
            free(cs);
            free(rs);
        }
        return tz_build((int64_t)mm, (int64_t)nn,
                        a0->data.function.args, a1->data.function.args, false);
    }

    /* Any other shape (e.g. a non-list, non-integer spec): leave the call
     * unevaluated so symbolic arguments flow through unchanged. */
    return NULL;
}
