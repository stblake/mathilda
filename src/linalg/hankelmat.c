/* HankelMatrix — a matrix that is constant along its antidiagonals.
 *
 *   HankelMatrix[n]            n x n Hankel matrix whose first row and first
 *                              column are the successive integers 1..n, with
 *                              zeros below the main antidiagonal.
 *   HankelMatrix[{c1,...,cm}]  m x m Hankel matrix whose first column is the
 *                              given list, with zeros below the antidiagonal.
 *   HankelMatrix[{c1,...,cm},  m x n Hankel matrix with the first list down
 *               {r1,...,rn}]   the first column and the second list across
 *                              the last row.
 *
 * The entry (i, j) is c_{i+j-1} when i+j-1 <= m, and r_{i+j-m} otherwise
 * (1-based indices).  The shared corner c_m and r_1 should be equal; if they
 * differ, the column element is used and a HankelMatrix::crs warning is
 * emitted.  Entries are copied verbatim, so symbolic, complex, exact and
 * inexact entries all flow through unchanged; arbitrary precision comes from
 * the entries themselves (e.g. `1`20`) or from wrapping in N.
 *
 * Diagnostics mirror Wolfram's surface text:
 *   - zero arguments  ->  HankelMatrix::argb
 *   - mismatched corner element  ->  HankelMatrix::crs (warning; still builds)
 */

#include "linalg.h"
#include "ndlinalg.h"
#include "sym_names.h"
#include "print.h"
#include <stdio.h>
#include <stdlib.h>

/* Recognise a `head[args...]` whose head is the symbol `sym`. */
static bool hk_is_call(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol == sym;
}

/* Parse a single positive machine integer (Integer > 0). */
static bool hk_positive_int(const Expr* e, int64_t* out) {
    if (e->type == EXPR_INTEGER && e->data.integer > 0) {
        *out = e->data.integer;
        return true;
    }
    return false;
}

/* Build the m x n Hankel matrix as a List of Lists.
 *
 * Entry (i, j) (0-based here; s = i + j + 1 is the 1-based antidiagonal
 * index i + j - 1) is:
 *   - the integer s (or 0 past the antidiagonal) when integer_form;
 *   - cvals[s - 1] when s <= m;
 *   - rvals[s - m] when rvals != NULL;
 *   - the integer 0 otherwise (square first-column form).
 * Source entries are deep-copied; the inputs keep their ownership. */
static Expr* hk_build(int64_t m, int64_t n,
                      Expr* const* cvals, Expr* const* rvals,
                      bool integer_form) {
    Expr** rows = malloc(sizeof(Expr*) * (size_t)m);
    for (int64_t i = 0; i < m; i++) {
        Expr** cells = malloc(sizeof(Expr*) * (size_t)n);
        for (int64_t j = 0; j < n; j++) {
            int64_t s = i + j + 1;  /* (i+1) + (j+1) - 1 */
            if (integer_form) {
                cells[j] = expr_new_integer(s <= m ? s : 0);
            } else if (s <= m) {
                cells[j] = expr_copy(cvals[s - 1]);
            } else if (rvals != NULL) {
                cells[j] = expr_copy(rvals[s - m]);
            } else {
                cells[j] = expr_new_integer(0);
            }
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), cells, (size_t)n);
        free(cells);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)m);
    free(rows);
    return result;
}

Expr* builtin_hankelmatrix(Expr* res) {
    if (linalg_call_has_ndarray(res)) return linalg_delist_and_reeval(res);
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 0) {
        fprintf(stderr,
                "HankelMatrix::argb: HankelMatrix called with 0 arguments; "
                "between 1 and 3 arguments are expected.\n");
        return NULL;
    }

    Expr* a0 = res->data.function.args[0];

    /* Form 1: HankelMatrix[n] — successive integers, zeros below. */
    int64_t nint;
    if (argc == 1 && hk_positive_int(a0, &nint)) {
        return hk_build(nint, nint, NULL, NULL, true);
    }

    /* Form 2: HankelMatrix[{c1,...,cm}] — square, first column given. */
    if (argc == 1 && hk_is_call(a0, SYM_List)) {
        size_t mm = a0->data.function.arg_count;
        if (mm == 0) return NULL;  /* empty list: leave unevaluated */
        return hk_build((int64_t)mm, (int64_t)mm,
                        a0->data.function.args, NULL, false);
    }

    /* Form 3: HankelMatrix[{c1,...,cm}, {r1,...,rn}] — first column and
     * last row given. */
    if (argc == 2 && hk_is_call(a0, SYM_List) &&
        hk_is_call(res->data.function.args[1], SYM_List)) {
        Expr* a1 = res->data.function.args[1];
        size_t mm = a0->data.function.arg_count;
        size_t nn = a1->data.function.arg_count;
        if (mm == 0 || nn == 0) return NULL;  /* leave unevaluated */

        /* The shared corner: c_m must match r_1.  If not, warn and use the
         * column element (the formula always reads c_m, never r_1). */
        Expr* clast = a0->data.function.args[mm - 1];
        Expr* rfirst = a1->data.function.args[0];
        if (!expr_eq(clast, rfirst)) {
            char* cs = expr_to_string(clast);
            char* rs = expr_to_string(rfirst);
            fprintf(stderr,
                    "HankelMatrix::crs: Warning: the column element %s and row "
                    "element %s at positions %lld and 1 are not the same. "
                    "Using column element.\n",
                    cs ? cs : "?", rs ? rs : "?", (long long)mm);
            free(cs);
            free(rs);
        }
        return hk_build((int64_t)mm, (int64_t)nn,
                        a0->data.function.args, a1->data.function.args, false);
    }

    /* Any other shape (e.g. a non-list, non-integer spec): leave the call
     * unevaluated so symbolic arguments flow through unchanged. */
    return NULL;
}
