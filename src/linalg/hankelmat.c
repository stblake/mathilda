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
#include "common.h"
#include "pack.h"
#include "ndlinalg.h"
#include "sym_names.h"
#include "print.h"
#include <stdio.h>
#include <stdlib.h>

/* Recognise a `head[args...]` whose head is the symbol `sym`. */
static bool hk_is_call(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym;
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
/* True when any of the `n` given values is a machine Real, so the pad zeros
 * this builder invents -- and the exact entries it was given -- are Real too.
 * See common.h on machine-real contagion; Mathematica's
 * HankelMatrix[{1, 2, 3.}] is {{1., 2., 3.}, {2., 3., 0.}, {3., 0., 0.}}, all
 * Real, where this used to write an exact 0 in the corner and leave a matrix of
 * two heads that no buffer can hold. */
static bool hk_any_real(Expr* const* vals, int64_t n) {
    if (!vals) return false;
    for (int64_t i = 0; i < n; i++)
        if (vals[i] && vals[i]->type == EXPR_REAL) return true;
    return false;
}

/* One entry, coerced when the data is machine-real. */
static Expr* hk_cell(const Expr* e, bool real) {
    double d;
    if (real && common_machine_real_value(e, &d)) return expr_new_real(d);
    return expr_copy((Expr*)e);
}

/* Every value a machine Integer (int64, not bignum) / machine-real-valued —
 * the two cases the buffer fast path in hk_build can write directly. */
static bool hk_all_machine_int(Expr* const* v, int64_t n) {
    for (int64_t i = 0; i < n; i++)
        if (!v[i] || v[i]->type != EXPR_INTEGER) return false;
    return true;
}
static bool hk_all_machine_real(Expr* const* v, int64_t n) {
    double d;
    for (int64_t i = 0; i < n; i++)
        if (!v[i] || !common_machine_real_value(v[i], &d)) return false;
    return true;
}

static Expr* hk_build(int64_t m, int64_t n,
                      Expr* const* cvals, Expr* const* rvals,
                      bool integer_form, bool real) {
    /* Buffer fast path — bit-identical to the nested build below (same cell
     * formula, same machine-real contagion), it just writes the flat row-major
     * machine buffer directly instead of allocating m*n Expr cells and offering
     * them to the packer.  An all-integer or all-machine-real Hankel matrix is a
     * uniform machine array; complex / symbolic / bignum entries fail both
     * guards and take the nested path (where pack_offer packs or declines as it
     * always has). */
    bool all_int = integer_form ||
        (!real && hk_all_machine_int(cvals, m) && (!rvals || hk_all_machine_int(rvals, n)));
    bool all_real = real && hk_all_machine_real(cvals, m) &&
        (!rvals || hk_all_machine_real(rvals, n));
    if (all_int || all_real) {
        int64_t dims[2] = { m, n };
        void* raw = NULL;
        Expr* packed = ndbuild_open(2, dims, all_int ? NDT_INT64 : NDT_FLOAT64, &raw);
        if (packed) {
            int64_t* bi = (int64_t*)raw;
            double*  bd = (double*)raw;
            for (int64_t i = 0; i < m; i++) {
                for (int64_t j = 0; j < n; j++) {
                    int64_t s = i + j + 1;  /* 1-based antidiagonal index */
                    const Expr* el = NULL;
                    if (!integer_form) {
                        if (s <= m)      el = cvals[s - 1];
                        else if (rvals)  el = rvals[s - m];
                    }
                    if (all_int) {
                        bi[i * n + j] = integer_form ? (s <= m ? s : 0)
                                      : el ? el->data.integer : 0;
                    } else {
                        double d = 0.0;
                        if (el) common_machine_real_value(el, &d);
                        bd[i * n + j] = d;
                    }
                }
            }
            return packed;
        }
    }

    Expr** rows = malloc(sizeof(Expr*) * (size_t)m);
    for (int64_t i = 0; i < m; i++) {
        Expr** cells = malloc(sizeof(Expr*) * (size_t)n);
        for (int64_t j = 0; j < n; j++) {
            int64_t s = i + j + 1;  /* (i+1) + (j+1) - 1 */
            if (integer_form) {
                cells[j] = expr_new_integer(s <= m ? s : 0);
            } else if (s <= m) {
                cells[j] = hk_cell(cvals[s - 1], real);
            } else if (rvals != NULL) {
                cells[j] = hk_cell(rvals[s - m], real);
            } else {
                cells[j] = real ? expr_new_real(0.0) : expr_new_integer(0);
            }
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), cells, (size_t)n);
        free(cells);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)m);
    free(rows);
    /* Uniform now, so it can be a buffer -- offered rather than built directly
     * because the entries may be symbolic, which pack_offer declines. */
    return pack_offer(result);
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
        return hk_build(nint, nint, NULL, NULL, true, false);
    }

    /* Form 2: HankelMatrix[{c1,...,cm}] — square, first column given. */
    if (argc == 1 && hk_is_call(a0, SYM_List)) {
        size_t mm = a0->data.function.arg_count;
        if (mm == 0) return NULL;  /* empty list: leave unevaluated */
        return hk_build((int64_t)mm, (int64_t)mm,
                        a0->data.function.args, NULL, false,
                        hk_any_real(a0->data.function.args, (int64_t)mm));
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
                        a0->data.function.args, a1->data.function.args, false,
                        hk_any_real(a0->data.function.args, (int64_t)mm) ||
                        hk_any_real(a1->data.function.args, (int64_t)nn));
    }

    /* Any other shape (e.g. a non-list, non-integer spec): leave the call
     * unevaluated so symbolic arguments flow through unchanged. */
    return NULL;
}
