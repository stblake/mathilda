/* nullspace.c
 *
 * NullSpace[m]                       -- list of basis vectors v such that
 *                                       m . v == 0.
 * NullSpace[m, Method -> "<name>"]   -- explicit RREF method dispatch.
 *
 * Implementation strategy:
 *
 *   1. Validate that `m` is a rank-2 non-empty matrix.  Wrong shape
 *      returns NULL so the call remains unevaluated (matching the
 *      RowReduce / LinearSolve convention).
 *
 *   2. Reduce `m` to reduced row echelon form (RREF) by calling
 *      RowReduce as a builtin.  Routing through the existing
 *      `builtin_rowreduce` dispatcher avoids duplicating the
 *      Bareiss-like / OneStep / Cofactor algorithms and means
 *      NullSpace inherits any future RREF improvements automatically.
 *
 *   3. Identify pivot columns from the RREF: the leftmost non-zero
 *      entry in each row marks one pivot column.
 *
 *   4. For each free (non-pivot) column f, build a basis vector v of
 *      length cols by setting:
 *        - v[f] = 1
 *        - v[p] = -RREF[row_of_p, f]  for each pivot column p
 *        - v[other free] = 0
 *      Iterate `f` from the rightmost free column to the leftmost so
 *      the basis ordering matches Mathematica.
 *
 *   5. For exact rational input, scale each vector by the LCM of its
 *      entries' integer denominators so the result is integer-valued
 *      whenever the input is integer-valued.  For symbolic / inexact
 *      input the vectors are left in their natural form.
 *
 * Memory ownership: standard builtin contract.  On success this file
 * owns `res` and frees it; on NULL return the caller retains
 * ownership of `res`.
 */

#include "nullspace.h"
#include "linalg.h"
#include "linsolve.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"
#include "poly.h"
#include "sym_names.h"

#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Build the empty result `List[]` used when `m` has full column rank. */
static Expr* empty_basis(void) {
    return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
}

/* Call `RowReduce[m, Method -> "<name>"]` via the evaluator.  When
 * `method` is MATSOL_AUTOMATIC we omit the option entirely so the
 * default (DivisionFreeRowReduction) is used.  Returned tree is owned
 * by the caller.  Returns NULL on evaluator failure (which shouldn't
 * happen for any well-formed matrix). */
static Expr* call_rowreduce(Expr* m, MatsolMethod method) {
    if (method == MATSOL_AUTOMATIC) {
        Expr** args = malloc(sizeof(Expr*) * 1);
        args[0] = expr_copy(m);
        Expr* call = expr_new_function(expr_new_symbol("RowReduce"), args, 1);
        free(args);
        return eval_and_free(call);
    }

    const char* name = NULL;
    switch (method) {
        case MATSOL_DIVFREE:  name = "DivisionFreeRowReduction"; break;
        case MATSOL_ONESTEP:  name = "OneStepRowReduction";      break;
        case MATSOL_COFACTOR: name = "CofactorExpansion";        break;
        default:              name = "Automatic";                break;
    }
    Expr** rule_args = malloc(sizeof(Expr*) * 2);
    rule_args[0] = expr_new_symbol(SYM_Method);
    rule_args[1] = expr_new_string(name);
    Expr* opt = expr_new_function(expr_new_symbol(SYM_Rule), rule_args, 2);
    free(rule_args);

    Expr** args = malloc(sizeof(Expr*) * 2);
    args[0] = expr_copy(m);
    args[1] = opt;
    Expr* call = expr_new_function(expr_new_symbol("RowReduce"), args, 2);
    free(args);
    return eval_and_free(call);
}

/* If every entry of v is an exact rational (Integer, Bignum, or
 * Rational[p, q]), multiply each entry by the LCM of all denominators
 * so the vector is integer-valued.  Otherwise leave the vector alone
 * (symbolic entries are returned as natural rationals).
 *
 * Mutates v[0..n-1] in place.  Takes ownership of each replaced entry. */
static void clear_int_denominators(Expr** v, int n) {
    mpz_t lcm;
    mpz_init_set_ui(lcm, 1);
    bool all_int_dens = true;

    for (int i = 0; i < n; i++) {
        Expr** den_args = malloc(sizeof(Expr*) * 1);
        den_args[0] = expr_copy(v[i]);
        Expr* d = eval_and_free(expr_new_function(
            expr_new_symbol("Denominator"), den_args, 1));
        free(den_args);

        if (d->type == EXPR_INTEGER) {
            mpz_t di;
            mpz_init(di);
            mpz_set_si(di, d->data.integer);
            mpz_lcm(lcm, lcm, di);
            mpz_clear(di);
        } else if (d->type == EXPR_BIGINT) {
            mpz_lcm(lcm, lcm, d->data.bigint);
        } else {
            all_int_dens = false;
            expr_free(d);
            break;
        }
        expr_free(d);
    }

    if (!all_int_dens || mpz_cmp_ui(lcm, 1) == 0) {
        mpz_clear(lcm);
        return;
    }

    Expr* lcm_expr = expr_new_bigint_from_mpz(lcm);
    mpz_clear(lcm);
    lcm_expr = expr_bigint_normalize(lcm_expr);

    for (int i = 0; i < n; i++) {
        Expr** mul_args = malloc(sizeof(Expr*) * 2);
        mul_args[0] = v[i];                       /* steal ownership */
        mul_args[1] = expr_copy(lcm_expr);
        Expr* new_val = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Times), mul_args, 2));
        free(mul_args);
        v[i] = new_val;
    }
    expr_free(lcm_expr);
}

/* The shared NullSpace core, parameterised by the RREF method. */
static Expr* nullspace_core(Expr* m, MatsolMethod method) {
    int64_t dims[64];
    int rank = get_tensor_dims(m, dims);

    if (rank != 2 || dims[0] == 0 || dims[1] == 0) {
        char* m_str = expr_to_string(m);
        fprintf(stderr,
                "NullSpace::matrix: Argument %s at position 1 is not a "
                "non-empty rectangular matrix.\n",
                m_str);
        free(m_str);
        return NULL;
    }

    int rows = (int)dims[0];
    int cols = (int)dims[1];

    Expr* rref = call_rowreduce(m, method);
    if (!rref) return NULL;

    /* RowReduce returns a {rows, cols} List of Lists.  Re-probe the
     * shape because future variants of RREF could (in principle)
     * return a rank-0 / non-list result on a degenerate input. */
    int64_t rref_dims[64];
    int rref_rank = get_tensor_dims(rref, rref_dims);
    if (rref_rank != 2 || rref_dims[0] != rows || rref_dims[1] != cols) {
        expr_free(rref);
        return NULL;
    }

    Expr** flat = malloc(sizeof(Expr*) * (size_t)rows * (size_t)cols);
    {
        size_t idx = 0;
        flatten_tensor(rref, flat, &idx);
    }
    expr_free(rref);

    /* Locate the pivot column of each row: the leftmost non-zero
     * entry.  After a successful RowReduce that entry is 1, but we
     * never depend on the value -- only on its position. */
    int* pivot_for_col = malloc(sizeof(int) * (size_t)cols);
    for (int c = 0; c < cols; c++) pivot_for_col[c] = -1;

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            if (!is_zero_poly(flat[r * cols + c])) {
                /* Only record the first pivot we see for this column;
                 * a well-formed RREF guarantees each column has at
                 * most one pivot row, but guard against degenerate
                 * inputs anyway. */
                if (pivot_for_col[c] == -1) pivot_for_col[c] = r;
                break;
            }
        }
    }

    /* Build basis vectors -- one per free column, rightmost first. */
    Expr** basis_rows = NULL;
    int basis_count = 0;
    int basis_cap = 0;

    for (int f = cols - 1; f >= 0; f--) {
        if (pivot_for_col[f] != -1) continue;

        Expr** v = malloc(sizeof(Expr*) * (size_t)cols);
        for (int j = 0; j < cols; j++) {
            if (j == f) {
                v[j] = expr_new_integer(1);
            } else if (pivot_for_col[j] != -1) {
                int pr = pivot_for_col[j];
                Expr** neg_args = malloc(sizeof(Expr*) * 2);
                neg_args[0] = expr_new_integer(-1);
                neg_args[1] = expr_copy(flat[pr * cols + f]);
                v[j] = eval_and_free(expr_new_function(
                    expr_new_symbol(SYM_Times), neg_args, 2));
                free(neg_args);
            } else {
                v[j] = expr_new_integer(0);
            }
        }

        clear_int_denominators(v, cols);

        Expr* row = expr_new_function(expr_new_symbol(SYM_List), v, (size_t)cols);
        free(v);

        if (basis_count >= basis_cap) {
            basis_cap = basis_cap == 0 ? 4 : basis_cap * 2;
            basis_rows = realloc(basis_rows, sizeof(Expr*) * (size_t)basis_cap);
        }
        basis_rows[basis_count++] = row;
    }

    for (int i = 0; i < rows * cols; i++) expr_free(flat[i]);
    free(flat);
    free(pivot_for_col);

    Expr* result;
    if (basis_count == 0) {
        free(basis_rows);
        result = empty_basis();
    } else {
        result = expr_new_function(expr_new_symbol(SYM_List),
                                    basis_rows, (size_t)basis_count);
        free(basis_rows);
    }
    return result;
}

Expr* builtin_nullspace(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    Expr* m = res->data.function.args[0];

    MatsolMethod method = MATSOL_AUTOMATIC;
    if (argc == 2) {
        method = matsol_parse_method_option(res->data.function.args[1]);
        if (method == MATSOL_INVALID) {
            static uint64_t last_warned = 0;
            matsol_warn_once(&last_warned, res,
                "NullSpace::method: Method option value is not one of "
                "\"Automatic\", \"DivisionFreeRowReduction\", "
                "\"OneStepRowReduction\", \"CofactorExpansion\".\n");
            return NULL;
        }
    }

    /* On success the evaluator (eval.c:813) frees `res` for us; on NULL
     * return it retains ownership.  Do not free `res` here. */
    return nullspace_core(m, method);
}

void matnull_init(void) {
    symtab_add_builtin("NullSpace", builtin_nullspace);
    symtab_get_def("NullSpace")->attributes |= ATTR_PROTECTED;
}
