/* VandermondeMatrix — a matrix whose rows are the successive powers of a
 * sequence of nodes.
 *
 *   VandermondeMatrix[{x1, ..., xn}]      n x n Vandermonde matrix on the
 *                                         nodes x_i.
 *   VandermondeMatrix[{x1, ..., xn}, k]   n x k Vandermonde matrix.
 *
 * The (1-based) entry (i, j) is x_i^(j-1), so the first column is all ones,
 * the second column is the nodes themselves, the third their squares, and so
 * on.  The nodes need not be numerical and need not be distinct: the entries
 * are built as Power[x_i, j-1] nodes (the first column emitted as the literal
 * integer 1, so 0^0 reads as 1 to match the interpolation semantics) and the
 * evaluator then simplifies them — numeric powers fold to their value,
 * Power[x, 1] folds to x, leaving symbolic nodes as clean Power expressions.
 *
 * Vandermonde matrices arise in polynomial interpolation and in computing
 * moments in the monomial basis: LinearSolve[V, b] recovers the coefficients
 * of the polynomial through the points {x_i, b_i}.
 *
 * The single-argument structured-array conversion form, VandermondeMatrix[vmat],
 * is not supported (Mathilda has no structured-array representation); a single
 * matrix argument (a list of lists) is therefore left unevaluated.
 *
 * Diagnostics mirror Wolfram's surface text:
 *   - zero arguments  ->  VandermondeMatrix::argt
 */

#include "linalg.h"
#include "ndlinalg.h"
#include "sym_names.h"
#include <stdio.h>
#include <stdlib.h>

/* Recognise a `head[args...]` whose head is the symbol `sym`. */
static bool vm_is_call(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym;
}

/* Parse a single positive machine integer (Integer > 0). */
static bool vm_positive_int(const Expr* e, int64_t* out) {
    if (e->type == EXPR_INTEGER && e->data.integer > 0) {
        *out = e->data.integer;
        return true;
    }
    return false;
}

/* True when every element of the (non-empty) List is itself a List — i.e. the
 * argument is a matrix.  Such an argument is the unsupported structured-array
 * conversion form, not a node list, and is left unevaluated. */
static bool vm_is_matrix(const Expr* list) {
    size_t n = list->data.function.arg_count;
    if (n == 0) return false;
    for (size_t i = 0; i < n; i++) {
        if (!vm_is_call(list->data.function.args[i], SYM_List)) return false;
    }
    return true;
}

/* Build a single entry x^e.  The exponent-0 column is the literal integer 1
 * (so 0^0 reads as 1); every other entry is a Power[x, e] node which the
 * evaluator simplifies (folding numeric powers and Power[x, 1] -> x). */
static Expr* vm_entry(Expr* node, int64_t exp) {
    if (exp == 0) {
        return expr_new_integer(1);
    }
    Expr** args = malloc(sizeof(Expr*) * 2);
    args[0] = expr_copy(node);
    args[1] = expr_new_integer(exp);
    Expr* p = expr_new_function(expr_new_symbol(SYM_Power), args, 2);
    free(args);
    return p;
}

/* Build the n x k Vandermonde matrix as a List of Lists; entry (i, j) is
 * nodes[i]^(j) for 0-based i in [0, n) and j in [0, k).  Source nodes are
 * deep-copied; the input keeps its ownership. */
static Expr* vm_build(int64_t n, int64_t k, Expr* const* nodes) {
    Expr** rows = malloc(sizeof(Expr*) * (size_t)n);
    for (int64_t i = 0; i < n; i++) {
        Expr** cells = malloc(sizeof(Expr*) * (size_t)k);
        for (int64_t j = 0; j < k; j++) {
            cells[j] = vm_entry(nodes[i], j);
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), cells, (size_t)k);
        free(cells);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    free(rows);
    return result;
}

Expr* builtin_vandermondematrix(Expr* res) {
    if (linalg_call_has_ndarray(res)) return linalg_delist_and_reeval(res);
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 0) {
        fprintf(stderr,
                "VandermondeMatrix::argt: VandermondeMatrix called with 0 "
                "arguments; 1 or 2 arguments are expected.\n");
        return NULL;
    }

    Expr* a0 = res->data.function.args[0];

    /* The node list must be a non-empty flat List.  A list of lists is the
     * unsupported matrix-conversion form; leave it unevaluated. */
    if (!vm_is_call(a0, SYM_List)) return NULL;
    size_t nn = a0->data.function.arg_count;
    if (nn == 0 || vm_is_matrix(a0)) return NULL;

    /* Form 1: VandermondeMatrix[{x1, ..., xn}] — square. */
    if (argc == 1) {
        return vm_build((int64_t)nn, (int64_t)nn, a0->data.function.args);
    }

    /* Form 2: VandermondeMatrix[{x1, ..., xn}, k] — n x k. */
    int64_t k;
    if (argc == 2 && vm_positive_int(res->data.function.args[1], &k)) {
        return vm_build((int64_t)nn, k, a0->data.function.args);
    }

    /* Any other shape: leave the call unevaluated. */
    return NULL;
}
