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
#include "pack.h"
#include "eval.h"
#include "ndlinalg.h"
#include "sym_names.h"
#include "checked_int.h"   /* ci_powi_i64 — exact integer powers, overflow flag */
#include <math.h>          /* pow — matches power.c's real^integer branch */
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

/* Build a single entry x^e.  The exponent-0 column is the literal 1 (so 0^0
 * reads as 1); every other entry is a Power[x, e] node which the evaluator
 * simplifies (folding numeric powers and Power[x, 1] -> x).
 *
 * That leading 1 is INVENTED, so it takes the node list's exactness rather than
 * a default of exact Integer -- Mathematica's VandermondeMatrix[{1., 2., 3.}]
 * is {{1., 1., 1.}, {1., 2., 4.}, {1., 3., 9.}}, all Real. Writing an exact 1
 * there left an entire column of a machine matrix with the wrong head, which no
 * buffer can hold. See common.h on machine-real contagion. */
static Expr* vm_entry(Expr* node, int64_t exp, bool real) {
    if (exp == 0) {
        return real ? expr_new_real(1.0) : expr_new_integer(1);
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
    bool real = false;
    for (int64_t i = 0; i < n && !real; i++)
        if (nodes[i] && nodes[i]->type == EXPR_REAL) real = true;

    /* Buffer fast path — bit-identical to the nested Power/eval/pack path below,
     * but only for the two UNIFORM cases that path actually packs: all nodes
     * Integer (exact integer powers -> an int64 buffer) or all nodes machine
     * Real (pow(), matching power.c's real^integer branch -> a float64 buffer).
     *
     * The triggers are STRICTER than Hankel/Toeplitz's on purpose.  vm_entry
     * copies each node into Power[node, e] WITHOUT the machine-real coercion
     * hk_cell/tz_cell apply, so a MIXED list like {1, 2, 3.} keeps its Integer
     * nodes' powers Integer beside the Real ones -- the nested path leaves that
     * a mixed-head, UNPACKED List (VandermondeMatrix[{1,2,3.}] is
     * {{1.,1,1},{1.,2,4},{1.,3.,9.}}), which no uniform buffer reproduces.
     * Requiring every node to be the same numeric head keeps this exact.
     * Rational / symbolic / complex nodes fail both guards -> nested path. */
    bool all_int = true, all_real = true;
    for (int64_t i = 0; i < n; i++) {
        if (nodes[i]->type != EXPR_INTEGER) all_int = false;
        if (nodes[i]->type != EXPR_REAL)    all_real = false;
    }
    if (all_int) {
        int64_t dims[2] = { n, k };
        void* raw = NULL;
        Expr* packed = ndbuild_open(2, dims, NDT_INT64, &raw);
        if (packed) {
            int64_t* bi = (int64_t*)raw;
            bool overflow = false;
            for (int64_t i = 0; i < n && !overflow; i++) {
                int64_t base = nodes[i]->data.integer;
                for (int64_t j = 0; j < k; j++)
                    if (ci_powi_i64(base, j, &bi[i * k + j])) { overflow = true; break; }
            }
            if (!overflow) return packed;
            /* A power exceeded int64: abandon the buffer and fall through to the
             * exact nested path, which promotes to a bignum (and then pack_offer
             * declines, so the whole matrix stays exact). */
            expr_free(packed);
        }
    } else if (all_real) {
        int64_t dims[2] = { n, k };
        void* raw = NULL;
        Expr* packed = ndbuild_open(2, dims, NDT_FLOAT64, &raw);
        if (packed) {
            double* bd = (double*)raw;
            for (int64_t i = 0; i < n; i++) {
                double base = nodes[i]->data.real;
                for (int64_t j = 0; j < k; j++)
                    bd[i * k + j] = pow(base, (double)j);
            }
            return packed;
        }
    }

    Expr** rows = malloc(sizeof(Expr*) * (size_t)n);
    for (int64_t i = 0; i < n; i++) {
        Expr** cells = malloc(sizeof(Expr*) * (size_t)k);
        for (int64_t j = 0; j < k; j++) {
            cells[j] = vm_entry(nodes[i], j, real);
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), cells, (size_t)k);
        free(cells);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    free(rows);
    /* The cells are Power[x, e] NODES; the evaluator folds them. Doing that
     * here rather than on the way out is what lets the finished matrix be
     * offered for packing -- pack_offer sees machine numbers instead of Power
     * expressions, and declines harmlessly for a symbolic node list. Evaluating
     * once inside is idempotent: the evaluator re-evaluates the returned value
     * to a fixed point anyway. */
    return pack_offer(eval_and_free(result));
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
