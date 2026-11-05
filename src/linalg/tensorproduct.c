/* TensorProduct[t1, t2, ...] -- the tensor product of arrays.
 *
 * Reference semantics (probed): Flat, OneIdentity, Protected. For arrays the
 * product is the full outer product, rank r1 + r2 + ...
 * (TensorProduct[{1, 2}, {3, 4}] is {{3, 4}, {6, 8}}); a NUMBER multiplies
 * (TensorProduct[2, {1, 2}] is {2, 4}); TensorProduct[t] is t and
 * TensorProduct[] is 1; a symbol stays a formal factor (TensorProduct[x, {1, 2}]
 * and TensorProduct[a, b] are left as they are). Any SparseArray factor makes
 * the answer a SparseArray.
 *
 * The outer product itself is Outer[Times, a, b] -- the buffer-aware,
 * sparse-native Outer -- so a packed or SparseArray factor never round-trips
 * through boxed Lists here. ndla_tensorproduct is the two-array entry used by
 * the NDArray surface and by Compile (ND_FN2S, rank rule R2_OUTER).
 */
#include "linalg.h"
#include "ndlinalg.h"
#include "sparse/sparse.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

static bool is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}
static bool is_array(const Expr* e) {
    return is_list(e) || e->type == EXPR_NDARRAY || e->type == EXPR_SPARSE;
}
static bool is_number(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        (e->data.function.head->data.symbol.name == SYM_Rational ||
         e->data.function.head->data.symbol.name == SYM_Complex)) return true;
    return false;
}

static Expr* eval_owned(Expr* e) {
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

/* a (x) b for two arrays: Outer[Times, a, b]. */
static Expr* outer2(const Expr* a, const Expr* b) {
    Expr* args[3] = { expr_new_symbol("Times"), expr_copy((Expr*)a), expr_copy((Expr*)b) };
    return eval_owned(expr_new_function(expr_new_symbol("Outer"), args, 3));
}

Expr* ndla_tensorproduct(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    const Expr* a = res->data.function.args[0];
    const Expr* b = res->data.function.args[1];
    if (!is_array(a) || !is_array(b)) return NULL;
    return outer2(a, b);
}

Expr* builtin_tensorproduct(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 0) return expr_new_integer(1);
    Expr* const* a = res->data.function.args;
    if (n == 1) return expr_copy(a[0]);
    bool any_sparse = false;
    size_t narr = 0;
    for (size_t i = 0; i < n; i++) {
        if (is_array(a[i])) { narr++; if (a[i]->type == EXPR_SPARSE) any_sparse = true; }
        else if (!is_number(a[i])) return NULL;      /* a formal factor: leave the product */
    }
    if (narr == 0) return NULL;                      /* numbers only: Times is not ours to do */
    /* the numeric scalar factor, then the arrays left to right */
    Expr* scal = NULL;
    Expr* acc = NULL;
    for (size_t i = 0; i < n; i++) {
        if (!is_array(a[i])) {
            if (!scal) scal = expr_copy(a[i]);
            else {
                Expr* t[2] = { scal, expr_copy(a[i]) };
                scal = eval_owned(expr_new_function(expr_new_symbol("Times"), t, 2));
            }
            continue;
        }
        if (!acc) { acc = expr_copy(a[i]); continue; }
        Expr* nx = outer2(acc, a[i]);
        expr_free(acc);
        if (!nx) { expr_free(scal); return NULL; }
        acc = nx;
    }
    if (scal) {
        Expr* t[2] = { scal, acc };
        acc = eval_owned(expr_new_function(expr_new_symbol("Times"), t, 2));
    }
    if (any_sparse && acc && acc->type != EXPR_SPARSE) {
        Expr* t[1] = { acc };
        acc = eval_owned(expr_new_function(expr_new_symbol("SparseArray"), t, 1));
    }
    return acc;
}
