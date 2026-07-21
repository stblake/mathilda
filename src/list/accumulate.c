#include "list_common.h"
#include "accumulate.h"
#include "assoc.h"
#include "ndarray.h"
#include "ndreduce.h"

/* Accumulate[list] returns the running cumulative totals of list, with the
 * same head and the same length as the input. The intermediate sums are
 * built as Plus[acc, next] and reduced by the evaluator so that integers,
 * arbitrary-precision bignums, machine doubles, lists (matrix columns via
 * Listable Plus), and symbolic expressions all combine correctly.
 *
 * An optional second argument of the form Method -> "CompensatedSummation"
 * triggers Kahan compensated summation when every element reduces to a
 * machine double (EXPR_REAL, EXPR_INTEGER, or EXPR_BIGINT). For other
 * inputs the option is silently ignored and the symbolic accumulation
 * is used. */
static bool accumulate_is_compensated_method(Expr* opt) {
    if (opt->type != EXPR_FUNCTION) return false;
    if (opt->data.function.head->type != EXPR_SYMBOL) return false;
    const char* hd = opt->data.function.head->data.symbol.name;
    if ((hd != SYM_Rule && hd != SYM_RuleDelayed) ||
        opt->data.function.arg_count != 2) return false;
    Expr* lhs = opt->data.function.args[0];
    Expr* rhs = opt->data.function.args[1];
    if (lhs->type != EXPR_SYMBOL || lhs->data.symbol.name != SYM_Method) return false;
    if (rhs->type != EXPR_STRING) return false;
    return strcmp(rhs->data.string, "CompensatedSummation") == 0;
}

static bool accumulate_to_double(Expr* e, double* out) {
    if (e->type == EXPR_REAL)    { *out = e->data.real;             return true; }
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer;  return true; }
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint); return true; }
    return false;
}

Expr* builtin_accumulate(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    Expr* lst = res->data.function.args[0];

    /* Accumulate[assoc] gives the running totals of the values, keeping every
     * key aligned with its cumulative sum. */
    if (is_association(lst)) { Expr* r = assoc_rekey_over_values(res); if (r) return r; }

    /* Accumulate[ndarray] prefix-sums the leading axis on the buffer (ndreduce.c). */
    if (argc == 1 && is_ndarray(lst)) return ndred_accumulate(res);

    if (lst->type != EXPR_FUNCTION) return NULL;

    bool kahan_requested = false;
    if (argc == 2) {
        if (!accumulate_is_compensated_method(res->data.function.args[1])) return NULL;
        kahan_requested = true;
    }

    Expr* head = lst->data.function.head;
    size_t n = lst->data.function.arg_count;

    if (n == 0) return expr_copy(lst);

    /* Kahan compensated summation in double precision when every element
     * is a machine number. */
    if (kahan_requested) {
        bool all_numeric = true;
        for (size_t i = 0; i < n; i++) {
            double tmp;
            if (!accumulate_to_double(lst->data.function.args[i], &tmp)) {
                all_numeric = false;
                break;
            }
        }
        if (all_numeric) {
            Expr** out = malloc(sizeof(Expr*) * n);
            if (!out) return NULL;
            double sum = 0.0;
            double c = 0.0;
            for (size_t i = 0; i < n; i++) {
                double x = 0.0;
                accumulate_to_double(lst->data.function.args[i], &x);
                double y = x - c;
                double t = sum + y;
                c = (t - sum) - y;
                sum = t;
                out[i] = expr_new_real(sum);
            }
            Expr* result = expr_new_function(expr_copy(head), out, n);
            free(out);
            return result;
        }
        /* Mixed/symbolic input: fall through to the standard accumulator. */
    }

    Expr** out = malloc(sizeof(Expr*) * n);
    if (!out) return NULL;

    out[0] = expr_copy(lst->data.function.args[0]);
    for (size_t i = 1; i < n; i++) {
        Expr** plus_args = malloc(sizeof(Expr*) * 2);
        if (!plus_args) {
            for (size_t j = 0; j < i; j++) expr_free(out[j]);
            free(out);
            return NULL;
        }
        plus_args[0] = expr_copy(out[i - 1]);
        plus_args[1] = expr_copy(lst->data.function.args[i]);
        Expr* plus_expr = expr_new_function(expr_new_symbol(SYM_Plus), plus_args, 2);
        free(plus_args);
        out[i] = evaluate(plus_expr);
        expr_free(plus_expr);
    }

    Expr* result = expr_new_function(expr_copy(head), out, n);
    free(out);
    return result;
}
