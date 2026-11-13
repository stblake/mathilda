#include "list_common.h"
#include "flatten.h"
#include "ndarray.h"
#include "ndstruct.h"
#include "message.h"

static void flatten_rec(Expr* e, const char* h, int64_t level, Expr*** results, size_t* count, size_t* cap) {
    /* A visible NDArray nested inside an ordinary List. `head_is` is false for
     * it -- an EXPR_NDARRAY is not an EXPR_FUNCTION headed List -- so it was
     * collected as one opaque element and `Flatten[{ndA, ndB}]` answered with
     * the list unchanged where `Flatten[{{1., 2.}, {3., 4.}}]` gives
     * `{1., 2., 3., 4.}`. It IS a list of values by every other measure
     * (`ArrayQ`, `Dimensions`, `Length` all say so), so descending is the
     * answer that agrees with the other surface.
     *
     * This is why the sweep read `Dot`/`Inverse`/`LinearSolve` on a 6x6 as
     * ND-UNSUPPORTED: all three answered correctly and it was the probe's own
     * checksum, `N[Total[Flatten[{r}]]]` over a List of 200 result arrays, that
     * could not reduce to a number.
     *
     * The PACKED form cannot reach here -- the no-nesting invariant means a
     * packed node never sits inside a plain List (pack.c) -- so this is a
     * visible-surface repair only, and Flatten[ndarray] at the top level still
     * takes ndstruct_flatten's reshape below. */
    if (level != 0 && is_ndarray(e) && h == SYM_List) {
        Expr* rows = ndarray_to_nested_list(e);
        for (size_t i = 0; i < rows->data.function.arg_count; i++)
            flatten_rec(rows->data.function.args[i], h,
                        level == -1 ? -1 : level - 1, results, count, cap);
        expr_free(rows);
        return;
    }
    if (level != 0 && head_is(e, intern_symbol(h))) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            flatten_rec(e->data.function.args[i], h, level == -1 ? -1 : level - 1, results, count, cap);
        }
    } else {
        if (*count == *cap) {
            *cap *= 2;
            *results = realloc(*results, sizeof(Expr*) * (*cap));
        }
        (*results)[(*count)++] = expr_copy(e);
    }
}

bool flatten_level_is_infinity(const Expr* lv) {
    if (lv->type == EXPR_SYMBOL) return lv->data.symbol.name == SYM_Infinity;
    return head_is(lv, SYM_DirectedInfinity) && lv->data.function.arg_count == 1 &&
           lv->data.function.args[0]->type == EXPR_INTEGER &&
           lv->data.function.args[0]->data.integer == 1;
}

void flatten_flev_message(const Expr* call) {
    char* ls = expr_to_string_outputform(call->data.function.args[1]);
    char* cs = expr_to_string_outputform((Expr*)call);
    mth_message("Flatten", "flev",
                "The level argument %s in position 2 of %s should be a non-negative integer or "
                "Infinity giving the levels to flatten through or a list of lists of levels to "
                "flatten together.", ls ? ls : "", cs ? cs : "");
    free(ls);
    free(cs);
}

Expr* builtin_flatten(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 3) return NULL;

    Expr* list = res->data.function.args[0];

    /* The level: a non-negative integer, or Infinity (every level). A negative
     * integer is Flatten::flev and the call stays as it is -- checked before
     * any fast path, so a buffer argument gets the same answer. */
    int64_t n = -1; // -1 means infinity
    if (res->data.function.arg_count >= 2) {
        Expr* lv = res->data.function.args[1];
        if (flatten_level_is_infinity(lv)) {
            n = -1;
        } else {
            if (lv->type != EXPR_INTEGER) return NULL;
            n = lv->data.integer;
            if (n < 0) {
                /* A buffer is reported on its List form, once: declining with
                 * the buffer in place would let the evaluator's post-gate
                 * materialise it and ask (and report) a second time. */
                if (is_ndarray(list)) return ndarray_delist_and_reeval(res);
                flatten_flev_message(res);
                return NULL;
            }
        }
    }

    /* Flatten[ndarray]: reshape the flat buffer to rank 1 (see ndstruct.c). */
    if (is_ndarray(list)) return ndstruct_flatten(res);

    /* a Rational/Complex number is an atom: ::normal, unevaluated. (Other atoms
     * keep the identity result, which the reference also rejects with ::normal:
     * internal code relies on Flatten[atom] === atom, so that stays.) */
    if (expr_is_number_atom(list)) return builtin_normal_error("Flatten", res, 1);
    if (list->type != EXPR_FUNCTION) return expr_copy(list);

    const char* h = SYM_List;
    if (res->data.function.arg_count == 3) {
        if (res->data.function.args[2]->type != EXPR_SYMBOL) return NULL;
        h = res->data.function.args[2]->data.symbol.name;
    }

    size_t cap = 16;
    size_t count = 0;
    Expr** results = malloc(sizeof(Expr*) * cap);

    // Initial call: we flatten children of the head if they also have head h.
    for (size_t i = 0; i < list->data.function.arg_count; i++) {
        flatten_rec(list->data.function.args[i], h, n, &results, &count, &cap);
    }

    Expr* result = expr_new_function(expr_copy(list->data.function.head), results, count);
    free(results);
    return result;
}
