#include "ds_ops.h"
#include "list_common.h"
#include "reverse.h"
#include "ndarray.h"
#include "ndstruct.h"
#include "assoc.h"
#include "linkedlist.h"
#include "message.h"
#include "print.h"

static bool should_reverse_at_level(Expr* level_spec, size_t current_level) {
    if (!level_spec) return current_level == 1;
    if (level_spec->type == EXPR_INTEGER) return (size_t)level_spec->data.integer == current_level;
    if (level_spec->type == EXPR_FUNCTION && level_spec->data.function.head->data.symbol.name == SYM_List) {
        for (size_t i = 0; i < level_spec->data.function.arg_count; i++) {
            if (level_spec->data.function.args[i]->type == EXPR_INTEGER && 
                (size_t)level_spec->data.function.args[i]->data.integer == current_level) return true;
        }
    }
    return false;
}

/* Rational[p, q] and Complex[a, b] are numbers -- atoms -- even though they are
 * stored as two-argument nodes: Reverse never reaches inside one, so
 * Reverse[{5/4, 7/2}, 2] is the list unchanged (not {4/5, 2/7}). */
static bool is_number_atom(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
           (head_is(e, SYM_Rational) || head_is(e, SYM_Complex));
}

static Expr* reverse_rec(Expr* expr, Expr* level_spec, size_t current_level) {
    if (expr->type != EXPR_FUNCTION || is_number_atom(expr)) return expr_copy(expr);

    size_t len = expr->data.function.arg_count;
    Expr** new_args = malloc(sizeof(Expr*) * len);
    bool do_rev = should_reverse_at_level(level_spec, current_level);

    /* An association's entries are its elements, and the next level down is
     * each entry's VALUE, not the Rule: Reverse[<|x -> {1, 2}, y -> {3}|>, 2]
     * is <|x -> {2, 1}, y -> {3}|> (never <|{1, 2} -> x, ...|>). */
    bool assoc = is_association(expr);
    for (size_t i = 0; i < len; i++) {
        size_t src_idx = do_rev ? (len - 1 - i) : i;
        Expr* child = expr->data.function.args[src_idx];
        if (assoc && is_rule2(child))
            new_args[i] = assoc_entry_with_value(child,
                reverse_rec(child->data.function.args[1], level_spec, current_level + 1));
        else
            new_args[i] = reverse_rec(child, level_spec, current_level + 1);
    }

    Expr* result = expr_new_function(expr_copy(expr->data.function.head), new_args, len);
    free(new_args);
    return result;
}

/* The level argument: a positive machine integer or a list of them.  Anything
 * else is Reverse::ilsmp and the call stays unevaluated (audit Reverse#2: 0,
 * -1, x, 1.5 and {1, x} were ignored or half-read). */
static bool reverse_levels_ok(const Expr* ls) {
    if (ls->type == EXPR_INTEGER) return ls->data.integer > 0;
    if (!head_is((Expr*)ls, SYM_List) || ls->data.function.arg_count == 0) return false;
    for (size_t i = 0; i < ls->data.function.arg_count; i++) {
        const Expr* a = ls->data.function.args[i];
        if (a->type != EXPR_INTEGER || a->data.integer <= 0) return false;
    }
    return true;
}

Expr* builtin_reverse(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* expr = res->data.function.args[0];

    if (res->data.function.arg_count == 2 && !reverse_levels_ok(res->data.function.args[1])) {
        char* cs = expr_to_string(res);
        mth_message("Reverse", "ilsmp",
                    "Single or list of positive machine-sized integers expected at position 2 of %s.",
                    cs ? cs : "?");
        free(cs);
        return NULL;
    }

    { bool dh; Expr* dr = ds_native(expr, SYM_Reverse, res, &dh); if (dh) return dr; }
    /* Reverse[ndarray]: block-reverse the leading axis on the buffer (ndstruct.c). */
    if (is_ndarray(expr)) return ndstruct_reverse(res);

    Expr* level_spec = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    /* a Rational/Complex number at the top is an atom: ::normal, unevaluated */
    if (expr_is_number_atom(expr)) return builtin_normal_error("Reverse", res, 1);

    return reverse_rec(expr, level_spec, 1);
}
