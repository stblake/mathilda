/* Pick[expr, sel] / Pick[expr, sel, patt] — select elements of `expr` whose
 * positionally-corresponding element of `sel` matches `patt` (literal `True`
 * for the two-argument form).
 *
 * The selector array mirrors the structure of `expr`, so the walk is a
 * simultaneous recursive descent over both trees. At each element:
 *
 *   - selector matches `patt`            -> keep the whole `expr` element,
 *   - selector is compound and no match  -> recurse into the pair,
 *   - selector is atomic and no match    -> drop the element.
 *
 * The head at every level comes from `expr`, never from `sel`, so
 * Pick[f[a, b, c], {True, False, True}] is f[a, c].
 *
 * Any structural disagreement between the two trees (differing lengths, or a
 * compound selector against an atomic expression element) is not an error:
 * the builtin returns NULL and the evaluator leaves `Pick[...]` unevaluated,
 * matching Mathematica's Pick::incomp behaviour of returning the input.
 * Detection is exact — a mismatch found arbitrarily deep aborts the whole
 * call rather than yielding a partially picked result. */

#include "list_common.h"
#include "match.h"
#include "pick.h"

/* True when `sel` matches `patt`. Bindings are irrelevant here (Pick has no
 * right-hand side to substitute into), so the environment is discarded. */
static bool pick_selects(Expr* sel, Expr* patt) {
    MatchEnv* env = env_new();
    bool matched = match(sel, patt, env);
    env_free(env);
    return matched;
}

/* Pick one level. On structure mismatch sets *ok to false and returns NULL;
 * callers must check *ok before using the result. */
static Expr* pick_rec(Expr* expr, Expr* sel, Expr* patt, bool* ok) {
    if (expr->type != EXPR_FUNCTION || sel->type != EXPR_FUNCTION ||
        expr->data.function.arg_count != sel->data.function.arg_count) {
        *ok = false;
        return NULL;
    }

    size_t n = expr->data.function.arg_count;
    Expr** kept = NULL;
    if (n > 0) {
        kept = malloc(sizeof(Expr*) * n);
        if (!kept) {
            *ok = false;
            return NULL;
        }
    }

    size_t kept_count = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* e_i = expr->data.function.args[i];
        Expr* s_i = sel->data.function.args[i];

        if (pick_selects(s_i, patt)) {
            kept[kept_count++] = expr_copy(e_i);
        } else if (s_i->type == EXPR_FUNCTION) {
            /* Compound selector that did not match as a whole: descend. The
             * corresponding expression element must be compound and of equal
             * length, which pick_rec verifies. */
            Expr* sub = pick_rec(e_i, s_i, patt, ok);
            if (!*ok) {
                for (size_t j = 0; j < kept_count; j++) expr_free(kept[j]);
                free(kept);
                return NULL;
            }
            kept[kept_count++] = sub;
        }
        /* else: atomic selector, no match -> element is simply dropped. */
    }

    Expr* result = expr_new_function(expr_copy(expr->data.function.head),
                                     kept, kept_count);
    free(kept);
    return result;
}

Expr* builtin_pick(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;

    Expr* expr = res->data.function.args[0];
    Expr* sel  = res->data.function.args[1];

    /* Two-argument form selects on literal True only; anything else (False,
     * 1, a symbol, ...) simply fails to match and is dropped. */
    Expr* implicit_true = (argc == 2) ? expr_new_symbol(SYM_True) : NULL;
    Expr* patt = implicit_true ? implicit_true : res->data.function.args[2];

    bool ok = true;
    Expr* result = pick_rec(expr, sel, patt, &ok);

    if (implicit_true) expr_free(implicit_true);
    if (!ok) {
        if (result) expr_free(result);
        return NULL;
    }
    return result;
}
