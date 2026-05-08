/* options.c — see options.h. */

#include "options.h"
#include "sym_names.h"
#include <stdbool.h>
#include <string.h>

const Expr* extract_extension_option(const Expr* res, size_t* new_argc) {
    if (!res || res->type != EXPR_FUNCTION) {
        if (new_argc) *new_argc = 0;
        return NULL;
    }
    size_t n = res->data.function.arg_count;
    if (new_argc) *new_argc = n;
    if (n == 0) return NULL;

    const Expr* alpha = NULL;
    /* Walk right-to-left.  Stop at the first non-option argument so we
     * never reorder `Foo[Rule[a,b], poly, Extension -> α]` (the middle
     * `poly` is not an option even though the leftmost arg looks like
     * one).  Multiple trailing Extension rules are all consumed; the
     * rightmost wins. */
    while (n > 0) {
        const Expr* opt = res->data.function.args[n - 1];
        if (opt && opt->type == EXPR_FUNCTION
            && opt->data.function.head
            && opt->data.function.head->type == EXPR_SYMBOL
            && (opt->data.function.head->data.symbol == SYM_Rule
                || opt->data.function.head->data.symbol == SYM_RuleDelayed)
            && opt->data.function.arg_count == 2) {
            const Expr* lhs = opt->data.function.args[0];
            const Expr* rhs = opt->data.function.args[1];
            if (lhs && lhs->type == EXPR_SYMBOL
                && lhs->data.symbol == SYM_Extension) {
                /* Treat None / Automatic as "no extension". */
                bool is_none = (rhs && rhs->type == EXPR_SYMBOL
                                && rhs->data.symbol
                                && (strcmp(rhs->data.symbol, "None") == 0
                                    || strcmp(rhs->data.symbol, "Automatic") == 0));
                if (alpha == NULL && !is_none) {
                    alpha = rhs;
                }
                /* Either way, this option arg is consumed. */
                n--;
                continue;
            }
        }
        break;
    }
    if (new_argc) *new_argc = n;
    return alpha;
}
