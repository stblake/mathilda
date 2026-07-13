/* options.c — see options.h. */

#include "options.h"
#include "sym_names.h"
#include <stdbool.h>
#include <string.h>

const Expr* extract_extension_option_full(const Expr* res, size_t* new_argc,
                                          bool* automatic_out) {
    if (automatic_out) *automatic_out = false;
    if (!res || res->type != EXPR_FUNCTION) {
        if (new_argc) *new_argc = 0;
        return NULL;
    }
    size_t n = res->data.function.arg_count;
    if (new_argc) *new_argc = n;
    if (n == 0) return NULL;

    const Expr* alpha = NULL;
    /* `seen_*` flags track whether we've already encountered an
     * Extension rule of that flavour to the *right* of the current
     * scan position.  Mathematica semantics: rightmost wins.  We're
     * walking right-to-left so the FIRST rule we see (most rightmost)
     * is authoritative; everything to its left of the same kind is
     * still consumed but doesn't override the rightmost value. */
    bool seen_rule = false;        /* any Extension rule yet? */
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
            && (opt->data.function.head->data.symbol.name == SYM_Rule
                || opt->data.function.head->data.symbol.name == SYM_RuleDelayed)
            && opt->data.function.arg_count == 2) {
            const Expr* lhs = opt->data.function.args[0];
            const Expr* rhs = opt->data.function.args[1];
            if (lhs && lhs->type == EXPR_SYMBOL
                && lhs->data.symbol.name == SYM_Extension) {
                bool is_none      = (rhs && rhs->type == EXPR_SYMBOL
                                     && rhs->data.symbol.name
                                     && strcmp(rhs->data.symbol.name, "None") == 0);
                bool is_automatic = (rhs && rhs->type == EXPR_SYMBOL
                                     && rhs->data.symbol.name
                                     && strcmp(rhs->data.symbol.name, "Automatic") == 0);
                if (!seen_rule) {
                    /* Rightmost setting wins. */
                    if (is_automatic) {
                        if (automatic_out) *automatic_out = true;
                        alpha = NULL;
                    } else if (is_none) {
                        alpha = NULL;
                    } else {
                        alpha = rhs;
                    }
                    seen_rule = true;
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

const Expr* extract_extension_option(const Expr* res, size_t* new_argc) {
    return extract_extension_option_full(res, new_argc, NULL);
}
