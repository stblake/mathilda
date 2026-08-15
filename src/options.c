/* options.c — see options.h. */

#include "options.h"
#include "symtab.h"
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

/* ---- generic named options (see options.h) --------------------------------- */

static bool opt_rule_parts(const Expr* e, const char** name, const Expr** value) {
    if (!e || e->type != EXPR_FUNCTION || !e->data.function.head) return false;
    const Expr* hd = e->data.function.head;
    if (hd->type != EXPR_SYMBOL) return false;
    const char* hn = hd->data.symbol.name;
    if (!hn || (strcmp(hn, "Rule") != 0 && strcmp(hn, "RuleDelayed") != 0)) return false;
    if (e->data.function.arg_count != 2) return false;
    const Expr* lhs = e->data.function.args[0];
    if (!lhs || lhs->type != EXPR_SYMBOL) return false;
    *name = lhs->data.symbol.name;
    *value = e->data.function.args[1];
    return true;
}

bool options_extract(const Expr* res, const char* head_name,
                     const OptEntry* entries, size_t n_entries, size_t* new_argc) {
    if (!res || res->type != EXPR_FUNCTION) return false;
    size_t argc = res->data.function.arg_count;

    /* Find where the trailing options begin: the first argument that is a Rule with a symbol on the
     * left. Scanning from the END rather than the front, so a positional argument that happens to be
     * a rule -- ReplaceAll's second argument, say -- cannot be mistaken for an option. */
    size_t first_opt = argc;
    while (first_opt > 0) {
        const char* nm = NULL; const Expr* vv = NULL;
        if (!opt_rule_parts(res->data.function.args[first_opt - 1], &nm, &vv)) break;
        first_opt--;
    }
    if (new_argc) *new_argc = first_opt;

    /* Apply the call's own options left to right, so the LAST setting wins. */
    for (size_t i = first_opt; i < argc; i++) {
        const char* nm = NULL; const Expr* vv = NULL;
        if (!opt_rule_parts(res->data.function.args[i], &nm, &vv)) return false;
        bool known = false;
        for (size_t j = 0; j < n_entries; j++) {
            if (strcmp(nm, entries[j].name) == 0) {
                if (entries[j].out) *entries[j].out = vv;
                if (entries[j].given) *entries[j].given = true;
                known = true;
                break;
            }
        }
        if (!known) return false;            /* refuse rather than silently ignore */
    }

    /* Anything the call did not set comes from Options[head], which keeps the registered defaults as
     * the single source of truth and makes SetOptions work with no further code here. */
    Expr* defs = head_name ? symtab_get_options(head_name) : NULL;
    if (defs && defs->type == EXPR_FUNCTION) {
        for (size_t j = 0; j < n_entries; j++) {
            if (!entries[j].out || *entries[j].out) continue;
            for (size_t i = 0; i < defs->data.function.arg_count; i++) {
                const char* nm = NULL; const Expr* vv = NULL;
                if (!opt_rule_parts(defs->data.function.args[i], &nm, &vv)) continue;
                if (strcmp(nm, entries[j].name) == 0) { *entries[j].out = vv; break; }
            }
        }
    }
    return true;
}
