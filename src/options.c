/* options.c — see options.h. */

#include "options.h"
#include "symtab.h"
#include "sym_names.h"
#include "message.h"
#include "print.h"
#include <stdbool.h>
#include <stdlib.h>
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
                                     && rhs->data.symbol.name == SYM_None);
                bool is_automatic = (rhs && rhs->type == EXPR_SYMBOL
                                     && rhs->data.symbol.name
                                     && rhs->data.symbol.name == SYM_Automatic);
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
    if (!hn || (hn != SYM_Rule && hn != SYM_RuleDelayed)) return false;
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

/* ---- validated Boolean options (see options.h) ----------------------------- */

/* Split a Rule/RuleDelayed argument into its key and value; false when `e` is
 * not a two-argument rule. */
static bool opt_bool_rule(const Expr* e, const Expr** key, const Expr** val) {
    if (!e || e->type != EXPR_FUNCTION || !e->data.function.head) return false;
    const Expr* hd = e->data.function.head;
    if (hd->type != EXPR_SYMBOL) return false;
    if (hd->data.symbol.name != SYM_Rule && hd->data.symbol.name != SYM_RuleDelayed)
        return false;
    if (e->data.function.arg_count != 2) return false;
    *key = e->data.function.args[0];
    *val = e->data.function.args[1];
    return true;
}

/* Index of `key` in names[], or -1. */
static int opt_bool_index(const Expr* key, const char* const* names, size_t n) {
    for (size_t j = 0; j < n; j++) {
        if (key->type == EXPR_SYMBOL && key->data.symbol.name == names[j]) return (int)j;
        if (key->type == EXPR_STRING && strcmp(key->data.string, names[j]) == 0) return (int)j;
    }
    return -1;
}

bool options_parse_bools(Expr* res, const char* head, size_t max_pos,
                         const char* const* names, bool* vals, size_t n_names,
                         size_t* npos) {
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;
    const Expr *key = NULL, *val = NULL;

    size_t first = argc;
    while (first > 0 && opt_bool_rule(args[first - 1], &key, &val)) first--;
    *npos = first;
    if (first > max_pos) {
        char* bad = expr_to_string(args[first - 1]);
        char* call = expr_to_string(res);
        mth_message(head, "nonopt",
                    "Options expected (instead of %s) beyond position %zu in %s. "
                    "An option must be a rule or a list of rules.",
                    bad ? bad : "?", max_pos, call ? call : "?");
        free(bad); free(call);
        return false;
    }

    /* Defaults from Options[head], so SetOptions takes effect. */
    Expr* defs = symtab_get_options(head);
    if (defs && defs->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < defs->data.function.arg_count; i++) {
            if (!opt_bool_rule(defs->data.function.args[i], &key, &val)) continue;
            int j = opt_bool_index(key, names, n_names);
            if (j >= 0 && val->type == EXPR_SYMBOL)
                vals[j] = (val->data.symbol.name == SYM_True);
        }
    }

    bool seen[8] = { false };             /* n_names <= 8 */
    for (size_t i = first; i < argc; i++) {
        opt_bool_rule(args[i], &key, &val);
        int j = opt_bool_index(key, names, n_names);
        if (j < 0) {
            char* kname = key->type == EXPR_STRING ? NULL : expr_to_string((Expr*)key);
            char* call = expr_to_string(res);
            mth_message(head, "optx", "Unknown option %s in %s.",
                        key->type == EXPR_STRING ? key->data.string : (kname ? kname : "?"),
                        call ? call : "?");
            free(kname); free(call);
            return false;
        }
        const char* v = val->type == EXPR_SYMBOL ? val->data.symbol.name : NULL;
        if (v != SYM_True && v != SYM_False) {
            char* bad = expr_to_string(args[i]);
            mth_message(head, "opttf", "Value of option %s should be True or False.",
                        bad ? bad : "?");
            free(bad);
            return false;
        }
        if ((size_t)j < 8 && !seen[j]) { vals[j] = (v == SYM_True); seen[j] = true; }
    }
    return true;
}
