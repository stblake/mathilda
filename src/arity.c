/* Mathilda — declared argument counts and option validation.  See arity.h. */

#include "arity.h"

#include "common.h"     /* builtin_arg_error */
#include "message.h"
#include "print.h"      /* expr_to_string */
#include "sym_names.h"

#include <stdlib.h>
#include <string.h>

static void arity_store(const char* name, int min, int max, int validate_from,
                        uint8_t flags) {
    SymbolDef* d = symtab_get_def(name);
    if (min < 0) min = 0;
    if (min > 254) min = 254;
    if (validate_from > 127) validate_from = 127;
    d->arity_min = (uint8_t)min;
    d->arity_max = (max < 0 || max >= (int)ARITY_UNBOUNDED) ? (uint8_t)ARITY_UNBOUNDED
                                                            : (uint8_t)max;
    d->arity_vfrom = (int8_t)(validate_from < 0 ? -1 : validate_from);
    d->arity_flags = flags;
    /* the inline nibbles: lo >= min, hi <= max; an unrepresentable minimum
     * (lo > hi) sends every call to the full check */
    unsigned lo = (unsigned)min, hi = (max < 0 || max > 15) ? 15u : (unsigned)max;
    if (lo > 15) { lo = 15; hi = 0; }
    d->arity_fast = (uint8_t)(lo | (hi << 4));
}

void symtab_set_arity(const char* name, int min, int max) {
    arity_store(name, min, max, -1, ARITY_F_SET);
}

void symtab_set_arity_opts(const char* name, int min, int max, int validate_from) {
    uint8_t f = ARITY_F_SET | ARITY_F_OPTS;
    if (validate_from >= 0) f |= ARITY_F_OPTX;
    arity_store(name, min, max, validate_from, f);
}

void symtab_set_arity_between(const char* name) {
    symtab_get_def(name)->arity_flags |= ARITY_F_BETWEEN;
}

void symtab_clear_arity(const char* name) {
    SymbolDef* d = symtab_get_def(name);
    d->arity_flags = 0;
    d->arity_min = 0;
    d->arity_max = (uint8_t)ARITY_UNBOUNDED;
    d->arity_vfrom = -1;
    d->arity_fast = 0;
}

static bool is_rule(const Expr* e) {
    if (e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    const Expr* h = e->data.function.head;
    return h->type == EXPR_SYMBOL &&
           (h->data.symbol.name == SYM_Rule || h->data.symbol.name == SYM_RuleDelayed);
}

bool arity_is_option_arg(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (is_rule(e)) return true;
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_List) return false;
    size_t n = e->data.function.arg_count;
    if (n == 0) return false;
    for (size_t i = 0; i < n; i++)
        if (!arity_is_option_arg(e->data.function.args[i])) return false;
    return true;
}

/* The name of an option key as the reference spells it in ::optx. */
static const char* key_name(const Expr* k) {
    if (k->type == EXPR_SYMBOL) return k->data.symbol.name;
    if (k->type == EXPR_STRING) return k->data.string;
    return NULL;
}

/* Is the option key `k` among the names of Options[head] (`defs`)?  A symbol
 * and a string of the same name are the same option. */
static bool known_option(const Expr* defs, const Expr* k) {
    const char* kn = key_name(k);
    if (!kn) return false;
    if (!defs || defs->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < defs->data.function.arg_count; i++) {
        const Expr* r = defs->data.function.args[i];
        if (!is_rule(r)) continue;
        const char* dn = key_name(r->data.function.args[0]);
        if (dn && (dn == kn || strcmp(dn, kn) == 0)) return true;
    }
    return false;
}

/* First unknown option key inside the option argument `o`, or NULL. */
static const Expr* first_unknown(const Expr* defs, const Expr* o) {
    if (is_rule(o))
        return known_option(defs, o->data.function.args[0]) ? NULL : o->data.function.args[0];
    for (size_t i = 0; i < o->data.function.arg_count; i++) {
        const Expr* u = first_unknown(defs, o->data.function.args[i]);
        if (u) return u;
    }
    return NULL;
}


/* Does Options[head] (`defs`) name every option inside the option argument
 * `o`?  False for a head without options. */
static bool all_known(const Expr* defs, const Expr* o) {
    if (!defs) return false;
    return first_unknown(defs, o) == NULL;
}

/* Count the rules inside the option argument `o`. */
static size_t option_rule_count(const Expr* o) {
    if (is_rule(o)) return 1;
    size_t k = 0;
    for (size_t i = 0; i < o->data.function.arg_count; i++)
        k += option_rule_count(o->data.function.args[i]);
    return k;
}

/* Append copies of the rules inside `o` to out[*k...]. */
static void option_rules_into(const Expr* o, Expr** out, size_t* k) {
    if (is_rule(o)) { out[(*k)++] = expr_copy((Expr*)o); return; }
    for (size_t i = 0; i < o->data.function.arg_count; i++)
        option_rules_into(o->data.function.args[i], out, k);
}

/* A list of options is the same options spliced in: f[x, {a -> 1, b -> 2}]
 * is f[x, a -> 1, b -> 2].  Rewrites the option arguments at index >= from
 * of `call` (owned by the evaluator, unshared) in place, so the head's C
 * code sees only rules. */
static void flatten_option_lists(Expr* call, size_t from) {
    size_t n = call->data.function.arg_count;
    Expr** args = call->data.function.args;
    size_t total = from, lists = 0;
    for (size_t i = from; i < n; i++) {
        if (is_rule(args[i])) { total++; continue; }
        lists++;
        total += option_rule_count(args[i]);
    }
    if (!lists) return;
    Expr** na = malloc(sizeof(Expr*) * (total ? total : 1));
    size_t k = 0;
    for (size_t i = 0; i < from; i++) na[k++] = args[i];
    for (size_t i = from; i < n; i++) {
        if (is_rule(args[i])) { na[k++] = args[i]; continue; }
        option_rules_into(args[i], na, &k);
        expr_free(args[i]);
    }
    free(args);
    call->data.function.args = na;
    call->data.function.arg_count = k;
    expr_invalidate_hash(call);
}

/* The argument-count message.  A head flagged ARITY_F_BETWEEN words a choice
 * of two counts as a range, as the reference does for For and While
 * ("between 1 and 2 arguments are expected"). */
static void count_error(const SymbolDef* d, size_t got, size_t min, size_t max) {
    const char* head = d->symbol_name;
    if ((d->arity_flags & ARITY_F_BETWEEN) && max == min + 1) {
        mth_message(head, (got == 1) ? "argbu" : "argb",
                    "%s called with %zu argument%s; between %zu and %zu arguments are expected.",
                    head, got, (got == 1) ? "" : "s", min, max);
        return;
    }
    builtin_arg_error(head, got, min, max);
}

bool arity_check_slow(const SymbolDef* d, Expr* call) {
    const char* head = d->symbol_name;
    size_t n = call->data.function.arg_count;
    Expr* const* args = call->data.function.args;
    size_t min = d->arity_min;
    size_t max = d->arity_max == ARITY_UNBOUNDED ? SIZE_MAX : d->arity_max;
    const Expr* defs = d->default_options;
    if (defs && (defs->type != EXPR_FUNCTION || defs->data.function.arg_count == 0))
        defs = NULL;

    if (!(d->arity_flags & ARITY_F_OPTS) && !defs) {
        if (n < min || n > max) {
            count_error(d, n, min, max);
            return true;
        }
        return false;
    }

    /* Trailing options.  One in the first `min` slots is positional unless
     * every name in it is a known option of the head. */
    size_t p = n;
    while (p > 0 && arity_is_option_arg(args[p - 1]) &&
           (p - 1 >= min || all_known(defs, args[p - 1])))
        p--;
    if (p < min) {
        count_error(d, p, min, max);
        return true;
    }
    if (p > max) {
        /* the first non-option argument beyond the last positional slot */
        size_t bad = max;
        while (bad < n && arity_is_option_arg(args[bad])) bad++;
        if (bad >= n) bad = max;
        char* bs = expr_to_string(args[bad]);
        char* cs = expr_to_string((Expr*)call);
        mth_message(head, "nonopt",
                    "Options expected (instead of %s) beyond position %zu in %s. "
                    "An option must be a rule or a list of rules.",
                    bs ? bs : "", max, cs ? cs : "");
        free(bs); free(cs);
        return true;
    }
    if (d->arity_flags & ARITY_F_OPTX) {
        size_t from = (size_t)(d->arity_vfrom < 0 ? 0 : d->arity_vfrom);
        if (from < p) from = p;
        for (size_t i = from; i < n; i++) {
            const Expr* u = first_unknown(defs, args[i]);   /* no Options: none known */
            if (!u) continue;
            char* us = expr_to_string((Expr*)u);
            char* cs = expr_to_string((Expr*)call);
            mth_message(head, "optx", "Unknown option %s in %s.",
                        u->type == EXPR_STRING ? u->data.string : (us ? us : ""),
                        cs ? cs : "");
            free(us); free(cs);
            return true;
        }
        flatten_option_lists(call, from);
    }
    return false;
}
