/*
 * stringcases.c - StringCases[subject, pattern]
 *
 * Returns a List of the non-overlapping substrings of `subject` that match
 * `pattern`, left to right.  With a rule pattern (patt -> rhs / patt :> rhs)
 * each match is replaced by the rhs, with $0/$1... expanded to the whole match
 * and capture groups.  The pattern may also be a List of alternatives/rules;
 * at each position the leftmost match wins, ties broken by rule order.  A list
 * of subjects threads.
 */

#include "picostrings.h"
#include "regex_common.h"
#include "sym_names.h"

#include <string.h>

/* New EXPR_STRING from subj[s..e). */
static Expr* substr_expr(const char* subj, size_t s, size_t e) {
    size_t n = e - s;
    char* buf = malloc(n + 1);
    if (!buf) return expr_new_string("");
    memcpy(buf, subj + s, n);
    buf[n] = '\0';
    Expr* r = expr_new_string(buf);
    free(buf);
    return r;
}

static Expr* sc_scalar(const char* subj, RegexRule* rules, int nr) {
    size_t len = strlen(subj);
    Expr** items = NULL;
    size_t count = 0, cap = 0;

    size_t pos = 0;
    while (pos <= len) {
        int found = 0, best_rule = 0;
        size_t best_ms = 0, best_me = 0, best_pairs = 0;
        size_t best_ov[REGEX_MAX_PAIRS * 2];

        for (int i = 0; i < nr; i++) {
            int gc = regex_group_count(rules[i].prog) + 1;
            size_t pairs = (gc > REGEX_MAX_PAIRS) ? REGEX_MAX_PAIRS : (size_t)gc;
            size_t ov[REGEX_MAX_PAIRS * 2];
            if (regex_match(rules[i].prog, subj, len, pos, ov, pairs) == 1) {
                if (!found || ov[0] < best_ms) {
                    found = 1;
                    best_rule = i;
                    best_ms = ov[0];
                    best_me = ov[1];
                    best_pairs = pairs;
                    memcpy(best_ov, ov, sizeof(size_t) * 2 * pairs);
                }
            }
        }
        if (!found) break;

        Expr* item;
        if (rules[best_rule].rhs) {
            char* rep = regex_rule_replacement(&rules[best_rule], subj,
                                               best_ov, best_pairs);
            if (rep) { item = expr_new_string(rep); free(rep); }
            else     { item = substr_expr(subj, best_ms, best_me); }
        } else {
            item = substr_expr(subj, best_ms, best_me);
        }

        if (count == cap) {
            cap = cap ? cap * 2 : 8;
            Expr** ni = realloc(items, cap * sizeof(Expr*));
            if (!ni) { expr_free(item); break; }
            items = ni;
        }
        items[count++] = item;

        pos = (best_me > best_ms) ? best_me : best_ms + 1;   /* progress */
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), items, count);
    free(items);
    return result;
}

Expr* builtin_stringcases(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;

    Expr* arg0 = res->data.function.args[0];
    Expr* patt = res->data.function.args[1];

    RegexRule* rules;
    int nr = regex_rules_build(patt, /*anchored=*/0, &rules, "StringCases");
    if (nr < 0) return NULL;

    if (arg0->type == EXPR_FUNCTION &&
        arg0->data.function.head->type == EXPR_SYMBOL &&
        arg0->data.function.head->data.symbol.name == SYM_List) {
        size_t n = arg0->data.function.arg_count;
        Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
        if (!out) { regex_rules_free(rules, nr); return NULL; }
        for (size_t i = 0; i < n; i++) {
            Expr* si = arg0->data.function.args[i];
            out[i] = (si->type == EXPR_STRING)
                         ? sc_scalar(si->data.string, rules, nr)
                         : expr_copy(si);
        }
        Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, n);
        free(out);
        regex_rules_free(rules, nr);
        return result;
    }

    if (arg0->type != EXPR_STRING) {
        regex_rules_free(rules, nr);
        return NULL;
    }

    Expr* result = sc_scalar(arg0->data.string, rules, nr);
    regex_rules_free(rules, nr);
    return result;
}
