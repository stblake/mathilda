/*
 * stringreplace.c - StringReplace[subject, rule | {rules...}]
 *
 * Replaces every non-overlapping match of a rule's pattern by its right-hand
 * side, scanning left to right; unmatched text is copied verbatim.  The RHS is
 * a string in which $0/$1... expand to the whole match and capture groups.
 * With several rules, at each position the leftmost match wins, ties broken by
 * rule order.  A list of subjects threads.
 */

#include "picostrings.h"
#include "regex_common.h"
#include "sym_names.h"

#include <string.h>

/* Assemble the replaced text for one subject. Returns a malloc'd string. */
static char* sr_scalar_str(const char* subj, RegexRule* rules, int nr) {
    size_t len = strlen(subj);
    RegexBuf b = {0};

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
        if (!found) {
            regexbuf_add(&b, subj + pos, len - pos);   /* tail */
            break;
        }

        regexbuf_add(&b, subj + pos, best_ms - pos);   /* literal before match */

        char* rep = regex_rule_replacement(&rules[best_rule], subj,
                                           best_ov, best_pairs);
        if (rep) { regexbuf_add(&b, rep, strlen(rep)); free(rep); }
        else     { regexbuf_add(&b, subj + best_ms, best_me - best_ms); }

        if (best_me > best_ms) {
            pos = best_me;
        } else {                                       /* zero-width match */
            if (best_ms < len) regexbuf_add(&b, subj + best_ms, 1);
            pos = best_ms + 1;
        }
    }

    if (!b.p) { b.p = malloc(1); if (b.p) b.p[0] = '\0'; }
    return b.p;
}

static Expr* sr_scalar(const char* subj, RegexRule* rules, int nr) {
    char* s = sr_scalar_str(subj, rules, nr);
    if (!s) return expr_new_string("");
    Expr* r = expr_new_string(s);
    free(s);
    return r;
}

Expr* builtin_stringreplace(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;

    Expr* arg0 = res->data.function.args[0];
    Expr* patt = res->data.function.args[1];

    RegexRule* rules;
    int nr = regex_rules_build(patt, /*anchored=*/0, &rules, "StringReplace");
    if (nr < 0) return NULL;

    /* StringReplace requires rule(s): every element must carry an RHS. */
    for (int i = 0; i < nr; i++) {
        if (!rules[i].rhs) { regex_rules_free(rules, nr); return NULL; }
    }

    if (arg0->type == EXPR_FUNCTION &&
        arg0->data.function.head->type == EXPR_SYMBOL &&
        arg0->data.function.head->data.symbol == SYM_List) {
        size_t n = arg0->data.function.arg_count;
        Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
        if (!out) { regex_rules_free(rules, nr); return NULL; }
        for (size_t i = 0; i < n; i++) {
            Expr* si = arg0->data.function.args[i];
            out[i] = (si->type == EXPR_STRING)
                         ? sr_scalar(si->data.string, rules, nr)
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

    Expr* result = sr_scalar(arg0->data.string, rules, nr);
    regex_rules_free(rules, nr);
    return result;
}
