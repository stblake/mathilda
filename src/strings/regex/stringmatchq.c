/*
 * stringmatchq.c - StringMatchQ[subject, pattern]
 *
 * Returns True if the WHOLE subject string matches the pattern, else False.
 * The pattern may be RegularExpression["re"], a literal string, or a List of
 * alternatives (matches if any one matches).  A list of subjects threads,
 * giving a list of True/False.
 */

#include "picostrings.h"
#include "regex_common.h"
#include "sym_names.h"

#include <string.h>

/* Whole-string match of one subject against a prebuilt (anchored) rule set. */
static Expr* smq_scalar(const char* subj, RegexRule* rules, int nr) {
    size_t len = strlen(subj);
    size_t ov[2];
    int matched = 0;
    for (int i = 0; i < nr && !matched; i++) {
        if (regex_match(rules[i].prog, subj, len, 0, ov, 1) == 1)
            matched = 1;
    }
    return expr_new_symbol(matched ? SYM_True : SYM_False);
}

Expr* builtin_stringmatchq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;

    Expr* arg0 = res->data.function.args[0];
    Expr* patt = res->data.function.args[1];

    RegexRule* rules;
    int nr = regex_rules_build(patt, /*anchored=*/1, &rules, "StringMatchQ");
    if (nr < 0) return NULL;

    /* StringMatchQ[{s1, s2, ...}, patt] threads over the subject list. */
    if (arg0->type == EXPR_FUNCTION &&
        arg0->data.function.head->type == EXPR_SYMBOL &&
        arg0->data.function.head->data.symbol == SYM_List) {
        size_t n = arg0->data.function.arg_count;
        Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
        if (!out) { regex_rules_free(rules, nr); return NULL; }
        for (size_t i = 0; i < n; i++) {
            Expr* si = arg0->data.function.args[i];
            out[i] = (si->type == EXPR_STRING)
                         ? smq_scalar(si->data.string, rules, nr)
                         : expr_copy(si);   /* non-string element: leave as-is */
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

    Expr* result = smq_scalar(arg0->data.string, rules, nr);
    regex_rules_free(rules, nr);
    return result;
}
