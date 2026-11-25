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
#include "linkedlist.h"   /* a LinkedList of subjects */
#include "symtab.h"
#include "sym_names.h"
#include "strings_internal.h"   /* strings_msg_strse */

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
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;

    /* IgnoreCase: registered default, overridden by a trailing option. */
    int caseless = 0;
    Expr* defs = symtab_get_options("StringMatchQ");   /* borrowed */
    if (defs && defs->type == EXPR_FUNCTION)
        for (size_t i = 0; i < defs->data.function.arg_count; i++)
            regex_match_opt(defs->data.function.args[i], SYM_IgnoreCase, &caseless, 0);
    for (;;) {
        int spell;
        if (argc > 2 && regex_match_opt(a[argc - 1], SYM_IgnoreCase, &caseless, 0)) argc--;
        else if (argc > 2 && regex_match_opt(a[argc - 1], SYM_SpellingCorrection, &spell, 0)) argc--;
        else break;
    }
    if (argc != 2) return NULL;

    Expr* arg0 = a[0];
    Expr* patt = a[1];

    /* A string or a list of strings; 5, {"a", 5} and {{"a"}} are
     * StringMatchQ::strse and stay unevaluated (a Q-function never answers
     * with a non-Boolean element). */
    if (arg0->type != EXPR_STRING &&
        !(strings_is_seq(arg0) && strings_seq_all_strings(arg0)))
        return strings_msg_strse("StringMatchQ", res, 1);

    RegexRule* rules;
    int nr = regex_rules_build_opts(patt, REGEX_RULE_ANCHORED | REGEX_RULE_WILDCARDS,
                                    caseless ? REGEX_CASELESS : 0u,
                                    &rules, "StringMatchQ");
    if (nr < 0) return NULL;

    /* StringMatchQ[{s1, s2, ...}, patt] threads over the subject list. */
    LLReader sv;
    if (ll_rd_open(arg0, &sv)) {          /* a List or a LinkedList of subjects */
        size_t n = sv.n;
        Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
        if (!out) { ll_rd_close(&sv); regex_rules_free(rules, nr); return NULL; }
        for (size_t i = 0; i < n; i++) {
            Expr* si = ll_rd_at(&sv, i);
            out[i] = (si->type == EXPR_STRING)
                         ? smq_scalar(si->data.string, rules, nr)
                         : expr_copy(si);   /* non-string element: leave as-is */
        }
        Expr* result = ll_seq_like(arg0, out, n);
        free(out);
        ll_rd_close(&sv);
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
