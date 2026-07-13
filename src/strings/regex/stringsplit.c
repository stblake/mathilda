/*
 * stringsplit.c - StringSplit[subject, delimiter]
 *
 * Splits `subject` into the substrings between non-overlapping matches of the
 * delimiter pattern.  Empty pieces are dropped (Wolfram default).  The
 * delimiter may be RegularExpression["re"], a literal string, or a List of
 * alternative delimiters; a list of subjects threads.
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

static Expr* ss_scalar(const char* subj, RegexRule* rules, int nr) {
    size_t len = strlen(subj);
    Expr** items = NULL;
    size_t count = 0, cap = 0;

    size_t last_end = 0, pos = 0;
    while (pos <= len) {
        int found = 0;
        size_t best_ms = 0, best_me = 0;
        for (int i = 0; i < nr; i++) {
            size_t ov[2];
            if (regex_match(rules[i].prog, subj, len, pos, ov, 1) == 1) {
                if (!found || ov[0] < best_ms) {
                    found = 1;
                    best_ms = ov[0];
                    best_me = ov[1];
                }
            }
        }
        if (!found) break;

        if (best_ms > last_end) {                 /* drop empty pieces */
            if (count == cap) {
                cap = cap ? cap * 2 : 8;
                Expr** ni = realloc(items, cap * sizeof(Expr*));
                if (!ni) break;
                items = ni;
            }
            items[count++] = substr_expr(subj, last_end, best_ms);
        }
        last_end = best_me;
        pos = (best_me > best_ms) ? best_me : best_ms + 1;
    }

    if (len > last_end) {                          /* trailing piece */
        if (count == cap) {
            cap = cap ? cap * 2 : 8;
            Expr** ni = realloc(items, cap * sizeof(Expr*));
            if (ni) items = ni; else len = last_end;   /* skip on OOM */
        }
        if (len > last_end)
            items[count++] = substr_expr(subj, last_end, len);
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), items, count);
    free(items);
    return result;
}

Expr* builtin_stringsplit(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;

    Expr* arg0 = res->data.function.args[0];
    Expr* patt = res->data.function.args[1];

    RegexRule* rules;
    int nr = regex_rules_build(patt, /*anchored=*/0, &rules, "StringSplit");
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
                         ? ss_scalar(si->data.string, rules, nr)
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

    Expr* result = ss_scalar(arg0->data.string, rules, nr);
    regex_rules_free(rules, nr);
    return result;
}
