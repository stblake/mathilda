/*
 * stringcases.c - StringCases[subject, pattern]
 *
 * Returns a List of the substrings of `subject` that match `pattern`, left to
 * right.  With a rule pattern (patt -> rhs / patt :> rhs) each match is replaced
 * by the rhs, with $0/$1... expanded to the whole match and capture groups.  The
 * pattern may also be a List of alternatives/rules; at each position the
 * leftmost match wins, ties broken by rule order.  A list of subjects threads.
 *
 * Options:
 *   Overlaps -> False (default) | True | All
 *     False - non-overlapping, greedy left-to-right.
 *     True  - overlapping substrings count separately, but only the first match
 *             starting at each position.
 *     All   - every matching substring at every start (all lengths).
 *   IgnoreCase -> True | False (default)
 *     Treat upper/lowercase as equivalent.
 *
 * The match enumeration itself is regex_scan() in regex_common.c, shared with
 * StringCount and StringPosition, so StringCount[s, p, opts] always equals
 * Length[StringCases[s, p, opts]].
 */

#include "picostrings.h"
#include "regex_common.h"
#include "sym_names.h"
#include "symtab.h"
#include "common.h"

#include <stdlib.h>
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

/*
 * Collect the matches for one subject string.  `want_captures` is set by the
 * caller when any rule carries a replacement RHS, so a pure-extraction call
 * never pays for the capture pool.  Returns a fresh List[...], or NULL on OOM.
 */
static Expr* sc_scalar(const char* subj, RegexRule* rules, int nr,
                       RegexOverlapMode mode, int want_captures) {
    RegexScan scan;
    if (regex_scan(subj, strlen(subj), rules, nr, mode, want_captures, &scan) < 0)
        return NULL;

    Expr** items = malloc(sizeof(Expr*) * (scan.count ? scan.count : 1));
    if (!items) { regex_scan_free(&scan); return NULL; }

    for (size_t k = 0; k < scan.count; k++) {
        RegexSpan* sp = &scan.spans[k];
        Expr* item = NULL;
        if (rules[sp->rule].rhs && scan.caps) {
            char* rep = regex_rule_replacement(&rules[sp->rule], subj,
                                               scan.caps + sp->caps_off,
                                               sp->npairs);
            if (rep) { item = expr_new_string(rep); free(rep); }
        }
        /* Bare pattern, or a rule whose RHS we cannot expand: the match itself. */
        items[k] = item ? item : substr_expr(subj, sp->ms, sp->me);
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), items, scan.count);
    free(items);
    regex_scan_free(&scan);
    return result;
}

Expr* builtin_stringcases(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;

    if (argc == 0) return builtin_arg_error("StringCases", 0, 2, 2);

    /* Seed option state from the registered defaults (so SetOptions[StringCases,
     * ...] takes effect), then let explicit trailing options override below.
     * Defaults are {IgnoreCase -> False, Overlaps -> False}. */
    int caseless = 0;
    RegexOverlapMode mode = REGEX_OV_FALSE;
    Expr* defs = symtab_get_options("StringCases");   /* borrowed */
    if (defs && defs->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < defs->data.function.arg_count; i++) {
            int v;
            if (regex_match_opt(defs->data.function.args[i], SYM_IgnoreCase, &v, 0))
                caseless = v;
            else if (regex_match_opt(defs->data.function.args[i], SYM_Overlaps, &v, 1))
                mode = (RegexOverlapMode)v;
        }
    }

    /* Strip trailing IgnoreCase / Overlaps options, leaving positional args. */
    size_t pargc = argc;
    while (pargc >= 2) {
        int v;
        if (regex_match_opt(a[pargc - 1], SYM_IgnoreCase, &v, 0)) {
            caseless = v; pargc--;
        } else if (regex_match_opt(a[pargc - 1], SYM_Overlaps, &v, 1)) {
            mode = (RegexOverlapMode)v; pargc--;
        } else {
            break;
        }
    }
    /* Too few positional args is a genuine arity error. Extra ones are not: the
     * WL occurrence-limit form StringCases["s", patt, n] is simply unsupported
     * here, and claiming "2 arguments are expected" would misdescribe it, so it
     * is left unevaluated silently (as StringCases has always done). */
    if (pargc < 2) return builtin_arg_error("StringCases", argc, 2, 2);
    if (pargc > 2) return NULL;

    Expr* subject = a[0];
    Expr* patt = a[1];

    /* All-mode enumerates exact-substring matches, so build anchored rules. */
    int anchored = (mode == REGEX_OV_ALL) ? 1 : 0;
    RegexRule* rules;
    int nr = regex_rules_build_ex(patt, anchored, caseless, &rules, "StringCases");
    if (nr < 0) return NULL;

    /* Capture offsets are only needed to expand a rule's $n template. */
    int want_captures = 0;
    for (int i = 0; i < nr; i++) if (rules[i].rhs) { want_captures = 1; break; }

    Expr* result;
    if (subject->type == EXPR_FUNCTION &&
        subject->data.function.head->type == EXPR_SYMBOL &&
        subject->data.function.head->data.symbol.name == SYM_List) {
        size_t m = subject->data.function.arg_count;
        Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
        if (!out) { regex_rules_free(rules, nr); return NULL; }
        size_t built = 0;
        for (; built < m; built++) {
            Expr* si = subject->data.function.args[built];
            Expr* one = (si->type == EXPR_STRING)
                            ? sc_scalar(si->data.string, rules, nr, mode, want_captures)
                            : expr_copy(si);
            if (!one) break;                     /* OOM: unwind and bail out */
            out[built] = one;
        }
        if (built < m) {                         /* partial: free and leave unevaluated */
            for (size_t k = 0; k < built; k++) expr_free(out[k]);
            free(out);
            regex_rules_free(rules, nr);
            return NULL;
        }
        result = expr_new_function(expr_new_symbol(SYM_List), out, m);
        free(out);
    } else if (subject->type == EXPR_STRING) {
        result = sc_scalar(subject->data.string, rules, nr, mode, want_captures);
    } else {
        result = NULL;   /* non-string subject: leave unevaluated */
    }

    regex_rules_free(rules, nr);
    return result;
}
