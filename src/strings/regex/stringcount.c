/*
 * stringcount.c - StringCount[subject, pattern]
 *
 * Counts the substrings of `subject` that match `pattern`.  This is the
 * counting-only companion of StringCases: it runs the same match enumeration
 * (regex_scan in regex_common.c) but never materialises a substring, so it costs
 * one small span record per match instead of a malloc + Expr + List element.
 *
 * The pattern may be a literal string, a general string expression
 * (Blank/Pattern/~~/RegularExpression/character classes), a Rule or RuleDelayed
 * (only the LHS matters -- the replacement is irrelevant to a count), or a List
 * of any of those.  A List of subjects threads, giving one count per subject.
 *
 * Options:
 *   Overlaps -> False (default) | True | All
 *     False - overlapping substrings are not counted separately.
 *     True  - overlapping substrings count separately, but only the first
 *             matching substring at a given position is counted.
 *     All   - every matching substring at every position is counted separately.
 *   IgnoreCase -> True | False (default)
 *     Treat upper/lowercase as equivalent.
 *
 * Because the scan is shared, StringCount[s, p, opts] is always exactly
 * Length[StringCases[s, p, opts]] and Length[StringPosition[s, p, opts]].
 *
 * Byte semantics: like the rest of src/strings, offsets are byte offsets (no
 * UTF-8 codepoint decoding), consistent with StringLength / StringPart.
 */

#include "picostrings.h"
#include "regex_common.h"
#include "sym_names.h"
#include "symtab.h"
#include "common.h"

#include <stdlib.h>
#include <string.h>

/* Count the matches in one subject string. Returns a fresh Integer, or NULL on
 * OOM (which leaves the whole call unevaluated rather than reporting a wrong
 * count). */
static Expr* sct_scalar(const char* subj, RegexRule* rules, int nr,
                        RegexOverlapMode mode) {
    RegexScan scan;
    long n = regex_scan(subj, strlen(subj), rules, nr, mode,
                        /*want_captures=*/0, &scan);
    if (n < 0) return NULL;
    regex_scan_free(&scan);
    return expr_new_integer((int64_t)n);
}

Expr* builtin_stringcount(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;

    if (argc == 0) return builtin_arg_error("StringCount", 0, 2, 2);

    /* Seed option state from the registered defaults (so SetOptions[StringCount,
     * ...] takes effect), then let explicit trailing options override below.
     * Defaults are {IgnoreCase -> False, Overlaps -> False}. */
    int caseless = 0;
    RegexOverlapMode mode = REGEX_OV_FALSE;
    Expr* defs = symtab_get_options("StringCount");   /* borrowed */
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
    if (pargc != 2) return builtin_arg_error("StringCount", argc, 2, 2);

    Expr* subject = a[0];
    Expr* patt = a[1];

    /* All-mode enumerates exact-substring matches, so build anchored rules. */
    int anchored = (mode == REGEX_OV_ALL) ? 1 : 0;
    RegexRule* rules;
    int nr = regex_rules_build_ex(patt, anchored, caseless, &rules, "StringCount");
    if (nr < 0) return NULL;

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
                            ? sct_scalar(si->data.string, rules, nr, mode)
                            : expr_copy(si);   /* non-string element passes through */
            if (!one) break;                   /* OOM: unwind and bail out */
            out[built] = one;
        }
        if (built < m) {                       /* partial: free and leave unevaluated */
            for (size_t k = 0; k < built; k++) expr_free(out[k]);
            free(out);
            regex_rules_free(rules, nr);
            return NULL;
        }
        result = expr_new_function(expr_new_symbol(SYM_List), out, m);
        free(out);
    } else if (subject->type == EXPR_STRING) {
        result = sct_scalar(subject->data.string, rules, nr, mode);
    } else {
        result = NULL;   /* non-string subject: leave unevaluated */
    }

    regex_rules_free(rules, nr);
    return result;
}
