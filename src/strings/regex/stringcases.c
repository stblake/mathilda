/*
 * stringcases.c - StringCases[subject, pattern]
 *
 * Returns a List of the substrings of `subject` that match `pattern`, left to
 * right.  With a rule pattern (patt -> rhs / patt :> rhs) each match is replaced
 * by the rhs with the pattern's names bound to the pieces they matched
 * (StringCases["abc123", x__ ~~ "123" :> x] = {"abc"}); for a RegularExpression
 * pattern a string rhs may use $0/$1... for the whole match and capture groups.  The
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
#include "linkedlist.h"   /* a LinkedList of subjects */
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
                       RegexOverlapMode mode, int want_captures, size_t limit) {
    RegexScan scan;
    if (regex_scan_limit(subj, strlen(subj), rules, nr, mode, want_captures,
                         limit, &scan) < 0)
        return NULL;

    Expr** items = malloc(sizeof(Expr*) * (scan.count ? scan.count : 1));
    if (!items) { regex_scan_free(&scan); return NULL; }

    for (size_t k = 0; k < scan.count; k++) {
        RegexSpan* sp = &scan.spans[k];
        /* A rule: its RHS instantiated with this match's named pieces. A bare
         * pattern: the matched text itself. */
        Expr* item = rules[sp->rule].rhs
            ? regex_rule_eval_rhs(&rules[sp->rule], subj,
                                  scan.caps ? scan.caps + sp->caps_off : NULL,
                                  scan.caps ? sp->npairs : 0)
            : NULL;
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

    if (argc == 0) return builtin_arg_error("StringCases", 0, 1, 3);

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
    /* StringCases[patt] is the operator form (applied by src/opform.c), so a
     * single positional argument stays unevaluated. StringCases[s, patt, n]
     * keeps the first n matches (n an integer >= 0, or Infinity). */
    if (pargc == 1 || pargc > 3) return NULL;
    size_t limit = 0;
    if (pargc == 3) {
        Expr* n = a[2];
        if (n->type == EXPR_INTEGER && n->data.integer >= 0) {
            if (n->data.integer == 0)
                return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
            limit = (size_t)n->data.integer;
        } else if (!(n->type == EXPR_SYMBOL && n->data.symbol.name == SYM_Infinity)) {
            return NULL;
        }
    }

    Expr* subject = a[0];
    Expr* patt = a[1];

    /* All-mode enumerates exact-substring matches, so build anchored rules. */
    int anchored = (mode == REGEX_OV_ALL) ? 1 : 0;
    RegexRule* rules;
    int nr = regex_rules_build_ex(patt, anchored, caseless, &rules, "StringCases");
    if (nr < 0) return NULL;

    /* Capture offsets are only needed to instantiate a rule's RHS. */
    int want_captures = 0;
    for (int i = 0; i < nr; i++) if (rules[i].rhs) { want_captures = 1; break; }

    Expr* result;
    LLReader sv;
    if (ll_rd_open(subject, &sv)) {       /* a List or a LinkedList of subjects */
        size_t m = sv.n;
        Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
        if (!out) { ll_rd_close(&sv); regex_rules_free(rules, nr); return NULL; }
        size_t built = 0;
        for (; built < m; built++) {
            Expr* si = ll_rd_at(&sv, built);
            Expr* one = (si->type == EXPR_STRING)
                            ? sc_scalar(si->data.string, rules, nr, mode, want_captures, limit)
                            : expr_copy(si);
            if (!one) break;                     /* OOM: unwind and bail out */
            out[built] = one;
        }
        if (built < m) {                         /* partial: free and leave unevaluated */
            for (size_t k = 0; k < built; k++) expr_free(out[k]);
            free(out);
            ll_rd_close(&sv); regex_rules_free(rules, nr);
            return NULL;
        }
        result = ll_seq_like(subject, out, m);
        free(out);
        ll_rd_close(&sv);
    } else if (subject->type == EXPR_STRING) {
        result = sc_scalar(subject->data.string, rules, nr, mode, want_captures, limit);
    } else {
        result = NULL;   /* non-string subject: leave unevaluated */
    }

    regex_rules_free(rules, nr);
    return result;
}
