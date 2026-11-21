/*
 * stringposition.c - StringPosition[subject, pattern, n]
 *
 * Returns a List of {start, end} character-position pairs at which substrings of
 * `subject` match the string pattern `pattern`, in the 1-based inclusive form
 * consumed by StringTake / StringDrop / StringReplacePart.  The pattern may be a
 * literal string, a general string expression (Blank/Pattern/~~/RegularExpression
 * /character classes), or a List of patterns; a List of subjects threads.
 *
 * Options:
 *   Overlaps -> True (default) | False | All
 *     True  - include overlapping substrings, but only the first (natural) match
 *             starting at each position.
 *     False - exclude overlapping substrings (greedy left-to-right, global).
 *     All   - include every matching substring at every start (all lengths).
 *   IgnoreCase -> True | False (default)
 *     Treat upper/lowercase as equivalent.
 *
 * A third positional integer argument n keeps only the first n matches.
 *
 * The match enumeration itself is regex_scan() in regex_common.c, shared with
 * StringCases and StringCount; this file only turns spans into position pairs.
 *
 * Positions are CHARACTER positions (UTF-8 code points), consistent with
 * StringLength / StringTake: the scan's byte offsets are mapped through the
 * subject's character-offset table (a binary search; identity for ASCII).
 */

#include "picostrings.h"
#include "regex_common.h"
#include "linkedlist.h"   /* a LinkedList of subjects */
#include "sym_names.h"
#include "symtab.h"
#include "common.h"
#include "utf8.h"

#include <string.h>
#include <stdlib.h>

/* Character index of byte offset `b` (a character boundary) in view `v`. */
static size_t sp_char_index(const U8View* v, size_t b) {
    if (!v->off) return b;                         /* ASCII: identity */
    size_t lo = 0, hi = v->nchars;                 /* off[0..nchars] ascending */
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (v->off[mid] < b) lo = mid + 1; else hi = mid;
    }
    return lo;
}

/* Build the 1-based inclusive position pair List[ms+1, me] for the match
 * covering characters [ms, me). */
static Expr* sp_pair(size_t ms, size_t me) {
    Expr* p[2] = { expr_new_integer((int64_t)(ms + 1)), expr_new_integer((int64_t)me) };
    return expr_new_function(expr_new_symbol(SYM_List), p, 2);
}

/*
 * Compute the ordered list of position pairs for one subject string.  The match
 * enumeration itself lives in regex_scan (regex_common.c), shared with
 * StringCases and StringCount so the three agree on every Overlaps policy.
 * Returns a freshly allocated List[...] (empty on no match), or NULL on OOM.
 */
static Expr* sp_scalar(const char* subj, RegexRule* rules, int nr,
                       RegexOverlapMode mode, long n_limit) {
    size_t len = strlen(subj);
    RegexScan scan;
    if (regex_scan_limit(subj, len, rules, nr, mode, /*want_captures=*/0,
                         n_limit > 0 ? (size_t)n_limit : 0, &scan) < 0)
        return NULL;

    U8View v;
    if (u8_view_init(&v, subj, len) != 0) { regex_scan_free(&scan); return NULL; }

    size_t out_n = scan.count;
    Expr** items = malloc(sizeof(Expr*) * (out_n ? out_n : 1));
    if (!items) { u8_view_free(&v); regex_scan_free(&scan); return NULL; }
    for (size_t k = 0; k < out_n; k++)
        items[k] = sp_pair(sp_char_index(&v, scan.spans[k].ms),
                           sp_char_index(&v, scan.spans[k].me));
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), items, out_n);
    free(items);
    u8_view_free(&v);
    regex_scan_free(&scan);
    return result;
}

Expr* builtin_stringposition(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;

    if (argc == 0) return builtin_arg_error("StringPosition", 0, 1, 3);

    /* Seed option state from the registered defaults (so SetOptions[
     * StringPosition, ...] takes effect), then let explicit trailing options
     * override below. Defaults are {IgnoreCase -> False, Overlaps -> True}. */
    int caseless = 0;
    RegexOverlapMode mode = REGEX_OV_TRUE;
    Expr* defs = symtab_get_options("StringPosition");   /* borrowed */
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
    /* StringPosition[patt] is the operator form (src/opform.c). */
    if (pargc == 1) return NULL;
    if (pargc > 3)
        return builtin_arg_error("StringPosition", argc, 1, 3);

    Expr* subject = a[0];
    Expr* patt = a[1];

    /* Optional third positional argument: an integer occurrence count. */
    long n_limit = 0;
    if (pargc == 3) {
        Expr* third = a[2];
        if (third->type == EXPR_INTEGER) n_limit = (long)third->data.integer;
        else return NULL;   /* unsupported third arg: leave unevaluated */
    }

    /* All-mode enumerates exact-substring matches, so build anchored rules. */
    int anchored = (mode == REGEX_OV_ALL) ? 1 : 0;
    RegexRule* rules;
    int nr = regex_rules_build_ex(patt, anchored, caseless, &rules, "StringPosition");
    if (nr < 0) return NULL;

    Expr* result;
    LLElems sv;
    if (ll_elems_open(subject, &sv)) {       /* a List or a LinkedList of subjects */
        size_t m = sv.n;
        Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
        size_t built = 0;
        for (; built < m; built++) {
            Expr* si = sv.v[built];
            Expr* one = (si->type == EXPR_STRING)
                            ? sp_scalar(si->data.string, rules, nr, mode, n_limit)
                            : expr_copy(si);
            if (!one) break;                     /* OOM: unwind and bail out */
            out[built] = one;
        }
        if (built < m) {                         /* partial: free and leave unevaluated */
            for (size_t k = 0; k < built; k++) expr_free(out[k]);
            free(out);
            ll_elems_close(&sv); regex_rules_free(rules, nr);
            return NULL;
        }
        result = ll_seq_like(subject, out, m);
        free(out);
        ll_elems_close(&sv);
    } else if (subject->type == EXPR_STRING) {
        result = sp_scalar(subject->data.string, rules, nr, mode, n_limit);
    } else {
        result = NULL;   /* non-string subject: leave unevaluated */
    }

    regex_rules_free(rules, nr);
    return result;
}
