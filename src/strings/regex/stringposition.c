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
 * Byte semantics: like the rest of src/strings, positions are byte offsets (no
 * UTF-8 codepoint decoding), consistent with StringLength / StringPart.
 */

#include "picostrings.h"
#include "regex_common.h"
#include "sym_names.h"
#include "symtab.h"
#include "common.h"

#include <string.h>
#include <stdlib.h>

/* Overlap handling mode for the third-per-position policy. */
typedef enum { OV_TRUE = 0, OV_FALSE = 1, OV_ALL = 2 } OverlapMode;

/* One recorded match: half-open byte offsets [ms, me) plus discovery rule index
 * (used only to keep a deterministic stable order for equal start positions). */
typedef struct { size_t ms, me; int rule; } SpMatch;

/* Build the 1-based inclusive position pair List[ms+1, me] for match [ms, me). */
static Expr* sp_pair(size_t ms, size_t me) {
    Expr** p = malloc(sizeof(Expr*) * 2);
    p[0] = expr_new_integer((int64_t)(ms + 1));
    p[1] = expr_new_integer((int64_t)me);
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), p, 2);
    free(p);   /* expr_new_function copies the args array; free our copy */
    return r;
}

/* Append a match to a growable SpMatch array. Returns 0, or -1 on OOM. */
static int sp_push(SpMatch** arr, size_t* count, size_t* cap,
                   size_t ms, size_t me, int rule) {
    if (*count == *cap) {
        size_t nc = *cap ? *cap * 2 : 16;
        SpMatch* na = realloc(*arr, nc * sizeof(SpMatch));
        if (!na) return -1;
        *arr = na;
        *cap = nc;
    }
    (*arr)[(*count)++] = (SpMatch){ ms, me, rule };
    return 0;
}

/*
 * Collect, per rule independently, one match per starting position (Overlaps->
 * True candidate set): the leftmost match at or after each position, advancing by
 * one past the match START so every overlapping match-start is enumerated.
 */
static int sp_collect_true(const char* subj, size_t len, RegexRule* rules, int nr,
                           SpMatch** arr, size_t* count, size_t* cap) {
    for (int i = 0; i < nr; i++) {
        int gc = regex_group_count(rules[i].prog) + 1;
        size_t pairs = (gc > REGEX_MAX_PAIRS) ? REGEX_MAX_PAIRS : (size_t)gc;
        size_t pos = 0;
        while (pos <= len) {
            size_t ov[REGEX_MAX_PAIRS * 2];
            if (regex_match(rules[i].prog, subj, len, pos, ov, pairs) != 1) break;
            if (sp_push(arr, count, cap, ov[0], ov[1], i) != 0) return -1;
            pos = ov[0] + 1;   /* advance past the match start (overlap-friendly) */
        }
    }
    return 0;
}

/* Stable insertion sort of the match array by start offset (ties keep the
 * existing order: rule index, then discovery order). qsort is not stable. */
static void sp_stable_sort(SpMatch* a, size_t n) {
    for (size_t i = 1; i < n; i++) {
        SpMatch key = a[i];
        size_t j = i;
        while (j > 0 && a[j - 1].ms > key.ms) { a[j] = a[j - 1]; j--; }
        a[j] = key;
    }
}

/*
 * Compute the ordered list of position pairs for one subject string.
 * Returns a freshly allocated List[...] (empty on no match).
 */
static Expr* sp_scalar(const char* subj, RegexRule* rules, int nr,
                       OverlapMode mode, long n_limit) {
    size_t len = strlen(subj);
    SpMatch* arr = NULL;
    size_t count = 0, cap = 0;

    if (mode == OV_ALL) {
        /* Every matching substring at every start: for each start p ascending
         * and end e descending, test an exact (anchored) match of subj[p..e).
         * `rules` are built anchored (\A(?:...)\z), so a whole-substring match
         * means the pattern matches [p, e) exactly. */
        for (size_t p = 0; p < len; p++) {
            for (size_t e = len; e > p; e--) {
                for (int i = 0; i < nr; i++) {
                    size_t ov[2];
                    if (regex_match(rules[i].prog, subj + p, e - p, 0, ov, 1) == 1) {
                        if (sp_push(&arr, &count, &cap, p, e, i) != 0) goto build;
                    }
                }
            }
        }
    } else {
        if (sp_collect_true(subj, len, rules, nr, &arr, &count, &cap) != 0) goto build;
        sp_stable_sort(arr, count);
        if (mode == OV_FALSE) {
            /* Greedy global non-overlap: keep a match only when it starts at or
             * after the end of the last kept match. */
            size_t kept = 0, last_me = 0;
            int have = 0;
            for (size_t k = 0; k < count; k++) {
                if (!have || arr[k].ms >= last_me) {
                    arr[kept++] = arr[k];
                    last_me = arr[k].me;
                    have = 1;
                }
            }
            count = kept;
        }
    }

build:;
    size_t out_n = count;
    if (n_limit > 0 && (size_t)n_limit < out_n) out_n = (size_t)n_limit;

    Expr** items = malloc(sizeof(Expr*) * (out_n ? out_n : 1));
    for (size_t k = 0; k < out_n; k++) items[k] = sp_pair(arr[k].ms, arr[k].me);
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), items, out_n);
    free(items);
    free(arr);
    return result;
}

/* Is e a supported option Rule/RuleDelayed[opt_sym, value]?  Sets *value from
 * the RHS: for IgnoreCase, 1 iff True; for Overlaps, the OverlapMode integer. */
static int sp_match_opt(const Expr* e, const char* opt_sym, int* value,
                        int overlaps) {
    if (e->type != EXPR_FUNCTION ||
        e->data.function.head->type != EXPR_SYMBOL ||
        (e->data.function.head->data.symbol.name != SYM_Rule &&
         e->data.function.head->data.symbol.name != SYM_RuleDelayed) ||
        e->data.function.arg_count != 2 ||
        e->data.function.args[0]->type != EXPR_SYMBOL ||
        e->data.function.args[0]->data.symbol.name != opt_sym)
        return 0;
    Expr* v = e->data.function.args[1];
    if (overlaps) {
        if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_All)
            *value = OV_ALL;
        else if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_False)
            *value = OV_FALSE;
        else
            *value = OV_TRUE;
    } else {
        *value = (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_True);
    }
    return 1;
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
    OverlapMode mode = OV_TRUE;
    Expr* defs = symtab_get_options("StringPosition");   /* borrowed */
    if (defs && defs->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < defs->data.function.arg_count; i++) {
            int v;
            if (sp_match_opt(defs->data.function.args[i], SYM_IgnoreCase, &v, 0))
                caseless = v;
            else if (sp_match_opt(defs->data.function.args[i], SYM_Overlaps, &v, 1))
                mode = (OverlapMode)v;
        }
    }

    /* Strip trailing IgnoreCase / Overlaps options, leaving positional args. */
    size_t pargc = argc;
    while (pargc >= 2) {
        int v;
        if (sp_match_opt(a[pargc - 1], SYM_IgnoreCase, &v, 0)) {
            caseless = v; pargc--;
        } else if (sp_match_opt(a[pargc - 1], SYM_Overlaps, &v, 1)) {
            mode = (OverlapMode)v; pargc--;
        } else {
            break;
        }
    }
    if (pargc < 2 || pargc > 3)
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
    int anchored = (mode == OV_ALL) ? 1 : 0;
    RegexRule* rules;
    int nr = regex_rules_build_ex(patt, anchored, caseless, &rules, "StringPosition");
    if (nr < 0) return NULL;

    Expr* result;
    if (subject->type == EXPR_FUNCTION &&
        subject->data.function.head->type == EXPR_SYMBOL &&
        subject->data.function.head->data.symbol.name == SYM_List) {
        size_t m = subject->data.function.arg_count;
        Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
        for (size_t i = 0; i < m; i++) {
            Expr* si = subject->data.function.args[i];
            out[i] = (si->type == EXPR_STRING)
                         ? sp_scalar(si->data.string, rules, nr, mode, n_limit)
                         : expr_copy(si);
        }
        result = expr_new_function(expr_new_symbol(SYM_List), out, m);
        free(out);
    } else if (subject->type == EXPR_STRING) {
        result = sp_scalar(subject->data.string, rules, nr, mode, n_limit);
    } else {
        result = NULL;   /* non-string subject: leave unevaluated */
    }

    regex_rules_free(rules, nr);
    return result;
}
