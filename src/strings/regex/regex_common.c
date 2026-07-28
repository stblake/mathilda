/*
 * regex_common.c - implementation of the Expr-aware regex helpers.
 */

#include "regex_common.h"
#include "sym_names.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* ------------------------------------------------------------------ */
/* Growable byte buffer                                               */
/* ------------------------------------------------------------------ */

int regexbuf_add(RegexBuf* b, const char* s, size_t n) {
    if (b->len + n + 1 > b->cap) {
        size_t nc = b->cap ? b->cap : 32;
        while (nc < b->len + n + 1) nc *= 2;
        char* np = realloc(b->p, nc);
        if (!np) return -1;
        b->p = np;
        b->cap = nc;
    }
    memcpy(b->p + b->len, s, n);
    b->len += n;
    b->p[b->len] = '\0';
    return 0;
}

/* ------------------------------------------------------------------ */
/* Small local helpers                                                */
/* ------------------------------------------------------------------ */

/* C99-safe strdup replacement. */
static char* rc_strdup(const char* s) {
    size_t n = strlen(s) + 1;
    char* d = malloc(n);
    if (d) memcpy(d, s, n);
    return d;
}

/* Wrap regex source as \A(?:src)\z for a whole-string (anchored) match. */
static char* wrap_anchored(const char* src) {
    size_t n = strlen(src);
    char* out = malloc(n + 9);              /* "\A(?:" (5) + src + ")\z" (3) + NUL */
    if (!out) return NULL;
    memcpy(out, "\\A(?:", 5);
    memcpy(out + 5, src, n);
    memcpy(out + 5 + n, ")\\z", 3);
    out[8 + n] = '\0';
    return out;
}

/* Regex source (malloc'd) for a single pattern element. Delegates to the shared
 * Wolfram-string-pattern translator (literal strings, RegularExpression, character
 * classes, StringExpression, Alternatives, Repeated, Except, ...). NULL if the
 * pattern is unsupported. */
static char* pattern_source(Expr* e) {
    return wl_pattern_to_regex(e, NULL);
}

/* Is e a Rule[...] or RuleDelayed[...] with two arguments? */
static int is_rule2(Expr* e, Expr** lhs, Expr** rhs) {
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        (e->data.function.head->data.symbol.name == SYM_Rule ||
         e->data.function.head->data.symbol.name == SYM_RuleDelayed) &&
        e->data.function.arg_count == 2) {
        *lhs = e->data.function.args[0];
        *rhs = e->data.function.args[1];
        return 1;
    }
    return 0;
}

/* Compile one element into `out`. Returns 0 on success, -1 on failure. */
static int build_one(Expr* e, int anchored, int caseless, RegexRule* out,
                     const char* head) {
    Expr* patt = e;
    Expr* rhs = NULL;
    Expr* l; Expr* r;
    if (is_rule2(e, &l, &r)) { patt = l; rhs = r; }

    char* src = pattern_source(patt);
    if (!src) return -1;                    /* unsupported pattern -> unevaluated */

    if (anchored) {
        char* a = wrap_anchored(src);
        free(src);
        if (!a) return -1;
        src = a;
    }
    if (caseless) {                          /* case-insensitive: prepend (?i) */
        size_t n = strlen(src);
        char* c = malloc(n + 5);            /* "(?i)" (4) + src + NUL */
        if (!c) { free(src); return -1; }
        memcpy(c, "(?i)", 4);
        memcpy(c + 4, src, n + 1);
        free(src);
        src = c;
    }

    char err[256];
    RegexProgram* prog = regex_compile(src, err, sizeof err);
    free(src);
    if (!prog) {
        fprintf(stderr, "%s::regex: %s\n", head, err);
        return -1;
    }
    out->prog = prog;
    out->rhs = rhs;
    out->lhs = patt;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Public: build / free a rule set                                    */
/* ------------------------------------------------------------------ */

int regex_rules_build(Expr* patt, int anchored, RegexRule** out, const char* head) {
    return regex_rules_build_ex(patt, anchored, /*caseless=*/0, out, head);
}

int regex_rules_build_ex(Expr* patt, int anchored, int caseless,
                         RegexRule** out, const char* head) {
    if (!regex_available()) {
        fprintf(stderr,
                "%s::regavail: regular-expression support is not available; "
                "rebuild Mathilda with PCRE2 (USE_REGEX=1).\n", head);
        return -1;
    }

    Expr** elems;
    int n;
    if (patt->type == EXPR_FUNCTION &&
        patt->data.function.head->type == EXPR_SYMBOL &&
        patt->data.function.head->data.symbol.name == SYM_List) {
        n = (int)patt->data.function.arg_count;
        elems = patt->data.function.args;
        if (n == 0) return -1;
    } else {
        n = 1;
        elems = &patt;
    }

    RegexRule* rules = calloc((size_t)n, sizeof(RegexRule));
    if (!rules) return -1;

    for (int i = 0; i < n; i++) {
        if (build_one(elems[i], anchored, caseless, &rules[i], head) != 0) {
            regex_rules_free(rules, i);     /* free the i already built */
            return -1;
        }
    }
    *out = rules;
    return n;
}

void regex_rules_free(RegexRule* rules, int n) {
    if (!rules) return;
    for (int i = 0; i < n; i++) regex_free(rules[i].prog);
    free(rules);
}

/* ------------------------------------------------------------------ */
/* Shared match scanner                                               */
/* ------------------------------------------------------------------ */

/*
 * Accumulator for one scan: the growable span array plus the optional flat
 * capture pool.  Spans index the pool by offset rather than pointer so that
 * growing the pool never invalidates already-recorded spans.
 */
typedef struct {
    RegexSpan* spans;
    size_t     count, cap;
    size_t*    caps;
    size_t     caps_len, caps_cap;
    int        want_captures;
} ScanState;

/*
 * Record one match.  `ov` holds `pairs` (start, end) offset pairs as written by
 * regex_match, relative to `base` (nonzero only in the All mode, which matches
 * against subj + p); they are rebased so the pool is always subject-relative.
 * Returns 0, or -1 on allocation failure.
 */
static int scan_push(ScanState* st, size_t ms, size_t me, int rule,
                     const size_t* ov, size_t pairs, size_t base) {
    if (st->count == st->cap) {
        size_t nc = st->cap ? st->cap * 2 : 16;
        RegexSpan* ns = realloc(st->spans, nc * sizeof(RegexSpan));
        if (!ns) return -1;
        st->spans = ns;
        st->cap = nc;
    }

    size_t off = 0, np = 0;
    if (st->want_captures && ov && pairs) {
        np = pairs;
        if (st->caps_len + 2 * np > st->caps_cap) {
            size_t nc = st->caps_cap ? st->caps_cap : 32;
            while (nc < st->caps_len + 2 * np) nc *= 2;
            size_t* ncaps = realloc(st->caps, nc * sizeof(size_t));
            if (!ncaps) return -1;
            st->caps = ncaps;
            st->caps_cap = nc;
        }
        off = st->caps_len;
        for (size_t k = 0; k < 2 * np; k++)
            st->caps[off + k] =
                (ov[k] == REGEX_UNSET) ? REGEX_UNSET : ov[k] + base;
        st->caps_len += 2 * np;
    }

    st->spans[st->count].ms       = ms;
    st->spans[st->count].me       = me;
    st->spans[st->count].rule     = rule;
    st->spans[st->count].caps_off = off;
    st->spans[st->count].npairs   = np;
    st->count++;
    return 0;
}

/*
 * Offset pairs to request per rule: the whole match plus every capture group,
 * clamped to the range we expose ($0..$63).  Computed once per scan -- All-mode
 * probes O(len^2 * nr) substrings, so querying the group count inside that loop
 * would be pure overhead.  Returns a malloc'd array, or NULL on failure.
 */
static size_t* scan_pairs_table(RegexRule* rules, int nr) {
    size_t* pairs = malloc((size_t)nr * sizeof(size_t));
    if (!pairs) return NULL;
    for (int i = 0; i < nr; i++) {
        int gc = regex_group_count(rules[i].prog) + 1;
        pairs[i] = (gc > REGEX_MAX_PAIRS) ? (size_t)REGEX_MAX_PAIRS : (size_t)gc;
    }
    return pairs;
}

/*
 * Overlaps -> False.  Streaming left-to-right: at each position take the match
 * with the smallest start across all rules (ties broken by rule order), then
 * resume at its end.  Zero-width matches advance by one so the loop always
 * makes progress.
 */
static int scan_false(const char* subj, size_t len, RegexRule* rules, int nr,
                      const size_t* pairs, ScanState* st) {
    size_t pos = 0;
    while (pos <= len) {
        int found = 0, best_rule = 0;
        size_t best_ms = 0, best_me = 0, best_pairs = 0;
        size_t best_ov[REGEX_MAX_PAIRS * 2];

        for (int i = 0; i < nr; i++) {
            size_t ov[REGEX_MAX_PAIRS * 2];
            if (regex_match(rules[i].prog, subj, len, pos, ov, pairs[i]) != 1)
                continue;
            if (!found || ov[0] < best_ms) {
                found = 1;
                best_rule = i;
                best_ms = ov[0];
                best_me = ov[1];
                best_pairs = pairs[i];
                memcpy(best_ov, ov, sizeof(size_t) * 2 * pairs[i]);
            }
        }
        if (!found) break;

        if (scan_push(st, best_ms, best_me, best_rule, best_ov, best_pairs, 0) != 0)
            return -1;
        pos = (best_me > best_ms) ? best_me : best_ms + 1;   /* progress */
    }
    return 0;
}

/*
 * Overlaps -> True.  Per rule independently, the leftmost match at or after each
 * position, advancing by one past the match START so every distinct match start
 * is enumerated.
 */
static int scan_true(const char* subj, size_t len, RegexRule* rules, int nr,
                     const size_t* pairs, ScanState* st) {
    for (int i = 0; i < nr; i++) {
        size_t pos = 0;
        while (pos <= len) {
            size_t ov[REGEX_MAX_PAIRS * 2];
            if (regex_match(rules[i].prog, subj, len, pos, ov, pairs[i]) != 1) break;
            if (scan_push(st, ov[0], ov[1], i, ov, pairs[i], 0) != 0) return -1;
            pos = ov[0] + 1;   /* past the match start (overlap-friendly) */
        }
    }
    return 0;
}

/*
 * Overlaps -> All.  Every matching substring at every start: for each start p
 * ascending and end e descending, test an exact match of subj[p, e).  The rules
 * are anchored (\A(?:...)\z) by the caller, so a match means the pattern covers
 * [p, e) exactly.
 */
static int scan_all(const char* subj, size_t len, RegexRule* rules, int nr,
                    const size_t* pairs, ScanState* st) {
    for (size_t p = 0; p < len; p++) {
        for (size_t e = len; e > p; e--) {
            for (int i = 0; i < nr; i++) {
                size_t ov[REGEX_MAX_PAIRS * 2];
                if (regex_match(rules[i].prog, subj + p, e - p, 0, ov, pairs[i]) != 1)
                    continue;
                if (scan_push(st, p, e, i, ov, pairs[i], p) != 0) return -1;
            }
        }
    }
    return 0;
}

/* Stable insertion sort by start offset (ties keep discovery order: rule index,
 * then match order). qsort is not stable, and the input is nearly sorted -- it
 * is a concatenation of per-rule ascending runs -- so this is close to linear. */
static void scan_stable_sort(RegexSpan* a, size_t n) {
    for (size_t i = 1; i < n; i++) {
        RegexSpan key = a[i];
        size_t j = i;
        while (j > 0 && a[j - 1].ms > key.ms) { a[j] = a[j - 1]; j--; }
        a[j] = key;
    }
}

long regex_scan(const char* subj, size_t len, RegexRule* rules, int nr,
                RegexOverlapMode mode, int want_captures, RegexScan* out) {
    ScanState st;
    memset(&st, 0, sizeof st);
    st.want_captures = want_captures;

    out->spans = NULL;
    out->count = 0;
    out->caps  = NULL;

    size_t* pairs = scan_pairs_table(rules, nr);
    if (!pairs) return -1;

    int rc;
    if (mode == REGEX_OV_ALL)        rc = scan_all(subj, len, rules, nr, pairs, &st);
    else if (mode == REGEX_OV_FALSE) rc = scan_false(subj, len, rules, nr, pairs, &st);
    else                             rc = scan_true(subj, len, rules, nr, pairs, &st);
    free(pairs);

    if (rc != 0) {                  /* OOM mid-scan: *out stays the zeroed state */
        free(st.spans);
        free(st.caps);
        return -1;
    }

    /* The False and All scans already emit in ascending start order; only the
     * True scan concatenates one run per rule and needs merging. */
    if (mode == REGEX_OV_TRUE && nr > 1) scan_stable_sort(st.spans, st.count);

    out->spans = st.spans;
    out->count = st.count;
    out->caps  = st.caps;
    return (long)st.count;
}

void regex_scan_free(RegexScan* s) {
    if (!s) return;
    free(s->spans);
    free(s->caps);
    s->spans = NULL;
    s->caps  = NULL;
    s->count = 0;
}

int regex_match_opt(const Expr* e, const char* opt_sym, int* value, int overlaps) {
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
            *value = REGEX_OV_ALL;
        else if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_False)
            *value = REGEX_OV_FALSE;
        else
            *value = REGEX_OV_TRUE;
    } else {
        *value = (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_True);
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* Public: $n template expansion                                      */
/* ------------------------------------------------------------------ */

char* regex_expand_template(const char* tpl, const char* subj,
                            const size_t* ov, size_t npairs) {
    RegexBuf b = {0};
    for (size_t i = 0; tpl[i];) {
        if (tpl[i] == '$') {
            char c = tpl[i + 1];
            if (c == '$') {                          /* $$ -> literal $ */
                if (regexbuf_add(&b, "$", 1)) goto oom;
                i += 2;
                continue;
            }
            if (c >= '0' && c <= '9') {              /* $n -> group n */
                size_t k = 0, j = i + 1;
                while (tpl[j] >= '0' && tpl[j] <= '9') {
                    k = k * 10 + (size_t)(tpl[j] - '0');
                    j++;
                }
                if (k < npairs && ov[2 * k] != REGEX_UNSET) {
                    size_t s = ov[2 * k], e = ov[2 * k + 1];
                    if (regexbuf_add(&b, subj + s, e - s)) goto oom;
                }
                i = j;
                continue;
            }
            if (regexbuf_add(&b, "$", 1)) goto oom;   /* lone $ */
            i += 1;
            continue;
        }
        if (regexbuf_add(&b, &tpl[i], 1)) goto oom;
        i += 1;
    }
    if (!b.p) return rc_strdup("");                   /* empty template */
    return b.p;
oom:
    free(b.p);
    return NULL;
}

char* regex_rule_replacement(const RegexRule* r, const char* subj,
                             const size_t* ov, size_t npairs) {
    if (!r->rhs) return NULL;
    if (r->rhs->type == EXPR_STRING)
        return regex_expand_template(r->rhs->data.string, subj, ov, npairs);
    return NULL;                                       /* non-string RHS: out of scope */
}
