/*
 * regex_common.c - implementation of the Expr-aware regex helpers.
 */

#include "regex_common.h"
#include "linkedlist.h"   /* a LinkedList of patterns */
#include "sym_names.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "match.h"     /* MatchEnv / replace_bindings: RHS instantiation */
#include "eval.h"
#include "utf8.h"

#include <stdint.h>
#include <stddef.h>

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

/* Length of the leading PCRE start-of-pattern verbs -- (*UTF), (*UCP), (*CR),
 * ... -- which must stay at the very front of the compiled source. */
static size_t leading_verbs_len(const char* src) {
    size_t i = 0;
    while (src[i] == '(' && src[i + 1] == '*') {
        const char* close = strchr(src + i, ')');
        if (!close) break;
        i = (size_t)(close - src) + 1;
    }
    return i;
}

/* Wrap regex source as \A(?:src)\z for a whole-string (anchored) match,
 * keeping any leading (*VERB)s in front. */
static char* wrap_anchored(const char* src) {
    size_t v = leading_verbs_len(src);
    size_t n = strlen(src);
    char* out = malloc(n + 9);              /* "\A(?:" (5) + src + ")\z" (3) + NUL */
    if (!out) return NULL;
    memcpy(out, src, v);
    memcpy(out + v, "\\A(?:", 5);
    memcpy(out + v + 5, src + v, n - v);
    memcpy(out + 5 + n, ")\\z", 3);
    out[8 + n] = '\0';
    return out;
}

/* New EXPR_STRING holding subj[s, e). */
static Expr* rc_substr(const char* subj, size_t s, size_t e) {
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
 * A MatchEnv binding each pattern name of rule `r` to the text of the first
 * participating group that carries it (a name repeated across Alternatives
 * branches owns a group in the branch where it first occurs). A name with no
 * participating group is bound to "" when `unset_empty` (a finished match:
 * Mathematica binds the untaken branch's name to ""), else left unbound (a
 * Condition mid-match: the name has not been reached yet).
 */
static MatchEnv* rule_bind_env(const RegexRule* r, const char* subj,
                               const size_t* ov, size_t npairs, int unset_empty) {
    MatchEnv* env = env_new();
    const char* const* gname = r->info.group_name;
    for (int g = 1; g <= r->info.ngroups; g++) {
        const char* nm = gname[g];
        if (!nm || env_get(env, nm)) continue;
        Expr* v = NULL;
        for (int h = g; h <= r->info.ngroups && !v; h++) {
            if (gname[h] == nm && (size_t)h < npairs && ov[2 * h] != REGEX_UNSET)
                v = rc_substr(subj, ov[2 * h], ov[2 * h + 1]);
        }
        if (!v && unset_empty) v = expr_new_string("");
        if (v) { env_set(env, nm, v); expr_free(v); }
    }
    return env;
}

/*
 * Callout hook for a translated pattern (see string_pattern.c): evaluate the
 * Condition / PatternTest it marks against the match in progress. Returns 0 to
 * let the match continue, 1 to make the engine backtrack.
 */
static int rule_callout(int number, const char* subj, size_t len,
                        const size_t* ov, size_t npairs, void* data) {
    (void)len;
    const RegexRule* r = (const RegexRule*)data;
    if (number < 1 || number > r->info.ncallouts) return 0;
    const SpCallout* c = &r->info.callouts[number - 1];

    Expr* test;
    if (c->kind == SP_CALLOUT_TEST) {
        Expr* text = ((size_t)c->group < npairs && ov[2 * c->group] != REGEX_UNSET)
                   ? rc_substr(subj, ov[2 * c->group], ov[2 * c->group + 1])
                   : expr_new_string("");
        Expr* targs[1] = { text };
        test = expr_new_function(expr_copy(c->expr), targs, 1);
    } else {
        MatchEnv* env = rule_bind_env(r, subj, ov, npairs, /*unset_empty=*/0);
        test = replace_bindings(c->expr, env);
        env_free(env);
    }
    Expr* v = eval_and_free(test);
    int pass = v && v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_True;
    if (v) expr_free(v);
    return pass ? 0 : 1;
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
static int build_one(Expr* e, unsigned anchored, unsigned regex_flags,
                     RegexRule* out, const char* head) {
    Expr* patt = e;
    Expr* rhs = NULL;
    Expr* l; Expr* r;
    if (is_rule2(e, &l, &r)) { patt = l; rhs = r; }

    /* Shared Wolfram-string-pattern translator (literal strings,
     * RegularExpression, character classes, StringExpression, Alternatives,
     * Repeated, Except, ...), keeping its capture-name table for the RHS. */
    char* src = wl_pattern_to_regex_wild(patt, &out->info,
                                         anchored & REGEX_RULE_WILDCARDS);
    if (!src) return -1;                    /* unsupported pattern -> unevaluated */

    if (anchored & REGEX_RULE_ANCHORED) {
        char* a = wrap_anchored(src);
        free(src);
        if (!a) { sp_info_free(&out->info); return -1; }
        src = a;
    }

    /* A translated pattern runs DOTALL (blanks match newlines); a bare
     * RegularExpression keeps PCRE's defaults. */
    unsigned flags = regex_flags;
    if (!out->info.raw_regex)  flags |= REGEX_DOTALL;

    char err[256];
    RegexProgram* prog = regex_compile_ex(src, flags,
                                          out->info.ncallouts ? rule_callout : NULL,
                                          out, err, sizeof err);
    free(src);
    if (!prog) {
        mth_message(head, "regex", "%s", err);
        sp_info_free(&out->info);
        return -1;
    }
    out->prog = prog;
    out->rhs = rhs;
    out->lhs = patt;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Public: build / free a rule set, with a compiled-pattern cache      */
/* ------------------------------------------------------------------ */

/*
 * A rule set lives in a refcounted block. Building is expensive -- pattern
 * translation, pcre2_compile and a JIT compile per rule -- and string code
 * routinely applies ONE pattern to many strings (Map[StringMatchQ[#, p] &, list],
 * Select[..., StringContainsQ[p]]), so built sets are kept in a small LRU cache
 * keyed by the pattern expression (structural hash + equality) and the build
 * flags. A hit costs a hash lookup instead of a compile.
 *
 * The cache holds a reference to the pattern (expr_copy is a refcount bump), so
 * the rules' borrowed lhs/rhs/callout pointers stay valid while it is cached.
 * Each regex_rules_build* result holds one more reference, released by
 * regex_rules_free; a block is destroyed when both are gone. The evaluator is
 * single-threaded, so the cache needs no lock; nested use of one cached set
 * (a callout evaluating the same pattern) is safe because regex_match gives a
 * nested match private match data.
 */
typedef struct RulesBlock {
    int        refs;
    int        nr;
    Expr*      key;            /* owned pattern reference; NULL if uncached */
    uint64_t   hash;
    unsigned   rule_flags, regex_flags;
    RegexRule  rules[];        /* C99 flexible array member */
} RulesBlock;

#define RULES_CACHE_SIZE 64
static RulesBlock* g_rules_cache[RULES_CACHE_SIZE];   /* [0] most recent */
static int         g_rules_cache_n = 0;

static RulesBlock* block_of(RegexRule* rules) {
    return (RulesBlock*)((char*)rules - offsetof(RulesBlock, rules));
}

static void block_destroy(RulesBlock* b) {
    for (int i = 0; i < b->nr; i++) {
        regex_free(b->rules[i].prog);
        sp_info_free(&b->rules[i].info);
    }
    if (b->key) expr_free(b->key);
    free(b);
}

static void block_release(RulesBlock* b) {
    if (--b->refs == 0) block_destroy(b);
}

/* Insert at the front, evicting the least recently used entry if full. */
static void cache_insert(RulesBlock* b) {
    if (g_rules_cache_n == RULES_CACHE_SIZE) {
        block_release(g_rules_cache[RULES_CACHE_SIZE - 1]);
        g_rules_cache_n--;
    }
    memmove(g_rules_cache + 1, g_rules_cache,
            (size_t)g_rules_cache_n * sizeof(RulesBlock*));
    g_rules_cache[0] = b;
    g_rules_cache_n++;
    b->refs++;
}

static RulesBlock* cache_lookup(Expr* patt, uint64_t h, unsigned rf, unsigned xf) {
    for (int i = 0; i < g_rules_cache_n; i++) {
        RulesBlock* b = g_rules_cache[i];
        if (b->hash != h || b->rule_flags != rf || b->regex_flags != xf ||
            !expr_eq(b->key, patt))
            continue;
        if (i > 0) {                                   /* move to front */
            memmove(g_rules_cache + 1, g_rules_cache, (size_t)i * sizeof(RulesBlock*));
            g_rules_cache[0] = b;
        }
        return b;
    }
    return NULL;
}

void regex_cache_clear(void) {
    for (int i = 0; i < g_rules_cache_n; i++) block_release(g_rules_cache[i]);
    g_rules_cache_n = 0;
}

int regex_rules_build(Expr* patt, int anchored, RegexRule** out, const char* head) {
    return regex_rules_build_ex(patt, anchored, /*caseless=*/0, out, head);
}

int regex_rules_build_ex(Expr* patt, int anchored, int caseless,
                         RegexRule** out, const char* head) {
    return regex_rules_build_opts(patt, (unsigned)anchored,
                                  caseless ? REGEX_CASELESS : 0u, out, head);
}

/* True when a LinkedList occurs anywhere in a (small) pattern expression. */
static bool rx_has_linkedlist(const Expr* e) {
    if (e->type == EXPR_LINKEDLIST) return true;
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rx_has_linkedlist(e->data.function.args[i])) return true;
    return false;
}

int regex_rules_build_opts(Expr* patt, unsigned rule_flags, unsigned regex_flags,
                           RegexRule** out, const char* head) {
    if (!regex_available()) {
        mth_message(head, "regavail",
                "regular-expression support is not available; "
                "rebuild Mathilda with PCRE2 (USE_REGEX=1).");
        return -1;
    }
    regex_flags |= REGEX_UTF;              /* strings are UTF-8, classes Unicode */

    /* A LinkedList anywhere in the pattern is a mutable object compared by
     * identity: a rule set keyed on it would go stale when it changed, so such
     * a pattern is built afresh every call and never cached. */
    bool mutable_patt = rx_has_linkedlist(patt);
    uint64_t h = mutable_patt ? 0 : expr_hash(patt);
    RulesBlock* hit = mutable_patt ? NULL : cache_lookup(patt, h, rule_flags, regex_flags);
    if (hit) {
        hit->refs++;
        *out = hit->rules;
        return hit->nr;
    }

    /* A top-level List -- or LinkedList -- of patterns is one rule per
     * element. The rules borrow into the elements, so a LinkedList is read
     * from a snapshot (ll_copy) that the block keeps as its key: its boxed
     * slots stay alive for the block's lifetime whatever happens to the list.
     * A typed (numeric) LinkedList holds no string pattern at all. */
    Expr** elems = NULL;
    int n;
    Expr* key;
    const LinkedList* pl = expr_linkedlist(patt);
    const LinkedList* snap = NULL;          /* elements read slot by slot */
    if (pl) {
        if (ll_kind(pl) != LL_BOXED) return -1;
        key = expr_new_linkedlist(ll_copy(pl));
        snap = key->data.llist;
        n = (int)ll_length(snap);
    } else if (patt->type == EXPR_FUNCTION &&
        patt->data.function.head->type == EXPR_SYMBOL &&
        patt->data.function.head->data.symbol.name == SYM_List) {
        n = (int)patt->data.function.arg_count;
        elems = patt->data.function.args;
        if (n == 0) return -1;
        key = NULL;
    } else {
        n = 1;
        elems = &patt;
        key = NULL;
    }
    if (n == 0) { expr_free(key); return -1; }

    RulesBlock* b = calloc(1, sizeof(RulesBlock) + (size_t)n * sizeof(RegexRule));
    if (!b) { expr_free(key); return -1; }
    b->refs = 1;
    /* The snapshot's elements are held by the snapshot itself (the key). */
    b->key = key ? key : expr_copy(patt);
    for (int i = 0; i < n; i++) {
        Expr* el = snap ? (Expr*)ll_peek_boxed(snap, i) : elems[i];
        if (build_one(el, rule_flags, regex_flags, &b->rules[i], head) != 0) {
            b->nr = i;                      /* free the i already built */
            block_destroy(b);
            return -1;
        }
        int gc = regex_group_count(b->rules[i].prog) + 1;
        b->rules[i].npairs = (gc > REGEX_MAX_PAIRS) ? (size_t)REGEX_MAX_PAIRS : (size_t)gc;
    }
    b->nr = n;
    b->hash = h;
    b->rule_flags = rule_flags;
    b->regex_flags = regex_flags;
    if (!mutable_patt) cache_insert(b);
    *out = b->rules;
    return n;
}

void regex_rules_free(RegexRule* rules, int n) {
    (void)n;
    if (rules) block_release(block_of(rules));
}

/* ------------------------------------------------------------------ */
/* Character stepping                                                 */
/* ------------------------------------------------------------------ */

size_t regex_next_char(const char* subj, size_t len, size_t pos) {
    if (pos >= len) return pos + 1;
    return pos + u8_char_len(subj, len, pos);
}

/* ------------------------------------------------------------------ */
/* Multi-rule cursor: leftmost match across a rule set                 */
/* ------------------------------------------------------------------ */

int regex_cursor_init(RegexCursor* c, RegexRule* rules, int nr,
                      const char* subj, size_t len) {
    c->rules = rules;
    c->nr = nr;
    c->subj = subj;
    c->len = len;
    c->from = malloc((size_t)nr * sizeof(size_t));
    c->found = malloc((size_t)nr * sizeof(int));
    c->ov = malloc((size_t)nr * REGEX_MAX_PAIRS * 2 * sizeof(size_t));
    if (!c->from || !c->found || !c->ov) { regex_cursor_free(c); return -1; }
    for (int i = 0; i < nr; i++) c->from[i] = REGEX_UNSET;   /* nothing cached */
    return 0;
}

void regex_cursor_free(RegexCursor* c) {
    free(c->from);
    free(c->found);
    free(c->ov);
    c->from = NULL;
    c->found = NULL;
    c->ov = NULL;
}

/*
 * The match a regex finds at a given START is independent of where the search
 * began, so rule i's leftmost match from position p is still its leftmost match
 * from any later q <= that match's start. Each rule therefore keeps its last
 * result and is re-searched only when the cursor has moved past it -- a scan is
 * O(n · rules) instead of re-running every rule at every step.
 */
int regex_cursor_next(RegexCursor* c, size_t pos, int* rule,
                      const size_t** ov, size_t* npairs) {
    int best = -1;
    for (int i = 0; i < c->nr; i++) {
        size_t* rov = c->ov + (size_t)i * REGEX_MAX_PAIRS * 2;
        int valid = c->from[i] != REGEX_UNSET && c->from[i] <= pos &&
                    (!c->found[i] || rov[0] >= pos);
        if (!valid) {
            int rc = regex_match(c->rules[i].prog, c->subj, c->len, pos,
                                 rov, c->rules[i].npairs);
            if (rc < 0) return -1;
            c->from[i] = pos;
            c->found[i] = rc;
        }
        if (c->found[i] &&
            (best < 0 || rov[0] < c->ov[(size_t)best * REGEX_MAX_PAIRS * 2]))
            best = i;
    }
    if (best < 0) return 0;
    *rule = best;
    *ov = c->ov + (size_t)best * REGEX_MAX_PAIRS * 2;
    *npairs = c->rules[best].npairs;
    return 1;
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
    size_t     limit;                  /* stop after this many spans (0: all) */
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

static int scan_full(const ScanState* st) {
    return st->limit && st->count >= st->limit;
}

/*
 * Overlaps -> False.  Streaming left-to-right: at each position take the match
 * with the smallest start across all rules (ties broken by rule order), then
 * resume at its end.  A zero-width match advances by one character so the loop
 * always makes progress.
 */
static int scan_false(const char* subj, size_t len, RegexRule* rules, int nr,
                      ScanState* st) {
    RegexCursor cur;
    if (regex_cursor_init(&cur, rules, nr, subj, len) != 0) return -1;
    size_t pos = 0;
    int rc = 0;
    while (pos <= len && !scan_full(st)) {
        int rule;
        const size_t* ov;
        size_t np;
        int m = regex_cursor_next(&cur, pos, &rule, &ov, &np);
        if (m < 0) { rc = -1; break; }
        if (m == 0) break;
        if (scan_push(st, ov[0], ov[1], rule, ov, np, 0) != 0) { rc = -1; break; }
        pos = (ov[1] > ov[0]) ? ov[1] : regex_next_char(subj, len, ov[0]);
    }
    regex_cursor_free(&cur);
    return rc;
}

/*
 * Overlaps -> True.  Per rule independently, the leftmost match at or after each
 * position, advancing one character past the match START so every distinct
 * match start is enumerated.
 */
static int scan_true(const char* subj, size_t len, RegexRule* rules, int nr,
                     ScanState* st) {
    for (int i = 0; i < nr; i++) {
        size_t pos = 0;
        while (pos <= len) {
            size_t ov[REGEX_MAX_PAIRS * 2];
            int m = regex_match(rules[i].prog, subj, len, pos, ov, rules[i].npairs);
            if (m < 0) return -1;
            if (m == 0) break;
            if (scan_push(st, ov[0], ov[1], i, ov, rules[i].npairs, 0) != 0) return -1;
            pos = regex_next_char(subj, len, ov[0]);
        }
    }
    return 0;
}

/*
 * Overlaps -> All.  Every matching substring at every start: for each character
 * start p ascending and character end e descending, test an exact match of
 * subj[p, e).  The rules are anchored (\A(?:...)\z) by the caller, so a match
 * means the pattern covers [p, e) exactly.
 */
static int scan_all(const char* subj, size_t len, RegexRule* rules, int nr,
                    ScanState* st) {
    U8View v;
    if (u8_view_init(&v, subj, len) != 0) return -1;
    int rc = 0;
    for (size_t pc = 0; pc < v.nchars && rc == 0; pc++) {
        size_t p = u8_view_off(&v, pc);
        for (size_t ec = v.nchars; ec > pc && rc == 0; ec--) {
            size_t e = u8_view_off(&v, ec);
            for (int i = 0; i < nr; i++) {
                size_t ov[REGEX_MAX_PAIRS * 2];
                int m = regex_match(rules[i].prog, subj + p, e - p, 0, ov, rules[i].npairs);
                if (m < 0) { rc = -1; break; }
                if (m == 0) continue;
                if (scan_push(st, p, e, i, ov, rules[i].npairs, p) != 0) { rc = -1; break; }
            }
        }
    }
    u8_view_free(&v);
    return rc;
}

/* Stable merge sort of spans by start offset (ties keep discovery order: rule
 * index, then match order). The True scan concatenates one ascending run per
 * rule, so this is O(n log n) worst case where the old insertion sort was
 * O(n^2) for two or more rules with many matches. */
static int scan_stable_sort(RegexSpan* a, size_t n) {
    if (n < 2) return 0;
    RegexSpan* tmp = malloc(n * sizeof(RegexSpan));
    if (!tmp) return -1;
    for (size_t width = 1; width < n; width *= 2) {
        for (size_t lo = 0; lo < n; lo += 2 * width) {
            size_t mid = lo + width < n ? lo + width : n;
            size_t hi = lo + 2 * width < n ? lo + 2 * width : n;
            size_t i = lo, j = mid, k = lo;
            while (i < mid && j < hi) tmp[k++] = (a[j].ms < a[i].ms) ? a[j++] : a[i++];
            while (i < mid) tmp[k++] = a[i++];
            while (j < hi)  tmp[k++] = a[j++];
        }
        memcpy(a, tmp, n * sizeof(RegexSpan));
    }
    free(tmp);
    return 0;
}

long regex_scan(const char* subj, size_t len, RegexRule* rules, int nr,
                RegexOverlapMode mode, int want_captures, RegexScan* out) {
    return regex_scan_limit(subj, len, rules, nr, mode, want_captures, 0, out);
}

long regex_scan_limit(const char* subj, size_t len, RegexRule* rules, int nr,
                      RegexOverlapMode mode, int want_captures, size_t limit,
                      RegexScan* out) {
    ScanState st;
    memset(&st, 0, sizeof st);
    st.want_captures = want_captures;
    /* Only the False scan emits in final order as it goes, so only it can stop
     * early; the others are cut after sorting. */
    st.limit = (mode == REGEX_OV_FALSE) ? limit : 0;

    out->spans = NULL;
    out->count = 0;
    out->caps  = NULL;

    int rc;
    if (mode == REGEX_OV_ALL)        rc = scan_all(subj, len, rules, nr, &st);
    else if (mode == REGEX_OV_FALSE) rc = scan_false(subj, len, rules, nr, &st);
    else                             rc = scan_true(subj, len, rules, nr, &st);

    /* The False and All scans already emit in ascending start order; only the
     * True scan concatenates one run per rule and needs merging. */
    if (rc == 0 && mode == REGEX_OV_TRUE && nr > 1)
        rc = scan_stable_sort(st.spans, st.count);

    if (rc != 0) {                  /* OOM / engine error: *out stays zeroed */
        free(st.spans);
        free(st.caps);
        return -1;
    }
    if (limit && st.count > limit) st.count = limit;

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
        /* Copy the literal run up to the next '$' in one append. */
        const char* d = strchr(tpl + i, '$');
        size_t run = d ? (size_t)(d - (tpl + i)) : strlen(tpl + i);
        if (regexbuf_add(&b, tpl + i, run)) goto oom;
        i += run;
    }
    if (!b.p) return rc_strdup("");                   /* empty template */
    return b.p;
oom:
    free(b.p);
    return NULL;
}

/* Expand $n in every string of an RHS tree, before it is evaluated.
 *
 * A RegularExpression rule's "$n" references are a property of the strings in
 * its RHS wherever they sit: {"$1", "$2"}, f["$1", {"$2"}] and
 * ToUpperCase["$1"] all see the capture -- the last because expansion precedes
 * evaluation. Copy-on-write: expr_copy is a refcount bump, so a subtree with no
 * '$' string is shared rather than rebuilt and the rule's RHS is never written
 * through. Returns an owned reference. */
static Expr* expand_rhs_tree(Expr* e, const char* subj,
                             const size_t* ov, size_t npairs) {
    if (e->type == EXPR_STRING) {
        if (!strchr(e->data.string, '$')) return expr_copy(e);
        char* t = regex_expand_template(e->data.string, subj, ov, npairs);
        if (!t) return expr_copy(e);
        Expr* s = expr_new_string(t);
        free(t);
        return s;
    }
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    size_t n = e->data.function.arg_count;
    Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
    if (!args) return expr_copy(e);
    Expr* head = expand_rhs_tree(e->data.function.head, subj, ov, npairs);
    int changed = (head != e->data.function.head);
    for (size_t i = 0; i < n; i++) {
        args[i] = expand_rhs_tree(e->data.function.args[i], subj, ov, npairs);
        if (args[i] != e->data.function.args[i]) changed = 1;
    }
    if (!changed) {
        for (size_t i = 0; i < n; i++) expr_free(args[i]);
        free(args);
        expr_free(head);
        return expr_copy(e);
    }
    Expr* out = expr_new_function(head, args, n);
    free(args);
    return out;
}

Expr* regex_rule_eval_rhs(const RegexRule* r, const char* subj,
                          const size_t* ov, size_t npairs) {
    if (!r->rhs) return NULL;

    if (r->rhs->type == EXPR_STRING) {
        if (!r->info.has_regex) return expr_copy(r->rhs);   /* literal text */
        char* t = regex_expand_template(r->rhs->data.string, subj, ov, npairs);
        if (!t) return expr_copy(r->rhs);
        Expr* e = expr_new_string(t);
        free(t);
        return e;
    }

    /* A regex rule's compound RHS: expand $n in each of its strings first. */
    Expr* rhs = r->info.has_regex ? expand_rhs_tree(r->rhs, subj, ov, npairs)
                                  : expr_copy(r->rhs);

    /* No named captures: nothing left to bind; just evaluate it. */
    if (r->info.nnamed == 0) return eval_and_free(rhs);

    MatchEnv* env = rule_bind_env(r, subj, ov, npairs, /*unset_empty=*/1);
    Expr* inst = replace_bindings(rhs, env);
    env_free(env);
    expr_free(rhs);
    return eval_and_free(inst);
}
