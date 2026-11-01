/*
 * string_pattern.c - translator from Wolfram string patterns to PCRE source.
 *
 * The regex string family (StringSplit / StringCases / StringReplace /
 * StringMatchQ / StringPosition / StringCount / StringContainsQ / ...) shares
 * this one translator to turn a Mathematica string pattern into a PCRE2
 * regular expression. The emitted source is compiled with PCRE2_DOTALL (so the
 * blanks match newlines, as Mathematica's do) -- RegularExpression[...] pieces
 * are wrapped (?-s:...) so a user regex keeps PCRE's own defaults.
 *
 *   "literal"                     literal text (PCRE-escaped)
 *   RegularExpression["re"]       raw PCRE source
 *   Whitespace                    \s+        (a run of whitespace)
 *   WhitespaceCharacter           \s
 *   LetterCharacter               [[:alpha:]]
 *   DigitCharacter                [[:digit:]]
 *   WordCharacter                 [[:alnum:]]   (letters or digits)
 *   HexadecimalCharacter          [[:xdigit:]]
 *   PunctuationCharacter          [[:punct:]]
 *   NumberString                  [+-]?(?:\d+\.?\d*|\.\d+)
 *   CharacterRange["a", "z"]      [a-z]
 *   StartOfString / EndOfString   \A / \z
 *   StartOfLine / EndOfLine       (?:\A|(?<=\n)) / (?=\n|\z)
 *   WordBoundary                  \b
 *   StringExpression[a, b, ...]   concatenation  (a ~~ b ~~ ...)
 *   Alternatives[a, b, ...]       (?:a|b|...)    (a | b | ...)
 *   {a, b, ...}                   alternatives (a nested list in a pattern)
 *   Repeated[p]                   (?:p)+         (p ..)
 *   Repeated[p, n | {n} | {m, n}] (?:p){1,n} / {n} / {m,n}
 *   RepeatedNull[p, ...]          as Repeated, minimum 0  (p ...)
 *   Shortest[p] / Longest[p]      p with lazy / greedy quantifiers
 *   Except[p]                     (?:(?!p).)     one char that does not start p
 *   Except[p, q]                  (?:(?!p)q)     q, unless it starts a p
 *   Blank[]                       .              (_)
 *   BlankSequence[]               .+             (__)
 *   BlankNullSequence[]           .*             (___)
 *   Pattern[x, p]                 (p)            capture group (x: p)
 *   PatternTest[p, f]             (p)(?C<k>)     p, then a callout testing
 *                                   f[text of p] (a known predicate on a blank
 *                                   -- LetterQ, DigitQ, UpperCaseQ, LowerCaseQ
 *                                   -- is a plain character class instead)
 *   Condition[p, c]               (?:p)(?C<k>)   p, then a callout testing c
 *                                   with the names bound so far
 *
 * Pattern names. A name's first occurrence is a capture group; every later
 * occurrence is a backreference \g{n} to that group, so x_ ~~ x_ means "two
 * equal characters". This is exactly Mathematica's semantics, including its
 * consequences: in (x:"a") | (x:"b") the second x is a backreference to the
 * first branch's group, which is unset on that branch, so only "a" can match.
 * A name inside a repeated pattern binds once and every repetition must equal
 * it -- (x:LetterCharacter).. matches "aaa" but not "abc" -- which is emitted as
 * the body followed by a repetition of the body with its names as
 * backreferences.
 *
 * Callouts. Condition and a general PatternTest need Mathilda evaluation in the
 * middle of a match; PCRE2 callouts provide exactly that, and a callout that
 * returns nonzero makes the engine backtrack and try the next alternative, which
 * is Mathematica's semantics for a failed test. The translator records each
 * callout's expression in SpInfo; regex_common.c evaluates it.
 *
 * Anything else yields NULL, so the calling builtin leaves the expression
 * unevaluated.
 */

#include "regex_common.h"
#include "sym_names.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Largest number of callouts one pattern may carry (PCRE2 numbers 1..255). */
#define SP_MAX_CALLOUTS 255

/*
 * Translation context. Source is appended to `out`. Capture groups are numbered
 * in the order their opening '(' is emitted; names[g] is the interned pattern
 * name bound by group g (NULL for an anonymous group).
 */
typedef struct {
    RegexBuf    out;
    const char* names[REGEX_MAX_PAIRS];   /* names[1..ngroups]; [0] unused */
    int         ngroups;
    int         lazy;                     /* inside Shortest[...] */
    int         has_regex;
    int         failed;                   /* unsupported form / budget exceeded */
    int         wildcards;                /* StringMatchQ: `*`, `@` in literals */
    SpCallout*  callouts;
    int         ncallouts, callouts_cap;
} SpCtx;

static void sp_emit(SpCtx* ctx, Expr* patt);

/* ------------------------------------------------------------------ */
/* Output helpers                                                     */
/* ------------------------------------------------------------------ */

static void put(SpCtx* ctx, const char* s) {
    if (regexbuf_add(&ctx->out, s, strlen(s))) ctx->failed = 1;
}

/* A quantifier, made lazy inside Shortest[...]. */
static void put_q(SpCtx* ctx, const char* q) {
    put(ctx, q);
    if (ctx->lazy) put(ctx, "?");
}

/* {m,n} quantifier; n < 0 means unbounded. */
static void put_range(SpCtx* ctx, long m, long n) {
    char buf[64];
    if (n < 0)       snprintf(buf, sizeof buf, "{%ld,}", m);
    else if (m == n) snprintf(buf, sizeof buf, "{%ld}", m);
    else             snprintf(buf, sizeof buf, "{%ld,%ld}", m, n);
    put_q(ctx, buf);
}

/* Escape PCRE metacharacters in `s` so it matches literally. With wildcards
 * (StringMatchQ), `*` is any run of characters and `@` one or more characters
 * that are not uppercase letters; a backslash makes either literal. */
static void put_literal(SpCtx* ctx, const char* s) {
    for (const char* p = s; *p; p++) {
        if (ctx->wildcards) {
            if (*p == '*') { put(ctx, "."); put_q(ctx, "*"); continue; }
            if (*p == '@') { put(ctx, "[^[:upper:]]"); put_q(ctx, "+"); continue; }
            if (*p == '\\' && (p[1] == '*' || p[1] == '@')) {
                put(ctx, p[1] == '*' ? "\\*" : "@");
                p++;
                continue;
            }
        }
        char c[3] = { '\\', *p, '\0' };
        put(ctx, strchr("\\.^$|?*+()[]{}#", *p) ? c : c + 1);
    }
}

/* Open a fresh capture group bound to `name` (NULL: anonymous). */
static int open_group(SpCtx* ctx, const char* name) {
    if (ctx->ngroups + 1 >= REGEX_MAX_PAIRS) { ctx->failed = 1; return 0; }
    int g = ++ctx->ngroups;
    ctx->names[g] = name;
    put(ctx, "(");
    return g;
}

/* Emit a callout (?C<k>) for `e` of the given kind. */
static void put_callout(SpCtx* ctx, int kind, Expr* e, int group) {
    if (ctx->ncallouts >= SP_MAX_CALLOUTS) { ctx->failed = 1; return; }
    if (ctx->ncallouts == ctx->callouts_cap) {
        int nc = ctx->callouts_cap ? 2 * ctx->callouts_cap : 4;
        SpCallout* na = realloc(ctx->callouts, (size_t)nc * sizeof(SpCallout));
        if (!na) { ctx->failed = 1; return; }
        ctx->callouts = na;
        ctx->callouts_cap = nc;
    }
    SpCallout* c = &ctx->callouts[ctx->ncallouts++];
    c->kind  = kind;
    c->expr  = e;
    c->group = group;
    char buf[16];
    snprintf(buf, sizeof buf, "(?C%d)", ctx->ncallouts);
    put(ctx, buf);
}

/* ------------------------------------------------------------------ */
/* Pattern inspection                                                 */
/* ------------------------------------------------------------------ */

/* True if `e` is a function with head symbol `sym`; writes arg count/args. */
static int is_fn(const Expr* e, const char* sym, size_t* argc, Expr*** args) {
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == sym) {
        if (argc) *argc = e->data.function.arg_count;
        if (args) *args = e->data.function.args;
        return 1;
    }
    return 0;
}

/* A machine integer >= 0 out of `e`, or -1. */
static long nonneg_int(const Expr* e) {
    return (e->type == EXPR_INTEGER && e->data.integer >= 0) ? (long)e->data.integer : -1;
}

/* Character class for a known predicate symbol, else NULL. Each class must
 * accept exactly the characters its predicate accepts on a one-character
 * string (src/strings/stringcase.c), or `_?f` in a string pattern disagrees
 * with f itself. The engine compiles in PCRE2_UCP mode, where the POSIX
 * classes are Unicode properties: [[:alpha:]] = \p{L} (LetterQ: L*),
 * [[:upper:]] = \p{Lu}, [[:lower:]] = \p{Ll} -- those agree. But
 * [[:digit:]] = \p{Nd}, every Unicode decimal digit, while DigitQ is
 * ASCII-only (DigitQ["\:0661"] is False, as in Mathematica), so DigitQ maps
 * to the literal range [0-9]. */
static const char* predicate_class(const Expr* f) {
    if (f->type != EXPR_SYMBOL) return NULL;
    const char* name = f->data.symbol.name;
    if (name == SYM_LetterQ)    return "[[:alpha:]]";
    if (name == SYM_DigitQ)     return "[0-9]";
    if (name == SYM_UpperCaseQ) return "[[:upper:]]";
    if (name == SYM_LowerCaseQ) return "[[:lower:]]";
    return NULL;
}

/* Emit one character of a CharacterRange endpoint inside [...]. */
static int put_class_char(SpCtx* ctx, const Expr* e) {
    if (e->type != EXPR_STRING || e->data.string[0] == '\0') return 0;
    /* One character (one UTF-8 sequence). */
    const unsigned char* s = (const unsigned char*)e->data.string;
    size_t n = (s[0] < 0x80) ? 1 : (s[0] >= 0xF0) ? 4 : (s[0] >= 0xE0) ? 3 : 2;
    if (strlen(e->data.string) != n) return 0;
    if (n == 1 && strchr("\\]^-[", s[0])) put(ctx, "\\");
    put(ctx, e->data.string);
    return 1;
}

/* ------------------------------------------------------------------ */
/* Translation                                                        */
/* ------------------------------------------------------------------ */

/* Children joined by `sep` ("|" or ""), each in its own non-capturing group. */
static void emit_join(SpCtx* ctx, Expr** args, size_t n, const char* sep) {
    put(ctx, "(?:");
    for (size_t i = 0; i < n; i++) {
        if (i) put(ctx, sep);
        put(ctx, "(?:");
        sp_emit(ctx, args[i]);
        put(ctx, ")");
    }
    put(ctx, ")");
}

/* Pattern[x, p]: a fresh name opens a capture group; a name already bound
 * becomes a backreference to its group. */
static void emit_pattern(SpCtx* ctx, Expr** args) {
    const char* nm = args[0]->type == EXPR_SYMBOL ? args[0]->data.symbol.name : NULL;
    if (nm) {
        for (int g = 1; g <= ctx->ngroups; g++) {
            if (ctx->names[g] == nm) {
                char buf[24];
                snprintf(buf, sizeof buf, "\\g{%d}", g);
                put(ctx, buf);
                return;
            }
        }
    }
    open_group(ctx, nm);
    sp_emit(ctx, args[1]);
    put(ctx, ")");
}

/* Is `nm` already bound to a capture group? */
static int name_bound(const SpCtx* ctx, const char* nm) {
    for (int g = 1; g <= ctx->ngroups; g++)
        if (ctx->names[g] == nm) return 1;
    return 0;
}

/* `blank`?f for a predicate f with character class `cls`: _ / __ / ___ become
 * cls, cls+, cls* (f[s] of a string tests every character). Returns 0, emitting
 * nothing, when `blank` is not an argument-free blank. */
static int emit_class_blank(SpCtx* ctx, const Expr* blank, const char* cls) {
    if (blank->type != EXPR_FUNCTION || blank->data.function.arg_count != 0 ||
        blank->data.function.head->type != EXPR_SYMBOL)
        return 0;
    const char* h = blank->data.function.head->data.symbol.name;
    if (h == SYM_Blank)             { put(ctx, cls); return 1; }
    if (h == SYM_BlankSequence)     { put(ctx, cls); put_q(ctx, "+"); return 1; }
    if (h == SYM_BlankNullSequence) { put(ctx, cls); put_q(ctx, "*"); return 1; }
    return 0;
}

/* Did translating the region since `mark` open a named group? */
static int binds_names_since(const SpCtx* ctx, int mark) {
    for (int g = mark + 1; g <= ctx->ngroups; g++)
        if (ctx->names[g]) return 1;
    return 0;
}

/* Repeated[p, ...] (min_floor 1) / RepeatedNull[p, ...] (min_floor 0). */
static void emit_repeated(SpCtx* ctx, Expr** args, size_t argc, long min_floor) {
    long m = min_floor, n = -1;                /* default: {min_floor,} */
    if (argc == 2) {
        Expr* spec = args[1];
        size_t sc; Expr** sa;
        if (nonneg_int(spec) >= 0) {
            n = nonneg_int(spec);              /* Repeated[p, n]: up to n */
        } else if (is_fn(spec, SYM_List, &sc, &sa) && sc == 1 && nonneg_int(sa[0]) >= 0) {
            m = n = nonneg_int(sa[0]);         /* exactly n */
        } else if (is_fn(spec, SYM_List, &sc, &sa) && sc == 2 && nonneg_int(sa[0]) >= 0 &&
                   (nonneg_int(sa[1]) >= 0 ||
                    (sa[1]->type == EXPR_SYMBOL && sa[1]->data.symbol.name == SYM_Infinity))) {
            m = nonneg_int(sa[0]);
            n = nonneg_int(sa[1]);             /* -1 for Infinity */
        } else {
            ctx->failed = 1;
            return;
        }
        if (n >= 0 && n < m) { ctx->failed = 1; return; }
    } else if (argc != 1) {
        ctx->failed = 1;
        return;
    }

    /* Translate the body once to learn whether it binds names. If it does,
     * every repetition must equal the first, so the source is
     *   body (?:body′){m-1, n-1}
     * where body′ re-translates the body with those names now backreferences;
     * wrapped (?: ... )? when zero repetitions are allowed. */
    size_t start = ctx->out.len;
    int mark = ctx->ngroups;
    put(ctx, "(?:");
    sp_emit(ctx, args[0]);
    put(ctx, ")");
    if (ctx->failed) return;

    if (!binds_names_since(ctx, mark)) {
        if (m == 1 && n < 0)      put_q(ctx, "+");
        else if (m == 0 && n < 0) put_q(ctx, "*");
        else                      put_range(ctx, m, n);
        return;
    }

    if (n == 0) {                              /* zero repetitions only */
        ctx->out.len = start;
        ctx->out.p[start] = '\0';
        return;
    }
    /* Splice "(?:" in front when the whole thing is optional. */
    if (m == 0) {
        size_t body_len = ctx->out.len - start;
        char* body = malloc(body_len + 1);
        if (!body) { ctx->failed = 1; return; }
        memcpy(body, ctx->out.p + start, body_len + 1);
        ctx->out.len = start;
        put(ctx, "(?:");
        put(ctx, body);
        free(body);
    }
    if (n != 1) {
        put(ctx, "(?:");
        sp_emit(ctx, args[0]);
        put(ctx, ")");
        long rm = m > 0 ? m - 1 : 0, rn = n < 0 ? -1 : n - 1;
        if (rm == 0 && rn < 0) put_q(ctx, "*");
        else                   put_range(ctx, rm, rn);
    }
    if (m == 0) put_q(ctx, ")?");
}

static void sp_emit(SpCtx* ctx, Expr* patt) {
    if (ctx->failed || !patt) { ctx->failed = 1; return; }

    if (patt->type == EXPR_STRING) { put_literal(ctx, patt->data.string); return; }

    if (patt->type == EXPR_SYMBOL) {
        const char* s = patt->data.symbol.name;
        if (s == SYM_Whitespace)            { put(ctx, "\\s"); put_q(ctx, "+"); return; }
        if (s == SYM_WhitespaceCharacter)   { put(ctx, "\\s"); return; }
        if (s == SYM_LetterCharacter)       { put(ctx, "[[:alpha:]]"); return; }
        if (s == SYM_DigitCharacter)        { put(ctx, "[[:digit:]]"); return; }
        if (s == SYM_WordCharacter)         { put(ctx, "[[:alnum:]]"); return; }
        if (s == SYM_HexadecimalCharacter)  { put(ctx, "[[:xdigit:]]"); return; }
        if (s == SYM_PunctuationCharacter)  { put(ctx, "[[:punct:]]"); return; }
        if (s == SYM_NumberString)          { put(ctx, "[+-]?(?:\\d+\\.?\\d*|\\.\\d+)"); return; }
        /* Zero-width anchors. \A / \z are absolute (unlike ^ / $), so
         * StartOfString ~~ p pins p to offset 0 even in an unanchored search. */
        if (s == SYM_StartOfString)         { put(ctx, "\\A"); return; }
        if (s == SYM_EndOfString)           { put(ctx, "\\z"); return; }
        if (s == SYM_StartOfLine)           { put(ctx, "(?:\\A|(?<=\\n))"); return; }
        if (s == SYM_EndOfLine)             { put(ctx, "(?=\\n|\\z)"); return; }
        if (s == SYM_WordBoundary)          { put(ctx, "\\b"); return; }
        ctx->failed = 1;
        return;
    }

    if (patt->type != EXPR_FUNCTION) { ctx->failed = 1; return; }

    size_t argc;
    Expr** args;

    /* RegularExpression["re"]: verbatim, under PCRE's default (non-DOTALL)
     * newline handling. */
    if (is_fn(patt, SYM_RegularExpression, &argc, &args)) {
        if (argc != 1 || args[0]->type != EXPR_STRING) { ctx->failed = 1; return; }
        ctx->has_regex = 1;
        put(ctx, "(?-s:");
        put(ctx, args[0]->data.string);
        put(ctx, ")");
        return;
    }

    if (is_fn(patt, SYM_StringExpression, &argc, &args)) { emit_join(ctx, args, argc, ""); return; }
    if (is_fn(patt, SYM_Alternatives, &argc, &args)) {
        if (argc == 0) { put(ctx, "(*FAIL)"); return; }   /* matches nothing */
        emit_join(ctx, args, argc, "|");
        return;
    }
    /* A nested list in a pattern acts as alternatives. Top-level lists are
     * split into separate rules before translation. */
    if (is_fn(patt, SYM_List, &argc, &args)) {
        if (argc == 0) { put(ctx, "(*FAIL)"); return; }
        emit_join(ctx, args, argc, "|");
        return;
    }

    if (is_fn(patt, SYM_Repeated, &argc, &args))     { emit_repeated(ctx, args, argc, 1); return; }
    if (is_fn(patt, SYM_RepeatedNull, &argc, &args)) { emit_repeated(ctx, args, argc, 0); return; }

    if ((is_fn(patt, SYM_Shortest, &argc, &args) || is_fn(patt, SYM_Longest, &argc, &args))
        && argc == 1) {
        int saved = ctx->lazy;
        ctx->lazy = (patt->data.function.head->data.symbol.name == SYM_Shortest);
        put(ctx, "(?:");
        sp_emit(ctx, args[0]);
        put(ctx, ")");
        ctx->lazy = saved;
        return;
    }

    /* Except[p]: one character that does not begin a p-match.
     * Except[p, q]: a q-match that does not begin with a p-match. */
    if (is_fn(patt, SYM_Except, &argc, &args) && (argc == 1 || argc == 2)) {
        put(ctx, "(?:(?!");
        sp_emit(ctx, args[0]);
        put(ctx, ")");
        if (argc == 2) { put(ctx, "(?:"); sp_emit(ctx, args[1]); put(ctx, ")"); }
        else           put(ctx, ".");
        put(ctx, ")");
        return;
    }

    if (is_fn(patt, SYM_Blank, &argc, &args) && argc == 0)             { put(ctx, "."); return; }
    if (is_fn(patt, SYM_BlankSequence, &argc, &args) && argc == 0)     { put(ctx, "."); put_q(ctx, "+"); return; }
    if (is_fn(patt, SYM_BlankNullSequence, &argc, &args) && argc == 0) { put(ctx, "."); put_q(ctx, "*"); return; }

    if (is_fn(patt, SYM_Pattern, &argc, &args) && argc == 2) { emit_pattern(ctx, args); return; }

    /* CharacterRange["a", "z"] -> [a-z]. */
    if (is_fn(patt, SYM_CharacterRange, &argc, &args) && argc == 2) {
        put(ctx, "[");
        if (!put_class_char(ctx, args[0])) { ctx->failed = 1; return; }
        put(ctx, "-");
        if (!put_class_char(ctx, args[1])) { ctx->failed = 1; return; }
        put(ctx, "]");
        return;
    }

    /* PatternTest[p, f] (p?f). A known character predicate on a blank is a
     * class (f[s] for a string s tests every character), otherwise the text p
     * matched is captured and handed to f by a callout. */
    if (is_fn(patt, SYM_PatternTest, &argc, &args) && argc == 2) {
        const char* cls = predicate_class(args[1]);
        if (cls && emit_class_blank(ctx, args[0], cls)) return;
        /* x__?LetterQ parses as PatternTest[Pattern[x, __], LetterQ]: with a
         * class predicate that is the capture x:(class+) -- unless x is
         * already bound, when the backreference must still pass the test. */
        size_t pc; Expr** pa;
        if (cls && is_fn(args[0], SYM_Pattern, &pc, &pa) && pc == 2 &&
            pa[0]->type == EXPR_SYMBOL && !name_bound(ctx, pa[0]->data.symbol.name)) {
            size_t mark = ctx->out.len;
            int g = open_group(ctx, pa[0]->data.symbol.name);
            if (emit_class_blank(ctx, pa[1], cls)) { put(ctx, ")"); return; }
            ctx->out.len = mark;               /* not a blank: general path */
            ctx->out.p[mark] = '\0';
            ctx->ngroups = g - 1;
        }
        int g = open_group(ctx, NULL);
        sp_emit(ctx, args[0]);
        put(ctx, ")");
        put_callout(ctx, SP_CALLOUT_TEST, args[1], g);
        return;
    }

    /* Condition[p, c] (p /; c): p, then c evaluated with the names bound so far. */
    if (is_fn(patt, SYM_Condition, &argc, &args) && argc == 2) {
        put(ctx, "(?:");
        sp_emit(ctx, args[0]);
        put(ctx, ")");
        put_callout(ctx, SP_CALLOUT_CONDITION, args[1], 0);
        return;
    }

    ctx->failed = 1;
}

/* ------------------------------------------------------------------ */
/* Public entry points                                                */
/* ------------------------------------------------------------------ */

char* wl_pattern_to_regex_ex(Expr* patt, SpInfo* info) {
    return wl_pattern_to_regex_wild(patt, info, 0);
}

char* wl_pattern_to_regex_wild(Expr* patt, SpInfo* info, int wildcards) {
    if (info) memset(info, 0, sizeof *info);
    if (!patt) return NULL;

    SpCtx ctx;
    memset(&ctx, 0, sizeof ctx);
    ctx.wildcards = wildcards;

    /* A bare RegularExpression is passed through untouched (no (?-s:) wrapper),
     * so leading PCRE verbs such as (*UCP) keep working. */
    size_t argc;
    Expr** args;
    int raw = is_fn(patt, SYM_RegularExpression, &argc, &args) &&
              argc == 1 && args[0]->type == EXPR_STRING;
    if (raw) {
        ctx.has_regex = 1;
        put(&ctx, args[0]->data.string);
    } else {
        sp_emit(&ctx, patt);
    }
    if (!ctx.failed && !ctx.out.p) put(&ctx, "");     /* "" -> empty source */

    if (ctx.failed) {
        free(ctx.out.p);
        free(ctx.callouts);
        return NULL;
    }

    if (!info) { free(ctx.callouts); return ctx.out.p; }

    info->is_null   = (patt->type == EXPR_STRING && patt->data.string[0] == '\0');
    info->has_regex = ctx.has_regex;
    info->raw_regex = raw;
    info->ngroups   = ctx.ngroups;
    info->callouts  = ctx.callouts;
    info->ncallouts = ctx.ncallouts;
    if (ctx.ngroups > 0) {
        info->group_name = calloc((size_t)ctx.ngroups + 1, sizeof(const char*));
        if (!info->group_name) { free(ctx.out.p); sp_info_free(info); return NULL; }
        for (int g = 1; g <= ctx.ngroups; g++) {
            info->group_name[g] = ctx.names[g];
            if (ctx.names[g]) info->nnamed++;
        }
    }
    return ctx.out.p;
}

void sp_info_free(SpInfo* info) {
    if (!info) return;
    free(info->group_name);
    free(info->callouts);
    info->group_name = NULL;
    info->callouts = NULL;
    info->ngroups = info->ncallouts = info->nnamed = 0;
}

char* wl_pattern_to_regex(Expr* patt, int* is_null) {
    if (is_null) *is_null = (patt && patt->type == EXPR_STRING && patt->data.string[0] == '\0');
    return wl_pattern_to_regex_ex(patt, NULL);
}
