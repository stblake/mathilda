/*
 * string_pattern.c - translator from Wolfram string patterns to PCRE source.
 *
 * The regex string family (StringSplit / StringCases / StringReplace /
 * StringMatchQ / StringPosition) shares this one function to turn a Mathematica
 * string pattern argument into a PCRE regular expression. Supported forms:
 *
 *   "literal"                     literal text (PCRE-escaped)
 *   RegularExpression["re"]       raw PCRE source
 *   Whitespace                    \s+        (a run of whitespace)
 *   WhitespaceCharacter           \s
 *   LetterCharacter               [[:alpha:]]
 *   DigitCharacter                [[:digit:]]
 *   WordCharacter                 [[:alnum:]]
 *   NumberString                  [+-]?(?:\d+\.?\d*|\.\d+)
 *   StringExpression[a, b, ...]   concatenation  (a ~~ b ~~ ...)
 *   Alternatives[a, b, ...]       (?:a|b|...)    (a | b | ...)
 *   {a, b, ...}                    alternatives (a nested list in a pattern)
 *   Repeated[p]                   (?:p)+         (p ..)
 *   RepeatedNull[p]               (?:p)*         (p ...)
 *   Except[p]                     (?:(?!p).)     one char that does not start p
 *   Blank[]                       .
 *   BlankSequence[]               .+             (__)
 *   BlankNullSequence[]           .*             (___)
 *   Pattern[x, p]                 (p)            capture group (x: p); a repeat
 *                                   of an already-bound name x becomes a
 *                                   backreference \g{n} to the earlier group,
 *                                   so e.g. x_ ~~ x_ matches two equal characters
 *   PatternTest[Blank[], f]       class for a known predicate  (_?f)
 *                                   LetterQ/DigitQ/UpperCaseQ/LowerCaseQ
 *
 * Anything else yields NULL, so the calling builtin leaves the expression
 * unevaluated (the historical behaviour when only literals / RegularExpression
 * were understood).
 */

#include "regex_common.h"
#include "sym_names.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/*
 * Translation context: named capture groups are numbered in the order their
 * opening '(' appears in the emitted regex. names[g] holds the interned symbol
 * name bound to group g (1-based); a repeated name emits a backreference to its
 * group instead of a fresh capture. This lets x_ ~~ x_ mean "two equal chars".
 */
typedef struct {
    const char* names[REGEX_MAX_PAIRS + 1];   /* names[1..ngroups]; [0] unused */
    int         ngroups;
} SpCtx;

static char* sp_to_regex(Expr* patt, SpCtx* ctx);

/* C99-safe strdup replacement. */
static char* sp_strdup(const char* s) {
    size_t n = strlen(s) + 1;
    char* d = malloc(n);
    if (d) memcpy(d, s, n);
    return d;
}

/* Escape PCRE metacharacters in `s` so it matches literally. */
static char* sp_escape_literal(const char* s) {
    size_t n = strlen(s);
    char* out = malloc(2 * n + 1);
    if (!out) return NULL;
    size_t j = 0;
    for (size_t i = 0; i < n; i++) {
        char c = s[i];
        if (c && strchr("\\.^$|?*+()[]{}", c)) out[j++] = '\\';
        out[j++] = c;
    }
    out[j] = '\0';
    return out;
}

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

/*
 * Join child patterns into a single group. `alternation` selects `|` vs
 * concatenation. Each child is wrapped `(?:...)` and the whole is wrapped in a
 * non-capturing group so the result composes safely in any context. Returns
 * malloc'd source, or NULL if any child is unsupported.
 */
static char* group_join(Expr** args, size_t n, int alternation, SpCtx* ctx) {
    if (n == 0) return sp_strdup("");
    char** parts = calloc(n, sizeof(char*));
    if (!parts) return NULL;

    size_t total = 8;                 /* outer "(?:" + ")" + slack */
    int ok = 1;
    for (size_t i = 0; i < n; i++) {
        parts[i] = sp_to_regex(args[i], ctx);
        if (!parts[i]) { ok = 0; break; }
        total += strlen(parts[i]) + 6;   /* "(?:" + part + ")" + "|" */
    }

    char* out = NULL;
    if (ok) {
        out = malloc(total + 1);
        if (out) {
            char* w = out;
            memcpy(w, "(?:", 3); w += 3;
            for (size_t i = 0; i < n; i++) {
                if (i && alternation) *w++ = '|';
                w += sprintf(w, "(?:%s)", parts[i]);
            }
            *w++ = ')';
            *w = '\0';
        }
    }

    for (size_t i = 0; i < n; i++) free(parts[i]);
    free(parts);
    return out;
}

/* Wrap `inner` as pre + inner + post; frees `inner`. Returns malloc'd. */
static char* wrap_free(const char* pre, char* inner, const char* post) {
    if (!inner) return NULL;
    size_t lp = strlen(pre), li = strlen(inner), lq = strlen(post);
    char* out = malloc(lp + li + lq + 1);
    if (out) {
        memcpy(out, pre, lp);
        memcpy(out + lp, inner, li);
        memcpy(out + lp + li, post, lq);
        out[lp + li + lq] = '\0';
    }
    free(inner);
    return out;
}

/* Map a known predicate symbol name to a POSIX character class, else NULL. */
static const char* predicate_class(const char* name) {
    if (strcmp(name, "LetterQ") == 0)    return "[[:alpha:]]";
    if (strcmp(name, "DigitQ") == 0)     return "[[:digit:]]";
    if (strcmp(name, "UpperCaseQ") == 0) return "[[:upper:]]";
    if (strcmp(name, "LowerCaseQ") == 0) return "[[:lower:]]";
    return NULL;
}

/* Translate Pattern[x, p]: a fresh name becomes a capture group; a repeated
 * name becomes a backreference \g{n} to the group it was first bound to. */
static char* translate_pattern(Expr** args, SpCtx* ctx) {
    if (args[0]->type == EXPR_SYMBOL) {
        const char* nm = args[0]->data.symbol.name;
        for (int g = 1; g <= ctx->ngroups; g++) {
            if (ctx->names[g] == nm) {
                char buf[24];
                snprintf(buf, sizeof buf, "\\g{%d}", g);
                return sp_strdup(buf);
            }
        }
        if (ctx->ngroups < REGEX_MAX_PAIRS) {
            int g = ++ctx->ngroups;      /* assign before recursing: outer '(' first */
            ctx->names[g] = nm;
            return wrap_free("(", sp_to_regex(args[1], ctx), ")");
        }
    }
    /* Anonymous name, or capture-group budget exhausted: a plain (non-recorded)
     * capture so nested numbering still advances consistently. */
    if (ctx->ngroups < REGEX_MAX_PAIRS) ctx->ngroups++;
    return wrap_free("(", sp_to_regex(args[1], ctx), ")");
}

static char* sp_to_regex(Expr* patt, SpCtx* ctx) {
    if (!patt) return NULL;

    /* Literal string. */
    if (patt->type == EXPR_STRING)
        return sp_escape_literal(patt->data.string);

    /* Character-class heads (bare symbols). */
    if (patt->type == EXPR_SYMBOL) {
        if (patt->data.symbol.name == SYM_Whitespace)          return sp_strdup("\\s+");
        if (patt->data.symbol.name == SYM_WhitespaceCharacter) return sp_strdup("\\s");
        if (patt->data.symbol.name == SYM_LetterCharacter)     return sp_strdup("[[:alpha:]]");
        if (patt->data.symbol.name == SYM_DigitCharacter)      return sp_strdup("[[:digit:]]");
        if (patt->data.symbol.name == SYM_WordCharacter)       return sp_strdup("[[:alnum:]]");
        if (patt->data.symbol.name == SYM_NumberString)
            return sp_strdup("[+-]?(?:\\d+\\.?\\d*|\\.\\d+)");
        return NULL;
    }

    if (patt->type != EXPR_FUNCTION) return NULL;

    size_t argc;
    Expr** args;

    /* RegularExpression["re"] -> verbatim source. */
    if (is_fn(patt, SYM_RegularExpression, &argc, &args) &&
        argc == 1 && args[0]->type == EXPR_STRING)
        return sp_strdup(args[0]->data.string);

    /* StringExpression[a, b, ...] -> concatenation. */
    if (is_fn(patt, SYM_StringExpression, &argc, &args))
        return group_join(args, argc, /*alternation=*/0, ctx);

    /* Alternatives[a, b, ...] -> (?:a|b|...). */
    if (is_fn(patt, SYM_Alternatives, &argc, &args))
        return group_join(args, argc, /*alternation=*/1, ctx);

    /* A nested list in a pattern acts as alternatives (e.g. inside
     * StringExpression). Top-level lists are split into separate rules earlier,
     * so this only sees genuinely nested lists. */
    if (is_fn(patt, SYM_List, &argc, &args))
        return group_join(args, argc, /*alternation=*/1, ctx);

    /* Repeated[p] / RepeatedNull[p] (single-argument forms). */
    if (is_fn(patt, SYM_Repeated, &argc, &args) && argc == 1)
        return wrap_free("(?:", sp_to_regex(args[0], ctx), ")+");
    if (is_fn(patt, SYM_RepeatedNull, &argc, &args) && argc == 1)
        return wrap_free("(?:", sp_to_regex(args[0], ctx), ")*");

    /* Except[p] -> a single char that does not begin a p-match. */
    if (is_fn(patt, SYM_Except, &argc, &args) && argc == 1)
        return wrap_free("(?:(?!", sp_to_regex(args[0], ctx), ").)");

    /* Blank[] -> any character; __ / ___ -> one-or-more / zero-or-more. */
    if (is_fn(patt, SYM_Blank, &argc, &args) && argc == 0)
        return sp_strdup(".");
    if (is_fn(patt, SYM_BlankSequence, &argc, &args) && argc == 0)
        return sp_strdup(".+");
    if (is_fn(patt, SYM_BlankNullSequence, &argc, &args) && argc == 0)
        return sp_strdup(".*");

    /* Pattern[x, p] -> capture group / backreference. */
    if (is_fn(patt, SYM_Pattern, &argc, &args) && argc == 2)
        return translate_pattern(args, ctx);

    /* PatternTest[Blank[], f] (i.e. _?f) for a known predicate f. */
    if (is_fn(patt, SYM_PatternTest, &argc, &args) && argc == 2 &&
        is_fn(args[0], SYM_Blank, NULL, NULL) &&
        args[1]->type == EXPR_SYMBOL) {
        const char* cls = predicate_class(args[1]->data.symbol.name);
        if (cls) return sp_strdup(cls);
    }

    return NULL;
}

char* wl_pattern_to_regex(Expr* patt, int* is_null) {
    if (is_null) *is_null = 0;
    if (!patt) return NULL;

    /* The null delimiter "" (split at every character) is flagged for callers. */
    if (patt->type == EXPR_STRING && patt->data.string[0] == '\0' && is_null)
        *is_null = 1;

    SpCtx ctx;
    ctx.ngroups = 0;
    return sp_to_regex(patt, &ctx);
}
