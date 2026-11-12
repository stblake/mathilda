/*
 * term_highlight.c -- see term_highlight.h.
 */
#include "term_highlight.h"
#include "term_lex.h"
#include "term_color.h"

#include "symtab.h"

#include <stdlib.h>
#include <string.h>

/* 16-colour SGR only, so every terminal that is not "dumb" shows them. Each
 * entry starts from a reset state (the renderer emits "\033[0m" first). */
static const char* const HL_SGR[TERM_HL_COUNT] = {
    "",                 /* PLAIN                                   */
    "\033[36m",         /* BUILTIN  cyan                           */
    "\033[34m",         /* USER     blue                           */
    "\033[35m",         /* NUMBER   magenta                        */
    "\033[32m",         /* STRING   green                          */
    "\033[90m",         /* COMMENT  bright black                   */
    "\033[33m",         /* PATTERN  yellow                         */
    "",                 /* BRACKET  plain                          */
    "\033[1;4m",        /* BRACKET_MATCH  bold underline           */
    "\033[1;31m",       /* BRACKET_BAD    bold red                 */
    NULL,               /* PROMPT: term_color's In style           */
    "\033[97;100m",     /* POPUP      white on grey                */
    "\033[1;97;44m",    /* POPUP_SEL  bold white on blue           */
    "\033[37;100m",     /* POPUP_DOC  light grey on grey           */
    "\033[2m",          /* DIM                                     */
};

const char* term_hl_sgr(int id) {
    if (id <= 0 || id >= (int)TERM_HL_COUNT) return "";
    if (id == TERM_HL_PROMPT) return term_style_on(TERM_STYLE_IN, true);
    return HL_SGR[id];
}

int term_hl_classify_symtab(const char* s, size_t len, void* ctx) {
    (void)ctx;
    char name[256];
    if (len == 0 || len >= sizeof(name)) return TERM_HL_PLAIN;
    memcpy(name, s, len);
    name[len] = '\0';
    SymbolDef* d = symtab_lookup(name);
    if (!d) return TERM_HL_PLAIN;
    if (d->is_system) return TERM_HL_BUILTIN;
    if (d->own_values || d->down_values || d->sub_values || d->builtin_func
#if UP_VALUES
        || d->up_values
#endif
        ) return TERM_HL_USER;
    return TERM_HL_PLAIN;
}

const char* term_hl_doc_symtab(const char* s, size_t len, void* ctx) {
    (void)ctx;
    char name[256];
    if (len == 0 || len >= sizeof(name)) return NULL;
    memcpy(name, s, len);
    name[len] = '\0';
    SymbolDef* d = symtab_lookup(name);
    if (!d || !d->docstring) return NULL;
    const char* doc = d->docstring;
    while (*doc == ' ' || *doc == '\n') doc++;
    return *doc ? doc : NULL;
}

static bool is_digit(char c) { return c >= '0' && c <= '9'; }
static bool is_alnum(char c) {
    return is_digit(c) || (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

/* Fill out[a, b) with style v. */
static void fill(uint8_t* out, size_t a, size_t b, int v) {
    memset(out + a, v, b - a);
}

TermHlState term_hl_line(const char* s, size_t n, TermHlState st, uint8_t* out,
                         TermHlClassify cls, void* ctx) {
    size_t i = 0;
    st.invalid = 0;
    while (i < n) {
        size_t start = i;
        if (st.in_string) {
            while (i < n) {
                if (s[i] == '\\') { i += (i + 1 < n) ? 2 : 1; continue; }
                if (s[i] == '"') { i++; st.in_string = 0; break; }
                i++;
            }
            fill(out, start, i, TERM_HL_STRING);
            continue;
        }
        if (st.comment > 0) {
            while (i < n && st.comment > 0) {
                if (s[i] == '(' && i + 1 < n && s[i + 1] == '*') { st.comment++; i += 2; }
                else if (s[i] == '*' && i + 1 < n && s[i + 1] == ')') { st.comment--; i += 2; }
                else i++;
            }
            fill(out, start, i, TERM_HL_COMMENT);
            continue;
        }
        char c = s[i];
        if (c == '(' && i + 1 < n && s[i + 1] == '*') {
            st.comment = 1;
            i += 2;
            fill(out, start, i, TERM_HL_COMMENT);
            continue;
        }
        if (c == '"') {
            st.in_string = 1;
            i++;
            fill(out, start, i, TERM_HL_STRING);
            continue;
        }
        if (is_digit(c) || (c == '.' && i + 1 < n && is_digit(s[i + 1])
                             && !(i > 0 && s[i - 1] == '.'))) {
            while (i < n && is_digit(s[i])) i++;
            if (i + 1 < n && s[i] == '^' && s[i + 1] == '^') {      /* 16^^ff */
                i += 2;
                while (i < n && (is_alnum(s[i]) || s[i] == '.')) i++;
            } else {
                if (i < n && s[i] == '.' && !(i + 1 < n && s[i + 1] == '.')) {
                    i++;
                    while (i < n && is_digit(s[i])) i++;
                }
                if (i < n && s[i] == '`') {                           /* 1.5`20 */
                    i++;
                    if (i < n && s[i] == '`') i++;
                    while (i < n && (is_digit(s[i]) || s[i] == '.')) i++;
                }
                if (i + 1 < n && s[i] == '*' && s[i + 1] == '^') {   /* 1*^-3 */
                    size_t j = i + 2;
                    if (j < n && (s[j] == '-' || s[j] == '+')) j++;
                    if (j < n && is_digit(s[j])) {
                        i = j;
                        while (i < n && is_digit(s[i])) i++;
                    }
                }
            }
            fill(out, start, i, TERM_HL_NUMBER);
            continue;
        }
        if (c == '#') {                                               /* slots */
            i++;
            if (i < n && s[i] == '#') i++;
            while (i < n && is_digit(s[i])) i++;
            fill(out, start, i, TERM_HL_PATTERN);
            continue;
        }
        if (c == '_' || (term_lex_is_word_char(c) && !is_digit(c))) {
            while (i < n && term_lex_is_word_char(s[i])) i++;
            size_t name_end = i;
            if (i < n && s[i] == '_') {                               /* a pattern */
                int k = 0;
                while (i < n && s[i] == '_' && k < 3) { i++; k++; }
                while (i < n && term_lex_is_word_char(s[i])) i++;    /* _Integer */
                if (i < n && s[i] == '.' && !(i + 1 < n && is_digit(s[i + 1]))) i++;
                fill(out, start, i, TERM_HL_PATTERN);
                continue;
            }
            int v = cls ? cls(s + start, name_end - start, ctx) : TERM_HL_PLAIN;
            fill(out, start, i, v);
            continue;
        }
        if (c == '(' || c == ')' || c == '[' || c == ']' || c == '{' || c == '}') {
            out[i++] = TERM_HL_BRACKET;
            continue;
        }
        out[i++] = TERM_HL_PLAIN;
    }
    return st;
}

bool term_hl_fputs(const char* s, size_t n, FILE* f, TermHlClassify cls, void* ctx) {
    uint8_t* sty = malloc(n ? n : 1);
    if (!sty) return false;
    TermHlState st = term_hl_state_initial();
    size_t i = 0;
    while (i < n) {
        const char* nl = memchr(s + i, '\n', n - i);
        size_t e = nl ? (size_t)(nl - s) : n;
        st = term_hl_line(s + i, e - i, st, sty + i, cls, ctx);
        if (e < n) sty[e] = TERM_HL_PLAIN;
        i = e + 1;
    }
    int cur = TERM_HL_PLAIN;
    size_t run = 0;
    for (i = 0; i < n; i++) {
        int v = sty[i];
        if (s[i] == '\n') v = TERM_HL_PLAIN;
        if (v == TERM_HL_BRACKET) v = TERM_HL_PLAIN;
        if (v != cur) {
            fwrite(s + run, 1, i - run, f);
            run = i;
            fputs("\033[0m", f);
            fputs(term_hl_sgr(v), f);
            cur = v;
        }
    }
    fwrite(s + run, 1, n - run, f);
    if (cur != TERM_HL_PLAIN) fputs("\033[0m", f);
    free(sty);
    return true;
}
