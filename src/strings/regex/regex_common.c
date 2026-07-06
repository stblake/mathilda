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

/* Escape PCRE metacharacters in `s` so it matches literally. */
static char* escape_literal(const char* s) {
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

/* Regex source (malloc'd) for a single pattern element that is either
 * RegularExpression["re"] or a literal string. Returns NULL if unsupported. */
static char* pattern_source(Expr* e) {
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_RegularExpression &&
        e->data.function.arg_count == 1 &&
        e->data.function.args[0]->type == EXPR_STRING) {
        return rc_strdup(e->data.function.args[0]->data.string);
    }
    if (e->type == EXPR_STRING) {
        return escape_literal(e->data.string);
    }
    return NULL;
}

/* Is e a Rule[...] or RuleDelayed[...] with two arguments? */
static int is_rule2(Expr* e, Expr** lhs, Expr** rhs) {
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        (e->data.function.head->data.symbol == SYM_Rule ||
         e->data.function.head->data.symbol == SYM_RuleDelayed) &&
        e->data.function.arg_count == 2) {
        *lhs = e->data.function.args[0];
        *rhs = e->data.function.args[1];
        return 1;
    }
    return 0;
}

/* Compile one element into `out`. Returns 0 on success, -1 on failure. */
static int build_one(Expr* e, int anchored, RegexRule* out, const char* head) {
    Expr* patt = e;
    Expr* rhs = NULL;
    Expr* l; Expr* r;
    if (is_rule2(e, &l, &r)) { patt = l; rhs = r; }

    char* src = pattern_source(patt);
    if (!src) return -1;                    /* unsupported pattern -> unevaluated */

    char* final;
    if (anchored) {
        final = wrap_anchored(src);
        free(src);
        if (!final) return -1;
    } else {
        final = src;
    }

    char err[256];
    RegexProgram* prog = regex_compile(final, err, sizeof err);
    free(final);
    if (!prog) {
        fprintf(stderr, "%s::regex: %s\n", head, err);
        return -1;
    }
    out->prog = prog;
    out->rhs = rhs;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Public: build / free a rule set                                    */
/* ------------------------------------------------------------------ */

int regex_rules_build(Expr* patt, int anchored, RegexRule** out, const char* head) {
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
        patt->data.function.head->data.symbol == SYM_List) {
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
        if (build_one(elems[i], anchored, &rules[i], head) != 0) {
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
