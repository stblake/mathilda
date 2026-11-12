/*
 * term_complete.c -- see term_complete.h.
 */
#include "term_complete.h"
#include "term_lex.h"

#include "expr.h"
#include "names.h"
#include "symtab.h"
#include "sym_names.h"

#include <stdlib.h>
#include <string.h>

#define TERM_COMPLETE_MAX_HEAD 256

typedef struct {
    char** items;
    size_t count, cap;
} StrList;

static void list_push(StrList* l, char* owned) {
    if (!owned) return;
    if (l->count == l->cap) {
        size_t nc = l->cap ? l->cap * 2 : 16;
        char** ni = realloc(l->items, nc * sizeof(char*));
        if (!ni) { free(owned); return; }
        l->items = ni;
        l->cap = nc;
    }
    l->items[l->count++] = owned;
}

/* name + suffix as a new string. */
static char* concat(const char* name, const char* suffix) {
    size_t a = strlen(name), b = strlen(suffix);
    char* r = malloc(a + b + 1);
    if (!r) return NULL;
    memcpy(r, name, a);
    memcpy(r + a, suffix, b + 1);
    return r;
}

static bool is_symbol_token(const char* s, size_t n) {
    if (n == 0 || (s[0] >= '0' && s[0] <= '9')) return false;
    for (size_t i = 0; i < n; i++) if (!term_lex_is_word_char(s[i])) return false;
    return true;
}

/* The registered options of `head` (a List of Rule / RuleDelayed), read
 * straight from the definition: symtab_get_options would note a dependency
 * in the middle of whatever is evaluating (a Dialog). */
static const Expr* head_options(const char* buf, const TermLexContext* ctx) {
    if (!ctx->has_call_head || ctx->head_len >= TERM_COMPLETE_MAX_HEAD) return NULL;
    if (!is_symbol_token(buf + ctx->head_start, ctx->head_len)) return NULL;
    char head[TERM_COMPLETE_MAX_HEAD];
    memcpy(head, buf + ctx->head_start, ctx->head_len);
    head[ctx->head_len] = '\0';
    SymbolDef* d = symtab_lookup(head);
    const Expr* o = d ? d->default_options : NULL;
    if (!o || o->type != EXPR_FUNCTION || o->data.function.head->type != EXPR_SYMBOL
        || o->data.function.head->data.symbol.name != SYM_List) return NULL;
    return o;
}

/* Option names of `opts` starting with word[0..wlen): symbol names when
 * !strings, string names when strings. */
static void add_options(StrList* l, const Expr* opts, const char* word, size_t wlen,
                        bool strings) {
    for (size_t i = 0; i < opts->data.function.arg_count; i++) {
        const Expr* r = opts->data.function.args[i];
        if (r->type != EXPR_FUNCTION || r->data.function.arg_count != 2
            || r->data.function.head->type != EXPR_SYMBOL) continue;
        const char* rh = r->data.function.head->data.symbol.name;
        if (rh != SYM_Rule && rh != SYM_RuleDelayed) continue;
        const Expr* lhs = r->data.function.args[0];
        const char* name = NULL;
        if (!strings && lhs->type == EXPR_SYMBOL) name = lhs->data.symbol.name;
        if (strings && lhs->type == EXPR_STRING) name = lhs->data.string;
        if (!name || strncmp(name, word, wlen) != 0) continue;
        list_push(l, concat(name, strings ? "\" -> " : " -> "));
    }
}

void term_complete(const char* buf, size_t cursor, TermCompletions* out) {
    memset(out, 0, sizeof(*out));
    if (!buf) return;
    TermLexContext ctx;
    term_lex_context(buf, cursor, &ctx);
    size_t n = strlen(buf);
    if (cursor > n) cursor = n;
    out->word_start = ctx.word_start;
    if (ctx.in_comment) return;

    size_t wlen = cursor - ctx.word_start;
    char* word = malloc(wlen + 1);
    if (!word) return;
    memcpy(word, buf + ctx.word_start, wlen);
    word[wlen] = '\0';

    StrList l = { NULL, 0, 0 };
    if (ctx.in_string) {
        /* Only a string that opens an argument: f["Ve|  */
        if (ctx.arg_position && ctx.string_start + 1 == ctx.word_start) {
            const Expr* opts = head_options(buf, &ctx);
            if (opts) add_options(&l, opts, word, wlen, true);
        }
    } else {
        if (ctx.arg_position) {
            const Expr* opts = head_options(buf, &ctx);
            if (opts) add_options(&l, opts, word, wlen, false);
        }
        if (wlen > 0 && !(word[0] >= '0' && word[0] <= '9')) {
            size_t n_opts = l.count;
            size_t k = 0;
            char** names = names_complete(word, &k);
            for (size_t i = 0; i < k; i++) {
                /* The option already offers this name, as "Name -> ". */
                bool dup = false;
                size_t nl = strlen(names[i]);
                for (size_t j = 0; j < n_opts && !dup; j++)
                    dup = strncmp(l.items[j], names[i], nl) == 0
                       && strcmp(l.items[j] + nl, " -> ") == 0;
                if (dup) free(names[i]);
                else list_push(&l, names[i]);
            }
            free(names);   /* the array only: the strings moved into l */
        }
    }
    free(word);
    out->items = l.items;
    out->count = l.count;
}

void term_completions_free(TermCompletions* c) {
    if (!c) return;
    for (size_t i = 0; i < c->count; i++) free(c->items[i]);
    free(c->items);
    c->items = NULL;
    c->count = 0;
}
