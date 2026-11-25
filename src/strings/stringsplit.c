/*
 * stringsplit.c - StringSplit[...], the full Wolfram surface.
 *
 *   StringSplit[s]                 split at runs of whitespace
 *   StringSplit[s, patt]           split at delimiters matching patt
 *   StringSplit[s, {p1, p2, ...}]  split at any of the pi
 *   StringSplit[s, patt -> val]    insert val at each delimiter
 *   StringSplit[s, patt, n]        at most n substrings
 *   StringSplit[s, patt, All]      keep leading/trailing empty substrings
 *   StringSplit[{s1, ...}, patt]   thread over a list of subjects
 *   IgnoreCase -> True             case-insensitive delimiters
 *
 * The delimiter pattern is translated to PCRE by the shared string-pattern
 * engine (regex_common.c / string_pattern.c), so it accepts literal strings,
 * RegularExpression["re"], the character-class heads (Whitespace, ...),
 * StringExpression (~~), Alternatives (|), Repeated (..), Except, and so on.
 *
 * Zero-length substrings between two adjacent interior delimiters are kept;
 * empty substrings at the very beginning or end are dropped unless All is given.
 * The empty-string delimiter "" splits at every character.
 */

#include "picostrings.h"
#include "regex_common.h"
#include "linkedlist.h"   /* a LinkedList of subjects */
#include "sym_names.h"
#include "common.h"
#include "eval.h"
#include "symtab.h"
#include "utf8.h"
#include "strings_internal.h"   /* strings_msg_strse */

#include <stdlib.h>
#include <string.h>

/* A result element: an Expr plus a flag marking it an (empty-trimmable)
 * substring piece rather than an inserted rule value. */
typedef struct { Expr* e; int is_sub; } Piece;
typedef struct { Piece* v; size_t n, cap; } PieceVec;

/* Append; takes ownership of e and frees it on OOM. Returns 0, or -1 on OOM. */
static int pv_push(PieceVec* pv, Expr* e, int is_sub) {
    if (pv->n == pv->cap) {
        size_t nc = pv->cap ? pv->cap * 2 : 8;
        Piece* nv = realloc(pv->v, nc * sizeof(Piece));
        if (!nv) { expr_free(e); return -1; }
        pv->v = nv;
        pv->cap = nc;
    }
    pv->v[pv->n].e = e;
    pv->v[pv->n].is_sub = is_sub;
    pv->n++;
    return 0;
}

/* New EXPR_STRING from subj[s..e). */
static Expr* substr_expr(const char* subj, size_t s, size_t e) {
    return expr_new_string_n(subj + s, e - s);
}

static int is_empty_string(const Expr* e) {
    return e->type == EXPR_STRING && e->data.string[0] == '\0';
}

/* The value inserted at a delimiter matched by `rule` (its RHS instantiated
 * with the delimiter's named pieces, as StringCases does), or NULL to insert
 * nothing (a bare delimiter pattern). */
static Expr* make_insertion(const RegexRule* rule, const char* subj,
                            const size_t* ov, size_t pairs) {
    return regex_rule_eval_rhs(rule, subj, ov, pairs);
}

/* Split one string subject. */
static Expr* ss_scalar(const char* subj, RegexRule* rules, int nr,
                       long n_limit, int keep_ends) {
    size_t len = strlen(subj);
    PieceVec pv = {0};

    size_t last_end = 0, pos = 0;
    long ns = 0;                        /* substrings emitted so far */
    RegexCursor cur;
    if (regex_cursor_init(&cur, rules, nr, subj, len) != 0) return NULL;

    while (pos <= len) {
        if (n_limit > 0 && ns >= n_limit - 1) break;   /* remainder is last piece */

        int rule;
        const size_t* ov;
        size_t np;
        int m = regex_cursor_next(&cur, pos, &rule, &ov, &np);
        if (m <= 0) break;
        size_t best_ms = ov[0], best_me = ov[1];

        pv_push(&pv, substr_expr(subj, last_end, best_ms), 1);   /* interior empties kept */
        ns++;

        Expr* ins = make_insertion(&rules[rule], subj, ov, np);
        if (ins) pv_push(&pv, ins, 0);

        last_end = best_me;
        pos = (best_me > best_ms) ? best_me
                                  : regex_next_char(subj, len, best_ms); /* zero-width guard */
    }

    regex_cursor_free(&cur);
    pv_push(&pv, substr_expr(subj, last_end, len), 1);           /* trailing piece */

    /* Drop the leading and trailing runs of empty substring pieces unless All. */
    size_t lo = 0, hi = pv.n;
    if (!keep_ends) {
        while (lo < hi && pv.v[lo].is_sub && is_empty_string(pv.v[lo].e)) lo++;
        while (hi > lo && pv.v[hi - 1].is_sub && is_empty_string(pv.v[hi - 1].e)) hi--;
    }

    size_t out_n = hi - lo;
    Expr** items = malloc(sizeof(Expr*) * (out_n ? out_n : 1));
    size_t k = 0;
    for (size_t i = 0; i < pv.n; i++) {
        if (items && i >= lo && i < hi) items[k++] = pv.v[i].e;
        else expr_free(pv.v[i].e);
    }
    free(pv.v);
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), items, out_n);
    free(items);
    return result;
}

/* Null delimiter "": split into individual characters (UTF-8 code points),
 * at most n_limit pieces when n_limit > 0 (the last one takes the rest:
 * StringSplit["abc", "", 2] is {"a", "bc"}). */
static Expr* ss_null(const char* subj, long n_limit) {
    size_t len = strlen(subj);
    Expr** items = malloc(sizeof(Expr*) * (len ? len : 1));
    if (!items) return NULL;
    size_t n = 0;
    for (size_t i = 0; i < len; ) {                /* one piece per character */
        if (n_limit > 0 && n + 1 == (size_t)n_limit) {
            items[n++] = expr_new_string_n(subj + i, len - i);
            break;
        }
        size_t k = u8_char_len(subj, len, i);
        items[n++] = expr_new_string_n(subj + i, k);
        i += k;
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), items, n);
    free(items);
    return result;
}

/* Is e a supported IgnoreCase option (Rule/RuleDelayed[IgnoreCase, bool])? */
static int is_ignorecase_opt(const Expr* e, int* value) {
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        (e->data.function.head->data.symbol.name == SYM_Rule ||
         e->data.function.head->data.symbol.name == SYM_RuleDelayed) &&
        e->data.function.arg_count == 2 &&
        e->data.function.args[0]->type == EXPR_SYMBOL &&
        e->data.function.args[0]->data.symbol.name == SYM_IgnoreCase) {
        Expr* v = e->data.function.args[1];
        *value = (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_True);
        return 1;
    }
    return 0;
}

Expr* builtin_stringsplit(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;

    if (argc == 0) return builtin_arg_error("StringSplit", 0, 1, 3);

    /* Strip trailing IgnoreCase options, leaving the positional arguments. */
    /* IgnoreCase: the registered default (SetOptions[StringSplit, ...]), then
     * a trailing option. */
    int caseless = 0;
    Expr* defs = symtab_get_options("StringSplit");   /* borrowed */
    if (defs && defs->type == EXPR_FUNCTION)
        for (size_t i = 0; i < defs->data.function.arg_count; i++)
            is_ignorecase_opt(defs->data.function.args[i], &caseless);
    size_t pargc = argc;
    while (pargc >= 2) {
        int v;
        if (is_ignorecase_opt(a[pargc - 1], &v)) { caseless = v; pargc--; }
        else break;
    }
    if (pargc < 1 || pargc > 3) return builtin_arg_error("StringSplit", argc, 1, 3);

    Expr* subject = a[0];
    /* A string or a list of strings, else StringSplit::strse (unevaluated). */
    if (subject->type != EXPR_STRING &&
        !(strings_is_seq(subject) && strings_seq_all_strings(subject)))
        return strings_msg_strse("StringSplit", res, 1);

    /* Third positional argument: an integer count n, or All. */
    long n_limit = 0;
    int keep_ends = 0;
    if (pargc == 3) {
        Expr* third = a[2];
        if (third->type == EXPR_SYMBOL && third->data.symbol.name == SYM_All)
            keep_ends = 1;
        else if (third->type == EXPR_INTEGER)
            n_limit = (long)third->data.integer;
        else
            return NULL;                        /* unsupported: leave unevaluated */
    }

    /* Delimiter pattern (default Whitespace). */
    int own_patt = 0;
    Expr* patt;
    if (pargc >= 2) {
        patt = a[1];
    } else {
        patt = expr_new_symbol(SYM_Whitespace);
        own_patt = 1;
    }

    int patt_is_null = (patt->type == EXPR_STRING && patt->data.string[0] == '\0');

    RegexRule* rules = NULL;
    int nr = 0;
    if (!patt_is_null) {
        nr = regex_rules_build_ex(patt, /*anchored=*/0, caseless, &rules, "StringSplit");
        if (nr < 0) {
            if (own_patt) expr_free(patt);
            return NULL;
        }
    }

    Expr* result;
    LLReader sv;
    if (ll_rd_open(subject, &sv)) {   /* a List or a LinkedList of subjects */
        size_t m = sv.n;
        Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
        for (size_t i = 0; i < m; i++) {
            Expr* si = ll_rd_at(&sv, i);
            if (si->type != EXPR_STRING)
                out[i] = expr_copy(si);
            else
                out[i] = patt_is_null ? ss_null(si->data.string, n_limit)
                                      : ss_scalar(si->data.string, rules, nr,
                                                  n_limit, keep_ends);
        }
        result = ll_seq_like(subject, out, m);
        free(out);
        ll_rd_close(&sv);
    } else if (subject->type == EXPR_STRING) {
        result = patt_is_null ? ss_null(subject->data.string, n_limit)
                              : ss_scalar(subject->data.string, rules, nr,
                                          n_limit, keep_ends);
    } else {
        result = NULL;                          /* non-string subject: unevaluated */
    }

    if (rules) regex_rules_free(rules, nr);
    if (own_patt) expr_free(patt);
    return result;
}
