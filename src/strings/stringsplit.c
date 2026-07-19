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
#include "sym_names.h"
#include "common.h"
#include "eval.h"

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
    size_t n = e - s;
    char* buf = malloc(n + 1);
    if (!buf) return expr_new_string("");
    memcpy(buf, subj + s, n);
    buf[n] = '\0';
    Expr* r = expr_new_string(buf);
    free(buf);
    return r;
}

static int is_empty_string(const Expr* e) {
    return e->type == EXPR_STRING && e->data.string[0] == '\0';
}

/*
 * Build the value inserted at a delimiter for a matched rule, or NULL to insert
 * nothing. String RHS uses $0/$n template expansion. For a non-string RHS whose
 * delimiter LHS is a named Pattern[x, ...], bind x to the matched text and
 * evaluate the RHS (handles `x:"--" :> x`); otherwise evaluate the RHS as a
 * constant.
 */
static Expr* make_insertion(const RegexRule* rule, const char* subj,
                            const size_t* ov, size_t pairs,
                            size_t ms, size_t me) {
    if (!rule->rhs) return NULL;

    if (rule->rhs->type == EXPR_STRING) {
        char* rep = regex_rule_replacement(rule, subj, ov, pairs);
        if (!rep) return NULL;
        Expr* r = expr_new_string(rep);
        free(rep);
        return r;
    }

    Expr* rhs_copy = expr_copy(rule->rhs);
    Expr* lhs = rule->lhs;
    if (lhs && lhs->type == EXPR_FUNCTION &&
        lhs->data.function.head->type == EXPR_SYMBOL &&
        lhs->data.function.head->data.symbol.name == SYM_Pattern &&
        lhs->data.function.arg_count == 2 &&
        lhs->data.function.args[0]->type == EXPR_SYMBOL) {
        Expr* var = expr_copy(lhs->data.function.args[0]);
        Expr* val = substr_expr(subj, ms, me);
        Expr* rule_args[2] = { var, val };
        Expr* r = expr_new_function(expr_new_symbol(SYM_Rule), rule_args, 2);
        Expr* ra_args[2] = { rhs_copy, r };
        Expr* ra = expr_new_function(expr_new_symbol(SYM_ReplaceAll), ra_args, 2);
        return eval_and_free(ra);
    }
    return eval_and_free(rhs_copy);
}

/* Split one string subject. */
static Expr* ss_scalar(const char* subj, RegexRule* rules, int nr,
                       long n_limit, int keep_ends) {
    size_t len = strlen(subj);
    PieceVec pv = {0};

    size_t last_end = 0, pos = 0;
    long ns = 0;                        /* substrings emitted so far */

    while (pos <= len) {
        if (n_limit > 0 && ns >= n_limit - 1) break;   /* remainder is last piece */

        int found = 0, best_rule = 0;
        size_t best_ms = 0, best_me = 0, best_pairs = 0;
        size_t best_ov[REGEX_MAX_PAIRS * 2];

        for (int i = 0; i < nr; i++) {
            int gc = regex_group_count(rules[i].prog) + 1;
            size_t pairs = (gc > REGEX_MAX_PAIRS) ? REGEX_MAX_PAIRS : (size_t)gc;
            size_t ov[REGEX_MAX_PAIRS * 2];
            if (regex_match(rules[i].prog, subj, len, pos, ov, pairs) == 1) {
                if (!found || ov[0] < best_ms) {
                    found = 1;
                    best_rule = i;
                    best_ms = ov[0];
                    best_me = ov[1];
                    best_pairs = pairs;
                    memcpy(best_ov, ov, sizeof(size_t) * 2 * pairs);
                }
            }
        }
        if (!found) break;

        pv_push(&pv, substr_expr(subj, last_end, best_ms), 1);   /* interior empties kept */
        ns++;

        Expr* ins = make_insertion(&rules[best_rule], subj, best_ov, best_pairs,
                                   best_ms, best_me);
        if (ins) pv_push(&pv, ins, 0);

        last_end = best_me;
        pos = (best_me > best_ms) ? best_me : best_ms + 1;       /* zero-width guard */
    }

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

/* Null delimiter "": split into individual (byte) characters. */
static Expr* ss_null(const char* subj) {
    size_t len = strlen(subj);
    Expr** items = malloc(sizeof(Expr*) * (len ? len : 1));
    for (size_t i = 0; i < len; i++) items[i] = substr_expr(subj, i, i + 1);
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), items, len);
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
    int caseless = 0;
    size_t pargc = argc;
    while (pargc >= 2) {
        int v;
        if (is_ignorecase_opt(a[pargc - 1], &v)) { caseless = v; pargc--; }
        else break;
    }
    if (pargc < 1 || pargc > 3) return builtin_arg_error("StringSplit", argc, 1, 3);

    Expr* subject = a[0];

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
    if (subject->type == EXPR_FUNCTION &&
        subject->data.function.head->type == EXPR_SYMBOL &&
        subject->data.function.head->data.symbol.name == SYM_List) {
        size_t m = subject->data.function.arg_count;
        Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
        for (size_t i = 0; i < m; i++) {
            Expr* si = subject->data.function.args[i];
            if (si->type != EXPR_STRING)
                out[i] = expr_copy(si);
            else
                out[i] = patt_is_null ? ss_null(si->data.string)
                                      : ss_scalar(si->data.string, rules, nr,
                                                  n_limit, keep_ends);
        }
        result = expr_new_function(expr_new_symbol(SYM_List), out, m);
        free(out);
    } else if (subject->type == EXPR_STRING) {
        result = patt_is_null ? ss_null(subject->data.string)
                              : ss_scalar(subject->data.string, rules, nr,
                                          n_limit, keep_ends);
    } else {
        result = NULL;                          /* non-string subject: unevaluated */
    }

    if (rules) regex_rules_free(rules, nr);
    if (own_patt) expr_free(patt);
    return result;
}
