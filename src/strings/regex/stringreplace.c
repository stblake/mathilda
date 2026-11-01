/*
 * stringreplace.c - StringReplace[subject, rule | {rules...}]
 *
 * Replaces every non-overlapping match of a rule's pattern by its right-hand
 * side, scanning left to right; unmatched text is copied verbatim.  Pattern
 * names in the LHS are bound to the text they matched and substituted into the
 * RHS (regex_rule_eval_rhs); a non-string value turns the result into a
 * StringExpression.  For a RegularExpression pattern a string RHS may use
 * $0/$1... for the whole match and capture groups.
 * With several rules, at each position the leftmost match wins, ties broken by
 * rule order.  A list of subjects threads.
 */

#include "picostrings.h"
#include "regex_common.h"
#include "sym_names.h"
#include "symtab.h"
#include "message.h"
#include "print.h"

#include <string.h>

/*
 * Result assembly. Unmatched text and string replacements accumulate in `buf`;
 * a replacement that is not a string (x:DigitCharacter :> ToExpression[x])
 * cannot be spliced into text, so -- as in Mathematica -- the result becomes a
 * StringExpression of the string runs and those values:
 *   StringReplace["abc123", x:DigitCharacter :> f[x]]
 *     = "abc" ~~ f["1"] ~~ f["2"] ~~ f["3"]
 * Sequence[...] replacements splice (Sequence[] deletes).
 */
typedef struct {
    RegexBuf buf;        /* current run of text */
    Expr**   items;      /* StringExpression pieces, once any non-string appears */
    size_t   n, cap;
    int      oom;
} SrOut;

static void sr_push_item(SrOut* o, Expr* e) {
    if (o->n == o->cap) {
        size_t nc = o->cap ? o->cap * 2 : 8;
        Expr** ni = realloc(o->items, nc * sizeof(Expr*));
        if (!ni) { expr_free(e); o->oom = 1; return; }
        o->items = ni;
        o->cap = nc;
    }
    o->items[o->n++] = e;
}

/* Close the current text run as a string piece (only once pieces exist). */
static void sr_flush(SrOut* o) {
    if (o->buf.len == 0) return;
    sr_push_item(o, expr_new_string(o->buf.p));
    o->buf.len = 0;
    o->buf.p[0] = '\0';
}

static void sr_text(SrOut* o, const char* p, size_t n) {
    if (n && regexbuf_add(&o->buf, p, n)) o->oom = 1;
}

/* Append one evaluated replacement value; takes ownership of v. */
static void sr_value(SrOut* o, Expr* v) {
    if (v->type == EXPR_STRING) {
        sr_text(o, v->data.string, strlen(v->data.string));
        expr_free(v);
        return;
    }
    if (v->type == EXPR_FUNCTION && v->data.function.head->type == EXPR_SYMBOL &&
        v->data.function.head->data.symbol.name == SYM_Sequence) {
        for (size_t i = 0; i < v->data.function.arg_count; i++)
            sr_value(o, expr_copy(v->data.function.args[i]));
        expr_free(v);
        return;
    }
    sr_flush(o);
    sr_push_item(o, v);
}

/* Replace in one subject. Returns a new String, or StringExpression[...] when a
 * replacement value is not a string; NULL on OOM. */
static Expr* sr_finish(SrOut* o);

/* Replace in one subject, at most `limit` matches (0: all). */
static Expr* sr_scalar(const char* subj, RegexRule* rules, int nr, size_t limit) {
    size_t len = strlen(subj);
    SrOut o;
    memset(&o, 0, sizeof o);

    RegexCursor cur;
    if (regex_cursor_init(&cur, rules, nr, subj, len) != 0) return NULL;

    size_t pos = 0, done = 0;
    while (pos <= len && !o.oom) {
        if (limit && done == limit) {
            sr_text(&o, subj + pos, len - pos);          /* rest untouched */
            break;
        }
        int rule;
        const size_t* ov;
        size_t np;
        int m = regex_cursor_next(&cur, pos, &rule, &ov, &np);
        if (m < 0) { o.oom = 1; break; }
        if (m == 0) {
            sr_text(&o, subj + pos, len - pos);          /* tail */
            break;
        }
        size_t ms = ov[0], me = ov[1];
        sr_text(&o, subj + pos, ms - pos);               /* literal before match */
        sr_value(&o, regex_rule_eval_rhs(&rules[rule], subj, ov, np));
        done++;

        if (me > ms) {
            pos = me;
        } else {                                         /* zero-width match */
            size_t next = regex_next_char(subj, len, ms);
            if (ms < len) sr_text(&o, subj + ms, next - ms);
            pos = next;
        }
    }
    regex_cursor_free(&cur);

    return sr_finish(&o);
}

/* The assembled result (a String, or a StringExpression when a value was not
 * a string); releases o. NULL on OOM. */
static Expr* sr_finish(SrOut* o) {
    Expr* result = NULL;
    if (!o->oom) {
        if (o->n == 0) {
            result = expr_new_string(o->buf.p ? o->buf.p : "");
        } else {
            sr_flush(o);
            if (!o->oom) {
                result = expr_new_function(expr_new_symbol(SYM_StringExpression),
                                           o->items, o->n);
                o->n = 0;
            }
        }
    }
    for (size_t i = 0; i < o->n; i++) expr_free(o->items[i]);
    free(o->items);
    free(o->buf.p);
    return result;
}

/* Is e a List[...] of Strings (or a String)? */
static int strings_or_list(const Expr* e) {
    if (e->type == EXPR_STRING) return 1;
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_List)
        return 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (e->data.function.args[i]->type != EXPR_STRING) return 0;
    return 1;
}

static Expr* strse_message(const char* head, Expr* res) {
    char* shown = expr_to_string_outputform(res);
    mth_message(head, "strse", "String or list of strings expected at position 1 in %s.",
                shown ? shown : head);
    free(shown);
    return NULL;
}

/*
 * Shared argument handling of StringReplace / StringDelete / StringReplaceList:
 * strips trailing IgnoreCase options (default from Options[head]), and decodes
 * the occurrence limit n (an integer >= 0, or Infinity). Returns the number of
 * positional arguments, or -1 if the limit is malformed.
 */
/* MetaCharacters setting: None (the default: string patterns are literal),
 * Automatic (a literal string is read as a regular expression), or a custom
 * metacharacter set (a 3-character string or a list of three strings), which
 * is not supported and leaves the call unevaluated. */
typedef enum { SR_META_NONE, SR_META_AUTO, SR_META_CUSTOM } SrMeta;

/* If `e` is MetaCharacters -> v (or :>), decode v into *meta and return 1;
 * an invalid v gives head::meta and sets *meta = -1.  Else return 0. */
static int sr_meta_opt(const Expr* e, const char* head, int* meta) {
    if (e->type != EXPR_FUNCTION || e->data.function.arg_count != 2 ||
        e->data.function.head->type != EXPR_SYMBOL ||
        (e->data.function.head->data.symbol.name != SYM_Rule &&
         e->data.function.head->data.symbol.name != SYM_RuleDelayed))
        return 0;
    const Expr* k = e->data.function.args[0];
    if (k->type != EXPR_SYMBOL || k->data.symbol.name != SYM_MetaCharacters) return 0;
    const Expr* v = e->data.function.args[1];
    if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_None)      *meta = SR_META_NONE;
    else if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_Automatic) *meta = SR_META_AUTO;
    else if (v->type == EXPR_STRING && strlen(v->data.string) == 3)     *meta = SR_META_CUSTOM;
    else if (v->type == EXPR_FUNCTION && v->data.function.head->type == EXPR_SYMBOL &&
             v->data.function.head->data.symbol.name == SYM_List && v->data.function.arg_count == 3 &&
             v->data.function.args[0]->type == EXPR_STRING &&
             v->data.function.args[1]->type == EXPR_STRING &&
             v->data.function.args[2]->type == EXPR_STRING)             *meta = SR_META_CUSTOM;
    else {
        char* shown = expr_to_string_outputform((Expr*)v);
        mth_message(head, "meta", "Value of option MetaCharacters -> %s is not None, Automatic, "
                    "a string of length 3, or a list of three strings.", shown ? shown : "");
        free(shown);
        *meta = -1;
    }
    return 1;
}

static int sr_args(Expr* res, const char* head, int* caseless, size_t* limit, int* has_limit,
                   int* meta) {
    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;
    *caseless = 0;
    *meta = SR_META_NONE;
    Expr* defs = symtab_get_options(head);                 /* borrowed */
    if (defs && defs->type == EXPR_FUNCTION)
        for (size_t i = 0; i < defs->data.function.arg_count; i++) {
            regex_match_opt(defs->data.function.args[i], SYM_IgnoreCase, caseless, 0);
            sr_meta_opt(defs->data.function.args[i], head, meta);
        }
    while (argc > 1 && (regex_match_opt(a[argc - 1], SYM_IgnoreCase, caseless, 0) ||
                        sr_meta_opt(a[argc - 1], head, meta))) {
        if (*meta < 0) return -1;
        argc--;
    }
    *limit = 0;
    *has_limit = 0;
    if (argc == 3) {
        Expr* n = a[2];
        *has_limit = 1;
        if (n->type == EXPR_SYMBOL && n->data.symbol.name == SYM_Infinity) *limit = 0;
        else if (n->type == EXPR_INTEGER && n->data.integer >= 0) *limit = (size_t)n->data.integer;
        else return -1;
    }
    return (int)argc;
}

/* MetaCharacters -> Automatic: a copy of the pattern argument in which every
 * literal string pattern (bare, or the LHS of a rule) is RegularExpression[s].
 * Returns a new Expr. */
static Expr* meta_auto_rules(const Expr* patt) {
    if (patt->type == EXPR_STRING) {
        Expr* a = expr_new_string(patt->data.string);
        return expr_new_function(expr_new_symbol(SYM_RegularExpression), &a, 1);
    }
    if (patt->type != EXPR_FUNCTION || patt->data.function.head->type != EXPR_SYMBOL)
        return expr_copy((Expr*)patt);
    const char* h = patt->data.function.head->data.symbol.name;
    if (h == SYM_List) {
        size_t n = patt->data.function.arg_count;
        Expr* out = expr_new_function(expr_new_symbol(SYM_List), NULL, n);
        for (size_t i = 0; i < n; i++) out->data.function.args[i] = meta_auto_rules(patt->data.function.args[i]);
        return out;
    }
    if ((h == SYM_Rule || h == SYM_RuleDelayed) && patt->data.function.arg_count == 2) {
        Expr* ra[2] = { meta_auto_rules(patt->data.function.args[0]),
                        expr_copy(patt->data.function.args[1]) };
        return expr_new_function(expr_new_symbol(h), ra, 2);
    }
    return expr_copy((Expr*)patt);
}

/* Is e the empty list {}? */
static int is_empty_list(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List && e->data.function.arg_count == 0;
}

/* rules for StringDelete: each pattern p becomes p -> "". Returns a new Expr. */
static Expr* delete_rules(Expr* patt) {
    int is_list = patt->type == EXPR_FUNCTION && patt->data.function.head->type == EXPR_SYMBOL &&
                  patt->data.function.head->data.symbol.name == SYM_List;
    size_t n = is_list ? patt->data.function.arg_count : 1;
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), NULL, n);
    for (size_t i = 0; i < n; i++) {
        Expr* p = is_list ? patt->data.function.args[i] : patt;
        Expr* rargs[2] = { expr_copy(p), expr_new_string("") };
        out->data.function.args[i] = expr_new_function(expr_new_symbol(SYM_Rule), rargs, 2);
    }
    return out;
}

typedef enum { SR_REPLACE, SR_DELETE, SR_LIST } SrKind;

static Expr* sr_list_scalar(const char* subj, RegexRule* rules, int nr, size_t limit,
                            int has_limit);

static Expr* sr_dispatch(Expr* res, const char* head, SrKind kind) {
    if (res->type != EXPR_FUNCTION) return NULL;
    int caseless, meta;
    size_t limit;
    int has_limit;
    int argc = sr_args(res, head, &caseless, &limit, &has_limit, &meta);
    /* One argument is the operator form head[rules][s] (see src/opform.c); a
     * malformed limit (-1) also stays unevaluated. */
    if (argc < 2 || argc > 3) return NULL;
    if (argc == 3 && kind == SR_DELETE) return NULL;
    Expr* arg0 = res->data.function.args[0];
    Expr* patt = res->data.function.args[1];

    if (!strings_or_list(arg0)) return strse_message(head, res);
    if (meta == SR_META_CUSTOM) return NULL;

    /* No rules: nothing is replaced or deleted. */
    if (is_empty_list(patt) && kind != SR_LIST) return expr_copy(arg0);

    Expr* meta_owned = (meta == SR_META_AUTO) ? meta_auto_rules(patt) : NULL;
    if (meta_owned) patt = meta_owned;
    Expr* owned = (kind == SR_DELETE) ? delete_rules(patt) : NULL;
    RegexRule* rules;
    int nr = regex_rules_build_opts(owned ? owned : patt, 0,
                                    caseless ? REGEX_CASELESS : 0u, &rules, head);
    if (nr < 0) {
        if (owned) expr_free(owned);
        if (meta_owned) expr_free(meta_owned);
        return NULL;
    }

    /* A replacement needs rule(s): every element must carry an RHS. */
    for (int i = 0; i < nr; i++) {
        if (!rules[i].rhs) {
            regex_rules_free(rules, nr);
            if (owned) expr_free(owned);
            if (meta_owned) expr_free(meta_owned);
            return NULL;
        }
    }

    Expr* result;
    if (arg0->type == EXPR_STRING) {
        result = (kind == SR_LIST)
            ? sr_list_scalar(arg0->data.string, rules, nr, limit, has_limit)
            : sr_scalar(arg0->data.string, rules, nr, limit);
    } else {
        size_t n = arg0->data.function.arg_count;
        result = expr_new_function(expr_new_symbol(SYM_List), NULL, n);
        for (size_t i = 0; result && i < n; i++) {
            const char* si = arg0->data.function.args[i]->data.string;
            result->data.function.args[i] = (kind == SR_LIST)
                ? sr_list_scalar(si, rules, nr, limit, has_limit)
                : sr_scalar(si, rules, nr, limit);
            if (!result->data.function.args[i]) { expr_free(result); result = NULL; }
        }
    }
    regex_rules_free(rules, nr);
    if (owned) expr_free(owned);          /* after the rules: they borrow into it */
    if (meta_owned) expr_free(meta_owned);
    return result;
}

/*
 * StringReplaceList: one result per single replacement -- every match start
 * of every rule (overlapping), in order of position then rule, each giving the
 * subject with just that match replaced. At most `limit` results when
 * has_limit (n = 0 gives {}).
 */
static Expr* sr_list_scalar(const char* subj, RegexRule* rules, int nr, size_t limit,
                            int has_limit) {
    size_t len = strlen(subj);
    RegexScan scan;
    if (has_limit && limit == 0) return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    if (regex_scan_limit(subj, len, rules, nr, REGEX_OV_TRUE, /*want_captures=*/1,
                         limit, &scan) < 0)
        return NULL;
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), NULL, scan.count);
    for (size_t k = 0; out && k < scan.count; k++) {
        const RegexSpan* sp = &scan.spans[k];
        SrOut o;
        memset(&o, 0, sizeof o);
        sr_text(&o, subj, sp->ms);
        sr_value(&o, regex_rule_eval_rhs(&rules[sp->rule], subj,
                                         scan.caps + sp->caps_off, sp->npairs));
        sr_text(&o, subj + sp->me, len - sp->me);
        out->data.function.args[k] = sr_finish(&o);
        if (!out->data.function.args[k]) { expr_free(out); out = NULL; }
    }
    regex_scan_free(&scan);
    return out;
}

Expr* builtin_stringreplace(Expr* res)     { return sr_dispatch(res, "StringReplace", SR_REPLACE); }
Expr* builtin_stringdelete(Expr* res)      { return sr_dispatch(res, "StringDelete", SR_DELETE); }
Expr* builtin_stringreplacelist(Expr* res) { return sr_dispatch(res, "StringReplaceList", SR_LIST); }
