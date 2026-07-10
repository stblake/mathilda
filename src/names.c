/*
 * names.c - Names[] and friends: enumerate symbol-table names by pattern.
 *
 *   Names["string"]            names matching a string pattern
 *   Names[patt]                names matching an arbitrary string pattern patt
 *   Names[{p1, p2, ...}]       names matching any of the p_i
 *   Names[]                    all names in the symbol table
 *
 * A string pattern is matched against the whole name (anchored) and supports
 * two metacharacters:
 *   *   matches zero or more characters
 *   @   matches one or more characters that are NOT uppercase letters
 * Every other character (including the ` used in context prefixes) is literal.
 *
 * A pattern element may instead be RegularExpression["re"], matched against the
 * whole name via the PCRE2 engine (src/strings/regex).  When the engine is
 * unavailable the call emits Names::regavail and stays unevaluated.
 *
 * The result is a List of Strings sorted ascending by byte value (strcmp),
 * matching Wolfram-Language ordering for the common ASCII case.  All symbols in
 * the table are candidates -- there is no filtering of internal helper symbols.
 */

#include "names.h"
#include "symtab.h"
#include "sym_names.h"
#include "attr.h"
#include "regex_engine.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* String-pattern (glob) matcher for the * and @ metacharacters.      */
/* ------------------------------------------------------------------ */

static int name_is_nonupper(char c) {
    return c != '\0' && !(c >= 'A' && c <= 'Z');
}

/* Recursive backtracking whole-string match of `p` against `s`. */
static int wl_glob_match(const char* p, const char* s) {
    if (*p == '\0') return *s == '\0';
    if (*p == '*') {
        /* zero or more of any character */
        if (wl_glob_match(p + 1, s)) return 1;          /* consume nothing */
        if (*s == '\0') return 0;
        return wl_glob_match(p, s + 1);                 /* consume one char */
    }
    if (*p == '@') {
        /* one or more non-uppercase characters */
        if (!name_is_nonupper(*s)) return 0;            /* need at least one */
        if (wl_glob_match(p + 1, s + 1)) return 1;      /* run ends here */
        return wl_glob_match(p, s + 1);                 /* run continues */
    }
    if (*p == *s) return wl_glob_match(p + 1, s + 1);   /* literal */
    return 0;
}

/* ------------------------------------------------------------------ */
/* Compiled pattern set                                               */
/* ------------------------------------------------------------------ */

typedef struct {
    int          is_regex;   /* 1 -> use prog; 0 -> use glob */
    const char*  glob;       /* borrowed from the input Expr (when !is_regex) */
    RegexProgram* prog;      /* owned (when is_regex); freed by free_pats */
} NamePat;

/* RegularExpression["re"] -> anchored PCRE program, or NULL if `e` is not one. */
static int is_regular_expression(Expr* e, const char** re_out) {
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_RegularExpression &&
        e->data.function.arg_count == 1 &&
        e->data.function.args[0]->type == EXPR_STRING) {
        *re_out = e->data.function.args[0]->data.string;
        return 1;
    }
    return 0;
}

/* Compile a single element into `out`. Returns 0 on success, -1 on an
 * unsupported element / regex compile failure / missing regex engine. */
static int build_pat(Expr* e, NamePat* out) {
    const char* re;
    if (e->type == EXPR_STRING) {
        out->is_regex = 0;
        out->glob = e->data.string;
        out->prog = NULL;
        return 0;
    }
    if (is_regular_expression(e, &re)) {
        if (!regex_available()) {
            fprintf(stderr, "Names::regavail: regular-expression support is not "
                            "available; rebuild Mathilda with PCRE2 "
                            "(USE_REGEX=1).\n");
            return -1;
        }
        /* Anchor for a whole-string match: \A(?:re)\z */
        size_t n = strlen(re);
        char* wrapped = malloc(n + 9);
        if (!wrapped) return -1;
        memcpy(wrapped, "\\A(?:", 5);
        memcpy(wrapped + 5, re, n);
        memcpy(wrapped + 5 + n, ")\\z", 3);
        wrapped[8 + n] = '\0';

        char err[256];
        RegexProgram* prog = regex_compile(wrapped, err, sizeof err);
        free(wrapped);
        if (!prog) {
            fprintf(stderr, "Names::regex: %s\n", err);
            return -1;
        }
        out->is_regex = 1;
        out->glob = NULL;
        out->prog = prog;
        return 0;
    }
    return -1;   /* unsupported pattern element */
}

static void free_pats(NamePat* pats, int n) {
    if (!pats) return;
    for (int i = 0; i < n; i++)
        if (pats[i].is_regex) regex_free(pats[i].prog);
    free(pats);
}

static int pat_matches(const NamePat* p, const char* name) {
    if (p->is_regex) {
        size_t ov[2];
        return regex_match(p->prog, name, strlen(name), 0, ov, 1) == 1;
    }
    return wl_glob_match(p->glob, name);
}

/* ------------------------------------------------------------------ */
/* Collector: gather matching names during symtab_for_each            */
/* ------------------------------------------------------------------ */

typedef struct {
    const NamePat* pats;
    int            npats;
    int            match_all;   /* Names[] with no argument: take every name */
    const char**   names;       /* borrowed interned names (not owned) */
    size_t         count, cap;
} NameCollect;

static void collect_visit(const char* name, SymbolDef* def, void* user) {
    (void)def;
    NameCollect* c = (NameCollect*)user;
    int keep = c->match_all;
    for (int i = 0; !keep && i < c->npats; i++)
        if (pat_matches(&c->pats[i], name)) keep = 1;
    if (!keep) return;

    if (c->count == c->cap) {
        size_t nc = c->cap ? c->cap * 2 : 64;
        const char** np = realloc(c->names, nc * sizeof(const char*));
        if (!np) return;   /* OOM: drop this name rather than crash */
        c->names = np;
        c->cap = nc;
    }
    c->names[c->count++] = name;
}

/* Order two String Exprs by Mathilda's canonical comparison, so that the
 * result of Names[patt] is identical to Sort[Names[patt]]. */
static int name_expr_cmp(const void* a, const void* b) {
    return expr_compare(*(Expr* const*)a, *(Expr* const*)b);
}

/* ------------------------------------------------------------------ */
/* Builtin                                                            */
/* ------------------------------------------------------------------ */

Expr* builtin_names(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc > 1) return NULL;

    NamePat* pats = NULL;
    int npats = 0;
    int match_all = 0;

    if (argc == 0) {
        match_all = 1;
    } else {
        Expr* arg = res->data.function.args[0];
        /* A List argument is a set of alternative patterns. */
        if (arg->type == EXPR_FUNCTION &&
            arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol == SYM_List) {
            size_t n = arg->data.function.arg_count;
            if (n > 0) {
                pats = calloc(n, sizeof(NamePat));
                if (!pats) return NULL;
                for (size_t i = 0; i < n; i++) {
                    if (build_pat(arg->data.function.args[i], &pats[npats]) != 0) {
                        free_pats(pats, npats);
                        return NULL;
                    }
                    npats++;
                }
            }
            /* n == 0 (empty list) -> zero patterns -> matches nothing -> {} */
        } else {
            pats = calloc(1, sizeof(NamePat));
            if (!pats) return NULL;
            if (build_pat(arg, &pats[0]) != 0) {
                free_pats(pats, 0);
                return NULL;
            }
            npats = 1;
        }
    }

    NameCollect c = { pats, npats, match_all, NULL, 0, 0 };
    symtab_for_each(collect_visit, &c);
    free_pats(pats, npats);

    Expr** items = NULL;
    if (c.count > 0) {
        items = malloc(c.count * sizeof(Expr*));
        if (!items) { free(c.names); return NULL; }
        for (size_t i = 0; i < c.count; i++)
            items[i] = expr_new_string(c.names[i]);
        /* Sort by canonical order so Names[p] === Sort[Names[p]]. */
        if (c.count > 1)
            qsort(items, c.count, sizeof(Expr*), name_expr_cmp);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, c.count);
    free(items);
    free(c.names);
    return out;
}

/* ------------------------------------------------------------------ */

void names_init(void) {
    symtab_add_builtin("Names", builtin_names);
    symtab_get_def("Names")->attributes |= ATTR_PROTECTED;
}
