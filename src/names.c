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
 *
 * Context handling: symbols are stored under bare names for the System` and
 * Global` contexts (builtins and unqualified user symbols) and under explicit
 * backtick-qualified names for other contexts.  A pattern element that itself
 * contains a backtick (e.g. "System`*") is matched against -- and returns --
 * each symbol's fully context-qualified name (System`Sin, Global`x, ...); a
 * plain pattern (no backtick) is matched against, and returns, the stored name
 * exactly as before.  This is what makes Names["System`*"] enumerate the
 * builtins instead of returning {} (nothing is stored with a literal "System`"
 * prefix).
 */

#include "names.h"
#include "symtab.h"
#include "sym_names.h"
#include "sym_intern.h"
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
/* Context-qualified name helpers                                     */
/* ------------------------------------------------------------------ */

/* C99-safe strdup replacement (strdup is POSIX, not C99). */
static char* name_strdup(const char* s) {
    size_t n = strlen(s) + 1;
    char* r = (char*)malloc(n);
    if (r) memcpy(r, s, n);
    return r;
}

static int name_has_backtick(const char* s) {
    for (; *s; s++) if (*s == '`') return 1;
    return 0;
}

/* Fully context-qualified name of a symbol (caller frees). Names that already
 * carry an explicit context prefix are returned verbatim; bare names are
 * qualified with "System`" when the symbol is a builtin (or a kernel-interned
 * System symbol) and "Global`" otherwise, mirroring Context[]'s home-context
 * rule. */
static char* full_qualified_name(const char* nm, SymbolDef* def) {
    if (name_has_backtick(nm)) return name_strdup(nm);
    const char* ctx =
        (def && (def->builtin_func != NULL || intern_is_system(nm)))
            ? "System`" : "Global`";
    size_t nc = strlen(ctx), nn = strlen(nm);
    char* r = (char*)malloc(nc + nn + 1);
    if (!r) return NULL;
    memcpy(r, ctx, nc);
    memcpy(r + nc, nm, nn + 1);
    return r;
}

/* ------------------------------------------------------------------ */
/* Collector: gather matching names during symtab_for_each            */
/* ------------------------------------------------------------------ */

typedef struct {
    const NamePat* pats;
    int            npats;
    int            match_all;   /* Names[] with no argument: take every name */
    char**         names;       /* owned emit strings (bare or context-qualified) */
    size_t         count, cap;
} NameCollect;

/* Determine the string to return for `name` under a single pattern, or NULL if
 * it does not match.  A glob containing a backtick matches (and yields) the
 * fully-qualified name; a plain glob matches (and yields) the stored name.  A
 * regex is tried against the stored name first, then the qualified name, so a
 * context-bearing regex works too.  `*full` caches the qualified name across
 * the caller's pattern loop; the caller frees it. */
static char* match_one(const char* name, SymbolDef* def,
                       const NamePat* p, char** full) {
    if (!p->is_regex && p->glob && name_has_backtick(p->glob)) {
        if (!*full) *full = full_qualified_name(name, def);
        return (*full && pat_matches(p, *full)) ? name_strdup(*full) : NULL;
    }
    if (!p->is_regex) {
        return pat_matches(p, name) ? name_strdup(name) : NULL;
    }
    /* regex */
    if (pat_matches(p, name)) return name_strdup(name);
    if (!*full) *full = full_qualified_name(name, def);
    return (*full && pat_matches(p, *full)) ? name_strdup(*full) : NULL;
}

static void collect_visit(const char* name, SymbolDef* def, void* user) {
    NameCollect* c = (NameCollect*)user;
    char* emit = NULL;

    if (c->match_all) {
        emit = name_strdup(name);
    } else {
        char* full = NULL;   /* lazily built, reused across patterns */
        for (int i = 0; !emit && i < c->npats; i++)
            emit = match_one(name, def, &c->pats[i], &full);
        free(full);
    }
    if (!emit) return;

    if (c->count == c->cap) {
        size_t nc = c->cap ? c->cap * 2 : 64;
        char** np = realloc(c->names, nc * sizeof(char*));
        if (!np) { free(emit); return; }  /* OOM: drop this name rather than crash */
        c->names = np;
        c->cap = nc;
    }
    c->names[c->count++] = emit;
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
        if (!items) {
            for (size_t i = 0; i < c.count; i++) free(c.names[i]);
            free(c.names);
            return NULL;
        }
        for (size_t i = 0; i < c.count; i++)
            items[i] = expr_new_string(c.names[i]);
        /* Sort by canonical order so Names[p] === Sort[Names[p]]. */
        if (c.count > 1)
            qsort(items, c.count, sizeof(Expr*), name_expr_cmp);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, c.count);
    free(items);
    for (size_t i = 0; i < c.count; i++) free(c.names[i]);
    free(c.names);
    return out;
}

/* ------------------------------------------------------------------ */

void names_init(void) {
    symtab_add_builtin("Names", builtin_names);
    symtab_get_def("Names")->attributes |= ATTR_PROTECTED;
}
