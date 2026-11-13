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
 * A string pattern containing a backtick is split at its LAST backtick into a
 * context part and a short-name part, matched separately against the symbol's
 * context and short name: in the context part `*` may span subcontexts
 * ("Tes*`f*" and "*`fun" both find Test`Private`fun), but the short-name part
 * never crosses a backtick, so "Test`*" lists only the symbols of Test` itself,
 * not those of Test`Private` (issue #109).
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
 * backtick-qualified names for other contexts, so each symbol's fully
 * context-qualified name (System`Sin, Global`x, Foo`bar) is reconstructed
 * before matching.  Each name is RETURNED the way it prints -- without its
 * context when that context is $Context or on $ContextPath, fully qualified
 * otherwise (context_display_name) -- and the pattern decides what it is
 * matched against:
 *   - a pattern containing a backtick ("System`*", "*`x", a regex such as
 *     "Foo`ba.") is matched against the fully-qualified name;
 *   - a pattern without one ("Sin", "F*") names symbols as they are typed, so
 *     it is matched against the short name of the symbols whose context is
 *     visible, and never sees a symbol in a context off $ContextPath.
 * So Names["System`File*"] gives {"File", ...}; Names["Foo`*"] gives
 * {"Foo`bar", ...} while Foo` is off the path and {"bar", ...} inside
 * Begin["Foo`"]; Names["bar"] finds Foo`bar only once Foo` is visible
 * (issue #105).  Names[] lists every symbol, by display name.
 */

#include "names.h"
#include "symtab.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "attr.h"
#include "regex_engine.h"
#include "context.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "message.h"   /* mth_message: Quiet/Check funnel */

/* ------------------------------------------------------------------ */
/* String-pattern (glob) matcher for the * and @ metacharacters.      */
/* ------------------------------------------------------------------ */

static int name_is_nonupper(char c) {
    return c != '\0' && !(c >= 'A' && c <= 'Z');
}

/* Recursive backtracking whole-string match of the pattern [p, pe) against
 * the subject [s, se).  Bounded by end pointers so a context-qualified name
 * can be matched piecewise (context part, then short-name part) in place. */
static int wl_glob_match_n(const char* p, const char* pe,
                           const char* s, const char* se) {
    if (p == pe) return s == se;
    if (*p == '*') {
        /* zero or more of any character */
        if (wl_glob_match_n(p + 1, pe, s, se)) return 1;  /* consume nothing */
        if (s == se) return 0;
        return wl_glob_match_n(p, pe, s + 1, se);         /* consume one char */
    }
    if (*p == '@') {
        /* one or more non-uppercase characters */
        if (s == se || !name_is_nonupper(*s)) return 0;   /* need at least one */
        if (wl_glob_match_n(p + 1, pe, s + 1, se)) return 1;  /* run ends here */
        return wl_glob_match_n(p, pe, s + 1, se);             /* run continues */
    }
    if (s != se && *p == *s) return wl_glob_match_n(p + 1, pe, s + 1, se);
    return 0;
}

static int wl_glob_match(const char* p, const char* s) {
    return wl_glob_match_n(p, p + strlen(p), s, s + strlen(s));
}

/* ------------------------------------------------------------------ */
/* Compiled pattern set                                               */
/* ------------------------------------------------------------------ */

typedef struct {
    int          is_regex;   /* 1 -> use prog; 0 -> use glob */
    int          has_ctx;    /* pattern text contains a backtick */
    const char*  glob;       /* borrowed from the input Expr (when !is_regex) */
    RegexProgram* prog;      /* owned (when is_regex); freed by free_pats */
} NamePat;

/* RegularExpression["re"] -> anchored PCRE program, or NULL if `e` is not one. */
static int is_regular_expression(Expr* e, const char** re_out) {
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_RegularExpression &&
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
        out->has_ctx = strchr(e->data.string, '`') != NULL;
        out->glob = e->data.string;
        out->prog = NULL;
        return 0;
    }
    if (is_regular_expression(e, &re)) {
        if (!regex_available()) {
            mth_message("Names", "regavail", "regular-expression support is not available; rebuild Mathilda with PCRE2 (USE_REGEX=1).");
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
            mth_message("Names", "regex", "%s", err);
            return -1;
        }
        out->is_regex = 1;
        out->has_ctx = strchr(re, '`') != NULL;
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

/* Context-qualified glob match: split both the pattern and the fully-qualified
 * name `full` at their last backtick and match context against context, short
 * name against short name (see the header comment; issue #109). */
static int glob_matches_qualified(const char* glob, const char* full) {
    const char* pb = strrchr(glob, '`');
    const char* nb = strrchr(full, '`');
    if (!pb || !nb) return 0;
    return wl_glob_match_n(glob, pb + 1, full, nb + 1) &&
           wl_glob_match(pb + 1, nb + 1);
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
 * qualified with their home context -- context_home, the same rule as
 * Context[]. */
static char* full_qualified_name(const char* nm, SymbolDef* def) {
    if (name_has_backtick(nm)) return name_strdup(nm);
    const char* ctx = context_home(nm, def);
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
    int            want_keys;   /* emit symbol-table keys, not display names */
    char**         names;       /* owned emit strings (display names) */
    size_t         count, cap;
} NameCollect;

/* Does a symbol match a single pattern?  `full` is its fully-qualified name
 * and `disp` its display name (see the header comment for the rule). */
static int match_one(const NamePat* p, const char* full, const char* disp) {
    if (p->has_ctx)
        return p->is_regex ? pat_matches(p, full)
                           : glob_matches_qualified(p->glob, full);
    return !name_has_backtick(disp) && pat_matches(p, disp);
}

static void collect_visit(const char* name, SymbolDef* def, void* user) {
    NameCollect* c = (NameCollect*)user;
    char* full = full_qualified_name(name, def);
    if (!full) return;   /* OOM: drop this name rather than crash */
    const char* disp = context_display_name(full);   /* points into full */

    int matched = c->match_all;
    for (int i = 0; !matched && i < c->npats; i++)
        matched = match_one(&c->pats[i], full, disp);
    char* emit = matched ? name_strdup(c->want_keys ? name : disp) : NULL;
    free(full);
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

/* ------------------------------------------------------------------ */
/* Completion candidates (the REPL's Tab key, src/term/term_complete.c) */
/* ------------------------------------------------------------------ */

typedef struct {
    const char* prefix;
    size_t      plen;
    int         qualified;   /* the prefix has a backtick: match full names */
    char**      names;
    size_t      count, cap;
} NameComplete;

/* A Module temporary: name$digits (x$123). */
static int is_module_temporary(const char* s) {
    const char* d = strrchr(s, '$');
    if (!d || d == s || !d[1]) return 0;
    for (const char* p = d + 1; *p; p++) if (*p < '0' || *p > '9') return 0;
    return 1;
}

/* Something a user could mean: a builtin, attributes, any rule or value,
 * options, or a System` / Mathilda` symbol. A name the parser merely
 * interned (a typo, a bare unbound symbol) is not offered. */
static int is_defined_for_completion(const SymbolDef* d) {
    return d->builtin_func || d->attributes || d->own_values || d->down_values
        || d->sub_values || d->up_values || d->default_options || d->is_system;
}

static void complete_visit(const char* name, SymbolDef* def, void* user) {
    NameComplete* c = (NameComplete*)user;
    if (!is_defined_for_completion(def)) return;
    char* full = full_qualified_name(name, def);
    if (!full) return;
    const char* disp = context_display_name(full);   /* points into full */
    const char* cand = c->qualified ? full : disp;
    if (strncmp(cand, c->prefix, c->plen) != 0
        || (!c->qualified && name_has_backtick(disp))
        || is_module_temporary(disp)) {
        free(full);
        return;
    }
    char* emit = name_strdup(cand);
    free(full);
    if (!emit) return;
    if (c->count == c->cap) {
        size_t nc = c->cap ? c->cap * 2 : 64;
        char** np = realloc(c->names, nc * sizeof(char*));
        if (!np) { free(emit); return; }
        c->names = np;
        c->cap = nc;
    }
    c->names[c->count++] = emit;
}

static int cstr_cmp(const void* a, const void* b) {
    return strcmp(*(char* const*)a, *(char* const*)b);
}

char** names_complete(const char* prefix, size_t* n) {
    *n = 0;
    if (!prefix) return NULL;
    NameComplete c = { prefix, strlen(prefix), name_has_backtick(prefix), NULL, 0, 0 };
    symtab_for_each(complete_visit, &c);
    if (c.count == 0) { free(c.names); return NULL; }
    qsort(c.names, c.count, sizeof(char*), cstr_cmp);
    size_t k = 0;
    for (size_t i = 0; i < c.count; i++) {
        if (k > 0 && strcmp(c.names[k - 1], c.names[i]) == 0) free(c.names[i]);
        else c.names[k++] = c.names[i];
    }
    *n = k;
    return c.names;
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
            arg->data.function.head->data.symbol.name == SYM_List) {
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

    NameCollect c = { pats, npats, match_all, 0, NULL, 0, 0 };
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

/* The symbol-table keys of every symbol a single Names-style pattern
 * (a string with * and @, or a RegularExpression) matches, resolved exactly
 * as Names[pattern] resolves it, in canonical order.  Returns a malloc'd
 * array of malloc'd strings (free with names_free_keys) and sets *n; NULL
 * with *n = 0 when nothing matches or the pattern is not one. */
char** names_match_keys(Expr* pattern, size_t* n) {
    *n = 0;
    NamePat pat;
    memset(&pat, 0, sizeof pat);
    if (build_pat(pattern, &pat) != 0) return NULL;
    NameCollect c = { &pat, 1, 0, 1, NULL, 0, 0 };
    symtab_for_each(collect_visit, &c);
    if (pat.is_regex) regex_free(pat.prog);
    if (c.count == 0) { free(c.names); return NULL; }
    /* Deterministic order for the messages a caller emits per name. */
    for (size_t i = 1; i < c.count; i++)
        for (size_t j = i; j > 0 && strcmp(c.names[j - 1], c.names[j]) > 0; j--) {
            char* t = c.names[j]; c.names[j] = c.names[j - 1]; c.names[j - 1] = t;
        }
    *n = c.count;
    return c.names;
}

void names_free_keys(char** keys, size_t n) {
    for (size_t i = 0; i < n; i++) free(keys[i]);
    free(keys);
}

/* ------------------------------------------------------------------ */

void names_init(void) {
    symtab_add_builtin("Names", builtin_names);
    symtab_get_def("Names")->attributes |= ATTR_PROTECTED;
}
