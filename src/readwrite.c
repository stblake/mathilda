/*
 * readwrite.c - File I/O builtins (Get, Put).
 *
 * Get reads Mathilda source from a file and evaluates each expression,
 * returning the last value (used by the REPL bootstrap to load
 * the internal .m initialization files).  Its argument is either a file
 * name or a context name ("Pkg`", "A`B`"); both are searched for along
 * $Path (see readwrite_resolve), which Needs relies on.
 *
 * Put writes one or more expressions to a file in InputForm so the
 * output can be read back with Get.  The parser also recognises the
 * infix shorthand `expr >> "file"` and lowers it to `Put[expr, "file"]`.
 */

#include "readwrite.h"
#include "sym_names.h"
#include "loadmodule.h"
#include "expr.h"
#include "symtab.h"
#include "print.h"
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------- Get path search ---------------------------- */

static int rw_readable(const char* path) {
    FILE* fp = fopen(path, "r");
    if (!fp) return 0;
    fclose(fp);
    return 1;
}

/* "<dir>/<rel>" in a fresh malloc'd string. */
static char* rw_join(const char* dir, const char* rel) {
    size_t nd = strlen(dir), nr = strlen(rel);
    int slash = nd > 0 && dir[nd - 1] != '/' && dir[nd - 1] != '\\';
    char* out = (char*)malloc(nd + (size_t)slash + nr + 1);
    memcpy(out, dir, nd);
    if (slash) out[nd] = '/';
    memcpy(out + nd + (size_t)slash, rel, nr + 1);
    return out;
}

static int rw_is_absolute(const char* p) {
    if (p[0] == '/' || p[0] == '\\') return 1;
#ifdef _WIN32
    if (((p[0] >= 'A' && p[0] <= 'Z') || (p[0] >= 'a' && p[0] <= 'z')) && p[1] == ':')
        return 1;
#endif
    return 0;
}

/* Try <dir>/<rel> for every string entry of $Path, in order.  A $Path that is
 * not a List (cleared, or assigned something odd) searches just ".". */
static char* rw_search_path(const char* rel) {
    Rule* r = symtab_get_own_values("$Path");
    Expr* path = r ? r->replacement : NULL;
    if (!path || path->type != EXPR_FUNCTION ||
        path->data.function.head->type != EXPR_SYMBOL ||
        path->data.function.head->data.symbol.name != SYM_List) {
        char* cand = rw_join(".", rel);
        if (rw_readable(cand)) return cand;
        free(cand);
        return NULL;
    }
    for (size_t i = 0; i < path->data.function.arg_count; i++) {
        Expr* d = path->data.function.args[i];
        if (d->type != EXPR_STRING || d->data.string[0] == '\0') continue;
        char* cand = rw_join(d->data.string, rel);
        if (rw_readable(cand)) return cand;
        free(cand);
    }
    return NULL;
}

char* readwrite_resolve(const char* spec) {
    if (!spec || !*spec) return NULL;
    size_t n = strlen(spec);

    if (spec[n - 1] == '`') {
        /* Context name: "A`B`" -> stem "A/B", tried as a file and as a
         * package directory, along $Path.  A leading backtick is not a
         * loadable context. */
        if (spec[0] == '`') return NULL;
        static const char* const suffixes[] = {
            ".m", ".wl", "/init.m", "/init.wl", "/Kernel/init.m", "/Kernel/init.wl"
        };
        char* stem = (char*)malloc(n);
        memcpy(stem, spec, n - 1);
        stem[n - 1] = '\0';
        for (char* p = stem; *p; p++) if (*p == '`') *p = '/';
        char* found = NULL;
        for (size_t k = 0; !found && k < sizeof(suffixes) / sizeof(suffixes[0]); k++) {
            size_t ns = strlen(stem), nx = strlen(suffixes[k]);
            char* rel = (char*)malloc(ns + nx + 1);
            memcpy(rel, stem, ns);
            memcpy(rel + ns, suffixes[k], nx + 1);
            found = rw_search_path(rel);
            free(rel);
        }
        free(stem);
        return found;
    }

    /* File name: as given first, then (when relative) along $Path. */
    if (rw_readable(spec)) {
        char* out = (char*)malloc(n + 1);
        memcpy(out, spec, n + 1);
        return out;
    }
    if (rw_is_absolute(spec)) return NULL;
    return rw_search_path(spec);
}

Expr* readwrite_get(const char* spec) {
    char* path = readwrite_resolve(spec);
    int opened = 0;
    Expr* last_eval = path ? mathilda_run_file(path, &opened) : NULL;
    free(path);
    if (!opened) {
        mth_message("Get", "noopen", "Cannot open %s.", spec);
        return expr_new_symbol(SYM_DollarFailed);
    }
    /* mathilda_run_file returns Null for an empty file; never NULL when
     * opened succeeded. */
    return last_eval ? last_eval : expr_new_symbol(SYM_Null);
}

Expr* builtin_get(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;

    Expr* file_arg = res->data.function.args[0];
    if (file_arg->type != EXPR_STRING) return NULL;
    return readwrite_get(file_arg->data.string);
}

/* Shared implementation for Put / PutAppend.  `mode` selects fopen
 * mode: "w" truncates ("Put"), "a" appends ("PutAppend").  `who` names
 * the calling builtin for diagnostics. */
static Expr* put_common(Expr* res, const char* mode, const char* who) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return NULL;

    Expr* last = res->data.function.args[argc - 1];
    if (last->type != EXPR_STRING) return NULL;

    const char* filename = last->data.string;
    FILE* fp = fopen(filename, mode);
    if (!fp) {
        printf("%s::noopen: Cannot open %s.\n", who, filename);
        return expr_new_symbol(SYM_DollarFailed);
    }

    /* When only the filename is supplied (Put["file"]) the loop runs
     * zero times — leaving an empty file in "w" mode and leaving the
     * file untouched (besides creating it if absent) in "a" mode. */
    for (size_t i = 0; i + 1 < argc; i++) {
        /* Put writes InputForm, as the reference: a/100 + 5*x, (x*y)/2 */
        Expr* inner = expr_copy(res->data.function.args[i]);
        Expr* w = expr_new_function(expr_new_symbol(SYM_InputForm), &inner, 1);
        char* s = expr_to_string(w);
        expr_free(w);
        if (s) {
            fputs(s, fp);
            free(s);
        }
        fputc('\n', fp);
    }

    fclose(fp);
    return expr_new_symbol(SYM_Null);
}

/* Put[expr_1, ..., expr_n, "file"] truncates "file" and writes each
 * expr_i as standard-form text followed by a newline.  The single-arg
 * form Put["file"] creates an empty file. */
Expr* builtin_put(Expr* res) {
    return put_common(res, "w", "Put");
}

/* PutAppend[expr_1, ..., expr_n, "file"] appends each expr_i to the
 * end of "file", creating the file if it does not exist.  Single-arg
 * form PutAppend["file"] creates an empty file if it doesn't exist
 * but does not truncate an existing file. */
Expr* builtin_putappend(Expr* res) {
    return put_common(res, "a", "PutAppend");
}

void readwrite_init(void) {
    /* Docstrings for Get / Put / PutAppend live in info_init() — keep
     * this module focused on builtin registration and attributes. */
    symtab_add_builtin("Get", builtin_get);
    symtab_get_def("Get")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("Put", builtin_put);
    symtab_get_def("Put")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("PutAppend", builtin_putappend);
    symtab_get_def("PutAppend")->attributes |= ATTR_PROTECTED;

    /* $Path: the directories Get (and so Needs) searches, in order.  An
     * ordinary assignable variable -- AppendTo[$Path, dir] extends it. */
    Expr* psym = expr_new_symbol(SYM_DollarPath);
    Expr** dirs = (Expr**)malloc(sizeof(Expr*));
    dirs[0] = expr_new_string(".");
    Expr* pval = expr_new_function(expr_new_symbol(SYM_List), dirs, 1);
    free(dirs);
    symtab_add_own_value(SYM_DollarPath, psym, pval);
    expr_free(psym);
    expr_free(pval);
}
