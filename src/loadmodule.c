/*
 * loadmodule.c -- runtime loading of Mathilda (.m) source modules.
 *
 * See loadmodule.h for the contract. The file-reading loop here is the same
 * one Get[] uses (builtin_get in readwrite.c is now a thin wrapper over
 * mathilda_run_file); mathilda_load_module adds working-directory-independent
 * path resolution and load-once bookkeeping on top.
 *
 * Path resolution (mathilda_resolve_internal) is deliberately CWD-independent:
 * a relocated or installed binary must still find its bundled src/internal
 * tree — and so must an installed shared libmathilda loaded into some other
 * program (a Python extension, a Jupyter kernel), which is located through
 * dladdr. On glibc, readlink("/proc/self/exe") needs the POSIX feature-test
 * macro and dladdr/Dl_info need _GNU_SOURCE to be visible under -std=c99;
 * define it before any system header.
 */

#if defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif

#include "loadmodule.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "parse.h"
#include "eval.h"
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "quit.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__linux__)
#include <unistd.h>          /* readlink */
#elif defined(__APPLE__)
#include <mach-o/dyld.h>     /* _NSGetExecutablePath */
#include <stdint.h>          /* uint32_t */
#elif defined(_WIN32)
#include <windows.h>         /* GetModuleFileNameA, GetModuleHandleExA */
#endif
#if !defined(_WIN32) && (defined(__linux__) || defined(__APPLE__))
#include <dlfcn.h>           /* dladdr (lm_self_dir) */
#endif

/* Longest module path we build; comfortably above any real filesystem path. */
#define LM_PATH_MAX 2048

/* Set $InputFileName's OwnValue to `val` (borrowed). */
static void input_file_set(Expr* val) {
    Expr* sym = expr_new_symbol("$InputFileName");
    symtab_add_own_value("$InputFileName", sym, val);
    expr_free(sym);
}

Expr* mth_input_file_enter(const char* path) {
    Rule* r = symtab_get_own_values("$InputFileName");
    Expr* prev = (r && r->replacement) ? expr_copy(r->replacement)
                                       : expr_new_string("");
#ifdef _WIN32
    char* abs = _fullpath(NULL, path, 0);
#else
    char* abs = realpath(path, NULL);
#endif
    Expr* val = expr_new_string(abs ? abs : path);   /* fall back to as given */
    free(abs);
    input_file_set(val);
    expr_free(val);
    return prev;
}

void mth_input_file_leave(Expr* prev) {
    input_file_set(prev);
    expr_free(prev);
}

/* Skip whitespace and (* comments *) -- what may legitimately trail the last
 * statement of a file. */
static const char* lm_skip_blanks(const char* p) {
    for (;;) {
        while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r' || *p == '\f' || *p == '\v') p++;
        if (p[0] == '(' && p[1] == '*') {
            int depth = 1;
            p += 2;
            while (*p && depth > 0) {
                if (p[0] == '(' && p[1] == '*') { depth++; p += 2; }
                else if (p[0] == '*' && p[1] == ')') { depth--; p += 2; }
                else p++;
            }
            continue;
        }
        return p;
    }
}

/* Syntax::sntx for the statement starting at err in buffer: the rest of its
 * line, and its 1-based line number in the file. */
static void lm_syntax_message(const char* buffer, const char* err, const char* path) {
    int line = 1;
    for (const char* q = buffer; q < err; q++) if (*q == '\n') line++;
    const char* eol = err;
    while (*eol && *eol != '\n' && *eol != '\r') eol++;
    int len = (int)(eol - err);
    if (len > 200) len = 200;
    mth_message("Syntax", "sntx", "Invalid syntax in or before \"%.*s\" (line %d of \"%s\").",
                len, err, line, path);
}

Expr* mathilda_run_file(const char* path, int* opened) {
    FILE* fp = fopen(path, "r");
    if (!fp) {
        if (opened) *opened = 0;
        return NULL;
    }
    if (opened) *opened = 1;

    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    if (fsize < 0) {
        fclose(fp);
        return expr_new_symbol(SYM_DollarFailed);
    }

    char* buffer = malloc((size_t)fsize + 1);
    if (!buffer) {
        fclose(fp);
        return expr_new_symbol(SYM_DollarFailed);
    }

    size_t read_len = fread(buffer, 1, (size_t)fsize, fp);
    buffer[read_len] = '\0';
    fclose(fp);

    Expr* last_eval = expr_new_symbol(SYM_Null);
    Expr* prev_input = mth_input_file_enter(path);
    const char* ptr = buffer;
    while (*ptr != '\0') {
        const char* stmt_start = ptr;
        int prev_quiet = parse_set_quiet(1);
        Expr* parsed = parse_next_expression(&ptr);
        parse_set_quiet(prev_quiet);
        if (!parsed) {
            /* End of input, or a syntax error: either way reading stops (as
             * the reference does), but an error is reported, through the
             * message funnel, and the result is Null. */
            const char* err = lm_skip_blanks(stmt_start);
            if (*err != '\0') {
                lm_syntax_message(buffer, err, path);
                expr_free(last_eval);
                last_eval = expr_new_symbol(SYM_Null);
            }
            break;
        }
        Expr* evaluated = evaluate(parsed);
        if (evaluated) {
            expr_free(last_eval);
            last_eval = evaluated;
        }
        expr_free(parsed);
        /* An abort stops the read and travels on out of Get (its value is
         * the in-flight Abort[] sentinel). */
        if (eval_is_inflight_abort(last_eval)) break;
        /* Quit[] inside the file ends the session: stop reading so the
         * top-level driver can shut down without running the rest. */
        if (mth_quit_requested()) break;
    }
    mth_input_file_leave(prev_input);

    free(buffer);
    return last_eval;
}

/* Modules loaded so far, keyed by the relpath requested. Bounded; a kernel
 * has a small fixed set of internal modules. */
#define LM_MAX_LOADED 256
#define LM_RELPATH_MAX 256
static char lm_loaded[LM_MAX_LOADED][LM_RELPATH_MAX];
static int  lm_loaded_count = 0;
static unsigned lm_loaded_gen = 0;   /* symtab_generation() the memo belongs to */

/* The memo describes definitions in the CURRENT symbol table; a symtab_init
 * since it was filled discarded them, so it is stale. */
static void lm_sync_generation(void) {
    if (lm_loaded_gen != symtab_generation()) {
        lm_loaded_count = 0;
        lm_loaded_gen = symtab_generation();
    }
}

static int lm_already_loaded(const char* relpath) {
    lm_sync_generation();
    for (int i = 0; i < lm_loaded_count; i++) {
        if (strcmp(lm_loaded[i], relpath) == 0) return 1;
    }
    return 0;
}

static void lm_mark_loaded(const char* relpath) {
    lm_sync_generation();
    if (lm_loaded_count >= LM_MAX_LOADED) return; /* silently cap */
    size_t n = strlen(relpath);
    if (n >= LM_RELPATH_MAX) return;
    memcpy(lm_loaded[lm_loaded_count], relpath, n + 1);
    lm_loaded_count++;
}

/* Directory (with trailing '/') of the winning base once one is found, so
 * subsequent lookups skip the search and go straight to the right tree. */
static char lm_base[LM_PATH_MAX];
static int  lm_base_known = 0;

/* Trim the file path `path` (modified in place) to its directory, keeping the
 * trailing '/', and copy that into `buf`. Returns 1 on success. */
static int lm_dir_of(char* path, char* buf, size_t bufsz) {
    char* slash = strrchr(path, '/');
    if (!slash) return 0;
    slash[1] = '\0';
    if (strlen(path) >= bufsz) return 0;
    memcpy(buf, path, strlen(path) + 1);
    return 1;
}

/* Write the directory containing the running executable (with a trailing '/')
 * into `buf`. Returns 1 on success, 0 where unsupported or on error. Purely a
 * hint for resolution; callers fall back to $MATHILDA_HOME and the CWD ladder
 * when this returns 0. */
static int lm_exe_dir(char* buf, size_t bufsz) {
    char path[LM_PATH_MAX];
    path[0] = '\0';

#if defined(__linux__)
    long n = readlink("/proc/self/exe", path, sizeof(path) - 1);
    if (n <= 0) return 0;
    path[n] = '\0';
#elif defined(__APPLE__)
    uint32_t sz = sizeof(path);
    if (_NSGetExecutablePath(path, &sz) != 0) return 0; /* buffer too small */
    path[sizeof(path) - 1] = '\0';
#elif defined(_WIN32)
    DWORD n = GetModuleFileNameA(NULL, path, (DWORD)sizeof(path));
    if (n == 0 || n >= sizeof(path)) return 0;
    path[n] = '\0';
    /* Normalise backslashes so the '/'-based split below works uniformly. */
    for (char* p = path; *p; p++) if (*p == '\\') *p = '/';
#else
    (void)path;
    return 0;
#endif

    return lm_dir_of(path, buf, bufsz);
}

/* Write the directory of the loaded image that contains THIS code (with a
 * trailing '/') into `buf`. For the Mathilda executable that is the executable
 * itself; for an installed shared libmathilda it is <prefix>/lib/, whatever
 * program loaded it — which is what lets the library find
 * <prefix>/share/mathilda/internal/ with no $MATHILDA_HOME. Statically linked
 * into some other image it names that image, and the lookup simply misses.
 * Returns 1 on success, 0 where unsupported or on error. */
static int lm_self_dir(char* buf, size_t bufsz) {
    char path[LM_PATH_MAX];
    path[0] = '\0';

#if !defined(_WIN32) && (defined(__linux__) || defined(__APPLE__))
    Dl_info info;
    if (!dladdr((void*)&mathilda_resolve_internal, &info) || !info.dli_fname)
        return 0;
    if (strlen(info.dli_fname) >= sizeof(path)) return 0;
    memcpy(path, info.dli_fname, strlen(info.dli_fname) + 1);
#elif defined(_WIN32)
    HMODULE mod = NULL;
    if (!GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                            GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                            (LPCSTR)(void*)&mathilda_resolve_internal, &mod))
        return 0;
    DWORD n = GetModuleFileNameA(mod, path, (DWORD)sizeof(path));
    if (n == 0 || n >= sizeof(path)) return 0;
    path[n] = '\0';
    for (char* p = path; *p; p++) if (*p == '\\') *p = '/';
#else
    (void)path;
    return 0;
#endif

    return lm_dir_of(path, buf, bufsz);
}

/* Try base directory `base` (which must end with '/'): if <base><relpath> is a
 * readable file, cache the base, write the full path into `out`, return 1. */
static int lm_try_base(const char* base, const char* relpath,
                       char* out, size_t outsz) {
    char full[LM_PATH_MAX];
    if ((size_t)snprintf(full, sizeof(full), "%s%s", base, relpath)
            >= sizeof(full))
        return 0; /* would truncate; not a real path */
    FILE* fp = fopen(full, "r");
    if (!fp) return 0;
    fclose(fp);
    snprintf(lm_base, sizeof(lm_base), "%s", base);
    lm_base_known = 1;
    if (strlen(full) < outsz) memcpy(out, full, strlen(full) + 1);
    return 1;
}

int mathilda_resolve_internal(const char* relpath, char* out, size_t outsz) {
    if (!relpath || !*relpath || !out || outsz == 0) return 0;

    /* Fast path: reuse the directory that satisfied an earlier lookup. */
    if (lm_base_known) {
        char full[LM_PATH_MAX];
        if ((size_t)snprintf(full, sizeof(full), "%s%s", lm_base, relpath)
                < sizeof(full)) {
            FILE* fp = fopen(full, "r");
            if (fp) {
                fclose(fp);
                if (strlen(full) < outsz) memcpy(out, full, strlen(full) + 1);
                return 1;
            }
        }
        /* Cached base no longer serves this file; fall through to a fresh
         * search rather than getting stuck on it. */
    }

    char cand[LM_PATH_MAX];

    /* 1. $MATHILDA_HOME points directly at the src/internal directory. Also
     * accept its parent (the directory CONTAINING internal/): that is what the
     * Jupyter kernelspec shipped before issue #96, and an easy mistake to make,
     * where it silently loaded init.m via the CWD ladder but missed the
     * tables Integrate loads later. */
    const char* home = getenv("MATHILDA_HOME");
    if (home && *home) {
        if ((size_t)snprintf(cand, sizeof(cand), "%s/", home) < sizeof(cand)
                && lm_try_base(cand, relpath, out, outsz))
            return 1;
        if ((size_t)snprintf(cand, sizeof(cand), "%s/internal/", home)
                < sizeof(cand)
                && lm_try_base(cand, relpath, out, outsz))
            return 1;
    }

    /* 2..3. Relative to the executable itself, so a binary that carries its
     * source tree (or an installed share/ layout) works from any CWD. */
    char exedir[LM_PATH_MAX];
    exedir[0] = '\0';
    if (lm_exe_dir(exedir, sizeof(exedir))) {
        /* binary sits in the repo root: <exe>/src/internal/ */
        if ((size_t)snprintf(cand, sizeof(cand), "%ssrc/internal/", exedir)
                < sizeof(cand)
                && lm_try_base(cand, relpath, out, outsz))
            return 1;
        /* installed layout: <prefix>/bin/<exe>, data at <prefix>/share/... */
        if ((size_t)snprintf(cand, sizeof(cand),
                             "%s../share/mathilda/internal/", exedir)
                < sizeof(cand)
                && lm_try_base(cand, relpath, out, outsz))
            return 1;
    }

    /* 3b. Relative to the image holding this code: an installed shared
     * libmathilda at <prefix>/lib/ finds <prefix>/share/mathilda/internal/
     * from inside any host program (issue #88). Identical to step 3 when the
     * code lives in the executable, so only retry when the directory differs. */
    char selfdir[LM_PATH_MAX];
    if (lm_self_dir(selfdir, sizeof(selfdir))
            && strcmp(selfdir, exedir) != 0) {
        if ((size_t)snprintf(cand, sizeof(cand),
                             "%s../share/mathilda/internal/", selfdir)
                < sizeof(cand)
                && lm_try_base(cand, relpath, out, outsz))
            return 1;
    }

    /* 4. Compile-time install prefix (make PREFIX=/usr/local). */
#ifdef MATHILDA_PREFIX
    if ((size_t)snprintf(cand, sizeof(cand),
                         "%s/share/mathilda/internal/", MATHILDA_PREFIX)
            < sizeof(cand)
            && lm_try_base(cand, relpath, out, outsz))
        return 1;
#endif

    /* 5. CWD ladder (repo root, tests/, tests/build/, one deeper). */
    static const char* const bases[] = {
        "src/internal/",
        "../src/internal/",
        "../../src/internal/",
        "../../../src/internal/",
        NULL
    };
    for (int i = 0; bases[i]; i++)
        if (lm_try_base(bases[i], relpath, out, outsz)) return 1;

    return 0;
}

int mathilda_load_module(const char* relpath) {
    if (!relpath || !*relpath) return 0;
    if (lm_already_loaded(relpath)) return 1;

    char full[LM_PATH_MAX];
    if (mathilda_resolve_internal(relpath, full, sizeof(full))) {
        int opened = 0;
        intern_system_scope_push();     /* heads this module defines are System */
        Expr* last = mathilda_run_file(full, &opened);
        intern_system_scope_pop();
        if (last) expr_free(last);
        if (opened) { lm_mark_loaded(relpath); return 1; }
    }

    mth_message("LoadModule", "nofile", "cannot locate src/internal/%s on disk.", relpath);
    return 0;
}

/* LoadModule["relpath"] -> True if the module was found and loaded (or had
 * already been loaded), False otherwise. */
Expr* builtin_loadmodule(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1)
        return NULL;
    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING) return NULL;
    int ok = mathilda_load_module(arg->data.string);
    return expr_new_symbol(ok ? "True" : "False");
}

void loadmodule_init(void) {
    symtab_add_builtin("LoadModule", builtin_loadmodule);
    symtab_get_def("LoadModule")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("LoadModule",
        "LoadModule[\"relpath\"]\n"
        "\tloads the internal Mathilda source module at relpath (relative to\n"
        "\tsrc/internal), resolving the location independently of the current\n"
        "\tworking directory. Each module is loaded at most once. Returns True\n"
        "\tif the module was located and loaded, False otherwise.");
}
