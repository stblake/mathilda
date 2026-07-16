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
 * tree. On glibc, readlink("/proc/self/exe") needs the POSIX feature-test
 * macro to be visible under -std=c99; define it before any system header.
 */

#if defined(__linux__) && !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200112L
#endif

#include "loadmodule.h"
#include "sym_names.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "parse.h"
#include "eval.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__linux__)
#include <unistd.h>          /* readlink */
#elif defined(__APPLE__)
#include <mach-o/dyld.h>     /* _NSGetExecutablePath */
#include <stdint.h>          /* uint32_t */
#elif defined(_WIN32)
#include <windows.h>         /* GetModuleFileNameA */
#endif

/* Longest module path we build; comfortably above any real filesystem path. */
#define LM_PATH_MAX 2048

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
    const char* ptr = buffer;
    while (*ptr != '\0') {
        Expr* parsed = parse_next_expression(&ptr);
        if (!parsed) break; /* EOF or unparseable: stop. */
        Expr* evaluated = evaluate(parsed);
        if (evaluated) {
            expr_free(last_eval);
            last_eval = evaluated;
        }
        expr_free(parsed);
    }

    free(buffer);
    return last_eval;
}

/* Modules loaded so far, keyed by the relpath requested. Bounded; a kernel
 * has a small fixed set of internal modules. */
#define LM_MAX_LOADED 256
#define LM_RELPATH_MAX 256
static char lm_loaded[LM_MAX_LOADED][LM_RELPATH_MAX];
static int  lm_loaded_count = 0;

static int lm_already_loaded(const char* relpath) {
    for (int i = 0; i < lm_loaded_count; i++) {
        if (strcmp(lm_loaded[i], relpath) == 0) return 1;
    }
    return 0;
}

static void lm_mark_loaded(const char* relpath) {
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

    /* Trim to the directory: keep everything up to and including the last '/'. */
    char* slash = strrchr(path, '/');
    if (!slash) return 0;
    slash[1] = '\0';
    if (strlen(path) >= bufsz) return 0;
    memcpy(buf, path, strlen(path) + 1);
    return 1;
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

    /* 1. $MATHILDA_HOME points directly at the src/internal directory. */
    const char* home = getenv("MATHILDA_HOME");
    if (home && *home) {
        if ((size_t)snprintf(cand, sizeof(cand), "%s/", home) < sizeof(cand)
                && lm_try_base(cand, relpath, out, outsz))
            return 1;
    }

    /* 2..3. Relative to the executable itself, so a binary that carries its
     * source tree (or an installed share/ layout) works from any CWD. */
    char exedir[LM_PATH_MAX];
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
        Expr* last = mathilda_run_file(full, &opened);
        if (last) expr_free(last);
        if (opened) { lm_mark_loaded(relpath); return 1; }
    }

    fprintf(stderr,
        "LoadModule::nofile: cannot locate src/internal/%s on disk.\n",
        relpath);
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
