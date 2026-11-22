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

#include "context.h"   /* restore $Context after an aborted read */
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
#include "core.h"      /* tc_register_unwind_counter: an abandoned module read */

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
    /* The context state at entry, put back if an abort stops the read: the
     * file's End[] / EndPackage[] never run, and the session would be left in
     * the package's private context (every later bare symbol a fresh
     * Private` one).  A lazily loaded package (Integrate's) is read inside the
     * caller's TimeConstrained, whose timeout is such an abort.  A package the
     * file loaded TO COMPLETION before the abort keeps its $ContextPath entry
     * (context_reapply_endpackages; see tc_run_guarded_mode). */
    size_t ctx_depth = context_stack_depth();
    size_t pkg_mark = context_endpackage_mark();
    Expr* saved_ctx; Expr* saved_path;
    context_capture(&saved_ctx, &saved_path);
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
    if (eval_is_inflight_abort(last_eval)) {
        context_unwind_to(ctx_depth);
        context_restore(saved_ctx, saved_path);
        context_reapply_endpackages(pkg_mark);
    } else {
        context_release(saved_ctx, saved_path);   /* a file may leave a Begin open */
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

/* Ends `s` with `suffix`? Returns the length of the prefix before it, or -1. */
static long lm_strip_suffix(const char* s, const char* suffix) {
    size_t n = strlen(s), m = strlen(suffix);
    if (m > n || strcmp(s + n - m, suffix) != 0) return -1;
    return (long)(n - m);
}

int mathilda_installation_directory(char* out, size_t outsz) {
    if (!out || outsz == 0) return 0;

    /* The top-level directory is read off the base that actually serves the
     * bootstrap tree, so $InstallationDirectory always names the tree this
     * process is running from (whichever resolution rule found it). */
    char top[LM_PATH_MAX];
    char probe[LM_PATH_MAX];
    if (mathilda_resolve_internal("init.m", probe, sizeof(probe)) && lm_base_known) {
        static const char* const suffixes[] = {
            "src/internal/",            /* source tree: <root>/src/internal/  */
            "share/mathilda/internal/", /* installed:   <prefix>/share/...    */
            "internal/",                /* $MATHILDA_HOME/<dir>/internal/     */
            NULL
        };
        snprintf(top, sizeof(top), "%s", lm_base);
        for (int i = 0; suffixes[i]; i++) {
            long k = lm_strip_suffix(lm_base, suffixes[i]);
            if (k >= 0) { top[k] = '\0'; break; }
        }
        if (!top[0]) snprintf(top, sizeof(top), "./"); /* CWD-relative base */
    } else if (!lm_exe_dir(top, sizeof(top))) {
        return 0;
    }

    /* Canonicalise: absolute, symlink-free, no trailing separator.
     *
     * Let the library allocate the result.  realpath(top, buf) into a fixed
     * buffer requires buf to hold PATH_MAX bytes, and glibc's fortified
     * realpath (_FORTIFY_SOURCE, on by default in distribution gcc at -O2)
     * aborts with "*** buffer overflow detected ***" when the destination is
     * smaller -- whatever the actual path length.  LM_PATH_MAX (2048) is
     * below Linux's PATH_MAX (4096), so every Linux build crashed at startup
     * while $InstallationDirectory was being set (issue #129).  macOS has
     * PATH_MAX 1024 and never saw it. */
#if defined(_WIN32)
    char* abs = _fullpath(NULL, top, 0);
#else
    char* abs = realpath(top, NULL);
#endif
    if (!abs) return 0;
    size_t n = strlen(abs);
    while (n > 1 && (abs[n - 1] == '/' || abs[n - 1] == '\\')) abs[--n] = '\0';
    int ok = n < outsz;
    if (ok) memcpy(out, abs, n + 1);
    free(abs);
    return ok;
}

/* ---- an interrupted read ------------------------------------------------- *
 * A module read can be cut short: by the cooperative Abort[] a caller's
 * TimeConstrained raises on timeout (since v0.6623; the read stops at the
 * statement that saw it), or by a timeout siglongjmp that leaves
 * mathilda_load_module altogether.  Either way the module is HALF defined.  It
 * is never marked loaded, so the next mathilda_load_module reads the whole
 * file -- but the lazy loaders (Integrate's cascade packages, the BVP and
 * Meijer-G modules) keep their own attempt-once flag, keyed on
 * mathilda_module_generation().  That generation moves on whenever a read was
 * interrupted, so a loader whose attempt was cut short asks again instead of
 * running its stage on a half-read package for the rest of the session.
 * (Measured: DSolve`Separable's budget expiring inside the lazy load of
 * ParallelMixedSpecial.m left that stage declining every later integrand at
 * once, and a later, unrelated separable equation shipped an implicit relation
 * where the full cascade would have handed it on to Linearizable.)
 *
 * The jump is detected afterwards: lm_depth is a TimeConstrained unwind
 * counter (restored to its pre-body value by the jump), lm_depth_raw is not,
 * so raw > depth means a read was abandoned mid-way. */
static int      lm_depth = 0;       /* reads in progress (unwind-restored) */
static int      lm_depth_raw = 0;   /* the same, never restored            */
static unsigned lm_interrupted = 0; /* reads cut short so far              */

static void lm_note_jumped_reads(void) {
    if (lm_depth_raw > lm_depth) { lm_interrupted++; lm_depth_raw = lm_depth; }
}

int mathilda_module_depth(void) { return lm_depth; }

unsigned mathilda_module_generation(void) {
    lm_note_jumped_reads();
    return symtab_generation() * 2654435761u + lm_interrupted;
}

int mathilda_load_module(const char* relpath) {
    if (!relpath || !*relpath) return 0;
    lm_note_jumped_reads();
    if (lm_already_loaded(relpath)) return 1;

    char full[LM_PATH_MAX];
    if (mathilda_resolve_internal(relpath, full, sizeof(full))) {
        int opened = 0;
        lm_depth++; lm_depth_raw++;
        intern_system_scope_push();     /* heads this module defines are System */
        Expr* last = mathilda_run_file(full, &opened);
        intern_system_scope_pop();
        lm_depth--; lm_depth_raw--;
        bool aborted = last && eval_is_inflight_abort(last);
        if (last) expr_free(last);
        if (aborted) {
            /* Half read: not marked loaded, and the loaders' generation moves
             * on.  Reported as found (1) -- the file is there, and the abort
             * itself travels on -- so no caller says "cannot locate". */
            lm_interrupted++;
            return 1;
        }
        if (opened) { lm_mark_loaded(relpath); return 1; }
    }

    mth_message("LoadModule", "nofile", "cannot locate src/internal/%s on disk.", relpath);
    return 0;
}

/* ---- lazy loaders -------------------------------------------------------- *
 * The rule every lazy loader keys its attempt on (see MthLazyLoad in the
 * header).  A COMPLETE read is final for its symtab generation: another
 * module's interrupted read moves mathilda_module_generation() on, and keying a
 * finished table on that number re-read it -- for the CRC table, whose last
 * statement Protects IntegrateTable, every one of its rules then failed with a
 * SetDelayed::write message.  Only an attempt that did not finish (file
 * missing) waits for the generation to move; an interrupted one is not recorded
 * at all, so the next call reads it again. */
int mathilda_lazy_due(const MthLazyLoad* st) {
    if (!st->attempted) return 1;
    if (st->complete && st->symgen == symtab_generation()) return 0;
    return st->gen != mathilda_module_generation();
}

void mathilda_lazy_record(MthLazyLoad* st, int complete) {
    st->attempted = 1;
    st->gen = mathilda_module_generation();
    st->complete = complete ? 1 : 0;
    st->symgen = symtab_generation();
}

int mathilda_lazy_load_module(MthLazyLoad* st, const char* relpath) {
    int ok = mathilda_load_module(relpath);
    if (!ok) mathilda_lazy_record(st, 0);            /* not found */
    else if (lm_already_loaded(relpath)) mathilda_lazy_record(st, 1);
    /* else: found but cut short -- unrecorded, so the next call retries */
    return ok;
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
    tc_register_unwind_counter(&lm_depth);   /* see mathilda_module_generation */
    symtab_add_builtin("LoadModule", builtin_loadmodule);
    symtab_get_def("LoadModule")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("LoadModule",
        "LoadModule[\"relpath\"]\n"
        "\tloads the internal Mathilda source module at relpath (relative to\n"
        "\tsrc/internal), resolving the location independently of the current\n"
        "\tworking directory. Each module is loaded at most once. Returns True\n"
        "\tif the module was located and loaded, False otherwise.");

    /* $InstallationDirectory -- fixed for the session (see
     * mathilda_installation_directory). Left unset if nothing resolves. */
    char dir[LM_PATH_MAX];
    if (mathilda_installation_directory(dir, sizeof(dir))) {
        Expr* sym = expr_new_symbol("$InstallationDirectory");
        Expr* val = expr_new_string(dir);
        symtab_add_own_value("$InstallationDirectory", sym, val);
        expr_free(sym);
        expr_free(val);
    }
    symtab_get_def("$InstallationDirectory")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("$InstallationDirectory",
        "$InstallationDirectory\n"
        "\tgives the top-level directory in which your Mathilda installation resides, as the full\n"
        "\tname of the directory as a string: the source-tree root for a binary run from its tree,\n"
        "\tor the install prefix (containing share/mathilda) for an installed build.");
}
