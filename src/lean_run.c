/*
 * lean_run.c -- the Lean 4 process runner and its builtins (see lean_run.h).
 *
 *   Lean`Internal`RunProcess[{exe, args...}, cwd | None]
 *       -> <|"ExitCode" -> n, "StandardOutput" -> s, "StandardError" -> s|>
 *       Internal (used for `lean --version` and by tests); there is no public
 *       RunProcess.
 *
 *   Lean`Internal`VerifySource[src, exe, project | None, keep, replay | None]
 *       Writes src to a fresh temp .lean file, runs `exe --json FILE` (with a
 *       project: `cd project && <lake> env exe --json FILE`, <lake> being the
 *       lake binary beside exe), parses the NDJSON message stream and the
 *       `#print axioms` reports, and returns
 *         <|"Success" -> b, "ExitCode" -> n, "Messages" -> {<|"Severity",
 *           "Kind", "Line", "Text"|>...}, "Axioms" -> <|thm -> {ax...}|>,
 *           "AxiomReports" -> {...}, "Timing" -> secs, "Command" -> "..."|>
 *       Success is the file-level verdict: exit code 0, no error message, no
 *       `sorry` (kind hasSorry), and every reported axiom list a subset of
 *       {propext, Classical.choice, Quot.sound} -- so sorryAx and
 *       Lean.ofReduceBool (native_decide) fail it.  Per-theorem attribution
 *       lives in Lean.m.  With `replay` (the path of Replay.lean) the source
 *       is also compiled to MathildaCert.olean and, after a clean run,
 *       replayed by `exe --run replay <olean>` in a fresh process; Success
 *       then also needs "Replay" -> True and standard axioms for every
 *       declaration of the module ("ReplayAxioms").
 *
 * Each run gets a private directory made by mkdtemp (mode 0700) under
 * $TemporaryDirectory, holding MathildaCert.lean (created with O_CREAT | O_EXCL |
 * O_NOFOLLOW, mode 0600), MathildaCert.out / .err and check.out / .err (created by the shell's
 * redirections under set -C).  No other user can plant a file or a symlink
 * there, so the .out parsed afterwards is the one Lean's run wrote (its
 * absence is reported as a failure).  The directory and its files are
 * removed unless keep is True (then the .lean file and the directory are kept and
 * the path is returned as "SourceFile").  _WIN32 has no mkdtemp in its C
 * runtime and keeps a best-effort probed fresh name in the temp directory.
 * mkdtemp / open / fdopen / rmdir need POSIX.1-2008, hence the feature-test
 * macro below (set before any header).
 *
 * "AxiomReports" lists every `#print axioms` report with its source line
 * (<|"Theorem", "Line", "Axioms"|>), so that Lean.m can accept only the
 * report printed by its own trusted footer line.
 */
#if defined(__APPLE__)
/* Darwin declares mkdtemp only in its full (non-strict) namespace, which
 * _POSIX_C_SOURCE alone would hide. */
#  ifndef _DARWIN_C_SOURCE
#  define _DARWIN_C_SOURCE
#  endif
#endif
#ifndef _WIN32
#  ifndef _POSIX_C_SOURCE
#  define _POSIX_C_SOURCE 200809L
#  endif
#endif

#include "lean_run.h"

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef _WIN32
#include <dirent.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

#include "assoc.h"
#include "attr.h"
#include "eval.h"
#include "message.h"
#include "mtime.h"
#include "print_lean.h"
#include "sym_names.h"
#include "symtab.h"

/* ------------------------------------------------------------------ */
/* Growable string                                                     */
/* ------------------------------------------------------------------ */

typedef struct { char* s; size_t n, cap; } RBuf;

static void rb_putn(RBuf* b, const char* s, size_t n) {
    if (b->n + n + 1 > b->cap) {
        size_t nc = b->cap ? b->cap : 128;
        while (nc < b->n + n + 1) nc *= 2;
        char* p = realloc(b->s, nc);
        if (!p) abort();
        b->s = p;
        b->cap = nc;
    }
    memcpy(b->s + b->n, s, n);
    b->n += n;
    b->s[b->n] = '\0';
}

static void rb_puts(RBuf* b, const char* s) { rb_putn(b, s, strlen(s)); }

static char* rb_take(RBuf* b) {
    if (!b->s) { char* e = malloc(1); if (!e) abort(); e[0] = '\0'; return e; }
    return b->s;
}

static char* str_copy(const char* s) {
    size_t n = strlen(s);
    char* r = malloc(n + 1);
    if (!r) abort();
    memcpy(r, s, n + 1);
    return r;
}

/* ------------------------------------------------------------------ */
/* Scratch space: a private per-run directory                          */
/* ------------------------------------------------------------------ */

#ifdef _WIN32
#define LEAN_SEP "\\"
#else
#define LEAN_SEP "/"
#endif

/* $TemporaryDirectory (its OwnValue), else "." -- caller frees. */
static char* temp_directory(void) {
    Expr* s = expr_new_symbol("$TemporaryDirectory");
    Expr* v = evaluate(s);
    expr_free(s);
    char* r = (v && v->type == EXPR_STRING) ? str_copy(v->data.string) : str_copy(".");
    if (v) expr_free(v);
    return r;
}

static bool file_exists(const char* path) {
    FILE* f = fopen(path, "r");
    if (!f) return false;
    fclose(f);
    return true;
}

/* Where one run keeps its files: <stem>.lean, <stem>.out, <stem>.err.
 *
 * POSIX: `dir` is a fresh directory made by mkdtemp (mode 0700, owned by this
 * process's user), and stem is <dir>/MathildaCert (the module name of the
 * compiled .olean).  Nobody else can create, replace or
 * symlink a file in it, so a file found there after the run was written by
 * this run.  _WIN32: there is no mkdtemp in the C runtime; `dir` is NULL and
 * stem is a fresh name in the temp directory, probed as absent (the old,
 * best-effort scheme -- this branch is compiled but not exercised). */
typedef struct { char* dir; char* stem; } Scratch;

static bool scratch_open(Scratch* sc) {
    sc->dir = sc->stem = NULL;
    char* tmp = temp_directory();
    RBuf b = {0};
    rb_puts(&b, tmp);
    free(tmp);
#ifndef _WIN32
    rb_puts(&b, LEAN_SEP "mathilda-lean-XXXXXX");
    if (!mkdtemp(b.s)) { free(b.s); return false; }
    sc->dir = b.s;
    RBuf st = {0};
    rb_puts(&st, sc->dir);
    rb_puts(&st, LEAN_SEP "MathildaCert");
    sc->stem = rb_take(&st);
    return true;
#else
    static unsigned long counter = 0;
    rb_puts(&b, LEAN_SEP);
    for (int tries = 0; tries < 1000; tries++) {
        RBuf c = {0};
        char t[96];
        snprintf(t, sizeof t, "mathilda-lean-%lu-%lu-%lu", (unsigned long)time(NULL),
                 (unsigned long)clock(), counter++);
        rb_puts(&c, b.s);
        rb_puts(&c, t);
        bool clash = false;
        const char* ext[3] = { ".lean", ".out", ".err" };
        for (int i = 0; i < 3 && !clash; i++) {
            RBuf p = {0};
            rb_puts(&p, c.s);
            rb_puts(&p, ext[i]);
            clash = file_exists(p.s);
            free(p.s);
        }
        if (!clash) { sc->stem = c.s; break; }
        free(c.s);
    }
    free(b.s);
    return sc->stem != NULL;
#endif
}

static char* with_ext(const char* stem, const char* ext) {
    RBuf b = {0};
    rb_puts(&b, stem);
    rb_puts(&b, ext);
    return rb_take(&b);
}

/* Remove the run's files (all but .lean when keep_source) and, when nothing
 * is kept, its directory.  Every file of the private directory goes, not
 * only the ones this run expects: first the known names (stem + suffix, the
 * only step on _WIN32), then, on POSIX, every remaining entry of OUR private
 * directory -- each is unlinked (a symlink is removed, never followed), a
 * subdirectory is only rmdir'ed (never recursed into), and nothing outside
 * the directory is touched. */
static void scratch_close(Scratch* sc, bool keep_source) {
    /* The source (.lean) last: it is the one entry kept by keep_source.  The
     * compiled module (.olean and its optional parts) never is. */
    const char* ext[] = { ".out", ".err", ".olean", ".olean.server", ".olean.private",
                          ".ilean", ".ir", ".ir.sig", ".c", ".hash", ".trace", ".forged",
                          "-check.out", "-check.err", ".lean" };
    const int n = (int)(sizeof ext / sizeof ext[0]);
    for (int i = 0; i < (keep_source ? n - 1 : n); i++) {
        char* p = with_ext(sc->stem, ext[i]);
        remove(p);
        free(p);
    }
#ifndef _WIN32
    if (sc->dir) {
        const char* checks[] = { LEAN_SEP "check.out", LEAN_SEP "check.err" };
        for (int i = 0; i < 2; i++) {
            char* p = with_ext(sc->dir, checks[i]);
            remove(p);
            free(p);
        }
        DIR* d = opendir(sc->dir);
        if (d) {
            const char* keep = keep_source ? "MathildaCert.lean" : NULL;
            struct dirent* de;
            while ((de = readdir(d)) != NULL) {
                const char* nm = de->d_name;
                if (strcmp(nm, ".") == 0 || strcmp(nm, "..") == 0) continue;
                if (keep && strcmp(nm, keep) == 0) continue;
                if (strchr(nm, '/')) continue;
                char* p = with_ext(sc->dir, LEAN_SEP);
                char* q = with_ext(p, nm);
                free(p);
                struct stat st;
                if (lstat(q, &st) == 0) {
                    if (S_ISDIR(st.st_mode)) rmdir(q);
                    else unlink(q);
                }
                free(q);
            }
            closedir(d);
        }
        if (!keep_source) rmdir(sc->dir);
    }
#endif
    free(sc->dir);
    free(sc->stem);
    sc->dir = sc->stem = NULL;
}

/* After the elaborating run, the private directory may hold exactly the
 * files this run created: MathildaCert.lean (written by us), .out / .err (the
 * shell's redirections) and .olean (Lean's -o).  Anything else (a forged
 * .olean.server / .olean.private next to the module, a planted helper) means
 * the run did something it should not have: NULL when the directory is as
 * expected, else the name of the first unexpected entry (owned). */
static char* unexpected_entry(const Scratch* sc) {
#ifndef _WIN32
    if (!sc->dir) return NULL;
    DIR* d = opendir(sc->dir);
    if (!d) return str_copy("(the private directory cannot be read)");
    const char* allowed[] = { ".", "..", "MathildaCert.lean", "MathildaCert.out",
                              "MathildaCert.err", "MathildaCert.olean", NULL };
    char* bad = NULL;
    struct dirent* de;
    while (!bad && (de = readdir(d)) != NULL) {
        bool ok = false;
        for (int i = 0; allowed[i]; i++)
            if (strcmp(de->d_name, allowed[i]) == 0) ok = true;
        if (!ok) bad = str_copy(de->d_name);
    }
    closedir(d);
    if (bad) return bad;
    /* The module itself must be a regular file, not a symlink planted in
     * its place. */
    char* olean = with_ext(sc->stem, ".olean");
    struct stat st;
    if (lstat(olean, &st) == 0 && !S_ISREG(st.st_mode)) bad = str_copy("MathildaCert.olean");
    free(olean);
    return bad;
#else
    (void)sc;
    return NULL;
#endif
}

/* Whole file as a NUL-terminated string ("" when unreadable).  Caller frees. */
static char* slurp(const char* path) {
    RBuf b = {0};
    FILE* f = fopen(path, "rb");
    if (f) {
        char chunk[8192];
        size_t k;
        while ((k = fread(chunk, 1, sizeof chunk, f)) > 0) rb_putn(&b, chunk, k);
        fclose(f);
    }
    return rb_take(&b);
}

/* Write a new file: POSIX creates it exclusively (O_EXCL, and O_NOFOLLOW
 * where the platform has it), mode 0600, so an existing file or a planted
 * symlink makes the write fail instead of being followed. */
static bool write_file(const char* path, const char* text) {
#ifndef _WIN32
    int flags = O_WRONLY | O_CREAT | O_EXCL;
#ifdef O_NOFOLLOW
    flags |= O_NOFOLLOW;
#endif
    int fd = open(path, flags, 0600);
    if (fd < 0) return false;
    FILE* f = fdopen(fd, "wb");
    if (!f) { close(fd); return false; }
#else
    FILE* f = fopen(path, "wb");
    if (!f) return false;
#endif
    size_t n = strlen(text);
    bool ok = fwrite(text, 1, n, f) == n;
    if (fclose(f) != 0) ok = false;
    return ok;
}

/* ------------------------------------------------------------------ */
/* Command construction and execution                                  */
/* ------------------------------------------------------------------ */

/* Append one shell-quoted word; false when it cannot be quoted safely. */
static bool quote_word(RBuf* b, const char* w) {
#ifdef _WIN32
    if (strpbrk(w, "\"%!^&|<>")) return false;
    rb_puts(b, "\"");
    rb_puts(b, w);
    rb_puts(b, "\"");
#else
    rb_puts(b, "'");
    for (const char* p = w; *p; p++) {
        if (*p == '\'') rb_puts(b, "'\\''");
        else rb_putn(b, p, 1);
    }
    rb_puts(b, "'");
#endif
    return true;
}

/* The full command line, or NULL (with *bad the refused word). */
static char* build_command(const char* const* argv, int argc, const char* cwd,
                           const char* out, const char* err, const char** bad) {
    RBuf b = {0};
    bool ok = true;
#ifndef _WIN32
    rb_puts(&b, "set -C; ");
#endif
    if (cwd) {
#ifdef _WIN32
        rb_puts(&b, "cd /d ");
#else
        rb_puts(&b, "cd ");
#endif
        if (!quote_word(&b, cwd)) { ok = false; *bad = cwd; }
        rb_puts(&b, " && ");
    }
    for (int i = 0; i < argc && ok; i++) {
        if (i) rb_puts(&b, " ");
        if (!quote_word(&b, argv[i])) { ok = false; *bad = argv[i]; }
    }
    rb_puts(&b, " > ");
    if (ok && !quote_word(&b, out)) { ok = false; *bad = out; }
    rb_puts(&b, " 2> ");
    if (ok && !quote_word(&b, err)) { ok = false; *bad = err; }
    if (!ok) { free(b.s); return NULL; }
    return rb_take(&b);
}

/* The shell-quoted command as it is shown to the user ("Command"). */
static char* display_command(const char* const* argv, int argc, const char* cwd) {
    RBuf b = {0};
    if (cwd) { rb_puts(&b, "cd "); quote_word(&b, cwd); rb_puts(&b, " && "); }
    for (int i = 0; i < argc; i++) {
        if (i) rb_puts(&b, " ");
        quote_word(&b, argv[i]);
    }
    return rb_take(&b);
}

static int decode_status(int st) {
#ifdef _WIN32
    return st;
#else
    if (st == -1) return -1;
    if (WIFEXITED(st)) return WEXITSTATUS(st);
    return 128 + (WIFSIGNALED(st) ? WTERMSIG(st) : 0);
#endif
}

/* Run argv with stdout / stderr redirected to <stem>.out / <stem>.err. */
static int run_in(const char* stem, const char* const* argv, int argc, const char* cwd,
                  LeanProcResult* r, const char** why) {
    memset(r, 0, sizeof *r);
    *why = NULL;
    char* out = with_ext(stem, ".out");
    char* err = with_ext(stem, ".err");
    const char* bad = NULL;
    char* cmd = build_command(argv, argc, cwd, out, err, &bad);
    if (!cmd) {
        mth_message("LeanVerify", "badpath",
                    "path %s cannot be passed safely to the shell on this platform.", bad);
        free(out); free(err);
        *why = "unquotable argument";
        return -1;
    }
    fflush(NULL);
    double t0 = mth_mono_seconds();
    int st = system(cmd);
    r->seconds = mth_mono_seconds() - t0;
    free(cmd);
    r->exit_code = decode_status(st);
    /* The redirection creates <stem>.out; in the private directory nobody
     * else can, so a missing file means the shell never ran the command. */
    bool created = file_exists(out);
    r->out = slurp(out);
    r->err = slurp(err);
    remove(out);
    remove(err);
    free(out);
    free(err);
    if (st == -1) { *why = "system() failed"; lean_proc_result_free(r); return -1; }
    if (!created) {
        *why = "the shell did not create the output file";
        lean_proc_result_free(r);
        return -1;
    }
    return 0;
}

int lean_run_process(const char* const* argv, int argc, const char* cwd, LeanProcResult* r,
                     const char** why) {
    Scratch sc;
    if (!scratch_open(&sc)) {
        memset(r, 0, sizeof *r);
        *why = "cannot create a private temp directory";
        return -1;
    }
    int rc = run_in(sc.stem, argv, argc, cwd, r, why);
    scratch_close(&sc, false);
    return rc;
}

void lean_proc_result_free(LeanProcResult* r) {
    free(r->out);
    free(r->err);
    r->out = r->err = NULL;
}

/* ------------------------------------------------------------------ */
/* NDJSON message parsing                                              */
/* ------------------------------------------------------------------ */

static const char* js_ws(const char* p) {
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
    return p;
}

static void utf8_put(RBuf* b, unsigned long cp) {
    char t[4];
    size_t n;
    if (cp < 0x80) { t[0] = (char)cp; n = 1; }
    else if (cp < 0x800) { t[0] = (char)(0xC0 | (cp >> 6)); t[1] = (char)(0x80 | (cp & 0x3F)); n = 2; }
    else if (cp < 0x10000) {
        t[0] = (char)(0xE0 | (cp >> 12)); t[1] = (char)(0x80 | ((cp >> 6) & 0x3F));
        t[2] = (char)(0x80 | (cp & 0x3F)); n = 3;
    } else {
        t[0] = (char)(0xF0 | (cp >> 18)); t[1] = (char)(0x80 | ((cp >> 12) & 0x3F));
        t[2] = (char)(0x80 | ((cp >> 6) & 0x3F)); t[3] = (char)(0x80 | (cp & 0x3F)); n = 4;
    }
    rb_putn(b, t, n);
}

static int hexval(char ch) {
    if (ch >= '0' && ch <= '9') return ch - '0';
    if (ch >= 'a' && ch <= 'f') return ch - 'a' + 10;
    if (ch >= 'A' && ch <= 'F') return ch - 'A' + 10;
    return -1;
}

static bool read_hex4(const char* p, unsigned long* v) {
    unsigned long x = 0;
    for (int i = 0; i < 4; i++) {
        int h = hexval(p[i]);
        if (h < 0) return false;
        x = x * 16 + (unsigned long)h;
    }
    *v = x;
    return true;
}

/* A JSON string at p (which points at '"'), decoded into out (may be NULL to
 * skip).  Returns the position after the closing quote, or NULL. */
static const char* js_string(const char* p, RBuf* out) {
    if (*p != '"') return NULL;
    p++;
    while (*p && *p != '"') {
        if (*p != '\\') { if (out) rb_putn(out, p, 1); p++; continue; }
        p++;
        char ch = *p++;
        const char* lit = NULL;
        switch (ch) {
        case '"': lit = "\""; break;
        case '\\': lit = "\\"; break;
        case '/': lit = "/"; break;
        case 'b': lit = "\b"; break;
        case 'f': lit = "\f"; break;
        case 'n': lit = "\n"; break;
        case 'r': lit = "\r"; break;
        case 't': lit = "\t"; break;
        case 'u': {
            unsigned long cp;
            if (!read_hex4(p, &cp)) return NULL;
            p += 4;
            if (cp >= 0xD800 && cp < 0xDC00 && p[0] == '\\' && p[1] == 'u') {
                unsigned long lo;
                if (read_hex4(p + 2, &lo) && lo >= 0xDC00 && lo < 0xE000) {
                    cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00);
                    p += 6;
                }
            }
            if (out) utf8_put(out, cp);
            continue;
        }
        default: return NULL;
        }
        if (out) rb_puts(out, lit);
    }
    return *p == '"' ? p + 1 : NULL;
}

static const char* js_skip(const char* p);

/* Skip an object or array starting at p. */
static const char* js_skip_compound(const char* p, char open, char close) {
    if (*p != open) return NULL;
    p = js_ws(p + 1);
    if (*p == close) return p + 1;
    for (;;) {
        if (open == '{') {
            p = js_string(js_ws(p), NULL);
            if (!p) return NULL;
            p = js_ws(p);
            if (*p != ':') return NULL;
            p++;
        }
        p = js_skip(js_ws(p));
        if (!p) return NULL;
        p = js_ws(p);
        if (*p == ',') { p++; continue; }
        if (*p == close) return p + 1;
        return NULL;
    }
}

static const char* js_skip(const char* p) {
    p = js_ws(p);
    if (*p == '"') return js_string(p, NULL);
    if (*p == '{') return js_skip_compound(p, '{', '}');
    if (*p == '[') return js_skip_compound(p, '[', ']');
    const char* q = p;
    while (*q && *q != ',' && *q != '}' && *q != ']' && *q != ' ' && *q != '\n') q++;
    return q > p ? q : NULL;
}

typedef struct {
    RBuf severity, kind, data;
    long line;
} LeanMsg;

/* The "line" member of a {"line": n, "column": m} object. */
static const char* js_pos_line(const char* p, long* line) {
    if (*p != '{') return js_skip(p);
    p = js_ws(p + 1);
    if (*p == '}') return p + 1;
    for (;;) {
        RBuf key = {0};
        p = js_string(js_ws(p), &key);
        if (!p) { free(key.s); return NULL; }
        p = js_ws(p);
        if (*p != ':') { free(key.s); return NULL; }
        p = js_ws(p + 1);
        if (key.s && strcmp(key.s, "line") == 0) *line = strtol(p, NULL, 10);
        free(key.s);
        p = js_skip(p);
        if (!p) return NULL;
        p = js_ws(p);
        if (*p == ',') { p++; continue; }
        if (*p == '}') return p + 1;
        return NULL;
    }
}

/* One NDJSON line -> m.  False when the line is not a JSON object. */
static bool parse_message_line(const char* p, LeanMsg* m) {
    memset(m, 0, sizeof *m);
    p = js_ws(p);
    if (*p != '{') return false;
    p = js_ws(p + 1);
    if (*p == '}') return true;
    for (;;) {
        RBuf key = {0};
        p = js_string(js_ws(p), &key);
        if (!p) { free(key.s); return false; }
        p = js_ws(p);
        if (*p != ':') { free(key.s); return false; }
        p = js_ws(p + 1);
        const char* k = key.s ? key.s : "";
        RBuf* target = NULL;
        if (strcmp(k, "severity") == 0) target = &m->severity;
        else if (strcmp(k, "kind") == 0) target = &m->kind;
        else if (strcmp(k, "data") == 0) target = &m->data;
        if (target && *p == '"') {
            free(target->s);
            memset(target, 0, sizeof *target);
            p = js_string(p, target);
        } else if (strcmp(k, "pos") == 0) {
            p = js_pos_line(p, &m->line);
        } else {
            p = js_skip(p);
        }
        free(key.s);
        if (!p) return false;
        p = js_ws(p);
        if (*p == ',') { p++; continue; }
        return *p == '}';
    }
}

static void lean_msg_free(LeanMsg* m) {
    free(m->severity.s);
    free(m->kind.s);
    free(m->data.s);
}

/* ------------------------------------------------------------------ */
/* Expr helpers                                                        */
/* ------------------------------------------------------------------ */

static Expr* mk_rule(const char* key, Expr* v) {
    Expr* a[2] = { expr_new_string(key), v };
    return expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
}

/* An Association from n Rule nodes (consumed). */
static Expr* mk_assoc(Expr** rules, size_t n) {
    Expr* r = assoc_from_rules(rules, n);
    for (size_t i = 0; i < n; i++) expr_free(rules[i]);
    return r;
}

static Expr* mk_bool(bool b) { return expr_new_symbol(b ? SYM_True : SYM_False); }

typedef struct { Expr** v; size_t n, cap; } EList;

static void el_push(EList* l, Expr* e) {
    if (l->n == l->cap) {
        size_t nc = l->cap ? l->cap * 2 : 8;
        Expr** p = realloc(l->v, nc * sizeof *p);
        if (!p) abort();
        l->v = p;
        l->cap = nc;
    }
    l->v[l->n++] = e;
}

static Expr* el_list(EList* l) {
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), l->v, l->n);
    free(l->v);
    l->v = NULL;
    l->n = l->cap = 0;
    return r;
}

static const char* const STANDARD_AXIOMS[] = { "propext", "Classical.choice", "Quot.sound", NULL };

static bool standard_axiom(const char* a) {
    for (size_t i = 0; STANDARD_AXIOMS[i]; i++)
        if (strcmp(a, STANDARD_AXIOMS[i]) == 0) return true;
    return false;
}

/* A `#print axioms` report: "'thm' depends on axioms: [a, b]" or
 * "'thm' does not depend on any axioms".  On a match, *name is the theorem
 * (fresh), the axioms are pushed onto ax, and *all_std says whether every one
 * is standard. */
static bool parse_axioms(const char* text, char** name, EList* ax, bool* all_std) {
    if (text[0] != '\'') return false;
    const char* deps = strstr(text, "' depends on axioms: [");
    const char* none = deps ? NULL : strstr(text, "' does not depend on any axioms");
    const char* end = deps ? deps : none;
    if (!end) return false;
    size_t nl = (size_t)(end - text - 1);
    char* nm = malloc(nl + 1);
    if (!nm) abort();
    memcpy(nm, text + 1, nl);
    nm[nl] = '\0';
    *name = nm;
    *all_std = true;
    if (end == none) return true;
    const char* p = end + strlen("' depends on axioms: [");
    while (*p && *p != ']') {
        while (*p == ' ' || *p == ',') p++;
        const char* q = p;
        while (*q && *q != ',' && *q != ']') q++;
        const char* e2 = q;
        while (e2 > p && e2[-1] == ' ') e2--;
        if (e2 > p) {
            Expr* s = expr_new_string_n(p, (size_t)(e2 - p));
            if (!standard_axiom(s->data.string)) *all_std = false;
            el_push(ax, s);
        }
        p = q;
    }
    return true;
}

/* ------------------------------------------------------------------ */
/* Builtins                                                            */
/* ------------------------------------------------------------------ */

static bool is_sym(const Expr* e, const char* s) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == s;
}

static bool is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/* Lean`Internal`RunProcess[{exe, args...}, cwd | None] */
static Expr* builtin_lean_runprocess(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* l = res->data.function.args[0];
    Expr* cw = res->data.function.args[1];
    if (!is_list(l) || l->data.function.arg_count == 0) return NULL;
    if (!(cw->type == EXPR_STRING || is_sym(cw, SYM_None))) return NULL;
    size_t n = l->data.function.arg_count;
    const char** argv = malloc(n * sizeof *argv);
    if (!argv) abort();
    for (size_t i = 0; i < n; i++) {
        if (l->data.function.args[i]->type != EXPR_STRING) { free((void*)argv); return NULL; }
        argv[i] = l->data.function.args[i]->data.string;
    }
    LeanProcResult r;
    const char* why = NULL;
    int rc = lean_run_process(argv, (int)n, cw->type == EXPR_STRING ? cw->data.string : NULL,
                              &r, &why);
    free((void*)argv);
    if (rc != 0) {
        mth_message("LeanVerify", "fail", "Lean process could not be started: %s.", why);
        return expr_new_symbol(SYM_DollarFailed);
    }
    Expr* rules[3] = {
        mk_rule("ExitCode", expr_new_integer(r.exit_code)),
        mk_rule("StandardOutput", expr_new_string(r.out)),
        mk_rule("StandardError", expr_new_string(r.err)),
    };
    lean_proc_result_free(&r);
    return mk_assoc(rules, 3);
}

/* The lake executable beside `exe` (bare "lake" when exe has no directory). */
static char* lake_beside(const char* exe) {
    const char* slash = strrchr(exe, '/');
#ifdef _WIN32
    const char* bs = strrchr(exe, '\\');
    if (!slash || (bs && bs > slash)) slash = bs;
#endif
    if (!slash) return str_copy("lake");
    char* dir = str_copy(exe);
    dir[slash - exe + 1] = '\0';      /* keep the directory and its separator */
    RBuf b = {0};
    rb_puts(&b, dir);
    rb_puts(&b, "lake");
    free(dir);
    return rb_take(&b);
}

/* Kernel replay (src/internal/lean/Replay.lean) of the compiled module.
 *
 * The checker prints "C <name>" per declaration of the module, followed by
 * "A <axiom>" per axiom it depends on, between MATHILDA-REPLAY-BEGIN and
 * MATHILDA-REPLAY-OK; MATHILDA-REPLAY-FAIL <why> on a kernel rejection. */
typedef struct {
    bool ran;           /* the checker was started */
    bool ok;            /* replay accepted every declaration */
    bool all_std;       /* every declaration depends on standard axioms only */
    Expr* axioms;       /* Association name -> {axioms...} (owned) */
    char* why;          /* a failure reason, or NULL (owned) */
    double seconds;
} ReplayResult;

/* The marker the checker reports for a declaration it could not replay
 * (unsafe / partial implementation detail). */
#define NOT_REPLAYED "mathilda.notReplayed."

static void replay_parse(const char* out, ReplayResult* rr) {
    EList rules = {0}, ax = {0};
    char* cur = NULL;          /* the current declaration */
    bool begun = false, ok = false;
    rr->all_std = true;
    const char* p = out;
    while (p && *p) {
        const char* nl = strchr(p, '\n');
        size_t n = nl ? (size_t)(nl - p) : strlen(p);
        if (n == strlen("MATHILDA-REPLAY-BEGIN") && strncmp(p, "MATHILDA-REPLAY-BEGIN", n) == 0) {
            begun = true;
        } else if (begun && n == strlen("MATHILDA-REPLAY-OK") &&
                   strncmp(p, "MATHILDA-REPLAY-OK", n) == 0) {
            ok = true;
        } else if (strncmp(p, "MATHILDA-REPLAY-FAIL", strlen("MATHILDA-REPLAY-FAIL")) == 0 &&
                   !rr->why) {
            const char* w = p + strlen("MATHILDA-REPLAY-FAIL");
            while (*w == ' ' && w < p + n) w++;
            rr->why = malloc((size_t)(p + n - w) + 1);
            if (!rr->why) abort();
            memcpy(rr->why, w, (size_t)(p + n - w));
            rr->why[p + n - w] = '\0';
        } else if (begun && !ok && n > 2 && (p[0] == 'C' || p[0] == 'A') && p[1] == ' ') {
            Expr* s = expr_new_string_n(p + 2, n - 2);
            if (p[0] == 'C') {
                if (cur) { el_push(&rules, mk_rule(cur, el_list(&ax))); free(cur); }
                cur = str_copy(s->data.string);
                expr_free(s);
            } else if (cur) {
                const char* a = s->data.string;
                /* A declaration's own not-replayed marker is no axiom: such a
                 * declaration is an unsafe implementation detail, and any
                 * safe declaration that uses it carries the marker too. */
                bool self = strncmp(a, NOT_REPLAYED, strlen(NOT_REPLAYED)) == 0 &&
                            strcmp(a + strlen(NOT_REPLAYED), cur) == 0;
                if (!standard_axiom(a) && !self) rr->all_std = false;
                el_push(&ax, s);
            } else {
                expr_free(s);
            }
        }
        p = nl ? nl + 1 : NULL;
    }
    if (cur) { el_push(&rules, mk_rule(cur, el_list(&ax))); free(cur); }
    free(ax.v);
    rr->ok = ok && !rr->why;
    if (!rr->ok && !rr->why) rr->why = str_copy("the replay checker gave no verdict");
    rr->axioms = mk_assoc(rules.v, rules.n);
    free(rules.v);
}

/* Replay <stem>.olean with `exe --run script` (through `lake env` in the
 * project); the checker's output goes to <dir>/check.out / .err. */
static void replay_module(const Scratch* sc, const char* exe, const char* lake,
                          const char* project, const char* script, ReplayResult* rr) {
    memset(rr, 0, sizeof *rr);
    char* olean = with_ext(sc->stem, ".olean");
    if (!file_exists(olean)) {
        rr->why = str_copy("Lean wrote no compiled module (the source has errors)");
        rr->axioms = mk_assoc(NULL, 0);
        free(olean);
        return;
    }
    char* odd = unexpected_entry(sc);
    if (odd) {
        RBuf w = {0};
        rb_puts(&w, "unexpected file in the private directory after elaboration: ");
        rb_puts(&w, odd);
        free(odd);
        rr->why = rb_take(&w);
        rr->axioms = mk_assoc(NULL, 0);
        free(olean);
        return;
    }
    RBuf cs = {0};
    if (sc->dir) { rb_puts(&cs, sc->dir); rb_puts(&cs, LEAN_SEP "check"); }
    else { rb_puts(&cs, sc->stem); rb_puts(&cs, "-check"); }
    const char* argv[6];
    int argc = 0;
    if (project) { argv[argc++] = lake; argv[argc++] = "env"; }
    argv[argc++] = exe;
    argv[argc++] = "--run";
    argv[argc++] = script;
    argv[argc++] = olean;
    LeanProcResult r;
    const char* why = NULL;
    rr->ran = run_in(cs.s, argv, argc, project, &r, &why) == 0;
    free(cs.s);
    free(olean);
    if (!rr->ran) {
        rr->why = str_copy(why ? why : "the replay checker could not be started");
        rr->axioms = mk_assoc(NULL, 0);
        return;
    }
    rr->seconds = r.seconds;
    replay_parse(r.out, rr);
    if (r.exit_code != 0) rr->ok = false;
    if (!rr->ok && !rr->why) rr->why = str_copy("the replay checker failed");
    lean_proc_result_free(&r);
}

/* Lean`Internal`VerifySource[src, exe, project | None, keep, replay | None]
 *
 * replay: the path of src/internal/lean/Replay.lean.  The source is then also
 * compiled to <dir>/MathildaCert.olean (`-R <dir> -o ...`), and after a clean
 * Lean run the checker replays it in a fresh `lean` process; "Success" then
 * also needs the replay to pass and every declaration of the module to depend
 * on the standard axioms only (as computed by the checker, not by the
 * checked file's own output). */
static Expr* builtin_lean_verifysource(Expr* res) {
    size_t nargs = res->data.function.arg_count;
    if (nargs != 4 && nargs != 5) return NULL;
    Expr* const* a = res->data.function.args;
    if (a[0]->type != EXPR_STRING || a[1]->type != EXPR_STRING) return NULL;
    if (!(a[2]->type == EXPR_STRING || is_sym(a[2], SYM_None))) return NULL;
    if (!is_sym(a[3], SYM_True) && !is_sym(a[3], SYM_False)) return NULL;
    if (nargs == 5 && !(a[4]->type == EXPR_STRING || is_sym(a[4], SYM_None))) return NULL;
    bool keep = is_sym(a[3], SYM_True);
    const char* exe = a[1]->data.string;
    const char* project = a[2]->type == EXPR_STRING ? a[2]->data.string : NULL;
    const char* script = (nargs == 5 && a[4]->type == EXPR_STRING) ? a[4]->data.string : NULL;

    Scratch sc;
    if (!scratch_open(&sc)) {
        mth_message("LeanVerify", "fail", "Lean process could not be started: %s.",
                    "cannot create a private temp directory");
        return expr_new_symbol(SYM_DollarFailed);
    }
    char* src = with_ext(sc.stem, ".lean");
    if (!write_file(src, a[0]->data.string)) {
        mth_message("LeanVerify", "fail", "Lean process could not be started: %s.",
                    "cannot write the temp source file");
        free(src);
        scratch_close(&sc, false);
        return expr_new_symbol(SYM_DollarFailed);
    }

    char* lake = project ? lake_beside(exe) : NULL;
    char* olean = with_ext(sc.stem, ".olean");
    const char* argv[9];
    int argc = 0;
    if (project) { argv[argc++] = lake; argv[argc++] = "env"; }
    argv[argc++] = exe;
    argv[argc++] = "--json";
    if (script) {
        if (sc.dir) { argv[argc++] = "-R"; argv[argc++] = sc.dir; }
        argv[argc++] = "-o";
        argv[argc++] = olean;
    }
    argv[argc++] = src;

    LeanProcResult r;
    const char* why = NULL;
    int rc = run_in(sc.stem, argv, argc, project, &r, &why);
    char* shown = display_command(argv, argc, project);
    ReplayResult rr;
    memset(&rr, 0, sizeof rr);
    if (rc == 0 && script && r.exit_code == 0)
        replay_module(&sc, exe, lake, project, script, &rr);
    else if (script)
        rr.axioms = mk_assoc(NULL, 0);
    free(lake);
    free(olean);
    scratch_close(&sc, keep && rc == 0);
    if (rc != 0) {
        mth_message("LeanVerify", "fail", "Lean process could not be started: %s.", why);
        if (rr.axioms) expr_free(rr.axioms);
        free(rr.why);
        free(shown);
        free(src);
        return expr_new_symbol(SYM_DollarFailed);
    }
    EList msgs = {0}, axr = {0}, reps = {0};
    bool any_error = false, any_sorry = false, axioms_ok = true;
    for (char* line = r.out; line && *line; ) {
        char* nl = strchr(line, '\n');
        if (nl) *nl = '\0';
        LeanMsg m;
        if (parse_message_line(line, &m)) {
            const char* sev = m.severity.s ? m.severity.s : "";
            const char* kind = m.kind.s ? m.kind.s : "";
            const char* text = m.data.s ? m.data.s : "";
            if (strcmp(sev, "error") == 0) any_error = true;
            if (strcmp(kind, "hasSorry") == 0) any_sorry = true;
            char* thm = NULL;
            EList ax = {0};
            bool std_only = true;
            if (strcmp(sev, "information") == 0 && parse_axioms(text, &thm, &ax, &std_only)) {
                if (!std_only) axioms_ok = false;
                Expr* axl = el_list(&ax);
                Expr* rep[3] = {
                    mk_rule("Theorem", expr_new_string(thm)),
                    mk_rule("Line", expr_new_integer(m.line)),
                    mk_rule("Axioms", expr_copy(axl)),
                };
                el_push(&reps, mk_assoc(rep, 3));
                el_push(&axr, mk_rule(thm, axl));
                free(thm);
            } else {
                free(ax.v);
            }
            Expr* fields[4] = {
                mk_rule("Severity", expr_new_string(sev)),
                mk_rule("Kind", expr_new_string(kind)),
                mk_rule("Line", expr_new_integer(m.line)),
                mk_rule("Text", expr_new_string(text)),
            };
            el_push(&msgs, mk_assoc(fields, 4));
        }
        lean_msg_free(&m);
        line = nl ? nl + 1 : NULL;
    }
    /* A process that wrote no message but failed (missing executable, a
     * crash): report its stderr as an error message. */
    if (r.exit_code != 0 && msgs.n == 0) {
        any_error = true;
        Expr* fields[4] = {
            mk_rule("Severity", expr_new_string("error")),
            mk_rule("Kind", expr_new_string("process")),
            mk_rule("Line", expr_new_integer(0)),
            mk_rule("Text", expr_new_string(*r.err ? r.err : r.out)),
        };
        el_push(&msgs, mk_assoc(fields, 4));
    }
    bool success = r.exit_code == 0 && !any_error && !any_sorry && axioms_ok;
    if (script) success = success && rr.ok && rr.all_std;

    size_t nax = axr.n;
    Expr* axioms = mk_assoc(axr.v, nax);
    free(axr.v);
    Expr* rules[12];
    size_t k = 0;
    rules[k++] = mk_rule("Success", mk_bool(success));
    rules[k++] = mk_rule("ExitCode", expr_new_integer(r.exit_code));
    rules[k++] = mk_rule("Messages", el_list(&msgs));
    rules[k++] = mk_rule("Axioms", axioms);
    rules[k++] = mk_rule("AxiomReports", el_list(&reps));
    rules[k++] = mk_rule("Timing", expr_new_real(r.seconds));
    rules[k++] = mk_rule("Command", expr_new_string(shown));
    if (keep) rules[k++] = mk_rule("SourceFile", expr_new_string(src));
    if (script) {
        /* "Replay": True, or the reason it failed; "ReplayAxioms": per
         * declaration of the module, the axioms the checker found. */
        rules[k++] = mk_rule("Replay", rr.ok ? mk_bool(true)
                                     : expr_new_string(rr.why ? rr.why : "not run"));
        rules[k++] = mk_rule("ReplayAxioms", rr.axioms);
        rules[k++] = mk_rule("ReplayTiming", expr_new_real(rr.seconds));
        rr.axioms = NULL;
    }
    free(rr.why);
    free(shown);
    free(src);
    lean_proc_result_free(&r);
    return mk_assoc(rules, k);
}

void lean_init(void) {
    lean_check_init();
    symtab_add_builtin("Lean`Internal`ToLean", builtin_lean_tolean);
    symtab_set_docstring("Lean`Internal`ToLean",
        "Lean`Internal`ToLean[expr, type, coreOnly, {heads...}, closed]\n"
        "\tgives the strict Lean 4 translation of expr as a String, with type \"Int\",\n"
        "\t\"Rat\", \"Real\" or None for rational literals and binders; coreOnly -> True\n"
        "\trefuses Mathlib-only heads; heads are the Lean names allowed as applied\n"
        "\tfunctions; closed -> True (optional, default False) refuses any symbol not\n"
        "\tbound by a quantifier. An untranslatable node gives $Failed and LeanForm::untrans.");
    symtab_get_def("Lean`Internal`ToLean")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("Lean`Internal`RunProcess", builtin_lean_runprocess);
    symtab_set_docstring("Lean`Internal`RunProcess",
        "Lean`Internal`RunProcess[{exe, args...}, cwd]\n"
        "\truns exe with the given arguments (each quoted for the shell) in directory\n"
        "\tcwd (None: the current one) and gives <|\"ExitCode\", \"StandardOutput\",\n"
        "\t\"StandardError\"|>. Internal to the Lean integration.");
    symtab_get_def("Lean`Internal`RunProcess")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("Lean`Internal`VerifySource", builtin_lean_verifysource);
    symtab_set_docstring("Lean`Internal`VerifySource",
        "Lean`Internal`VerifySource[src, exe, project, keep, replay]\n"
        "\twrites the Lean 4 source src to a temp file, checks it with exe --json (in\n"
        "\tthe Lake project directory project via lake env, unless project is None) and\n"
        "\tgives <|\"Success\", \"ExitCode\", \"Messages\", \"Axioms\", \"Timing\",\n"
        "\t\"Command\"|>. keep -> True keeps the source file (\"SourceFile\").\n"
        "\treplay (optional): the kernel-replay script; the module is compiled and\n"
        "\treplayed in a fresh lean process (\"Replay\", \"ReplayAxioms\").");
    symtab_get_def("Lean`Internal`VerifySource")->attributes |= ATTR_PROTECTED;
}
