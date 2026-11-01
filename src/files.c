/*
 * files.c — Filesystem-predicate and path-string builtins.
 *
 * Builtins:
 *   FileExistsQ["name"]    — touches the filesystem (lstat).
 *   FileSize["name"]       — touches the filesystem (stat); byte count.
 *   FileExtension["name"]  — pure string manipulation.
 *   FileBaseName["name"]   — pure string manipulation.
 *   FilePrint["name", ...] — streams file contents to stdout.
 *   FileNameJoin[{...}]    — assembles a file name from parts.
 *   FileNameSplit["name"]  — splits a file name into its parts.
 *   DirectoryName["name"]  — directory portion of a path (pure string).
 *   FileNames[forms, dirs, n] — files matching forms, searched to depth n.
 *   Directory[], SetDirectory["dir"], ResetDirectory[], DirectoryStack[]
 *                          — the working directory and its push/pop stack.
 *   $HomeDirectory         — the user's home directory (from the environment).
 *   $TemporaryDirectory    — the host temporary directory (TMPDIR / TEMP / TMP).
 *   Environment["var"]     — an environment variable's value, or $Failed.
 *   ParentDirectory, FileNameTake, AbsoluteFileName, FileType, DirectoryQ, ToFileName,
 *   DeleteFile             — path manipulation and deletion.
 *   $PathnameSeparator     — the host pathname separator ("/" or "\\").
 *
 * Also initialises $InputFileName to "" (loadmodule.c swaps it in and out
 * around every file read).
 *
 * Both string builtins follow the POSIX convention: the path
 * separator is `/`, the directory portion is everything up to and
 * including the final `/`, and the file-name component is everything
 * after.  The "extension" is the suffix following the LAST `.` in
 * the file-name component, but a leading `.` (e.g. ".bashrc") is
 * treated as part of the base name rather than an extension marker.
 *
 * No filesystem access happens for FileExtension / FileBaseName; they
 * operate purely on the given string.
 */

/* lstat is POSIX, not C99.  glibc hides it under -std=c99 unless we
 * declare the feature-test macro before any include.  Darwin exposes
 * it implicitly but accepts the macro.  realpath is X/Open: glibc's
 * <stdlib.h> declares it only under __USE_XOPEN_EXTENDED or __USE_MISC,
 * which _POSIX_C_SOURCE alone never sets (issue #110; Ubuntu's gcc hid
 * this because its default _FORTIFY_SOURCE declares realpath unguarded). */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include "files.h"
#include "expr.h"
#include "symtab.h"
#include "sym_names.h"
#include "attr.h"
#include "common.h"
#include "message.h"
#include "print.h"
#include "regex_common.h"

#include <dirent.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>     /* _getcwd / _chdir */
#define getcwd _getcwd
#define chdir  _chdir
#else
#include <unistd.h>     /* getcwd / chdir (POSIX, see _POSIX_C_SOURCE above) */
#endif
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Pathname separator for the host operating system.  FileNameJoin uses
 * this as its default when no OperatingSystem->"..." option is given. */
#ifdef _WIN32
#define HOST_SEP '\\'
#else
#define HOST_SEP '/'
#endif

/* Native Windows (MinGW) has no lstat.  FileExistsQ uses lstat so a dangling
 * symlink still counts as "exists"; on Windows there is no POSIX symlink to
 * distinguish, so stat is the exact equivalent. */
#ifdef _WIN32
#define lstat stat
#endif

/* Emit a Wolfram-style diagnostic ("Head::tag: text") unless Quiet is active.
 * The firing is always noted so Check[] sees it. Mirrors dt_msg in
 * src/datetime.c. */
static void fs_msg(const char* fmt, ...) {
    /* Through the funnel by name, so Quiet[e, Head::tag] and
     * Check[e, f, Head::tag] see which message this is. */
    va_list ap;
    va_start(ap, fmt);
    mth_message_linev(0, fmt, ap);
    va_end(ap);
}

static char* home_from_env(void);

/* The path a file argument names: a String or File["s"], with a leading `~`
 * (alone or before a separator) expanded to the home directory.  Caller
 * frees; NULL when the argument is neither form. */
static char* file_arg_path(const Expr* a) {
    const char* p = NULL;
    if (a->type == EXPR_STRING) p = a->data.string;
    else if (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL &&
             a->data.function.head->data.symbol.name == SYM_File &&
             a->data.function.arg_count == 1 && a->data.function.args[0]->type == EXPR_STRING)
        p = a->data.function.args[0]->data.string;
    if (!p) return NULL;
    if (p[0] == '~' && (p[1] == '\0' || p[1] == '/' || p[1] == '\\')) {
        char* home = home_from_env();
        if (home) {
            size_t hl = strlen(home), rl = strlen(p + 1);
            char* r = malloc(hl + rl + 1);
            if (r) { memcpy(r, home, hl); memcpy(r + hl, p + 1, rl + 1); }
            free(home);
            return r;
        }
    }
    char* r = malloc(strlen(p) + 1);
    if (r) strcpy(r, p);
    return r;
}

/* FileExistsQ["name"]
 *
 * Returns True if anything exists at the given path, False otherwise.
 * Uses lstat() rather than stat() so dangling symlinks (which are
 * themselves filesystem objects) are still reported as existing.
 *
 * Anything other than a single EXPR_STRING argument leaves the call
 * unevaluated (NULL return).  The evaluator owns `res`: on a
 * successful rewrite it frees the input itself after this builtin
 * returns, so we must NOT free it here. */
Expr* builtin_fileexistsq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1) return NULL;

    char* path = file_arg_path(res->data.function.args[0]);   /* String / File, `~` */
    if (!path) return NULL;

    struct stat st;
    int exists = (lstat(path, &st) == 0);
    free(path);
    return expr_new_symbol(exists ? "True" : "False");
}

/* FileSize["name"]
 *
 * Returns the size of the file at the given path, in bytes, as an Integer
 * (not a Quantity).  Uses stat() rather than lstat() so a symbolic link is
 * followed to its target — the size reported is the target file's size.
 *
 * If the path cannot be stat'd (most commonly because no such file exists),
 * a FileSize::nffil message is emitted through the Quiet[]/Check[]-aware
 * helper and $Failed is returned.
 *
 * Anything other than a single EXPR_STRING argument leaves the call
 * unevaluated (NULL return), so symbolic arguments flow through unchanged.
 * The evaluator owns `res`: return a fresh Expr* on rewrite and it frees the
 * input itself; we must NOT free it here. */
Expr* builtin_filesize(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1) return NULL;

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING) return NULL;

    struct stat st;
    if (stat(arg->data.string, &st) != 0) {
        fs_msg("FileSize::nffil: Cannot find file %s.", arg->data.string);
        return expr_new_symbol(SYM_DollarFailed);
    }
    return expr_new_integer((int64_t)st.st_size);
}

/* Split "path" into (directory-prefix-length, filename-component).
 * The filename component is everything after the last '/'; it may be
 * empty (when "path" is empty or ends with '/').  Pure string view —
 * no allocation. */
static const char* filename_component(const char* path) {
    const char* slash = strrchr(path, '/');
    return slash ? slash + 1 : path;
}

/* Find the extension start within `name` (the filename component
 * only — pass the output of filename_component, not a raw path).
 *
 * Returns the offset of the character AFTER the last '.' that
 * qualifies as an extension separator, or strlen(name) if there is
 * no extension.  By the Mathematica contract, a '.' qualifies only
 * when it is neither the first nor the last character of `name`:
 *
 *   "file.txt"  -> 5  ("txt")
 *   "a.b.c"     -> 4  ("c")
 *   "file"      -> 4  ("")
 *   "file."     -> 5  ("")  trailing dot is not an extension marker
 *   ".bashrc"   -> 7  ("")  leading dot is not an extension marker
 *   ""          -> 0  ("")
 */
static size_t extension_offset(const char* name) {
    size_t len = strlen(name);
    if (len == 0) return 0;
    /* Walk from the end looking for a qualifying '.'.  Stop before
     * index 0 so a lone leading dot is never picked up. */
    for (size_t i = len; i > 1; i--) {
        if (name[i - 1] == '.') {
            /* A '.' at the very end of `name` means no extension. */
            if (i == len) return len;
            return i;
        }
    }
    return len;
}

/* FileExtension["name"]
 *
 * Returns the substring after the last '.' in the filename component
 * of "name", excluding the dot.  Returns "" when no extension is
 * present.  Pure string manipulation — does not touch the
 * filesystem. */
Expr* builtin_fileextension(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1) return NULL;

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING) return NULL;

    const char* name = filename_component(arg->data.string);
    size_t off = extension_offset(name);
    /* extension_offset returns len(name) when there is none; passing
     * name + len yields an empty string, which is the right answer. */
    return expr_new_string(name + off);
}

/* FileBaseName["name"]
 *
 * Returns the filename component of "name" with its trailing
 * extension removed.  When there is no extension the filename
 * component is returned verbatim.  Directory specifications are
 * always dropped.  Pure string manipulation — does not touch the
 * filesystem. */
Expr* builtin_filebasename(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1) return NULL;

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING) return NULL;

    const char* name = filename_component(arg->data.string);
    size_t off = extension_offset(name);
    size_t base_len;
    if (off == strlen(name)) {
        /* No extension — keep the whole filename component, including
         * any trailing '.' (matches "file." -> "file."). */
        base_len = off;
    } else {
        /* off points just past the '.'; exclude the '.' from the
         * base name. */
        base_len = off - 1;
    }

    /* Build a NUL-terminated copy of the prefix so expr_new_string
     * (which strdup's its input) sees exactly the bytes we want.
     * Reusing arg->data.string is not safe — the evaluator will free
     * res (and with it arg) after we return. */
    char* buf = (char*)malloc(base_len + 1);
    if (!buf) return NULL;  /* leave res untouched on OOM */
    memcpy(buf, name, base_len);
    buf[base_len] = '\0';
    Expr* result = expr_new_string(buf);
    free(buf);
    return result;
}

/* === FilePrint =========================================================
 *
 * Reads "filename" into a single buffer, walks the buffer once to record
 * the offset and length of each line (the trailing '\n' is included in
 * the line's length when present), then writes the selected lines back
 * out via fwrite so the bytes — including any embedded NULs or non-UTF-8
 * sequences — are passed through verbatim.
 *
 * Selection semantics:
 *
 *   No selector             →  all lines.
 *   Integer n > 0           →  first n lines (clamped to total).
 *   Integer n < 0           →  last |n| lines (clamped).
 *   Integer 0               →  no output (still a successful call).
 *   Span[m, n]              →  lines m through n, step +1.
 *   Span[m, n, s]           →  lines m through n, step s (s ≠ 0; s may
 *                              be negative to walk backwards).
 *
 * Within a Span, negative index `-k` means `total + k + 1` (so -1 is the
 * last line).  `All` in any slot means the natural endpoint for that
 * slot (1 for start, total for end, 1 for step).  Out-of-range indices
 * after normalisation cause the call to be unevaluated rather than
 * producing a partial print — matching Mathematica's behaviour of
 * signalling a Part::partw and leaving the call alone.
 */

typedef struct {
    const char* start;  /* points into `buf` */
    size_t      len;    /* includes trailing '\n' when the line has one */
} LineSpan;

/* Read the whole file into a freshly-malloced buffer.  Returns NULL on
 * any I/O failure (the caller distinguishes "couldn't open" from
 * "couldn't read" only insofar as both yield $Failed).  *out_size gets
 * the byte count on success. */
static char* fileprint_slurp(const char* filename, size_t* out_size) {
    FILE* fp = fopen(filename, "rb");
    if (!fp) return NULL;

    if (fseek(fp, 0, SEEK_END) != 0) { fclose(fp); return NULL; }
    long sz = ftell(fp);
    if (sz < 0)                      { fclose(fp); return NULL; }
    if (fseek(fp, 0, SEEK_SET) != 0) { fclose(fp); return NULL; }

    /* +1 so a zero-byte file still gets a valid allocation. */
    char* buf = (char*)malloc((size_t)sz + 1);
    if (!buf) { fclose(fp); return NULL; }

    size_t n = fread(buf, 1, (size_t)sz, fp);
    fclose(fp);
    if (n != (size_t)sz) { free(buf); return NULL; }
    buf[n] = '\0';
    *out_size = n;
    return buf;
}

/* Normalise line ends in place -- CRLF and a lone CR each become one LF --
 * and return the new length.  Import "Text"/"Lines" see the file this way;
 * "String" and ReadString keep the bytes as they are. */
size_t files_crlf_normalize(char* buf, size_t len) {
    size_t w = 0;
    for (size_t r = 0; r < len; r++) {
        if (buf[r] == '\r') {
            if (r + 1 < len && buf[r + 1] == '\n') continue;   /* CRLF: keep the LF */
            buf[w++] = '\n';                                   /* lone CR */
        } else {
            buf[w++] = buf[r];
        }
    }
    buf[w] = '\0';
    return w;
}

/* Whole file, bytes untouched; NULL if it cannot be read. */
char* files_slurp_raw(const char* path, size_t* len) {
    return fileprint_slurp(path, len);
}

/* Whole file as text with line ends normalised; NULL if it cannot be read. */
char* files_slurp_text(const char* path, size_t* len) {
    char* buf = fileprint_slurp(path, len);
    if (buf) *len = files_crlf_normalize(buf, *len);
    return buf;
}

/*
 * Import[path, fmt] for the text formats, fmt one of:
 *   "String"  the contents, byte for byte
 *   "Text"    the contents with line ends (CRLF, CR) normalised to "\n",
 *             less ONE trailing line break
 *   "Lines"   a list of lines, leading and trailing empty lines dropped
 * A file that cannot be read gives Import::nffil and $Failed.  Returns NULL
 * for any other fmt.
 */
Expr* files_import_text(const char* path, const char* fmt) {
    int kind = strcmp(fmt, "Text") == 0 ? 0 : strcmp(fmt, "String") == 0 ? 1
             : strcmp(fmt, "Lines") == 0 ? 2 : -1;
    if (kind < 0) return NULL;
    size_t len = 0;
    char* buf = kind == 1 ? files_slurp_raw(path, &len) : files_slurp_text(path, &len);
    if (!buf) {
        mth_message("Import", "nffil", "File %s not found during Import.", path);
        return expr_new_symbol(SYM_DollarFailed);
    }
    Expr* out;
    if (kind == 0) {
        if (len > 0 && buf[len - 1] == '\n') buf[--len] = '\0';
        out = expr_new_string(buf);
    } else if (kind == 1) {
        out = expr_new_string(buf);
    } else {
        while (len > 0 && buf[len - 1] == '\n') buf[--len] = '\0';
        size_t start = 0;                      /* leading empty lines go too */
        while (start < len && buf[start] == '\n') start++;
        size_t n = 0, cap = 16;
        Expr** lines = malloc(cap * sizeof(Expr*));
        for (size_t i = start; start < len && i <= len; i++) {
            if (i == len || buf[i] == '\n') {
                if (n == cap) { cap *= 2; lines = realloc(lines, cap * sizeof(Expr*)); }
                buf[i] = '\0';
                lines[n++] = expr_new_string(buf + start);
                start = i + 1;
            }
        }
        out = expr_new_function(expr_new_symbol(SYM_List), lines, n);
        free(lines);
    }
    free(buf);
    return out;
}

/* Walk `buf` once and record one LineSpan per line.  A line ends at the
 * next '\n' (inclusive) or at end-of-buffer.  An empty buffer yields
 * zero lines; a buffer that is exactly "\n" yields one (empty-with-LF)
 * line.  Caller frees *out_lines. */
static int fileprint_split_lines(const char* buf, size_t size,
                                 LineSpan** out_lines, size_t* out_count) {
    *out_lines = NULL;
    *out_count = 0;
    if (size == 0) return 1;

    /* First pass: count newlines so we can allocate exactly once. */
    size_t lines = 0;
    for (size_t i = 0; i < size; i++) {
        if (buf[i] == '\n') lines++;
    }
    if (buf[size - 1] != '\n') lines++;  /* trailing unterminated line */

    LineSpan* spans = (LineSpan*)malloc(sizeof(LineSpan) * lines);
    if (!spans) return 0;

    size_t idx = 0;
    size_t start = 0;
    for (size_t i = 0; i < size; i++) {
        if (buf[i] == '\n') {
            spans[idx].start = buf + start;
            spans[idx].len   = i - start + 1;  /* include the '\n' */
            idx++;
            start = i + 1;
        }
    }
    if (start < size) {
        spans[idx].start = buf + start;
        spans[idx].len   = size - start;
        idx++;
    }

    *out_lines = spans;
    *out_count = idx;
    return 1;
}

/* Resolve one Span endpoint.  `slot` is 0 for start, 1 for end, 2 for
 * step — used so `All` picks the right default.  Returns 1 on success
 * (writes *out), 0 if the endpoint isn't a recognised form.  Negative
 * integers other than the step are normalised against `total` here so
 * the caller can treat the result as a plain 1-based index. */
static int fileprint_resolve_span_slot(const Expr* slot_expr, int which,
                                       int64_t total, int64_t* out) {
    if (slot_expr->type == EXPR_SYMBOL && slot_expr->data.symbol.name == SYM_All) {
        if (which == 0) *out = 1;
        else if (which == 1) *out = total;
        else *out = 1;
        return 1;
    }
    if (slot_expr->type != EXPR_INTEGER) return 0;
    int64_t v = slot_expr->data.integer;
    if (which == 2) { *out = v; return 1; }  /* step keeps sign as-is */
    if (v < 0) v = total + v + 1;
    *out = v;
    return 1;
}

/* Decode the optional second argument into a (start, end, step) triple
 * over 1-based line indices.  Returns 1 on success, 0 if the argument
 * is ill-formed (caller leaves the call unevaluated). */
static int fileprint_decode_selector(const Expr* selector, int64_t total,
                                     int64_t* out_start, int64_t* out_end,
                                     int64_t* out_step) {
    if (selector == NULL) {
        *out_start = 1; *out_end = total; *out_step = 1;
        return 1;
    }
    if (selector->type == EXPR_INTEGER) {
        int64_t n = selector->data.integer;
        if (n >= 0) {
            *out_start = 1;
            *out_end   = (n < total) ? n : total;
            *out_step  = 1;
        } else {
            int64_t want = -n;
            int64_t s = total - want + 1;
            if (s < 1) s = 1;
            *out_start = s;
            *out_end   = total;
            *out_step  = 1;
        }
        return 1;
    }
    if (selector->type == EXPR_FUNCTION
        && selector->data.function.head->type == EXPR_SYMBOL
        && selector->data.function.head->data.symbol.name == SYM_Span)
    {
        size_t argc = selector->data.function.arg_count;
        if (argc < 2 || argc > 3) return 0;

        int64_t start, end, step;
        if (!fileprint_resolve_span_slot(selector->data.function.args[0], 0, total, &start)) return 0;
        if (!fileprint_resolve_span_slot(selector->data.function.args[1], 1, total, &end))   return 0;
        if (argc == 3) {
            if (!fileprint_resolve_span_slot(selector->data.function.args[2], 2, total, &step)) return 0;
            if (step == 0) return 0;
        } else {
            step = 1;
        }
        *out_start = start; *out_end = end; *out_step = step;
        return 1;
    }
    return 0;
}

/* Print one selected range, taking care never to walk off the line
 * array.  When start (or end, in a reverse walk) is out of bounds we
 * silently skip — for the positive-integer cases this lets `n` larger
 * than the file just print everything; for Span we've already validated
 * the endpoints, but defensive bounds checks keep us safe against an
 * empty file with a Span selector. */
static void fileprint_emit(const LineSpan* lines, size_t total,
                           int64_t start, int64_t end, int64_t step) {
    int64_t t = (int64_t)total;
    int64_t last_emitted = 0;  /* 0 = nothing emitted yet */
    if (step > 0) {
        for (int64_t i = start; i <= end; i += step) {
            if (i < 1 || i > t) continue;
            fwrite(lines[i - 1].start, 1, lines[i - 1].len, stdout);
            last_emitted = i;
        }
    } else {
        for (int64_t i = start; i >= end; i += step) {
            if (i < 1 || i > t) continue;
            fwrite(lines[i - 1].start, 1, lines[i - 1].len, stdout);
            last_emitted = i;
        }
    }
    /* If the line we just wrote doesn't end with '\n', flush a newline
     * so the REPL prompt (or the next test assertion) doesn't share a
     * row with file content.  This can only happen when the file's
     * final line is itself unterminated AND we actually emitted it. */
    if (last_emitted > 0) {
        const LineSpan* ls = &lines[last_emitted - 1];
        if (ls->len == 0 || ls->start[ls->len - 1] != '\n') putchar('\n');
    }
    fflush(stdout);
}

/* FilePrint["file"]                 — entire file
 * FilePrint["file", n]              — first/last |n| lines
 * FilePrint["file", m;;n]           — explicit range, step +1
 * FilePrint["file", m;;n;;s]        — explicit range with signed step
 *
 * The evaluator owns `res`; on NULL return we leave it alone, on any
 * other return we hand back ownership.  Returning expr_new_symbol(SYM_Null)
 * mirrors Print's contract. */
Expr* builtin_fileprint(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    Expr* file_arg = res->data.function.args[0];
    if (file_arg->type != EXPR_STRING) return NULL;

    /* Sniff the selector type up front so a bad selector keeps the call
     * unevaluated instead of half-printing the file. */
    Expr* selector = (argc == 2) ? res->data.function.args[1] : NULL;
    if (selector) {
        if (selector->type != EXPR_INTEGER
            && !(selector->type == EXPR_FUNCTION
                 && selector->data.function.head->type == EXPR_SYMBOL
                 && selector->data.function.head->data.symbol.name == SYM_Span))
        {
            return NULL;
        }
    }

    size_t fsize = 0;
    char* buf = fileprint_slurp(file_arg->data.string, &fsize);
    if (!buf) {
        printf("FilePrint::noopen: Cannot open %s.\n", file_arg->data.string);
        return expr_new_symbol(SYM_DollarFailed);
    }

    LineSpan* lines = NULL;
    size_t total = 0;
    if (!fileprint_split_lines(buf, fsize, &lines, &total)) {
        free(buf);
        return expr_new_symbol(SYM_DollarFailed);
    }

    int64_t start = 1, end = (int64_t)total, step = 1;
    if (!fileprint_decode_selector(selector, (int64_t)total, &start, &end, &step)) {
        free(lines);
        free(buf);
        return NULL;  /* unevaluated — bad selector */
    }

    if (total > 0) {
        fileprint_emit(lines, total, start, end, step);
    }

    free(lines);
    free(buf);
    return expr_new_symbol(SYM_Null);
}

/* === FileNameJoin ======================================================
 *
 * Assembles a file name from a list of path components (or canonicalizes a
 * single name).  Pure string manipulation — never touches the filesystem.
 *
 *   FileNameJoin[{"a","b","c"}]         -> "a/b/c"
 *   FileNameJoin[{"a/b","c"}]           -> "a/b/c"   (components may
 *                                          themselves contain separators;
 *                                          they are split and rejoined)
 *   FileNameJoin[{"","a","b"}]          -> "/a/b"    (empty leading
 *                                          component yields an absolute path)
 *   FileNameJoin["a/b/c"]               -> "a/b/c"   (canonicalize one name)
 *   FileNameJoin[{}]                    -> ""
 *   FileNameJoin[..., OperatingSystem->"Windows"|"MacOSX"|"Unix"]
 *
 * "Windows" uses '\\' as the separator and treats a leading "\\server\share"
 * (UNC) prefix as a single unit; "MacOSX"/"Unix" use '/'.  The default is the
 * host operating system's separator.
 */

/* Is `c` a pathname separator?  '/' is always one; '\\' is a separator only
 * when the target operating system is Windows. */
static int fnj_is_sep(char c, int windows) {
    return c == '/' || (windows && c == '\\');
}

/* Join `ncomp` path components into a freshly-malloced NUL-terminated string
 * using separator `sep`.  Each component is split into maximal non-separator
 * runs; empty runs (from leading/trailing/duplicate separators) are dropped,
 * so "a//b" and "a/" collapse cleanly.  A leading empty (or separator-led)
 * first component makes the path absolute.  On Windows, a first component
 * beginning with two separators is preserved verbatim as a UNC prefix.
 * Returns NULL only on allocation failure. */
static char* fnj_build(const char* const* comps, size_t ncomp,
                       char sep, int windows) {
    size_t cap = 4;
    for (size_t i = 0; i < ncomp; i++) {
        cap += strlen(comps[i]) + 1;
    }
    char* out = (char*)malloc(cap);
    if (!out) return NULL;

    size_t pos = 0;
    int wrote = 0;  /* have we emitted at least one segment yet? */

    if (ncomp > 0) {
        const char* c0 = comps[0];
        if (windows && fnj_is_sep(c0[0], windows) && fnj_is_sep(c0[1], windows)) {
            /* UNC: emit the "\\" prefix; the segment walk below picks up the
             * share/host names that follow it. */
            out[pos++] = sep;
            out[pos++] = sep;
        } else if (c0[0] == '\0' || fnj_is_sep(c0[0], windows)) {
            /* Absolute path: a single leading separator.  Later segments are
             * joined without an extra separator because `wrote` stays 0. */
            out[pos++] = sep;
        }
    }

    for (size_t i = 0; i < ncomp; i++) {
        const char* p = comps[i];
        while (*p) {
            while (*p && fnj_is_sep(*p, windows)) p++;   /* skip separators */
            const char* seg = p;
            while (*p && !fnj_is_sep(*p, windows)) p++;  /* one segment */
            size_t seg_len = (size_t)(p - seg);
            if (seg_len == 0) continue;
            if (wrote) out[pos++] = sep;
            memcpy(out + pos, seg, seg_len);
            pos += seg_len;
            wrote = 1;
        }
    }

    out[pos] = '\0';
    return out;
}

/* FileNameJoin[spec]
 * FileNameJoin[spec, OperatingSystem->"os"]
 *
 * `spec` is either a single string (canonicalized) or a List of strings
 * (joined).  Options are trailing Rule[OperatingSystem, "..."] arguments.
 * Leaves the call unevaluated (NULL) on any malformed argument; prints the
 * standard argx message when called with zero arguments. */
Expr* builtin_filenamejoin(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return builtin_arg_error("FileNameJoin", 0, 1, 1);

    /* Determine the target separator from any OperatingSystem option. */
    char sep = HOST_SEP;
    int windows = (HOST_SEP == '\\');
    for (size_t i = 1; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (a->type != EXPR_FUNCTION ||
            a->data.function.arg_count != 2 ||
            a->data.function.head->type != EXPR_SYMBOL ||
            (a->data.function.head->data.symbol.name != SYM_Rule &&
             a->data.function.head->data.symbol.name != SYM_RuleDelayed)) {
            return NULL;  /* an unexpected non-option argument */
        }
        Expr* key = a->data.function.args[0];
        Expr* val = a->data.function.args[1];
        if (key->type != EXPR_SYMBOL ||
            key->data.symbol.name != SYM_OperatingSystem) return NULL;
        if (val->type != EXPR_STRING) return NULL;
        if (strcmp(val->data.string, "Windows") == 0) {
            sep = '\\'; windows = 1;
        } else if (strcmp(val->data.string, "Unix") == 0 ||
                   strcmp(val->data.string, "MacOSX") == 0) {
            sep = '/'; windows = 0;
        } else {
            return NULL;  /* unknown operating system */
        }
    }

    /* Gather the path components: a lone string, or a List of strings. */
    Expr* spec = res->data.function.args[0];
    const char* stackbuf[1];
    const char** comps;
    const char** heapbuf = NULL;
    size_t ncomp;

    if (spec->type == EXPR_STRING) {
        stackbuf[0] = spec->data.string;
        comps = stackbuf;
        ncomp = 1;
    } else if (spec->type == EXPR_FUNCTION &&
               spec->data.function.head->type == EXPR_SYMBOL &&
               spec->data.function.head->data.symbol.name == SYM_List) {
        size_t n = spec->data.function.arg_count;
        for (size_t i = 0; i < n; i++) {
            if (spec->data.function.args[i]->type != EXPR_STRING) return NULL;
        }
        if (n == 0) {
            comps = NULL;
            ncomp = 0;
        } else {
            heapbuf = (const char**)malloc(sizeof(const char*) * n);
            if (!heapbuf) return NULL;
            for (size_t i = 0; i < n; i++) {
                heapbuf[i] = spec->data.function.args[i]->data.string;
            }
            comps = heapbuf;
            ncomp = n;
        }
    } else {
        return NULL;
    }

    char* out = fnj_build(comps, ncomp, sep, windows);
    free(heapbuf);
    if (!out) return NULL;  /* allocation failure — leave call unevaluated */

    Expr* result = expr_new_string(out);
    free(out);
    return result;
}

/* === FileNameSplit =====================================================
 *
 * Splits a file name into its List of path components — the structural
 * inverse of FileNameJoin.  Pure string manipulation; never touches the
 * filesystem.
 *
 *   FileNameSplit["a/b/c"]              -> {"a", "b", "c"}
 *   FileNameSplit["/home/sb/x/"]        -> {"", "home", "sb", "x"}
 *                                          (leading separator => absolute
 *                                          => leading ""; trailing and
 *                                          duplicate separators dropped)
 *   FileNameSplit[""]                   -> {}
 *   FileNameSplit["/"]                  -> {""}
 *   FileNameSplit[..., OperatingSystem->"Windows"|"MacOSX"|"Unix"]
 *
 * On Windows the separator set is {'/','\\'}; a leading "\\host\share"
 * (UNC) prefix is kept as a single part, and a drive like "C:" falls out
 * naturally as an ordinary first part (it contains no separator).  On
 * "Unix"/"MacOSX" the separator is '/'.  The default is the host OS.
 *
 * FileNameJoin[FileNameSplit[x]] reconstructs a canonicalized x for the
 * common cases: both share fnj_is_sep and the same absolute/UNC rules.
 */

/* Duplicate the substring [start, start+len) as a fresh NUL-terminated
 * string.  Returns NULL on allocation failure. */
static char* fns_dup(const char* start, size_t len) {
    char* s = (char*)malloc(len + 1);
    if (!s) return NULL;
    memcpy(s, start, len);
    s[len] = '\0';
    return s;
}

/* Split `path` into its path components under the given separator
 * convention.  Returns a freshly-malloced array of `*out_n` NUL-terminated
 * strings (each malloced); the caller frees every element and the array.
 * An empty path yields zero parts (returns a valid 0-length allocation).
 * Returns NULL only on allocation failure. */
static char** fns_split_build(const char* path, int windows, size_t* out_n) {
    /* Upper bound on the number of parts: each part past the first needs at
     * least one separator, and the optional leading "" adds one.  strlen+2
     * is comfortably safe. */
    size_t cap = strlen(path) + 2;
    char** parts = (char**)malloc(sizeof(char*) * cap);
    if (!parts) return NULL;
    size_t np = 0;

    const char* p = path;

    if (windows && fnj_is_sep(p[0], windows) && fnj_is_sep(p[1], windows)) {
        /* UNC: the "\\host\share" prefix is a single part.  Consume the two
         * leading separators, the host segment, and (if present) a
         * separator plus the share segment. */
        const char* q = p + 2;
        while (*q && !fnj_is_sep(*q, windows)) q++;      /* host segment */
        if (*q) {
            q++;                                         /* skip separator */
            while (*q && !fnj_is_sep(*q, windows)) q++;  /* share segment */
        }
        parts[np] = fns_dup(p, (size_t)(q - p));
        if (!parts[np]) goto oom;
        np++;
        p = q;
    } else if (*p != '\0' && fnj_is_sep(*p, windows)) {
        /* Absolute path (non-UNC): a leading "" part. */
        parts[np] = fns_dup("", 0);
        if (!parts[np]) goto oom;
        np++;
    }

    /* Remaining components: maximal non-separator runs; empty runs from
     * leading/trailing/duplicate separators are dropped. */
    while (*p) {
        while (*p && fnj_is_sep(*p, windows)) p++;   /* skip separators */
        const char* seg = p;
        while (*p && !fnj_is_sep(*p, windows)) p++;  /* one segment */
        size_t seg_len = (size_t)(p - seg);
        if (seg_len == 0) continue;
        parts[np] = fns_dup(seg, seg_len);
        if (!parts[np]) goto oom;
        np++;
    }

    *out_n = np;
    return parts;

oom:
    for (size_t i = 0; i < np; i++) free(parts[i]);
    free(parts);
    return NULL;
}

/* FileNameSplit[spec]
 * FileNameSplit[spec, OperatingSystem->"os"]
 *
 * `spec` must be a single string.  Options are trailing
 * Rule[OperatingSystem, "..."] arguments, decoded exactly as FileNameJoin
 * does.  Leaves the call unevaluated (NULL) on any malformed argument;
 * prints the standard argx message when called with zero arguments. */
Expr* builtin_filenamesplit(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return builtin_arg_error("FileNameSplit", 0, 1, 1);

    /* Determine the target separator from any OperatingSystem option — same
     * decoding as builtin_filenamejoin. */
    int windows = (HOST_SEP == '\\');
    for (size_t i = 1; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (a->type != EXPR_FUNCTION ||
            a->data.function.arg_count != 2 ||
            a->data.function.head->type != EXPR_SYMBOL ||
            (a->data.function.head->data.symbol.name != SYM_Rule &&
             a->data.function.head->data.symbol.name != SYM_RuleDelayed)) {
            return NULL;  /* an unexpected non-option argument */
        }
        Expr* key = a->data.function.args[0];
        Expr* val = a->data.function.args[1];
        if (key->type != EXPR_SYMBOL ||
            key->data.symbol.name != SYM_OperatingSystem) return NULL;
        if (val->type != EXPR_STRING) return NULL;
        if (strcmp(val->data.string, "Windows") == 0) {
            windows = 1;
        } else if (strcmp(val->data.string, "Unix") == 0 ||
                   strcmp(val->data.string, "MacOSX") == 0) {
            windows = 0;
        } else {
            return NULL;  /* unknown operating system */
        }
    }

    Expr* spec = res->data.function.args[0];
    if (spec->type != EXPR_STRING) return NULL;

    size_t np = 0;
    char** parts = fns_split_build(spec->data.string, windows, &np);
    if (!parts) return NULL;  /* allocation failure — leave call unevaluated */

    /* Assemble the List.  expr_new_function memcpy-copies the args array and
     * adopts the element Expr* pointers, so we free our args array (and the
     * transient char* parts) but not the elements. */
    Expr** args = NULL;
    if (np > 0) {
        args = (Expr**)malloc(sizeof(Expr*) * np);
        if (!args) {
            for (size_t i = 0; i < np; i++) free(parts[i]);
            free(parts);
            return NULL;
        }
        for (size_t i = 0; i < np; i++) {
            args[i] = expr_new_string(parts[i]);
        }
    }
    for (size_t i = 0; i < np; i++) free(parts[i]);
    free(parts);

    Expr* list = expr_new_function(expr_new_symbol("List"), args, np);
    free(args);
    return list;
}

/* ------------------------------------------------------------------ */
/* DirectoryName                                                      */
/* ------------------------------------------------------------------ */

static int dn_is_sep(char c) {
    return c == '/' || (HOST_SEP == '\\' && c == '\\');
}

/* Emit Head::tag naming the whole call, e.g. "... in DirectoryName[3]." */
static void files_arg_message(Expr* call, const char* head, const char* tag,
                              const char* what) {
    char* txt = expr_to_string_outputform(call);
    mth_message(head, tag, "%s in %s.", what, txt ? txt : head);
    free(txt);
}

/*
 * DirectoryName["name"] / DirectoryName["name", n]
 *
 * Runs of separators are first collapsed to one ("a//b" -> "a/b"); each step
 * then drops trailing separators and the final component, keeping the
 * separator before it: "/a/b/c.txt" -> "/a/b/", "/a/b/" -> "/a/",
 * "c.txt" -> "", "/" -> "".  n repeats the step.
 */
Expr* builtin_directoryname(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    Expr* a0 = res->data.function.args[0];
    if (a0->type != EXPR_STRING) {
        files_arg_message(res, "DirectoryName", "string",
                          "String expected at position 1");
        return NULL;
    }
    int64_t n = 1;
    if (argc == 2) {
        Expr* a1 = res->data.function.args[1];
        if (a1->type != EXPR_INTEGER || a1->data.integer < 1) {
            files_arg_message(res, "DirectoryName", "intpm",
                "Positive machine-sized integer expected at position 2");
            return NULL;
        }
        n = a1->data.integer;
    }

    const char* src = a0->data.string;
    char* buf = malloc(strlen(src) + 1);
    if (!buf) return NULL;
    size_t len = 0;
    for (const char* p = src; *p; p++)
        if (!(dn_is_sep(*p) && len > 0 && dn_is_sep(buf[len - 1])))
            buf[len++] = *p;

    for (int64_t i = 0; i < n && len > 0; i++) {
        while (len > 0 && dn_is_sep(buf[len - 1])) len--;   /* trailing seps */
        while (len > 0 && !dn_is_sep(buf[len - 1])) len--;  /* last component */
    }
    buf[len] = '\0';
    Expr* out = expr_new_string(buf);
    free(buf);
    return out;
}

/* ------------------------------------------------------------------ */
/* FileNames                                                          */
/* ------------------------------------------------------------------ */

typedef struct {
    Expr** items;
    size_t n, cap;
} FnList;

static void fn_push(FnList* l, const char* s) {
    if (l->n == l->cap) {
        size_t nc = l->cap ? 2 * l->cap : 32;
        Expr** ni = realloc(l->items, nc * sizeof(Expr*));
        if (!ni) return;                    /* OOM: drop the name */
        l->items = ni;
        l->cap = nc;
    }
    l->items[l->n++] = expr_new_string(s);
}

/* dir + separator + name; an empty dir yields name alone, and no separator
 * is doubled when dir already ends in one.  Caller frees. */
static char* fn_join(const char* dir, const char* name) {
    size_t dl = strlen(dir), nl = strlen(name);
    int sep = dl > 0 && !dn_is_sep(dir[dl - 1]);
    char* r = malloc(dl + (size_t)sep + nl + 1);
    if (!r) return NULL;
    memcpy(r, dir, dl);
    if (sep) r[dl] = HOST_SEP;
    memcpy(r + dl + sep, name, nl + 1);
    return r;
}

static int fn_matches(const char* name, RegexRule* rules, int nr) {
    size_t len = strlen(name), ov[2];
    for (int i = 0; i < nr; i++)
        if (regex_match(rules[i].prog, name, len, 0, ov, 1) == 1) return 1;
    return 0;
}

/* List `dir` ("" = the working directory), appending dir/name for every entry
 * whose name matches, and descending into subdirectories while depth > 1
 * (depth < 0 = unlimited).  Subdirectories are found with lstat, so a
 * symlinked directory is listed but not followed (no cycles).  A directory
 * that cannot be opened contributes nothing. */
static void fn_walk(const char* dir, RegexRule* rules, int nr, int64_t depth,
                    FnList* out) {
    DIR* d = opendir(*dir ? dir : ".");
    if (!d) return;
    struct dirent* de;
    while ((de = readdir(d)) != NULL) {
        const char* nm = de->d_name;
        if (strcmp(nm, ".") == 0 || strcmp(nm, "..") == 0) continue;
        char* path = fn_join(dir, nm);
        if (!path) continue;
        if (fn_matches(nm, rules, nr)) fn_push(out, path);
        struct stat st;
        if (depth != 1 && lstat(path, &st) == 0 && S_ISDIR(st.st_mode))
            fn_walk(path, rules, nr, depth < 0 ? depth : depth - 1, out);
        free(path);
    }
    closedir(d);
}

static int fn_str_cmp(const void* a, const void* b) {
    return expr_compare(*(Expr* const*)a, *(Expr* const*)b);
}

/* IgnoreCase -> True / False / Automatic.  Automatic follows the host
 * filesystem: case-insensitive on macOS and Windows.  Returns 1 if `e` is an
 * IgnoreCase rule (and sets *caseless), else 0. */
static int fn_ignorecase_opt(const Expr* e, int* caseless) {
    if (e->type != EXPR_FUNCTION || e->data.function.arg_count != 2 ||
        e->data.function.head->type != EXPR_SYMBOL ||
        (e->data.function.head->data.symbol.name != SYM_Rule &&
         e->data.function.head->data.symbol.name != SYM_RuleDelayed) ||
        e->data.function.args[0]->type != EXPR_SYMBOL ||
        e->data.function.args[0]->data.symbol.name != SYM_IgnoreCase)
        return 0;
    const Expr* v = e->data.function.args[1];
    if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_True) *caseless = 1;
    else if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_False) *caseless = 0;
    else {
#if defined(__APPLE__) || defined(_WIN32)
        *caseless = 1;
#else
        *caseless = 0;
#endif
    }
    return 1;
}

static int fn_is_list(const Expr* e) {
    return e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/* Search every dir for one form.  A string form containing a separator has
 * its directory part moved onto the dirs (a form "sub/x*" searches sub/ for "x*").
 * Returns 0, or -1 when the form cannot be compiled. */
static int fn_search_form(Expr* form, Expr** dirs, size_t ndirs, int64_t depth,
                          int caseless, FnList* out) {
    const char* dpart = "";
    char* dbuf = NULL;
    Expr* patt = form;
    Expr* tail = NULL;
    if (form->type == EXPR_STRING) {
        const char* s = form->data.string;
        const char* last = NULL;
        for (const char* p = s; *p; p++) if (dn_is_sep(*p)) last = p;
        if (last) {
            size_t dl = (size_t)(last - s) + 1;
            dbuf = malloc(dl + 1);
            if (!dbuf) return -1;
            memcpy(dbuf, s, dl);
            dbuf[dl] = '\0';
            dpart = dbuf;
            patt = tail = expr_new_string(last + 1);
        }
    }

    RegexRule* rules;
    int nr = regex_rules_build_opts(patt, REGEX_RULE_ANCHORED | REGEX_RULE_WILDCARDS,
                                    caseless ? REGEX_CASELESS : 0u,
                                    &rules, "FileNames");
    if (tail) expr_free(tail);
    if (nr < 0) { free(dbuf); return -1; }

    for (size_t i = 0; i < ndirs; i++) {
        const char* base = dirs ? dirs[i]->data.string : "";
        char* dir = *dpart ? fn_join(base, dpart) : NULL;
        fn_walk(dir ? dir : base, rules, nr, depth, out);
        free(dir);
    }
    regex_rules_free(rules, nr);
    free(dbuf);
    return 0;
}

/*
 * FileNames[] / FileNames[forms] / FileNames[forms, dirs] /
 * FileNames[forms, dirs, n]
 *
 * forms is a string pattern ("*.m"), a StringExpression, a RegularExpression
 * or a list of them, matched against each entry's name; dirs is a directory
 * or list of them (default: the working directory, with names given bare);
 * n is the search depth, a positive integer or Infinity (default 1; 0 lists
 * the working directory, ignoring dirs).
 * Directories count as entries.  The result is sorted canonically, without
 * duplicates.  Option: IgnoreCase -> Automatic.
 */
Expr* builtin_filenames(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;

    int caseless = 0;
    Expr* defs = symtab_get_options("FileNames");   /* borrowed */
    if (defs && fn_is_list(defs))
        for (size_t i = 0; i < defs->data.function.arg_count; i++)
            fn_ignorecase_opt(defs->data.function.args[i], &caseless);
    while (argc > 0 && fn_ignorecase_opt(a[argc - 1], &caseless)) argc--;
    if (argc > 3) return NULL;

    /* dirs: a string or a list of strings */
    Expr** dirs = NULL;
    size_t ndirs = 1;
    if (argc >= 2) {
        if (a[1]->type == EXPR_STRING) {
            dirs = &a[1];
        } else if (fn_is_list(a[1])) {
            dirs = a[1]->data.function.args;
            ndirs = a[1]->data.function.arg_count;
            for (size_t i = 0; i < ndirs; i++)
                if (dirs[i]->type != EXPR_STRING) return NULL;
        } else {
            return NULL;
        }
    }

    /* depth: n >= 0 or Infinity (< 0 here) */
    int64_t depth = 1;
    if (argc == 3) {
        if (a[2]->type == EXPR_SYMBOL && a[2]->data.symbol.name == SYM_Infinity) {
            depth = -1;
        } else if (a[2]->type == EXPR_INTEGER && a[2]->data.integer >= 0) {
            depth = a[2]->data.integer;
        } else {
            files_arg_message(res, "FileNames", "innf",
                "Non-negative integer or Infinity expected at position 3");
            return NULL;
        }
    }

    /* Depth 0 is WL's odd corner: dirs are ignored and the working directory
     * is listed one level deep, with bare names. */
    if (depth == 0) { dirs = NULL; ndirs = 1; depth = 1; }

    FnList out = { NULL, 0, 0 };
    {
        Expr* star = NULL;
        Expr* forms = argc >= 1 ? a[0] : (star = expr_new_string("*"));
        Expr** fl = fn_is_list(forms) ? forms->data.function.args : &forms;
        size_t nf = fn_is_list(forms) ? forms->data.function.arg_count : 1;
        int bad = 0;
        for (size_t i = 0; i < nf && !bad; i++)
            bad = fn_search_form(fl[i], dirs, ndirs, depth, caseless, &out) != 0;
        if (star) expr_free(star);
        if (bad) {
            for (size_t i = 0; i < out.n; i++) expr_free(out.items[i]);
            free(out.items);
            return NULL;
        }
    }

    /* Sort canonically and drop the duplicates overlapping forms produce. */
    size_t k = 0;
    if (out.n > 1) qsort(out.items, out.n, sizeof(Expr*), fn_str_cmp);
    for (size_t i = 0; i < out.n; i++) {
        if (k > 0 && expr_eq(out.items[k - 1], out.items[i])) {
            expr_free(out.items[i]);
            continue;
        }
        out.items[k++] = out.items[i];
    }
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), out.items, k);
    free(out.items);
    return r;
}

/* ------------------------------------------------------------------ */
/* Directory, SetDirectory, ResetDirectory, DirectoryStack             */
/* ------------------------------------------------------------------ */

/* Directories SetDirectory left, most recent last.  Reachable until exit. */
static char** g_dir_stack = NULL;
static size_t g_dir_count = 0, g_dir_cap = 0;

/* The working directory as a malloc'd string, or NULL.  Grows the buffer
 * rather than relying on getcwd(NULL, 0), which POSIX leaves unspecified. */
static char* cwd_dup(void) {
    size_t cap = 256;
    for (;;) {
        char* buf = malloc(cap);
        if (!buf) return NULL;
        if (getcwd(buf, cap)) return buf;
        free(buf);
        if (cap > (size_t)1 << 20) return NULL;   /* not ERANGE: give up */
        cap *= 2;
    }
}

/* String Expr for the working directory, or $Failed. */
static Expr* cwd_expr(void) {
    char* d = cwd_dup();
    if (!d) return expr_new_symbol(SYM_DollarFailed);
    Expr* r = expr_new_string(d);
    free(d);
    return r;
}

Expr* builtin_directory(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 0) return builtin_arg_error("Directory", argc, 0, 0);
    return cwd_expr();
}

/* $HomeDirectory's current value (borrowed string), or NULL. */
static const char* home_directory(void) {
    Rule* r = symtab_get_own_values("$HomeDirectory");
    return (r && r->replacement && r->replacement->type == EXPR_STRING)
               ? r->replacement->data.string : NULL;
}

/*
 * SetDirectory["dir"] -- make dir the working directory, pushing the old one
 * onto the directory stack; returns the new (canonical) working directory.
 * SetDirectory[] goes to $HomeDirectory.  An unusable dir emits
 * SetDirectory::cdir and returns $Failed, leaving the stack unchanged.
 */
Expr* builtin_setdirectory(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc > 1) return builtin_arg_error("SetDirectory", argc, 0, 1);
    const char* target;
    if (argc == 0) {
        target = home_directory();
        if (!target) return expr_new_symbol(SYM_DollarFailed);
    } else if (res->data.function.args[0]->type == EXPR_STRING) {
        target = res->data.function.args[0]->data.string;
    } else {
        char* txt = expr_to_string_outputform(res->data.function.args[0]);
        mth_message("SetDirectory", "badfile",
                    "The specified argument, %s, should be a valid string or File.",
                    txt ? txt : "?");
        free(txt);
        return NULL;
    }

    char* prev = cwd_dup();
    if (!prev) return expr_new_symbol(SYM_DollarFailed);
    if (chdir(target) != 0) {
        mth_message("SetDirectory", "cdir", "Cannot set current directory to %s.", target);
        free(prev);
        return expr_new_symbol(SYM_DollarFailed);
    }
    if (g_dir_count == g_dir_cap) {
        size_t nc = g_dir_cap ? 2 * g_dir_cap : 8;
        char** ns = realloc(g_dir_stack, nc * sizeof(char*));
        if (!ns) { free(prev); return cwd_expr(); }   /* OOM: moved, not pushed */
        g_dir_stack = ns;
        g_dir_cap = nc;
    }
    g_dir_stack[g_dir_count++] = prev;
    return cwd_expr();
}

/* ResetDirectory[] -- pop the directory stack and return to the directory it
 * held.  On an empty stack emits ResetDirectory::dtop and returns the current
 * directory. */
Expr* builtin_resetdirectory(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 0) return builtin_arg_error("ResetDirectory", argc, 0, 0);
    if (g_dir_count == 0) {
        mth_message("ResetDirectory", "dtop", "Directory stack is empty.");
        return cwd_expr();
    }
    char* dir = g_dir_stack[--g_dir_count];
    if (chdir(dir) != 0)
        mth_message("ResetDirectory", "cdir", "Cannot set current directory to %s.", dir);
    free(dir);
    return cwd_expr();
}

/* DirectoryStack[] -- the directories SetDirectory left, most recent first.
 * Arguments are ignored, as in WL (DirectoryStack[1] is the stack too). */
Expr* builtin_directorystack(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    Expr** items = g_dir_count ? malloc(g_dir_count * sizeof(Expr*)) : NULL;
    if (g_dir_count && !items) return NULL;
    for (size_t i = 0; i < g_dir_count; i++)
        items[i] = expr_new_string(g_dir_stack[g_dir_count - 1 - i]);
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), items, g_dir_count);
    free(items);
    return r;
}

/*
 * Environment["var"] -- the value of the operating-system environment variable
 * var as a string, or $Failed when it is not set.  A non-string argument gives
 * Environment::string and leaves the call unevaluated; a wrong argument count
 * gives Environment::argx.
 */
Expr* builtin_environment(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("Environment", argc, 1, 1);
    Expr* a0 = res->data.function.args[0];
    if (a0->type != EXPR_STRING) {
        files_arg_message(res, "Environment", "string",
                          "String expected at position 1");
        return NULL;
    }
    const char* v = getenv(a0->data.string);
    return v ? expr_new_string(v) : expr_new_symbol(SYM_DollarFailed);
}

/* ---------------- path builtins: ParentDirectory, FileNameTake, ... ---------------- */

static bool is_sep(char c) {
    return c == '/' || c == HOST_SEP;
}

/* Absolute, canonical form of an EXISTING path (symlinks and ".." resolved),
 * malloc'd; NULL when the path does not exist. */
static char* abs_existing(const char* path) {
#ifdef _WIN32
    struct stat st;
    if (stat(path, &st) != 0) return NULL;
    return _fullpath(NULL, path, 0);
#else
    return realpath(path, NULL);
#endif
}

/* "The specified argument, arg, should be a valid string or File." */
static void badfile_message(const char* head, Expr* arg) {
    char* txt = expr_to_string_outputform(arg);
    mth_message(head, "badfile", "The specified argument, %s, should be a valid string or File.",
                txt ? txt : "");
    free(txt);
}

static void fdnfnd_message(const char* head, const char* path) {
    mth_message(head, "fdnfnd", "Directory or file \"%s\" not found.", path);
}

/*
 * ParentDirectory[] / ParentDirectory["dir"] / ParentDirectory["dir", n]
 *
 * The parent of an existing directory, as an absolute path; n repeats the
 * step, and the root is its own parent.  ParentDirectory[] is the parent of
 * the working directory.  A path that is not an existing directory gives
 * ParentDirectory::dirnex and leaves the call unevaluated.
 */
Expr* builtin_parentdirectory(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc > 2) return builtin_arg_error("ParentDirectory", argc, 0, 2);
    int64_t n = 1;
    if (argc == 2) {
        Expr* a1 = res->data.function.args[1];
        if (a1->type != EXPR_INTEGER || a1->data.integer < 1) {
            files_arg_message(res, "ParentDirectory", "intpm",
                "Positive machine-sized integer expected at position 2");
            return NULL;
        }
        n = a1->data.integer;
    }
    char* dir;
    if (argc == 0) {
        dir = cwd_dup();
        if (!dir) return NULL;
    } else {
        Expr* a0 = res->data.function.args[0];
        if (a0->type != EXPR_STRING) {
            files_arg_message(res, "ParentDirectory", "string", "String expected at position 1");
            return NULL;
        }
        dir = abs_existing(a0->data.string);
        struct stat st;
        if (!dir || stat(dir, &st) != 0 || !S_ISDIR(st.st_mode)) {
            free(dir);
            mth_message("ParentDirectory", "dirnex", "Directory %s does not exist.",
                        a0->data.string);
            return NULL;
        }
    }
    for (int64_t i = 0; i < n; i++) {
        size_t len = strlen(dir);
        while (len > 1 && is_sep(dir[len - 1])) dir[--len] = '\0';
        char* last = NULL;
        for (char* c = dir; *c; c++) if (is_sep(*c)) last = c;
        if (!last) break;
        /* Keep the separator of a root: "/" or (Windows) "C:\". */
        size_t keep = (size_t)(last - dir);
        if (keep == 0 || (keep == 2 && dir[1] == ':')) keep++;
        if (keep >= len) break;          /* already the root */
        dir[keep] = '\0';
    }
    Expr* r = expr_new_string(dir);
    free(dir);
    return r;
}

/* Split a file name into components.  A leading separator gives a leading ""
 * component (the root); empty components and trailing separators are dropped.
 * Returns the count; comps[] are malloc'd strings, freed by fnt_free. */
static size_t fnt_split(const char* s, char*** out) {
    size_t n = 0, cap = 8;
    char** comps = malloc(cap * sizeof(char*));
    if (is_sep(*s)) comps[n++] = fns_dup("", 0);
    const char* p = s;
    while (*p) {
        while (*p && is_sep(*p)) p++;
        if (!*p) break;
        const char* q = p;
        while (*q && !is_sep(*q)) q++;
        if (n == cap) { cap *= 2; comps = realloc(comps, cap * sizeof(char*)); }
        comps[n++] = fns_dup(p, (size_t)(q - p));
        p = q;
    }
    *out = comps;
    return n;
}

static void fnt_free(char** comps, size_t n) {
    for (size_t i = 0; i < n; i++) free(comps[i]);
    free(comps);
}

/* Join comps[lo..hi] (0-based, inclusive) with "/"; a leading root component
 * "" gives a leading "/".  An empty range is "". */
static Expr* fnt_join(char** comps, size_t n, int64_t lo, int64_t hi) {
    if (lo < 0) lo = 0;
    if (hi >= (int64_t)n) hi = (int64_t)n - 1;
    if (lo > hi) return expr_new_string("");
    size_t len = 2;
    for (int64_t i = lo; i <= hi; i++) len += strlen(comps[i]) + 1;
    char* buf = malloc(len);
    buf[0] = '\0';
    for (int64_t i = lo; i <= hi; i++) {
        /* The root component is "/" itself; others are "/"-separated. */
        if (comps[i][0] == '\0') strcat(buf, "/");
        else {
            if (i > lo && comps[i - 1][0] != '\0') strcat(buf, "/");
            strcat(buf, comps[i]);
        }
    }
    Expr* r = expr_new_string(buf);
    free(buf);
    return r;
}

/* Resolve a 1-based (negative: from the end) index into 0-based; out of
 * range gives -1 or n. */
static int64_t fnt_index(int64_t k, size_t n) {
    return k > 0 ? k - 1 : (int64_t)n + k;
}

/*
 * FileNameTake["name"] / FileNameTake["name", spec]
 *
 * Components of name as FileNameSplit gives them, rejoined:
 *   (none) the last component       n / -n  the first / last n (all if fewer)
 *   {k}    the k-th (negative: from the end; "" if out of range)
 *   {m, k} components m through k   0       ""
 * Any other spec gives FileNameTake::seqso and leaves the call unevaluated.
 */
Expr* builtin_filenametake(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return builtin_arg_error("FileNameTake", argc, 1, 2);
    Expr* a0 = res->data.function.args[0];
    if (a0->type != EXPR_STRING) {
        files_arg_message(res, "FileNameTake", "string", "String expected at position 1");
        return NULL;
    }
    int64_t lo, hi;
    char** comps;
    size_t n = fnt_split(a0->data.string, &comps);
    Expr* spec = argc == 2 ? res->data.function.args[1] : NULL;
    if (!spec) {
        lo = hi = (int64_t)n - 1;
    } else if (spec->type == EXPR_INTEGER) {
        int64_t k = spec->data.integer;
        if (k >= 0) { lo = 0; hi = k - 1; }
        else        { lo = (int64_t)n + k; hi = (int64_t)n - 1; }
    } else if (spec->type == EXPR_FUNCTION && spec->data.function.head->type == EXPR_SYMBOL
               && spec->data.function.head->data.symbol.name == SYM_List
               && (spec->data.function.arg_count == 1 || spec->data.function.arg_count == 2)
               && spec->data.function.args[0]->type == EXPR_INTEGER
               && spec->data.function.args[spec->data.function.arg_count - 1]->type == EXPR_INTEGER
               && spec->data.function.args[0]->data.integer != 0
               && spec->data.function.args[spec->data.function.arg_count - 1]->data.integer != 0) {
        lo = fnt_index(spec->data.function.args[0]->data.integer, n);
        hi = fnt_index(spec->data.function.args[spec->data.function.arg_count - 1]->data.integer, n);
        /* An index outside the name selects nothing rather than clamping. */
        if (lo < 0 || hi < 0 || lo >= (int64_t)n || hi >= (int64_t)n) lo = 1, hi = 0;
    } else {
        fnt_free(comps, n);
        char* txt = expr_to_string_outputform(spec);
        mth_message("FileNameTake", "seqso",
                    "Sequence specification (+n, -n, {+n}, {-n}, or {m, n}) expected at position 2 in %s.",
                    txt ? txt : "");
        free(txt);
        return NULL;
    }
    Expr* r = fnt_join(comps, n, lo, hi);
    fnt_free(comps, n);
    return r;
}

/*
 * AbsoluteFileName["name"] -- the absolute, canonical path of an existing file
 * or directory.  A missing one gives AbsoluteFileName::fdnfnd and $Failed.
 */
Expr* builtin_absolutefilename(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("AbsoluteFileName", argc, 1, 1);
    Expr* a0 = res->data.function.args[0];
    if (a0->type != EXPR_STRING) {
        badfile_message("AbsoluteFileName", a0);
        return NULL;
    }
    char* abs = abs_existing(a0->data.string);
    if (!abs) {
        fdnfnd_message("AbsoluteFileName", a0->data.string);
        return expr_new_symbol(SYM_DollarFailed);
    }
    Expr* r = expr_new_string(abs);
    free(abs);
    return r;
}

/*
 * FileType["name"] -- File, Directory, or None when nothing exists there.
 * Symbolic links are followed; File["name"] and a leading `~` are accepted.
 */
Expr* builtin_filetype(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("FileType", argc, 1, 1);
    Expr* a0 = res->data.function.args[0];
    char* path = file_arg_path(a0);
    if (!path) {
        badfile_message("FileType", a0);
        return NULL;
    }
    struct stat st;
    int rc = stat(path, &st);
    free(path);
    if (rc != 0) return expr_new_symbol(SYM_None);
    return expr_new_symbol(S_ISDIR(st.st_mode) ? "Directory" : "File");
}

/*
 * DirectoryQ["name"] -- True when name is an existing directory (links
 * followed; File["name"] and `~` accepted), else False.  Any other argument,
 * a list included, is simply False.
 */
Expr* builtin_directoryq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("DirectoryQ", argc, 1, 1);
    char* path = file_arg_path(res->data.function.args[0]);
    bool dir = false;
    if (path && path[0] != '\0') {
        struct stat st;
        dir = stat(path, &st) == 0 && S_ISDIR(st.st_mode);
    }
    free(path);
    return expr_new_symbol(dir ? SYM_True : SYM_False);
}

/* Append path piece `p` to buf (capacity grows), separator-joined: trailing
 * separators of the piece are dropped and exactly one is written after it.
 * A piece made only of separators is the root. */
static void tfn_append(char** buf, size_t* len, size_t* cap, const char* p, bool trailing_sep) {
    size_t pl = strlen(p);
    while (pl > 0 && is_sep(p[pl - 1])) pl--;
    bool root = pl == 0 && p[0] != '\0';
    if (pl == 0 && !root) return;
    while (*len + pl + 2 > *cap) { *cap *= 2; *buf = realloc(*buf, *cap); }
    if (!root) { memcpy(*buf + *len, p, pl); *len += pl; }
    if (trailing_sep || root) (*buf)[(*len)++] = '/';
    (*buf)[*len] = '\0';
}

/*
 * ToFileName[dirs, "name"] / ToFileName[dirs] -- the legacy file-name builder.
 * dirs is a string or a list of strings; each directory is followed by one
 * "/", so ToFileName[{"a", "b"}] is "a/b/" and ToFileName[{"a", "b"}, "c.m"]
 * is "a/b/c.m".  Empty pieces are skipped.
 */
Expr* builtin_tofilename(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return builtin_arg_error("ToFileName", argc, 1, 2);
    Expr* dirs = res->data.function.args[0];
    Expr* name = argc == 2 ? res->data.function.args[1] : NULL;
    if (name && name->type != EXPR_STRING) return NULL;
    bool list = dirs->type == EXPR_FUNCTION && dirs->data.function.head->type == EXPR_SYMBOL
        && dirs->data.function.head->data.symbol.name == SYM_List;
    if (!list && dirs->type != EXPR_STRING) return NULL;
    size_t nd = list ? dirs->data.function.arg_count : 1;
    for (size_t i = 0; list && i < nd; i++)
        if (dirs->data.function.args[i]->type != EXPR_STRING) return NULL;
    size_t len = 0, cap = 64;
    char* buf = malloc(cap);
    buf[0] = '\0';
    for (size_t i = 0; i < nd; i++)
        tfn_append(&buf, &len, &cap, list ? dirs->data.function.args[i]->data.string
                                          : dirs->data.string, true);
    if (name) tfn_append(&buf, &len, &cap, name->data.string, false);
    Expr* r = expr_new_string(buf);
    free(buf);
    return r;
}

/* Delete one file; false (with a message) when it is missing, is a
 * directory, or cannot be removed. */
static bool delete_one(const char* path) {
    struct stat st;
    if (lstat(path, &st) != 0) {
        fdnfnd_message("DeleteFile", path);
        return false;
    }
    if (S_ISDIR(st.st_mode)) {
        mth_message("DeleteFile", "fdir",
                    "The specified path %s refers to a directory; a file path was expected.", path);
        return false;
    }
    if (remove(path) != 0) {
        mth_message("DeleteFile", "nffil", "File %s could not be deleted.", path);
        return false;
    }
    return true;
}

/*
 * DeleteFile["file"] / DeleteFile[{"file1", ...}] -- delete files; Null on
 * success, $Failed (with a message) when any file is missing, is a directory,
 * or cannot be removed.  The remaining files of a list are still deleted.
 */
Expr* builtin_deletefile(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("DeleteFile", argc, 1, 1);
    Expr* a0 = res->data.function.args[0];
    bool list = a0->type == EXPR_FUNCTION && a0->data.function.head->type == EXPR_SYMBOL
        && a0->data.function.head->data.symbol.name == SYM_List;
    size_t n = list ? a0->data.function.arg_count : 1;
    for (size_t i = 0; i < n; i++) {
        Expr* f = list ? a0->data.function.args[i] : a0;
        if (f->type != EXPR_STRING) {
            badfile_message("DeleteFile", f);
            return NULL;
        }
    }
    bool ok = true;
    for (size_t i = 0; i < n; i++)
        ok &= delete_one((list ? a0->data.function.args[i] : a0)->data.string);
    return expr_new_symbol(ok ? SYM_Null : SYM_DollarFailed);
}

/* $HomeDirectory from the environment: HOME, else (Windows) USERPROFILE or
 * HOMEDRIVE + HOMEPATH.  Caller frees; NULL when none is set. */
static char* home_from_env(void) {
    const char* h = getenv("HOME");
    if (!h || !*h) h = getenv("USERPROFILE");
    if (h && *h) {
        char* r = malloc(strlen(h) + 1);
        if (r) strcpy(r, h);
        return r;
    }
    const char* drive = getenv("HOMEDRIVE");
    const char* path = getenv("HOMEPATH");
    if (!drive || !path) return NULL;
    char* r = malloc(strlen(drive) + strlen(path) + 1);
    if (r) { strcpy(r, drive); strcat(r, path); }
    return r;
}

/* $TemporaryDirectory: TMPDIR (POSIX), else TEMP / TMP (Windows), else the
 * host default.  A trailing separator is dropped, so FileNameJoin[{
 * $TemporaryDirectory, "f"}] never doubles it.  Caller frees. */
static char* tempdir_from_env(void) {
    const char* names[] = { "TMPDIR", "TEMP", "TMP" };
    const char* t = NULL;
    for (size_t i = 0; i < sizeof names / sizeof names[0] && !(t && *t); i++)
        t = getenv(names[i]);
#ifdef _WIN32
    if (!t || !*t) t = "C:\\Windows\\Temp";
#else
    if (!t || !*t) t = "/tmp";
#endif
    size_t n = strlen(t);
    while (n > 1 && (t[n - 1] == '/' || t[n - 1] == '\\')) n--;
    char* r = malloc(n + 1);
    if (r) { memcpy(r, t, n); r[n] = '\0'; }
    return r;
}

void files_init(void) {
    /* Docstrings live in info.c alongside the other File-I/O entries.
     * All builtins are Protected with no special evaluation behaviour —
     * argument validation handles non-string inputs by leaving the call
     * unevaluated. */
    symtab_add_builtin("FileExistsQ", builtin_fileexistsq);
    symtab_get_def("FileExistsQ")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("FileSize", builtin_filesize);
    symtab_get_def("FileSize")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("FileExtension", builtin_fileextension);
    symtab_get_def("FileExtension")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("FileBaseName", builtin_filebasename);
    symtab_get_def("FileBaseName")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("FileNameJoin", builtin_filenamejoin);
    symtab_get_def("FileNameJoin")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("FileNameSplit", builtin_filenamesplit);
    symtab_get_def("FileNameSplit")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("FilePrint", builtin_fileprint);
    symtab_get_def("FilePrint")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("DirectoryName", builtin_directoryname);
    symtab_get_def("DirectoryName")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("FileNames", builtin_filenames);
    symtab_get_def("FileNames")->attributes |= ATTR_PROTECTED;
    Expr* ic = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_new_symbol(SYM_IgnoreCase), expr_new_symbol(SYM_Automatic) }, 2);
    symtab_set_options("FileNames",
        expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ ic }, 1));

    symtab_add_builtin("Directory", builtin_directory);
    symtab_get_def("Directory")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("SetDirectory", builtin_setdirectory);
    symtab_get_def("SetDirectory")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("ResetDirectory", builtin_resetdirectory);
    symtab_get_def("ResetDirectory")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("DirectoryStack", builtin_directorystack);
    symtab_get_def("DirectoryStack")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Environment", builtin_environment);
    symtab_get_def("Environment")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("ParentDirectory", builtin_parentdirectory);
    symtab_get_def("ParentDirectory")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("FileNameTake", builtin_filenametake);
    symtab_get_def("FileNameTake")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("AbsoluteFileName", builtin_absolutefilename);
    symtab_get_def("AbsoluteFileName")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("FileType", builtin_filetype);
    symtab_get_def("FileType")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("DirectoryQ", builtin_directoryq);
    symtab_get_def("DirectoryQ")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("ToFileName", builtin_tofilename);
    symtab_get_def("ToFileName")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("DeleteFile", builtin_deletefile);
    symtab_get_def("DeleteFile")->attributes |= ATTR_PROTECTED;

    char* home = home_from_env();
    if (home) {
        Expr* hsym = expr_new_symbol("$HomeDirectory");
        Expr* hval = expr_new_string(home);
        symtab_add_own_value("$HomeDirectory", hsym, hval);
        expr_free(hsym);
        expr_free(hval);
        free(home);
    }
    symtab_get_def("$HomeDirectory")->attributes |= ATTR_PROTECTED;

    char* tmpdir = tempdir_from_env();
    if (tmpdir) {
        Expr* tsym = expr_new_symbol("$TemporaryDirectory");
        Expr* tval = expr_new_string(tmpdir);
        symtab_add_own_value("$TemporaryDirectory", tsym, tval);
        expr_free(tsym);
        expr_free(tval);
        free(tmpdir);
    }
    symtab_get_def("$TemporaryDirectory")->attributes |= ATTR_PROTECTED;

    /* $PathnameSeparator: the host separator FileNameJoin uses by default. */
    {
        const char sep[2] = { HOST_SEP, '\0' };
        Expr* psym = expr_new_symbol("$PathnameSeparator");
        Expr* pval = expr_new_string(sep);
        symtab_add_own_value("$PathnameSeparator", psym, pval);
        expr_free(psym);
        expr_free(pval);
        symtab_get_def("$PathnameSeparator")->attributes |= ATTR_PROTECTED;
    }

    /* $InputFileName: "" outside any file read (see mth_input_file_enter). */
    Expr* ifn_sym = expr_new_symbol("$InputFileName");
    Expr* ifn_val = expr_new_string("");
    symtab_add_own_value("$InputFileName", ifn_sym, ifn_val);
    expr_free(ifn_sym);
    expr_free(ifn_val);
}
