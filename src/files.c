/*
 * files.c — Filesystem-predicate and path-string builtins.
 *
 * Builtins:
 *   FileExistsQ["name"]    — touches the filesystem (lstat).
 *   FileExtension["name"]  — pure string manipulation.
 *   FileBaseName["name"]   — pure string manipulation.
 *   FilePrint["name", ...] — streams file contents to stdout.
 *   FileNameJoin[{...}]    — assembles a file name from parts.
 *   FileNameSplit["name"]  — splits a file name into its parts.
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
 * it implicitly but accepts the macro. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "files.h"
#include "expr.h"
#include "symtab.h"
#include "sym_names.h"
#include "attr.h"
#include "common.h"

#include <sys/stat.h>
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

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING) return NULL;

    struct stat st;
    int exists = (lstat(arg->data.string, &st) == 0);
    return expr_new_symbol(exists ? "True" : "False");
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
            strcmp(key->data.symbol.name, "OperatingSystem") != 0) return NULL;
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
            strcmp(key->data.symbol.name, "OperatingSystem") != 0) return NULL;
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

void files_init(void) {
    /* Docstrings live in info.c alongside the other File-I/O entries.
     * All builtins are Protected with no special evaluation behaviour —
     * argument validation handles non-string inputs by leaving the call
     * unevaluated. */
    symtab_add_builtin("FileExistsQ", builtin_fileexistsq);
    symtab_get_def("FileExistsQ")->attributes |= ATTR_PROTECTED;

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
}
