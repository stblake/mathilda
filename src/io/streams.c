/*
 * streams.c — the input/output stream registry and stream-management builtins.
 *
 * A process-global registry of open streams (see streams.h).  Input streams
 * hold the whole file in a buffer plus a moving read point; output streams hold
 * an open FILE*.  Each is handed to the language as an inert InputStream[name,
 * id] / OutputStream[name, id] object whose integer id addresses the registry
 * entry across calls.  Entries are freed by Close or, for anything still open at
 * exit, by an atexit hook, so nothing leaks under valgrind.
 *
 * Builtins: OpenRead, OpenWrite, OpenAppend, Close, Streams, StreamPosition,
 * SetStreamPosition, Write, WriteString.
 *
 * ANSI C99 only (fopen/fread/fseek/ftell/fclose/fflush/realloc/atexit); no
 * POSIX symbols, so no feature-test guards are needed (SPEC.md §10).
 */

#include "streams.h"
#include "message.h"   /* mth_message */
#include "sym_names.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Registry                                                            */
/* ------------------------------------------------------------------ */

static Stream* g_streams = NULL;   /* dynamically grown; id == 0 marks free */
static size_t  g_cap = 0;
static int     g_next_id = 1;      /* monotonic handle allocator */
static int     g_atexit_done = 0;

static void streams_shutdown(void) {
    for (size_t i = 0; i < g_cap; i++) {
        Stream* s = &g_streams[i];
        if (s->id == 0) continue;
        free(s->name);
        free(s->buf);
        if (s->fp && !s->borrowed) fclose(s->fp);
    }
    free(g_streams);
    g_streams = NULL;
    g_cap = 0;
}

/* Return a free slot (id == 0), growing the registry if needed. */
static Stream* alloc_slot(void) {
    if (!g_atexit_done) { atexit(streams_shutdown); g_atexit_done = 1; }
    for (size_t i = 0; i < g_cap; i++)
        if (g_streams[i].id == 0) return &g_streams[i];
    size_t cap = g_cap ? g_cap * 2 : 4;
    g_streams = realloc(g_streams, cap * sizeof(Stream));
    for (size_t i = g_cap; i < cap; i++) memset(&g_streams[i], 0, sizeof(Stream));
    Stream* slot = &g_streams[g_cap];
    g_cap = cap;
    return slot;
}

/* malloc a NUL-terminated copy of s. */
static char* dup_str(const char* s) {
    size_t n = strlen(s);
    char* r = malloc(n + 1);
    memcpy(r, s, n + 1);
    return r;
}

/* Slurp a whole file into a malloc'd buffer. Returns 1 and fills the buffer and
 * length outputs on success (caller frees the buffer), 0 on failure. */
static int slurp_file(const char* path, char** buf, size_t* len) {
    FILE* fp = fopen(path, "rb");
    if (!fp) return 0;
    if (fseek(fp, 0, SEEK_END) != 0) { fclose(fp); return 0; }
    long sz = ftell(fp);
    if (sz < 0) { fclose(fp); return 0; }
    if (fseek(fp, 0, SEEK_SET) != 0) { fclose(fp); return 0; }
    char* b = malloc((size_t)sz + 1);
    if (!b) { fclose(fp); return 0; }
    size_t n = fread(b, 1, (size_t)sz, fp);
    b[n] = '\0';
    fclose(fp);
    *buf = b;
    *len = n;
    return 1;
}

int stream_open_input(const char* name) {
    char* buf = NULL;
    size_t len = 0;
    if (!slurp_file(name, &buf, &len)) return -1;
    Stream* s = alloc_slot();
    memset(s, 0, sizeof(Stream));
    s->id = g_next_id++;
    s->is_output = 0;
    s->name = dup_str(name);
    s->buf = buf;
    s->len = len;
    s->pos = 0;
    return s->id;
}

/* "stdout" and "stderr" name the process's own streams, not files: Write /
 * WriteString["stdout", ...] print (to the same FILE* Print uses), and never
 * create a file called "stdout" in the working directory. */
static FILE* standard_stream(const char* name) {
    if (strcmp(name, "stdout") == 0) return stdout;
    if (strcmp(name, "stderr") == 0) return stderr;
    return NULL;
}

int stream_open_output(const char* name, int append) {
    FILE* std = standard_stream(name);
    FILE* fp = std ? std : fopen(name, append ? "ab" : "wb");
    if (!fp) return -1;
    Stream* s = alloc_slot();
    memset(s, 0, sizeof(Stream));
    s->id = g_next_id++;
    s->is_output = 1;
    s->name = dup_str(name);
    s->fp = fp;
    s->borrowed = std != NULL;
    return s->id;
}

Stream* stream_by_id(int id) {
    if (id <= 0) return NULL;
    for (size_t i = 0; i < g_cap; i++)
        if (g_streams[i].id == id) return &g_streams[i];
    return NULL;
}

Stream* stream_find_by_name(const char* name, int want_output) {
    Stream* best = NULL;
    for (size_t i = 0; i < g_cap; i++) {
        Stream* s = &g_streams[i];
        if (s->id == 0) continue;
        if (want_output >= 0 && s->is_output != want_output) continue;
        if (strcmp(s->name, name) != 0) continue;
        if (!best || s->id > best->id) best = s;   /* most recently opened */
    }
    return best;
}

int stream_close_id(int id) {
    Stream* s = stream_by_id(id);
    if (!s) return 0;
    free(s->name);
    free(s->buf);
    if (s->fp && s->borrowed) fflush(s->fp);
    else if (s->fp) fclose(s->fp);
    memset(s, 0, sizeof(Stream));   /* id -> 0: slot free */
    return 1;
}

size_t stream_capacity(void) { return g_cap; }
Stream* stream_slot(size_t i) { return (i < g_cap) ? &g_streams[i] : NULL; }

/* ------------------------------------------------------------------ */
/* Argument helpers                                                    */
/* ------------------------------------------------------------------ */

int stream_handle(const Expr* a, int* id, int* is_output) {
    if (!a || a->type != EXPR_FUNCTION) return 0;
    Expr* h = a->data.function.head;
    if (h->type != EXPR_SYMBOL) return 0;
    int outp;
    if (h->data.symbol.name == SYM_InputStream) outp = 0;
    else if (h->data.symbol.name == SYM_OutputStream) outp = 1;
    else return 0;
    if (a->data.function.arg_count < 2) return 0;
    Expr* idn = a->data.function.args[1];
    if (idn->type != EXPR_INTEGER) return 0;
    if (id) *id = (int)idn->data.integer;
    if (is_output) *is_output = outp;
    return 1;
}

const char* stream_filename_arg(const Expr* a) {
    if (!a) return NULL;
    if (a->type == EXPR_STRING) return a->data.string;
    if (a->type == EXPR_FUNCTION
        && a->data.function.head->type == EXPR_SYMBOL
        && a->data.function.head->data.symbol.name == SYM_File
        && a->data.function.arg_count == 1
        && a->data.function.args[0]->type == EXPR_STRING)
        return a->data.function.args[0]->data.string;
    return NULL;
}

Expr* stream_make_object(const Stream* s) {
    Expr* a[2] = { expr_new_string(s->name), expr_new_integer(s->id) };
    return expr_new_function(
        expr_new_symbol(s->is_output ? SYM_OutputStream : SYM_InputStream), a, 2);
}

/* ------------------------------------------------------------------ */
/* Open / Close / Streams                                              */
/* ------------------------------------------------------------------ */

static Expr* open_common(Expr* res, int output, int append, const char* who) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    const char* name = stream_filename_arg(res->data.function.args[0]);
    if (!name) return NULL;
    int id = output ? stream_open_output(name, append) : stream_open_input(name);
    if (id < 0) {
        mth_message(who, "noopen", "Cannot open %s.", name);
        return expr_new_symbol(SYM_DollarFailed);
    }
    return stream_make_object(stream_by_id(id));
}

/* StringToStream["string"] -- an input stream over a copy of the string, as
 * InputStream["String", id].  Read, ReadList, ReadString, Skip,
 * StreamPosition and Close all work on it exactly as on a file stream (an
 * input stream is an in-memory buffer either way).  A non-string argument
 * gives StringToStream::string and $Failed. */
Expr* builtin_stringtostream(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* a = res->data.function.args[0];
    if (a->type != EXPR_STRING) {
        char* txt = expr_to_string_outputform(res);
        mth_message("StringToStream", "string", "String expected at position 1 in %s.", txt ? txt : "");
        free(txt);
        return expr_new_symbol(SYM_DollarFailed);
    }
    size_t len = strlen(a->data.string);
    char* buf = malloc(len + 1);
    if (!buf) return NULL;
    memcpy(buf, a->data.string, len + 1);
    Stream* st = alloc_slot();
    memset(st, 0, sizeof(Stream));
    st->id = g_next_id++;
    st->is_output = 0;
    st->name = dup_str("String");
    st->buf = buf;
    st->len = len;
    st->pos = 0;
    return stream_make_object(st);
}

Expr* builtin_openread(Expr* res)   { return open_common(res, 0, 0, "OpenRead"); }
Expr* builtin_openwrite(Expr* res)  { return open_common(res, 1, 0, "OpenWrite"); }
Expr* builtin_openappend(Expr* res) { return open_common(res, 1, 1, "OpenAppend"); }

Expr* builtin_close(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* a = res->data.function.args[0];

    int id, is_output;
    Stream* s = NULL;
    if (stream_handle(a, &id, &is_output)) {
        s = stream_by_id(id);
    } else {
        const char* name = stream_filename_arg(a);
        if (!name) return NULL;             /* not a stream/filename: decline */
        s = stream_find_by_name(name, -1);
    }
    if (!s) {
        mth_message("Close", "stream", "the stream is not open.");
        return expr_new_symbol(SYM_DollarFailed);
    }
    Expr* name_result = expr_new_string(s->name);
    stream_close_id(s->id);
    return name_result;                     /* Close returns the file name */
}

Expr* builtin_streams(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc > 1) return NULL;

    const char* filter = NULL;
    if (argc == 1) {
        filter = stream_filename_arg(res->data.function.args[0]);
        if (!filter) return NULL;
    }

    size_t cap = 8, cnt = 0;
    Expr** arr = malloc(cap * sizeof(Expr*));
    for (size_t i = 0; i < g_cap; i++) {
        Stream* s = &g_streams[i];
        if (s->id == 0) continue;
        if (filter && strcmp(s->name, filter) != 0) continue;
        if (cnt == cap) { cap *= 2; arr = realloc(arr, cap * sizeof(Expr*)); }
        arr[cnt++] = stream_make_object(s);
    }
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), arr, cnt);
    free(arr);
    return list;
}

/* ------------------------------------------------------------------ */
/* StreamPosition / SetStreamPosition                                  */
/* ------------------------------------------------------------------ */

/* Resolve a stream argument that must already be open (no auto-open). */
static Stream* resolve_open(const Expr* a, const char* who) {
    int id, is_output;
    Stream* s = NULL;
    if (stream_handle(a, &id, &is_output)) s = stream_by_id(id);
    else {
        const char* name = stream_filename_arg(a);
        if (name) s = stream_find_by_name(name, -1);
    }
    if (!s) mth_message(who, "stream", "the stream is not open.");
    return s;
}

Expr* builtin_streamposition(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Stream* s = resolve_open(res->data.function.args[0], "StreamPosition");
    if (!s) return expr_new_symbol(SYM_DollarFailed);
    long pos = s->is_output ? ftell(s->fp) : (long)s->pos;
    if (pos < 0) pos = 0;
    return expr_new_integer(pos);
}

Expr* builtin_setstreamposition(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Stream* s = resolve_open(res->data.function.args[0], "SetStreamPosition");
    if (!s) return expr_new_symbol(SYM_DollarFailed);

    Expr* n = res->data.function.args[1];
    int to_end = (n->type == EXPR_SYMBOL && n->data.symbol.name == SYM_Infinity);
    long want = 0;
    if (!to_end) {
        if (n->type != EXPR_INTEGER) return NULL;
        want = (long)n->data.integer;
        if (want < 0) want = 0;
    }

    if (s->is_output) {
        if (to_end) fseek(s->fp, 0, SEEK_END);
        else        fseek(s->fp, want, SEEK_SET);
        long pos = ftell(s->fp);
        if (pos < 0) pos = 0;
        return expr_new_integer(pos);
    }
    /* input: clamp into the buffer */
    size_t p = to_end ? s->len : (size_t)want;
    if (p > s->len) p = s->len;
    s->pos = p;
    return expr_new_integer((long)p);
}

/* ------------------------------------------------------------------ */
/* Write / WriteString                                                 */
/* ------------------------------------------------------------------ */

/* Resolve a write target: an OutputStream object, or a filename/File[...] that
 * is auto-opened (truncating) and kept open.  *fail as in resolve_input. */
static Stream* resolve_output(const Expr* a, const char* who, int* fail) {
    *fail = 0;
    int id, is_output;
    if (stream_handle(a, &id, &is_output)) {
        if (!is_output) {
            mth_message(who, "openx", "an InputStream is not open for writing.");
            *fail = 1;
            return NULL;
        }
        Stream* s = stream_by_id(id);
        if (!s) {
            mth_message(who, "openx", "the stream is not open.");
            *fail = 1;
            return NULL;
        }
        return s;
    }
    const char* name = stream_filename_arg(a);
    if (!name) return NULL;                 /* decline */
    Stream* s = stream_find_by_name(name, 1);
    if (s) return s;
    int nid = stream_open_output(name, 0);
    if (nid < 0) {
        mth_message(who, "noopen", "Cannot open %s.", name);
        *fail = 1;
        return NULL;
    }
    return stream_by_id(nid);
}

/* Write[out, e1, e2, ...] — each eᵢ in re-readable text, then one newline.
 * On "stdout"/"stderr" the text is OutputForm (strings unquoted), as Print. */
Expr* builtin_write(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;
    int fail = 0;
    Stream* s = resolve_output(res->data.function.args[0], "Write", &fail);
    if (!s) return fail ? expr_new_symbol(SYM_DollarFailed) : NULL;

    size_t argc = res->data.function.arg_count;
    for (size_t i = 1; i < argc; i++) {
        Expr* e = res->data.function.args[i];
        /* A file stream writes InputForm, its default FormatType in the
         * reference (Write[f, 1/100 a + 5 x] writes "a/100 + 5*x"). */
        char* str;
        if (s->borrowed) {
            str = expr_to_string_outputform(e);
        } else {
            Expr* inner = expr_copy(e);
            Expr* w = expr_new_function(expr_new_symbol(SYM_InputForm), &inner, 1);
            str = expr_to_string(w);
            expr_free(w);
        }
        if (str) { fputs(str, s->fp); free(str); }
    }
    fputc('\n', s->fp);
    fflush(s->fp);
    return expr_new_symbol(SYM_Null);
}

/* WriteString[out, s1, ...] — strings written raw, no quotes, no newline. */
Expr* builtin_writestring(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;
    int fail = 0;
    Stream* s = resolve_output(res->data.function.args[0], "WriteString", &fail);
    if (!s) return fail ? expr_new_symbol(SYM_DollarFailed) : NULL;

    size_t argc = res->data.function.arg_count;
    for (size_t i = 1; i < argc; i++) {
        Expr* e = res->data.function.args[i];
        if (e->type == EXPR_STRING) {
            fputs(e->data.string, s->fp);   /* raw, no surrounding quotes */
        } else {
            char* str = expr_to_string(e);
            if (str) { fputs(str, s->fp); free(str); }
        }
    }
    fflush(s->fp);
    return expr_new_symbol(SYM_Null);
}

/* ------------------------------------------------------------------ */
/* Registration                                                        */
/* ------------------------------------------------------------------ */

static void protect(const char* name) {
    symtab_get_def(name)->attributes |= ATTR_PROTECTED;
}

void streams_init(void) {
    symtab_add_builtin("OpenRead", builtin_openread);           protect("OpenRead");
    symtab_add_builtin("OpenWrite", builtin_openwrite);         protect("OpenWrite");
    symtab_add_builtin("OpenAppend", builtin_openappend);       protect("OpenAppend");
    symtab_add_builtin("Close", builtin_close);                 protect("Close");
    symtab_add_builtin("StringToStream", builtin_stringtostream); protect("StringToStream");
    symtab_set_docstring("StringToStream",
        "StringToStream[\"string\"]\n"
        "\topens an input stream for reading from a string, as InputStream[\"String\", id].\n"
        "Read, ReadList, ReadString, Skip, StreamPosition and Close work on it as on a file stream.");
    symtab_add_builtin("Streams", builtin_streams);             protect("Streams");
    symtab_add_builtin("StreamPosition", builtin_streamposition);
    protect("StreamPosition");
    symtab_add_builtin("SetStreamPosition", builtin_setstreamposition);
    protect("SetStreamPosition");
    symtab_add_builtin("Write", builtin_write);                 protect("Write");
    symtab_add_builtin("WriteString", builtin_writestring);     protect("WriteString");

    /* Inert normal-form objects: no builtin, just Protected so the head is a
     * stable head the evaluator leaves alone. */
    protect("InputStream");
    protect("OutputStream");
    protect("File");
}
