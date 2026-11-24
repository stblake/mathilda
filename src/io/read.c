/*
 * read.c — the shared file-reading engine and the Read[] builtin.
 *
 * Read[stream]                 read one Expression from a stream
 * Read[stream, type]           read one object of the given type
 * Read[stream, {t1, t2, ...}]  read one object of each type into a list
 * Read[stream, structure]      any nested type-structure (Hold[...], f[...],
 *                              {{Number,Number},...}); filled depth-first
 *
 * `stream` is an InputStream[...] object (from OpenRead), a bare "file" string,
 * or File["file"].  A filename that is not already open is opened and left open,
 * so successive Read calls advance the same current point (Close finishes it).
 *
 * Read returns EndOfFile once past end of file, and $Failed for a malformed
 * numeric token or a stream that is not open.  The per-type readers and the
 * separator handling are shared with ReadList (readlist.c), which loops
 * read_structure() to end of file.
 *
 * ANSI C99 only (fopen/fread via streams.c, strtod via the parser); no POSIX
 * symbols, so no feature-test guards are needed (SPEC.md §10).
 */

#include "read.h"
#include "streams.h"
#include "sym_names.h"
#include "symtab.h"
#include "attr.h"
#include "parse.h"
#include "eval.h"
#include "options.h"
#include "files.h"     /* files_slurp_raw */
#include "message.h"   /* mth_message */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Small utilities                                                     */
/* ------------------------------------------------------------------ */

/* malloc a NUL-terminated copy of s[0..n). Avoids POSIX strndup, which glibc
 * hides under -std=c99 (SPEC.md §10). */
static char* dup_n(const char* s, size_t n) {
    char* out = malloc(n + 1);
    if (!out) return NULL;
    if (n) memcpy(out, s, n);
    out[n] = '\0';
    return out;
}

static int expr_is_true(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_True;
}

static int is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List;
}

/* ------------------------------------------------------------------ */
/* Separator configuration                                             */
/* ------------------------------------------------------------------ */

/* Extract a String or List-of-Strings option value into a SepStr array
 * (borrowing the underlying char* pointers). *out is malloc'd (caller frees the
 * array via read_cfg_free, never the borrowed strings). */
static void extract_seps(const Expr* v, SepStr** out, size_t* count) {
    *out = NULL;
    *count = 0;
    if (!v) return;
    if (v->type == EXPR_STRING) {
        SepStr* a = malloc(sizeof(SepStr));
        a[0].s = v->data.string;
        a[0].n = strlen(v->data.string);
        *out = a;
        *count = 1;
        return;
    }
    if (is_list(v)) {
        size_t m = v->data.function.arg_count;
        SepStr* a = malloc((m ? m : 1) * sizeof(SepStr));
        size_t cc = 0;
        for (size_t i = 0; i < m; i++) {
            Expr* e = v->data.function.args[i];
            if (e->type == EXPR_STRING) {
                a[cc].s = e->data.string;
                a[cc].n = strlen(e->data.string);
                cc++;
            }
        }
        *out = a;
        *count = cc;
    }
}

void read_cfg_build(ReadCfg* cfg, const char* msg_head,
                    const Expr* rec, const Expr* word, const Expr* tok,
                    const Expr* nullrec, const Expr* nullword) {
    memset(cfg, 0, sizeof *cfg);
    cfg->msg_head = msg_head ? msg_head : "Read";
    extract_seps(rec, &cfg->rec, &cfg->rec_n);
    extract_seps(word, &cfg->word, &cfg->word_n);
    extract_seps(tok, &cfg->token, &cfg->token_n);
    cfg->null_records = expr_is_true(nullrec);
    cfg->null_words = expr_is_true(nullword);
}

void read_cfg_free(ReadCfg* cfg) {
    free(cfg->rec);
    free(cfg->word);
    free(cfg->token);
    cfg->rec = cfg->word = cfg->token = NULL;
}

/* Length of the longest separator in `arr` matching buf at pos, else 0. */
static size_t sep_len_at(const char* buf, size_t len, size_t pos,
                         const SepStr* arr, size_t n) {
    size_t best = 0;
    for (size_t i = 0; i < n; i++) {
        size_t sl = arr[i].n;
        if (sl > best && pos + sl <= len && memcmp(buf + pos, arr[i].s, sl) == 0)
            best = sl;
    }
    return best;
}

/* Words are delimited by word separators AND record separators. */
static size_t word_sep_len(const char* buf, size_t len, size_t pos, const ReadCfg* c) {
    size_t a = sep_len_at(buf, len, pos, c->word, c->word_n);
    size_t b = sep_len_at(buf, len, pos, c->rec, c->rec_n);
    return a > b ? a : b;
}
static size_t rec_sep_len(const char* buf, size_t len, size_t pos, const ReadCfg* c) {
    return sep_len_at(buf, len, pos, c->rec, c->rec_n);
}
static size_t token_len(const char* buf, size_t len, size_t pos, const ReadCfg* c) {
    return sep_len_at(buf, len, pos, c->token, c->token_n);
}

/* ------------------------------------------------------------------ */
/* Per-type readers                                                    */
/* ------------------------------------------------------------------ */

typedef enum {
    RT_BYTE, RT_CHARACTER, RT_EXPRESSION, RT_NUMBER,
    RT_REAL, RT_RECORD, RT_STRING, RT_WORD
} ReadType;

static int type_from_symbol(const Expr* e, ReadType* out) {
    if (!e || e->type != EXPR_SYMBOL) return 0;
    const char* n = e->data.symbol.name;
    if (n == SYM_Byte)            *out = RT_BYTE;
    else if (n == SYM_Character)  *out = RT_CHARACTER;
    else if (n == SYM_Expression) *out = RT_EXPRESSION;
    else if (n == SYM_Number)     *out = RT_NUMBER;
    else if (n == SYM_Real)       *out = RT_REAL;
    else if (n == SYM_Record)     *out = RT_RECORD;
    else if (n == SYM_String)     *out = RT_STRING;
    else if (n == SYM_Word)       *out = RT_WORD;
    else return 0;
    return 1;
}

/* Next word field: text between word/record separators, with TokenWords emitted
 * as standalone fields. NullWords keeps empty fields between adjacent
 * separators. Returns a malloc'd string (caller frees) or NULL at EOF. */
static char* next_word(ReadCursor* c, const ReadCfg* cfg) {
    for (;;) {
        if (c->pos >= c->len) return NULL;
        /* A TokenWord at the current position is its own field. */
        size_t tw = token_len(c->buf, c->len, c->pos, cfg);
        if (tw > 0) {
            char* w = dup_n(c->buf + c->pos, tw);
            c->pos += tw;
            return w;
        }
        size_t start = c->pos;
        while (c->pos < c->len) {
            if (word_sep_len(c->buf, c->len, c->pos, cfg) > 0) break;
            if (token_len(c->buf, c->len, c->pos, cfg) > 0) break;
            c->pos++;
        }
        if (c->pos > start) {
            char* w = dup_n(c->buf + start, c->pos - start);
            /* consume exactly one following word separator, if present */
            size_t m = word_sep_len(c->buf, c->len, c->pos, cfg);
            if (m > 0) c->pos += m;
            return w;
        }
        /* empty field: sitting on a word separator (token handled above) */
        size_t m = word_sep_len(c->buf, c->len, c->pos, cfg);
        c->pos += (m > 0) ? m : 1;
        if (cfg->null_words) return dup_n("", 0);
        /* NullWords off: skip the empty field and continue */
    }
}

/* Next record: text up to a record separator. NullRecords keeps empty records.
 * Returns malloc'd string or NULL at EOF. */
static char* next_record(ReadCursor* c, const ReadCfg* cfg) {
    for (;;) {
        if (c->pos >= c->len) return NULL;
        size_t start = c->pos;
        while (c->pos < c->len && rec_sep_len(c->buf, c->len, c->pos, cfg) == 0)
            c->pos++;
        if (c->pos > start) {
            char* w = dup_n(c->buf + start, c->pos - start);
            size_t m = rec_sep_len(c->buf, c->len, c->pos, cfg);
            if (m > 0) c->pos += m;
            return w;
        }
        size_t m = rec_sep_len(c->buf, c->len, c->pos, cfg);
        c->pos += (m > 0) ? m : 1;
        if (cfg->null_records) return dup_n("", 0);
    }
}

/* Next String: a line terminated by LF, CRLF or a lone CR. Empty lines are
 * preserved. Returns malloc'd string or NULL at EOF. */
static char* next_string_line(ReadCursor* c) {
    if (c->pos >= c->len) return NULL;
    size_t start = c->pos;
    while (c->pos < c->len && c->buf[c->pos] != '\n' && c->buf[c->pos] != '\r') c->pos++;
    size_t end = c->pos;
    if (c->pos < c->len) {                   /* consume the line break */
        if (c->buf[c->pos] == '\r' && c->pos + 1 < c->len && c->buf[c->pos + 1] == '\n') c->pos++;
        c->pos++;
    }
    return dup_n(c->buf + start, end - start);
}

static void readn_message(const ReadCfg* cfg, const char* tok) {
    mth_message(cfg->msg_head, "readn", "Invalid input `%s` found when a number was expected.", tok ? tok : "");
}

/* Is `tok` a number LITERAL: [+-] digits [. digits] | [+-] . digits, then an optional exponent
 * written e/E, *^ or *10^ with an optional sign? On success writes into `out` (capacity >= 2*strlen+1)
 * the same literal with its exponent respelled as *^, which the expression parser reads as one
 * number. Anything else -- 1/2, 2^3, Print["x"], a symbol -- is not a number in the file, and must
 * never reach the evaluator: Read[s, Number] reads data, it does not run it. */
static int number_literal(const char* tok, char* out) {
    const char* p = tok;
    char* o = out;
    int digits = 0, point = 0;
    if (*p == '+' || *p == '-') *o++ = *p++;
    while (*p >= '0' && *p <= '9') { *o++ = *p++; digits++; }
    if (*p == '.') {
        point = 1;
        *o++ = *p++;
        while (*p >= '0' && *p <= '9') { *o++ = *p++; digits++; }
    }
    if (!digits) return 0;
    int sci = (p[0] == '*' && p[1] == '^') ? 2
            : (p[0] == '*' && p[1] == '1' && p[2] == '0' && p[3] == '^') ? 4
            : (*p == 'e' || *p == 'E') ? 1 : 0;
    if (sci) {
        /* An e/E exponent is the C/Fortran spelling of a REAL (1e2 is 100.); *^ and *10^
         * keep the language's own meaning (exact when the mantissa is). */
        if (sci == 1 && !point) *o++ = '.';
        p += sci;
        *o++ = '*'; *o++ = '^';
        if (*p == '+' || *p == '-') *o++ = *p++;
        if (!(*p >= '0' && *p <= '9')) return 0;
        while (*p >= '0' && *p <= '9') *o++ = *p++;
    }
    *o = '\0';
    return *p == '\0';
}

/* Read one numeric token. `as_real` forces an approximate result. Returns a
 * numeric Expr, or $Failed (with a message) for a malformed token, or NULL at
 * EOF. Owns the returned Expr. */
static Expr* read_number(ReadCursor* c, const ReadCfg* cfg, int as_real) {
    char* tok = NULL;
    do {
        free(tok);
        tok = next_word(c, cfg);
        if (!tok) return NULL;              /* EOF */
    } while (tok[0] == '\0');               /* skip empty (never a number) */

    char* lit = malloc(2 * strlen(tok) + 3);
    if (!lit) { free(tok); return expr_new_symbol(SYM_DollarFailed); }
    Expr* parsed = number_literal(tok, lit) ? parse_expression(lit) : NULL;
    free(lit);
    if (!parsed) {
        readn_message(cfg, tok);
        free(tok);
        return expr_new_symbol(SYM_DollarFailed);
    }
    Expr* v = evaluate(parsed);
    expr_free(parsed);
    if (!v || !expr_is_numeric_like(v)) {
        readn_message(cfg, tok);
        free(tok);
        if (v) expr_free(v);
        return expr_new_symbol(SYM_DollarFailed);
    }
    free(tok);
    if (as_real) {
        Expr* args[1] = { v };              /* expr_new_function adopts v */
        Expr* napp = expr_new_function(expr_new_symbol("N"), args, 1);
        Expr* r = evaluate(napp);
        expr_free(napp);
        return r ? r : expr_new_symbol(SYM_DollarFailed);
    }
    return v;
}

/* Read one complete top-level expression, UNEVALUATED. The caller / evaluator
 * evaluates the constructed result (so Hold[Expression] stays held). NULL at
 * EOF or when no further complete expression can be parsed. */
static Expr* read_expression(ReadCursor* c) {
    const char* p = c->buf + c->pos;
    Expr* parsed = parse_next_expression(&p);
    c->pos = (size_t)(p - c->buf);
    return parsed;                          /* NULL == EOF; unevaluated otherwise */
}

/* Dispatch one read of the given type. NULL means EOF (stop). */
static Expr* read_one(ReadType t, ReadCursor* c, const ReadCfg* cfg) {
    switch (t) {
        case RT_BYTE:
            if (c->pos >= c->len) return NULL;
            return expr_new_integer((unsigned char)c->buf[c->pos++]);
        case RT_CHARACTER: {
            if (c->pos >= c->len) return NULL;
            char s[2] = { c->buf[c->pos++], '\0' };
            return expr_new_string(s);
        }
        case RT_WORD: {
            char* w = next_word(c, cfg);
            if (!w) return NULL;
            Expr* e = expr_new_string(w);
            free(w);
            return e;
        }
        case RT_RECORD: {
            char* w = next_record(c, cfg);
            if (!w) return NULL;
            Expr* e = expr_new_string(w);
            free(w);
            return e;
        }
        case RT_STRING: {
            char* w = next_string_line(c);
            if (!w) return NULL;
            Expr* e = expr_new_string(w);
            free(w);
            return e;
        }
        case RT_NUMBER: return read_number(c, cfg, 0);
        case RT_REAL:   return read_number(c, cfg, 1);
        case RT_EXPRESSION: return read_expression(c);
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
/* The one-object / nested-structure reader                            */
/* ------------------------------------------------------------------ */

/* Scan a type spec: returns 0 if any leaf is not a read type; counts read-type
 * leaves into *leaves. */
static int spec_scan(const Expr* spec, int* leaves) {
    ReadType t;
    if (type_from_symbol(spec, &t)) { (*leaves)++; return 1; }
    if (spec && spec->type == EXPR_FUNCTION) {
        size_t k = spec->data.function.arg_count;
        for (size_t j = 0; j < k; j++)
            if (!spec_scan(spec->data.function.args[j], leaves)) return 0;
        return 1;
    }
    return 0;
}

int read_spec_valid(const Expr* spec) {
    int leaves = 0;
    return spec_scan(spec, &leaves) && leaves > 0;
}

Expr* read_structure(ReadCursor* c, const Expr* spec, const ReadCfg* cfg,
                     int* read_any, int* bad_spec) {
    ReadType t;
    if (type_from_symbol(spec, &t)) {
        Expr* v = read_one(t, c, cfg);
        if (v) { *read_any = 1; return v; }
        /* Leaf at EOF: pad with EndOfFile once anything has been read; otherwise
         * bubble a bare EOF so the caller stops. */
        return *read_any ? expr_new_symbol(SYM_EndOfFile) : NULL;
    }
    if (spec->type == EXPR_FUNCTION) {
        size_t k = spec->data.function.arg_count;
        Expr** out = malloc((k ? k : 1) * sizeof(Expr*));
        for (size_t j = 0; j < k; j++) {
            Expr* child = read_structure(c, spec->data.function.args[j], cfg,
                                         read_any, bad_spec);
            if (*bad_spec || !child) {
                /* bad spec, or bare EOF at the very first leaf (read_any == 0):
                 * unwind and propagate. */
                for (size_t r = 0; r < j; r++) expr_free(out[r]);
                free(out);
                return NULL;
            }
            out[j] = child;
        }
        Expr* r = expr_new_function(expr_copy(spec->data.function.head), out, k);
        free(out);
        return r;
    }
    /* Not a read-type symbol and not a function: an invalid type spec. */
    *bad_spec = 1;
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Read[] builtin                                                      */
/* ------------------------------------------------------------------ */

/* Resolve arg 0 of Read to a live input stream.  On a hard failure a message is
 * printed and *fail is set (return $Failed); a return of NULL with *fail == 0
 * means "not a stream/filename form" (decline). */
static Stream* resolve_input(const Expr* a, int* fail) {
    *fail = 0;
    int id, is_output;
    if (stream_handle(a, &id, &is_output)) {
        if (is_output) {
            mth_message("Read", "openx", "%s is not open for reading.", "OutputStream");
            *fail = 1;
            return NULL;
        }
        Stream* s = stream_by_id(id);
        if (!s) {
            mth_message("Read", "openx", "the stream is not open.");
            *fail = 1;
            return NULL;
        }
        return s;
    }
    const char* name = stream_filename_arg(a);
    if (!name) return NULL;                 /* decline */
    Stream* s = stream_find_by_name(name, 0);
    if (s) return s;
    int nid = stream_open_input(name);
    if (nid < 0) {
        mth_message("Read", "noopen", "Cannot open %s.", name);
        *fail = 1;
        return NULL;
    }
    return stream_by_id(nid);
}

Expr* builtin_read(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    /* Trailing separator options; values default from Options[Read]. */
    const Expr* o_rec = NULL, *o_word = NULL, *o_tok = NULL;
    const Expr* o_nrec = NULL, *o_nword = NULL;
    OptEntry entries[5] = {
        { "RecordSeparators", &o_rec,   NULL },
        { "WordSeparators",   &o_word,  NULL },
        { "TokenWords",       &o_tok,   NULL },
        { "NullRecords",      &o_nrec,  NULL },
        { "NullWords",        &o_nword, NULL },
    };
    size_t new_argc = argc;
    if (!options_extract(res, "Read", entries, 5, &new_argc)) return NULL;
    if (new_argc < 1 || new_argc > 2) return NULL;

    /* Type specification: default is a single Expression. Validate before
     * touching the file, so a malformed spec never opens a stream. */
    Expr* default_spec = NULL;
    const Expr* spec;
    if (new_argc >= 2) {
        spec = res->data.function.args[1];
        if (!read_spec_valid(spec)) return NULL;
    } else {
        default_spec = expr_new_symbol(SYM_Expression);
        spec = default_spec;
    }

    int fail = 0;
    Stream* s = resolve_input(res->data.function.args[0], &fail);
    if (!s) {
        if (default_spec) expr_free(default_spec);
        return fail ? expr_new_symbol(SYM_DollarFailed) : NULL;
    }

    ReadCfg cfg;
    read_cfg_build(&cfg, "Read", o_rec, o_word, o_tok, o_nrec, o_nword);

    ReadCursor c = { s->buf, s->len, s->pos };
    int read_any = 0, bad_spec = 0;
    Expr* result = read_structure(&c, spec, &cfg, &read_any, &bad_spec);
    s->pos = c.pos;                         /* advance the stream's current point */

    read_cfg_free(&cfg);
    if (default_spec) expr_free(default_spec);

    if (bad_spec) {
        if (result) expr_free(result);
        return NULL;                        /* invalid type spec: leave unevaluated */
    }
    return result ? result : expr_new_symbol(SYM_EndOfFile);
}

/* ------------------------------------------------------------------ */
/* Skip                                                                */
/* ------------------------------------------------------------------ */

/* Skip[src, type] / Skip[src, type, n] -- read and discard one (or n) objects
 * of the given type, exactly as Read would read them (same separator options).
 * Returns Null, or EndOfFile if the source ran out first. */
Expr* builtin_skip(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;
    /* Trailing options pass through to Read; an optional count sits at 3. */
    size_t pos_args = 2;
    int64_t n = 1;
    if (argc >= 3 && res->data.function.args[2]->type == EXPR_INTEGER) {
        n = res->data.function.args[2]->data.integer;
        if (n < 0) return NULL;
        pos_args = 3;
    }
    size_t nopts = argc - pos_args;
    Expr** ra = malloc((2 + nopts) * sizeof(Expr*));
    ra[0] = expr_copy(res->data.function.args[0]);
    ra[1] = expr_copy(res->data.function.args[1]);
    for (size_t i = 0; i < nopts; i++) ra[2 + i] = expr_copy(res->data.function.args[pos_args + i]);
    Expr* call = expr_new_function(expr_new_symbol(SYM_Read), ra, 2 + nopts);
    free(ra);
    Expr* out = expr_new_symbol(SYM_Null);
    for (int64_t k = 0; k < n; k++) {
        Expr* r = builtin_read(call);
        if (!r) { expr_free(out); out = NULL; break; }       /* not a valid Read */
        bool stop = (r->type == EXPR_SYMBOL &&
                     (r->data.symbol.name == SYM_EndOfFile || r->data.symbol.name == SYM_DollarFailed));
        if (stop) { expr_free(out); out = r; break; }
        expr_free(r);
    }
    expr_free(call);
    return out;
}

/* ------------------------------------------------------------------ */
/* ReadString                                                          */
/* ------------------------------------------------------------------ */

/* Copy buf[from, to) into a new string Expr, bytes untouched: ReadString
 * returns the file as it is, CRLF line ends included. */
static Expr* rs_chunk(const char* buf, size_t from, size_t to) {
    char* t = malloc(to - from + 1);
    memcpy(t, buf + from, to - from);
    t[to - from] = '\0';
    Expr* e = expr_new_string(t);
    free(t);
    return e;
}

/* First position >= pos where any terminator string starts, or len. */
static size_t rs_find(const char* buf, size_t len, size_t pos, const Expr* term) {
    size_t best = len;
    size_t nt = term->type == EXPR_STRING ? 1 : term->data.function.arg_count;
    for (size_t k = 0; k < nt; k++) {
        const char* ts = (term->type == EXPR_STRING ? term : term->data.function.args[k])->data.string;
        size_t tl = strlen(ts);
        if (tl == 0) continue;
        for (size_t i = pos; i + tl <= len && i < best; i++)
            if (memcmp(buf + i, ts, tl) == 0) { best = i; break; }
    }
    return best;
}

/*
 * ReadString["file"] / ReadString[stream] / ReadString[src, term]
 *
 * The text from the source's current point to its end -- or up to, but not
 * including, the first occurrence of the string term (or of any string in a
 * list of them), which is left unread.  Line ends are normalised to "\n".  At
 * end of file the result is EndOfFile.  A file name is read whole each time
 * without registering a stream; a stream advances.  A file that cannot be
 * opened gives OpenRead::noopen and $Failed.
 */
Expr* builtin_readstring(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    const Expr* term = NULL;
    if (argc == 2) {
        term = res->data.function.args[1];
        bool ok = term->type == EXPR_STRING;
        if (!ok && term->type == EXPR_FUNCTION && term->data.function.head->type == EXPR_SYMBOL
            && term->data.function.head->data.symbol.name == SYM_List) {
            ok = true;
            for (size_t k = 0; k < term->data.function.arg_count && ok; k++)
                ok = term->data.function.args[k]->type == EXPR_STRING;
        }
        if (!ok) return NULL;
    }
    Expr* src = res->data.function.args[0];
    int id, is_output;
    if (stream_handle(src, &id, &is_output)) {
        Stream* s = is_output ? NULL : stream_by_id(id);
        if (!s) return NULL;
        if (s->pos >= s->len) return expr_new_symbol(SYM_EndOfFile);
        size_t end = term ? rs_find(s->buf, s->len, s->pos, term) : s->len;
        Expr* out = rs_chunk(s->buf, s->pos, end);
        s->pos = end;
        return out;
    }
    const char* name = stream_filename_arg(src);
    if (!name) return NULL;
    size_t len = 0;
    char* buf = files_slurp_raw(name, &len);
    if (!buf) {
        mth_message("OpenRead", "noopen", "Cannot open %s.", name);
        return expr_new_symbol(SYM_DollarFailed);
    }
    Expr* out;
    if (len == 0) out = expr_new_symbol(SYM_EndOfFile);
    else out = rs_chunk(buf, 0, term ? rs_find(buf, len, 0, term) : len);
    free(buf);
    return out;
}

/* ------------------------------------------------------------------ */
/* Registration                                                        */
/* ------------------------------------------------------------------ */

/* Build a List[String, ...] from a NULL-terminated C string array. */
static Expr* string_list(const char** items) {
    size_t n = 0;
    while (items[n]) n++;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = expr_new_string(items[i]);
    Expr* l = expr_new_function(expr_new_symbol(SYM_List), a, n);
    free(a);
    return l;
}

/* Build Rule[name, value]; adopts value. */
static Expr* rule_of(const char* name_sym, Expr* value) {
    Expr* a[2] = { expr_new_symbol(name_sym), value };
    return expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
}

void read_init(void) {
    symtab_add_builtin("Read", builtin_read);
    symtab_get_def("Read")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("Skip", builtin_skip);
    symtab_get_def("Skip")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Skip",
        "Skip[stream, type]\n"
        "\tskips one object of the given type (as Read would read it) in an input stream.\n"
        "Skip[stream, type, n] skips n objects.  Returns Null, or EndOfFile if the stream ends first.");
    symtab_add_builtin("ReadString", builtin_readstring);
    symtab_get_def("ReadString")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ReadString",
        "ReadString[\"file\"]\n"
        "\treads the complete contents of a file and returns them as a string.\n"
        "ReadString[stream] reads from the stream's current point to the end.\n"
        "ReadString[src, term] reads up to, but not including, the first occurrence of the string\n"
        "term (or any of a list of strings), which is left unread.  Line ends are normalised to\n"
        "\"\\n\".  At end of file the result is EndOfFile.");

    /* Same default separator settings as ReadList; SetOptions[Read, ...] works
     * with no extra code because Options[Read] is the single source of truth. */
    static const char* rec_defaults[]  = { "\r\n", "\n", "\r", NULL };
    static const char* word_defaults[] = { " ", "\t", NULL };
    static const char* tok_defaults[]  = { NULL };

    Expr* opts_args[5] = {
        rule_of(SYM_RecordSeparators, string_list(rec_defaults)),
        rule_of(SYM_WordSeparators,   string_list(word_defaults)),
        rule_of(SYM_TokenWords,       string_list(tok_defaults)),
        rule_of(SYM_NullRecords,      expr_new_symbol(SYM_False)),
        rule_of(SYM_NullWords,        expr_new_symbol(SYM_False)),
    };
    Expr* opts = expr_new_function(expr_new_symbol(SYM_List), opts_args, 5);
    symtab_set_options("Read", opts);
}
