/*
 * readlist.c — ReadList[] file-reading builtin.
 *
 * ReadList["file"]                 read all remaining Expressions -> List
 * ReadList["file", type]           read `type` until EOF -> flat List
 * ReadList["file", {t1,...,tk}]    read one of each type per pass -> List of
 *                                  k-element sublists (EndOfFile-padded when
 *                                  EOF arrives partway through a pass)
 * ReadList["file", types, n]       stop after n top-level objects / passes
 *
 * Types: Byte, Character, Expression, Number, Real, Record, String, Word.
 * Options (trailing rules, defaults from Options[ReadList]): RecordSeparators,
 * WordSeparators, TokenWords, NullRecords, NullWords.
 *
 * The whole file is slurped into a buffer and consumed by a ReadCursor.  No
 * stream objects exist, so a file is always opened and closed per call; the
 * Mathematica "already-open stream" behaviour does not apply.
 *
 * Uses only ANSI C99 (fopen/fread/strtod via the parser/isspace); no POSIX
 * symbols, so no feature-test guards are required (see SPEC.md §10).
 */

#include "readlist.h"
#include "sym_names.h"
#include "symtab.h"
#include "attr.h"
#include "parse.h"
#include "eval.h"
#include "options.h"

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

/* Slurp a whole file into a malloc'd buffer. Returns 1 and fills the buffer
 * and length outputs on success (caller frees the buffer), 0 if the file
 * cannot be opened or read. */
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

/* One separator string, borrowed from an option Expr (valid for the call). */
typedef struct { const char* s; size_t n; } SepStr;

typedef struct {
    SepStr* rec;   size_t rec_n;    /* RecordSeparators */
    SepStr* word;  size_t word_n;   /* WordSeparators   */
    SepStr* token; size_t token_n;  /* TokenWords       */
    int null_records;
    int null_words;
} ReadCfg;

/* Extract a String or List-of-Strings option value into a SepStr array
 * (borrowing the underlying char* pointers). *out is malloc'd (caller frees
 * the array, never the borrowed strings). */
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

static void free_cfg(ReadCfg* cfg) {
    free(cfg->rec);
    free(cfg->word);
    free(cfg->token);
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
/* The read cursor and per-type readers                                */
/* ------------------------------------------------------------------ */

typedef struct { const char* buf; size_t len; size_t pos; } ReadCursor;

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

/* Next word field: text between word/record separators, with TokenWords
 * emitted as standalone fields. NullWords keeps empty fields between adjacent
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

/* Next record: text up to a record separator. NullRecords keeps empty
 * records. Returns malloc'd string or NULL at EOF. */
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

/* Next String: a line terminated by '\n', trailing '\r' stripped. Empty
 * lines are preserved. Returns malloc'd string or NULL at EOF. */
static char* next_string_line(ReadCursor* c) {
    if (c->pos >= c->len) return NULL;
    size_t start = c->pos;
    while (c->pos < c->len && c->buf[c->pos] != '\n') c->pos++;
    size_t end = c->pos;
    if (c->pos < c->len) c->pos++;           /* consume the newline */
    if (end > start && c->buf[end - 1] == '\r') end--;
    return dup_n(c->buf + start, end - start);
}

static void readn_message(const char* tok) {
    printf("ReadList::readn: Invalid input `%s` found when a number was expected.\n",
           tok ? tok : "");
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

    Expr* parsed = parse_expression(tok);
    if (!parsed) {
        readn_message(tok);
        free(tok);
        return expr_new_symbol(SYM_DollarFailed);
    }
    Expr* v = evaluate(parsed);
    expr_free(parsed);
    if (!v || !expr_is_numeric_like(v)) {
        readn_message(tok);
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

/* Read and evaluate one complete top-level expression. NULL at EOF or when no
 * further complete expression can be parsed (matching Get's semantics). */
static Expr* read_expression(ReadCursor* c) {
    const char* p = c->buf + c->pos;
    Expr* parsed = parse_next_expression(&p);
    c->pos = (size_t)(p - c->buf);
    if (!parsed) return NULL;
    Expr* v = evaluate(parsed);
    expr_free(parsed);
    return v ? v : expr_new_symbol(SYM_Null);
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
/* Drivers                                                             */
/* ------------------------------------------------------------------ */

/* Single-type mode: read `type` until EOF (or n objects) -> flat List. */
static Expr* read_flat(ReadCursor* c, ReadType type, const ReadCfg* cfg,
                       int have_n, long n) {
    size_t cap = 8, cnt = 0;
    Expr** arr = malloc(cap * sizeof(Expr*));
    while (!have_n || (long)cnt < n) {
        Expr* v = read_one(type, c, cfg);
        if (!v) break;
        if (cnt == cap) { cap *= 2; arr = realloc(arr, cap * sizeof(Expr*)); }
        arr[cnt++] = v;
    }
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), arr, cnt);
    free(arr);
    return list;
}

/* Group mode: read one of each type per pass -> List of k-element sublists.
 * On EOF partway through a pass, remaining slots are EndOfFile and that pass
 * is the last one. */
static Expr* read_groups(ReadCursor* c, const ReadType* types, size_t k,
                         const ReadCfg* cfg, int have_n, long n) {
    size_t cap = 8, cnt = 0;
    Expr** arr = malloc(cap * sizeof(Expr*));
    while (!have_n || (long)cnt < n) {
        Expr* first = read_one(types[0], c, cfg);
        if (!first) break;                  /* clean end at a pass boundary */
        Expr** grp = malloc(k * sizeof(Expr*));
        grp[0] = first;
        int truncated = 0;
        for (size_t j = 1; j < k; j++) {
            Expr* v = read_one(types[j], c, cfg);
            if (!v) {
                for (size_t r = j; r < k; r++) grp[r] = expr_new_symbol(SYM_EndOfFile);
                truncated = 1;
                break;
            }
            grp[j] = v;
        }
        Expr* sub = expr_new_function(expr_new_symbol(SYM_List), grp, k);
        free(grp);
        if (cnt == cap) { cap *= 2; arr = realloc(arr, cap * sizeof(Expr*)); }
        arr[cnt++] = sub;
        if (truncated) break;
    }
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), arr, cnt);
    free(arr);
    return list;
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                 */
/* ------------------------------------------------------------------ */

Expr* builtin_readlist(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    /* Options (trailing rules); values default from Options[ReadList]. */
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
    if (!options_extract(res, "ReadList", entries, 5, &new_argc)) return NULL;
    if (new_argc < 1 || new_argc > 3) return NULL;

    Expr* file_arg = res->data.function.args[0];
    if (file_arg->type != EXPR_STRING) return NULL;

    /* Type specification. */
    ReadType single = RT_EXPRESSION;
    ReadType* types = NULL;
    size_t k = 0;
    int group = 0;
    if (new_argc >= 2) {
        Expr* ta = res->data.function.args[1];
        if (is_list(ta)) {
            group = 1;
            k = ta->data.function.arg_count;
            if (k == 0) return NULL;
            types = malloc(k * sizeof(ReadType));
            for (size_t j = 0; j < k; j++) {
                if (!type_from_symbol(ta->data.function.args[j], &types[j])) {
                    free(types);
                    return NULL;
                }
            }
        } else if (!type_from_symbol(ta, &single)) {
            return NULL;
        }
    }

    /* Optional count n. */
    int have_n = 0;
    long n = 0;
    if (new_argc == 3) {
        Expr* na = res->data.function.args[2];
        if (na->type != EXPR_INTEGER || na->data.integer < 0) {
            free(types);
            return NULL;
        }
        have_n = 1;
        n = (long)na->data.integer;
    }

    /* Build separator configuration from the (possibly defaulted) options. */
    ReadCfg cfg;
    memset(&cfg, 0, sizeof cfg);
    extract_seps(o_rec, &cfg.rec, &cfg.rec_n);
    extract_seps(o_word, &cfg.word, &cfg.word_n);
    extract_seps(o_tok, &cfg.token, &cfg.token_n);
    cfg.null_records = expr_is_true(o_nrec);
    cfg.null_words = expr_is_true(o_nword);

    /* Slurp the file. */
    char* buffer = NULL;
    size_t len = 0;
    if (!slurp_file(file_arg->data.string, &buffer, &len)) {
        printf("ReadList::noopen: Cannot open %s.\n", file_arg->data.string);
        free(types);
        free_cfg(&cfg);
        return expr_new_symbol(SYM_DollarFailed);
    }

    ReadCursor c = { buffer, len, 0 };
    Expr* result = group
        ? read_groups(&c, types, k, &cfg, have_n, n)
        : read_flat(&c, single, &cfg, have_n, n);

    free(buffer);
    free(types);
    free_cfg(&cfg);
    return result;
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

void readlist_init(void) {
    symtab_add_builtin("ReadList", builtin_readlist);
    symtab_get_def("ReadList")->attributes |= ATTR_PROTECTED;

    /* Default options are the single source of truth (SetOptions works with
     * no extra code). Record separators are the actual newline byte
     * sequences; word separators are space and tab (record separators also
     * end words). */
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
    symtab_set_options("ReadList", opts);
}
