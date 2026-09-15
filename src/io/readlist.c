/*
 * readlist.c — ReadList[] file-reading builtin.
 *
 * ReadList["file"]                 read all remaining Expressions -> List
 * ReadList["file", type]           read `type` until EOF -> flat List
 * ReadList["file", {t1,...,tk}]    read one of each type per pass -> List of
 *                                  k-element sublists (EndOfFile-padded when
 *                                  EOF arrives partway through a pass)
 * ReadList["file", types, n]       stop after n top-level objects / passes
 * ReadList[InputStream[...], ...]  read from an already-open input stream
 *
 * ReadList is a sequence of calls to the shared single-object reader: it loops
 * read_structure() (read.c) to end of file, collecting each result.  The read
 * types, separator options, and EndOfFile padding all come from that one
 * engine, which Read[] also uses.  A filename that is not already open is opened
 * and closed by ReadList (there is no user-visible stream in that case).
 *
 * Options (trailing rules, defaults from Options[ReadList]): RecordSeparators,
 * WordSeparators, TokenWords, NullRecords, NullWords.
 *
 * ANSI C99 only; no POSIX symbols (SPEC.md §10).
 */

#include "readlist.h"
#include "read.h"
#include "streams.h"
#include "sym_names.h"
#include "symtab.h"
#include "attr.h"
#include "options.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Resolve arg 0 to a live input stream.  *opened_here is set when ReadList
 * itself opened the file (and must close it afterwards).  On a hard failure a
 * message is printed and *fail is set (return $Failed); NULL with *fail == 0
 * means "not a stream/filename form" (decline). */
static Stream* resolve_input(const Expr* a, int* opened_here, int* fail) {
    *opened_here = 0;
    *fail = 0;
    int id, is_output;
    if (stream_handle(a, &id, &is_output)) {
        if (is_output) {
            printf("ReadList::openx: an OutputStream is not open for reading.\n");
            *fail = 1;
            return NULL;
        }
        Stream* s = stream_by_id(id);
        if (!s) {
            printf("ReadList::openx: the stream is not open.\n");
            *fail = 1;
            return NULL;
        }
        return s;
    }
    const char* name = stream_filename_arg(a);
    if (!name) return NULL;                 /* decline */
    Stream* s = stream_find_by_name(name, 0);
    if (s) return s;                        /* already open: use it, leave open */
    int nid = stream_open_input(name);
    if (nid < 0) {
        printf("ReadList::noopen: Cannot open %s.\n", name);
        *fail = 1;
        return NULL;
    }
    *opened_here = 1;
    return stream_by_id(nid);
}

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

    /* Type specification (default: one Expression per pass). Validate before
     * opening the file. */
    Expr* default_spec = NULL;
    const Expr* spec;
    if (new_argc >= 2) {
        spec = res->data.function.args[1];
        if (!read_spec_valid(spec)) return NULL;
    } else {
        default_spec = expr_new_symbol(SYM_Expression);
        spec = default_spec;
    }

    /* Optional count n. */
    int have_n = 0;
    long n = 0;
    if (new_argc == 3) {
        Expr* na = res->data.function.args[2];
        if (na->type != EXPR_INTEGER || na->data.integer < 0) {
            if (default_spec) expr_free(default_spec);
            return NULL;
        }
        have_n = 1;
        n = (long)na->data.integer;
    }

    int opened_here = 0, fail = 0;
    Stream* s = resolve_input(res->data.function.args[0], &opened_here, &fail);
    if (!s) {
        if (default_spec) expr_free(default_spec);
        return fail ? expr_new_symbol(SYM_DollarFailed) : NULL;
    }

    ReadCfg cfg;
    read_cfg_build(&cfg, "ReadList", o_rec, o_word, o_tok, o_nrec, o_nword);

    /* Loop the single-object reader to EOF (or n passes), collecting results. */
    ReadCursor c = { s->buf, s->len, s->pos };
    size_t cap = 8, cnt = 0;
    Expr** arr = malloc(cap * sizeof(Expr*));
    while (!have_n || (long)cnt < n) {
        size_t prev = c.pos;
        int read_any = 0, bad_spec = 0;
        Expr* v = read_structure(&c, spec, &cfg, &read_any, &bad_spec);
        if (bad_spec) { if (v) expr_free(v); break; }
        if (!v) break;                                  /* EOF at a pass boundary */
        if (c.pos == prev) { expr_free(v); break; }     /* defensive: no progress */
        if (cnt == cap) { cap *= 2; arr = realloc(arr, cap * sizeof(Expr*)); }
        arr[cnt++] = v;
    }
    s->pos = c.pos;

    Expr* list = expr_new_function(expr_new_symbol(SYM_List), arr, cnt);
    free(arr);

    read_cfg_free(&cfg);
    if (default_spec) expr_free(default_spec);
    if (opened_here) stream_close_id(s->id);
    return list;
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

    /* Default options are the single source of truth (SetOptions works with no
     * extra code). Record separators are the actual newline byte sequences;
     * word separators are space and tab (record separators also end words). */
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
