/*
 * read.h — the shared file-reading engine and the Read[] builtin (src/io).
 *
 * Read[] and ReadList[] both read typed objects from a stream through ONE
 * engine defined here: a ReadCursor over a byte buffer, a per-call ReadCfg of
 * separator settings, and read_structure(), which reads exactly one object or
 * one (possibly nested) type-structure and advances the cursor.  Read[] calls
 * it once; ReadList[] loops it to end of file.  See readlist.c and read.c.
 */
#ifndef MATHILDA_IO_READ_H
#define MATHILDA_IO_READ_H

#include "expr.h"
#include <stddef.h>

/* One separator string, borrowed from an option Expr (valid for the call). */
typedef struct { const char* s; size_t n; } SepStr;

/* Separator / tokenisation settings for one read call. The SepStr arrays are
 * malloc'd by read_cfg_build and freed by read_cfg_free; the char* they point
 * at are borrowed from the option Exprs and must outlive the read. */
typedef struct {
    SepStr* rec;   size_t rec_n;    /* RecordSeparators */
    SepStr* word;  size_t word_n;   /* WordSeparators   */
    SepStr* token; size_t token_n;  /* TokenWords       */
    int null_records;
    int null_words;
    const char* msg_head;           /* "Read" / "ReadList" for messages */
} ReadCfg;

/* A read position over a fixed byte buffer. */
typedef struct { const char* buf; size_t len; size_t pos; } ReadCursor;

/* Build a ReadCfg from option value Exprs (any of which may be NULL to leave
 * that field empty). msg_head names the caller for diagnostics. */
void read_cfg_build(ReadCfg* cfg, const char* msg_head,
                    const Expr* rec, const Expr* word, const Expr* tok,
                    const Expr* nullrec, const Expr* nullword);

/* Free the arrays owned by a ReadCfg (not the borrowed strings). */
void read_cfg_free(ReadCfg* cfg);

/* Read ONE object or nested type-structure from the cursor.
 *
 * `spec` is either a read-type symbol (Byte, Character, Expression, Number,
 * Real, Record, String, Word) or a function whose leaves are read-type symbols;
 * in the latter case the result mirrors the spec's shape (List, Hold, any head),
 * filled depth-first.
 *
 * Returns a new Expr (which the caller/evaluator owns), or NULL for a bare
 * end-of-file — i.e. EOF reached before any leaf could be read.  Once any leaf
 * has been read, later leaves at EOF are filled with EndOfFile instead.
 * `*read_any` is set to 1 as soon as a leaf reads a value.  If `spec` contains a
 * leaf that is not a valid read type, `*bad_spec` is set to 1 and NULL returned.
 *
 * The Expression leaf is returned UNEVALUATED; the surrounding evaluator
 * evaluates the constructed result, so Read[s, Hold[Expression]] stays held
 * while Read[s, Expression] evaluates.  Number/Real leaves are validated
 * (and folded) internally.
 */
Expr* read_structure(ReadCursor* c, const Expr* spec, const ReadCfg* cfg,
                     int* read_any, int* bad_spec);

/* 1 iff `spec` is a valid type structure: every leaf is a read-type symbol and
 * there is at least one leaf (so an empty {} is rejected — it would read nothing
 * and loop forever in ReadList). */
int read_spec_valid(const Expr* spec);

/* Read[stream|"file"|File["file"], spec?] with trailing separator options. */
Expr* builtin_read(Expr* res);

/* Registers Read, its Protected attribute, and Options[Read]. */
void read_init(void);

#endif /* MATHILDA_IO_READ_H */
