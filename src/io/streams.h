/*
 * streams.h — the input/output stream registry and stream-management builtins.
 *
 * Mathilda's first stream layer.  A stream is a file opened for reading (a
 * whole-file buffer plus a moving read point) or for writing (an open FILE*).
 * Streams live in a process-global registry and are handed to the language as
 * inert InputStream["name", id] / OutputStream["name", id] objects; the small
 * integer id addresses the registry entry across successive Read/Write calls.
 *
 * Builtins registered here: OpenRead, OpenWrite, OpenAppend, Close, Streams,
 * StreamPosition, SetStreamPosition, Write, WriteString.  The Read[] builtin and
 * the shared reading engine live in read.c/read.h.
 */
#ifndef MATHILDA_IO_STREAMS_H
#define MATHILDA_IO_STREAMS_H

#include "expr.h"
#include <stddef.h>
#include <stdio.h>

/* A registered stream. `id == 0` marks a free slot. */
typedef struct {
    int    id;          /* stable handle (> 0 when in use) */
    int    is_output;   /* 0 = input, 1 = output */
    char*  name;        /* file name (owned) */
    /* input streams */
    char*  buf;         /* slurped file contents (owned, NUL-terminated) */
    size_t len;         /* bytes in buf */
    size_t pos;         /* current read point */
    /* output streams */
    FILE*  fp;          /* open output handle */
} Stream;

/* Open `name` for reading: slurp it and register an input stream.  Returns the
 * new id (> 0), or -1 if the file cannot be opened/read. */
int stream_open_input(const char* name);

/* Open `name` for writing (append ? "ab" : "wb") and register an output
 * stream.  Returns the new id (> 0), or -1 on failure. */
int stream_open_output(const char* name, int append);

/* Live stream with this id, or NULL. */
Stream* stream_by_id(int id);

/* Live stream with this name and kind (want_output: 0 input, 1 output, -1
 * either), most recently opened first, or NULL. */
Stream* stream_find_by_name(const char* name, int want_output);

/* Close and free the live stream with this id.  Returns 1 if it was open. */
int stream_close_id(int id);

/* Registry slot iteration (for Streams[]): capacity and slot access; a slot
 * with id == 0 is free. */
size_t stream_capacity(void);
Stream* stream_slot(size_t i);

/* If `a` is an InputStream[...] / OutputStream[...] object, fill *id / *is_output
 * from it (id may name a slot that is no longer open) and return 1; else 0. */
int stream_handle(const Expr* a, int* id, int* is_output);

/* Borrowed filename from a bare String or File["name"]; NULL if `a` is neither. */
const char* stream_filename_arg(const Expr* a);

/* Build the inert InputStream["name", id] / OutputStream["name", id] object. */
Expr* stream_make_object(const Stream* s);

/* Register the stream-management builtins, the inert InputStream/OutputStream/
 * File heads, and their attributes. */
void streams_init(void);

#endif /* MATHILDA_IO_STREAMS_H */
