/*
 * term_input.h -- the interactive REPL's line-reader seam.
 *
 * The REPL loop and Dialog read their lines through a TermReader, never
 * through readline or fgets directly, so the backend can be swapped: the
 * stdio backend here (the NO_READLINE fgets loop), the readline adapter
 * (term_readline.h), and later an in-tree editor. Pure C99, no readline.
 */
#ifndef MATHILDA_TERM_INPUT_H
#define MATHILDA_TERM_INPUT_H

#include <stdbool.h>

typedef enum { TERM_READ_OK, TERM_READ_EOF, TERM_READ_INTERRUPTED } TermReadStatus;
typedef enum { TERM_BACKEND_STDIO, TERM_BACKEND_READLINE } TermBackendKind;

typedef struct TermReader TermReader;

/* Backend vtable. `read` shows `prompt` and stores a malloc'd line (no
 * trailing newline) in *out on TERM_READ_OK; *out is NULL otherwise. */
typedef struct {
    TermReadStatus (*read)(TermReader*, const char* prompt, char** out);
    void           (*history_add)(TermReader*, const char* entry);
    void           (*close)(TermReader*);
} TermReaderOps;

struct TermReader {
    const TermReaderOps* ops;
    TermBackendKind      kind;
    void*                impl;
};

/* Lines from stdin with fgets, exactly as the NO_READLINE loop always read
 * them: the prompt on stdout, trailing CR/LF stripped, no history. */
TermReader*     term_reader_open_stdio(void);

/* Read one line. The interrupt state is READING for the duration of the call
 * (a Ctrl-C clears the line instead of aborting anything) and then restored,
 * so a Dialog reader inside an evaluation returns to EVALUATING. */
TermReadStatus  term_read(TermReader* rd, const char* prompt, char** out);
/* Record an accepted entry in the backend's history (a no-op for stdio). */
void            term_reader_history_add(TermReader* rd, const char* entry);
TermBackendKind term_reader_kind(const TermReader* rd);
/* Close and free the reader; NULL is a no-op. */
void            term_reader_close(TermReader* rd);

#endif /* MATHILDA_TERM_INPUT_H */
