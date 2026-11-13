/*
 * term_readline.h -- the readline backend of the line-reader seam.
 *
 * The ONLY file in Mathilda that includes <readline/...>. It is kept out of
 * libmathilda.a (makefile LIB_OBJ) and out of the CMake test build, and its
 * whole body compiles away under NO_READLINE.
 */
#ifndef MATHILDA_TERM_READLINE_H
#define MATHILDA_TERM_READLINE_H

#include <stdbool.h>
#include "term_input.h"

struct TermHistory;

typedef struct {
    /* The interactive session: smart Return, the Ctrl-C line clear, and (in
     * later phases) completion, the history file and the coloured prompt.
     * Off, it is the plain reader -file's Dialog has always used: readline()
     * plus in-memory add_history of non-empty lines, nothing else. */
    bool                interactive_features;
    struct TermHistory* history;   /* persistent history store, or NULL */
} TermReadlineOptions;

#ifndef NO_READLINE
TermReader* term_reader_open_readline(const TermReadlineOptions* o);
#endif

#endif /* MATHILDA_TERM_READLINE_H */
