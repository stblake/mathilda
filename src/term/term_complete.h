/*
 * term_complete.h -- Tab completion candidates for the interactive REPL.
 *
 * Pure C99, no readline: the readline adapter turns the candidates into
 * readline's match list; a later in-tree editor can use them directly.
 *
 * Rules (word characters are [A-Za-z0-9$`]):
 *   - nothing inside a comment;
 *   - inside a string, nothing -- unless the string opens an argument of a
 *     call whose head has string-named options; the candidates are then
 *     Name" -> ;
 *   - in argument position (just after the '[' of a call, or a ',' at that
 *     depth), the head's symbol options Name ->  come first, in their declared
 *     order, then symbols (a symbol an option already offers is dropped);
 *   - symbols: defined ones whose display name starts with the word (fully
 *     qualified names when the word has a backtick), sorted, with Module
 *     temporaries left out (names_complete in src/names.c);
 *   - an empty word completes nothing outside argument position, and only
 *     option names inside it: a bare Tab never lists thousands of names.
 * It never evaluates, never interns a symbol and makes no dependency note.
 */
#ifndef MATHILDA_TERM_COMPLETE_H
#define MATHILDA_TERM_COMPLETE_H

#include <stddef.h>

typedef struct {
    char** items;        /* malloc'd strings, the full text replacing the word */
    size_t count;
    size_t word_start;   /* byte offset of the word being completed            */
} TermCompletions;

/* Candidates for the word ending at `cursor` in `buf`. Never evaluates. */
void term_complete(const char* buf, size_t cursor, TermCompletions* out);
void term_completions_free(TermCompletions* c);

#endif /* MATHILDA_TERM_COMPLETE_H */
