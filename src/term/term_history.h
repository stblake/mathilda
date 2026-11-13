/*
 * term_history.h -- the interactive REPL's persistent input history.
 *
 * Pure C99 (plus a guarded POSIX open() so a new file is created 0600); no
 * readline. The readline adapter feeds the loaded entries to readline's own
 * history list; a later in-tree editor can use the store directly.
 *
 * File format (one entry per line, so multi-line inputs survive):
 *
 *     #mathilda-history v1
 *     f[\n1]
 *     a = "x\\y"
 *
 * The first line is the header. Each following line is one entry with '\'
 * written as "\\", LF as "\n" and CR as "\r"; any other "\x" decodes to "\x"
 * literally. Later lines starting with '#' are entries (#^2 & is valid input).
 * A non-empty file without the header is never written to: history is then
 * kept for the session only, with one "Mathilda: ..." note on stderr.
 *
 * Every accepted entry is appended at once (open, write one line, close), so
 * a crash or a forced quit loses nothing. At load, a file holding more than
 * the size limit is rewritten with its newest entries (via <path>.tmp and
 * rename). Blank and comment-only entries, and an entry equal to the one
 * before it, are not recorded.
 */
#ifndef MATHILDA_TERM_HISTORY_H
#define MATHILDA_TERM_HISTORY_H

#include <stdbool.h>
#include <stddef.h>

#define TERM_HISTORY_DEFAULT_SIZE 1000u
#define TERM_HISTORY_MAX_SIZE     1000000u

typedef struct TermHistory TermHistory;

/* Path resolution: MATHILDA_HISTFILE (set but empty = disabled), else
 * $HOME/.mathilda_history. Returns a malloc'd path, or NULL when history is
 * disabled or HOME is unset or empty. The environment values are passed in. */
char*  term_history_resolve_path(const char* env_histfile, const char* env_home);

/* MATHILDA_HISTSIZE: unset or empty gives the default (1000); "0" disables
 * history; anything that is not a plain decimal number gives the default;
 * values above 1,000,000 are capped. */
size_t term_history_size_from_env(const char* env_histsize);

/* One entry as one file line, and back. Both return malloc'd strings. */
char*  term_history_encode(const char* entry);
char*  term_history_decode(const char* line);

/* Open a store holding at most `max` entries. With a path, load (and if
 * needed truncate) the file, creating it with the header when absent; with a
 * NULL path, an in-memory store for the session. max == 0 returns NULL. */
TermHistory* term_history_open(const char* path, size_t max);

/* Filter, remember and append one accepted entry. Returns true when the entry
 * was kept (so the caller adds it to its own history list too). */
bool   term_history_add(TermHistory* h, const char* entry);

size_t      term_history_count(const TermHistory* h);
const char* term_history_get(const TermHistory* h, size_t i);   /* 0 = oldest */
size_t      term_history_max(const TermHistory* h);
/* True when entries reach the file (a path, a Mathilda history file, no write
 * failure so far). */
bool        term_history_persistent(const TermHistory* h);

/* Free the store; NULL is a no-op. Nothing is written at close. */
void   term_history_close(TermHistory* h);

#endif /* MATHILDA_TERM_HISTORY_H */
