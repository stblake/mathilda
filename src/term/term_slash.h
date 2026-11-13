/*
 * term_slash.h -- the interactive REPL's slash commands: the registry, the
 * menu's candidates and the parsing of a command line.
 *
 * A line whose first character is '/' is a command, never Mathilda input
 * (no Mathilda expression can start with '/': /. // /@ /; are all infix or
 * postfix). "/name rest" runs command `name` with the argument text "rest"
 * (trimmed). Names are a letter followed by letters and digits.
 *
 * Two kinds of command share the registry:
 *   built-in   the REPL's own (/help, /theme, ...), identified by an int id
 *              the REPL dispatches on;
 *   mod        defined from Mathilda code with SlashCommand["name", f, ...];
 *              the payload is the function (an Expr*) owned by the registry,
 *              released through the free function given at registration.
 * Built-ins come first in the menu, in registration order, then mod commands
 * by name. A mod command cannot replace a built-in one.
 *
 * Pure C99, single-threaded (the main thread only).
 */
#ifndef MATHILDA_TERM_SLASH_H
#define MATHILDA_TERM_SLASH_H

#include <stdbool.h>
#include <stddef.h>

#include "term_complete.h"

typedef enum { TERM_SLASH_BUILTIN = 0, TERM_SLASH_MOD } TermSlashKind;

typedef struct {
    char*         name;      /* without the '/'                         */
    char*         summary;   /* one line for the menu; "" when none     */
    char*         owner;     /* the mod that defined it, or NULL        */
    TermSlashKind kind;
    int           id;        /* built-in: the REPL's command id         */
    void*         payload;   /* mod: the function                       */
    void        (*payload_free)(void*);
} TermSlashCmd;

/* A letter followed by letters and digits, at most 32 characters. */
bool  term_slash_valid_name(const char* s, size_t n);

/* Register a built-in command. False when the name is invalid or taken. */
bool  term_slash_add_builtin(const char* name, const char* summary, int id);

/* Register (or replace) a mod command. Takes ownership of `payload` on
 * success; on failure (invalid name, a built-in of that name, no memory) the
 * caller keeps it. `owner` may be NULL. */
bool  term_slash_add_mod(const char* name, const char* summary, const char* owner,
                         void* payload, void (*payload_free)(void*));

/* Remove a mod command; false when there is none of that name. */
bool  term_slash_remove(const char* name);

/* Remove everything (start-up tests, shut-down). */
void  term_slash_clear(void);

const TermSlashCmd* term_slash_find(const char* name, size_t n);
size_t              term_slash_count(void);
const TermSlashCmd* term_slash_at(size_t i);   /* menu order */

/* The menu's candidates for a buffer starting with '/': "/name" for every
 * command whose name starts with buf[1, cursor) (case-insensitive), in menu
 * order. word_start is 0. Nothing when the buffer holds a space before the
 * cursor or does not start with '/'. */
void  term_slash_complete(const char* buf, size_t cursor, TermCompletions* out);

/* The summary of "/name" (or "name"), or NULL. */
const char* term_slash_summary(const char* item, size_t n);

/* Split a command line: true when line[0] == '/'. *name and *nlen get the name
 * (possibly empty or invalid), *args the trimmed rest (never NULL; points
 * into `line`, its length in *alen). */
bool  term_slash_parse(const char* line, const char** name, size_t* nlen,
                       const char** args, size_t* alen);

#endif /* MATHILDA_TERM_SLASH_H */
