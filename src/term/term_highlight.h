/*
 * term_highlight.h -- incremental syntax colouring for the in-tree editor.
 *
 * Pure C99. A line-at-a-time tokenizer: term_hl_line() colours one line given
 * the lexer state at its start (inside a "..." string, inside (* ... *)
 * comments and how deeply) and returns the state at its end, so the editor
 * relexes only from an edited line until the carried state stops changing.
 * The string and comment rules are term_lex's (and the real lexer's): a '"'
 * inside a comment and a "(*" inside a string are plain text, comments nest,
 * a backslash escapes the next character inside a string.
 *
 * The style ids double as the editor's render styles: the lexer classes come
 * first, then the editor's own (bracket match, popup, ...). The SGR table is
 * private to this module so term_color.{c,h} stay untouched.
 */
#ifndef MATHILDA_TERM_HIGHLIGHT_H
#define MATHILDA_TERM_HIGHLIGHT_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

typedef enum {
    TERM_HL_PLAIN = 0,
    TERM_HL_BUILTIN,        /* a System` / Mathilda` symbol            */
    TERM_HL_USER,           /* a symbol with a value or rules          */
    TERM_HL_NUMBER,
    TERM_HL_STRING,
    TERM_HL_COMMENT,
    TERM_HL_PATTERN,        /* x_  _h  x__  x___  x_.  #  #2  ##        */
    TERM_HL_BRACKET,        /* ( ) [ ] { } outside strings and comments */
    TERM_HL_LEX_COUNT,      /* -- editor styles below ------------------ */
    TERM_HL_BRACKET_MATCH = TERM_HL_LEX_COUNT,
    TERM_HL_BRACKET_BAD,
    TERM_HL_PROMPT,
    TERM_HL_POPUP,
    TERM_HL_POPUP_SEL,
    TERM_HL_POPUP_DOC,
    TERM_HL_DIM,
    TERM_HL_COUNT
} TermHlStyle;

typedef struct {
    uint8_t  in_string;     /* inside "..."                 */
    uint8_t  invalid;       /* editor bookkeeping: not lexed */
    uint16_t comment;       /* (* ... *) nesting depth       */
} TermHlState;

/* Classify an identifier token s[0, len): TERM_HL_BUILTIN, TERM_HL_USER or
 * TERM_HL_PLAIN. */
typedef int (*TermHlClassify)(const char* s, size_t len, void* ctx);

/* The default classifier: symtab_lookup (which never creates a symbol);
 * is_system -> builtin; values, rules or a C builtin -> user. */
int term_hl_classify_symtab(const char* s, size_t len, void* ctx);

/* The docstring of symbol s[0, len) (from symtab_lookup: never creates a
 * symbol), leading blank lines skipped, or NULL. The editor shows its first
 * line as the completion popup's description. */
const char* term_hl_doc_symtab(const char* s, size_t len, void* ctx);

static inline TermHlState term_hl_state_initial(void) {
    TermHlState st = { 0, 0, 0 };
    return st;
}
static inline bool term_hl_state_eq(TermHlState a, TermHlState b) {
    return a.in_string == b.in_string && a.comment == b.comment && !a.invalid && !b.invalid;
}

/* Colour one line s[0, n) (no '\n' inside) into out[0, n), one style id per
 * byte, starting in state `st`; return the state after the line (a line break
 * does not end a string or a comment). cls may be NULL (every symbol plain).
 * Never allocates. */
TermHlState term_hl_line(const char* s, size_t n, TermHlState st, uint8_t* out,
                         TermHlClassify cls, void* ctx);

/* The SGR sequence that selects style `id` from a reset state ("" for
 * TERM_HL_PLAIN and for an id out of range). */
const char* term_hl_sgr(int id);

/* Write s[0, n) to f with SGR colouring (every line lexed, styles reset at
 * the end). Returns false, writing nothing, when it cannot allocate. */
bool term_hl_fputs(const char* s, size_t n, FILE* f, TermHlClassify cls, void* ctx);

/* Print a result as expr_print would, with syntax colouring: the same bytes
 * plus SGR sequences, then a newline. Returns false, having printed nothing,
 * when the text is over 64 KiB or cannot be built; the caller then uses
 * expr_print. The interactive editor session's results only. */
struct Expr;
bool term_hl_print_result(struct Expr* e);

#endif /* MATHILDA_TERM_HIGHLIGHT_H */
