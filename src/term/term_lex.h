/*
 * term_lex.h -- the interactive terminal's lightweight input lexer.
 *
 * Pure C99, no readline: shared by the readline adapter (smart Return,
 * Tab completion) and, later, the in-tree line editor. It answers two
 * questions about a partially typed input buffer:
 *
 *   - Is it a complete expression, so Return should submit it?
 *     (term_lex_input_complete)
 *   - What surrounds the cursor: a string, a comment, the word being typed,
 *     the head of the innermost enclosing call, and whether the cursor is in
 *     that call's argument position? (term_lex_context, for completion)
 *
 * Both follow the real lexer's string and comment rules: a '"' inside a
 * comment and a "(*" inside a string are plain text, comments nest, and a
 * backslash escapes the next character inside a string.
 */
#ifndef MATHILDA_TERM_LEX_H
#define MATHILDA_TERM_LEX_H

#include <stdbool.h>
#include <stddef.h>

/* True when `s` is a syntactically complete Mathilda input: no unterminated
 * "..." string or (* ... *) comment, and every '(', '[' or '{' has a matching
 * closer. An over-closed buffer counts as complete (appending text cannot
 * repair it, so it should submit and let the parser report the error); the
 * empty string is complete. */
bool term_lex_input_complete(const char* s);

/* What surrounds the cursor in a partially typed buffer. Offsets are byte
 * offsets into the buffer passed to term_lex_context. */
typedef struct {
    bool   in_string, in_comment;
    size_t string_start;        /* offset of the opening '"' when in_string        */
    size_t word_start;          /* start of the [A-Za-z0-9$`]* run ending at cursor */
    bool   has_call_head;       /* innermost unclosed '[' is a call  f[  (not [[ )  */
    size_t head_start, head_len;/* the symbol token before that '['                 */
    bool   arg_position;        /* last significant char at that depth is '[' or ','
                                   (for a string: the char before its opening '"') */
} TermLexContext;

/* Fill `out` for the first `cursor` bytes of `buf` (cursor is clamped to the
 * buffer length). Never allocates. */
void term_lex_context(const char* buf, size_t cursor, TermLexContext* out);

/* True for the characters a completion word is made of: [A-Za-z0-9$`]. */
bool term_lex_is_word_char(char c);

#endif /* MATHILDA_TERM_LEX_H */
