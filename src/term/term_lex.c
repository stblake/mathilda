/*
 * term_lex.c -- see term_lex.h. Pure C99; no allocation, no readline.
 */
#include "term_lex.h"

#include <string.h>

/* True when `s` is a syntactically complete Mathilda input: no unterminated
 * "..." string or (* ... *) comment, and every '(', '[' or '{' has a matching
 * closer. Newlines are whitespace to the lexer, so a balanced multi-line
 * buffer parses exactly as its single-line spelling. An *over*-closed buffer
 * (a stray ')') reports complete on purpose: appending text cannot repair it,
 * so it should submit now and let the parser flag the error rather than trap
 * the user on an endlessly growing line. The empty string is complete —
 * submitting it is a harmless no-op the caller skips. This is the predicate
 * behind the smart Return key: complete -> evaluate, incomplete -> open a
 * fresh continuation line.
 *
 * The string/comment lexing here mirrors find_unterminated() exactly (a '"'
 * inside a comment and a "(*" inside a string are both just text, and a
 * backslash escapes the next character inside a string) so the completeness
 * check never disagrees with the real lexer on where a token ends. */
bool term_lex_input_complete(const char* s) {
    int depth = 0;      /* net (), [], {} nesting outside strings/comments   */
    int comment = 0;    /* (* ... *) nesting depth                           */
    int in_string = 0;  /* inside a "..." literal                            */
    for (const char* p = s; *p; ) {
        if (in_string) {
            if (*p == '\\' && p[1]) { p += 2; continue; }
            if (*p == '"') in_string = 0;
            p++;
        } else if (comment > 0) {
            if (p[0] == '(' && p[1] == '*') { comment++; p += 2; }
            else if (p[0] == '*' && p[1] == ')') { comment--; p += 2; }
            else p++;
        } else if (p[0] == '(' && p[1] == '*') {
            comment++; p += 2;
        } else if (*p == '"') {
            in_string = 1; p++;
        } else {
            if (*p == '(' || *p == '[' || *p == '{') depth++;
            else if (*p == ')' || *p == ']' || *p == '}') depth--;
            p++;
        }
    }
    return !in_string && comment == 0 && depth <= 0;
}

bool term_lex_is_word_char(char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')
        || (c >= '0' && c <= '9') || c == '$' || c == '`';
}

static bool is_space(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v';
}

/* ---- context scanner ----------------------------------------------------- */

#define TERM_LEX_MAX_FRAMES 64

typedef enum { FRAME_PAREN, FRAME_BRACE, FRAME_CALL, FRAME_BRACKET, FRAME_PART } FrameKind;

typedef struct {
    FrameKind kind;
    size_t    head_start, head_len;   /* FRAME_CALL only                         */
    char      last_sig;               /* last significant char seen at this depth */
} Frame;

/* The word run ending at `cursor` never contains a bracket, quote or comment
 * delimiter, so the lexer state at its start equals the state at the cursor.
 * Scanning only the text before the word therefore yields both the string /
 * comment state and the "previous significant character" without the word's
 * own letters getting in the way. */
void term_lex_context(const char* buf, size_t cursor, TermLexContext* out) {
    memset(out, 0, sizeof(*out));
    if (!buf) return;
    size_t n = strlen(buf);
    if (cursor > n) cursor = n;

    size_t word_start = cursor;
    while (word_start > 0 && term_lex_is_word_char(buf[word_start - 1])) word_start--;
    out->word_start = word_start;

    Frame  frames[TERM_LEX_MAX_FRAMES];
    size_t depth = 0;                 /* may exceed TERM_LEX_MAX_FRAMES          */
    int    comment = 0;
    bool   in_string = false;
    size_t string_start = 0;
    char   sig_before_string = 0;     /* last_sig of the frame when '"' opened    */
    /* The most recent symbol token, valid while only whitespace or comments
     * have followed it. */
    bool   sym_valid = false;
    size_t sym_start = 0, sym_len = 0;
    char   prev_sig = 0;              /* previous significant char (any depth)    */

    size_t i = 0;
    while (i < word_start) {
        char c = buf[i];
        if (in_string) {
            if (c == '\\' && i + 1 < word_start) { i += 2; continue; }
            if (c == '"') in_string = false;
            i++;
            continue;
        }
        if (comment > 0) {
            if (c == '(' && i + 1 < word_start && buf[i + 1] == '*') { comment++; i += 2; }
            else if (c == '*' && i + 1 < word_start && buf[i + 1] == ')') { comment--; i += 2; }
            else i++;
            continue;
        }
        if (c == '(' && i + 1 < word_start && buf[i + 1] == '*') { comment++; i += 2; continue; }
        if (is_space(c)) { i++; continue; }

        if (term_lex_is_word_char(c)) {
            size_t s = i;
            while (i < word_start && term_lex_is_word_char(buf[i])) i++;
            bool symbol = !(buf[s] >= '0' && buf[s] <= '9');
            sym_valid = symbol;
            sym_start = s; sym_len = i - s;
            if (depth > 0 && depth <= TERM_LEX_MAX_FRAMES) frames[depth - 1].last_sig = buf[i - 1];
            prev_sig = buf[i - 1];
            continue;
        }

        /* Any other significant character. */
        if (c == '"') {
            in_string = true;
            string_start = i;
            sig_before_string = (depth > 0 && depth <= TERM_LEX_MAX_FRAMES)
                              ? frames[depth - 1].last_sig : 0;
            if (depth > 0 && depth <= TERM_LEX_MAX_FRAMES) frames[depth - 1].last_sig = '"';
        } else if (c == '(' || c == '{' || c == '[') {
            if (depth < TERM_LEX_MAX_FRAMES) {
                Frame* f = &frames[depth];
                f->head_start = f->head_len = 0;
                f->last_sig = c;
                if (c == '(') f->kind = FRAME_PAREN;
                else if (c == '{') f->kind = FRAME_BRACE;
                else if (prev_sig == '[' && i > 0 && buf[i - 1] == '['
                         && depth > 0 && frames[depth - 1].last_sig == '['
                         && frames[depth - 1].kind != FRAME_PART) {
                    /* "[[" -- the Part form; it and its partner are not calls. */
                    frames[depth - 1].kind = FRAME_PART;
                    f->kind = FRAME_PART;
                } else if (sym_valid) {
                    f->kind = FRAME_CALL;
                    f->head_start = sym_start;
                    f->head_len = sym_len;
                } else {
                    f->kind = FRAME_BRACKET;
                }
            }
            depth++;
        } else if (c == ')' || c == '}' || c == ']') {
            if (depth > 0) {
                depth--;
                if (depth > 0 && depth <= TERM_LEX_MAX_FRAMES) frames[depth - 1].last_sig = c;
            }
        } else {
            if (depth > 0 && depth <= TERM_LEX_MAX_FRAMES) frames[depth - 1].last_sig = c;
        }
        sym_valid = false;
        prev_sig = c;
        i++;
    }

    out->in_string = in_string;
    out->in_comment = comment > 0;
    if (in_string) out->string_start = string_start;
    if (depth == 0 || depth > TERM_LEX_MAX_FRAMES) return;
    const Frame* top = &frames[depth - 1];
    if (top->kind != FRAME_CALL) return;
    out->has_call_head = true;
    out->head_start = top->head_start;
    out->head_len = top->head_len;
    char sig = in_string ? sig_before_string : top->last_sig;
    out->arg_position = (sig == '[' || sig == ',');
}
