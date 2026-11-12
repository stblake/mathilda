/*
 * term_width.h -- display width of text for the in-tree line editor.
 *
 * Pure C99, locale-independent: a terminal column count per code point from
 * a compact table (East Asian Wide/Fullwidth ranges and emoji are 2, combining
 * marks, joiners and variation selectors are 0), so the editor knows where its
 * cursor is without asking the C library's locale-dependent wcwidth().
 * Control characters (C0, DEL, C1) are drawn by the editor as ^X, so
 * term_cp_width reports them as -1 and leaves the choice to the caller.
 */
#ifndef MATHILDA_TERM_WIDTH_H
#define MATHILDA_TERM_WIDTH_H

#include <stddef.h>
#include <stdint.h>

/* Columns a code point occupies: 0, 1 or 2; -1 for a control character. */
int    term_cp_width(uint32_t cp);

/* Decode the UTF-8 character at s[i] (i < n), store its code point in *cp and
 * return the offset just past it. A byte that does not start a well-formed
 * sequence is one character on its own, code point = the byte value. */
size_t term_utf8_next(const char* s, size_t n, size_t i, uint32_t* cp);

/* Offset of the start of the character before offset i (i > 0): steps back
 * over continuation bytes, at most three. */
size_t term_utf8_prev(const char* s, size_t i);

/* Display columns of s[0, n): controls count 2 (^X), zero-width count 0. */
size_t term_text_width(const char* s, size_t n);

#endif /* MATHILDA_TERM_WIDTH_H */
