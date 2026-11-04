#ifndef MATHILDA_UTF8_H
#define MATHILDA_UTF8_H

/*
 * utf8.h - character (code point) view of Mathilda's UTF-8 strings.
 *
 * Mathilda stores a String as NUL-terminated UTF-8; Mathematica indexes strings
 * by character. Every String* builtin that counts, indexes, slices or reverses
 * therefore works in code points through these helpers:
 *
 *   - A "character" is one well-formed UTF-8 sequence (1-4 bytes, no overlong
 *     forms, no surrogates, <= U+10FFFF). A byte that does not start a valid
 *     sequence is one character on its own (a Latin-1-style fallback), so
 *     malformed input is never split mid-character, never crashes, and
 *     round-trips byte-for-byte.
 *   - Pure-ASCII strings (the overwhelmingly common case) are detected with a
 *     word-at-a-time scan, and callers then use byte arithmetic directly:
 *     character index == byte index costs nothing extra.
 *
 * Positions in this API are 0-based; builtins convert from Mathematica's
 * 1-based / negative positions themselves.
 */

#include <stddef.h>
#include <stdint.h>

/* True iff s[0, n) is pure 7-bit ASCII. */
int u8_is_ascii(const char* s, size_t n);

/* Byte length of the character starting at s[i] (1..4); s[i, n) nonempty. */
size_t u8_char_len(const char* s, size_t n, size_t i);

/* Number of characters in s[0, n). */
size_t u8_length(const char* s, size_t n);

/* Byte offset of character `k` (0-based) in s[0, n); k == length gives n.
 * k beyond the end also gives n. O(k) (O(1) for ASCII). */
size_t u8_offset(const char* s, size_t n, size_t k);

/*
 * Character view of a string: `nchars` characters whose byte offsets are
 * off[0..nchars] (off[nchars] == nbytes). For an ASCII string `off` is NULL
 * and the offset of character k is k itself. Build with u8_view_init, release
 * with u8_view_free. Use when a builtin needs many positions of one string.
 */
typedef struct {
    const char* s;
    size_t      nbytes, nchars;
    size_t*     off;       /* NULL when ASCII */
} U8View;

int    u8_view_init(U8View* v, const char* s, size_t n);   /* 0, or -1 on OOM */
void   u8_view_free(U8View* v);
static inline size_t u8_view_off(const U8View* v, size_t k) {
    return v->off ? v->off[k] : k;
}

/* Decode the character at s[*i] and advance *i past it. A malformed byte
 * decodes to its own value (0x80..0xFF). */
uint32_t u8_decode(const char* s, size_t n, size_t* i);

/* Encode code point `cp` into out (>= 4 bytes); returns the byte count, or 0
 * for an invalid code point (surrogate or > U+10FFFF). */
size_t u8_encode(uint32_t cp, char* out);

/* Code point of the named character \[name] (name without brackets, length
 * `len`), or -1 if unknown. Covers Mathematica's named-character table. */
long u8_named_char(const char* name, size_t len);

/* How the parser treats a non-ASCII code point in source code (outside a
 * string literal).  A code point not in the table is a LETTER: it forms
 * identifiers, as Greek letters do (a\[Mu]1 is one symbol). */
typedef enum {
    U8C_LETTER,   /* identifier character; `head` names the System symbol a
                     one-character identifier aliases (\[Pi] alone is Pi) */
    U8C_SPACE,    /* whitespace (\[InvisibleSpace], \[VeryThinSpace], ...) */
    U8C_NEWLINE,  /* a line break (\[IndentingNewLine]) */
    U8C_INFIX,    /* binary operator `head` (\[CenterDot], \[Element], ...) */
    U8C_PREFIX,   /* prefix operator `head` (\[Not], \[Del], ...) */
    U8C_POSTFIX,  /* postfix operator `head` (\[Transpose], ...) */
    U8C_OTHER     /* cannot appear in source (brackets, \[PartialD], ...) */
} U8CharClass;

typedef enum { U8A_NONE, U8A_FLAT, U8A_LEFT, U8A_RIGHT } U8Assoc;

typedef struct {
    uint32_t cp;
    U8CharClass cls;
    U8Assoc assoc;     /* INFIX: a<c>b<c>d is h[a,b,d] / h[h[a,b],d] / h[a,h[b,d]] */
    int prec;          /* operator precedence on the parser's 1..10000 ladder */
    const char* head;  /* operator head, letter alias, or NULL */
} U8CharInfo;

/* The source-code class of non-ASCII code point `cp`, or NULL when it is a
 * plain letter.  Generated from the reference parser's behaviour by
 * tools/classify_named_chars.wls + tools/gen_named_chars.py --classes. */
const U8CharInfo* u8_char_class(uint32_t cp);

/* Call fn(head) for every operator head and letter alias in the class table
 * (Wedge, CircleTimes, Element, Degree, ...). Every one of them is a System`
 * symbol the parser can emit, so sym_names_init interns them all up front:
 * otherwise a package that defines rules for, say, Wedge inside BeginPackage
 * mints a private DifferentialForms`Wedge while `a \[Wedge] b` still parses to
 * System`Wedge, and the infix form never reaches the package's rules. */
void u8_for_each_class_head(void (*fn)(const char* head));

#endif /* MATHILDA_UTF8_H */
