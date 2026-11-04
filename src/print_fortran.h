/*
 * print_fortran.h -- FortranForm and CForm text renderers.
 *
 * FortranForm[expr] / CForm[expr] print an expression as Fortran / C source
 * text: x**2 / Power(x,2), f(x,y), (0,1) / Complex(0,1), a.eq.b / a == b, ...
 * The renderer first builds a small layout tree whose rows mirror the
 * reference's box nesting (a parenthesised group is its own row, a function
 * call is one flat row of head, brackets and arguments), then lays it out.
 *
 * With a finite page width (an OpenWrite[..., FormatType -> FortranForm,
 * PageWidth -> n] stream, or ToString[e, FormatType -> FortranForm,
 * PageWidth -> n]) an expression that does not fit on one line is written as
 * fixed-form Fortran: an 8-space first line and continuation lines
 * "     -" followed by an indent of one more space than the depth of the row
 * being broken.  Sums break after their + / - operators, argument lists
 * after their commas, and an over-long atom is split with a backslash.
 * CForm breaks the same way with no first-line prefix and a continuation
 * indent of depth + 2 spaces.
 *
 * The break points are chosen greedily (a sum term moves to the next line
 * when it does not fit and would fit there; an argument when it does not
 * fit).  That reproduces the reference on sums at any depth -- the shape
 * generated code has -- but the reference optimises over the whole layout:
 * at narrow widths it may break before a parenthesised factor
 * (m**8* / (v1 + ...)) or keep a call on the line after a comma.
 *
 * Format[f[...], FortranForm] := rhs (and CForm) display rules apply; plain
 * Format rules do not, as in the reference.
 */
#ifndef MATHILDA_PRINT_FORTRAN_H
#define MATHILDA_PRINT_FORTRAN_H

#include <stddef.h>
#include "expr.h"

typedef enum { FORM_LANG_FORTRAN = 0, FORM_LANG_C = 1 } FormLang;

/* The text of e in the given language.  page_width <= 0 means unbounded (one
 * line, no prefix).  Caller frees.  Never NULL. */
char* fortran_form_string(const Expr* e, FormLang lang, int page_width);

/* As fortran_form_string, for the n expressions written by one
 * Write[stream, e1, ..., en] call: they form one row, printed back to back. */
char* fortran_form_string_seq(Expr* const* es, size_t n, FormLang lang, int page_width);

#endif /* MATHILDA_PRINT_FORTRAN_H */
