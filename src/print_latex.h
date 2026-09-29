#ifndef MATHILDA_PRINT_LATEX_H
#define MATHILDA_PRINT_LATEX_H

#include "expr.h"

/*
 * expr_to_latex(e)
 *
 * Convert an expression to a KaTeX-compatible LaTeX string in StandardForm
 * style: fractions, roots, superscripts, Greek letters, trig functions, etc.
 *
 * Returns a heap-allocated string — caller must free().
 * Returns NULL only on allocation failure.
 *
 * Returns the EMPTY string for an expression carrying a printer directive
 * (InputForm, FullForm, TeXForm, NumberForm anywhere in the tree): those ask
 * for a notation typesetting cannot express, and the plain-text form from
 * expr_to_string is already the answer. Every caller treats an empty result as
 * "no LaTeX" and falls back to that text, so a directive is honoured rather
 * than overridden by the typeset StandardForm the user asked not to see.
 *
 * Examples:
 *   Times[Pi, Power[E, -2]]  →  "\frac{\pi}{e^{2}}"
 *   Power[x, Rational[1,2]]  →  "\sqrt{x}"
 *   Rational[1,6]            →  "\frac{1}{6}"
 *   Sin[x]                   →  "\sin(x)"
 *   InputForm[Rational[1,6]] →  ""   (use expr_to_string's "1/6")
 */
char* expr_to_latex(const Expr* e);

#endif /* MATHILDA_PRINT_LATEX_H */
