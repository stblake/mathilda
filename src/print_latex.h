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
 * Examples:
 *   Times[Pi, Power[E, -2]]  →  "\frac{\pi}{e^{2}}"
 *   Power[x, Rational[1,2]]  →  "\sqrt{x}"
 *   Rational[1,6]            →  "\frac{1}{6}"
 *   Sin[x]                   →  "\sin(x)"
 */
char* expr_to_latex(const Expr* e);

#endif /* MATHILDA_PRINT_LATEX_H */
