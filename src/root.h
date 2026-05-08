/* root.h
 *
 * Mathematica-style symbolic Root and RootSum.  These are held forms
 * used by integration and polynomial-algorithm internals when no
 * closed-form algebraic expression is available.
 *
 *   Root[Function[t, p[t]], k]
 *       The k-th real root of the univariate polynomial p in t.
 *       Symbolic — the system makes no attempt to express the root
 *       as a radical.
 *
 *   RootSum[Function[t, p[t]], Function[t, body[t, ...]]]
 *       The formal sum of body[α] over the roots α of p[α] == 0.
 *       Used by NaiveLogPart in src/intrat.c when LogToReal cannot
 *       close a logarithmic part to a real elementary form.
 *
 * Differentiation is wired in src/deriv.c:
 *     D[RootSum[f1, Function[t, body]], x]
 *         -> RootSum[f1, Function[t, D[body, x]]]
 *
 * Both heads carry HoldAll + Protected so the Function arguments
 * stay structurally intact and cannot be redefined by user code.
 */
#ifndef ROOT_H
#define ROOT_H

#include "expr.h"

/* Module init: registers Root and RootSum and assigns attributes. */
void root_init(void);

/* Build a RootSum[Function[bvar, poly], Function[bvar, body]] tree.
 * Takes ownership of every argument.  The bound variable is used in
 * both Function nodes; the caller must ensure poly and body are
 * expressed in terms of bvar. */
Expr* root_make_rootsum(Expr* bvar, Expr* poly, Expr* body);

#endif
