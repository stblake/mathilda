/* Mathilda — Precision / Accuracy introspection and coercion.
 *
 * Exposes four builtins:
 *   Precision[x]        — returns Infinity, MachinePrecision, or decimal digits.
 *   Accuracy[x]         — same family, using digits-after-decimal semantics.
 *   SetPrecision[x, n]  — produces an equivalent expression at precision n.
 *   SetAccuracy[x, n]   — same with accuracy n.
 *
 * All four are ATTR_PROTECTED | ATTR_LISTABLE. SetPrecision/SetAccuracy
 * walk the expression tree similarly to N[], so operating on compound
 * expressions (e.g. `SetPrecision[Sin[1] + 1/3, 30]`) does the right thing.
 *
 * Without USE_MPFR, the 2-arg Set* forms emit a one-shot warning and
 * fall back to machine precision, mirroring N[…, prec]'s behavior.
 */
#ifndef PRECISION_H
#define PRECISION_H

#include "expr.h"

Expr* builtin_precision(Expr* res);
Expr* builtin_accuracy(Expr* res);
Expr* builtin_set_precision(Expr* res);
Expr* builtin_set_accuracy(Expr* res);

void precision_init(void);

#endif /* PRECISION_H */
