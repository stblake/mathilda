#ifndef SOWREAP_H
#define SOWREAP_H

#include "expr.h"

/* Sow / Reap — Mathematica's dynamic-scope result-accumulation mechanism.
 *
 * Reap[expr] evaluates expr and returns {value, collected}, where every
 * Sow[e] executed during that evaluation contributes e. Sow routes each
 * value to the innermost enclosing Reap whose pattern(s) match its tag.
 *
 * Reap has attribute HoldFirst (it drives evaluation of its body itself so
 * the intervening Sow side effects are captured); Sow is non-held. Both are
 * Protected. See src/core.c for registration and docstrings, and the design
 * notes in sowreap.c for the collector-stack data structures. */

Expr* builtin_sow(Expr* res);
Expr* builtin_reap(Expr* res);

#endif /* SOWREAP_H */
