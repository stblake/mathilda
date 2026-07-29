#ifndef RIFFLE_H
#define RIFFLE_H

#include "expr.h"

/* Riffle[list, x] / Riffle[list, {x1, x2, ...}]
 *
 * Interleaves separators into the gaps BETWEEN successive elements of `list`,
 * never before the first element and never after the last. A list of n
 * elements has exactly n - 1 gaps, so a list of length 0 or 1 comes back
 * unchanged. A List second argument supplies its elements cyclically, filling
 * the gaps left to right; any other expression is a single separator used in
 * every gap. See riffle.c for the full semantics. */
Expr* builtin_riffle(Expr* res);

#endif /* RIFFLE_H */
