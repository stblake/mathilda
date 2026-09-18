/*
 * bitwise.h
 *
 * Bit-level integer builtins. Currently:
 *   - BitLength[n]   number of binary bits needed to represent the integer n.
 *
 * This is the seed of the Bit* family (BitAnd / BitOr / BitXor / BitNot /
 * BitShiftLeft / BitShiftRight, ...). Each builtin lives in its own <name>.c
 * file and is registered by bitwise_init(); docstrings live in src/info.c.
 */

#ifndef BITWISE_H
#define BITWISE_H

#include "expr.h"

Expr* builtin_bitlength(Expr* res);

void bitwise_init(void);

#endif /* BITWISE_H */
