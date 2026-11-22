/*
 * bitwise.h
 *
 * Bit-level integer builtins. Currently:
 *   - BitLength[n]          number of binary bits needed to represent n.
 *   - BitXor[n1, n2, ...]   bitwise XOR (two's complement for negatives).
 *   - BitAnd[n1, n2, ...]   bitwise AND (two's complement for negatives).
 *   - BitOr[n1, n2, ...]    bitwise OR (two's complement for negatives).
 *   - BitNot[n]             bitwise NOT, -1 - n (two's complement).
 *   - BitSet[n, k]          n with bit k (the 2^k coefficient) set to 1.
 *   - BitClear[n, k]        n with bit k (the 2^k coefficient) set to 0.
 *   - BitGet[n, k]          bit k (the 2^k coefficient) of n, as 0 or 1.
 *   - BitShiftLeft[n, k]    n * 2^k; a negative k shifts right.
 *   - BitShiftRight[n, k]   IntegerPart[n / 2^k]; a negative k shifts left.
 *
 * Each builtin lives in its own <name>.c
 * file and is registered by bitwise_init(); docstrings live in src/info.c.
 */

#ifndef BITWISE_H
#define BITWISE_H

#include "expr.h"

Expr* builtin_bitlength(Expr* res);
Expr* builtin_bitxor(Expr* res);
Expr* builtin_bitand(Expr* res);
Expr* builtin_bitor(Expr* res);
Expr* builtin_bitnot(Expr* res);
Expr* builtin_bitset(Expr* res);
Expr* builtin_bitclear(Expr* res);
Expr* builtin_bitget(Expr* res);
Expr* builtin_bitshiftleft(Expr* res);
Expr* builtin_bitshiftright(Expr* res);

/* The shared core of BitShiftLeft (dir +1) and BitShiftRight (dir -1), in
 * bitshiftleft.c. `head` names the caller in its messages. */
Expr* bitshift_eval(Expr* res, const char* head, int dir);

void bitwise_init(void);

#endif /* BITWISE_H */
