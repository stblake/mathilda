#ifndef READWRITE_H
#define READWRITE_H

#include "expr.h"

/* File I/O builtins.
 *
 *   Get["file"]                           reads and evaluates every
 *                                         expression in "file", returning
 *                                         the value of the last one.
 *
 *   Put[expr_1, ..., expr_n, "file"]      writes the InputForm of each
 *                                         expr_i to "file", one per line,
 *                                         replacing any existing content.
 *   Put["file"]                           creates an empty "file".
 *   expr >> "file"                        parser shorthand for
 *                                         Put[expr, "file"].
 *
 *   PutAppend[expr_1, ..., expr_n, "file"]
 *                                         like Put but appends to the
 *                                         existing contents of "file".
 *   expr >>> "file"                       parser shorthand for
 *                                         PutAppend[expr, "file"].
 *
 * All three builtins return Null on success and $Failed on I/O error.
 */

void readwrite_init(void);

Expr* builtin_get(Expr* res);
Expr* builtin_put(Expr* res);
Expr* builtin_putappend(Expr* res);

#endif /* READWRITE_H */
