#ifndef NAMES_H
#define NAMES_H

#include "expr.h"

/* Names[] / Names[patt] / Names[{p1, p2, ...}] -- list symbol names matching a
 * string pattern (with the * and @ metacharacters), a RegularExpression["re"],
 * or any of a list of such patterns. See src/names.c for details. */
Expr* builtin_names(Expr* res);

void names_init(void);

#endif /* NAMES_H */
