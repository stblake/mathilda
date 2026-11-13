#ifndef NAMES_H
#define NAMES_H

#include "expr.h"

/* Names[] / Names[patt] / Names[{p1, p2, ...}] -- list symbol names matching a
 * string pattern (with the * and @ metacharacters), a RegularExpression["re"],
 * or any of a list of such patterns. See src/names.c for details. */
Expr* builtin_names(Expr* res);

/* Symbol-table keys matched by one Names-style pattern (see src/names.c);
 * used by Remove / ClearAll / Clear to expand "ctx`x*" arguments. */
#include <stddef.h>
char** names_match_keys(Expr* pattern, size_t* n);
void names_free_keys(char** keys, size_t n);

/* Completion candidates for `prefix`: the display names of defined symbols
 * that start with it (fully qualified names when the prefix has a backtick),
 * sorted and de-duplicated, Module temporaries and context-hidden names left
 * out. Never interns and never evaluates. Free with names_free_keys; NULL
 * with *n = 0 when nothing matches. */
char** names_complete(const char* prefix, size_t* n);

void names_init(void);

#endif /* NAMES_H */
