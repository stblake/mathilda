#ifndef PICOSTRINGS_H
#define PICOSTRINGS_H

#include "expr.h"

/*
 * String manipulation builtins for Mathilda. Each builtin lives in its own
 * translation unit under src/strings/; strings_init (strings_init.c) registers
 * them all with their attributes and docstrings.
 */
Expr* builtin_stringlength(Expr* res);
Expr* builtin_characters(Expr* res);
Expr* builtin_stringjoin(Expr* res);
Expr* builtin_stringpart(Expr* res);
Expr* builtin_stringtake(Expr* res);
Expr* builtin_stringdrop(Expr* res);
Expr* builtin_stringreverse(Expr* res);
Expr* builtin_stringinsert(Expr* res);
Expr* builtin_stringreplacepart(Expr* res);
Expr* builtin_stringpartition(Expr* res);

/*
 * Regular-expression string builtins (src/strings/regex/, backed by PCRE2).
 * Registered by regex_init(), which core_init() calls after strings_init().
 */
Expr* builtin_regularexpression(Expr* res);
Expr* builtin_stringmatchq(Expr* res);
Expr* builtin_stringcases(Expr* res);
Expr* builtin_stringreplace(Expr* res);
Expr* builtin_stringsplit(Expr* res);
Expr* builtin_stringtrim(Expr* res);
Expr* builtin_stringextract(Expr* res);

void strings_init(void);
void regex_init(void);

#endif // PICOSTRINGS_H
