#ifndef PICOSTRINGS_H
#define PICOSTRINGS_H

#include "expr.h"

Expr* builtin_stringlength(Expr* res);
Expr* builtin_characters(Expr* res);
Expr* builtin_stringjoin(Expr* res);
Expr* builtin_stringpart(Expr* res);
Expr* builtin_stringtake(Expr* res);
Expr* builtin_stringdrop(Expr* res);
Expr* builtin_stringreverse(Expr* res);
Expr* builtin_stringinsert(Expr* res);
Expr* builtin_stringreplacepart(Expr* res);

void strings_init(void);

#endif // PICOSTRINGS_H
