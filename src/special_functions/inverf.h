#ifndef INVERF_H
#define INVERF_H

#include "expr.h"

Expr* builtin_inverf(Expr* res);   /* InverseErf[s], InverseErf[z0, s] */
void  inverf_init(void);

#endif /* INVERF_H */
