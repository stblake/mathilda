/*
 * findmin.h
 *
 * Numerical local optimization: FindMinimum / FindMaximum.
 * Mirror of findroot.h. See findmin.c for the full surface and option set.
 */
#ifndef FINDMIN_H
#define FINDMIN_H

#include "expr.h"

Expr* builtin_findminimum(Expr* res);
Expr* builtin_findmaximum(Expr* res);
void  findmin_init(void);

#endif /* FINDMIN_H */
