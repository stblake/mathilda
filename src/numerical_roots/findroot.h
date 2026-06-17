/*
 * findroot.h
 *
 * Numerical root finding: FindRoot[f, {x, x0}] and related forms.
 * See findroot.c for the full surface and option set.
 */
#ifndef FINDROOT_H
#define FINDROOT_H

#include "expr.h"

Expr* builtin_findroot(Expr* res);
void  findroot_init(void);

#endif /* FINDROOT_H */
