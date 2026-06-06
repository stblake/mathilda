#ifndef EXPAND_POWER_H
#define EXPAND_POWER_H

#include "expr.h"

/*
 * PowerExpand[expr]                 expand powers of products and nested
 *                                   powers, and logs of products/powers.
 * PowerExpand[expr, {x1, x2, ...}]  expand only with respect to the listed
 *                                   variables.
 * PowerExpand[expr, Assumptions->a] generate branch-correct results under
 *                                   the assumption a (Assumptions->True gives
 *                                   the universally-correct formula).
 *
 * See src/expand_power.c for the transformation rules and correction-term
 * formulas.
 */
Expr* builtin_powerexpand(Expr* res);
void  expand_power_init(void);

#endif /* EXPAND_POWER_H */
