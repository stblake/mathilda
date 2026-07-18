/*
 * complex_expand.h
 *
 * ComplexExpand[expr] and ComplexExpand[expr, {vars}].
 *
 * Expands an expression into explicit real and imaginary parts, assuming
 * every free symbol is real unless it matches one of the (possibly
 * pattern) variables named in the optional second argument.  Propagates a
 * (real, imaginary) decomposition recursively through Plus, Times, Power,
 * Exp, Log, all circular and hyperbolic functions and their inverses, and
 * the Re/Im/Abs/Arg/Conjugate/Sign/ReIm wrapper heads.
 *
 * The output basis is controlled by the option
 *   TargetFunctions -> {Re, Im}          (default; Cartesian form)
 *   TargetFunctions -> {Abs, Arg}        (polar form)
 *   TargetFunctions -> Conjugate         (conjugate form)
 *
 * ComplexExpand threads over List, and over equations, inequalities and
 * logic heads (Equal, Less, And, Or, ...).
 */
#ifndef COMPLEX_EXPAND_H
#define COMPLEX_EXPAND_H

#include "expr.h"

/* Builtin entry point.  Owns `res`; returns a new Expr on success or NULL
 * to leave the call unevaluated (bad argument count). */
Expr* builtin_complex_expand(Expr* res);

/* Register ComplexExpand (attributes + builtin). */
void complex_expand_init(void);

#endif /* COMPLEX_EXPAND_H */
