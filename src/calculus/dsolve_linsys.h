/*
 * dsolve_linsys.h — shared fundamental-matrix machinery for linear ODE systems.
 *
 * Factored out of dsolve_linsys.c so that both the constant-coefficient method
 * (DSolve`LinearFirstOrderSystem) and the variable-coefficient method
 * (DSolve`LinearSystemVarCoeff) share one Jordan -> e^{Mt} -> variation-of-
 * parameters -> realify pipeline.  Only the two system-solver entry points are
 * intended for dsolve.c; the helpers are for sibling method files.
 */
#ifndef MATHILDA_DSOLVE_LINSYS_H
#define MATHILDA_DSOLVE_LINSYS_H

#include "../expr.h"
#include "dsolve_common.h"
#include <stddef.h>
#include <stdbool.h>

/* Realify + tidy a solution body; consumes `body`, returns owned.
 * = Simplify[ComplexExpand[body]] //. {Cosh[a] + Sinh[a] :> E^a}. */
Expr* dsolve_linsys_tidy(Expr* body);

/* e^{M t} for a CONSTANT matrix M given by its Jordan factors S, J (M = S J S^{-1});
 * S, J, t are borrowed and the returned n x n matrix (List of Lists) is owned.
 * `t` may be any expression (a plain symbol x for the constant method, or an
 * antiderivative Integrate[f, x] for the scalar-factor variable-coefficient one). */
Expr* dsolve_linsys_matexp(Expr* S, Expr* J, Expr* t, size_t n);

/* Extract the explicit first-order-linear structure  Y' == A Y + b  from the parsed
 * square system.  Returns false on a STRUCTURAL mismatch only (neq != nfun, some
 * order != 1, two derivatives in one equation, a non-constant / Y-dependent
 * leading-derivative coefficient, two equations for one function, or an equation not
 * solvable for its leading derivative).  On true, *A (n x n List of Lists) and *b
 * (n-vector List) are owned and MAY still depend on x (no constant-A guard here —
 * that is each caller's decision); *b_zero is set when b == 0.  n = P->nfun. */
bool dsolve_linsys_extract_Ab(DSolveProblem* P, Expr** A, Expr** b, bool* b_zero);

/* General solution of a system whose fundamental matrix is Phi = e^{M t} for the
 * CONSTANT matrix M and exponent argument `t`:
 *     Y = tidy( Phi . (C + Integrate[Phi^{-1} b, x]) ).
 * Constant-coefficient case: t = the symbol x.  Scalar-factor variable-coefficient
 * case: t = Integrate[f, x] (so Phi = e^{M Integrate[f,x]}).  `M`, `t`, `b` borrowed;
 * the variation-of-parameters integral is taken over `xvar`.  Returns a malloc'd
 * array of n owned bodies, or NULL if Jordan fails, the Jordan shape is unexpected,
 * or a VoP integral is not elementary. */
Expr** dsolve_linsys_assemble(Expr* M, Expr* t, const char* xvar,
                              Expr* b, bool b_zero, size_t n);

/* The two system solvers (DSolveSysFn).  Forward-declared here so dsolve.c and the
 * variable-coefficient file can reach them without ad-hoc externs. */
Expr** dsolve_linsys_solve(DSolveProblem* P);
Expr** dsolve_linsys_varcoeff_solve(DSolveProblem* P);

#endif /* MATHILDA_DSOLVE_LINSYS_H */
