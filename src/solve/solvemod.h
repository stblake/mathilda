/*
 * solvemod.h
 *
 * Modular solving pre-pass for Solve[..., Modulus -> p].
 *
 * Solve[poly == 0, x, Modulus -> p] solves a single-variable polynomial
 * equation over the finite ring Z/pZ.  This is a fundamentally different
 * value domain from the characteristic-0 dispatch (radicals, Cardano,
 * Root[] objects are all meaningless mod p), so it runs as a pre-pass in
 * builtin_solve before the ordinary specialists.
 *
 * The implementation is residue enumeration: reduce the polynomial over
 * each residue r in [0, p) and keep those where it vanishes mod p.  This
 * is O(p * deg) and correct for prime OR composite p (and for rational
 * coefficients, via internal_mod's modular inverse), covering the common
 * cases; it is capped at SOLVEMOD_P_MAX.  A modulus outside the
 * supported range, a system, a multivariable spec, or a non-polynomial
 * residual all cause a clean REFUSAL (NULL -> Solve stays unevaluated)
 * rather than silently ignoring the Modulus option.
 */
#ifndef SOLVEMOD_H
#define SOLVEMOD_H

#include "expr.h"

/* Largest modulus for which residue enumeration is attempted. */
#define SOLVEMOD_P_MAX 100000

/* Solve `equation` (an Equal[lhs, rhs]) for `vars` over Z/pZ where p is
 * given by the borrowed `modulus` expression (the rhs of Modulus -> p).
 * `dom` is borrowed and ignored (Modulus overrides any domain).  Returns
 * a solution List {{var -> r1}, ...} with r ascending in [0, p), or NULL
 * to refuse (leaving Solve unevaluated).  All arguments are borrowed. */
Expr* solvemod_solve_modular(Expr* equation, Expr* vars, Expr* dom,
                             Expr* modulus);

/* Registers the qualified builtin Solve`SolveModular. */
void solvemod_init(void);

#endif /* SOLVEMOD_H */
