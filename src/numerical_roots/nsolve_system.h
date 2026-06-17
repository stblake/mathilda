/* nsolve_system.h — numerical solver for square polynomial systems.
 *
 * Two methods, both returning every isolated complex solution of a
 * zero-dimensional square system  f_1 == ... == f_n == 0  in variables
 * x_1 .. x_n:
 *
 *   NSYS_HOMOTOPY      total-degree polynomial homotopy continuation
 *                      (the gamma-trick start system + predictor/corrector
 *                      path tracking, then Newton refinement).
 *   NSYS_ENDOMORPHISM  the eigenvalue / multiplication-matrix method:
 *                      Groebner basis over Q, quotient-ring monomial basis,
 *                      eigenstructure of multiplication by a generic linear
 *                      form.
 *
 * Both build their numeric work from exact Expr coefficients, so integer,
 * rational, real, and complex coefficients are all supported, with solutions
 * over the complexes (the Reals domain just discards the non-real ones).
 *
 * Returns a freshly owned List of solution rule-lists
 *   {{x -> r, y -> s, ...}, ...}
 * or NULL when the system is outside the supported envelope (non-square,
 * positive-dimensional, too large, or not a numeric polynomial system in the
 * given variables), letting the caller fall back or leave NSolve unevaluated.
 */
#ifndef MATHILDA_NSOLVE_SYSTEM_H
#define MATHILDA_NSOLVE_SYSTEM_H

#include "expr.h"
#include <stdbool.h>

typedef enum { NSYS_HOMOTOPY = 0, NSYS_ENDOMORPHISM } NSysMethod;

Expr* nsolve_polynomial_system(Expr** polys, int npoly, Expr** vars, int nvar,
                               NSysMethod method, bool reals_only,
                               bool want_machine, long bits, long max_roots,
                               int verify, unsigned long seed);

/* Elimination / triangular solver (Method -> "Symbolic", and the fallback
 * for nsolve_polynomial_system): a lexicographic Gröbner basis solved by
 * NRoots back-substitution.  Same return contract as above. */
Expr* nsolve_system_eliminate(Expr** polys, int npoly, Expr** vars, int nvar,
                              bool reals_only, bool want_machine, long bits,
                              long max_roots, int verify);

#endif /* MATHILDA_NSOLVE_SYSTEM_H */
