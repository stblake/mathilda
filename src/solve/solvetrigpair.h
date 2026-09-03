/*
 * solvetrigpair.h
 *
 * Argument-pair reducer for the public `Solve` router (src/solve.c).
 * Solves generalized two-argument trigonometric / hyperbolic equations
 *
 *     f[A(var)] (+/-) g[B(var)] == c
 *
 * where f, g are among Sin/Cos/Tan/Cot/Sec/Csc and the hyperbolic
 * analogues, A and B are affine / polynomial in `var` with var-free
 * (possibly symbolic) coefficients, and c is a constant.  These are the
 * cases the single-head inverse-function isolator (solveinv) cannot peel
 * because the variable appears in more than one additive trig term, and
 * that the TrigToExp polynomial-in-u path (solvetrig) rejects because the
 * arguments are not integer multiples of a single variable.
 *
 * Algorithm:
 *
 *   1. residual = lhs - rhs.
 *   2. Reciprocal-normalize: Tan/Cot/Sec/Csc -> Sin/Cos ratios (and the
 *      hyperbolic analogues), then Together; split into a numerator N and
 *      denominator D.
 *   3. Factor N into single-argument atoms via explicit sum-to-product /
 *      reverse-angle-addition identities (a private, deterministic rule
 *      table -- TrigFactor/TrigReduce are unreliable on the bare Sin/Cos
 *      differences this engine must handle).
 *   4. Degenerate cases: N identically 0 -> {{}} (True); N a nonzero
 *      constant -> {} (False).
 *   5. Solve each var-bearing atom H[arg] == 0 through the existing
 *      solveinv specialist (which supplies ArcXxx + period, the C[k]
 *      parameter, and the ConditionalExpression wrapping), and union the
 *      results.
 *   6. Pole gate: drop any solution family on which a denominator (pole)
 *      factor vanishes identically -- this turns Sec[A]-Tan[A]==0 and
 *      Csch[A]-Coth[A]==0 into no-solution.
 *
 * Memory contract: the caller retains ownership of every argument.
 * Returns a freshly-owned `List[List[Rule[var, val]], ...]` on success,
 * or NULL when the equation does not fit this shape (in which case the
 * caller falls through to the solvetrig fallback / decline).
 */

#ifndef SOLVETRIGPAIR_H
#define SOLVETRIGPAIR_H

#include <stdbool.h>
#include "expr.h"
#include "solveinv.h"  /* for SolveInvOpts */

/* Internal C-callable entry.  Borrowed args. */
Expr* solvetrigpair_solve(Expr* equation,
                          Expr* var,
                          Expr* dom,
                          const SolveInvOpts* opts);

/* Cheap structural probe: true iff `expr` carries at least one trig /
 * hyperbolic head over `var`, so the engine is worth attempting.  (Same
 * fast-fail guard the solvetrig pre-pass uses.) */
bool solvetrigpair_looks_applicable(const Expr* expr, const Expr* var);

/* Builtin entry for the context-qualified specialist
 *   `Solve`SolveTrigPair[lhs == rhs, var]`
 *   `Solve`SolveTrigPair[lhs == rhs, var, dom]`. */
Expr* builtin_solve_trig_pair(Expr* res);

/* Registration hook -- parses and caches the rule tables and registers
 * the `Solve`SolveTrigPair` qualified builtin. */
void solvetrigpair_init(void);

#endif /* SOLVETRIGPAIR_H */
