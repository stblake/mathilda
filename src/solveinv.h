/*
 * solveinv.h
 *
 * Inverse-function specialist for the public `Solve` router
 * (src/solve.c).  Handles single-variable equations whose outermost
 * dependence on the variable is an elementary invertible head:
 * Log, Exp, Sin/Cos/Tan/Cot/Sec/Csc, the corresponding hyperbolics
 * and inverse trig/hyperbolic, and Power[g, n] for integer n >= 2.
 *
 * Algorithm:
 *
 *   1. Form residual = lhs - rhs.
 *   2. Additive-shift isolation: split residual into "with var" and
 *      "free of var" parts; if the with-var portion is c * head[g(x)]
 *      for an invertible head, rewrite as head[g(x)] == new_rhs.
 *   3. Look up the head in an inverse table that produces one or more
 *      branches (inner_eq, condition).  Multi-branch heads (Sin, Cos,
 *      Exp, ...) introduce a fresh integer parameter C[k] and a
 *      ConditionalExpression wrapper carrying Element[C[k], Integers]
 *      or a principal-strip predicate on Re[a]/Im[a].
 *   4. Recursively solve each inner_eq through solvepoly -> solveinv
 *      -> solverad (depth-capped).  Wrap each {var -> val} solution in
 *      ConditionalExpression[val, condition].
 *   5. Fallback: if head[var] == rhs but the inner solve declines,
 *      return {{var -> InverseFunction[head][rhs]}} and raise
 *      Solve::ifun.
 *
 * Memory contract: the caller retains ownership of `equation`, `var`,
 * `dom`, and `opts`.  Returns a freshly-owned solution
 * List[List[Rule[var, val]], ...] on success, or NULL when the input
 * is not invertible at the top level (in which case the caller hands
 * off to solverad).
 */

#ifndef SOLVEINV_H
#define SOLVEINV_H

#include <stdbool.h>
#include "expr.h"

typedef struct {
    bool        enabled;     /* InverseFunctions; default true (Automatic)  */
    const char* param_head;  /* GeneratedParameters; interned, default "C"  */
} SolveInvOpts;

/* Internal C-callable entry point.  Borrowed args.  Returns a
 * freshly-owned `List[List[Rule[var, val]], ...]` on success, or NULL
 * when the equation is not invertible by this module. */
Expr* solveinv_solve_inverse_equality(Expr* equation,
                                      Expr* var,
                                      Expr* dom,
                                      const SolveInvOpts* opts);

/* Cheap structural probe: returns true iff `expr` has at least one
 * peelable inverse-function head applied to a subterm that contains
 * `var`.  Used as a fast-fail guard in solve.c's dispatch so we do
 * not pay any setup cost on polynomial-only inputs. */
bool solveinv_looks_invertible(const Expr* expr, const Expr* var);

/* Builtin entry for the context-qualified specialist
 *   `Solve`SolveInverseFunctions[lhs == rhs, var]`
 *   `Solve`SolveInverseFunctions[lhs == rhs, var, dom]`. */
Expr* builtin_solve_inverse_functions(Expr* res);

/* Registration hook -- interns the inverse-table head pointers and
 * registers the `Solve`SolveInverseFunctions` qualified builtin. */
void solveinv_init(void);

#endif /* SOLVEINV_H */
