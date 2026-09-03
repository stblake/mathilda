/*
 * reduce_trigregion.h
 *
 * Reduce over the Reals for a univariate CONJUNCTION that pairs one periodic
 * (circular-trig / hyperbolic) equation `f[A(x)] == c` with inequality
 * constraints that bound the variable to a finite region, e.g.
 *
 *     Reduce[Sin[x] == 1/2 && 0 < x < 2 Pi, x]
 *       -> x == Pi/6 || x == 5 Pi/6
 *
 * The single-equation dispatch in reduce.c solves the periodic equation into
 * its infinite families `a + p C[k]`; a bounding region selects the finitely
 * many members inside it.  The sign-diagram engines cannot represent a periodic
 * family, so before this they returned a WRONG `False` for such conjunctions.
 *
 * Method (sound by construction): solve the equation into its real families,
 * substitute each family `a + p k` into the remaining constraints, `Reduce`
 * that over the *real* parameter `k` to a bounded interval, enumerate the
 * integers in it (with slack), and keep only the candidates that make the FULL
 * original statement evaluate to `True`.  The numeric interval only bounds the
 * search; the exact re-check is what makes the result sound.
 */
#ifndef REDUCE_TRIGREGION_H
#define REDUCE_TRIGREGION_H

#include "expr.h"
#include "reduce_opts.h"

/* Solve `peq && region_stmt` over the Reals for `var`.
 *   peq         a periodic equation `Equal[lhs, rhs]` over `var`
 *   region_stmt the conjunction of the remaining constraints (an Expr in `var`)
 *   var         the solve variable
 *
 * Returns (owned) an `Or` of `var == value` for the finitely many real
 * solutions, `False` if there are none, or NULL to decline (unrecognised
 * region, an unbounded feasible set, a non-integer / non-real family, ...).
 * The caller must NOT fall through to the sign diagram on NULL: it is unsound
 * for this shape, so a decline means "leave the statement unevaluated". */
Expr* reduce_periodic_region(const Expr* peq, const Expr* region_stmt,
                             const Expr* var, const ReduceOpts* opts);

/* Solve a univariate conjunction over the Reals that contains at least one
 * trig/hyperbolic INEQUALITY over `var` together with inequalities that bound
 * `var` to a finite region, e.g.
 *
 *     Reduce[Sin[x] > 1/2 && 0 < x < 2 Pi, x]  ->  Pi/6 < x < 5 Pi/6
 *
 * The zeros of the trig atoms partition each bounded feasible interval into
 * sign-constant cells; the FULL statement is evaluated at a sample point of
 * each cell (numerically in the interior, exactly at each boundary), and the
 * true cells are merged into a union of intervals / points.  Only pole-free
 * real heads (Sin, Cos, Sinh, Cosh, Tanh, Sech) are handled; a pole-bearing
 * head (Tan, Cot, Sec, Csc, Coth, Csch), an unbounded region, or a shape it
 * cannot read makes it return NULL (the caller then leaves the statement
 * unevaluated -- the sign diagram is unsound for this shape). */
Expr* reduce_trig_ineq_region(const Expr* conj, const Expr* var,
                              const ReduceOpts* opts);

#endif /* REDUCE_TRIGREGION_H */
