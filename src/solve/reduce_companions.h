/*
 * reduce_companions.h
 *
 * Companion builtins for `Reduce` (REDUCE_PLAN.md, Phase 8).
 *
 * v1 ships `LogicalExpand`: it expands a logical combination of equations,
 * inequalities and Boolean atoms into disjunctive normal form (an Or of Ands)
 * with idempotence / complementation / absorption contractions, collapsing to
 * True (tautology) or False (contradiction) when the statement decides.
 *
 * Unlike the rest of Reduce, LogicalExpand treats every non-connective
 * subexpression as an OPAQUE Boolean atom (a symbol p, a relation x==a, an
 * Element[..] membership) -- it does NO domain reasoning, exactly as
 * Mathematica's LogicalExpand does.  It therefore cannot reuse the RForm/RAtom
 * layer (which canonicalises leaves into polynomial relations); the engine here
 * is a standalone Expr-level DNF distributor.
 *
 * A minimal `NotElement` head ships alongside so negated membership prints as
 * `x` \[NotElement] `dom` (Mathematica-faithful) rather than Not[Element[..]].
 *
 * `FindInstance[expr, vars, dom, n]` finds up to n witness points satisfying
 * expr over Complexes / Reals / Integers / Rationals / Booleans, returned in
 * Solve's form.  Rather than reaching into the three static CAD paths, it reads
 * witnesses off the PUBLIC cylindrical outputs of Reduce and Solve, samples
 * intervals with rru_rational_between, and VERIFIES every candidate against the
 * original expr -- so a returned instance is always a true solution, {} means
 * provably empty, and it stays unevaluated when it can neither exhibit an
 * instance nor prove emptiness.  The Booleans domain reuses the DNF engine here.
 *
 * Four extensions reach witnesses the raw Reduce/Solve outputs do not surface:
 *   - generated parameters (`x -> ConditionalExpression[.., C[1]>=1]`) are
 *     instantiated over a small integer grid -- reaching the fundamental Pell
 *     solution x^2 - 61 y^2 == 1 at C[1] == 1;
 *   - a single generated parameter the grid misses is solved against the
 *     remaining constraints over the Reals and rounded to an integer -- reaching
 *     the periodic Sin[1/x] == 0 && 0 < x < 10^-5 at C[1] == 15916;
 *   - indexed variables c[i] are matched structurally, so systems in c[1..n]
 *     (e.g. 0/1 knapsacks) are accepted;
 *   - over the Integers, a bounded integer-box search gives a Diophantine
 *     witness, or a finite-domain emptiness proof ({} when the box is exhausted
 *     or a linear reach-range excludes the target);
 *   - for transcendental / inexact Real systems (which Reduce cannot soundly
 *     decide -- it wrongly reports False for 0<x<0.001 && Sin[1/x]>0.999) a
 *     numerical feasibility search (NMinimize[{0, expr}, vars]) supplies a
 *     verified inexact witness, and Reduce's False is not trusted as {};
 *   - a declined positive-dimensional polynomial system is proven empty via a
 *     Rabinowitsch Groebner certificate (basis {1} over C ⊇ R ⊇ Z).
 *
 * The Booleans domain also expands `Equivalent` (via the DNF engine here), which
 * additionally now evaluates as a builtin (src/boolean.c).
 *
 * CylindricalDecomposition[expr, vars] gives a cylindrical algebraic
 * decomposition of the real solution set of `expr`.  It is Reals-only (the sole
 * semantic difference from Reduce over an ordered field) and is implemented as a
 * thin front-end that forces the Reals domain and delegates to Reduce, whose
 * Reals engine (Fourier-Motzkin / CAD / sign diagram) already emits the merged
 * cylindrical formula.
 */
#ifndef REDUCE_COMPANIONS_H
#define REDUCE_COMPANIONS_H

#include "expr.h"

/* LogicalExpand[expr] -- distribute `expr` to disjunctive normal form (Or of
 * Ands of literals), with contractions; True / False when it decides.  Always
 * returns a freshly-owned Expr for a single argument (a bare atom returns
 * itself); NULL only for the wrong arity, leaving the call unevaluated. */
Expr* builtin_logical_expand(Expr* res);

/* NotElement[x, dom] -- the negation of Element.  Decides to True / False when
 * Element[x, dom] decides, else stays symbolic (NULL). */
Expr* builtin_not_element(Expr* res);

/* FindInstance[expr, vars, dom, n] -- up to n verified witness points as a List
 * of Solve-style rule-lists; {} if provably empty; NULL (unevaluated) when it
 * can neither exhibit an instance nor prove emptiness. */
Expr* builtin_find_instance(Expr* res);

/* CylindricalDecomposition[expr, vars] -- the cylindrical algebraic
 * decomposition of the real solution set of `expr`, a quantifier-free And/Or
 * formula in which each variable is bounded cylindrically in terms of the
 * earlier ones.  Reals-only; forces the Reals domain and delegates to Reduce.
 * Returns True / False / the formula, or NULL (unevaluated) when the engine
 * cannot decide.  A trailing Reals domain arg is accepted and ignored; option
 * Rules are forwarded to Reduce. */
Expr* builtin_cylindrical_decomposition(Expr* res);

/* Register LogicalExpand, NotElement, FindInstance and
 * CylindricalDecomposition.  Called from reduce_init. */
void reduce_companions_init(void);

#endif /* REDUCE_COMPANIONS_H */
