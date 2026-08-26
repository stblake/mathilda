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
 * FindInstance / CylindricalDecomposition (the remaining Phase-8 companions)
 * are not implemented here yet.
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

/* Register LogicalExpand and NotElement.  Called from reduce_init. */
void reduce_companions_init(void);

#endif /* REDUCE_COMPANIONS_H */
