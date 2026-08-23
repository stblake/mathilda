/*
 * reduce.h
 *
 * `Reduce` -- the complete equation/inequality reducer.  Unlike `Solve`, which
 * returns the generic solution of equations as a list of rules, `Reduce`
 * returns a complete, quantifier-free LOGICAL description (an And/Or tree of
 * equations and inequalities) of the solution set over a domain (Complexes,
 * Reals, Integers, Rationals).
 *
 * Architecture and phased roadmap: REDUCE_PLAN.md.  This is the Phase-0
 * front-end: it parses arguments like Solve, normalises the input into the
 * internal DNF layer (reduce_form.h), and returns True/False when the statement
 * decides; the per-domain solving engines are wired in later phases.
 */
#ifndef REDUCE_H
#define REDUCE_H

#include "expr.h"

Expr* builtin_reduce(Expr* res);

/* Register Reduce (Protected) and its docstring.  Idempotent. */
void reduce_init(void);

#endif /* REDUCE_H */
