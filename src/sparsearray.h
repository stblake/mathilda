#ifndef SPARSEARRAY_H
#define SPARSEARRAY_H

#include "expr.h"
#include <stdbool.h>

/*
 * SparseArray support.
 *
 * Mathilda has no sparse storage: SparseArray[...] is an inert expression, and
 * every consumer that meets one sees the specification exactly as written. The
 * one conversion Mathematica code relies on everywhere is Normal, so that is
 * what this module provides -- the dense nested List a specification denotes.
 *
 * Accepted specifications (the documented constructor forms, plus the
 * InputForm of Mathematica's own internal representation so a value pasted out
 * of Mathematica round-trips):
 *
 *   SparseArray[{pos1 -> v1, ...}]                 dims inferred from positions
 *   SparseArray[{pos1 -> v1, ...}, dims]           dims an Integer or a List
 *   SparseArray[{pos1 -> v1, ...}, dims, default]
 *   SparseArray[{pos1, pos2, ...} -> {v1, v2, ...}, ...]
 *   SparseArray[{i_, j_} :> f[i, j], dims]         pattern rules (dims needed)
 *   SparseArray[Band[start] -> v, dims]            Band[start], Band[start, end],
 *                                                  Band[start, end, step]; v a
 *                                                  scalar or a List of values
 *   SparseArray[list], SparseArray[list, dims, default]    a dense array
 *   SparseArray[Automatic, dims, default, {1, {rowptr, colidx}, vals}]
 *
 * A position is a List of Integers (an Integer for a vector); negative indices
 * count from the end once dims are known. Where several rules name one
 * position the FIRST wins, as in Mathematica. Unfilled positions take the
 * default (0 unless given).
 */

/* True when `e` is a SparseArray[...] expression (any arity). */
bool is_sparse_array(const Expr* e);

/* The dense nested List that SparseArray expression `sa` denotes, or NULL when
 * `sa` is not a SparseArray or its specification is not one of the forms above
 * (malformed, pattern rules without dims, a result too large to allocate).
 * Borrows `sa`; the result is owned by the caller and is unevaluated -- pattern
 * right-hand sides are substituted, and the caller's evaluator finishes them. */
Expr* sparse_array_to_dense(const Expr* sa);

#endif /* SPARSEARRAY_H */
