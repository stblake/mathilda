#ifndef MATHILDA_ASSOC_PACKED_H
#define MATHILDA_ASSOC_PACKED_H

/* ---------------------------------------------------------------------------
 * assoc_packed — the machine-buffer paths of the Association builtins.
 *
 * Every function here takes a rank-1 machine buffer (a packed List or a visible
 * NDArray[...]) where the builtin takes a list, and answers EXACTLY what the
 * List path answers for the materialised list -- same keys, same order, same
 * element heads -- without first boxing the input into n Exprs. Each returns
 * NULL when the buffer is outside what it handles (complex or narrow dtype,
 * rank > 1, a non-finite double, a key function outside the compilable subset,
 * ...), and the caller then takes its ordinary List path over a materialised
 * copy (assoc_list_view). NULL never means "wrong", only "not here".
 *
 * WHY THE OUTPUTS ARE BOXED. The groups, positions and values these build sit
 * inside Rule[] nodes inside an Association, and a packed node may never sit
 * inside a plain EXPR_FUNCTION tree (the NO-NESTING INVARIANT, src/pack.c and
 * docs/design/packed_arrays.md D3): the transparency gate scans only the top
 * level of a call, so a buffer nested in an association would reach unaware
 * heads (ReplaceAll, Cases, MatchQ, Level, ...) that walk args[] and treat it as
 * an atom. The win here is therefore on the INPUT side -- no materialisation of
 * the argument, machine-word hashing instead of expr_hash/expr_eq, one compiled
 * pass for the key function -- and on the top-level OUTPUTS that may pack
 * (Values, Keys, and a reducer's group argument).
 * -------------------------------------------------------------------------- */

#include "expr.h"
#include <stdbool.h>

/* A List view of `a`: a visible NDArray or packed list is materialised into
 * *owned (which the caller frees) and returned; anything else is returned as-is
 * with *owned = NULL. Borrows `a`. */
Expr* assoc_list_view(Expr* a, Expr** owned);

/* PositionIndex[buffer]. */
Expr* assoc_packed_positionindex(const Expr* arr);

/* GroupBy[buffer, f] and GroupBy[buffer, f, reducer]: the key function is run
 * over the whole buffer as one compiled loop (autocompile_map_unary) and the
 * machine keys are grouped by word. `reducer` may be NULL. */
Expr* assoc_packed_groupby(const Expr* arr, const Expr* f, const Expr* reducer);

/* GatherBy[buffer, f]: the GroupBy grouping without the keys. */
Expr* assoc_packed_gatherby(const Expr* arr, const Expr* f);

/* CountsBy[buffer, f]. */
Expr* assoc_packed_countsby(const Expr* arr, const Expr* f);

/* AssociationThread[keys, vals] where `keys` is a rank-1 machine buffer and
 * `vals` a buffer or a List of the same length. */
Expr* assoc_packed_thread(const Expr* keys, const Expr* vals);

/* <|word -> count|> over a rank-1 int64 / float64 / bool buffer, first
 * appearance order: the tail of CountsBy once its key function has run (the
 * compiled CountsBy hands its computed key buffer here). */
Expr* assoc_packed_counts_words(const Expr* keys);

/* AssociationMap[f, buffer]. */
Expr* assoc_packed_map(const Expr* f, const Expr* arr);

/* Lookup[assoc, buffer] / Lookup[assoc, buffer, default]: each machine key is
 * probed through the association's persistent index with no boxed key, so a
 * bulk lookup costs O(m) after the one-off index build rather than a fresh
 * O(n) table per call. Absent keys give `deflt` (borrowed; may be NULL for
 * Missing["KeyAbsent", k]). The answer is packed when every value found is a
 * machine Integer (or every one a machine Real) and it clears the threshold. */
Expr* assoc_packed_lookup(const Expr* assoc, const Expr* keys, const Expr* deflt);

/* Keys[assoc] / Values[assoc] as a PACKED list when every key (value) is a
 * machine Integer, every one a machine Real, or every one True/False, and the
 * association is at least the packing threshold long. The result is a fresh
 * top-level List, so packing it cannot put a buffer inside another expression.
 * NULL: build the ordinary boxed List. Borrows `assoc` (an Association or a
 * List of rules, every element a 2-arg rule). */
Expr* assoc_pack_column(const Expr* assoc, bool want_keys);

#endif /* MATHILDA_ASSOC_PACKED_H */
