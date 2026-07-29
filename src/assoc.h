#ifndef MATHILDA_ASSOC_H
#define MATHILDA_ASSOC_H

/* ---------------------------------------------------------------------------
 * Association — the Wolfram-Language <| key -> value, ... |> data structure.
 *
 * Representation.  An association is stored as an ordinary expression:
 *
 *     <|a -> 1, b -> 2|>   ==   Association[Rule[a, 1], Rule[b, 2]]
 *
 * i.e. an EXPR_FUNCTION node with head `Association` whose arguments are
 * two-argument `Rule[k, v]` nodes.  Keeping associations as plain expressions
 * means the whole "everything is an expression" toolchain — expr_copy,
 * expr_free, expr_eq, expr_hash, pattern matching, ReplaceAll, FullForm — all
 * keep working with no special cases.
 *
 * Keys are unique and insertion-ordered: the first occurrence fixes a key's
 * position, the last occurrence its value (matching Wolfram semantics).
 *
 * Performance.  The naive way to build/query an association is an O(n) linear
 * scan per key, giving O(n^2) construction and O(n*m) bulk lookups.  Instead
 * this module builds a transient open-addressing hash index (keyed by
 * expr_hash / expr_eq) so that construction, de-duplication, Merge, KeyDrop,
 * bulk Lookup, Counts and GroupBy are all amortised O(n) — see assoc.c.
 * -------------------------------------------------------------------------- */

#include "expr.h"
#include <stdbool.h>
#include <stddef.h>

/* True when `e` is an Association[...] node (head == SYM_Association). */
bool is_association(const Expr* e);

/* Build a canonical Association from `count` Rule[k,v] nodes.  The rules are
 * copied (the caller keeps ownership of `rules`).  Duplicate keys collapse
 * with last-value-wins while preserving first-occurrence order.  O(count). */
Expr* assoc_from_rules(Expr** rules, size_t count);

/* Thread `f` over the values of an association, preserving keys:
 * <|k -> v|>  ->  <|k -> f[v]|>. The f[v] applications are left for the
 * evaluator to reduce. Used by Map so associations flow through the functional
 * toolchain the Wolfram way. Caller owns the result. */
Expr* assoc_map_values(Expr* f, const Expr* assoc);

/* Keep the entries whose value satisfies the predicate `pred` (pred[v] must
 * evaluate to True), up to `max` of them (max < 0 means no limit). Preserves
 * order. Caller owns the result. */
Expr* assoc_select_values(Expr* pred, const Expr* assoc, int64_t max);

/* The values of an association as a List (deep copies). Caller owns it. */
Expr* assoc_values_list(const Expr* assoc);

/* If res is `head[assoc, rest...]`, rebuild it as `head[Values[assoc], rest...]`
 * and evaluate — so aggregations (Total, Min, Max, Mean, ...) act on the values
 * of an association. Returns NULL (leave to the caller) when the first argument
 * is not an association. Caller owns any returned result. */
Expr* assoc_apply_over_values(Expr* res);

/* An association reordered by its values in canonical order (keys follow their
 * values). Caller owns the result. */
Expr* assoc_sort_by_value(const Expr* assoc);

/* For windowed value transforms (Accumulate, Differences) whose result is a
 * List parallel to the entries: rewrite head[assoc, rest...] as
 * head[Values[assoc], rest...], evaluate, then re-pair the resulting values
 * with the trailing keys (so a dropped leading value drops the leading key).
 * Returns NULL when the first argument is not an association or the transform
 * did not yield a same-or-shorter List. Caller owns any returned result. */
Expr* assoc_rekey_over_values(Expr* res);

/* Re-pair the elements of a values List with the trailing keys of assoc,
 * giving an association (used when the collection is not the first argument,
 * e.g. FoldList[f, assoc]). Returns NULL when `values` is not a List or is
 * longer than assoc. Caller owns any returned result. */
Expr* assoc_rekey_from_list(const Expr* assoc, const Expr* values);

/* DeleteCases[assoc, pattern]: keep the entries whose value does NOT match the
 * pattern (matching is done via MatchQ on each value). Caller owns the result. */
Expr* assoc_delete_cases(const Expr* assoc, Expr* pattern);

/* DeleteDuplicates[assoc]: keep the first entry for each distinct value
 * (hash-indexed, O(n)). Caller owns the result. */
Expr* assoc_delete_duplicate_values(const Expr* assoc);

/* GatherBy[list, f]: group elements with equal f[element] into a list of
 * sublists (first-appearance order). Caller owns the result. */
Expr* builtin_gatherby(Expr* res);

/* Shared grouping engine behind GatherBy[list, f] and Gather[list]: groups the
 * elements of `list` (a List, or an Association whose entries group by value)
 * into a list of sublists in first-appearance order, hash-indexed in O(n).
 *
 * `f` is the key function; passing NULL selects the identity key — the element
 * itself is the group key and no function application is evaluated, which is
 * the Gather[list] path. Returns NULL if `list` is neither a List nor an
 * Association. Caller owns the result. */
Expr* assoc_gather_core(Expr* list, Expr* f);

/* KeyUnion[{assoc1, ...}]: pad every association to the union of all keys,
 * filling absent values with Missing["KeyAbsent", key]. Caller owns result. */
Expr* builtin_keyunion(Expr* res);

/* A FRESH association entry carrying `entry`'s key and head but the value
 * `newval` (adopted). Caller owns the result; `entry` is left untouched.
 *
 * Use this for EVERY value rewrite on an association. expr_copy is a refcount
 * bump (src/expr.c), so an entry reached from a caller's input is the very node
 * every other holder of that association sees -- assigning into its args[1] is
 * how `u2 = u1; u2[["x"]] = 9` silently corrupted `u1`. Copying the head keeps
 * a RuleDelayed entry delayed.
 *
 * `entry` MUST be a two-argument rule: is_association() only checks the head,
 * so Association[1, 2] can exist and callers must guard with is_rule2 first. */
Expr* assoc_entry_with_value(const Expr* entry, Expr* newval);

/* True for a two-argument Rule/RuleDelayed, i.e. a well-formed entry. */
bool is_rule2(const Expr* e);

/* Register every Association-family builtin, with attributes and docstrings. */
void assoc_init(void);

#endif /* MATHILDA_ASSOC_H */
