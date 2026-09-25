#ifndef MATHILDA_ASSOC_STRUCT_H
#define MATHILDA_ASSOC_STRUCT_H

/* ---------------------------------------------------------------------------
 * Association atomicity for structural traversal.
 *
 * Physically an association is the ordinary expression
 *
 *     <|a -> 1, b -> 2|>   ==   Association[Rule[a, 1], Rule[b, 2]]
 *
 * but in the Wolfram Language it is an ATOM: AtomQ gives True, and every
 * level-based function sees only its VALUES. The Rule wrappers and the keys
 * are not parts of the expression:
 *
 *     Level[<|a -> 1, b -> 2|>, {1}]     {1, 2}
 *     Depth[<|a -> 1|>]                  2        (not 3)
 *     FreeQ[<|a -> 1|>, a]               True
 *     <|a -> 1, b -> 2|> /. b -> a       <|a -> 1, b -> 2|>   (keys untouched)
 *     Position[{<|a -> 1|>}, 1]          {{1, Key[a]}}
 *     Map[f, <|a -> {1}|>, {2}]          <|a -> {f[1]}|>
 *
 * A walker that descends through data.function.args would instead visit the
 * Rule nodes (one level too many) and the keys (which must never be matched
 * or rewritten), and a rebuild of the args could produce a malformed
 * Association[f[a -> 1]]. Every level walker -- Level, Depth, LeafCount,
 * Map/MapAll/Apply/Scan/FreeQ, Replace/ReplaceAll, Cases/Position/Count/
 * MemberQ/DeleteCases -- therefore reaches an expression's parts through the
 * three helpers below instead of through args[] directly:
 *
 *   assoc_is_wellformed(e)       is e an association in the atomic sense?
 *   struct_part(e, i, assoc)     part i: args[i], or the VALUE of entry i
 *   struct_rebuild(...)          reassemble e from new parts, keeping keys
 *
 * Only a WELL-FORMED association (every argument a two-argument Rule or
 * RuleDelayed) is atomic. A malformed Association[1, 2] -- which the
 * constructor leaves unevaluated, exactly as Mathematica does -- is an ordinary
 * compound expression: AtomQ[Association[1, 2]] is False and Level sees 1, 2.
 *
 * Deliberate difference: Mathematica distinguishes an evaluated association
 * (atomic) from an unevaluated Association[...] expression held inside Hold
 * (not atomic), so Hold[<|a -> 1|>] /. a -> b gives Hold[Association[b -> 1]].
 * Mathilda has one representation for both and treats every well-formed node
 * as atomic, so the key is left alone there too.
 *
 * Header-only (static inline): the helpers are one-liners called in the hot
 * loop of every walker, and keeping them here needs no build-list changes.
 * assoc_is_wellformed is O(n) in the entry count, and each walker calls it
 * once per visited node, so a traversal stays linear in the tree size.
 * -------------------------------------------------------------------------- */

#include "expr.h"
#include "sym_names.h"
#include <stdbool.h>
#include <stdlib.h>

/* A two-argument Rule[k, v] or RuleDelayed[k, v]: a well-formed entry. */
static inline bool assoc_struct_entry_ok(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2)
        return false;
    const Expr* h = e->data.function.head;
    return h->type == EXPR_SYMBOL &&
           (h->data.symbol.name == SYM_Rule || h->data.symbol.name == SYM_RuleDelayed);
}

/* True when e is Association[entries...] and every entry is a two-argument
 * Rule/RuleDelayed -- the only shape that is atomic, that AssociationQ
 * accepts, and whose entries' args[0]/args[1] may be read. */
static inline bool assoc_is_wellformed(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_Association)
        return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!assoc_struct_entry_ok(e->data.function.args[i])) return false;
    return true;
}

/* Structural part i (0-based) of the compound expression e. For an atomic
 * association (`assoc` == assoc_is_wellformed(e), computed once by the
 * caller) that is the VALUE of entry i; otherwise the plain argument.
 * Borrowed pointer. The part count is data.function.arg_count either way. */
static inline Expr* struct_part(const Expr* e, size_t i, bool assoc) {
    Expr* a = e->data.function.args[i];
    return assoc ? a->data.function.args[1] : a;
}

/* The position component of part i: Key[k] for an association entry, the
 * Integer i+1 otherwise. Owned result. */
static inline Expr* struct_part_component(const Expr* e, size_t i, bool assoc) {
    if (!assoc) return expr_new_integer((int64_t)(i + 1));
    Expr* k = expr_copy(e->data.function.args[i]->data.function.args[0]);
    return expr_new_function(expr_new_symbol(SYM_Key), &k, 1);
}

/* Reassemble e from a (possibly rewritten) head and n rewritten parts.
 * ADOPTS `head` and every parts[i]; the caller frees only the parts array.
 *
 * Not an association: head[parts...].
 *
 * Atomic association (n == its entry count): while the head is still the
 * symbol Association, each entry keeps its key and its Rule/RuleDelayed head
 * and takes the new value -- keys are never rewritten, so no duplicates can
 * appear and no re-canonicalisation is needed. When nothing changed (every
 * part is the very node it replaced, which is what an untouched subtree's
 * refcount-bumping expr_copy returns) the original node is shared instead,
 * keeping its key index.
 *
 * If a Heads->True walker rewrote the Association head itself the result is
 * no longer an association. Mathematica is not uniform here: Map passes the
 * values (Map[f, <|a -> 1|>, Heads -> True] is f[Association][f[1]]) while
 * the replace family keeps the entries (<|a -> 1|> /. Association -> g is
 * g[a -> 1]); `keep_entries_on_new_head` selects the latter. */
static inline Expr* struct_rebuild(const Expr* e, bool assoc, Expr* head,
                                   Expr** parts, size_t n,
                                   bool keep_entries_on_new_head) {
    if (!assoc) return expr_new_function(head, parts, n);

    bool same_head = head->type == EXPR_SYMBOL &&
                     head->data.symbol.name == SYM_Association;
    if (!same_head && !keep_entries_on_new_head)
        return expr_new_function(head, parts, n);

    if (same_head) {
        bool unchanged = true;
        for (size_t i = 0; i < n && unchanged; i++)
            if (parts[i] != e->data.function.args[i]->data.function.args[1])
                unchanged = false;
        if (unchanged) {
            for (size_t i = 0; i < n; i++) expr_free(parts[i]);
            expr_free(head);
            return expr_copy((Expr*)e);
        }
    }

    Expr** entries = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        const Expr* old = e->data.function.args[i];
        Expr* rargs[2] = { expr_copy(old->data.function.args[0]), parts[i] };
        entries[i] = expr_new_function(expr_copy(old->data.function.head), rargs, 2);
    }
    Expr* out = expr_new_function(head, entries, n);
    free(entries);
    return out;
}

#endif /* MATHILDA_ASSOC_STRUCT_H */
