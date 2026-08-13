#include "print.h"

#include "funcprog.h"
#include "eval.h"
#include "arithmetic.h"
#include "match.h"
#include "sym_names.h"
#include "assoc.h"
#include "common.h"
#include "numloop.h"
#include "ndarray.h"
#include "ndstruct.h"   /* ndstruct_delist_repack — packed-argument fallback */
#include "ndreduce.h"   /* ndred_scan — Fold/FoldList over an associative operator */
#include "pack.h"       /* pack_offer — Select narrows one numeric list to another */
#include "part.h"
#include "checked_int.h"  /* ci_add_i64 / ci_mul_i64 — the int64 MapThread arm */
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>       /* isfinite — Outer's buffer path abandons on overflow */

typedef struct {
    int64_t min;
    int64_t max;
    bool heads;
} LevelSpec;

static int64_t get_depth(Expr* e) {
    if (e->type != EXPR_FUNCTION) return 1;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_Rational || h == SYM_Complex) return 1;
    }
    int64_t max_d = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        int64_t d = get_depth(e->data.function.args[i]);
        if (d > max_d) max_d = d;
    }
    return max_d + 1;
}

/* The Infinity sentinel used for an unbounded upper level. */
#define LEVEL_SPEC_INF 1000000

/* is_rule2 — true for a two-argument Rule[a,b] or RuleDelayed[a,b], the shape
 * of an option and of an association entry — now lives in assoc.c, shared with
 * the association machinery that needs the same guard. See assoc.h. */

/* Lenient level-spec parser: anything it does not recognise (including
 * {n, Infinity}) silently leaves the caller's defaults in place. Used by
 * Map / Apply / MapAll, which have always behaved this way; Scan and
 * MapIndexed use parse_level_spec_strict below instead. */
static LevelSpec parse_level_spec(Expr* ls, int64_t default_min, int64_t default_max) {
    LevelSpec spec = {default_min, default_max, false};
    if (!ls) return spec;

    if (ls->type == EXPR_INTEGER) {
        spec.min = 1;
        spec.max = ls->data.integer;
    } else if (head_is(ls, SYM_List)) {
        if (ls->data.function.arg_count == 1 && ls->data.function.args[0]->type == EXPR_INTEGER) {
            spec.min = spec.max = ls->data.function.args[0]->data.integer;
        } else if (ls->data.function.arg_count == 2 &&
                   ls->data.function.args[0]->type == EXPR_INTEGER &&
                   ls->data.function.args[1]->type == EXPR_INTEGER) {
            spec.min = ls->data.function.args[0]->data.integer;
            spec.max = ls->data.function.args[1]->data.integer;
        }
    } else if (ls->type == EXPR_SYMBOL && ls->data.symbol.name == SYM_Infinity) {
        spec.min = 1;
        spec.max = LEVEL_SPEC_INF;
    }
    return spec;
}

/* Parse one level bound: an integer yields its value; Infinity yields
 * +LEVEL_SPEC_INF. Returns false for anything else. */
static bool level_bound(Expr* e, int64_t* out) {
    if (e->type == EXPR_INTEGER) { *out = e->data.integer; return true; }
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Infinity) {
        *out = LEVEL_SPEC_INF;
        return true;
    }
    return false;
}

/* Strict level-spec parser, shared by Scan and MapIndexed. Unlike
 * parse_level_spec it accepts Infinity in *either* slot of a two-element spec
 * (e.g. {0, Infinity}) and, crucially, reports failure instead of silently
 * falling back to a default — so an unrecognised spec leaves the call
 * unevaluated rather than quietly mapping somewhere else. Map / Apply /
 * MapAll deliberately still use the lenient parser above.
 *
 * Accepts n, {n}, {n1,n2} (each bound integer or Infinity), and Infinity.
 * Default (ls == NULL) is {1, 1}; `heads` is always reset here so the
 * subsequent parse_options call is the single source of that flag. */
static bool parse_level_spec_strict(Expr* ls, LevelSpec* spec) {
    spec->min = 1;
    spec->max = 1;
    spec->heads = false;
    if (!ls) return true;

    if (ls->type == EXPR_INTEGER) {          /* n -> {1, n} */
        spec->min = 1;
        spec->max = ls->data.integer;
        return true;
    }
    if (ls->type == EXPR_SYMBOL && ls->data.symbol.name == SYM_Infinity) {
        spec->min = 1;                        /* Infinity -> {1, Infinity} */
        spec->max = LEVEL_SPEC_INF;
        return true;
    }
    if (head_is(ls, SYM_List)) {
        size_t n = ls->data.function.arg_count;
        if (n == 1) {                         /* {n} -> {n, n} */
            int64_t v;
            if (level_bound(ls->data.function.args[0], &v)) {
                spec->min = spec->max = v;
                return true;
            }
        } else if (n == 2) {                  /* {n1, n2} */
            int64_t a, b;
            if (level_bound(ls->data.function.args[0], &a) &&
                level_bound(ls->data.function.args[1], &b)) {
                spec->min = a;
                spec->max = b;
                return true;
            }
        }
    }
    return false;
}

static void parse_options(Expr* res, size_t start_idx, LevelSpec* spec) {
    for (size_t i = start_idx; i < res->data.function.arg_count; i++) {
        Expr* opt = res->data.function.args[i];
        if (is_rule2(opt)) {
            if (opt->data.function.args[0]->type == EXPR_SYMBOL &&
                opt->data.function.args[0]->data.symbol.name == SYM_Heads) {
                if (opt->data.function.args[1]->type == EXPR_SYMBOL &&
                    opt->data.function.args[1]->data.symbol.name == SYM_True) {
                    spec->heads = true;
                }
            }
        }
    }
}

/* ------------------- Apply ------------------- */

static Expr* apply_at_level(Expr* f, Expr* expr, int64_t current_level, LevelSpec spec);

/* Recurse into one child at the next level down, but short-circuit the common
 * case where a non-negative spec has already passed its maximum: no application
 * can occur in the child or any deeper node, so refcount-share it (O(1)) instead
 * of paying a recursive apply_at_level that would do exactly that after
 * re-checking. This is what makes a wide MapApply (`f @@@ list`) cheap -- every
 * element's arguments sit one level past max and would otherwise each cost a
 * recursive call that only shares. Inlined by the compiler at -O3. */
static Expr* apply_child(Expr* f, Expr* child, int64_t child_level, LevelSpec spec) {
    if (spec.min >= 0 && child_level > spec.max)
        return expr_copy(child);
    return apply_at_level(f, child, child_level, spec);
}

static Expr* apply_at_level(Expr* f, Expr* expr, int64_t current_level, LevelSpec spec) {
    if (expr->type != EXPR_FUNCTION) {
        return expr_copy(expr);
    }

    /* With a non-negative level spec the negative-depth branch of should_apply
     * can never fire, so get_depth's full-subtree traversal is dead weight.
     * Moreover, once we are past the maximum level no application can occur here
     * or at any deeper node (levels only increase with depth), so the whole
     * subtree passes through unchanged -- an O(1) refcount bump instead of a
     * structural clone. Only the negative-spec case still needs get_depth. */
    bool should_apply;
    if (spec.min >= 0) {
        if (current_level > spec.max) return expr_copy(expr);
        should_apply = (current_level >= spec.min && current_level <= spec.max);
    } else {
        int64_t d = get_depth(expr);
        should_apply = (current_level >= spec.min && current_level <= spec.max) ||
                       (-d >= spec.min && -d <= spec.max);
    }

    if (should_apply) {
        size_t count = expr->data.function.arg_count;
        Expr** args = malloc(sizeof(Expr*) * count);
        for (size_t i = 0; i < count; i++) {
            args[i] = apply_child(f, expr->data.function.args[i], current_level + 1, spec);
        }
        Expr* new_func = expr_new_function(expr_copy(f), args, count);
        free(args);
        Expr* eval_res = evaluate(new_func);
        expr_free(new_func);
        return eval_res;
    }

    size_t count = expr->data.function.arg_count;
    Expr** new_args = malloc(sizeof(Expr*) * count);
    for (size_t i = 0; i < count; i++) {
        new_args[i] = apply_child(f, expr->data.function.args[i], current_level + 1, spec);
    }
    
    Expr* new_head = NULL;
    if (spec.heads) {
        new_head = apply_at_level(f, expr->data.function.head, current_level + 1, spec);
    } else {
        new_head = expr_copy(expr->data.function.head);
    }
    
    Expr* result = expr_new_function(new_head, new_args, count);
    free(new_args);
    return result;
}

Expr* builtin_apply(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    
    Expr* f = res->data.function.args[0];
    Expr* expr = res->data.function.args[1];

    Expr* ls = (res->data.function.arg_count >= 3) ? res->data.function.args[2] : NULL;
    if (is_rule2(ls)) ls = NULL;            /* option in the level-spec slot */

    /* Apply[f, assoc] uses the association's values as f's arguments:
     * f @@ <|k1 -> v1, ...|> is f[v1, ...] (matching Wolfram, and consistent
     * with Total = Plus @@ assoc). Only the default level applies here; an
     * explicit level spec falls through to the generic traversal. */
    if (ls == NULL && is_association(expr)) {
        size_t n = expr->data.function.arg_count;
        Expr** vargs = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t i = 0; i < n; i++)
            vargs[i] = expr_copy(expr->data.function.args[i]->data.function.args[1]);
        Expr* call = expr_new_function(expr_copy(f), vargs, n);
        free(vargs);
        Expr* r = evaluate(call);
        expr_free(call);
        return r;
    }

    /* Specialised top-level fast path for `f @@ expr` (the default, no explicit
     * level spec). Applying at level 0 only rewrites the head, so `f @@ g[a,b,c]`
     * is exactly `evaluate(f[a,b,c])` with every element untouched. The general
     * apply_at_level below would instead call get_depth over the whole tree (its
     * only consumer is the negative-level branch a default spec never takes) and
     * structurally clone every element subtree. Here we build f[elements...] by
     * refcount-sharing each element (Expr trees are immutable; expr_copy is an
     * O(1) refcount bump) and evaluate once -- no traversal, no per-subtree clone.
     * Provably equivalent to apply_at_level(f, expr, 0, {0,0}). Mirrors the
     * association fast path above. Explicit `0`/`{0}`/`{0,0}` parse to non-default
     * specs and keep the general path, so only the ls == NULL form lands here. */
    if (ls == NULL && expr->type == EXPR_FUNCTION) {
        size_t count = expr->data.function.arg_count;
        Expr** args = malloc(sizeof(Expr*) * (count ? count : 1));
        for (size_t i = 0; i < count; i++)
            args[i] = expr_copy(expr->data.function.args[i]);
        Expr* call = expr_new_function(expr_copy(f), args, count);
        free(args);
        Expr* r = evaluate(call);
        expr_free(call);
        return r;
    }

    LevelSpec spec = parse_level_spec(ls, 0, 0);
    parse_options(res, ls ? 3 : 2, &spec);

    return apply_at_level(f, expr, 0, spec);
}

/* ------------------- Map ------------------- */

static Expr* map_at_level(Expr* f, Expr* expr, int64_t current_level, LevelSpec spec) {
    /* Non-negative spec: once we are past the maximum level, nothing here or at
     * any deeper node can be mapped, so the whole subtree passes through by an
     * O(1) refcount bump instead of a full bottom-up descent + rebuild. This is
     * what stops the common `f /@ list` from walking (and rebuilding) the entire
     * tree when f only touches level 1 -- cost was growing with element depth.
     * Only the negative-spec case still needs get_depth below. */
    if (spec.min >= 0 && current_level > spec.max)
        return expr_copy(expr);

    Expr* intermediate = NULL;
    /* Descend (bottom-up) only while a deeper level can still be mapped. When the
     * children are all one level past max, none of them maps, so the node is
     * shared wholesale as `intermediate` -- no per-child recursion and no
     * structural rebuild -- and f (if it applies at this level) wraps the shared
     * subtree. */
    if (expr->type == EXPR_FUNCTION && !(spec.min >= 0 && current_level + 1 > spec.max)) {
        // Recurse first (Bottom-up)
        size_t count = expr->data.function.arg_count;
        Expr** new_args = malloc(sizeof(Expr*) * count);
        for (size_t i = 0; i < count; i++) {
            new_args[i] = map_at_level(f, expr->data.function.args[i], current_level + 1, spec);
        }

        Expr* new_head = NULL;
        if (spec.heads) {
            new_head = map_at_level(f, expr->data.function.head, current_level + 1, spec);
        } else {
            new_head = expr_copy(expr->data.function.head);
        }

        intermediate = expr_new_function(new_head, new_args, count);
        free(new_args);
    } else {
        intermediate = expr_copy(expr);
    }

    bool should_map;
    if (spec.min >= 0) {
        should_map = (current_level >= spec.min && current_level <= spec.max);
    } else {
        int64_t d = get_depth(intermediate);
        should_map = (current_level >= spec.min && current_level <= spec.max) ||
                     (-d >= spec.min && -d <= spec.max);
    }

    if (should_map) {
        /* The f[...] application is left for the evaluator to reduce, so a
         * surrounding Hold suppresses it: Map[f, Hold[1+1]] is Hold[f[1+1]],
         * not Hold[f[2]]. Callers that need the reduced form now -- only the
         * NDArray repack paths -- evaluate the finished tree themselves. */
        Expr* f_copy = expr_copy(f);
        return expr_new_function(f_copy, &intermediate, 1);   /* adopts intermediate */
    }

    return intermediate;
}

/* If `list` (owned) is a rectangular List of machine numbers, repack it into an
 * NDArray of dtype `dt` (freeing `list`), so mapping a numeric function over a
 * packed array yields a packed array again — matching Sin[arr], arr + 1, etc.
 * A symbolic or non-rectangular result is returned unchanged as a List. */
/* Repack a result whose input was the array `src`, inheriting src's dtype AND
 * its presentation -- packed in, packed out; visible in, visible out. Consumes
 * `list` on success, returns it untouched when the result is not packable. */
static Expr* map_try_repack(const Expr* src, Expr* list, NDType dt) {
    if (!list || list->type != EXPR_FUNCTION) return list;
    /* A packed source re-sniffs, for the same reason ndstruct_delist_repack does:
     * f can return a different exactness than it was given, and repacking at the
     * SOURCE dtype would coerce it -- Map[If[# > 150., #, 1] &, Range[1., 300.]]
     * would turn every exact 1 into 1.. A visible NDArray[...] keeps coercing,
     * which is its documented machine-buffer behaviour. */
    if (is_packed_list(src)) return pack_repack_like(src, list);
    Expr* packed = ndarray_from_nested_list_like(src, list, dt);
    if (packed) { expr_free(list); return packed; }
    return list;
}

/* Map f over the leading axis of NDArray `a` at the default level {1}: apply f
 * to each part (a scalar leaf for rank 1, a sub-NDArray row for rank >= 2),
 * collect the evaluated results, and repack them when they are numeric. */
static Expr* map_ndarray_axis(Expr* f, Expr* a) {
    int64_t n = a->data.ndarray.dims[0];
    size_t count = (n > 0) ? (size_t)n : 0;
    Expr** results = malloc(sizeof(Expr*) * (count ? count : 1));
    for (int64_t i = 0; i < n; i++) {
        Expr* idx = expr_new_integer(i + 1);
        Expr* iargs[1] = { idx };
        bool degrade = false;
        Expr* part = ndarray_part(a, iargs, 1, &degrade);   /* owned */
        expr_free(idx);
        if (!part) part = expr_new_symbol(SYM_Null);
        Expr* call = expr_new_function(expr_copy(f), &part, 1);  /* adopts part */
        results[i] = evaluate(call);
        expr_free(call);
    }
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), results, count);
    free(results);
    return map_try_repack(a, list, a->data.ndarray.dtype);
}

Expr* builtin_map(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;

    Expr* f = res->data.function.args[0];
    Expr* expr = res->data.function.args[1];

    /* Map over an association threads over its values, preserving keys:
     * Map[f, <|k -> v|>] -> <|k -> f[v]|>. Only the default level (no explicit
     * level spec) is special-cased; deeper level specs fall through. */
    if (res->data.function.arg_count == 2 && is_association(expr))
        return assoc_map_values(f, expr);

    /* Automatic numeric fast-path: Map[f, {machine numbers}] with f a numeric
     * pure function runs the compiled body per element (no per-element evaluate).
     * NULL means "not numeric-closed" -> the general level-mapper below. */
    if (res->data.function.arg_count == 2) {
        Expr* fast = numloop_map(f, expr);
        if (fast) return fast;
    }

    Expr* ls = (res->data.function.arg_count >= 3) ? res->data.function.args[2] : NULL;
    if (is_rule2(ls)) ls = NULL;            /* option in the level-spec slot */

    LevelSpec spec = parse_level_spec(ls, 1, 1);
    parse_options(res, ls ? 3 : 2, &spec);

    /* NDArray is atomic, so the generic level-mapper would never descend into
     * it. At the default level iterate the leading axis directly; otherwise
     * materialize once and reuse the generic path. Numeric results repack. */
    if (is_ndarray(expr)) {
        if (spec.min == 1 && spec.max == 1 && !spec.heads)
            return map_ndarray_axis(f, expr);
        Expr* nested = ndarray_to_nested_list(expr);
        Expr* out = map_at_level(f, nested, 0, spec);
        expr_free(nested);
        /* map_at_level leaves the f[...] applications unreduced, but repacking
         * needs actual machine numbers -- evaluate before offering it back. */
        return map_try_repack(expr, eval_and_free(out), expr->data.ndarray.dtype);
    }

    return map_at_level(f, expr, 0, spec);
}

/* ------------------- MapIndexed -------------------
 *
 * MapIndexed[f, expr]             applies f to the parts of expr at level 1,
 *                                 handing each part's position to f as a
 *                                 second argument: {f[e1,{1}], f[e2,{2}], ...}
 * MapIndexed[f, expr, levelspec]  the same, over the parts selected by
 *                                 levelspec (default {1}).
 *
 * The second argument is a full position of the shape Part / Extract /
 * Position use: {1,2} for the second part of the first part, {0} for a head
 * under Heads->True, {Key[k]} for an association value, {} for level 0. So
 * Extract[expr, #2] is #1 at every visited node.
 *
 * Traversal is bottom-up: a node's parts are rebuilt before f is applied to
 * the node itself, so MapIndexed[f, {{a}}, 2] gives f[{f[a, {1,1}]}, {1}].
 *
 * The f[...] applications are NOT evaluated here. Wolfram specifies that
 * MapIndexed "always effectively constructs a complete new expression and then
 * evaluates it", and the difference is observable: deferring keeps
 * MapIndexed[f, Hold[1+1]] as Hold[f[1+1, {1}]] instead of collapsing it to
 * Hold[f[2]]. Returning the assembled tree makes the evaluator do exactly what
 * that sentence describes. */

/* One position component per recursion level, threaded through the C stack.
 * `comp` is what this frame contributes to the position (a part index, the
 * Integer 0 of a head, or Key[k] for an association value) and lives in the
 * caller's stack frame; `len` is the total path length at this depth. Keeping
 * the path on the C stack rather than in a heap stack means there is no
 * push/pop bookkeeping for a future early return to get wrong. */
typedef struct MIPath {
    const struct MIPath* prev;
    Expr*  comp;
    size_t len;
} MIPath;

/* The accumulated position as a List. An empty path (level 0) gives {}. */
static Expr* mi_position(const MIPath* path) {
    size_t n = path ? path->len : 0;
    Expr** items = malloc(sizeof(Expr*) * (n ? n : 1));
    for (const MIPath* p = path; p; p = p->prev)
        items[p->len - 1] = expr_copy(p->comp);
    Expr* pos = expr_new_function(expr_new_symbol(SYM_List), items, n);
    free(items);
    return pos;
}

/* True for an Association whose every entry is a two-argument rule.
 * is_association() only checks the head, and builtin_association leaves a
 * malformed Association[a, b] intact, so an entry's args[0]/args[1] must not
 * be read until this has passed — doing so type-puns a symbol's SymbolDef*. */
static bool mi_is_assoc(const Expr* e) {
    if (!is_association(e)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!is_rule2(e->data.function.args[i])) return false;
    return true;
}

/* Rebuild `expr` bottom-up, wrapping each part selected by `spec` in
 * f[part, position].
 *
 * `*out_depth` receives the depth of the ORIGINAL subtree, never the rebuilt
 * one — the rebuilt tree carries f[...] wrappers that would shift every
 * negative level bound. Accumulating it as the recursion unwinds also keeps
 * the whole traversal O(n); recomputing a depth at each visited node (as
 * map_at_level and scan_at_level do) is O(n*d). It counts the head when
 * Heads->True is set, matching Depth[expr, Heads->True]. */
static Expr* mi_at_level(Expr* f, Expr* expr, int64_t level, LevelSpec spec,
                         const MIPath* path, int64_t* out_depth) {
    /* Non-negative spec fast path (mirrors map_at_level/scan_at_level). The depth
     * this function accumulates into *out_depth is consulted ONLY by a negative
     * level bound (see the membership test below), so with both bounds >= 0 it is
     * never used; and once past the maximum level nothing here or deeper is
     * wrapped. So the subtree passes through by an O(1) refcount bump -- no
     * descent, no rebuild, no depth accumulation. This is what stops
     * MapIndexed[f, list] (default level 1) from walking every element's whole
     * subtree; its cost had grown with element depth. Negative-spec calls keep
     * the full traversal, which is the only case that reads the depth. */
    if (spec.min >= 0 && spec.max >= 0 && level > spec.max) {
        *out_depth = 1;                /* unused under a non-negative spec */
        return expr_copy(expr);
    }

    Expr*   rebuilt;
    int64_t maxpart = 0;               /* depth of the deepest original part */
    int64_t depth   = 1;               /* atoms, NDArray, Rational, Complex */
    size_t  plen = path ? path->len : 0;

    if (spec.min >= 0 && spec.max >= 0 && level + 1 > spec.max) {
        /* Deepest level a non-negative spec wraps at: every child (and, with
         * Heads->True, the head) is one level past max, so none is wrapped and a
         * rebuild would just reproduce the original node. Share it instead -- no
         * descent, no rebuild -- and let the membership test below wrap it. This
         * removes the per-element node rebuild on top of the descent already cut
         * by the early return above; depth stays 1 (unused under this spec). */
        rebuilt = expr_copy(expr);
    } else if (mi_is_assoc(expr)) {
        /* An association's parts are its VALUES, positioned by Key[k]: the
         * Rule wrapper is not a level of its own, which is what makes a
         * two-deep association position {Key[a], Key[b]} rather than four
         * components long. The Association head is never mapped. */
        size_t n = expr->data.function.arg_count;
        Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t i = 0; i < n; i++) {
            Expr* entry = expr->data.function.args[i];
            Expr* key   = entry->data.function.args[0];
            Expr* karg[1] = { expr_copy(key) };
            Expr* comp = expr_new_function(expr_new_symbol(SYM_Key), karg, 1);
            MIPath frame = { path, comp, plen + 1 };
            int64_t vd = 1;
            Expr* newv = mi_at_level(f, entry->data.function.args[1],
                                     level + 1, spec, &frame, &vd);
            expr_free(comp);
            if (vd > maxpart) maxpart = vd;
            /* Preserve Rule vs RuleDelayed: <|a :> 1|> survives evaluation as
             * Association[RuleDelayed[a, 1]] and must stay delayed. */
            Expr* rargs[2] = { expr_copy(key), newv };
            out[i] = expr_new_function(expr_copy(entry->data.function.head), rargs, 2);
        }
        rebuilt = expr_new_function(expr_new_symbol(SYM_Association), out, n);
        free(out);
        depth = maxpart + 1;
    } else if (expr->type == EXPR_FUNCTION &&
               !head_is(expr, SYM_Rational) && !head_is(expr, SYM_Complex)) {
        /* Rational and Complex are atomic here, as they are for get_depth and
         * for Level[]: MapIndexed[f, 1/2] is 1/2, not f[1,{1}]/f[2,{2}]. */
        size_t n = expr->data.function.arg_count;
        Expr* newhead;
        if (spec.heads) {
            Expr* comp = expr_new_integer(0);          /* a head is part 0 */
            MIPath frame = { path, comp, plen + 1 };
            int64_t hd = 1;
            newhead = mi_at_level(f, expr->data.function.head,
                                  level + 1, spec, &frame, &hd);
            expr_free(comp);
            if (hd > maxpart) maxpart = hd;
        } else {
            newhead = expr_copy(expr->data.function.head);
        }
        Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t i = 0; i < n; i++) {
            Expr* comp = expr_new_integer((int64_t)(i + 1));
            MIPath frame = { path, comp, plen + 1 };
            int64_t cd = 1;
            out[i] = mi_at_level(f, expr->data.function.args[i],
                                 level + 1, spec, &frame, &cd);
            expr_free(comp);
            if (cd > maxpart) maxpart = cd;
        }
        rebuilt = expr_new_function(newhead, out, n);
        free(out);
        depth = maxpart + 1;
    } else {
        rebuilt = expr_copy(expr);      /* atom, NDArray, Rational, Complex */
    }

    *out_depth = depth;

    /* Standard mixed positive/negative level membership: a bound >= 0 is a
     * level counted from the root; a bound < 0 is a negative depth (leaves at
     * depth 1 sit at level -1). */
    bool lo = (spec.min >= 0) ? (level >= spec.min) : (-depth >= spec.min);
    bool hi = (spec.max >= 0) ? (level <= spec.max) : (-depth <= spec.max);
    if (lo && hi) {
        Expr* fargs[2] = { rebuilt, mi_position(path) };   /* both adopted */
        return expr_new_function(expr_copy(f), fargs, 2);
    }
    return rebuilt;
}

/* MapIndexed over an NDArray at the default level pairs each leading-axis part
 * with its position {i}: {f[part1, {1}], f[part2, {2}], ...}. Parts come from
 * ndarray_part (a scalar for rank 1, a sub-NDArray row for rank >= 2). The
 * f[...] nodes are left for the evaluator to reduce, as in the list case. */
static Expr* mi_ndarray_axis(Expr* f, Expr* arr) {
    int64_t nd = arr->data.ndarray.dims[0];
    size_t cnt = (nd > 0) ? (size_t)nd : 0;
    Expr** out = malloc(sizeof(Expr*) * (cnt ? cnt : 1));
    for (int64_t i = 0; i < nd; i++) {
        Expr* idx = expr_new_integer(i + 1);
        Expr* iargs[1] = { idx };
        bool degrade = false;
        Expr* part = ndarray_part(arr, iargs, 1, &degrade);   /* owned */
        expr_free(idx);
        if (!part) part = expr_new_symbol(SYM_Null);
        Expr* pos_inner = expr_new_integer(i + 1);
        Expr* pos = expr_new_function(expr_new_symbol(SYM_List), &pos_inner, 1);
        Expr* fargs[2] = { part, pos };   /* both adopted */
        out[i] = expr_new_function(expr_copy(f), fargs, 2);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, cnt);
    free(out);
    return result;
}

/* Validate the trailing arguments: everything from `start` on must be an
 * option rule naming a known option. Reports the rightmost offender via
 * ::nonopt / ::optx and returns false so the call is left unevaluated. */
static bool mi_check_options(Expr* res, size_t start) {
    size_t argc = res->data.function.arg_count;
    Expr* bad_shape = NULL;
    Expr* bad_name  = NULL;
    for (size_t i = start; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (!is_rule2(opt)) { bad_shape = opt; continue; }
        Expr* name = opt->data.function.args[0];
        if (!(name->type == EXPR_SYMBOL && name->data.symbol.name == SYM_Heads))
            bad_name = name;
    }
    if (!bad_shape && !bad_name) return true;

    char* call = expr_to_string(res);
    char* bad  = expr_to_string(bad_shape ? bad_shape : bad_name);
    if (bad_shape)
        fprintf(stderr,
                "MapIndexed::nonopt: Options expected (instead of %s) beyond "
                "position %zu in %s. An option must be a rule or a list of "
                "rules.\n", bad ? bad : "?", start, call ? call : "?");
    else
        fprintf(stderr, "MapIndexed::optx: Unknown option %s in %s.\n",
                bad ? bad : "?", call ? call : "?");
    free(bad);
    free(call);
    return false;
}

Expr* builtin_mapindexed(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return builtin_arg_error("MapIndexed", 0, 1, 3); /* ::argb */
    if (argc == 1) return NULL;                  /* MapIndexed[f]: leave alone */

    Expr* f    = res->data.function.args[0];
    Expr* expr = res->data.function.args[1];

    Expr* ls = (argc >= 3) ? res->data.function.args[2] : NULL;
    if (is_rule2(ls)) ls = NULL;            /* option in the level-spec slot */
    size_t opt_start = ls ? 3 : 2;

    if (!mi_check_options(res, opt_start)) return NULL;

    LevelSpec spec;
    if (!parse_level_spec_strict(ls, &spec)) return NULL;  /* bad level spec */
    parse_options(res, opt_start, &spec);                  /* Heads->True */

    /* An empty level range (0, {3,1}, ...) selects nothing. Returning the
     * expression itself keeps a packed NDArray packed and skips a full
     * rebuild that could only produce an identical tree.
     *
     * Only two non-negative bounds are comparable this way: a negative bound
     * is a depth, not a level, so {1,-1} (the spelling of the level spec -1)
     * spans everything at level 1 or below despite reading as min > max. */
    if (spec.min >= 0 && spec.max >= 0 && spec.min > spec.max) {
        res->data.function.args[1] = NULL;   /* transfer ownership; SPEC.md §4 */
        return expr;
    }

    /* NDArray is an atomic value, so the generic traversal would never descend
     * into it. Iterate the leading axis directly at the default level;
     * otherwise materialize once and reuse the generic path. */
    if (is_ndarray(expr)) {
        if (spec.min == 1 && spec.max == 1 && !spec.heads)
            return mi_ndarray_axis(f, expr);
        Expr* nested = ndarray_to_nested_list(expr);
        int64_t depth = 1;
        Expr* out = mi_at_level(f, nested, 0, spec, NULL, &depth);
        expr_free(nested);
        return out;
    }

    int64_t depth = 1;
    return mi_at_level(f, expr, 0, spec, NULL, &depth);
}

/* ------------------- MapAt ------------------- */

/* Leaf action for MapAt: build f[part]. The application is left for the
 * evaluator to reduce, so a surrounding Hold suppresses it the Wolfram way
 * (MapAt[f, Hold[1+1], 1] is Hold[f[1+1]], not Hold[f[2]]). `ctx` is f. */
static Expr* mapat_leaf(void* ctx, Expr* leaf) {
    Expr* arg_copy = expr_copy(leaf);
    return expr_new_function(expr_copy((Expr*)ctx), &arg_copy, 1);
}

/*
 * MapAt[f, expr, pos]        apply f at one position, or at each of several.
 * MapAt[f, pos]              operator form: MapAt[f, pos][expr].
 *
 * Position resolution -- integers (negatives from the end, 0 for the head),
 * All, Span, and Key[k]/literal keys on associations -- lives in the shared
 * walker expr_apply_at_path (src/part.c), which ReplaceAt uses too. A position
 * that does not exist makes the walker return NULL and MapAt stay unevaluated,
 * matching Part[{a, b, c}, 5].
 */
Expr* builtin_map_at(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 2) {
        /* Operator form: MapAt[f, pos] -> Function[MapAt[f, #1, pos]]. */
        Expr* slot_args[1] = { expr_new_integer(1) };
        Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot), slot_args, 1);
        Expr* inner_args[3] = { expr_copy(res->data.function.args[0]), slot,
                                expr_copy(res->data.function.args[1]) };
        Expr* inner = expr_new_function(expr_new_symbol(SYM_MapAt), inner_args, 3);
        Expr* func_args[1] = { inner };
        return expr_new_function(expr_new_symbol(SYM_Function), func_args, 1);
    }
    if (argc != 3) return NULL;

    Expr* f    = res->data.function.args[0];
    Expr* expr = res->data.function.args[1];
    Expr* pos  = res->data.function.args[2];

    /* An NDArray is atomic, so the walker cannot descend into it. Materialize
     * to a nested list, map, and repack if the result is still numeric --
     * applying a symbolic f leaves a plain List, as in MapAll. */
    if (is_ndarray(expr)) {
        Expr* nested = ndarray_to_nested_list(expr);
        Expr* out = expr_apply_at_positions(nested, pos, mapat_leaf, f);
        expr_free(nested);
        if (!out) return NULL;
        /* mapat_leaf leaves f[part] unreduced, but repacking needs actual
         * machine numbers -- evaluate before offering the result back. */
        return map_try_repack(expr, eval_and_free(out), expr->data.ndarray.dtype);
    }

    return expr_apply_at_positions(expr, pos, mapat_leaf, f);
}

/* ------------------- MapAll ------------------- */

Expr* builtin_map_all(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;

    Expr* f = res->data.function.args[0];
    Expr* expr = res->data.function.args[1];
    
    LevelSpec spec = {0, LEVEL_SPEC_INF, false}; // {0, Infinity}
    parse_options(res, 2, &spec);

    /* NDArray is atomic; materialize to a nested list so f is applied at every
     * level (including the whole array at level 0), as MapAll does for lists. */
    if (is_ndarray(expr)) {
        Expr* nested = ndarray_to_nested_list(expr);
        Expr* out = map_at_level(f, nested, 0, spec);
        expr_free(nested);
        return out;
    }

    return map_at_level(f, expr, 0, spec);
}

/* ------------------- compiled predicates -------------------
 *
 * Seven heads apply a boolean test to every element -- Select, AllTrue,
 * AnyTrue, NoneTrue, TakeWhile, LengthWhile, SelectFirst -- and every one of
 * them opened with ndstruct_delist_repack: materialise the buffer into one Expr
 * per element, then call the test through the interpreter for each. Map has had
 * a compiled path since the auto-compilation experiment (numloop_map); the
 * boolean half was never written, and it left these the slowest numeric heads
 * in the system. AllTrue[v, # > 0 &] over 10^6 measured 416 ms against
 * np.all(v > 0)'s 319 us.
 *
 * numloop_pred_compile turns an ordinary comparison predicate into bytecode;
 * these three helpers are the loop shapes the seven heads need. Everything
 * outside the compilable subset returns NULL here and the caller runs exactly
 * the code it ran before.
 *
 * FLOAT64 ONLY, and that is a correctness bound rather than an omission: the
 * compiled comparison reproduces the interpreter's INEXACT rule, which carries
 * a relative tolerance of 2^-46. Two exact Integers compare through GMP with no
 * tolerance at all, so an int64 buffer must not come down here. See the
 * contract note in numloop.h. */

/* The common opening: `v` a rank-1 float64 buffer and `f` a compilable
 * predicate. On success the caller owns *pred and must numloop_pred_free it. */
static bool pred_open(const Expr* f, const Expr* v, NumPred** pred,
                      const double** data, size_t* n) {
    if (!is_packed_list(v) && !is_ndarray(v)) return false;
    if (v->data.ndarray.rank != 1) return false;
    if (v->data.ndarray.dtype != NDT_FLOAT64) return false;
    size_t len = (size_t)v->data.ndarray.dims[0];
    if (len == 0) return false;            /* trivial; let the interpreter do it */
    NumPred* p = numloop_pred_compile(f);
    if (!p) return false;
    *pred = p;
    *data = (const double*)v->data.ndarray.data;
    *n = len;
    return true;
}

/* AllTrue / AnyTrue / NoneTrue. `want` is the value that short-circuits.
 * Returns NULL to fall back; otherwise True or False. */
static Expr* pred_quantify(const Expr* f, const Expr* v, int mode) {
    NumPred* p; const double* d; size_t n;
    if (!pred_open(f, v, &p, &d, &n)) return NULL;
    bool fired = false, ok = true;
    for (size_t i = 0; i < n; i++) {
        bool t;
        if (!numloop_pred_run(p, d[i], &t)) { ok = false; break; }
        /* mode 0 AllTrue: a False decides it. mode 1 AnyTrue and mode 2
         * NoneTrue: a True decides it. */
        if (mode == 0 ? !t : t) { fired = true; break; }
    }
    numloop_pred_free(p);
    if (!ok) return NULL;                  /* non-finite -> interpreter answers */
    bool val;
    if (mode == 0)      val = !fired;      /* All: True unless a False appeared */
    else if (mode == 1) val = fired;       /* Any: True iff a True appeared */
    else                val = !fired;      /* None: True unless a True appeared */
    return expr_new_symbol(val ? SYM_True : SYM_False);
}

/* Select[list, crit] with no count limit. */
static Expr* pred_select(const Expr* f, const Expr* v) {
    NumPred* p; const double* d; size_t n;
    if (!pred_open(f, v, &p, &d, &n)) return NULL;

    unsigned char* keep = malloc(n);
    if (!keep) { numloop_pred_free(p); return NULL; }
    size_t kept = 0;
    bool ok = true;
    for (size_t i = 0; i < n; i++) {
        bool t;
        if (!numloop_pred_run(p, d[i], &t)) { ok = false; break; }
        keep[i] = (unsigned char)t;
        if (t) kept++;
    }
    numloop_pred_free(p);
    if (!ok) { free(keep); return NULL; }

    /* DERIVED from a buffer, so no threshold: a filtered packed list stays
     * packed however few elements survive, matching Sin[packedList]. If
     * ndbuild_open_like declines (rank/dtype it will not carry) the plain List
     * below is built instead, which is what Select produced before. */
    int64_t dk = (int64_t)kept;
    void* vb = NULL;
    Expr* packed = kept ? ndbuild_open_like(v, 1, &dk, NDT_FLOAT64, &vb) : NULL;
    if (packed) {
        double* out = (double*)vb;
        size_t k = 0;
        for (size_t i = 0; i < n; i++) if (keep[i]) out[k++] = d[i];
        free(keep);
        return packed;
    }
    Expr** args = kept ? malloc(kept * sizeof(Expr*)) : NULL;
    if (kept && !args) { free(keep); return NULL; }
    size_t k = 0;
    for (size_t i = 0; i < n; i++) if (keep[i]) args[k++] = expr_new_real(d[i]);
    free(keep);
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), args, kept);
    free(args);
    return r;
}

/* The first `k` elements of a rank-1 float64 buffer, as a buffer. A leading
 * slice is contiguous, so this is one memcpy -- the whole reason TakeWhile can
 * answer without materialising. */
static Expr* pred_leading_slice(const Expr* v, size_t k) {
    const double* d = (const double*)v->data.ndarray.data;
    int64_t dk = (int64_t)k;
    void* vb = NULL;
    Expr* packed = k ? ndbuild_open_like(v, 1, &dk, NDT_FLOAT64, &vb) : NULL;
    if (packed) {
        memcpy(vb, d, k * sizeof(double));
        return packed;
    }
    Expr** args = k ? malloc(k * sizeof(Expr*)) : NULL;
    if (k && !args) return NULL;
    for (size_t i = 0; i < k; i++) args[i] = expr_new_real(d[i]);
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), args, k);
    free(args);
    return r;
}

/* pred_find_first's two non-index answers. PRED_FALL_BACK keeps the caller's
 * old path; PRED_NO_MATCH means the scan ran and nothing satisfied it, which is
 * a real answer and not a decline. */
#define PRED_FALL_BACK (-1)
#define PRED_NO_MATCH  (-2)

/* Index of the first element satisfying the predicate, or one of the two
 * sentinels above. */
static int64_t pred_find_first(const Expr* f, const Expr* v) {
    NumPred* p; const double* d; size_t n;
    if (!pred_open(f, v, &p, &d, &n)) return PRED_FALL_BACK;
    int64_t at = PRED_NO_MATCH;
    bool ok = true;
    for (size_t i = 0; i < n; i++) {
        bool t;
        if (!numloop_pred_run(p, d[i], &t)) { ok = false; break; }
        if (t) { at = (int64_t)i; break; }
    }
    numloop_pred_free(p);
    return ok ? at : PRED_FALL_BACK;
}

/* The length of the leading run satisfying the predicate. Returns -1 to fall
 * back, which TakeWhile and LengthWhile both read as "run the old path". */
static int64_t pred_run_length(const Expr* f, const Expr* v) {
    NumPred* p; const double* d; size_t n;
    if (!pred_open(f, v, &p, &d, &n)) return -1;
    size_t i = 0;
    bool ok = true;
    for (; i < n; i++) {
        bool t;
        if (!numloop_pred_run(p, d[i], &t)) { ok = false; break; }
        if (!t) break;
    }
    numloop_pred_free(p);
    return ok ? (int64_t)i : -1;
}

/* ------------------- Select ------------------- */

Expr* builtin_select(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count != 2 && res->data.function.arg_count != 3)) {
        return NULL;
    }
    
    Expr* list = res->data.function.args[0];
    Expr* crit = res->data.function.args[1];

    int64_t n_max = -1; // -1 means all
    if (res->data.function.arg_count == 3) {
        Expr* n_expr = res->data.function.args[2];
        if (n_expr->type == EXPR_INTEGER) {
            n_max = n_expr->data.integer;
        } else {
            return NULL; // Invalid n
        }
    }

    /* Select over an association filters by value, preserving keys:
     * Select[<|k -> v|>, p] keeps the entries where p[v] is True;
     * Select[<|k -> v|>, p, n] keeps the first n such entries. */
    if (is_association(list))
        return assoc_select_values(crit, list, n_max);

    /* Compiled-predicate fast path (numloop.h). Falls through to the
     * materialising path below for any predicate outside the subset. */
    if (n_max < 0) {
        Expr* fast = pred_select(crit, list);
        if (fast) return fast;
    }

    /* An NDArray is atomic, so the element walk below looks straight past one
     * and the call comes back UNEVALUATED, while the identical List call works.
     * Materialise, reuse the List implementation, repack — see ndstruct.h. */
    if (is_ndarray(list)) return ndstruct_delist_repack(res, list);

    if (list->type != EXPR_FUNCTION) return NULL; // Can only select from compound expressions
    
    size_t count = list->data.function.arg_count;
    Expr** kept_args = NULL;
    if (count > 0) kept_args = malloc(sizeof(Expr*) * count);
    size_t kept_count = 0;
    
    for (size_t i = 0; i < count; i++) {
        if (n_max >= 0 && (int64_t)kept_count >= n_max) break;
        
        Expr* elem = list->data.function.args[i];
        
        Expr* call_args[1] = { expr_copy(elem) };
        Expr* call = expr_new_function(expr_copy(crit), call_args, 1);
        
        Expr* eval_res = evaluate(call);
        
        if (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol.name == SYM_True) {
            kept_args[kept_count++] = expr_copy(elem);
        }
        
        expr_free(eval_res);
        expr_free(call);
    }
    
    Expr* result = expr_new_function(expr_copy(list->data.function.head), kept_args, kept_count);
    if (kept_args) free(kept_args);

    /* The kept elements are copies of the input's, already evaluated, so
     * offering the filtered list is safe -- and it is the one Select case worth
     * packing, where a large numeric list is being narrowed to another one. */
    return pack_offer(result);
}

/* ------------------- TakeWhile / LengthWhile -------------------
 *
 * TakeWhile[list, crit]    the longest leading run of elements e with crit[e] === True
 * LengthWhile[list, crit]  the length of that run
 *
 * Over an association the criterion is applied to the values; TakeWhile keeps
 * the matching leading entries as an association (keys preserved), LengthWhile
 * counts them. */

/* Length of the longest leading run for which crit[e] evaluates to True. When
 * `over_values` the elements are Rule[k, v] nodes and the test uses v. */
static size_t leading_run_length(Expr* crit, Expr** elems, size_t n, bool over_values) {
    size_t k = 0;
    for (; k < n; k++) {
        Expr* e = over_values ? elems[k]->data.function.args[1] : elems[k];
        Expr* call_args[1] = { expr_copy(e) };
        Expr* call = expr_new_function(expr_copy(crit), call_args, 1);
        Expr* v = evaluate(call);
        expr_free(call);
        bool ok = (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_True);
        expr_free(v);
        if (!ok) break;
    }
    return k;
}

Expr* builtin_takewhile(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* coll = res->data.function.args[0];
    Expr* crit = res->data.function.args[1];
    /* Compiled-predicate fast path: the leading run is found on the buffer and
     * Take does the copy, so neither the scan nor the result materialises. */
    {
        int64_t k = pred_run_length(crit, coll);
        if (k >= 0) {
            Expr* r = pred_leading_slice(coll, (size_t)k);
            if (r) return r;
        }
    }
    /* An NDArray is atomic, so the element walk below looks straight past one
     * and the call comes back UNEVALUATED, while the identical List call works.
     * Materialise, reuse the List implementation, repack — see ndstruct.h. */
    if (is_ndarray(coll)) return ndstruct_delist_repack(res, coll);
    bool assoc = is_association(coll);
    if (!assoc && coll->type != EXPR_FUNCTION) return NULL;

    size_t n = coll->data.function.arg_count;
    Expr** elems = coll->data.function.args;
    size_t k = leading_run_length(crit, elems, n, assoc);

    Expr** out = malloc(sizeof(Expr*) * (k ? k : 1));
    for (size_t i = 0; i < k; i++) out[i] = expr_copy(elems[i]);
    Expr* head = assoc ? expr_new_symbol(SYM_Association)
                       : expr_copy(coll->data.function.head);
    Expr* result = expr_new_function(head, out, k);
    free(out);
    return result;
}

Expr* builtin_lengthwhile(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* coll = res->data.function.args[0];
    Expr* crit = res->data.function.args[1];
    /* Compiled-predicate fast path: LengthWhile needs only the count, so
     * nothing is built at all. */
    {
        int64_t k = pred_run_length(crit, coll);
        if (k >= 0) return expr_new_integer(k);
    }
    /* An NDArray is atomic, so the element walk below looks straight past one
     * and the call comes back UNEVALUATED, while the identical List call works.
     * Materialise, reuse the List implementation, repack — see ndstruct.h. */
    if (is_ndarray(coll)) return ndstruct_delist_repack(res, coll);
    bool assoc = is_association(coll);
    if (!assoc && coll->type != EXPR_FUNCTION) return NULL;

    size_t k = leading_run_length(crit, coll->data.function.args,
                                  coll->data.function.arg_count, assoc);
    return expr_new_integer((int64_t)k);
}

/* ------------------- SelectFirst -------------------
 *
 * SelectFirst[list, pred]           first e with pred[e] === True, else Missing["NotFound"]
 * SelectFirst[list, pred, default]  default when none match
 * Over an association, tests the values and returns the first matching value. */
Expr* builtin_select_first(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2 && argc != 3) return NULL;

    Expr* coll = res->data.function.args[0];
    Expr* pred = res->data.function.args[1];
    Expr* deflt = (argc == 3) ? res->data.function.args[2] : NULL;

    /* Compiled-predicate fast path. Cheap when a match comes early, but a
     * predicate nothing satisfies is a full scan and paid the interpreted
     * per-element cost like the rest of the family. */
    {
        int64_t at = pred_find_first(pred, coll);
        if (at >= 0)
            return expr_new_real(
                ((const double*)coll->data.ndarray.data)[at]);
        if (at == PRED_NO_MATCH) {
            if (deflt) return expr_copy(deflt);
            Expr* ma[1] = { expr_new_string("NotFound") };
            return expr_new_function(expr_new_symbol(SYM_Missing), ma, 1);
        }
    }

    /* An NDArray is atomic, so the element walk below looks straight past one
     * and the call comes back UNEVALUATED, while the identical List call works.
     * Materialise, reuse the List implementation, repack — see ndstruct.h. */
    if (is_ndarray(coll)) return ndstruct_delist_repack(res, coll);
    if (is_association(coll)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }
    if (coll->type != EXPR_FUNCTION) return NULL;

    size_t n = coll->data.function.arg_count;
    for (size_t i = 0; i < n; i++) {
        Expr* elem = coll->data.function.args[i];
        Expr* call_args[1] = { expr_copy(elem) };
        Expr* call = expr_new_function(expr_copy(pred), call_args, 1);
        Expr* v = evaluate(call);
        expr_free(call);
        /* A Throw from the predicate propagates. */
        if (eval_is_inflight_throw(v)) return v;
        bool is_true = (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_True);
        expr_free(v);
        if (is_true) return expr_copy(elem);
    }
    if (deflt) return expr_copy(deflt);
    Expr* margs[1] = { expr_new_string("NotFound") };
    return expr_new_function(expr_new_symbol(SYM_Missing), margs, 1);
}

/* ------------------- Catch / Throw -------------------
 *
 * Throw[value] / Throw[value, tag] / Throw[value, tag, f] stop evaluation and
 * hand the Throw[...] node up to the nearest enclosing Catch. Throw is
 * non-held: value/tag/f are evaluated by the arg loop before the throw
 * propagates. The plain Throw[...] node *is* the in-flight sentinel
 * (eval_is_inflight_throw); evaluate_step's argument loop short-circuits on it,
 * so this builtin only validates arity and declines. */
Expr* builtin_throw(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n < 1 || n > 3) return NULL;   /* wrong arity: leave unevaluated */
    return NULL;                        /* Throw[...] stands as the sentinel */
}

/* ------------------- Goto / Label -------------------
 *
 * Goto[tag] transfers control to Label[tag] in the CompoundExpression it
 * appears in directly, then in enclosing ones. Goto is non-held: tag is
 * evaluated by the arg loop, and the plain Goto[tag] node *is* the in-flight
 * sentinel (eval_is_inflight_goto). evaluate_step's argument loop
 * short-circuits on it so a Goto fired inside a nested call bubbles up to the
 * enclosing CompoundExpression, which scans for the matching Label and jumps
 * (builtin_compoundexpression). Like builtin_throw, this only validates arity
 * and declines. */
Expr* builtin_goto(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1) return NULL;  /* wrong arity: leave alone */
    return NULL;                                          /* Goto[tag] stands as the sentinel */
}

/* Label[tag] marks a jump target inside a CompoundExpression. As a statement it
 * is a no-op that evaluates to Null; the raw held Label[tag] node is what
 * builtin_compoundexpression scans for when consuming a Goto[tag] sentinel. */
Expr* builtin_label(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1) return NULL;  /* wrong arity: leave alone */
    return expr_new_symbol(SYM_Null);
}

/* Catch[expr] returns the argument of the first Throw generated while
 * evaluating expr, or expr itself if none. Catch[expr, form] catches only a
 * Throw[v, tag] whose tag matches form (tag is re-evaluated before each
 * comparison); other throws propagate to an outer Catch. Catch[expr, form, f]
 * returns f[value, tag]. Attribute HoldFirst: the body is held so we drive its
 * evaluation ourselves and consume any sentinel; form and f (args 1,2) are
 * evaluated normally by the arg loop. */
Expr* builtin_catch(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return NULL;

    Expr* result = evaluate(res->data.function.args[0]);   /* HoldFirst body */

    if (!eval_is_inflight_throw(result))
        return result;   /* no throw -> Catch yields the evaluated body */

    size_t targc = result->data.function.arg_count;
    Expr* value = result->data.function.args[0];
    Expr* tag   = (targc >= 2) ? result->data.function.args[1] : NULL;

    /* 1-arg Catch catches ANY throw; tag/handler are irrelevant. */
    if (argc == 1) {
        Expr* out = expr_copy(value);
        expr_free(result);
        return out;
    }

    /* 2/3-arg Catch: only a *tagged* throw is eligible, and only if its tag
     * matches form. A tagless Throw[value] is never caught here. */
    bool matched = false;
    if (targc >= 2) {
        Expr* tag_eval = evaluate(expr_copy(tag));   /* WL: tag re-evaluated */
        MatchEnv* env = env_new();
        matched = match(tag_eval, res->data.function.args[1], env);
        env_free(env);
        expr_free(tag_eval);
    }

    if (!matched)
        return result;   /* propagate the same sentinel to an outer Catch */

    if (argc == 2) {
        Expr* out = expr_copy(value);
        expr_free(result);
        return out;
    }

    /* argc == 3: return f[value, tag] (evaluated). tag is non-NULL because a
     * match required targc >= 2. */
    Expr* handler = res->data.function.args[2];   /* evaluated by arg loop */
    Expr* fa[2] = { expr_copy(value), expr_copy(tag) };
    Expr* call = expr_new_function(expr_copy(handler), fa, 2);
    expr_free(result);
    Expr* out = evaluate(call);
    expr_free(call);
    return out;
}

/* ------------------- Scan -------------------
 *
 * Scan[f, expr] / Scan[f, expr, levelspec] applies f to each part of expr
 * selected by levelspec (default {1}) for its side effects and returns Null,
 * discarding f's results. Traversal is depth-first, leaves before roots
 * (post-order). Supports the Heads->True option, early exit via Throw
 * (propagates to an enclosing Catch) and Return (Return[ret] makes the final
 * Scan value ret). Over an association at the default level, scans the values.
 *
 * These helpers mirror Map's LevelSpec/get_depth/parse_options machinery above.
 * The level spec goes through parse_level_spec_strict (shared with MapIndexed)
 * rather than the lenient parse_level_spec, so Scan accepts Infinity in either
 * slot of a two-element spec (e.g. {0, Infinity}) and leaves the call
 * unevaluated on a spec it does not recognise. */

/* Apply f to `part` for its side effects. CONSUMES `part` (adopts ownership).
 * Returns NULL to continue scanning, or non-NULL as the value builtin_scan
 * must return: an in-flight Throw sentinel, a propagating Return[v, h], or the
 * payload of a Return consumed at the Scan boundary. */
static Expr* scan_apply(Expr* f, Expr* part) {
    Expr* call = expr_new_function(expr_copy(f), &part, 1);  /* adopts part */
    Expr* r = evaluate(call);
    expr_free(call);                                         /* frees f-copy + part */
    if (eval_is_inflight_throw(r)) return r;                 /* propagate Throw */
    Expr* rv = NULL;
    EvalReturnAction ra = eval_classify_return(r, SYM_Scan, &rv);
    if (ra == EVAL_RETURN_CONSUME)   { expr_free(r); return rv; }  /* Return[ret] */
    if (ra == EVAL_RETURN_PROPAGATE) return r;                     /* Return[v, h] */
    if (r) expr_free(r);                                    /* result discarded */
    return NULL;
}

/* Depth-first, leaves-before-roots traversal. Visits sub-parts (and, with
 * Heads->True, the head) before the node itself. Returns NULL to continue, or
 * the sentinel/return value to stop and hand back to builtin_scan. */
static Expr* scan_at_level(Expr* f, Expr* expr, int64_t current_level, LevelSpec spec) {
    /* Non-negative upper bound: once past it, neither this node nor any deeper
     * one can be in range (a visit needs current_level <= max, and the level only
     * rises with depth), so stop descending. This is what stops the common
     * Scan[f, list] (default level 1) from walking every element's whole subtree
     * -- the traversal + per-node get_depth cost was growing with element depth,
     * though f only ever fires at level 1. Mirrors the Map/Apply short-circuits.
     * A negative max keeps the full descent (a deeper leaf may still match). */
    if (spec.max >= 0 && current_level > spec.max)
        return NULL;

    if (expr->type == EXPR_FUNCTION) {
        if (spec.heads) {
            Expr* s = scan_at_level(f, expr->data.function.head, current_level + 1, spec);
            if (s) return s;
        }
        size_t count = expr->data.function.arg_count;
        for (size_t i = 0; i < count; i++) {
            Expr* s = scan_at_level(f, expr->data.function.args[i], current_level + 1, spec);
            if (s) return s;
        }
    }

    /* Standard mixed positive/negative level membership: a bound >= 0 is a
     * level counted from the root; a bound < 0 is a negative depth (leaves at
     * depth 1 have negative level -1). get_depth is only needed when a bound is
     * negative; a pure non-negative spec (the common case) skips its full-subtree
     * traversal entirely. */
    bool lo, hi;
    if (spec.min >= 0 && spec.max >= 0) {
        lo = current_level >= spec.min;
        hi = current_level <= spec.max;
    } else {
        int64_t d = get_depth(expr);
        lo = (spec.min >= 0) ? (current_level >= spec.min) : (-d >= spec.min);
        hi = (spec.max >= 0) ? (current_level <= spec.max) : (-d <= spec.max);
    }
    if (lo && hi)
        return scan_apply(f, expr_copy(expr));   /* visit this node last */
    return NULL;
}

Expr* builtin_scan(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return builtin_arg_error("Scan", 0, 1, 3);  /* Scan::argb */
    if (argc < 2) return NULL;                                  /* Scan[f]: leave alone */

    Expr* f = res->data.function.args[0];
    Expr* expr = res->data.function.args[1];

    Expr* ls = (argc >= 3) ? res->data.function.args[2] : NULL;
    if (is_rule2(ls)) ls = NULL;                               /* option in slot 2 */

    LevelSpec spec;
    if (!parse_level_spec_strict(ls, &spec)) return NULL;      /* unrecognized levelspec */
    parse_options(res, ls ? 3 : 2, &spec);                     /* Heads->True */

    /* Association at the default level: scan the values (mirrors Map). */
    if (ls == NULL && is_association(expr)) {
        size_t n = expr->data.function.arg_count;
        for (size_t i = 0; i < n; i++) {
            Expr* elem = expr->data.function.args[i];
            Expr* val = (elem->type == EXPR_FUNCTION && elem->data.function.arg_count == 2)
                        ? elem->data.function.args[1] : elem;
            Expr* s = scan_apply(f, expr_copy(val));
            if (s) return s;
        }
        return expr_new_symbol(SYM_Null);
    }

    /* Automatic numeric fast-path at the default level, mirroring Map's. A
     * numeric-closed body has no side effects for Scan to produce, so all this
     * removes is the per-element substitute+evaluate; it still runs the body so
     * that a non-finite element hands the call back to the interpreter. */
    if (ls == NULL && spec.min == 1 && spec.max == 1 && !spec.heads &&
        !is_ndarray(expr)) {
        Expr* fast = numloop_scan(f, expr);
        if (fast) return fast;
    }

    /* NDArray fast path: NDArray is an atomic value, so the generic traversal
     * would never descend into it. At the default level iterate the leading
     * axis directly (scalar leaf for rank 1, sub-NDArray row for rank >= 2);
     * for any other level spec, materialize once and reuse the generic path. */
    if (is_ndarray(expr)) {
        if (spec.min == 1 && spec.max == 1 && !spec.heads) {
            int64_t n = expr->data.ndarray.dims[0];
            for (int64_t i = 1; i <= n; i++) {
                Expr* idx = expr_new_integer(i);
                Expr* iargs[1] = { idx };
                bool degrade = false;
                Expr* part = ndarray_part(expr, iargs, 1, &degrade);  /* owned */
                expr_free(idx);
                if (!part) continue;
                Expr* s = scan_apply(f, part);   /* consumes part */
                if (s) return s;
            }
            return expr_new_symbol(SYM_Null);
        }
        Expr* nested = ndarray_to_nested_list(expr);
        Expr* s = scan_at_level(f, nested, 0, spec);
        expr_free(nested);
        return s ? s : expr_new_symbol(SYM_Null);
    }

    Expr* s = scan_at_level(f, expr, 0, spec);
    return s ? s : expr_new_symbol(SYM_Null);
}

/* ------------------- AllTrue / AnyTrue / NoneTrue -------------------
 *
 * mode 0 = AllTrue  (True iff test[e] is True for every element)
 * mode 1 = AnyTrue  (True iff test[e] is True for some element)
 * mode 2 = NoneTrue (True iff test[e] is True for no element)
 *
 * Over an association these test the values (via assoc_apply_over_values).
 * If any test result is neither True nor False the whole call is left
 * unevaluated (return NULL), matching Wolfram; short-circuits otherwise. */
static Expr* all_any_none_true(Expr* res, int mode) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* coll = res->data.function.args[0];
    Expr* test = res->data.function.args[1];

    /* Bool-buffer fast path: AllTrue/AnyTrue/NoneTrue over a bool array with the
     * identity-shaped predicate (TrueQ or Identity, both of which return a bool
     * element unchanged) is a single early-exit scan of the byte buffer -- np.all/
     * np.any, and the reason a bool dtype closes the sign-predicate -> quantifier
     * pipeline without materialising 10^6 True/False symbols. Any other predicate
     * or dtype falls through to the compiled/materialising paths below. */
    if (is_ndarray(coll) && coll->data.ndarray.dtype == NDT_BOOL &&
        test->type == EXPR_SYMBOL &&
        (test->data.symbol.name == SYM_Identity ||
         strcmp(test->data.symbol.name, "TrueQ") == 0)) {
        size_t n = ndarray_size(coll);
        const unsigned char* b = (const unsigned char*)coll->data.ndarray.data;
        if (mode == 0) {                    /* AllTrue: every byte set? */
            for (size_t i = 0; i < n; i++)
                if (!b[i]) return expr_new_symbol(SYM_False);
            return expr_new_symbol(SYM_True);
        }
        bool any = false;                   /* AnyTrue / NoneTrue: any byte set? */
        for (size_t i = 0; i < n; i++) if (b[i]) { any = true; break; }
        if (mode == 1) return expr_new_symbol(any ? SYM_True : SYM_False);
        return expr_new_symbol(any ? SYM_False : SYM_True);      /* NoneTrue */
    }

    /* Compiled-predicate fast path (numloop.h) before the materialising one:
     * this is the head the sweep measured at 416 ms against np.all's 319 us. */
    { Expr* fast = pred_quantify(test, coll, mode); if (fast) return fast; }

    /* An NDArray is atomic, so the element walk below looks straight past one
     * and the call comes back UNEVALUATED, while the identical List call works.
     * Materialise, reuse the List implementation, repack — see ndstruct.h. */
    if (is_ndarray(coll)) return ndstruct_delist_repack(res, coll);
    if (is_association(coll)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }
    if (coll->type != EXPR_FUNCTION) return NULL;

    size_t n = coll->data.function.arg_count;
    bool indeterminate = false;
    for (size_t i = 0; i < n; i++) {
        Expr* call_args[1] = { expr_copy(coll->data.function.args[i]) };
        Expr* call = expr_new_function(expr_copy(test), call_args, 1);
        Expr* v = evaluate(call);
        expr_free(call);
        /* A Throw from the test propagates. */
        if (eval_is_inflight_throw(v)) return v;
        bool is_true  = (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_True);
        bool is_false = (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_False);
        expr_free(v);

        if (is_true) {
            if (mode == 1) return expr_new_symbol(SYM_True);   /* AnyTrue short-circuit */
            if (mode == 2) return expr_new_symbol(SYM_False);  /* NoneTrue short-circuit */
        } else if (is_false) {
            if (mode == 0) return expr_new_symbol(SYM_False);  /* AllTrue short-circuit */
        } else {
            indeterminate = true;
        }
    }
    if (indeterminate) return NULL;  /* leave unevaluated */
    /* No short-circuit fired: All -> True, Any -> False, None -> True. */
    return expr_new_symbol(mode == 1 ? SYM_False : SYM_True);
}

Expr* builtin_all_true(Expr* res)  { return all_any_none_true(res, 0); }
Expr* builtin_any_true(Expr* res)  { return all_any_none_true(res, 1); }
Expr* builtin_none_true(Expr* res) { return all_any_none_true(res, 2); }

/* ------------------- Through ------------------- */

static Expr* transform_head(Expr* head, Expr* h_spec, Expr** args, size_t arg_count, bool* transformed) {
    if (head->type != EXPR_FUNCTION) return expr_copy(head);
    
    Expr* P = head->data.function.head;
    bool match = false;
    
    if (h_spec == NULL) {
        match = true;
    } else {
        match = expr_eq(P, h_spec);
    }
    
    if (match) {
        *transformed = true;
        size_t n = head->data.function.arg_count;
        Expr** new_args = malloc(sizeof(Expr*) * n);
        
        for (size_t i = 0; i < n; i++) {
            Expr* fi = head->data.function.args[i];
            
            Expr** call_args = malloc(sizeof(Expr*) * arg_count);
            for (size_t j = 0; j < arg_count; j++) {
                call_args[j] = expr_copy(args[j]);
            }
            Expr* call = expr_new_function(expr_copy(fi), call_args, arg_count);
            free(call_args);
            
            new_args[i] = evaluate(call);
            expr_free(call);
        }
        
        Expr* ret = expr_new_function(expr_copy(P), new_args, n);
        free(new_args);
        return ret;
    } else {
        size_t n = head->data.function.arg_count;
        Expr** new_args = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            new_args[i] = transform_head(head->data.function.args[i], h_spec, args, arg_count, transformed);
        }
        Expr* ret = expr_new_function(expr_copy(P), new_args, n);
        free(new_args);
        return ret;
    }
}

Expr* builtin_through(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    
    Expr* expr = res->data.function.args[0];
    Expr* h_spec = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    
    if (expr->type != EXPR_FUNCTION) return expr_copy(expr);
    
    Expr* head = expr->data.function.head;
    Expr** args = expr->data.function.args;
    size_t arg_count = expr->data.function.arg_count;
    
    bool transformed = false;
    Expr* new_expr = NULL;
    
    if (h_spec == NULL) {
        if (head->type == EXPR_FUNCTION) {
            Expr* P = head->data.function.head;
            size_t n = head->data.function.arg_count;
            Expr** new_H_args = malloc(sizeof(Expr*) * n);
            for (size_t i = 0; i < n; i++) {
                Expr* fi = head->data.function.args[i];
                
                Expr** call_args = malloc(sizeof(Expr*) * arg_count);
                for (size_t j = 0; j < arg_count; j++) {
                    call_args[j] = expr_copy(args[j]);
                }
                Expr* call = expr_new_function(expr_copy(fi), call_args, arg_count);
                free(call_args);
                
                new_H_args[i] = evaluate(call);
                expr_free(call);
            }
            new_expr = expr_new_function(expr_copy(P), new_H_args, n);
            free(new_H_args);
            transformed = true;
        }
    } else {
        new_expr = transform_head(head, h_spec, args, arg_count, &transformed);
    }
    
    if (transformed) {
        Expr* eval_res = evaluate(new_expr);
        expr_free(new_expr);
        return eval_res;
    }
    
    if (new_expr) expr_free(new_expr);
    return expr_copy(expr);
}

/* ------------------- FreeQ ------------------- */

static bool freeq_at_level(Expr* expr, Expr* form, int64_t current_level, LevelSpec spec) {
    int64_t d = get_depth(expr);
    bool should_check = false;
    
    if (spec.min == -1 && spec.max == -1) {
        if (expr->type != EXPR_FUNCTION || (expr->type == EXPR_FUNCTION && expr->data.function.head->type == EXPR_SYMBOL && (expr->data.function.head->data.symbol.name == SYM_Rational || expr->data.function.head->data.symbol.name == SYM_Complex))) {
            should_check = true;
        }
    } else {
        should_check = (current_level >= spec.min && current_level <= spec.max) ||
                       (-d >= spec.min && -d <= spec.max);
    }

    if (should_check) {
        MatchEnv* env = env_new();
        if (match(expr, form, env)) {
            env_free(env);
            return false; // Not free!
        }
        env_free(env);
    }

    if (expr->type == EXPR_FUNCTION && !(expr->data.function.head->type == EXPR_SYMBOL && (expr->data.function.head->data.symbol.name == SYM_Rational || expr->data.function.head->data.symbol.name == SYM_Complex))) {
        if (spec.heads) {
            if (!freeq_at_level(expr->data.function.head, form, current_level + 1, spec)) {
                return false;
            }
        }
        for (size_t i = 0; i < expr->data.function.arg_count; i++) {
            if (!freeq_at_level(expr->data.function.args[i], form, current_level + 1, spec)) {
                return false;
            }
        }
    }

    return true;
}

Expr* builtin_freeq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    /* A visible NDArray has no args[] for the level walk below to descend into,
     * so FreeQ[NDArray[Range[300]], 5] answered True where FreeQ[Range[300], 5]
     * answers False. Materialise and re-evaluate, exactly as the pattern family
     * in patterns.c does — see patterns_delist_visible there for why this is a
     * correctness fix and not a missing fast path. A packed argument has already
     * been materialised by the gate, so this fires only for a visible array.
     *
     * args[0] only, never the form being searched FOR: a visible NDArray[...]
     * and the plain List of the same values are distinct expressions, so
     * materialising the form would change what FreeQ is looking for. */
    if (is_ndarray(res->data.function.args[0])) {
        size_t n_ = res->data.function.arg_count;
        Expr** a_ = malloc(sizeof(Expr*) * n_);
        if (a_) {
            a_[0] = ndarray_to_nested_list(res->data.function.args[0]);
            for (size_t i = 1; i < n_; i++)
                a_[i] = expr_copy(res->data.function.args[i]);
            Expr* rb_ = expr_new_function(expr_copy(res->data.function.head), a_, n_);
            free(a_);
            return eval_and_free(rb_);
        }
    }
    Expr* expr = res->data.function.args[0];
    Expr* form = res->data.function.args[1];

    LevelSpec spec = {0, 1000000, true}; // Default: {0, Infinity}, Heads -> True
    if (res->data.function.arg_count >= 3) {
        Expr* ls = res->data.function.args[2];
        spec = parse_level_spec(ls, 1, 1);
        if (ls->type == EXPR_FUNCTION && ls->data.function.head->data.symbol.name == SYM_List) {
            // Already parsed by parse_level_spec
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol.name == SYM_All) {
            spec.min = 0; spec.max = 1000000;
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol.name == SYM_Infinity) {
            spec.min = 1; spec.max = 1000000;
        } else if (ls->type == EXPR_INTEGER) {
            spec.min = 1; spec.max = ls->data.integer;
        }
    }

    parse_options(res, res->data.function.arg_count >= 3 ? 3 : 2, &spec);
    
    // Check options anywhere for FreeQ
    for (size_t i = 2; i < res->data.function.arg_count; i++) {
        Expr* opt = res->data.function.args[i];
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->data.symbol.name == SYM_Rule) {
            if (opt->data.function.arg_count == 2 && 
                opt->data.function.args[0]->type == EXPR_SYMBOL &&
                opt->data.function.args[0]->data.symbol.name == SYM_Heads) {
                if (opt->data.function.args[1]->type == EXPR_SYMBOL &&
                    opt->data.function.args[1]->data.symbol.name == SYM_False) {
                    spec.heads = false;
                }
            }
        }
    }

    if (freeq_at_level(expr, form, 0, spec)) {
        return expr_new_symbol(SYM_True);
    } else {
        return expr_new_symbol(SYM_False);
    }
}

/* ------------------- Distribute ------------------- */

static void distribute_recursive(Expr*** components, size_t* component_counts, size_t n_args, size_t current_arg, Expr** current_tuple, Expr*** results, size_t* res_count, size_t* res_cap, Expr* fp_head) {
    if (current_arg == n_args) {
        Expr** tuple_copy = malloc(sizeof(Expr*) * n_args);
        for (size_t i = 0; i < n_args; i++) tuple_copy[i] = expr_copy(current_tuple[i]);
        Expr* term = expr_new_function(expr_copy(fp_head), tuple_copy, n_args);
        free(tuple_copy);
        
        if (*res_count >= *res_cap) {
            *res_cap = (*res_cap == 0) ? 16 : (*res_cap * 2);
            *results = realloc(*results, sizeof(Expr*) * (*res_cap));
        }
        (*results)[(*res_count)++] = term;
        return;
    }
    
    for (size_t i = 0; i < component_counts[current_arg]; i++) {
        current_tuple[current_arg] = components[current_arg][i];
        distribute_recursive(components, component_counts, n_args, current_arg + 1, current_tuple, results, res_count, res_cap, fp_head);
    }
}

Expr* builtin_distribute(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 5) {
        return NULL;
    }
    
    Expr* expr = res->data.function.args[0];
    if (expr->type != EXPR_FUNCTION) return expr_copy(expr);

    Expr* g = (res->data.function.arg_count >= 2) ? res->data.function.args[1] : expr_new_symbol(SYM_Plus);
    Expr* f = (res->data.function.arg_count >= 3) ? res->data.function.args[2] : expr_copy(expr->data.function.head);
    Expr* gp = (res->data.function.arg_count >= 4) ? res->data.function.args[3] : expr_copy(g);
    Expr* fp = (res->data.function.arg_count >= 5) ? res->data.function.args[4] : expr_copy(f);

    // If head(expr) != f, return expr
    if (!expr_eq(expr->data.function.head, f)) {
        if (res->data.function.arg_count < 2) expr_free(g);
        if (res->data.function.arg_count < 3) expr_free(f);
        if (res->data.function.arg_count < 4) expr_free(gp);
        if (res->data.function.arg_count < 5) expr_free(fp);
        return expr_copy(expr);
    }

    size_t n_args = expr->data.function.arg_count;
    Expr*** components = malloc(sizeof(Expr**) * n_args);
    size_t* component_counts = malloc(sizeof(size_t) * n_args);
    bool any_g = false;

    for (size_t i = 0; i < n_args; i++) {
        Expr* arg = expr->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && expr_eq(arg->data.function.head, g)) {
            components[i] = arg->data.function.args;
            component_counts[i] = arg->data.function.arg_count;
            any_g = true;
        } else {
            components[i] = &expr->data.function.args[i];
            component_counts[i] = 1;
        }
    }

    if (!any_g) {
        // No distribution possible, but we still apply fp if it's different?
        // Actually Distribute[f[x], g, f, gp, fp] -> gp[fp[x]] if no g is found?
        // Let's re-read: "Distribute explicitly constructs the complete result of a distribution"
        // If no g is found, it usually returns fp[args] wrapped in gp? No, wait.
        // Distribute[f[a], g] -> f[a]
        // Distribute[f[g[a]], g] -> g[f[a]]
        
        // Cleanup and return original if no g found
        free(components);
        free(component_counts);
        if (res->data.function.arg_count < 2) expr_free(g);
        if (res->data.function.arg_count < 3) expr_free(f);
        if (res->data.function.arg_count < 4) expr_free(gp);
        if (res->data.function.arg_count < 5) expr_free(fp);
        return expr_copy(expr);
    }

    Expr** results = NULL;
    size_t res_count = 0;
    size_t res_cap = 0;
    Expr** current_tuple = malloc(sizeof(Expr*) * n_args);

    distribute_recursive(components, component_counts, n_args, 0, current_tuple, &results, &res_count, &res_cap, fp);

    Expr* final_res = expr_new_function(expr_copy(gp), results, res_count);
    if (results) free(results);
    
    // results array itself is consumed by expr_new_function, but only if we don't free it.
    // wait, expr_new_function usually COPIES the array if we pass it, OR it takes ownership?
    // In Mathilda, expr_new_function takes Expr** args.
    // Let's check expr_new_function implementation.
    
    free(current_tuple);
    free(components);
    free(component_counts);
    if (res->data.function.arg_count < 2) expr_free(g);
    if (res->data.function.arg_count < 3) expr_free(f);
    if (res->data.function.arg_count < 4) expr_free(gp);
    if (res->data.function.arg_count < 5) expr_free(fp);
    
    // We need to evaluate the result because gp or fp might be builtins
    Expr* eval_res = evaluate(final_res);
    expr_free(final_res);
    return eval_res;
}

/* ------------------- Inner ------------------- */

/* --- Inner with n=1: contract first index of A with first index of B --- */

static Expr* inner_n1_B(Expr* f, Expr** A_leaves, Expr** B_slices, size_t L, Expr* g, Expr* head) {
    Expr* first_b = B_slices[0];

    if (first_b->type == EXPR_FUNCTION && expr_eq(first_b->data.function.head, head)) {
        /* B has remaining structure — iterate over its sub-elements ("columns") */
        size_t M = first_b->data.function.arg_count;
        Expr** res_args = malloc(sizeof(Expr*) * M);
        Expr** new_B = malloc(sizeof(Expr*) * L);
        for (size_t j = 0; j < M; j++) {
            for (size_t k = 0; k < L; k++) {
                Expr* bk = B_slices[k];
                if (bk->type == EXPR_FUNCTION && j < bk->data.function.arg_count) {
                    new_B[k] = bk->data.function.args[j];
                } else {
                    new_B[k] = bk; /* fallback for ragged */
                }
            }
            res_args[j] = inner_n1_B(f, A_leaves, new_B, L, g, head);
        }
        free(new_B);
        Expr* ret = expr_new_function(expr_copy(head), res_args, M);
        free(res_args);
        return ret;
    } else {
        /* Base case: both A and B at leaf level — form g[f[a0,b0], f[a1,b1], ...] */
        Expr** g_args = malloc(sizeof(Expr*) * L);
        for (size_t k = 0; k < L; k++) {
            Expr* f_args_arr[2];
            f_args_arr[0] = expr_copy(A_leaves[k]);
            f_args_arr[1] = expr_copy(B_slices[k]);
            g_args[k] = expr_new_function(expr_copy(f), f_args_arr, 2);
        }
        Expr* ret = expr_new_function(expr_copy(g), g_args, L);
        free(g_args);
        return ret;
    }
}

static Expr* inner_n1_A(Expr* f, Expr** A_slices, Expr** B_slices, size_t L, Expr* g, Expr* head) {
    Expr* first_a = A_slices[0];

    /* If first A slice is not a list, A is at leaf level — go to B-side */
    if (first_a->type != EXPR_FUNCTION || !expr_eq(first_a->data.function.head, head)) {
        return inner_n1_B(f, A_slices, B_slices, L, g, head);
    }

    /* A still has remaining structure — descend one level along A's next free index */
    size_t M = first_a->data.function.arg_count;
    Expr** res_args = malloc(sizeof(Expr*) * M);
    Expr** new_slices = malloc(sizeof(Expr*) * L);
    for (size_t i = 0; i < M; i++) {
        for (size_t k = 0; k < L; k++) {
            Expr* ak = A_slices[k];
            if (ak->type == EXPR_FUNCTION && i < ak->data.function.arg_count) {
                new_slices[k] = ak->data.function.args[i];
            } else {
                new_slices[k] = ak; /* fallback for ragged */
            }
        }
        res_args[i] = inner_n1_A(f, new_slices, B_slices, L, g, head);
    }
    free(new_slices);
    Expr* ret = expr_new_function(expr_copy(head), res_args, M);
    free(res_args);
    return ret;
}

/* --- Standard Inner: contract last index of A with first index of B --- */

static Expr* contract_V_B(Expr* f, Expr* V, Expr* B, Expr* g, Expr* head) {
    if (V->type != EXPR_FUNCTION || B->type != EXPR_FUNCTION) return NULL;
    size_t N = V->data.function.arg_count;
    
    bool b_is_matrix = false;
    if (B->data.function.arg_count > 0 && B->data.function.args[0]->type == EXPR_FUNCTION &&
        expr_eq(B->data.function.args[0]->data.function.head, head)) {
        b_is_matrix = true;
    }
    
    if (b_is_matrix) {
        size_t M = B->data.function.args[0]->data.function.arg_count;
        Expr** res_args = malloc(sizeof(Expr*) * M);
        for (size_t j = 0; j < M; j++) {
            Expr** col_args = malloc(sizeof(Expr*) * N);
            for (size_t i = 0; i < N; i++) {
                if (i < B->data.function.arg_count) {
                    Expr* B_i = B->data.function.args[i];
                    if (B_i->type == EXPR_FUNCTION && j < B_i->data.function.arg_count) {
                        col_args[i] = expr_copy(B_i->data.function.args[j]);
                    } else {
                        col_args[i] = expr_new_symbol(SYM_Null);
                    }
                } else {
                    col_args[i] = expr_new_symbol(SYM_Null);
                }
            }
            Expr* B_col = expr_new_function(expr_copy(head), col_args, N);
            Expr* contracted = contract_V_B(f, V, B_col, g, head);
            expr_free(B_col);
            if (!contracted) {
                for (size_t k = 0; k < j; k++) expr_free(res_args[k]);
                free(res_args);
                return NULL;
            }
            res_args[j] = contracted;
        }
        Expr* ret = expr_new_function(expr_copy(head), res_args, M);
        free(res_args);
        return ret;
    } else {
        size_t B_len = B->data.function.arg_count;
        size_t min_len = N < B_len ? N : B_len;
        if (N != B_len) {
            // Technically lengths should match. We'll evaluate up to the min_len.
        }
        Expr** g_args = malloc(sizeof(Expr*) * min_len);
        for (size_t i = 0; i < min_len; i++) {
            Expr* f_args[2] = { expr_copy(V->data.function.args[i]), expr_copy(B->data.function.args[i]) };
            g_args[i] = expr_new_function(expr_copy(f), f_args, 2);
        }
        Expr* ret = expr_new_function(expr_copy(g), g_args, min_len);
        free(g_args);
        return ret;
    }
}

static Expr* inner_A(Expr* f, Expr* A, Expr* B, Expr* g, Expr* head) {
    if (A->type != EXPR_FUNCTION) return NULL;
    bool a_is_matrix = false;
    if (A->data.function.arg_count > 0 && A->data.function.args[0]->type == EXPR_FUNCTION &&
        expr_eq(A->data.function.args[0]->data.function.head, head)) {
        a_is_matrix = true;
    }
    
    if (a_is_matrix) {
        size_t N = A->data.function.arg_count;
        Expr** res_args = malloc(sizeof(Expr*) * N);
        for (size_t i = 0; i < N; i++) {
            Expr* inner_res = inner_A(f, A->data.function.args[i], B, g, head);
            if (!inner_res) {
                for (size_t j = 0; j < i; j++) expr_free(res_args[j]);
                free(res_args);
                return NULL;
            }
            res_args[i] = inner_res;
        }
        Expr* ret = expr_new_function(expr_copy(head), res_args, N);
        free(res_args);
        return ret;
    } else {
        return contract_V_B(f, A, B, g, head);
    }
}

Expr* builtin_inner(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 3) return NULL;
    
    Expr* f = res->data.function.args[0];
    Expr* A = res->data.function.args[1];
    Expr* B = res->data.function.args[2];
    Expr* g = (res->data.function.arg_count >= 4) ? res->data.function.args[3] : expr_new_symbol(SYM_Plus);
    
    Expr* n_expr = (res->data.function.arg_count >= 5) ? res->data.function.args[4] : NULL;

    /* Inner[Times, a, b, Plus] IS Dot -- that is the definition of a Dot, and
     * the pair is how a contraction gets written when the reader wants the
     * operators named. With two buffers that is ndarray_dot2, i.e. BLAS, where
     * the generic path below builds one Times Expr per element pair and folds
     * them through the evaluator. Only this operator pair qualifies: Inner's
     * whole point is that f and g may be anything, and anything else has no
     * kernel. */
    if (is_ndarray(A) && is_ndarray(B) && !n_expr &&
        f->type == EXPR_SYMBOL && f->data.symbol.name == SYM_Times &&
        g->type == EXPR_SYMBOL && g->data.symbol.name == SYM_Plus) {
        bool shape_error = false;
        Expr* dot = ndarray_dot2(A, B, &shape_error);
        if (dot) {
            if (res->data.function.arg_count < 4) expr_free(g);
            return dot;
        }
    }

    /* Any other operator pair, or a shape ndarray_dot2 declines: materialise
     * and let the generic contraction below answer. Inner is on pack.c's AWARE
     * list for the Dot case above, so without this an f other than Times would
     * reach the EXPR_FUNCTION test with a buffer and fall out UNEVALUATED --
     * the exact regression that marking a head aware invites. */
    if (is_ndarray(A) || is_ndarray(B)) {
        if (res->data.function.arg_count < 4) expr_free(g);
        return ndarray_delist_and_reeval(res);
    }

    if (A->type != EXPR_FUNCTION || B->type != EXPR_FUNCTION) {
        if (res->data.function.arg_count < 4) expr_free(g);
        return NULL;
    }
    
    Expr* head = A->data.function.head;

    Expr* inner_res;
    if (n_expr && n_expr->type == EXPR_INTEGER && n_expr->data.integer == 1) {
        /* Contract first index of A with first index of B directly */
        size_t L = A->data.function.arg_count;
        if (L != B->data.function.arg_count) {
            if (res->data.function.arg_count < 4) expr_free(g);
            return NULL;
        }
        inner_res = inner_n1_A(f, A->data.function.args, B->data.function.args, L, g, head);
    } else {
        Expr* A_used = expr_copy(A);
        inner_res = inner_A(f, A_used, B, g, head);
        expr_free(A_used);
    }

    if (res->data.function.arg_count < 4) expr_free(g);

    if (!inner_res) return NULL;

    Expr* final_eval = evaluate(inner_res);
    expr_free(inner_res);
    return final_eval;
}

/* ------------------- Outer ------------------- */

static Expr* outer_rec(Expr* f, Expr** orig_tensors, size_t num_tensors, int64_t* target_depths, 
                size_t arg_idx, Expr* curr_subtensor, int64_t curr_depth, Expr** current_atoms, Expr* head) {
    
    bool treat_as_atom = false;
    if (curr_depth >= target_depths[arg_idx]) {
        treat_as_atom = true;
    } else if (curr_subtensor->type != EXPR_FUNCTION || (head && !expr_eq(curr_subtensor->data.function.head, head))) {
        treat_as_atom = true;
    }
    
    if (treat_as_atom) {
        current_atoms[arg_idx] = curr_subtensor;
if (arg_idx + 1 == num_tensors) {
            Expr** f_args = malloc(sizeof(Expr*) * num_tensors);
            for(size_t i=0; i<num_tensors; i++) f_args[i] = expr_copy(current_atoms[i]);
            Expr* ret = expr_new_function(expr_copy(f), f_args, num_tensors);
            free(f_args);
            return ret;
        } else {
            return outer_rec(f, orig_tensors, num_tensors, target_depths, 
                             arg_idx + 1, orig_tensors[arg_idx + 1], 0, current_atoms, head);
        }
    }
    
    size_t count = curr_subtensor->data.function.arg_count;
    Expr** new_args = malloc(sizeof(Expr*) * count);
    for (size_t i = 0; i < count; i++) {
        new_args[i] = outer_rec(f, orig_tensors, num_tensors, target_depths,
                                arg_idx, curr_subtensor->data.function.args[i], curr_depth + 1, current_atoms, head);
    }
    Expr* ret = expr_new_function(expr_copy(head), new_args, count);
    free(new_args);
    return ret;
}

/* Outer[f, a, b] over two rank-1 float64 buffers, for the arithmetic heads.
 *
 * The generic path builds the ENTIRE result as a symbolic tree -- one
 * f[a_i, b_j] node per output element -- and then calls evaluate() on it, so a
 * 1000x1000 outer product allocates and evaluates a million nodes: 753 ms, or
 * ~750 ns per element, against ~3 ms for the same product written into a
 * buffer. That is the shape an all-pairs kernel is written in (displacement
 * matrices in an N-body step, a Gram matrix, a distance matrix), and it was the
 * single slowest thing in the third benchmark sweep.
 *
 * Only float64 in, and only total functions of two finite doubles. Integer
 * operands keep the generic path, where exactness and bignum promotion are
 * already handled; a non-finite result anywhere abandons the buffer and lets
 * the generic path produce whatever symbolic form it produces (DirectedInfinity
 * and friends), so this can only ever be faster, never different. */
/* A plain numeric List lifted to a float64 buffer so it can pair with one.
 * Threshold-ignoring (pack_force): the decision here is about the PAIR, and the
 * pair is already large. NULL when the list is not uniformly machine-numeric. */
static Expr* outer_pack_operand(Expr* list) {
    if (!list || list->type != EXPR_FUNCTION) return NULL;
    Expr* cand = pack_force(expr_copy(list), false, NDT_FLOAT64);
    if (cand && cand->type == EXPR_NDARRAY) return cand;
    if (cand) expr_free(cand);
    return NULL;
}

static Expr* nd_outer2(const Expr* f, const Expr* a, const Expr* b) {
    if (!f || f->type != EXPR_SYMBOL) return NULL;
    if (a->type != EXPR_NDARRAY || b->type != EXPR_NDARRAY) return NULL;
    if (a->data.ndarray.rank != 1 || b->data.ndarray.rank != 1) return NULL;
    if (a->data.ndarray.dtype != NDT_FLOAT64 ||
        b->data.ndarray.dtype != NDT_FLOAT64) return NULL;

    const char* op = f->data.symbol.name;
    enum { OP_PLUS, OP_SUB, OP_TIMES, OP_MIN, OP_MAX } k;
    if      (op == SYM_Plus)     k = OP_PLUS;
    else if (op == SYM_Subtract) k = OP_SUB;
    else if (op == SYM_Times)    k = OP_TIMES;
    else if (op == SYM_Min)      k = OP_MIN;
    else if (op == SYM_Max)      k = OP_MAX;
    else return NULL;

    int64_t m = a->data.ndarray.dims[0], n = b->data.ndarray.dims[0];
    const double* A = (const double*)a->data.ndarray.data;
    const double* B = (const double*)b->data.ndarray.data;
    size_t total = (size_t)m * (size_t)n;
    double* out = malloc(sizeof(double) * (total ? total : 1));
    if (!out) return NULL;

    bool bad = false;
    for (int64_t i = 0; i < m; i++) {
        double ai = A[i];
        double* row = out + (size_t)i * (size_t)n;
        switch (k) {
        case OP_PLUS:  for (int64_t j = 0; j < n; j++) row[j] = ai + B[j]; break;
        case OP_SUB:   for (int64_t j = 0; j < n; j++) row[j] = ai - B[j]; break;
        case OP_TIMES: for (int64_t j = 0; j < n; j++) row[j] = ai * B[j]; break;
        case OP_MIN:   for (int64_t j = 0; j < n; j++) row[j] = ai < B[j] ? ai : B[j]; break;
        case OP_MAX:   for (int64_t j = 0; j < n; j++) row[j] = ai > B[j] ? ai : B[j]; break;
        }
        for (int64_t j = 0; j < n; j++) if (!isfinite(row[j])) { bad = true; break; }
        if (bad) break;
    }
    if (bad) { free(out); return NULL; }

    int64_t dims[2] = { m, n };
    return expr_new_ndarray_like(nd_present_src2(a, b), 2, dims, out, NDT_FLOAT64);
}

Expr* builtin_outer(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;

    Expr* f = res->data.function.args[0];

    if (res->data.function.arg_count == 3) {
        /* One operand packed and the other an ordinary List: pack the plain one
         * so the PAIR can take the buffer path — the same rule dot2 keeps, for
         * the same reason. Outer's cost is the PRODUCT of its operands' lengths,
         * but PACK_MIN_ELEMENTS judges each one alone, so a k-means with 16
         * centroids had its 16-element operand correctly left unpacked and paid
         * for it across 1.6 million output cells: 849 ms, against 4 ms once both
         * sides are buffers. Packing is value-preserving, so lifting the small
         * side costs nothing and there is nothing to weigh. */
        Expr *ta = NULL, *tb = NULL;
        Expr *a = res->data.function.args[1], *b = res->data.function.args[2];
        if (pack_any_created()) {
            if (a->type == EXPR_NDARRAY && b->type != EXPR_NDARRAY) {
                tb = outer_pack_operand(b);
                if (tb) b = tb;
            } else if (b->type == EXPR_NDARRAY && a->type != EXPR_NDARRAY) {
                ta = outer_pack_operand(a);
                if (ta) a = ta;
            }
        }
        Expr* fast = nd_outer2(f, a, b);
        expr_free(ta);
        expr_free(tb);
        if (fast) return fast;
    }
    /* Outer is packed-aware ONLY so nd_outer2 can see a buffer. Everything below
     * this point walks the operands as nested Lists (outer_rec indexes
     * data.function.args), which an NDArray does not have, so anything the fast
     * path declined has to be handed back materialised. */
    for (size_t i = 1; i < res->data.function.arg_count; i++) {
        if (is_ndarray(res->data.function.args[i]))
            return ndarray_delist_and_reeval(res);
    }
    
    size_t num_depths = 0;
    for (int64_t i = res->data.function.arg_count - 1; i >= 1; i--) {
        Expr* a = res->data.function.args[i];
        if (a->type == EXPR_INTEGER || (a->type == EXPR_SYMBOL && a->data.symbol.name == SYM_Infinity)) {
            num_depths++;
        } else {
            break;
        }
    }
    
    size_t num_tensors = res->data.function.arg_count - 1 - num_depths;
    if (num_tensors == 0) {
        Expr* ret = expr_new_function(expr_copy(f), NULL, 0);
        Expr* evaluated = evaluate(ret);
        expr_free(ret);
        return evaluated;
    }
    
    Expr** tensors = &res->data.function.args[1];
    
    int64_t* target_depths = malloc(sizeof(int64_t) * num_tensors);
    for (size_t i = 0; i < num_tensors; i++) target_depths[i] = INT64_MAX;

    if (num_depths == 1) {
        Expr* d = res->data.function.args[1 + num_tensors];
        int64_t val = (d->type == EXPR_INTEGER) ? d->data.integer : INT64_MAX;
        for (size_t i = 0; i < num_tensors; i++) target_depths[i] = val;
    } else if (num_depths > 1) {
        for (size_t i = 0; i < num_depths && i < num_tensors; i++) {
            Expr* d = res->data.function.args[1 + num_tensors + i];
            target_depths[i] = (d->type == EXPR_INTEGER) ? d->data.integer : INT64_MAX;
        }
    }
    
    Expr* head = NULL;
    for (size_t i = 0; i < num_tensors; i++) {
        if (tensors[i]->type == EXPR_FUNCTION) {
            head = tensors[i]->data.function.head;
            break;
        }
    }
    
    Expr** current_atoms = malloc(sizeof(Expr*) * num_tensors);
    Expr* raw_res = outer_rec(f, tensors, num_tensors, target_depths, 0, tensors[0], 0, current_atoms, head);
    
    free(current_atoms);
    free(target_depths);
    
    Expr* final_eval = evaluate(raw_res);
    expr_free(raw_res);
    return final_eval;
}

/* ------------------- Tuples ------------------- */

static void tuples_rec(Expr** lists, size_t num_lists, size_t curr_list, Expr** current_tuple, Expr* head, Expr*** results, size_t* res_count, size_t* res_cap) {
    if (curr_list == num_lists) {
        Expr** t_args = malloc(sizeof(Expr*) * num_lists);
        for(size_t i=0; i<num_lists; i++) t_args[i] = expr_copy(current_tuple[i]);
        Expr* t = expr_new_function(expr_copy(head), t_args, num_lists);
        free(t_args);  /* expr_new_function copies the array; free the container */
        if (*res_count >= *res_cap) {
            *res_cap = (*res_cap == 0) ? 16 : (*res_cap * 2);
            *results = realloc(*results, sizeof(Expr*) * (*res_cap));
        }
        (*results)[(*res_count)++] = t;
        return;
    }
    
    Expr* lst = lists[curr_list];
    size_t count = lst->type == EXPR_FUNCTION ? lst->data.function.arg_count : 0;
    if (count == 0) return;
    
    for (size_t i = 0; i < count; i++) {
        current_tuple[curr_list] = lst->data.function.args[i];
        tuples_rec(lists, num_lists, curr_list + 1, current_tuple, head, results, res_count, res_cap);
    }
}

static Expr* reshape_rec(Expr** flat_args, size_t* offset, int64_t* dims, size_t dims_count, size_t current_dim, Expr* head) {
    if (current_dim == dims_count - 1) {
        size_t n = dims[current_dim];
        Expr** args = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            args[i] = expr_copy(flat_args[(*offset)++]);
        }
        Expr* ret = expr_new_function(expr_copy(head), args, n);
        free(args);
        return ret;
    } else {
        size_t n = dims[current_dim];
        Expr** args = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            args[i] = reshape_rec(flat_args, offset, dims, dims_count, current_dim + 1, head);
        }
        Expr* ret = expr_new_function(expr_copy(head), args, n);
        free(args);
        return ret;
    }
}

Expr* builtin_tuples(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    
    if (res->data.function.arg_count == 1) {
        Expr* arg = res->data.function.args[0];
        if (arg->type != EXPR_FUNCTION) return expr_copy(res);
        Expr* head = arg->data.function.head;
        
        size_t num_lists = arg->data.function.arg_count;
        Expr** lists = arg->data.function.args;
        
        if (num_lists == 0) {
            Expr* t = expr_new_function(expr_copy(head), NULL, 0);
            Expr* final_res = expr_new_function(expr_copy(head), (Expr*[]){t}, 1);
            return final_res;
        }

        Expr** results = NULL;
        size_t res_count = 0, res_cap = 0;
        Expr** current_tuple = malloc(sizeof(Expr*) * num_lists);
        
        tuples_rec(lists, num_lists, 0, current_tuple, head, &results, &res_count, &res_cap);
        free(current_tuple);

        Expr* final_res = expr_new_function(expr_copy(head), results, res_count);
        free(results);  /* expr_new_function copied the pointers into its own storage */
        return final_res;
    } else {
        Expr* list = res->data.function.args[0];
        Expr* n_expr = res->data.function.args[1];
        if (list->type != EXPR_FUNCTION) return expr_copy(res);
        Expr* head = list->data.function.head;
        
        if (n_expr->type == EXPR_INTEGER) {
            int64_t n = n_expr->data.integer;
            if (n < 0) n = 0;
            if (n == 0) {
                Expr* t = expr_new_function(expr_copy(head), NULL, 0);
                Expr* final_res = expr_new_function(expr_copy(head), (Expr*[]){t}, 1);
                return final_res;
            }

            Expr** lists = malloc(sizeof(Expr*) * n);
            for (int64_t i = 0; i < n; i++) lists[i] = list;
            
            Expr** results = NULL;
            size_t res_count = 0, res_cap = 0;
            Expr** current_tuple = malloc(sizeof(Expr*) * n);
            
            tuples_rec(lists, n, 0, current_tuple, head, &results, &res_count, &res_cap);
            free(current_tuple);
            free(lists);

            Expr* out = expr_new_function(expr_new_symbol(SYM_List), results, res_count);
            free(results);
            return out;
        } else if (n_expr->type == EXPR_FUNCTION && expr_eq(n_expr->data.function.head, head)) {
            size_t dims_count = n_expr->data.function.arg_count;
            int64_t total_elements = 1;
            int64_t* dims = malloc(sizeof(int64_t) * dims_count);
            for (size_t i = 0; i < dims_count; i++) {
                if (n_expr->data.function.args[i]->type != EXPR_INTEGER) {
                    free(dims);
                    return expr_copy(res);
                }
                dims[i] = n_expr->data.function.args[i]->data.integer;
                if (dims[i] < 0) dims[i] = 0;
                total_elements *= dims[i];
            }
            
            if (total_elements == 0) {
                free(dims);
                return expr_new_function(expr_copy(head), NULL, 0);
            }

            Expr** lists = malloc(sizeof(Expr*) * total_elements);
            for (int64_t i = 0; i < total_elements; i++) lists[i] = list;
            
            Expr** results = NULL;
            size_t res_count = 0, res_cap = 0;
            Expr** current_tuple = malloc(sizeof(Expr*) * total_elements);
            
            tuples_rec(lists, total_elements, 0, current_tuple, head, &results, &res_count, &res_cap);
            free(current_tuple);
            free(lists);
            
            for (size_t i = 0; i < res_count; i++) {
                size_t offset = 0;
                Expr* reshaped = reshape_rec(results[i]->data.function.args, &offset, dims, dims_count, 0, head);
                expr_free(results[i]);
                results[i] = reshaped;
            }
            
            free(dims);
            Expr* out = expr_new_function(expr_new_symbol(SYM_List), results, res_count);
            free(results);
            return out;
        }
    }
    return expr_copy(res);
}

/* ------------------- Permutations ------------------- */

typedef struct {
    Expr* expr;
    int count;
} UniqueElement;

static void permutations_rec(UniqueElement* elements, size_t num_unique, size_t target_len, size_t current_len, 
                             Expr** current_perm, Expr* head, Expr*** results, size_t* res_count, size_t* res_cap) {
    if (current_len == target_len) {
        Expr** p_args = NULL;
        if (target_len > 0) {
            p_args = malloc(sizeof(Expr*) * target_len);
            for (size_t i = 0; i < target_len; i++) p_args[i] = expr_copy(current_perm[i]);
        }
        Expr* p = expr_new_function(expr_copy(head), p_args, target_len);
        if (p_args) free(p_args);
        if (*res_count >= *res_cap) {
            *res_cap = (*res_cap == 0) ? 16 : (*res_cap * 2);
            *results = realloc(*results, sizeof(Expr*) * (*res_cap));
        }
        (*results)[(*res_count)++] = p;
        return;
    }
    
    for (size_t i = 0; i < num_unique; i++) {
        if (elements[i].count > 0) {
            elements[i].count--;
            current_perm[current_len] = elements[i].expr;
            permutations_rec(elements, num_unique, target_len, current_len + 1, current_perm, head, results, res_count, res_cap);
            elements[i].count++;
        }
    }
}

Expr* builtin_permutations(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    
    Expr* list = res->data.function.args[0];
    if (list->type != EXPR_FUNCTION) return expr_copy(res); 
    
    Expr* head = list->data.function.head;
    size_t len = list->data.function.arg_count;
    
    int64_t n_min = len;
    int64_t n_max = len;
    int64_t dn = 1;
    
    if (res->data.function.arg_count == 2) {
        Expr* spec = res->data.function.args[1];
        if (spec->type == EXPR_INTEGER) {
            int64_t n = spec->data.integer;
            if (n > (int64_t)len) n = len;
            if (n < 0) n = 0;
            n_min = n;
            n_max = n;
            dn = 1;
        } else if (spec->type == EXPR_SYMBOL && spec->data.symbol.name == SYM_All) {
            n_min = 0;
            n_max = len;
            dn = 1;
        } else if (spec->type == EXPR_FUNCTION && spec->data.function.head->data.symbol.name == SYM_List) {
            size_t s_len = spec->data.function.arg_count;
            if (s_len == 1 && spec->data.function.args[0]->type == EXPR_INTEGER) {
                n_min = n_max = spec->data.function.args[0]->data.integer;
            } else if (s_len == 2 && spec->data.function.args[0]->type == EXPR_INTEGER && spec->data.function.args[1]->type == EXPR_INTEGER) {
                n_min = spec->data.function.args[0]->data.integer;
                n_max = spec->data.function.args[1]->data.integer;
                dn = (n_min <= n_max) ? 1 : -1;
            } else if (s_len == 3 && spec->data.function.args[0]->type == EXPR_INTEGER && spec->data.function.args[1]->type == EXPR_INTEGER && spec->data.function.args[2]->type == EXPR_INTEGER) {
                n_min = spec->data.function.args[0]->data.integer;
                n_max = spec->data.function.args[1]->data.integer;
                dn = spec->data.function.args[2]->data.integer;
            }
        }
    }
    
    UniqueElement* elements = malloc(sizeof(UniqueElement) * len);
    size_t num_unique = 0;
    for (size_t i = 0; i < len; i++) {
        Expr* e = list->data.function.args[i];
        bool found = false;
        for (size_t j = 0; j < num_unique; j++) {
            if (expr_eq(e, elements[j].expr)) {
                elements[j].count++;
                found = true;
                break;
            }
        }
        if (!found) {
            elements[num_unique].expr = e;
            elements[num_unique].count = 1;
            num_unique++;
        }
    }
    
    Expr** results = NULL;
    size_t res_count = 0, res_cap = 0;
    Expr** current_perm = malloc(sizeof(Expr*) * (n_min > n_max ? n_min : n_max) + sizeof(Expr*) * len); 
    
    if (dn != 0) {
        if (dn > 0) {
            for (int64_t l = n_min; l <= n_max; l += dn) {
                if (l >= 0 && l <= (int64_t)len) {
                    permutations_rec(elements, num_unique, l, 0, current_perm, head, &results, &res_count, &res_cap);
                }
            }
        } else {
            for (int64_t l = n_min; l >= n_max; l += dn) {
                if (l >= 0 && l <= (int64_t)len) {
                    permutations_rec(elements, num_unique, l, 0, current_perm, head, &results, &res_count, &res_cap);
                }
            }
        }
    }
    
    free(current_perm);
    free(elements);
    
    Expr* final_res = expr_new_function(expr_new_symbol(SYM_List), results, res_count);
    if (results) free(results);
    return final_res;
}

/* ================================================================= *
 * Shared iteration machinery for Nest/NestList, Fold/FoldList,
 * NestWhile/NestWhileList, and FixedPoint/FixedPointList.
 *
 * Each of those builtin pairs is a thin wrapper around a single
 * implementation function that takes a boolean "as_list" flag. The
 * implementations differ only in argument parsing, the per-iteration
 * step callback, and any post-processing of the history buffer. The
 * generic runner iter_run() and the finalization helper ebuf_finalize()
 * handle the common concerns: history growth, safety capping, early
 * return, and packaging the result as either a full list of iterates
 * (the *List variant) or just the last one (the scalar variant).
 * ================================================================= */

/* --- Dynamic buffer of Expr* used as iteration history. --- */

typedef struct {
    Expr** items;
    size_t count;
    size_t capacity;
} ExprBuf;

static void ebuf_init(ExprBuf* b) {
    b->capacity = 16;
    b->count = 0;
    b->items = malloc(sizeof(Expr*) * b->capacity);
}

static void ebuf_push(ExprBuf* b, Expr* e) {
    if (b->count >= b->capacity) {
        b->capacity *= 2;
        b->items = realloc(b->items, sizeof(Expr*) * b->capacity);
    }
    b->items[b->count++] = e;
}

/* Free every stored expression and the backing array, leaving buf empty. */
static void ebuf_free_all(ExprBuf* b) {
    for (size_t i = 0; i < b->count; i++) expr_free(b->items[i]);
    free(b->items);
    b->items = NULL;
    b->count = b->capacity = 0;
}

/* Drop the last `drop` items (freeing them). Safe if drop >= count. */
static void ebuf_truncate(ExprBuf* b, size_t drop) {
    if (drop >= b->count) drop = b->count;
    for (size_t i = b->count - drop; i < b->count; i++) expr_free(b->items[i]);
    b->count -= drop;
}

/* Drop the first `drop` items (freeing them) and shift the rest down. Used by
 * the bounded-history window (see iter_run's max_history): scalar Nest/Fold/
 * FixedPoint only ever read the most recent iterate, so superseded ones are
 * freed mid-loop instead of piling up to O(n) live nodes (which would defeat
 * the Expr node pool). Safe if drop >= count. */
static void ebuf_drop_front(ExprBuf* b, size_t drop) {
    if (drop >= b->count) drop = b->count;
    for (size_t i = 0; i < drop; i++) expr_free(b->items[i]);
    if (drop < b->count)
        memmove(b->items, b->items + drop, (b->count - drop) * sizeof(Expr*));
    b->count -= drop;
}

/*
 * Convert an iteration history into a result expression.
 *   as_list=true  -> returns out_head[items...], taking ownership of out_head
 *                    and of every Expr* in buf.
 *   as_list=false -> returns the last element of buf, freeing the others and
 *                    freeing out_head (which is unused in this form).
 * In either case the backing array is freed. If buf is empty in the scalar
 * form, out_head is freed and NULL is returned.
 */
static Expr* ebuf_finalize(ExprBuf* b, bool as_list, Expr* out_head) {
    Expr* result;
    if (as_list) {
        result = expr_new_function(out_head, b->items, b->count);
        free(b->items);
        /* The iterate history of NestList / FoldList / NestWhileList /
         * FixedPointList. Offered here rather than at each of the four call
         * sites, and it also keeps this path in step with numloop's compiled one
         * (reals_to_list packs), so whether the body happened to compile does not
         * decide the result's representation. */
        result = pack_offer(result);
    } else {
        expr_free(out_head);
        if (b->count == 0) {
            free(b->items);
            return NULL;
        }
        result = b->items[b->count - 1];
        for (size_t i = 0; i + 1 < b->count; i++) expr_free(b->items[i]);
        free(b->items);
    }
    b->items = NULL;
    b->count = b->capacity = 0;
    return result;
}

/* --- Small helpers for building and evaluating calls. --- */

static inline Expr* apply_unary(Expr* f, Expr* arg) {
    Expr* a = expr_copy(arg);
    Expr* call = expr_new_function(expr_copy(f), &a, 1);
    return eval_and_free(call);
}

static inline Expr* apply_binary(Expr* f, Expr* a, Expr* b) {
    Expr* args[2] = { expr_copy(a), expr_copy(b) };
    Expr* call = expr_new_function(expr_copy(f), args, 2);
    return eval_and_free(call);
}

static inline bool sym_is_true(Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_True;
}

/* --- Generic iteration runner. --- */

typedef enum {
    ITER_STEP_CONT,       /* push *out and keep going */
    ITER_STEP_HALT,       /* stop; *out is ignored */
    ITER_STEP_HALT_ADD,   /* push *out and stop */
    ITER_STEP_RETURN      /* free history and hand *out to the caller */
} IterStep;

typedef IterStep (*IterStepFn)(ExprBuf* hist, void* ctx, Expr** out_next);

typedef enum {
    ITER_RUN_OK,       /* normal completion; buf contains the history */
    ITER_RUN_EARLY,    /* early return from step; buf freed, value in *out_early */
    ITER_RUN_SAFETY    /* hit safety cap; buf freed */
} IterRunResult;

/* Hard cap for "infinite" runs so a runaway iteration can't consume all memory. */
#define ITER_SAFETY_CAP ((int64_t)1000000)

/*
 * Iterate until step_fn halts or `limit` applications have been performed.
 *
 * limit_is_safety=false  -> `limit` is a user-facing bound (the iteration
 *                           simply stops on reaching it; whatever is in buf
 *                           is the result).
 * limit_is_safety=true   -> `limit` is a safety cap for an otherwise unbounded
 *                           run; reaching it is treated as failure (buf is
 *                           freed and ITER_RUN_SAFETY is returned).
 *
 * max_history bounds how many trailing iterates the history buffer retains.
 * The *List variants keep every iterate (pass SIZE_MAX). The scalar variants
 * that only ever read the most recent value(s) pass a small window (1 for
 * Nest/Fold/FixedPoint): once the buffer exceeds the window the oldest entries
 * are freed mid-loop, so a 10^6-step Nest holds O(window) live nodes instead of
 * O(n) — which keeps the bounded Expr node pool effective (an O(n) live set
 * empties the free-list and forces every allocation to fresh malloc). Steps
 * whose logic depends on buf->count (NestWhile's m_min/m_max gating) MUST pass
 * SIZE_MAX so the count keeps reflecting the true application total.
 */
static IterRunResult iter_run(ExprBuf* buf, IterStepFn step_fn, void* ctx,
                              int64_t limit, bool limit_is_safety,
                              size_t max_history, Expr** out_early) {
    *out_early = NULL;
    int64_t apps = 0;
    while (1) {
        if (apps >= limit) {
            if (limit_is_safety) {
                ebuf_free_all(buf);
                return ITER_RUN_SAFETY;
            }
            return ITER_RUN_OK;
        }
        Expr* next = NULL;
        IterStep s = step_fn(buf, ctx, &next);
        /* Catch/Throw: an in-flight Throw produced by the step (via f applied
         * to the running value) short-circuits the whole iteration and
         * propagates as the result. Covers Nest/NestList/Fold/FoldList/
         * NestWhile(List)/FixedPoint(List) in one place. */
        if (next && eval_is_inflight_throw(next)) {
            ebuf_free_all(buf);
            *out_early = next;
            return ITER_RUN_EARLY;
        }
        if (s == ITER_STEP_CONT) {
            ebuf_push(buf, next);
            apps++;
            /* Bounded-history window: drop superseded iterates so scalar runs
             * stay O(window) in memory. Trimming every push means the buffer
             * only ever overshoots by one, so this is O(window) per step. */
            if (buf->count > max_history)
                ebuf_drop_front(buf, buf->count - max_history);
        } else if (s == ITER_STEP_HALT_ADD) {
            ebuf_push(buf, next);
            return ITER_RUN_OK;
        } else if (s == ITER_STEP_RETURN) {
            ebuf_free_all(buf);
            *out_early = next;
            return ITER_RUN_EARLY;
        } else { /* ITER_STEP_HALT */
            return ITER_RUN_OK;
        }
    }
}

/* ------------------- Nest / NestList ------------------- */

/*
 * Nest[f, expr, n]       -- apply f to expr n times and return the final value.
 * NestList[f, expr, n]   -- return the list {expr, f[expr], ..., f^n[expr]}.
 *
 * n must be a non-negative integer. Returns NULL (unevaluated) otherwise.
 */

typedef struct { Expr* f; } NestCtx;

static IterStep nest_step(ExprBuf* hist, void* vctx, Expr** out) {
    NestCtx* c = (NestCtx*)vctx;
    *out = apply_unary(c->f, hist->items[hist->count - 1]);
    return ITER_STEP_CONT;
}

static Expr* nest_impl(Expr* res, bool as_list) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;

    Expr* f = res->data.function.args[0];
    Expr* expr = res->data.function.args[1];
    Expr* n_expr = res->data.function.args[2];
    if (n_expr->type != EXPR_INTEGER) return NULL;
    int64_t n = n_expr->data.integer;
    if (n < 0) return NULL;

    /* Automatic numeric fast-path: a machine-real-arithmetic pure function
     * iterated over a machine real runs in compiled doubles (no per-step Expr
     * allocation). NULL means "not numeric-closed", fall through to the
     * interpreted loop. NestList keeps every iterate, so it still allocates the
     * result -- what it saves is the substitute+evaluate of the body per step. */
    {
        Expr* fast = as_list ? numloop_nestlist(f, expr, n) : numloop_nest(f, expr, n);
        if (fast) return fast;
    }

    ExprBuf buf;
    ebuf_init(&buf);
    ebuf_push(&buf, expr_copy(expr));

    NestCtx ctx = { .f = f };
    Expr* early = NULL;
    /* Scalar Nest reads only the latest iterate -> window 1; NestList keeps all. */
    IterRunResult r = iter_run(&buf, nest_step, &ctx, n, false,
                               as_list ? SIZE_MAX : 1, &early);
    if (r == ITER_RUN_SAFETY) return NULL;
    if (r == ITER_RUN_EARLY) return early;
    return ebuf_finalize(&buf, as_list, expr_new_symbol(SYM_List));
}

Expr* builtin_nest(Expr* res)     { return nest_impl(res, false); }
Expr* builtin_nestlist(Expr* res) { return nest_impl(res, true); }

/* ------------------- Fold / FoldList ------------------- */

/*
 * Fold[f, x, list]       -- last element of FoldList[f, x, list]
 * Fold[f, list]          -- Fold[f, First[list], Rest[list]]
 * FoldList[f, x, list]   -- {x, f[x, list[[1]]], f[f[x, list[[1]]], list[[2]]], ...}
 * FoldList[f, list]      -- {list[[1]], f[list[[1]], list[[2]]], ...}
 *
 * The head of list is preserved in the FoldList output. Empty-list
 * behaviour: FoldList[f, {}] -> {}, Fold[f, {}] is unevaluated.
 */

typedef struct {
    Expr* f;
    Expr** elems;    /* borrowed; elements to consume */
    size_t total;
    size_t idx;
} FoldCtx;

static IterStep fold_step(ExprBuf* hist, void* vctx, Expr** out) {
    FoldCtx* c = (FoldCtx*)vctx;
    *out = apply_binary(c->f, hist->items[hist->count - 1], c->elems[c->idx++]);
    return ITER_STEP_CONT;
}

static Expr* fold_impl(Expr* res, bool as_list) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2 && argc != 3) return NULL;

    /* Folding over an association folds over its values: rebuild the call with
     * the association replaced by Values[assoc] and re-evaluate. Fold collapses
     * to a scalar; FoldList[f, assoc] keeps the keys, pairing each with its
     * running result (n -> n). */
    size_t coll_idx = (argc == 3) ? 2 : 1;
    if (is_association(res->data.function.args[coll_idx])) {
        Expr* assoc = res->data.function.args[coll_idx];
        Expr** na = malloc(sizeof(Expr*) * argc);
        for (size_t i = 0; i < argc; i++)
            na[i] = (i == coll_idx) ? assoc_values_list(res->data.function.args[i])
                                    : expr_copy(res->data.function.args[i]);
        Expr* call = expr_new_function(expr_copy(res->data.function.head), na, argc);
        free(na);
        Expr* r = evaluate(call);
        expr_free(call);
        if (as_list && argc == 2) {              /* FoldList[f, assoc] -> keyed */
            Expr* keyed = assoc_rekey_from_list(assoc, r);
            if (keyed) { expr_free(r); return keyed; }
        }
        return r;
    }

    /* An associative machine operator over a packed argument is a two-line C
     * loop — the running max/min/sum/product. Taken before the compiled VM
     * below because the VM cannot express Max/Min at all (there is no OP_MAX)
     * and pays per-element dispatch for the two it can. See ndreduce.h. */
    {
        Expr* coll = res->data.function.args[coll_idx];
        NDScanOp sop;
        if (is_ndarray(coll) && ndred_scan_op_for(res->data.function.args[0], &sop)) {
            Expr* fast = ndred_scan(coll, sop,
                                    argc == 3 ? res->data.function.args[1] : NULL,
                                    as_list);
            if (fast) return fast;
        }
    }

    /* Automatic numeric fast-path for the seeded Fold[f, x0, list] over machine
     * numbers: run the binary reduction in compiled doubles. FoldList takes the
     * same path and keeps every partial result. */
    if (argc == 3) {
        Expr* fast = as_list
            ? numloop_foldlist(res->data.function.args[0],
                               res->data.function.args[1],
                               res->data.function.args[2])
            : numloop_fold(res->data.function.args[0],
                           res->data.function.args[1],
                           res->data.function.args[2]);
        if (fast) return fast;
    }

    Expr* f = res->data.function.args[0];
    Expr* seed_src;
    Expr* list;
    bool seed_from_list;

    if (argc == 3) {
        seed_src = res->data.function.args[1];
        list = res->data.function.args[2];
        seed_from_list = false;
    } else {
        seed_src = NULL;
        list = res->data.function.args[1];
        seed_from_list = true;
    }

    /* An NDArray is atomic, so the element walk below would never see into it
     * and Fold[f, x, NDArray[...]] came back UNEVALUATED while the same call on
     * the equivalent List folded fine.  Materialise once and reuse the general
     * path: the answer is then identical to the List call by construction,
     * which is the contract the rest of the ND layer keeps (ndstruct.h:14).
     * `materialized` is freed on every exit below; nothing returned aliases it,
     * because both the seed and each element are copied on the way into the
     * history buffer. */
    Expr* materialized = NULL;
    const Expr* list_src = NULL;      /* the array we materialised, for the repack */
    NDType list_src_dtype = NDT_FLOAT64;
    if (is_ndarray(list)) {
        list_src = list;
        list_src_dtype = list->data.ndarray.dtype;
        materialized = ndarray_to_nested_list(list);
        if (!materialized) return NULL;
        list = materialized;
    }

    if (list->type != EXPR_FUNCTION) { expr_free(materialized); return NULL; }

    Expr* list_head = list->data.function.head;
    Expr** elems = list->data.function.args;
    size_t n = list->data.function.arg_count;
    size_t start = 0;

    if (seed_from_list) {
        if (n == 0) {
            Expr* empty = as_list ? expr_new_function(expr_copy(list_head), NULL, 0)
                                  : NULL;               /* Fold[f, {}] stays unevaluated */
            expr_free(materialized);
            return empty;
        }
        seed_src = elems[0];
        start = 1;
    }

    size_t m = n - start;

    ExprBuf buf;
    ebuf_init(&buf);
    ebuf_push(&buf, expr_copy(seed_src));

    FoldCtx ctx = { .f = f, .elems = elems + start, .total = m, .idx = 0 };
    Expr* early = NULL;
    /* Scalar Fold reads only the latest accumulator -> window 1; FoldList keeps all. */
    IterRunResult r = iter_run(&buf, fold_step, &ctx, (int64_t)m, false,
                               as_list ? SIZE_MAX : 1, &early);
    Expr* out;
    if (r == ITER_RUN_SAFETY)     out = NULL;
    else if (r == ITER_RUN_EARLY) out = early;
    else                          out = ebuf_finalize(&buf, as_list, expr_copy(list_head));

    /* A packed argument gives a packed history back, matching Map[f, NDArray[...]]. */
    if (materialized && as_list && out && r == ITER_RUN_OK)
        out = map_try_repack(list_src, out, list_src_dtype);
    expr_free(materialized);
    return out;
}

Expr* builtin_fold(Expr* res)     { return fold_impl(res, false); }
Expr* builtin_foldlist(Expr* res) { return fold_impl(res, true); }

/* ------------------- NestWhile / NestWhileList ------------------- */

/*
 * NestWhile / NestWhileList accept up to six arguments:
 *   (f, expr, test[, m[, max[, n]]])
 *
 * m    : positive integer, All, or {mmin, mmax|Infinity} -- how many recent
 *        results to pass to test (default 1).
 * max  : non-negative integer or Infinity -- cap on f-applications.
 * n    : integer -- positive means apply f n extra times after the while
 *        loop ends; negative means drop the last |n| iterates.
 *
 * Returns NULL (unevaluated) on malformed argument specs.
 */

typedef struct {
    Expr* f;
    Expr* test;
    int64_t m_min, m_max;
    bool m_max_inf;
} NestWhileCtx;

/*
 * Mathematica-style NestWhile semantics: test is called with the most recent
 * min(count, m_max) history entries, starting only once we have at least
 * m_min entries. A False result halts before applying f again.
 */
static IterStep nestwhile_step(ExprBuf* hist, void* vctx, Expr** out) {
    NestWhileCtx* c = (NestWhileCtx*)vctx;

    if ((int64_t)hist->count >= c->m_min) {
        int64_t k = c->m_max_inf ? (int64_t)hist->count
                                 : ((int64_t)hist->count < c->m_max
                                        ? (int64_t)hist->count : c->m_max);
        Expr** ta = malloc(sizeof(Expr*) * (size_t)k);
        for (int64_t i = 0; i < k; i++) {
            ta[i] = expr_copy(hist->items[hist->count - (size_t)k + (size_t)i]);
        }
        Expr* test_call = expr_new_function(expr_copy(c->test), ta, (size_t)k);
        free(ta);
        Expr* tr = eval_and_free(test_call);
        /* A Throw from the while-test propagates as the result. */
        if (eval_is_inflight_throw(tr)) { *out = tr; return ITER_STEP_RETURN; }
        bool keep = sym_is_true(tr);
        expr_free(tr);
        if (!keep) return ITER_STEP_HALT;
    }

    *out = apply_unary(c->f, hist->items[hist->count - 1]);
    return ITER_STEP_CONT;
}

static Expr* nestwhile_impl(Expr* res, bool as_list) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 3 || argc > 6) return NULL;

    Expr* f = res->data.function.args[0];
    Expr* expr = res->data.function.args[1];
    Expr* test = res->data.function.args[2];

    /* Automatic numeric fast-path for the default NestWhile[f, x0, test]
     * (m = 1, no max, no extra n): iterate in compiled doubles. */
    if (argc == 3) {
        Expr* fast = as_list ? numloop_nestwhilelist(f, expr, test)
                             : numloop_nestwhile(f, expr, test);
        if (fast) return fast;
    }

    /* Parse optional m. */
    int64_t m_min = 1, m_max = 1;
    bool m_max_inf = false;
    if (argc >= 4) {
        Expr* m_arg = res->data.function.args[3];
        if (m_arg->type == EXPR_INTEGER) {
            if (m_arg->data.integer < 1) return NULL;
            m_min = m_max = m_arg->data.integer;
        } else if (m_arg->type == EXPR_SYMBOL && m_arg->data.symbol.name == SYM_All) {
            m_min = 1;
            m_max_inf = true;
        } else if (m_arg->type == EXPR_FUNCTION &&
                   m_arg->data.function.head->type == EXPR_SYMBOL &&
                   m_arg->data.function.head->data.symbol.name == SYM_List &&
                   m_arg->data.function.arg_count == 2) {
            Expr* a0 = m_arg->data.function.args[0];
            Expr* a1 = m_arg->data.function.args[1];
            if (a0->type != EXPR_INTEGER || a0->data.integer < 1) return NULL;
            m_min = a0->data.integer;
            if (a1->type == EXPR_INTEGER) {
                if (a1->data.integer < m_min) return NULL;
                m_max = a1->data.integer;
            } else if (a1->type == EXPR_SYMBOL && a1->data.symbol.name == SYM_Infinity) {
                m_max_inf = true;
            } else {
                return NULL;
            }
        } else {
            return NULL;
        }
    }

    /* Parse optional max. */
    int64_t max_apps = 0;
    bool max_inf = true;
    if (argc >= 5) {
        Expr* max_arg = res->data.function.args[4];
        if (max_arg->type == EXPR_INTEGER) {
            if (max_arg->data.integer < 0) return NULL;
            max_apps = max_arg->data.integer;
            max_inf = false;
        } else if (max_arg->type == EXPR_SYMBOL && max_arg->data.symbol.name == SYM_Infinity) {
            max_inf = true;
        } else {
            return NULL;
        }
    }

    /* Parse optional n_extra. */
    int64_t n_extra = 0;
    if (argc >= 6) {
        Expr* n_arg = res->data.function.args[5];
        if (n_arg->type != EXPR_INTEGER) return NULL;
        n_extra = n_arg->data.integer;
    }

    ExprBuf buf;
    ebuf_init(&buf);
    ebuf_push(&buf, expr_copy(expr));

    NestWhileCtx ctx = { .f = f, .test = test,
                         .m_min = m_min, .m_max = m_max, .m_max_inf = m_max_inf };
    Expr* early = NULL;
    int64_t limit = max_inf ? ITER_SAFETY_CAP : max_apps;
    /* NestWhile's step gates on buf->count (m_min/m_max) and n_extra<0 trims the
     * tail, so the full history must be retained -> SIZE_MAX (no windowing). */
    IterRunResult r = iter_run(&buf, nestwhile_step, &ctx, limit, max_inf,
                               SIZE_MAX, &early);
    if (r == ITER_RUN_SAFETY) return NULL;
    if (r == ITER_RUN_EARLY) return early;

    /* Post-processing: positive n_extra appends more applications;
     * negative n_extra trims iterates from the end. */
    if (n_extra > 0) {
        NestCtx nctx = { .f = f };
        r = iter_run(&buf, nest_step, &nctx, n_extra, false, SIZE_MAX, &early);
        if (r == ITER_RUN_SAFETY) return NULL;
        if (r == ITER_RUN_EARLY) return early;
    } else if (n_extra < 0) {
        ebuf_truncate(&buf, (size_t)(-n_extra));
    }

    return ebuf_finalize(&buf, as_list, expr_new_symbol(SYM_List));
}

Expr* builtin_nestwhile(Expr* res)     { return nestwhile_impl(res, false); }
Expr* builtin_nestwhilelist(Expr* res) { return nestwhile_impl(res, true); }

/* ------------------- FixedPoint / FixedPointList ------------------- */

/*
 * FixedPoint / FixedPointList iterate f starting from expr until two
 * successive results are SameQ (or, when given SameTest -> s, until s
 * yields True on the consecutive pair). The result list always begins
 * with expr and ends with the fixed point; in the *List* form the fixed
 * point appears as the last two elements.
 *
 * Inside FixedPoint, a Throw/Abort/Quit/Return emitted by f propagates
 * as the result (mirrors Mathematica behaviour). The *List* form does
 * not intercept control-flow heads and will include them in the history
 * as ordinary values.
 *
 * Supported argument shapes (after f, expr):
 *     -- nothing (iterate until fixed point / safety cap)
 *     -- n                                (bound on applications)
 *     -- SameTest -> s
 *     -- n, SameTest -> s
 */

typedef struct {
    Expr* f;
    Expr* same_test;       /* NULL => use SameQ (expr_eq) */
    bool propagate_throw;  /* FixedPoint only */
} FixedPointCtx;

static IterStep fixedpoint_step(ExprBuf* hist, void* vctx, Expr** out) {
    FixedPointCtx* c = (FixedPointCtx*)vctx;
    Expr* last = hist->items[hist->count - 1];
    Expr* next = apply_unary(c->f, last);

    if (c->propagate_throw && next && next->type == EXPR_FUNCTION &&
        next->data.function.head->type == EXPR_SYMBOL) {
        const char* h = next->data.function.head->data.symbol.name;
        if (h == SYM_Throw || h == SYM_Abort ||
            h == SYM_Quit || h == SYM_Return) {
            *out = next;
            return ITER_STEP_RETURN;
        }
    }

    bool same;
    if (c->same_test == NULL) {
        same = expr_eq(last, next);
    } else {
        Expr* tr = apply_binary(c->same_test, last, next);
        /* A Throw from SameTest propagates as the result. */
        if (eval_is_inflight_throw(tr)) { expr_free(next); *out = tr; return ITER_STEP_RETURN; }
        same = sym_is_true(tr);
        expr_free(tr);
    }

    *out = next;
    return same ? ITER_STEP_HALT_ADD : ITER_STEP_CONT;
}

/*
 * Parse optional n + SameTest -> s arguments, starting at index `start`.
 * Returns true on success, false (leaving the builtin unevaluated) on
 * malformed or duplicated specs.
 */
static bool parse_fp_opts(Expr* res, size_t start,
                          Expr** out_same_test,
                          int64_t* out_max, bool* out_max_inf) {
    Expr* same_test = NULL;
    bool has_max = false;
    int64_t max_apps = 0;
    bool max_inf = true;
    size_t argc = res->data.function.arg_count;

    for (size_t i = start; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (a->type == EXPR_FUNCTION &&
            a->data.function.head->type == EXPR_SYMBOL &&
            (a->data.function.head->data.symbol.name == SYM_Rule ||
             a->data.function.head->data.symbol.name == SYM_RuleDelayed) &&
            a->data.function.arg_count == 2 &&
            a->data.function.args[0]->type == EXPR_SYMBOL &&
            a->data.function.args[0]->data.symbol.name == SYM_SameTest) {
            if (same_test != NULL) return false;
            same_test = a->data.function.args[1];
        } else if (!has_max && a->type == EXPR_INTEGER) {
            if (a->data.integer < 0) return false;
            max_apps = a->data.integer;
            max_inf = false;
            has_max = true;
        } else if (!has_max && a->type == EXPR_SYMBOL &&
                   a->data.symbol.name == SYM_Infinity) {
            has_max = true;
            max_inf = true;
        } else {
            return false;
        }
    }

    *out_same_test = same_test;
    *out_max = max_apps;
    *out_max_inf = max_inf;
    return true;
}

static Expr* fixedpoint_impl(Expr* res, bool as_list) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    Expr* f = res->data.function.args[0];
    Expr* expr = res->data.function.args[1];

    Expr* same_test = NULL;
    int64_t max_apps = 0;
    bool max_inf = true;
    if (!parse_fp_opts(res, 2, &same_test, &max_apps, &max_inf)) return NULL;

    /* Automatic numeric fast-path for the plain FixedPoint[f, x0] (no SameTest,
     * no application cap): iterate in compiled doubles. */
    if (same_test == NULL && max_inf) {
        Expr* fast = as_list ? numloop_fixedpointlist(f, expr)
                             : numloop_fixedpoint(f, expr);
        if (fast) return fast;
    }

    ExprBuf buf;
    ebuf_init(&buf);
    ebuf_push(&buf, expr_copy(expr));

    FixedPointCtx ctx = { .f = f, .same_test = same_test,
                          .propagate_throw = !as_list };
    Expr* early = NULL;
    int64_t limit = max_inf ? ITER_SAFETY_CAP : max_apps;
    /* Scalar FixedPoint compares only the latest pair -> window 1; the *List
     * form keeps every iterate. */
    IterRunResult r = iter_run(&buf, fixedpoint_step, &ctx, limit, max_inf,
                               as_list ? SIZE_MAX : 1, &early);
    if (r == ITER_RUN_SAFETY) return NULL;
    if (r == ITER_RUN_EARLY) return early;

    return ebuf_finalize(&buf, as_list, expr_new_symbol(SYM_List));
}

Expr* builtin_fixedpoint(Expr* res)     { return fixedpoint_impl(res, false); }
Expr* builtin_fixedpointlist(Expr* res) { return fixedpoint_impl(res, true); }

/* ------------------- Thread ------------------- */

/*
 * Parse a Thread position specification into a boolean mask of length K.
 *
 *   spec    | meaning
 *   --------+--------------------------------------------------------------
 *   NULL    | All (default when omitted)
 *   All     | all elements
 *   None    | no elements
 *   n>=0    | first n elements                       (n>K is clamped to K)
 *   n<0     | last -n elements                       (-n>K is clamped to K)
 *   {n}     | only element n   (1-based, n<0 counts from the end)
 *   {m,n}   | elements m..n inclusive
 *   {m,n,s} | elements m..n in steps of s            (s != 0)
 *
 * Returns true on success and fills mask[0..K-1]. Returns false on a
 * malformed spec; mask contents are undefined in that case.
 */
static bool thread_parse_spec(Expr* spec, size_t K, bool* mask) {
    /* Default: All. */
    if (!spec) {
        for (size_t i = 0; i < K; i++) mask[i] = true;
        return true;
    }

    if (spec->type == EXPR_SYMBOL) {
        if (spec->data.symbol.name == SYM_All) {
            for (size_t i = 0; i < K; i++) mask[i] = true;
            return true;
        }
        if (spec->data.symbol.name == SYM_None) {
            for (size_t i = 0; i < K; i++) mask[i] = false;
            return true;
        }
        return false;
    }

    if (spec->type == EXPR_INTEGER) {
        int64_t n = spec->data.integer;
        for (size_t i = 0; i < K; i++) mask[i] = false;
        if (n >= 0) {
            size_t lim = ((size_t)n > K) ? K : (size_t)n;
            for (size_t i = 0; i < lim; i++) mask[i] = true;
        } else {
            int64_t cnt = -n;
            size_t lim = ((size_t)cnt > K) ? K : (size_t)cnt;
            for (size_t i = K - lim; i < K; i++) mask[i] = true;
        }
        return true;
    }

    if (spec->type == EXPR_FUNCTION &&
        spec->data.function.head->type == EXPR_SYMBOL &&
        spec->data.function.head->data.symbol.name == SYM_List) {
        size_t na = spec->data.function.arg_count;
        int64_t m_idx, n_idx, s_step;
        if (na == 1) {
            Expr* a0 = spec->data.function.args[0];
            if (a0->type != EXPR_INTEGER) return false;
            m_idx = n_idx = a0->data.integer;
            s_step = 1;
        } else if (na == 2) {
            Expr* a0 = spec->data.function.args[0];
            Expr* a1 = spec->data.function.args[1];
            if (a0->type != EXPR_INTEGER || a1->type != EXPR_INTEGER) return false;
            m_idx = a0->data.integer;
            n_idx = a1->data.integer;
            s_step = 1;
        } else if (na == 3) {
            Expr* a0 = spec->data.function.args[0];
            Expr* a1 = spec->data.function.args[1];
            Expr* a2 = spec->data.function.args[2];
            if (a0->type != EXPR_INTEGER || a1->type != EXPR_INTEGER || a2->type != EXPR_INTEGER) return false;
            m_idx = a0->data.integer;
            n_idx = a1->data.integer;
            s_step = a2->data.integer;
        } else {
            return false;
        }
        if (s_step == 0) return false;

        /* Negative indices count from the end: -1 -> K, -2 -> K-1, ... */
        if (m_idx < 0) m_idx = (int64_t)K + 1 + m_idx;
        if (n_idx < 0) n_idx = (int64_t)K + 1 + n_idx;

        for (size_t i = 0; i < K; i++) mask[i] = false;

        if (s_step > 0) {
            for (int64_t i = m_idx; i <= n_idx; i += s_step) {
                if (i >= 1 && i <= (int64_t)K) mask[i - 1] = true;
            }
        } else {
            for (int64_t i = m_idx; i >= n_idx; i += s_step) {
                if (i >= 1 && i <= (int64_t)K) mask[i - 1] = true;
            }
        }
        return true;
    }

    return false;
}

Expr* builtin_thread(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return NULL;

    Expr* expr = res->data.function.args[0];

    /* Thread on an atom -- nothing to thread, return the atom. */
    if (expr->type != EXPR_FUNCTION) return expr_copy(expr);

    Expr* h_spec = (argc >= 2) ? res->data.function.args[1] : NULL;
    Expr* n_spec = (argc >= 3) ? res->data.function.args[2] : NULL;

    /* Default threading head is List. We use the interned symbol so
     * pointer comparisons in expr_eq work. */
    Expr* h = h_spec ? expr_ref(h_spec) : expr_new_symbol(SYM_List);

    size_t K = expr->data.function.arg_count;
    bool* mask = (K > 0) ? (bool*)malloc(sizeof(bool) * K) : NULL;
    if (K > 0 && !mask) { expr_free(h); return NULL; }

    if (!thread_parse_spec(n_spec, K, mask)) {
        if (mask) free(mask);
        expr_free(h);
        return NULL;
    }

    /* Determine threading length L from threadable args whose head is h.
     * All such args must have the same length, otherwise we leave the
     * expression unchanged (matches Mathematica's "unequal length" path
     * after issuing the message; we simply skip the message). */
    int64_t L = -1;
    bool has_threadable_h = false;
    for (size_t i = 0; i < K; i++) {
        if (!mask[i]) continue;
        Expr* a = expr->data.function.args[i];
        if (a->type == EXPR_FUNCTION && expr_eq(a->data.function.head, h)) {
            has_threadable_h = true;
            int64_t len = (int64_t)a->data.function.arg_count;
            if (L == -1) {
                L = len;
            } else if (len != L) {
                if (mask) free(mask);
                expr_free(h);
                return expr_copy(expr);
            }
        }
    }

    if (!has_threadable_h) {
        /* Nothing in the selected positions matches h; return expr unchanged. */
        if (mask) free(mask);
        expr_free(h);
        return expr_copy(expr);
    }

    Expr* f = expr->data.function.head;

    /* Build h[ f[...], f[...], ..., f[...] ] with L copies. */
    Expr** wrap_args = (L > 0) ? (Expr**)malloc(sizeof(Expr*) * (size_t)L) : NULL;
    if (L > 0 && !wrap_args) {
        if (mask) free(mask);
        expr_free(h);
        return NULL;
    }

    for (int64_t k = 0; k < L; k++) {
        Expr** new_args = (K > 0) ? (Expr**)malloc(sizeof(Expr*) * K) : NULL;
        for (size_t j = 0; j < K; j++) {
            Expr* aj = expr->data.function.args[j];
            if (mask[j] && aj->type == EXPR_FUNCTION &&
                expr_eq(aj->data.function.head, h)) {
                new_args[j] = expr_copy(aj->data.function.args[(size_t)k]);
            } else {
                new_args[j] = expr_copy(aj);
            }
        }
        wrap_args[k] = expr_new_function(expr_copy(f), new_args, K);
        if (new_args) free(new_args);
    }

    Expr* wrapped = expr_new_function(expr_copy(h), wrap_args, (size_t)L);
    if (wrap_args) free(wrap_args);
    if (mask) free(mask);
    expr_free(h);

    /* Evaluate so f attributes (Listable, OneIdentity, ...) take effect. */
    Expr* eval_res = evaluate(wrapped);
    expr_free(wrapped);
    return eval_res;
}

/* ---------------------------------------------------------------------------
 * MapThread — the "function and arguments separately" sibling of Thread.
 *
 *   MapThread[f, {l1, ..., lk}]      -> {f[l1[[1]],...], ..., f[l1[[L]],...]}
 *   MapThread[f, {e1, ..., ek}, n]   -> thread f over the parts at level n
 *
 * mapthread_rec threads f in parallel through the k BORROWED sub-expressions in
 * `exprs`, descending `level` dimensions. At level 0 it builds the leaf
 * f[exprs[0], ..., exprs[k-1]] (deep-copying each part). At level > 0 all k
 * sub-expressions must have the same shape: either all Lists of equal length
 * (thread over corresponding elements -> List[...]) or all Associations with
 * identical key sequences (thread over values, preserving keys -> Association).
 * Any structural mismatch returns NULL so the caller leaves MapThread[...]
 * unevaluated (matching builtin_thread's silent unequal-length path).
 * -------------------------------------------------------------------------- */
static Expr* mapthread_rec(Expr* f, Expr** exprs, size_t k, int64_t level) {
    /* Leaf: apply f to the k corresponding parts. */
    if (level <= 0) {
        Expr** fargs = (Expr**)malloc(sizeof(Expr*) * (k ? k : 1));
        if (!fargs) return NULL;
        for (size_t j = 0; j < k; j++) fargs[j] = expr_copy(exprs[j]);
        Expr* leaf = expr_new_function(expr_copy(f), fargs, k);
        free(fargs);
        return leaf;
    }

    if (k == 0) return NULL;   /* cannot determine shape with no sub-exprs */

    /* Classify the k sub-expressions: all associations, or all lists? */
    bool all_assoc = true, all_list = true;
    for (size_t j = 0; j < k; j++) {
        if (!is_association(exprs[j])) all_assoc = false;
        if (!(exprs[j]->type == EXPR_FUNCTION &&
              exprs[j]->data.function.head->type == EXPR_SYMBOL &&
              exprs[j]->data.function.head->data.symbol.name == SYM_List))
            all_list = false;
    }

    /* Association threading: identical key sequences, thread over values. */
    if (all_assoc) {
        size_t L = exprs[0]->data.function.arg_count;
        for (size_t j = 1; j < k; j++)
            if (exprs[j]->data.function.arg_count != L) return NULL;
        for (size_t p = 0; p < L; p++) {
            Expr* key0 = exprs[0]->data.function.args[p]->data.function.args[0];
            for (size_t j = 1; j < k; j++) {
                Expr* keyj = exprs[j]->data.function.args[p]->data.function.args[0];
                if (!expr_eq(key0, keyj)) return NULL;
            }
        }
        Expr** rules = (Expr**)malloc(sizeof(Expr*) * (L ? L : 1));
        if (!rules) return NULL;
        Expr** col = (Expr**)malloc(sizeof(Expr*) * k);
        if (!col) { free(rules); return NULL; }
        for (size_t p = 0; p < L; p++) {
            for (size_t j = 0; j < k; j++)
                col[j] = exprs[j]->data.function.args[p]->data.function.args[1];
            Expr* child = mapthread_rec(f, col, k, level - 1);
            if (!child) {
                for (size_t q = 0; q < p; q++) expr_free(rules[q]);
                free(rules); free(col);
                return NULL;
            }
            Expr* key = expr_copy(exprs[0]->data.function.args[p]->data.function.args[0]);
            Expr* rargs[2] = { key, child };
            rules[p] = expr_new_function(expr_new_symbol(SYM_Rule), rargs, 2);
        }
        Expr* result = expr_new_function(expr_new_symbol(SYM_Association), rules, L);
        free(rules); free(col);
        return result;
    }

    /* List threading: all lists of equal length, thread over elements. */
    if (all_list) {
        size_t L = exprs[0]->data.function.arg_count;
        for (size_t j = 1; j < k; j++)
            if (exprs[j]->data.function.arg_count != L) return NULL;
        Expr** out = (Expr**)malloc(sizeof(Expr*) * (L ? L : 1));
        if (!out) return NULL;
        Expr** col = (Expr**)malloc(sizeof(Expr*) * k);
        if (!col) { free(out); return NULL; }
        for (size_t i = 0; i < L; i++) {
            for (size_t j = 0; j < k; j++)
                col[j] = exprs[j]->data.function.args[i];
            Expr* child = mapthread_rec(f, col, k, level - 1);
            if (!child) {
                for (size_t q = 0; q < i; q++) expr_free(out[q]);
                free(out); free(col);
                return NULL;
            }
            out[i] = child;
        }
        Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, L);
        free(out); free(col);
        return result;
    }

    return NULL;  /* mixed / non-threadable shapes -> leave unevaluated */
}

/* MapThread[f, arr] over a rank-2 float64 buffer: fold f down each column.
 *
 * MapThread[f, {l1, ..., lk}] pairs the j-th element of every list, so over a
 * k x n buffer it is exactly a reduction along the LEADING axis -- the same
 * shape of work Total[a, {1}] does, for a different operator. Restricted to the
 * associative n-ary heads, which is what "fold down the column" means at all.
 *
 * The generic path builds one f[x1, ..., xk] node per column and evaluates the
 * lot: 365 ms for a 16 x 100000 buffer, which was more than half of the k-means
 * benchmark (the assignment step takes an elementwise minimum across k centroid
 * distance rows, once per iteration).
 *
 * Accumulation is in row order, matching the left-to-right fold builtin_plus
 * does on the same arguments, so Plus and Times are bit-identical to the List
 * path rather than merely close. */
static Expr* nd_mapthread2(const Expr* f, const Expr* a) {
    if (!f || f->type != EXPR_SYMBOL) return NULL;
    if (a->type != EXPR_NDARRAY || a->data.ndarray.rank != 2 ||
        (a->data.ndarray.dtype != NDT_FLOAT64 &&
         a->data.ndarray.dtype != NDT_INT64)) return NULL;

    const char* op = f->data.symbol.name;
    enum { M_MIN, M_MAX, M_PLUS, M_TIMES } k;
    if      (op == SYM_Min)   k = M_MIN;
    else if (op == SYM_Max)   k = M_MAX;
    else if (op == SYM_Plus)  k = M_PLUS;
    else if (op == SYM_Times) k = M_TIMES;
    else return NULL;

    int64_t rows = a->data.ndarray.dims[0], n = a->data.ndarray.dims[1];
    if (rows < 1) return NULL;

    /* int64 arm. Min and Max hand back one of their arguments, so they are
     * exact on an integer buffer by construction; Plus and Times compute, so
     * they go through the checked helpers and ABANDON the whole result on the
     * first overflow -- the List path then re-runs it and GMP answers exactly,
     * rather than a wrapped sum reaching the user.
     *
     * Without this arm every integer MapThread materialised: the two rows of a
     * Needleman-Wunsch DP -- MapThread[Max, {diag, up}] over 2000 int64 columns,
     * once per row -- cost 853 ms per 10^6 columns where the identical Real
     * pair cost 12 ms. */
    if (a->data.ndarray.dtype == NDT_INT64) {
        const int64_t* A = (const int64_t*)a->data.ndarray.data;
        int64_t* out = malloc(sizeof(int64_t) * (size_t)(n > 0 ? n : 1));
        if (!out) return NULL;
        for (int64_t j = 0; j < n; j++) out[j] = A[j];
        for (int64_t i = 1; i < rows; i++) {
            const int64_t* row = A + (size_t)i * (size_t)n;
            switch (k) {
            case M_MIN: for (int64_t j = 0; j < n; j++) if (row[j] < out[j]) out[j] = row[j]; break;
            case M_MAX: for (int64_t j = 0; j < n; j++) if (row[j] > out[j]) out[j] = row[j]; break;
            case M_PLUS:
                for (int64_t j = 0; j < n; j++)
                    if (ci_add_i64(out[j], row[j], &out[j])) { free(out); return NULL; }
                break;
            case M_TIMES:
                for (int64_t j = 0; j < n; j++)
                    if (ci_mul_i64(out[j], row[j], &out[j])) { free(out); return NULL; }
                break;
            }
        }
        int64_t idims[1] = { n };
        return expr_new_ndarray_like(a, 1, idims, out, NDT_INT64);
    }

    const double* A = (const double*)a->data.ndarray.data;
    double* out = malloc(sizeof(double) * (size_t)(n > 0 ? n : 1));
    if (!out) return NULL;
    for (int64_t j = 0; j < n; j++) out[j] = A[j];
    for (int64_t i = 1; i < rows; i++) {
        const double* row = A + (size_t)i * (size_t)n;
        switch (k) {
        case M_MIN:   for (int64_t j = 0; j < n; j++) if (row[j] < out[j]) out[j] = row[j]; break;
        case M_MAX:   for (int64_t j = 0; j < n; j++) if (row[j] > out[j]) out[j] = row[j]; break;
        case M_PLUS:  for (int64_t j = 0; j < n; j++) out[j] += row[j]; break;
        case M_TIMES: for (int64_t j = 0; j < n; j++) out[j] *= row[j]; break;
        }
    }
    for (int64_t j = 0; j < n; j++)
        if (!isfinite(out[j])) { free(out); return NULL; }  /* let the List path form it */

    int64_t dims[1] = { n };
    return expr_new_ndarray_like(a, 1, dims, out, NDT_FLOAT64);
}

Expr* builtin_mapthread(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3)
        return builtin_arg_error("MapThread", argc, 2, 3);

    Expr* f     = res->data.function.args[0];
    Expr* lists = res->data.function.args[1];

    if (argc == 2) {                          /* level 1, the default */
        Expr* fast = nd_mapthread2(f, lists);
        if (fast) return fast;
    }

    int64_t level = 1;                       /* default threading level */
    if (argc == 3) {
        Expr* nspec = res->data.function.args[2];
        if (nspec->type != EXPR_INTEGER || nspec->data.integer < 0) return NULL;
        level = nspec->data.integer;
    }

    /* Outer container must be a List of the k expressions to thread. */
    if (lists->type != EXPR_FUNCTION ||
        lists->data.function.head->type != EXPR_SYMBOL ||
        lists->data.function.head->data.symbol.name != SYM_List)
        return NULL;

    size_t k = lists->data.function.arg_count;
    if (k == 0)                              /* MapThread[f, {}] -> {} */
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);

    /* Materialize any NDArray entry to a nested list so it threads like a list:
     * MapThread[f, {NDArray[{1,2}], NDArray[{3,4}]}] -> {f[1,3], f[2,4]}. */
    Expr** entries = lists->data.function.args;
    Expr* const* entries0 = lists->data.function.args;   /* pre-delisting, for the dtype */
    Expr** delisted = NULL;
    bool any_nd = false;
    for (size_t j = 0; j < k; j++)
        if (is_ndarray(entries[j])) { any_nd = true; break; }
    if (any_nd) {
        delisted = malloc(sizeof(Expr*) * k);
        for (size_t j = 0; j < k; j++)
            delisted[j] = is_ndarray(entries[j]) ? ndarray_to_nested_list(entries[j])
                                                 : expr_copy(entries[j]);
        entries = delisted;
    }

    Expr* built = mapthread_rec(f, entries, k, level);
    if (delisted) {
        for (size_t j = 0; j < k; j++) expr_free(delisted[j]);
        free(delisted);
    }
    if (!built) return NULL;                 /* structural mismatch -> unchanged */

    /* Evaluate so f's attributes fire and the f[...] leaves reduce
     * (mirrors builtin_thread's final evaluate). */
    Expr* out = evaluate(built);
    expr_free(built);
    /* Packed in, packed out — the same convention Map keeps, so threading over
     * NDArrays does not silently unpack the data. */
    if (any_nd && out) {
        NDType dt = NDT_FLOAT64;
        const Expr* src_nd = NULL;
        for (size_t j = 0; j < k; j++)
            if (is_ndarray(entries0[j])) { src_nd = entries0[j]; dt = src_nd->data.ndarray.dtype; break; }
        out = map_try_repack(src_nd, out, dt);
    }
    return out;
}
