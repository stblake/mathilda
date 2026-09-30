/* ===========================================================================
 * Query and Dataset -- the Wolfram Language query language (see assoc_query.h)
 *
 * Query[op1, op2, ..., opn][data] applies op_i at level i of data.  The rule
 * that makes this more than "Map at successive levels" is Mathematica's split
 * of operators into DESCENDING and ASCENDING ones:
 *
 *   * A descending operator is applied to the ORIGINAL data at its level,
 *     before the deeper operators run.  Part specifications ("key", Key[k],
 *     i, i;;j, {p1, p2, ...}, All) and filtering / reordering operators
 *     (Select, SortBy, GroupBy, MaximalBy, ... Reverse, Sort, Keys, Values,
 *     DeleteMissing) are descending: they leave the structure below them
 *     intact, so the deeper operators still see records.
 *
 *   * An ascending operator is applied AFTER every deeper operator has run,
 *     to the results.  Everything else -- aggregations (Total, Length, Max,
 *     Mean, Counts, First, ...) and arbitrary functions -- is ascending.
 *
 * So Query[Total, "a"] maps "a" over the records first and totals the result,
 * whereas Query[Select[p], "a"] filters the records first and then extracts
 * "a" from the survivors.  The classification table below was read off
 * Mathematica 15's own compiled form: Normal[Query[op, h]] gives
 * RightComposition[op, Map[h]] for a descending op and
 * RightComposition[Map[h], op] for an ascending one.  q_normal reproduces
 * that compiled form.
 *
 * A single-part descending operator ("key", Key[k], integer i, SelectFirst,
 * Lookup[k]) consumes its level: the next operator is applied to its result
 * directly rather than mapped over it (RightComposition[Slice["a"], h]).
 *
 * Missing propagation: a part that does not exist gives Missing["KeyAbsent",
 * k] (absent key), Missing["PartAbsent", i] (index out of range) or
 * Missing["PartInvalid", p] (part of an atom / key of a list); a Missing flows
 * through later part specifications unchanged, and the numeric aggregations
 * Total/Mean/Max/Min/Median/Variance/StandardDeviation drop Missing values
 * (Mathematica's MissingBehavior -> Automatic).
 *
 * Query applies operator forms ITSELF (Select[p] on x becomes Select[x, p],
 * Map[f] becomes Map[f, x], ...), so it does not depend on the evaluator
 * supporting curried forms.  An operator it does not know is applied as op[x]
 * and left to the evaluator.
 *
 * Dataset[data] is a thin wrapper: ds[ops...] runs Query[ops...] on the data
 * and re-wraps a List or Association result as a Dataset; a scalar (or
 * Missing) comes back bare, as in Mathematica.  A handful of list builtins
 * are made Dataset-transparent by wrapping their registered C functions
 * (install_dataset_wrappers), so no other module needs to know Dataset exists.
 * ======================================================================== */

#include "assoc_query.h"
#include "assoc.h"
#include "common.h"
#include "eval.h"
#include "expr.h"
#include "ndarray.h"
#include "print.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "symtab.h"
#include "attr.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ names */

/* Heads the classifier compares against.  Interned once at init (pointer
 * identity), local to this module: they are operator names Query recognises,
 * not internal symbols of their own. */
static const char *Q_Select, *Q_SortBy, *Q_ReverseSortBy, *Q_GroupBy,
    *Q_MaximalBy, *Q_MinimalBy, *Q_DeleteDuplicatesBy, *Q_KeySortBy,
    *Q_KeySelect, *Q_KeyMap, *Q_TakeLargestBy, *Q_TakeSmallestBy,
    *Q_TakeLargest, *Q_TakeSmallest, *Q_SelectFirst, *Q_Lookup, *Q_KeyTake,
    *Q_KeyDrop, *Q_CountsBy, *Q_Merge, *Q_DeleteCases, *Q_Cases, *Q_Count,
    *Q_Position, *Q_MemberQ, *Q_FreeQ, *Q_KeyExistsQ, *Q_KeyMemberQ,
    *Q_KeyFreeQ, *Q_ReplaceAll, *Q_Replace, *Q_Append, *Q_Prepend,
    *Q_Extract, *Q_Delete, *Q_ReplacePart, *Q_Insert, *Q_AllTrue,
    *Q_AnyTrue, *Q_NoneTrue, *Q_FirstCase, *Q_FirstPosition,
    *Q_Map, *Q_Apply, *Q_MapIndexed, *Q_KeyValueMap, *Q_AssociationMap,
    *Q_Scan, *Q_MapAt, *Q_Reverse, *Q_Sort, *Q_ReverseSort, *Q_DeleteMissing,
    *Q_Keys, *Q_Values, *Q_Transpose, *Q_Total, *Q_Mean, *Q_Max, *Q_Min,
    *Q_Median, *Q_Variance, *Q_StandardDeviation, *Q_RightComposition,
    *Q_Composition, *Q_Slice, *Q_ApplyThrough, *Q_AssocTranspose,
    *Q_Identity, *Q_Normal, *Q_Length, *Q_Dimensions, *Q_First, *Q_Last,
    *Q_Take, *Q_Part;

/* ---------------------------------------------------------------- helpers */

static bool is_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == name;
}

static const char* head_name(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || !e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL) return NULL;
    return e->data.function.head->data.symbol.name;
}

static size_t nargs(const Expr* e) { return e->data.function.arg_count; }
static Expr* arg(const Expr* e, size_t i) { return e->data.function.args[i]; }

static bool is_list(const Expr* e) { return head_name(e) == SYM_List; }
static bool is_missing(const Expr* e) { return head_name(e) == SYM_Missing; }
static bool is_query(const Expr* e) { return head_name(e) == SYM_Query; }
static bool is_dataset(const Expr* e) {
    return head_name(e) == SYM_Dataset && nargs(e) >= 1;
}

static Expr* mk1(const char* h, Expr* a) {
    Expr* args[1] = { a };
    return expr_new_function(expr_new_symbol(h), args, 1);
}
static Expr* mk2(const char* h, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(h), args, 2);
}

static Expr* missing2(const char* why, Expr* what /* adopted */) {
    return mk2(SYM_Missing, expr_new_string(why), what);
}

/* A packed list presents as a List everywhere; the query walker reads
 * args[], so materialise it.  Returns an owned reference either way. */
static Expr* unpacked(Expr* e) {
    if (e && is_packed_list(e)) return ndarray_to_nested_list(e);
    return expr_copy(e);
}

/* "key" or Key[k] -- a key specification. */
static bool is_key_spec(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_STRING) return true;
    return head_name(e) == SYM_Key && nargs(e) == 1;
}
static const Expr* key_of(const Expr* e) {
    return e->type == EXPR_STRING ? e : arg(e, 0);
}
static bool is_index(const Expr* e) {
    return e && (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT);
}
static bool is_part_elem(const Expr* e) { return is_key_spec(e) || is_index(e); }

/* A List of part elements (possibly empty) -- a multi-part specification. */
static bool is_part_list(const Expr* e) {
    if (!is_list(e)) return false;
    for (size_t i = 0; i < nargs(e); i++)
        if (!is_part_elem(arg(e, i))) return false;
    return true;
}

static Expr* new_list(Expr** items, size_t n) {
    return expr_new_function(expr_new_symbol(SYM_List), items, n);
}

/* Evaluate an owned expression, freeing it. */
static Expr* ev(Expr* e) { return eval_and_free(e); }

/* ------------------------------------------------------- classification */

typedef enum {
    QK_SKIP,    /* Nothing: removed from the operator chain            */
    QK_ALL,     /* All / Identity: map the rest over this level         */
    QK_PART1,   /* "key", Key[k], i: take one part, rest applies to it  */
    QK_PARTN,   /* i;;j, {p1, ...}: take parts, map the rest over them  */
    QK_DESC,    /* descending operator: apply, then map the rest        */
    QK_DESC1,   /* descending, single result (SelectFirst, Lookup[k])   */
    QK_ASC      /* ascending: map the rest, then apply to the results   */
} QKind;

static QKind classify(const Expr* op) {
    if (!op) return QK_ASC;
    if (op->type == EXPR_SYMBOL) {
        const char* s = op->data.symbol.name;
        if (s == SYM_Nothing) return QK_SKIP;
        if (s == SYM_All || s == Q_Identity) return QK_ALL;
        if (s == Q_Reverse || s == Q_Sort || s == Q_ReverseSort ||
            s == Q_DeleteMissing || s == Q_Keys || s == Q_Values)
            return QK_DESC;
        return QK_ASC;
    }
    if (is_part_elem(op)) return QK_PART1;
    if (head_name(op) == SYM_Span) return QK_PARTN;
    if (is_part_list(op)) return QK_PARTN;
    {
        const char* h = head_name(op);
        if (h && nargs(op) >= 1) {
            if (h == Q_Select || h == Q_SortBy || h == Q_ReverseSortBy ||
                h == Q_GroupBy || h == Q_MaximalBy || h == Q_MinimalBy ||
                h == Q_DeleteDuplicatesBy || h == Q_KeySortBy ||
                h == Q_KeySelect || h == Q_KeyMap || h == Q_TakeLargestBy ||
                h == Q_TakeSmallestBy)
                return QK_DESC;
            if (h == Q_SelectFirst || h == Q_Lookup) return QK_DESC1;
        }
    }
    return QK_ASC;
}

/* ------------------------------------------------------------- slicing */

/* Resolve a 1-based (negative = from the end) index against length n.
 * Returns the 0-based position, or SIZE_MAX when out of range. */
static size_t resolve_index(const Expr* ie, size_t n) {
    int64_t i;
    if (ie->type != EXPR_INTEGER) return SIZE_MAX;
    i = ie->data.integer;
    if (i > 0 && (uint64_t)i <= n) return (size_t)(i - 1);
    if (i < 0 && (uint64_t)(-(i + 1)) < n) return (size_t)((int64_t)n + i);
    return SIZE_MAX;
}

/* The i-th element (List) or value (Association) of a container. */
static Expr* container_item(const Expr* c, size_t i) {
    Expr* a = arg(c, i);
    if (head_name(c) == SYM_Association) return arg(a, 1);
    return a;
}

/* One part: "key", Key[k] or an integer.  data is borrowed. */
static Expr* q_slice1(const Expr* op, Expr* data) {
    Expr* d;
    Expr* out;
    if (is_missing(data)) return expr_copy(data);
    d = unpacked(data);
    if (is_index(op)) {
        if (is_list(d) || is_association(d)) {
            size_t p = resolve_index(op, nargs(d));
            out = (p == SIZE_MAX) ? missing2("PartAbsent", expr_copy((Expr*)op))
                                  : expr_copy(container_item(d, p));
        } else {
            out = missing2("PartInvalid", expr_copy((Expr*)op));
        }
    } else {
        const Expr* k = key_of(op);
        if (is_association(d)) {
            Expr* v = assoc_lookup_value(d, k);
            out = v ? expr_copy(v) : missing2("KeyAbsent", expr_copy((Expr*)k));
        } else {
            out = missing2("PartInvalid", expr_copy((Expr*)k));
        }
    }
    expr_free(d);
    return out;
}

/* Resolve Span[a, b(, c)] against length n into start/stop/step (0-based,
 * inclusive).  Returns false when the span is malformed or out of range;
 * *empty is set for the legitimate empty span a;;a-1. */
static bool resolve_span(const Expr* sp, size_t n, int64_t* start,
                         int64_t* stop, int64_t* step, bool* empty) {
    int64_t N = (int64_t)n, a = 1, b = N, c = 1;
    size_t k = nargs(sp);
    *empty = false;
    if (k < 2 || k > 3) return false;
    if (!is_sym(arg(sp, 0), SYM_All)) {
        if (arg(sp, 0)->type != EXPR_INTEGER) return false;
        a = arg(sp, 0)->data.integer;
    }
    if (!is_sym(arg(sp, 1), SYM_All)) {
        if (arg(sp, 1)->type != EXPR_INTEGER) return false;
        b = arg(sp, 1)->data.integer;
    }
    if (k == 3) {
        if (arg(sp, 2)->type != EXPR_INTEGER || arg(sp, 2)->data.integer == 0)
            return false;
        c = arg(sp, 2)->data.integer;
    }
    if (a < 0) a = N + a + 1;
    if (b < 0) b = N + b + 1;
    if (c > 0 && a == b + 1) { *empty = true; return a >= 1 && a <= N + 1; }
    if (c < 0 && b == a + 1) { *empty = true; return b >= 1 && b <= N + 1; }
    if (a < 1 || a > N || b < 1 || b > N) return false;
    if ((c > 0 && a > b) || (c < 0 && a < b)) return false;
    *start = a - 1; *stop = b - 1; *step = c;
    return true;
}

/* Several parts: a Span or a List of part elements.  data is borrowed.
 * A List gives a List; an Association gives an Association (keys kept). */
static Expr* q_slicen(const Expr* op, Expr* data) {
    Expr* d;
    Expr* out = NULL;
    if (is_missing(data)) return expr_copy(data);
    d = unpacked(data);
    if (!is_list(d) && !is_association(d)) {
        expr_free(d);
        return missing2("PartInvalid", expr_copy((Expr*)op));
    }
    if (head_name(op) == SYM_Span) {
        int64_t s = 0, e = 0, st = 1, i;
        bool empty;
        size_t cnt = 0;
        Expr** items;
        if (!resolve_span(op, nargs(d), &s, &e, &st, &empty)) {
            expr_free(d);
            return missing2("PartAbsent", expr_copy((Expr*)op));
        }
        items = malloc(sizeof(Expr*) * (nargs(d) ? nargs(d) : 1));
        if (!empty) {
            for (i = s; st > 0 ? i <= e : i >= e; i += st)
                items[cnt++] = expr_copy(arg(d, (size_t)i));   /* List elem or Rule */
        }
        out = expr_new_function(expr_new_symbol(is_list(d) ? SYM_List : SYM_Association),
                                items, cnt);
        free(items);
    } else if (is_list(d)) {
        size_t n = nargs(op), i;
        Expr** items = malloc(sizeof(Expr*) * (n ? n : 1));
        for (i = 0; i < n; i++) items[i] = q_slice1(arg(op, i), d);
        out = new_list(items, n);
        free(items);
    } else {
        /* Association: key specs keep (or report Missing for) their key;
         * integers take the i-th entry and drop an out-of-range one. */
        size_t n = nargs(op), i, cnt = 0;
        Expr** rules = malloc(sizeof(Expr*) * (n ? n : 1));
        for (i = 0; i < n; i++) {
            const Expr* p = arg(op, i);
            if (is_index(p)) {
                size_t pos = resolve_index(p, nargs(d));
                if (pos != SIZE_MAX) rules[cnt++] = expr_copy(arg(d, pos));
            } else {
                const Expr* k = key_of(p);
                Expr* v = assoc_lookup_value(d, k);
                rules[cnt++] = mk2(SYM_Rule, expr_copy((Expr*)k),
                                   v ? expr_copy(v)
                                     : missing2("KeyAbsent", expr_copy((Expr*)k)));
            }
        }
        out = assoc_from_rules(rules, cnt);
        for (i = 0; i < cnt; i++) expr_free(rules[i]);
        free(rules);
    }
    expr_free(d);
    return out;
}

/* ------------------------------------------------------- the query walker */

static Expr* q_run(Expr** ops, size_t n, Expr* data);
static Expr* q_apply(Expr* op, Expr* x);
static Expr* q_normal(Expr** ops, size_t n);
static Expr* q_normal_step(Expr* op, Expr** ops, size_t n);

/* Map the operator chain ops[0..n) over the elements of a List or the values
 * of an Association (keys kept).  Anything else is returned unchanged, as
 * Map leaves an atom alone. */
static Expr* q_map(Expr** ops, size_t n, Expr* data) {
    Expr* d;
    Expr* out;
    size_t m, i;
    if (n == 0) return expr_copy(data);
    if (!is_list(data) && !is_association(data) && !is_packed_list(data))
        return expr_copy(data);
    d = unpacked(data);
    m = nargs(d);
    {
        Expr** items = malloc(sizeof(Expr*) * (m ? m : 1));
        if (is_list(d)) {
            for (i = 0; i < m; i++) items[i] = q_run(ops, n, arg(d, i));
            out = new_list(items, m);
        } else {
            for (i = 0; i < m; i++) {
                Expr* r = arg(d, i);
                items[i] = mk2(SYM_Rule, expr_copy(arg(r, 0)), q_run(ops, n, arg(r, 1)));
            }
            /* The keys are already unique, so the node is canonical as built. */
            out = expr_new_function(expr_new_symbol(SYM_Association), items, m);
        }
        free(items);
    }
    expr_free(d);
    return out;
}

/* Run the operator chain ops[0..n) on data (borrowed); owned result. */
static bool data_first_head(const char* h);

/* Mathilda evaluates a few operator forms (SortBy[f], MaximalBy[f], ...) to
 * the pure function Function[H[#, args...]].  Recognise that shape and hand
 * back the operator form H[args...] (owned), so the operator keeps its
 * descending classification; NULL when op is not such a function. */
static Expr* q_decurry(const Expr* op) {
    const Expr* body;
    const char* h;
    size_t n, i;
    Expr** a;
    Expr* out;
    if (head_name(op) != SYM_Function || nargs(op) != 1) return NULL;
    body = arg(op, 0);
    h = head_name(body);
    if (!h || !data_first_head(h) || nargs(body) < 2) return NULL;
    if (head_name(arg(body, 0)) != SYM_Slot || nargs(arg(body, 0)) != 1 ||
        arg(arg(body, 0), 0)->type != EXPR_INTEGER ||
        arg(arg(body, 0), 0)->data.integer != 1) return NULL;
    n = nargs(body) - 1;
    a = malloc(sizeof(Expr*) * n);
    for (i = 0; i < n; i++) a[i] = expr_copy(arg(body, i + 1));
    out = expr_new_function(expr_copy(body->data.function.head), a, n);
    free(a);
    return out;
}

/* Run the operator chain ops[0..n) on data (borrowed); owned result. */
static Expr* q_run(Expr** ops, size_t n, Expr* data) {
    Expr* dc;
    Expr* op;
    Expr* r;
    Expr* out;
    if (n == 0) return expr_copy(data);
    dc = q_decurry(ops[0]);
    op = dc ? dc : ops[0];
    switch (classify(op)) {
    case QK_SKIP:
        out = q_run(ops + 1, n - 1, data);
        break;
    case QK_ALL:
        out = q_map(ops + 1, n - 1, data);
        break;
    case QK_PART1:
        r = q_slice1(op, data);
        if (n == 1) { out = r; break; }
        out = q_run(ops + 1, n - 1, r);
        expr_free(r);
        break;
    case QK_PARTN:
        r = q_slicen(op, data);
        out = q_map(ops + 1, n - 1, r);
        expr_free(r);
        break;
    case QK_DESC:
        if (is_missing(data)) { out = expr_copy(data); break; }
        r = q_apply(op, data);
        out = q_map(ops + 1, n - 1, r);
        expr_free(r);
        break;
    case QK_DESC1:
        if (is_missing(data)) { out = expr_copy(data); break; }
        r = q_apply(op, data);
        out = q_run(ops + 1, n - 1, r);
        expr_free(r);
        break;
    case QK_ASC:
    default:
        r = q_map(ops + 1, n - 1, data);
        out = q_apply(op, r);
        expr_free(r);
        break;
    }
    if (dc) expr_free(dc);
    return out;
}

/* ------------------------------------------------ operator application */

/* Is `h` an operator-form head whose operand goes FIRST (Select[p][x] is
 * Select[x, p])? */
static bool data_first_head(const char* h) {
    return h == Q_Select || h == Q_SortBy || h == Q_ReverseSortBy ||
           h == Q_GroupBy || h == Q_MaximalBy || h == Q_MinimalBy ||
           h == Q_DeleteDuplicatesBy || h == Q_KeySortBy || h == Q_KeySelect ||
           h == Q_TakeLargestBy || h == Q_TakeSmallestBy || h == Q_TakeLargest ||
           h == Q_TakeSmallest || h == Q_SelectFirst || h == Q_Lookup ||
           h == Q_KeyTake || h == Q_KeyDrop || h == Q_CountsBy || h == Q_Merge ||
           h == Q_DeleteCases || h == Q_Cases || h == Q_Count ||
           h == Q_Position || h == Q_MemberQ || h == Q_FreeQ ||
           h == Q_KeyExistsQ || h == Q_KeyMemberQ || h == Q_KeyFreeQ ||
           h == Q_ReplaceAll || h == Q_Replace || h == Q_Append ||
           h == Q_Prepend || h == Q_Extract || h == Q_Delete ||
           h == Q_ReplacePart || h == Q_Insert || h == Q_AllTrue ||
           h == Q_AnyTrue || h == Q_NoneTrue || h == Q_FirstCase ||
           h == Q_FirstPosition;
}

/* Is `h` an operator-form head whose operand goes LAST (Map[f][x] is
 * Map[f, x])? */
static bool data_last_head(const char* h) {
    return h == Q_Map || h == Q_Apply || h == Q_MapIndexed ||
           h == Q_KeyValueMap || h == Q_AssociationMap || h == Q_KeyMap ||
           h == Q_Scan;
}

/* Operators whose first argument is a criterion applied to each element; a
 * key specification there ("c" or Key["c"]) means "the value at that key". */
static bool criterion_head(const char* h) {
    return h == Q_Select || h == Q_SelectFirst || h == Q_SortBy ||
           h == Q_ReverseSortBy || h == Q_GroupBy || h == Q_MaximalBy ||
           h == Q_MinimalBy || h == Q_DeleteDuplicatesBy || h == Q_CountsBy ||
           h == Q_TakeLargestBy || h == Q_TakeSmallestBy;
}

/* The numeric aggregations that skip Missing values inside a query. */
static bool missing_aware(const char* s) {
    return s == Q_Total || s == Q_Mean || s == Q_Max || s == Q_Min ||
           s == Q_Median || s == Q_Variance || s == Q_StandardDeviation;
}

/* The values of x (a List or Association) with every Missing[...] removed, as
 * a List -- or NULL when x holds no Missing (so the caller can keep x as is,
 * packed fast paths and all). */
static Expr* drop_missing_values(Expr* x) {
    size_t n, i, cnt = 0;
    bool any = false;
    Expr** items;
    Expr* out;
    if (!is_list(x) && !is_association(x)) return NULL;
    n = nargs(x);
    for (i = 0; i < n && !any; i++) any = is_missing(container_item(x, i));
    if (!any) return NULL;
    items = malloc(sizeof(Expr*) * (n ? n : 1));
    for (i = 0; i < n; i++) {
        Expr* v = container_item(x, i);
        if (!is_missing(v)) items[cnt++] = expr_copy(v);
    }
    out = new_list(items, cnt);
    free(items);
    return out;
}

/* Transpose a List of Associations into an Association of Lists (Query's
 * GeneralUtilities`AssociationTranspose).  Keys are taken in first-seen
 * order; a record lacking a key contributes Missing["KeyAbsent", k].
 * Returns NULL when x is not a non-empty List of Associations. */
static Expr* assoc_transpose(Expr* x) {
    size_t n, i, j, nk = 0, cap = 8, n0;
    Expr** keys;
    Expr** rules;
    Expr* out;
    if (!is_list(x) || nargs(x) == 0) return NULL;
    n = nargs(x);
    for (i = 0; i < n; i++) if (!is_association(arg(x, i))) return NULL;
    /* Union of keys in first-seen order.  The first record's keys are unique
     * already; a later record's key is new when the first record (O(1) index
     * probe) and the few extra keys seen so far both lack it -- records
     * normally share their keys, so the extra scan is short. */
    n0 = nargs(arg(x, 0));
    cap = n0 + 8;
    keys = malloc(sizeof(Expr*) * cap);
    for (j = 0; j < n0; j++) keys[nk++] = arg(arg(arg(x, 0), j), 0);
    for (i = 1; i < n; i++) {
        Expr* a = arg(x, i);
        for (j = 0; j < nargs(a); j++) {
            Expr* k = arg(arg(a, j), 0);
            bool dup = assoc_lookup_value(arg(x, 0), k) != NULL;
            size_t q;
            for (q = n0; q < nk && !dup; q++) dup = expr_eq(keys[q], k);
            if (dup) continue;
            if (nk == cap) { cap *= 2; keys = realloc(keys, sizeof(Expr*) * cap); }
            keys[nk++] = k;
        }
    }
    rules = malloc(sizeof(Expr*) * (nk ? nk : 1));
    for (j = 0; j < nk; j++) {
        Expr** col = malloc(sizeof(Expr*) * n);
        for (i = 0; i < n; i++) {
            Expr* v = assoc_lookup_value(arg(x, i), keys[j]);
            col[i] = v ? expr_copy(v) : missing2("KeyAbsent", expr_copy(keys[j]));
        }
        rules[j] = mk2(SYM_Rule, expr_copy(keys[j]), new_list(col, n));
        free(col);
    }
    out = expr_new_function(expr_new_symbol(SYM_Association), rules, nk);
    free(rules);
    free(keys);
    return out;
}

/* MapAt-style update inside a query: {"a" -> f} applies f to the part "a". */
static Expr* q_mapat(Expr* x, const Expr* part, Expr* f) {
    Expr* d = unpacked(x);
    Expr* out;
    size_t pos = SIZE_MAX, i, n;
    if (!is_list(d) && !is_association(d)) return d;
    n = nargs(d);
    if (is_index(part)) {
        pos = resolve_index(part, n);
    } else if (is_key_spec(part) && is_association(d)) {
        const Expr* k = key_of(part);
        for (i = 0; i < n; i++)
            if (expr_eq(arg(arg(d, i), 0), k)) { pos = i; break; }
    }
    if (pos == SIZE_MAX) return d;
    {
        Expr** items = malloc(sizeof(Expr*) * n);
        for (i = 0; i < n; i++) {
            if (i != pos) { items[i] = expr_copy(arg(d, i)); continue; }
            if (is_list(d)) {
                items[i] = q_apply(f, arg(d, i));
            } else {
                Expr* r = arg(d, i);
                items[i] = mk2(SYM_Rule, expr_copy(arg(r, 0)), q_apply(f, arg(r, 1)));
            }
        }
        out = expr_new_function(expr_copy(d->data.function.head), items, n);
        free(items);
    }
    expr_free(d);
    return out;
}

/* A criterion given as a key specification becomes the callable Query[k]. */
static Expr* criterion(Expr* c) {
    if (is_key_spec(c)) return mk1(SYM_Query, expr_copy(c));
    return expr_copy(c);
}

/* Apply one operator to x (borrowed); owned, evaluated result. */
static Expr* q_apply(Expr* op, Expr* x) {
    const char* h;
    if (op->type == EXPR_SYMBOL) {
        const char* s = op->data.symbol.name;
        if ((s == Q_Keys || s == Q_Values) && (is_list(x) || is_packed_list(x))) {
            /* Keys/Values of a list of records: one key/value list per record. */
            Expr* d = unpacked(x);
            size_t n = nargs(d), i;
            Expr** items = malloc(sizeof(Expr*) * (n ? n : 1));
            Expr* out;
            for (i = 0; i < n; i++) items[i] = q_apply(op, arg(d, i));
            out = new_list(items, n);
            free(items);
            expr_free(d);
            return out;
        }
        if (s == Q_Transpose) {
            Expr* t = assoc_transpose(x);
            if (t) return t;
        }
        if (missing_aware(s)) {
            Expr* clean = drop_missing_values(x);
            if (clean) return ev(mk1(s, clean));
        }
        if (s == Q_Identity) return expr_copy(x);
        return ev(mk1(s, expr_copy(x)));
    }
    if (is_part_elem(op)) return q_slice1(op, x);
    if (head_name(op) == SYM_Span || is_part_list(op)) return q_slicen(op, x);

    h = head_name(op);
    if (h == SYM_Query) {
        return q_run(op->data.function.args, nargs(op), x);
    }
    if (h == Q_RightComposition || h == SYM_Composition) {
        /* RightComposition[f1, f2, ...] applies f1 first; Composition the last. */
        size_t n = nargs(op), i;
        Expr* cur = expr_copy(x);
        for (i = 0; i < n; i++) {
            Expr* f = arg(op, h == Q_RightComposition ? i : n - 1 - i);
            Expr* next = q_apply(f, cur);
            expr_free(cur);
            cur = next;
        }
        return cur;
    }
    if (h == SYM_List) {
        size_t n = nargs(op), i;
        bool rules = n > 0;
        for (i = 0; i < n && rules; i++) {
            const char* rh = head_name(arg(op, i));
            rules = (rh == SYM_Rule || rh == SYM_RuleDelayed) && nargs(arg(op, i)) == 2;
        }
        if (rules) {
            /* {part -> f, ...}: update those parts in place. */
            Expr* cur = expr_copy(x);
            for (i = 0; i < n; i++) {
                Expr* r = arg(op, i);
                Expr* next = q_mapat(cur, arg(r, 0), arg(r, 1));
                expr_free(cur);
                cur = next;
            }
            return cur;
        } else {
            /* {f, g, ...}: apply each (as a query) to the same data. */
            Expr** items = malloc(sizeof(Expr*) * (n ? n : 1));
            Expr* out;
            for (i = 0; i < n; i++) items[i] = q_run(&op->data.function.args[i], 1, x);
            out = new_list(items, n);
            free(items);
            return out;
        }
    }
    if (h == SYM_Association) {
        /* <|k1 -> q1, ...|>: an association of query results. */
        size_t n = nargs(op), i;
        Expr** rules = malloc(sizeof(Expr*) * (n ? n : 1));
        Expr* out;
        for (i = 0; i < n; i++) {
            Expr* r = arg(op, i);
            rules[i] = mk2(SYM_Rule, expr_copy(arg(r, 0)),
                           q_run(&r->data.function.args[1], 1, x));
        }
        out = expr_new_function(expr_new_symbol(SYM_Association), rules, n);
        free(rules);
        return out;
    }
    if (h && data_first_head(h)) {
        size_t n = nargs(op), i;
        Expr** a = malloc(sizeof(Expr*) * (n + 1));
        Expr* call;
        a[0] = expr_copy(x);
        for (i = 0; i < n; i++)
            a[i + 1] = (i == 0 && criterion_head(h)) ? criterion(arg(op, 0))
                                                     : expr_copy(arg(op, i));
        call = expr_new_function(expr_copy(op->data.function.head), a, n + 1);
        free(a);
        return ev(call);
    }
    if (h == Q_MapAt && nargs(op) == 2) {
        Expr* a[3] = { expr_copy(arg(op, 0)), expr_copy(x), expr_copy(arg(op, 1)) };
        return ev(expr_new_function(expr_copy(op->data.function.head), a, 3));
    }
    if (h && data_last_head(h)) {
        size_t n = nargs(op), i;
        Expr** a = malloc(sizeof(Expr*) * (n + 1));
        Expr* call;
        for (i = 0; i < n; i++) a[i] = expr_copy(arg(op, i));
        a[n] = expr_copy(x);
        call = expr_new_function(expr_copy(op->data.function.head), a, n + 1);
        free(a);
        return ev(call);
    }
    /* Any other function: op[x], left to the evaluator. */
    {
        Expr* a[1] = { expr_copy(x) };
        return ev(expr_new_function(expr_copy(op), a, 1));
    }
}

/* ------------------------------------------------ Normal[Query[...]] */

static Expr* slice_of(Expr* spec) { return mk1(Q_Slice, spec); }

/* RightComposition[f, g], nested RightCompositions flattened, Identity dropped. */
static Expr* rcompose(Expr* f, Expr* g) {
    Expr** items;
    size_t nf, ng, i, k = 0;
    Expr* out;
    if (is_sym(f, Q_Identity)) { expr_free(f); return g; }
    if (is_sym(g, Q_Identity)) { expr_free(g); return f; }
    nf = head_name(f) == Q_RightComposition ? nargs(f) : 1;
    ng = head_name(g) == Q_RightComposition ? nargs(g) : 1;
    items = malloc(sizeof(Expr*) * (nf + ng));
    if (nf == 1 && head_name(f) != Q_RightComposition) items[k++] = expr_copy(f);
    else for (i = 0; i < nf; i++) items[k++] = expr_copy(arg(f, i));
    if (ng == 1 && head_name(g) != Q_RightComposition) items[k++] = expr_copy(g);
    else for (i = 0; i < ng; i++) items[k++] = expr_copy(arg(g, i));
    out = expr_new_function(expr_new_symbol(Q_RightComposition), items, k);
    free(items);
    expr_free(f);
    expr_free(g);
    return out;
}

/* Map[f], folding Map[Slice[specs]] into Slice[All, specs] as Mathematica
 * does. */
static Expr* map_of(Expr* f) {
    if (is_sym(f, Q_Identity)) return f;
    if (head_name(f) == Q_Slice) {
        size_t n = nargs(f), i;
        Expr** a = malloc(sizeof(Expr*) * (n + 1));
        Expr* out;
        a[0] = expr_new_symbol(SYM_All);
        for (i = 0; i < n; i++) a[i + 1] = expr_copy(arg(f, i));
        out = expr_new_function(expr_new_symbol(Q_Slice), a, n + 1);
        free(a);
        expr_free(f);
        return out;
    }
    return mk1(Q_Map, f);
}

/* The function an ascending/descending operator compiles to. */
static Expr* op_form(Expr* op) {
    if (is_sym(op, Q_Transpose)) return expr_new_symbol(Q_AssocTranspose);
    if (is_query(op)) {
        return q_normal(op->data.function.args, nargs(op));
    }
    if (head_name(op) == SYM_List) {
        size_t n = nargs(op), i;
        bool rules = n > 0;
        for (i = 0; i < n && rules; i++) {
            const char* rh = head_name(arg(op, i));
            rules = (rh == SYM_Rule || rh == SYM_RuleDelayed);
        }
        if (rules) {
            Expr* acc = expr_new_symbol(Q_Identity);
            for (i = 0; i < n; i++)
                acc = rcompose(acc, mk2(Q_MapAt, expr_copy(arg(arg(op, i), 1)),
                                        expr_copy(arg(arg(op, i), 0))));
            return acc;
        }
        return mk1(Q_ApplyThrough, expr_copy(op));
    }
    if (head_name(op) == SYM_Association) return mk1(Q_ApplyThrough, expr_copy(op));
    return expr_copy(op);
}

/* The compiled form of Query[ops...] -- what Mathematica's Normal gives:
 * a RightComposition of Slice / Map / operator pieces. */
static Expr* q_normal(Expr** ops, size_t n) {
    Expr* op;
    Expr* dc;
    Expr* out;
    if (n == 0) return expr_new_symbol(Q_Identity);
    dc = q_decurry(ops[0]);
    op = dc ? dc : ops[0];
    out = q_normal_step(op, ops, n);
    if (dc) expr_free(dc);
    return out;
}

/* One step of q_normal with the (de-curried) first operator `op`. */
static Expr* q_normal_step(Expr* op, Expr** ops, size_t n) {
    Expr* rest;
    switch (classify(op)) {
    case QK_SKIP:  return q_normal(ops + 1, n - 1);
    case QK_ALL:   return map_of(q_normal(ops + 1, n - 1));
    default: break;
    }
    rest = q_normal(ops + 1, n - 1);
    switch (classify(op)) {
    case QK_PART1:
        if (head_name(rest) == Q_Slice) {
            /* RightComposition[Slice[a], Slice[b...]] -> Slice[a, b...] */
            size_t m = nargs(rest), i;
            Expr** a = malloc(sizeof(Expr*) * (m + 1));
            Expr* out;
            a[0] = expr_copy(op);
            for (i = 0; i < m; i++) a[i + 1] = expr_copy(arg(rest, i));
            out = expr_new_function(expr_new_symbol(Q_Slice), a, m + 1);
            free(a);
            expr_free(rest);
            return out;
        }
        return rcompose(slice_of(expr_copy(op)), rest);
    case QK_PARTN:
        return rcompose(slice_of(expr_copy(op)), map_of(rest));
    case QK_DESC:
        return rcompose(op_form(op), map_of(rest));
    case QK_DESC1:
        return rcompose(op_form(op), rest);
    case QK_ASC:
    default:
        return rcompose(map_of(rest), op_form(op));
    }
}

/* ------------------------------------------------------ callable heads */

bool query_callable_probe(const Expr* head) {
    const char* h = head_name(head);
    return h && (h == SYM_Query || h == SYM_Dataset || h == Q_Slice);
}

/* A query result stays a Dataset when it is a List or an Association. */
static Expr* rewrap(Expr* r) {
    if (is_list(r) || is_packed_list(r) || is_association(r))
        return mk1(SYM_Dataset, r);
    return r;
}

Expr* query_callable_apply(const Expr* head, Expr** args, size_t n) {
    const char* h = head_name(head);
    if (h == SYM_Query || h == Q_Slice) {
        Expr* data;
        if (n != 1) return NULL;
        data = args[0];
        if (is_dataset(data))
            return rewrap(q_run(head->data.function.args, nargs(head), arg(data, 0)));
        return q_run(head->data.function.args, nargs(head), data);
    }
    if (h == SYM_Dataset && nargs(head) >= 1)
        return rewrap(q_run(args, n, arg(head, 0)));
    return NULL;
}

/* -------------------------------------------------------- Dimensions */

/* Dimensions of data with associations counted like lists of their values
 * (Mathematica: Dimensions[Dataset[{<|a->1,b->2|>, ...}]] is {n, 2}).
 * Writes up to `max` dims into out and returns the count. */
static size_t ds_dims(const Expr* e, int64_t* out, size_t max) {
    size_t n, i, k, best;
    int64_t sub[64];
    if (max == 0) return 0;
    if (is_packed_list(e)) {
        Expr* d = ndarray_to_nested_list(e);
        k = ds_dims(d, out, max);
        expr_free(d);
        return k;
    }
    if (!is_list(e) && !is_association(e)) return 0;
    n = nargs(e);
    out[0] = (int64_t)n;
    if (n == 0 || max == 1) return 1;
    best = ds_dims(container_item(e, 0), out + 1, max - 1 > 64 ? 64 : max - 1);
    for (i = 1; i < n && best > 0; i++) {
        size_t kk = ds_dims(container_item(e, i), sub, best);
        size_t j = 0;
        while (j < kk && j < best && sub[j] == out[1 + j]) j++;
        best = j;
    }
    return 1 + best;
}

/* ------------------------------------------------- Dataset print form */

typedef struct { char* s; size_t len; size_t cap; } SBuf;

static void sb_put(SBuf* b, const char* s, size_t n) {
    if (b->len + n + 1 > b->cap) {
        size_t nc = b->cap ? b->cap * 2 : 256;
        while (nc < b->len + n + 1) nc *= 2;
        b->s = realloc(b->s, nc);
        b->cap = nc;
    }
    memcpy(b->s + b->len, s, n);
    b->len += n;
    b->s[b->len] = '\0';
}
static void sb_puts(SBuf* b, const char* s) { sb_put(b, s, strlen(s)); }

/* Display width of a UTF-8 string (one column per code point). */
static size_t u8_width(const char* s) {
    size_t w = 0;
    for (; *s; s++) if (((unsigned char)*s & 0xC0) != 0x80) w++;
    return w;
}

#define DS_MAX_ROWS 20
#define DS_MAX_COLS 12
#define DS_CELL_W   24

/* One cell's text: strings raw, everything else in the standard print form,
 * clipped to DS_CELL_W columns. */
static char* cell_text(Expr* e) {
    char* s;
    size_t w, i, cut;
    if (!e) return mathilda_strdup("");
    if (e->type == EXPR_STRING) s = mathilda_strdup(e->data.string);
    else s = expr_to_string(e);
    if (!s) return mathilda_strdup("?");
    for (i = 0; s[i]; i++) if (s[i] == '\n' || s[i] == '\t') s[i] = ' ';
    w = u8_width(s);
    if (w > DS_CELL_W) {
        /* Cut at a code-point boundary after DS_CELL_W - 3 columns. */
        size_t cols = 0;
        for (cut = 0; s[cut]; cut++) {
            if (((unsigned char)s[cut] & 0xC0) != 0x80) {
                if (cols == DS_CELL_W - 3) break;
                cols++;
            }
        }
        s[cut] = '\0';
        s = realloc(s, cut + 4);
        strcat(s, "...");
    }
    return s;
}

static void sb_pad(SBuf* b, const char* s, size_t width) {
    size_t w = u8_width(s);
    sb_puts(b, s);
    while (w++ < width) sb_put(b, " ", 1);
}

/* Render a rows x cols grid of cell strings (row 0 may be a header) as an
 * aligned table.  `header` draws a rule under row 0. */
static void render_grid(SBuf* b, char** cells, size_t rows, size_t cols, bool header) {
    size_t* widths = calloc(cols ? cols : 1, sizeof(size_t));
    size_t r, c;
    for (r = 0; r < rows; r++)
        for (c = 0; c < cols; c++) {
            size_t w = u8_width(cells[r * cols + c]);
            if (w > widths[c]) widths[c] = w;
        }
    for (r = 0; r < rows; r++) {
        for (c = 0; c < cols; c++) {
            if (c) sb_puts(b, " | ");
            if (c + 1 == cols) sb_puts(b, cells[r * cols + c]);
            else sb_pad(b, cells[r * cols + c], widths[c]);
        }
        sb_puts(b, "\n");
        if (r == 0 && header) {
            for (c = 0; c < cols; c++) {
                size_t k;
                if (c) sb_puts(b, "-+-");
                for (k = 0; k < widths[c]; k++) sb_put(b, "-", 1);
            }
            sb_puts(b, "\n");
        }
    }
    free(widths);
}

/* Union of the keys of the associations among items[0..n), first-seen order,
 * capped at DS_MAX_COLS.  Borrowed pointers. */
static size_t key_union(Expr** items, size_t n, Expr** keys, size_t max, bool* more) {
    size_t nk = 0, i, j, q;
    *more = false;
    for (i = 0; i < n; i++) {
        Expr* a = items[i];
        if (!is_association(a)) continue;
        for (j = 0; j < nargs(a); j++) {
            Expr* k = arg(arg(a, j), 0);
            bool dup = false;
            for (q = 0; q < nk; q++) if (expr_eq(keys[q], k)) { dup = true; break; }
            if (dup) continue;
            if (nk == max) { *more = true; continue; }
            keys[nk++] = k;
        }
    }
    return nk;
}

static bool all_assoc(Expr** items, size_t n) {
    size_t i;
    if (n == 0) return false;
    for (i = 0; i < n; i++) if (!is_association(items[i])) return false;
    return true;
}
static bool none_container(Expr** items, size_t n) {
    size_t i;
    for (i = 0; i < n; i++)
        if (is_list(items[i]) || is_association(items[i]) || is_packed_list(items[i]))
            return false;
    return true;
}

char* dataset_format(const Expr* ds) {
    Expr* data;
    Expr* d;
    SBuf b = { NULL, 0, 0 };
    size_t n, shown, i, c;
    bool is_a;
    Expr** items;
    if (!is_dataset(ds)) return NULL;
    data = arg(ds, 0);
    if (!is_list(data) && !is_association(data) && !is_packed_list(data)) {
        char* inner = expr_to_string(data);
        sb_puts(&b, "Dataset[");
        sb_puts(&b, inner ? inner : "");
        sb_puts(&b, "]");
        free(inner);
        return b.s;
    }
    d = unpacked(data);
    is_a = is_association(d);
    n = nargs(d);
    shown = n > DS_MAX_ROWS ? DS_MAX_ROWS : n;
    items = malloc(sizeof(Expr*) * (n ? n : 1));
    for (i = 0; i < n; i++) items[i] = container_item(d, i);

    if (n == 0) {
        sb_puts(&b, is_a ? "Dataset[<||>]" : "Dataset[{}]");
        free(items);
        expr_free(d);
        return b.s;
    }
    {
        char hdr[96];
        int64_t dims[2];
        size_t nd = ds_dims(d, dims, 2);
        if (nd >= 2)
            snprintf(hdr, sizeof hdr, "Dataset <%lld x %lld>\n",
                     (long long)dims[0], (long long)dims[1]);
        else
            snprintf(hdr, sizeof hdr, "Dataset <%lld>\n", (long long)dims[0]);
        sb_puts(&b, hdr);
    }

    if (all_assoc(items, n)) {
        /* Records: header of keys, one row per record (row labels for an
         * association of records). */
        Expr* keys[DS_MAX_COLS];
        bool more;
        size_t nk = key_union(items, shown, keys, DS_MAX_COLS, &more);
        size_t lab = is_a ? 1 : 0;
        size_t cols = nk + lab + (more ? 1 : 0);
        size_t rows = shown + 1;
        char** cells = malloc(sizeof(char*) * rows * (cols ? cols : 1));
        if (lab) cells[0] = mathilda_strdup("");
        for (c = 0; c < nk; c++) cells[lab + c] = cell_text(keys[c]);
        if (more) cells[cols - 1] = mathilda_strdup("...");
        for (i = 0; i < shown; i++) {
            Expr* rec = items[i];
            char** row = cells + (i + 1) * cols;
            if (lab) row[0] = cell_text(arg(arg(d, i), 0));
            for (c = 0; c < nk; c++) {
                Expr* v = assoc_lookup_value(rec, keys[c]);
                row[lab + c] = v ? cell_text(v) : mathilda_strdup("-");
            }
            if (more) row[cols - 1] = mathilda_strdup("...");
        }
        render_grid(&b, cells, rows, cols, true);
        for (i = 0; i < rows * cols; i++) free(cells[i]);
        free(cells);
    } else if (none_container(items, n)) {
        /* Scalars: one column, labelled by key for an association. */
        size_t cols = is_a ? 2 : 1;
        char** cells = malloc(sizeof(char*) * shown * cols);
        for (i = 0; i < shown; i++) {
            if (is_a) cells[i * 2] = cell_text(arg(arg(d, i), 0));
            cells[i * cols + (is_a ? 1 : 0)] = cell_text(items[i]);
        }
        render_grid(&b, cells, shown, cols, false);
        for (i = 0; i < shown * cols; i++) free(cells[i]);
        free(cells);
    } else {
        /* Rows of lists (a matrix) or mixed rows: one cell per element. */
        size_t cols = 0, lab = is_a ? 1 : 0;
        char** cells;
        for (i = 0; i < shown; i++) {
            size_t w = is_list(items[i]) ? nargs(items[i]) : 1;
            if (w > cols) cols = w;
        }
        if (cols > DS_MAX_COLS) cols = DS_MAX_COLS;
        cols += lab;
        cells = malloc(sizeof(char*) * shown * cols);
        for (i = 0; i < shown; i++) {
            char** row = cells + i * cols;
            Expr* it = items[i];
            if (lab) row[0] = cell_text(arg(arg(d, i), 0));
            for (c = lab; c < cols; c++) {
                size_t k = c - lab;
                if (is_list(it)) row[c] = k < nargs(it) ? cell_text(arg(it, k)) : mathilda_strdup("");
                else row[c] = k == 0 ? cell_text(it) : mathilda_strdup("");
            }
        }
        render_grid(&b, cells, shown, cols, false);
        for (i = 0; i < shown * cols; i++) free(cells[i]);
        free(cells);
    }
    if (shown < n) {
        char foot[96];
        snprintf(foot, sizeof foot, "rows 1-%lu of %lu\n",
                 (unsigned long)shown, (unsigned long)n);
        sb_puts(&b, foot);
    }
    /* Drop the trailing newline: the REPL adds its own. */
    if (b.len && b.s[b.len - 1] == '\n') b.s[--b.len] = '\0';
    free(items);
    expr_free(d);
    return b.s;
}

/* ------------------------------------------------------------ builtins */

/* Query[ops...] is inert: it only does something when applied. */
static Expr* builtin_query(Expr* res) { (void)res; return NULL; }

/* Dataset[data] is inert; Dataset[Dataset[d]] flattens, and the
 * type/metadata arguments of Mathematica's Dataset[data, type, meta] are
 * dropped (Mathilda keeps only the data). */
static Expr* builtin_dataset(Expr* res) {
    if (res->type != EXPR_FUNCTION || nargs(res) == 0) return NULL;
    if (is_dataset(arg(res, 0)) && nargs(res) == 1) {
        Expr* inner = arg(res, 0);
        res->data.function.args[0] = NULL;
        return inner;
    }
    if (nargs(res) > 1) return mk1(SYM_Dataset, expr_copy(arg(res, 0)));
    return NULL;
}

/* ------------------------------------- Dataset-transparent list builtins */

typedef enum { DSW_NORMAL, DSW_LENGTH, DSW_DIMS, DSW_QUERY, DSW_CALL } DSMode;

typedef struct {
    const char** name;      /* interned head name (points at a Q_* slot) */
    BuiltinFunc orig;       /* the builtin being wrapped (may be NULL)   */
    size_t pos;             /* argument position holding the Dataset     */
    DSMode mode;
} DSWrap;

static DSWrap g_wraps[] = {
    { &Q_Normal,     NULL, 0, DSW_NORMAL },
    { &Q_Length,     NULL, 0, DSW_LENGTH },
    { &Q_Dimensions, NULL, 0, DSW_DIMS   },
    { &Q_Keys,       NULL, 0, DSW_QUERY  },   /* ds[Keys]            */
    { &Q_Values,     NULL, 0, DSW_QUERY  },   /* ds[Values]          */
    { &Q_Part,       NULL, 0, DSW_QUERY  },   /* ds[[specs]] = ds[specs] */
    { &Q_Select,     NULL, 0, DSW_CALL   },
    { &Q_SortBy,     NULL, 0, DSW_CALL   },
    { &Q_First,      NULL, 0, DSW_CALL   },
    { &Q_Last,       NULL, 0, DSW_CALL   },
    { &Q_Take,       NULL, 0, DSW_CALL   },
    { &Q_Reverse,    NULL, 0, DSW_CALL   },
    { &Q_Map,        NULL, 1, DSW_CALL   },
};
#define N_WRAPS (sizeof(g_wraps) / sizeof(g_wraps[0]))

/* The shared body of every wrapper: handle a Dataset (or, for Normal, a
 * Query) argument, otherwise hand the call to the original builtin.  The
 * common case -- no Dataset -- costs one head-pointer compare. */
static Expr* ds_dispatch(DSWrap* w, Expr* res) {
    const char* h = head_name(res);
    if (res->type != EXPR_FUNCTION) return w->orig ? w->orig(res) : NULL;

    /* Normal[Query[...]] -- the compiled operator form. */
    if (w->mode == DSW_NORMAL && nargs(res) == 1 && is_query(arg(res, 0))) {
        Expr* q = arg(res, 0);
        return q_normal(q->data.function.args, nargs(q));
    }
    if (nargs(res) > w->pos && is_dataset(arg(res, w->pos))) {
        Expr* data = arg(arg(res, w->pos), 0);
        switch (w->mode) {
        case DSW_NORMAL:
            if (nargs(res) == 1) return expr_copy(data);
            break;
        case DSW_LENGTH:
            if (nargs(res) == 1) {
                if (is_list(data) || is_association(data))
                    return expr_new_integer((int64_t)nargs(data));
                return ev(mk1(Q_Length, expr_copy(data)));
            }
            break;
        case DSW_DIMS:
            if (nargs(res) == 1) {
                int64_t dims[64];
                size_t nd = ds_dims(data, dims, 64), k;
                Expr** items = malloc(sizeof(Expr*) * (nd ? nd : 1));
                Expr* out;
                for (k = 0; k < nd; k++) items[k] = expr_new_integer(dims[k]);
                out = new_list(items, nd);
                free(items);
                return out;
            }
            break;
        case DSW_QUERY:
            if (w->name == &Q_Part)
                return rewrap(q_run(res->data.function.args + 1, nargs(res) - 1, data));
            if (nargs(res) == 1) {
                Expr* op = expr_new_symbol(h);
                Expr* out = rewrap(q_run(&op, 1, data));
                expr_free(op);
                return out;
            }
            break;
        case DSW_CALL: {
            /* Re-issue the call on the bare data and re-wrap. */
            size_t n = nargs(res), k;
            Expr** a = malloc(sizeof(Expr*) * n);
            Expr* call;
            for (k = 0; k < n; k++)
                a[k] = expr_copy(k == w->pos ? data : arg(res, k));
            call = expr_new_function(expr_copy(res->data.function.head), a, n);
            free(a);
            return rewrap(ev(call));
        }
        }
    }
    return w->orig ? w->orig(res) : NULL;
}

/* One entry point per wrapped head (C has no closures), each bound to its
 * g_wraps slot. */
#define DS_ENTRY(i) static Expr* ds_entry_##i(Expr* res) { return ds_dispatch(&g_wraps[i], res); }
DS_ENTRY(0) DS_ENTRY(1) DS_ENTRY(2) DS_ENTRY(3) DS_ENTRY(4) DS_ENTRY(5) DS_ENTRY(6)
DS_ENTRY(7) DS_ENTRY(8) DS_ENTRY(9) DS_ENTRY(10) DS_ENTRY(11) DS_ENTRY(12)
#undef DS_ENTRY
static const BuiltinFunc g_entries[] = {
    ds_entry_0, ds_entry_1, ds_entry_2, ds_entry_3, ds_entry_4, ds_entry_5, ds_entry_6,
    ds_entry_7, ds_entry_8, ds_entry_9, ds_entry_10, ds_entry_11, ds_entry_12,
};
typedef char ds_entries_match_wraps[(sizeof g_entries / sizeof g_entries[0]) == N_WRAPS ? 1 : -1];

static void install_dataset_wrappers(void) {
    size_t i;
    for (i = 0; i < N_WRAPS; i++) {
        SymbolDef* def = symtab_get_def(*g_wraps[i].name);
        if (def->builtin_func == g_entries[i]) continue;   /* idempotent */
        g_wraps[i].orig = def->builtin_func;
        def->builtin_func = g_entries[i];
    }
}

/* ---------------------------------------------------------------- init */

void assoc_query_init(void) {
#define I(v, s) v = intern_symbol(s)
    I(Q_Select, "Select"); I(Q_SortBy, "SortBy"); I(Q_ReverseSortBy, "ReverseSortBy");
    I(Q_GroupBy, "GroupBy"); I(Q_MaximalBy, "MaximalBy"); I(Q_MinimalBy, "MinimalBy");
    I(Q_DeleteDuplicatesBy, "DeleteDuplicatesBy"); I(Q_KeySortBy, "KeySortBy");
    I(Q_KeySelect, "KeySelect"); I(Q_KeyMap, "KeyMap");
    I(Q_TakeLargestBy, "TakeLargestBy"); I(Q_TakeSmallestBy, "TakeSmallestBy");
    I(Q_TakeLargest, "TakeLargest"); I(Q_TakeSmallest, "TakeSmallest");
    I(Q_SelectFirst, "SelectFirst"); I(Q_Lookup, "Lookup"); I(Q_KeyTake, "KeyTake");
    I(Q_KeyDrop, "KeyDrop"); I(Q_CountsBy, "CountsBy"); I(Q_Merge, "Merge");
    I(Q_DeleteCases, "DeleteCases"); I(Q_Cases, "Cases"); I(Q_Count, "Count");
    I(Q_Position, "Position"); I(Q_MemberQ, "MemberQ"); I(Q_FreeQ, "FreeQ");
    I(Q_KeyExistsQ, "KeyExistsQ"); I(Q_KeyMemberQ, "KeyMemberQ"); I(Q_KeyFreeQ, "KeyFreeQ");
    I(Q_ReplaceAll, "ReplaceAll"); I(Q_Replace, "Replace"); I(Q_Append, "Append");
    I(Q_Prepend, "Prepend"); I(Q_Extract, "Extract"); I(Q_Delete, "Delete");
    I(Q_ReplacePart, "ReplacePart"); I(Q_Insert, "Insert"); I(Q_AllTrue, "AllTrue");
    I(Q_AnyTrue, "AnyTrue"); I(Q_NoneTrue, "NoneTrue"); I(Q_FirstCase, "FirstCase");
    I(Q_FirstPosition, "FirstPosition");
    I(Q_Map, "Map"); I(Q_Apply, "Apply"); I(Q_MapIndexed, "MapIndexed");
    I(Q_KeyValueMap, "KeyValueMap"); I(Q_AssociationMap, "AssociationMap");
    I(Q_Scan, "Scan"); I(Q_MapAt, "MapAt"); I(Q_Reverse, "Reverse"); I(Q_Sort, "Sort");
    I(Q_ReverseSort, "ReverseSort"); I(Q_DeleteMissing, "DeleteMissing");
    I(Q_Keys, "Keys"); I(Q_Values, "Values"); I(Q_Transpose, "Transpose");
    I(Q_Total, "Total"); I(Q_Mean, "Mean"); I(Q_Max, "Max"); I(Q_Min, "Min");
    I(Q_Median, "Median"); I(Q_Variance, "Variance");
    I(Q_StandardDeviation, "StandardDeviation");
    I(Q_RightComposition, "RightComposition"); I(Q_Composition, "Composition");
    I(Q_Slice, "GeneralUtilities`Slice");
    I(Q_ApplyThrough, "GeneralUtilities`ApplyThrough");
    I(Q_AssocTranspose, "GeneralUtilities`AssociationTranspose");
    I(Q_Identity, "Identity"); I(Q_Normal, "Normal"); I(Q_Length, "Length");
    I(Q_Dimensions, "Dimensions"); I(Q_First, "First"); I(Q_Last, "Last");
    I(Q_Take, "Take"); I(Q_Part, "Part");
#undef I

    symtab_add_builtin("Query", builtin_query);
    symtab_get_def("Query")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Query",
        "Query[op1, op2, ...] is a query that, applied to data, applies op1 at level 1, op2 at "
        "level 2, and so on. Part specifications (\"key\", Key[k], i, i;;j, {p1, p2, ...}, All) "
        "and the filtering/reordering operators Select, SortBy, ReverseSortBy, GroupBy, "
        "MaximalBy, MinimalBy, DeleteDuplicatesBy, KeySortBy, KeySelect, KeyMap, TakeLargestBy, "
        "TakeSmallestBy, SelectFirst, Lookup, Reverse, Sort, ReverseSort, Keys, Values and "
        "DeleteMissing are descending: they act on the data before deeper levels are "
        "processed. Every other function (Total, Length, Max, Mean, Counts, First, f, ...) is "
        "ascending: it acts on the results of the deeper levels. {f, g} applies each operator, "
        "<|k -> f, ...|> builds an association of results, {\"key\" -> f} updates a part, "
        "and f /* g composes. Absent parts give Missing[\"KeyAbsent\", k], "
        "Missing[\"PartAbsent\", i] or Missing[\"PartInvalid\", p], which flow through later "
        "part specifications; Total, Mean, Max, Min, Median, Variance and StandardDeviation "
        "skip Missing values. Query[...][Dataset[d]] gives a Dataset. Normal[Query[...]] gives "
        "the equivalent composition of GeneralUtilities`Slice, Map and operator pieces.");

    symtab_add_builtin("Dataset", builtin_dataset);
    symtab_get_def("Dataset")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Dataset",
        "Dataset[data] represents structured data -- a list, an association, a list of "
        "associations (records) or an association of associations. ds[op1, op2, ...] applies "
        "Query[op1, op2, ...] to the data; a List or Association result is returned as a "
        "Dataset, anything else (a number, a string, Missing[...]) bare. ds[[parts]] is the same "
        "as ds[parts]. Normal[ds] gives the data; Length, Dimensions, Keys, Values, Select, "
        "SortBy, Map, First, Last, Take and Reverse accept a Dataset. Prints as a text table in "
        "the REPL; the type and metadata arguments of Mathematica's Dataset[data, type, meta] "
        "are not kept.");

    install_dataset_wrappers();
}
