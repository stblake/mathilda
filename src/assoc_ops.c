/* ---------------------------------------------------------------------------
 * assoc_ops.c — second-tier Association heads, extended forms, Splice, and
 * Listable threading over association values.
 *
 * See assoc_ops.h for the inventory.  Two design points:
 *
 *   1. Extended forms of an existing head (Keys, Values, KeyUnion, Lookup,
 *      Merge, AssociationMap, PositionIndex, Association, Normal, Transpose,
 *      DeleteMissing) are installed as WRAPPERS: assoc_ops_init reads the
 *      builtin already registered for the head and installs a function that
 *      handles only the new forms and otherwise calls the original.  Every form
 *      the first-tier implementation handled keeps its exact code path (and its
 *      speed), and the first tier can evolve independently.
 *
 *   2. Every bulk operation hashes (ExprSet below, keyed by expr_hash /
 *      expr_eq), so KeyIntersection, KeyComplement, JoinAcross, CountDistinct,
 *      SubsetQ and the KeyUnion fill are all O(total size) rather than
 *      O(n * m).
 *
 * Messages go to stderr through ops_msg, which honours Quiet (the global
 * suppression depth) and notes the firing for Check.
 * -------------------------------------------------------------------------- */

#include "assoc_ops.h"
#include "assoc.h"
#include "attr.h"
#include "common.h"
#include "eval.h"
#include "message.h"
#include "ndarray.h"
#include "print.h"
#include "sym_names.h"
#include "symtab.h"

#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ======================================================================
 * Small utilities.
 * ====================================================================== */

/* Emit a Wolfram-style diagnostic ("Head::tag: text") unless Quiet is active.
 * The firing is always noted so Check[] sees it. */
static void ops_msg(const char* fmt, ...) {
    mth_msg_note_fired();
    if (mth_msg_suppressed()) return;
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

#define ARGC(e)   ((e)->data.function.arg_count)
#define ARG(e, i) ((e)->data.function.args[(i)])
#define RKEY(r)   ((r)->data.function.args[0])
#define RVAL(r)   ((r)->data.function.args[1])

static bool is_fn(const Expr* e) { return e && e->type == EXPR_FUNCTION; }
static bool is_list(const Expr* e) { return head_is(e, SYM_List); }

/* An association whose every entry is a two-argument Rule/RuleDelayed.
 * is_association() only checks the head, so Association[f[1]] can exist; no
 * routine here reads args[0]/args[1] of an entry without this guarantee. */
static bool assoc_ok(const Expr* e) {
    if (!is_association(e)) return false;
    for (size_t i = 0; i < ARGC(e); i++)
        if (!is_rule2(ARG(e, i))) return false;
    return true;
}

/* A List all of whose elements are two-argument rules ({} qualifies). */
static bool is_rule_list(const Expr* e) {
    if (!is_list(e)) return false;
    for (size_t i = 0; i < ARGC(e); i++)
        if (!is_rule2(ARG(e, i))) return false;
    return true;
}

static bool is_true_sym(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_True;
}

static bool is_infinity_expr(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == SYM_Infinity;
    return head_is(e, SYM_DirectedInfinity) && ARGC(e) == 1 &&
           ARG(e, 0)->type == EXPR_INTEGER && ARG(e, 0)->data.integer == 1;
}

static Expr* mk_fn(const char* head, Expr** args, size_t n) {
    return expr_new_function(expr_new_symbol(head), args, n);
}

static Expr* mk_rule(Expr* k, Expr* v) {            /* adopts both */
    Expr* a[2] = { k, v };
    return mk_fn(SYM_Rule, a, 2);
}

static Expr* mk_missing(const char* why) {
    Expr* a[1] = { expr_new_string(why) };
    return mk_fn(SYM_Missing, a, 1);
}

/* f[arg] with `arg` adopted; unevaluated. */
static Expr* mk_call1(const Expr* f, Expr* arg) {
    return expr_new_function(expr_copy((Expr*)f), &arg, 1);
}

/* evaluate(f[arg]); borrows both, returns an owned result. */
static Expr* eval_call1(const Expr* f, const Expr* arg) {
    Expr* call = mk_call1(f, expr_copy((Expr*)arg));
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* A growable array of owned Expr*. */
typedef struct { Expr** v; size_t n, cap; } ExprBuf;

static void eb_push(ExprBuf* b, Expr* e) {
    if (b->n == b->cap) {
        b->cap = b->cap ? 2 * b->cap : 8;
        b->v = realloc(b->v, sizeof(Expr*) * b->cap);
    }
    b->v[b->n++] = e;
}

static void eb_free(ExprBuf* b) {
    for (size_t i = 0; i < b->n; i++) if (b->v[i]) expr_free(b->v[i]);
    free(b->v);
    b->v = NULL; b->n = b->cap = 0;
}

/* Wrap the buffer's contents in head[...] (adopting them) and reset it. */
static Expr* eb_take(ExprBuf* b, const char* head) {
    Expr* r = mk_fn(head, b->v, b->n);
    free(b->v);
    b->v = NULL; b->n = b->cap = 0;
    return r;
}

/* The elements a structural head iterates: an association's values, else the
 * arguments.  Borrowed. */
static Expr* elem_at(const Expr* e, size_t i) {
    Expr* a = ARG(e, i);
    return is_association(e) && is_rule2(a) ? RVAL(a) : a;
}

/* A visible NDArray is an atom to these structural heads, so materialise the
 * first `upto` arguments and hand back the rebuilt call for the evaluator to
 * re-run: the answer must not depend on the representation.  (None of these
 * heads is packed-aware, so the evaluator's gate already materialises a packed
 * List; this covers the visible NDArray the gate never touches.)  NULL when no
 * such argument is present. */
static Expr* ops_delist_visible(Expr* res, size_t upto) {
    size_t n = ARGC(res);
    bool any = false;
    for (size_t i = 0; i < n && i < upto; i++) any = any || is_ndarray(ARG(res, i));
    if (!any) return NULL;
    Expr** args = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++)
        args[i] = (i < upto && is_ndarray(ARG(res, i))) ? ndarray_to_nested_list(ARG(res, i))
                                                        : expr_copy(ARG(res, i));
    Expr* r = expr_new_function(expr_copy(res->data.function.head), args, n);
    free(args);
    return r;
}

/* ======================================================================
 * ExprSet — fixed-capacity open-addressing hash set over BORROWED Expr*
 * (expr_hash / expr_eq).  Sized once for the maximum number of items, load
 * factor < 1/2, so it never rehashes.  Items keep insertion order in `items`.
 * ====================================================================== */
typedef struct {
    size_t* pos;     /* pos[slot] = item index + 1; 0 = empty */
    size_t  mask;
    Expr**  items;
    size_t  n;
} ExprSet;

static void es_init(ExprSet* s, size_t maxn) {
    size_t cap = 8;
    while (cap < 2 * maxn + 1) cap <<= 1;
    s->pos = calloc(cap, sizeof(size_t));
    s->mask = cap - 1;
    s->items = malloc(sizeof(Expr*) * (maxn ? maxn : 1));
    s->n = 0;
}

static void es_free(ExprSet* s) {
    free(s->pos); free(s->items);
    s->pos = NULL; s->items = NULL; s->n = 0;
}

/* Index of `e` in the set, or SIZE_MAX. */
static size_t es_find(const ExprSet* s, const Expr* e) {
    size_t slot = (size_t)expr_hash(e) & s->mask;
    while (s->pos[slot]) {
        size_t i = s->pos[slot] - 1;
        if (expr_eq(s->items[i], (Expr*)e)) return i;
        slot = (slot + 1) & s->mask;
    }
    return SIZE_MAX;
}

/* Index of `e`, inserting it if absent (then *added = true). */
static size_t es_add(ExprSet* s, Expr* e, bool* added) {
    size_t slot = (size_t)expr_hash(e) & s->mask;
    while (s->pos[slot]) {
        size_t i = s->pos[slot] - 1;
        if (expr_eq(s->items[i], e)) { if (added) *added = false; return i; }
        slot = (slot + 1) & s->mask;
    }
    s->items[s->n] = e;
    s->pos[slot] = s->n + 1;
    if (added) *added = true;
    return s->n++;
}

/* Normalise one element of an "associations or rules" collection into an
 * owned, well-formed association: an association is copied, a rule becomes a
 * one-entry association, a list of rules is canonicalised.  NULL otherwise. */
static Expr* as_assoc(const Expr* e) {
    if (assoc_ok(e)) return expr_copy((Expr*)e);
    if (is_rule2(e)) return assoc_from_rules((Expr**)&e, 1);
    if (is_rule_list(e)) return assoc_from_rules(e->data.function.args, ARGC(e));
    return NULL;
}

/* Normalise every element of the List `lst` with as_assoc.  Returns an owned
 * array (*n elements) or NULL if `lst` is not a List or an element is invalid. */
static Expr** as_assoc_array(const Expr* lst, size_t* n) {
    if (!is_list(lst)) return NULL;
    size_t m = ARGC(lst);
    Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
    for (size_t i = 0; i < m; i++) {
        out[i] = as_assoc(ARG(lst, i));
        if (!out[i]) {
            while (i--) expr_free(out[i]);
            free(out);
            return NULL;
        }
    }
    *n = m;
    return out;
}

static void free_array(Expr** a, size_t n) {
    for (size_t i = 0; i < n; i++) if (a[i]) expr_free(a[i]);
    free(a);
}

/* ======================================================================
 * Listable threading over association values (evaluator hook).
 *
 * Rule (Mathematica 15): a Listable head applied to arguments of which at
 * least one is an association, and none is a List, maps over the association
 * values.  All association arguments must carry the same keys in the same
 * order, otherwise Association::incmp and the call stays unevaluated.  Other
 * arguments are repeated for every key:
 *     <|a -> 1, b -> 2|> + 1               -> <|a -> 2, b -> 3|>
 *     <|a -> 1|> + <|a -> 10|>             -> <|a -> 11|>
 *     <|a -> 1, b -> 2|> + <|b -> 1, a -> 2|>   Association::incmp
 * When a List argument is also present, List threading runs first (the List
 * is the outer structure), so this hook only sees List-free calls.  A
 * RuleDelayed entry of the first association stays delayed, as in
 * Mathematica (<|a :> 1 + 1|> + 1 gives <|a :> 1 + (1 + 1)|>).
 * ====================================================================== */
static bool same_keys(const Expr* a, const Expr* b) {
    if (a == b) return true;
    size_t n = ARGC(a);
    if (ARGC(b) != n) return false;
    for (size_t i = 0; i < n; i++)
        if (!expr_eq(RKEY(ARG(a, i)), RKEY(ARG(b, i)))) return false;
    return true;
}

Expr* assoc_thread_listable(Expr* e) {
    size_t argc = ARGC(e);
    const Expr* tmpl = NULL;
    for (size_t i = 0; i < argc; i++) {
        const Expr* a = ARG(e, i);
        if (!is_association(a)) continue;
        if (!assoc_ok(a)) return NULL;
        if (!tmpl) { tmpl = a; continue; }
        if (!same_keys(tmpl, a)) {
            char* s1 = expr_to_string((Expr*)tmpl);
            char* s2 = expr_to_string((Expr*)a);
            char* s3 = expr_to_string(e);
            ops_msg("Association::incmp: The arguments %s and %s in %s are incompatible.",
                    s1, s2, s3);
            free(s1); free(s2); free(s3);
            return NULL;
        }
    }
    if (!tmpl) return NULL;

    size_t n = ARGC(tmpl);
    Expr** rules = malloc(sizeof(Expr*) * (n ? n : 1));
    Expr** fargs = malloc(sizeof(Expr*) * (argc ? argc : 1));
    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < argc; i++) {
            Expr* a = ARG(e, i);
            fargs[i] = is_association(a) ? expr_copy(RVAL(ARG(a, j))) : expr_copy(a);
        }
        Expr* call = expr_new_function(expr_copy(e->data.function.head), fargs, argc);
        rules[j] = assoc_entry_with_value(ARG(tmpl, j), call);
    }
    free(fargs);
    Expr* r = mk_fn(SYM_Association, rules, n);
    free(rules);
    return r;
}

/* ======================================================================
 * Splice[list] / Splice[list, h] (evaluator hook).
 *
 * Splice[{e1, e2, ...}] is replaced by e1, e2, ... when it appears as an
 * argument of a List or an Association; Splice[list, h] splices into any head
 * matching the pattern h (Splice[{1, 2}, _] splices everywhere).  Anywhere
 * else it stays inert.  The check runs in the evaluator's Sequence-flattening
 * step, i.e. only for heads without SequenceHold / HoldAllComplete, and only
 * once a Splice argument has actually been seen, so no other call pays for it.
 * ====================================================================== */
static bool splice_applies(const Expr* e, const Expr* s) {
    size_t n = ARGC(s);
    if (n != 1 && n != 2) return false;
    const Expr* inner = ARG(s, 0);
    if (!is_list(inner) && !is_packed_list(inner)) return false;
    const Expr* h = e->data.function.head;
    if (n == 1)
        return head_is(e, SYM_List) || head_is(e, SYM_Association);
    /* An explicit head spec never splices into an Association (Mathematica
     * 15 leaves <|Splice[{a -> 1}, Association]|> alone). */
    if (head_is(e, SYM_Association)) return false;
    Expr* margs[2] = { expr_copy((Expr*)h), expr_copy(ARG(s, 1)) };
    Expr* mq = mk_fn(SYM_MatchQ, margs, 2);
    Expr* v = evaluate(mq);
    expr_free(mq);
    bool ok = is_true_sym(v);
    expr_free(v);
    return ok;
}

bool eval_splice_args(Expr* e) {
    if (!is_fn(e)) return false;
    size_t argc = ARGC(e);
    bool* hit = calloc(argc ? argc : 1, sizeof(bool));
    size_t newn = 0;
    bool any = false;
    for (size_t i = 0; i < argc; i++) {
        Expr* a = ARG(e, i);
        if (head_is(a, SYM_Splice) && splice_applies(e, a)) {
            Expr* inner = ARG(a, 0);
            if (is_packed_list(inner)) {           /* materialise a buffer */
                Expr* lst = ndarray_to_nested_list(inner);
                expr_free(inner);
                a->data.function.args[0] = lst;
                expr_invalidate_hash(a);
                inner = lst;
            }
            hit[i] = true; any = true;
            newn += ARGC(inner);
        } else {
            newn++;
        }
    }
    if (!any) { free(hit); return false; }

    Expr** na = malloc(sizeof(Expr*) * (newn ? newn : 1));
    size_t k = 0;
    for (size_t i = 0; i < argc; i++) {
        Expr* a = ARG(e, i);
        if (hit[i]) {
            Expr* inner = ARG(a, 0);
            for (size_t j = 0; j < ARGC(inner); j++) na[k++] = expr_copy(ARG(inner, j));
            expr_free(a);
        } else {
            na[k++] = a;
        }
    }
    free(hit);
    free(e->data.function.args);
    e->data.function.args = na;
    e->data.function.arg_count = newn;
    expr_invalidate_hash(e);
    return true;
}

/* ======================================================================
 * KeyIntersection[{assoc1, assoc2, ...}]
 *   Each association restricted to the keys common to all of them, in the
 *   key order of the first.  Elements may also be rules or lists of rules.
 * KeyComplement[{assoc1, assoc2, ...}]
 *   The entries of assoc1 whose keys occur in none of the others.
 * ====================================================================== */
static Expr* builtin_keyintersection(Expr* res) {
    if (ARGC(res) != 1) return NULL;
    size_t m;
    Expr** as = as_assoc_array(ARG(res, 0), &m);
    if (!as) {
        char* s = expr_to_string(ARG(res, 0));
        ops_msg("KeyIntersection::invar: The argument %s is not a valid list of Associations or rules.", s);
        free(s);
        return NULL;
    }
    if (m == 0) { free(as); return mk_fn(SYM_List, NULL, 0); }

    /* Common keys, in the first association's order. */
    ExprBuf common = {0};
    const Expr* a0 = as[0];
    for (size_t i = 0; i < ARGC(a0); i++) {
        Expr* k = RKEY(ARG(a0, i));
        bool everywhere = true;
        for (size_t j = 1; j < m && everywhere; j++)
            everywhere = assoc_lookup_value(as[j], k) != NULL;
        if (everywhere) eb_push(&common, expr_copy(k));
    }

    ExprBuf out = {0};
    for (size_t j = 0; j < m; j++) {
        ExprBuf rules = {0};
        for (size_t c = 0; c < common.n; c++) {
            Expr* v = assoc_lookup_value(as[j], common.v[c]);
            eb_push(&rules, mk_rule(expr_copy(common.v[c]), expr_copy(v)));
        }
        eb_push(&out, eb_take(&rules, SYM_Association));
    }
    eb_free(&common);
    free_array(as, m);
    return eb_take(&out, SYM_List);
}

static Expr* builtin_keycomplement(Expr* res) {
    if (ARGC(res) != 1) return NULL;
    Expr* lst = ARG(res, 0);
    if (is_list(lst) && ARGC(lst) == 0) {
        ops_msg("KeyComplement::empt: Argument {} should be a nonempty list.");
        return NULL;
    }
    size_t m;
    Expr** as = as_assoc_array(lst, &m);
    if (!as) {
        char* s = expr_to_string(lst);
        ops_msg("KeyComplement::invar: The argument %s is not a valid list of Associations or rules.", s);
        free(s);
        return NULL;
    }
    const Expr* a0 = as[0];
    ExprBuf rules = {0};
    for (size_t i = 0; i < ARGC(a0); i++) {
        Expr* ent = ARG(a0, i);
        bool elsewhere = false;
        for (size_t j = 1; j < m && !elsewhere; j++)
            elsewhere = assoc_lookup_value(as[j], RKEY(ent)) != NULL;
        if (!elsewhere) eb_push(&rules, expr_copy(ent));
    }
    free_array(as, m);
    return eb_take(&rules, SYM_Association);
}

/* ======================================================================
 * MissingQ[expr] — True iff expr has head Missing (Missing[], Missing[r, ...]).
 * ====================================================================== */
static Expr* builtin_missingq(Expr* res) {
    if (ARGC(res) != 1) {
        ops_msg("MissingQ::argx: MissingQ called with %zu arguments; 1 argument is expected.",
                ARGC(res));
        return NULL;
    }
    return expr_new_symbol(head_is(ARG(res, 0), SYM_Missing) ? SYM_True : SYM_False);
}

/* ======================================================================
 * Discard[expr, crit] / Discard[expr, crit, n]
 *   The complement of Select: drop the elements for which crit gives True
 *   (at most n of them, first first).  Works on any non-atomic expression and
 *   on associations (tested on the values, keys kept).
 * ====================================================================== */
static Expr* builtin_discard(Expr* res) {
    size_t argc = ARGC(res);
    if (argc != 2 && argc != 3) return NULL;
    { Expr* d = ops_delist_visible(res, 1); if (d) return d; }
    Expr* x = ARG(res, 0);
    Expr* crit = ARG(res, 1);
    if (!is_fn(x)) {
        char* s = expr_to_string(res);
        ops_msg("Discard::normal: Nonatomic expression expected at position 1 in %s.", s);
        free(s);
        return NULL;
    }
    int64_t limit = INT64_MAX;
    if (argc == 3) {
        Expr* n = ARG(res, 2);
        if (n->type == EXPR_INTEGER && n->data.integer >= 0) limit = n->data.integer;
        else if (!is_infinity_expr(n)) {
            char* s = expr_to_string(res);
            ops_msg("Discard::innf: Non-negative integer or Infinity expected at position 3 in %s.", s);
            free(s);
            return NULL;
        }
    }
    size_t n = ARGC(x);
    Expr** keep = malloc(sizeof(Expr*) * (n ? n : 1));
    size_t nk = 0;
    int64_t dropped = 0;
    for (size_t i = 0; i < n; i++) {
        bool drop = false;
        if (dropped < limit) {
            Expr* v = eval_call1(crit, elem_at(x, i));
            drop = is_true_sym(v);
            expr_free(v);
        }
        if (drop) dropped++;
        else keep[nk++] = expr_copy(ARG(x, i));
    }
    Expr* r = expr_new_function(expr_copy(x->data.function.head), keep, nk);
    free(keep);
    return r;
}

/* ======================================================================
 * CountDistinct[expr] / CountDistinctBy[expr, f]
 *   The number of distinct elements (values, for an association), or of
 *   distinct f[element] values.  One hash pass.
 * ====================================================================== */
static Expr* count_distinct(Expr* res, const Expr* f, const char* name) {
    { Expr* d = ops_delist_visible(res, 1); if (d) return d; }
    Expr* x = ARG(res, 0);
    if (!is_fn(x)) {
        char* s = expr_to_string(res);
        ops_msg("%s::normal: Nonatomic expression expected at position 1 in %s.", name, s);
        free(s);
        return NULL;
    }
    size_t n = ARGC(x);
    Expr** owned = f ? malloc(sizeof(Expr*) * (n ? n : 1)) : NULL;
    ExprSet s;
    es_init(&s, n);
    for (size_t i = 0; i < n; i++) {
        Expr* el = elem_at(x, i);
        if (f) { owned[i] = eval_call1(f, el); el = owned[i]; }
        es_add(&s, el, NULL);
    }
    Expr* r = expr_new_integer((int64_t)s.n);
    es_free(&s);
    if (owned) free_array(owned, n);
    return r;
}

static Expr* builtin_countdistinct(Expr* res) {
    if (ARGC(res) != 1) return NULL;
    return count_distinct(res, NULL, "CountDistinct");
}

static Expr* builtin_countdistinctby(Expr* res) {
    if (ARGC(res) != 2) return NULL;
    return count_distinct(res, ARG(res, 1), "CountDistinctBy");
}

/* ======================================================================
 * SubsetQ[a, b] — True iff every element of b occurs in a (multiplicity is
 * ignored).  Lists and associations (by value) mix freely; any other pair of
 * expressions must share a head.
 * ====================================================================== */
static char* head_string(const Expr* e) {
    const char* atom = NULL;
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: atom = "Integer"; break;
        case EXPR_REAL:    atom = "Real"; break;
        case EXPR_STRING:  atom = "String"; break;
        case EXPR_SYMBOL:  atom = "Symbol"; break;
        case EXPR_FUNCTION: return expr_to_string(e->data.function.head);
        default:           atom = "List"; break;     /* packed / NDArray */
    }
    size_t len = strlen(atom) + 1;
    char* s = malloc(len);
    memcpy(s, atom, len);
    return s;
}

static Expr* builtin_subsetq(Expr* res) {
    if (ARGC(res) != 2) return NULL;
    { Expr* d = ops_delist_visible(res, 2); if (d) return d; }
    Expr* a = ARG(res, 0);
    Expr* b = ARG(res, 1);
    bool seqlike_a = is_list(a) || is_association(a);
    bool seqlike_b = is_list(b) || is_association(b);
    bool ok = is_fn(a) && is_fn(b) &&
              ((seqlike_a && seqlike_b) ||
               expr_eq(a->data.function.head, b->data.function.head));
    if (!ok) {
        char* ha = head_string(a);
        char* hb = head_string(b);
        ops_msg("SubsetQ::heads: Heads %s and %s at positions 1 and 2 are expected to be the same.",
                ha, hb);
        free(ha); free(hb);
        return NULL;
    }
    ExprSet s;
    es_init(&s, ARGC(a));
    for (size_t i = 0; i < ARGC(a); i++) es_add(&s, elem_at(a, i), NULL);
    bool subset = true;
    for (size_t i = 0; i < ARGC(b) && subset; i++)
        subset = es_find(&s, elem_at(b, i)) != SIZE_MAX;
    es_free(&s);
    return expr_new_symbol(subset ? SYM_True : SYM_False);
}

/* ======================================================================
 * AssociationComap[{f1, f2, ...}, x] — <|f1 -> f1[x], f2 -> f2[x], ...|>.
 * ====================================================================== */
static Expr* builtin_associationcomap(Expr* res) {
    if (ARGC(res) != 2) return NULL;
    Expr* fs = ARG(res, 0);
    Expr* x = ARG(res, 1);
    if (!is_list(fs)) {
        char* s = expr_to_string(fs);
        ops_msg("AssociationComap::invl: The argument %s is not a list.", s);
        free(s);
        return NULL;
    }
    size_t n = ARGC(fs);
    Expr** rules = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        rules[i] = mk_rule(expr_copy(ARG(fs, i)), mk_call1(ARG(fs, i), expr_copy(x)));
    Expr* r = assoc_from_rules(rules, n);
    free_array(rules, n);
    return r;
}

/* ======================================================================
 * ApplyTo[x, f] (HoldFirst) — x = f[x], returning the new value.  x may be a
 * symbol with a value, a Part of one (s[[i]]), or an association entry
 * (s[key]); the write-back goes through Set, which knows all three shapes.
 * ====================================================================== */
static const char* lvalue_root(const Expr* lhs) {
    while (lhs) {
        if (lhs->type == EXPR_SYMBOL) return lhs->data.symbol.name;
        if (!is_fn(lhs)) return NULL;
        if (head_is(lhs, SYM_Part) && ARGC(lhs) >= 1) { lhs = ARG(lhs, 0); continue; }
        if (lhs->data.function.head->type == EXPR_SYMBOL && ARGC(lhs) >= 1) {
            lhs = lhs->data.function.head;          /* s[key] */
            continue;
        }
        return NULL;
    }
    return NULL;
}

static Expr* builtin_applyto(Expr* res) {
    if (ARGC(res) != 2) {
        ops_msg("ApplyTo::argrx: ApplyTo called with %zu arguments; 2 arguments are expected.",
                ARGC(res));
        return NULL;
    }
    Expr* lhs = ARG(res, 0);
    Expr* f = ARG(res, 1);
    const char* root = lvalue_root(lhs);
    if (!root || symtab_get_own_values(root) == NULL) {
        char* s = expr_to_string(lhs);
        ops_msg("ApplyTo::rvalue: %s is not a variable with a value, so its value cannot be changed.", s);
        free(s);
        return NULL;
    }
    Expr* cur = evaluate(lhs);
    Expr* call = mk_call1(f, cur);                   /* adopts cur */
    Expr* nv = evaluate(call);
    expr_free(call);
    Expr* sargs[2] = { expr_copy(lhs), expr_copy(nv) };
    Expr* set = mk_fn(SYM_Set, sargs, 2);
    Expr* sr = evaluate(set);
    expr_free(set);
    if (sr) expr_free(sr);
    return nv;
}

/* ======================================================================
 * JoinAcross[{a1, ...}, {b1, ...}, spec]            inner join
 * JoinAcross[{a1, ...}, {b1, ...}, spec, type]      type = "Inner" | "Left" |
 *                                                    "Right" | "Outer"
 * spec: a key k, Key[k], k1 -> k2 (different key names on the two sides), or
 * a List of these (join on several keys at once).  Option
 * KeyCollisionFunction -> Left (default) | Right | f, with f[k] giving the
 * pair of keys {kl, kr} under which the two colliding values are kept.
 *
 * Output order (Mathematica 15): matched pairs, left-major (each left row with
 * its matching right rows in right order); then unmatched left rows; then
 * unmatched right rows.  An unmatched row carries the other side's common
 * keys as Missing["Unmatched"]; if the rows then disagree on their key sets,
 * every row is padded to the union (left keys, then right keys) with
 * Missing["NotAvailable"].  The right rows are hash-indexed on the join-key
 * tuple, so the join is O(|left| + |right| + |output|).
 * ====================================================================== */
typedef struct {
    Expr** lk;          /* borrowed left key per join column */
    Expr** rk;          /* borrowed right key per join column */
    size_t nk;
} JoinSpec;

static Expr* unkey(Expr* k) {
    return head_is(k, SYM_Key) && ARGC(k) == 1 ? ARG(k, 0) : k;
}

static void js_add(JoinSpec* js, Expr* item) {
    js->lk = realloc(js->lk, sizeof(Expr*) * (js->nk + 1));
    js->rk = realloc(js->rk, sizeof(Expr*) * (js->nk + 1));
    if (is_rule2(item) && head_is(item, SYM_Rule)) {
        js->lk[js->nk] = unkey(RKEY(item));
        js->rk[js->nk] = unkey(RVAL(item));
    } else {
        js->lk[js->nk] = js->rk[js->nk] = unkey(item);
    }
    js->nk++;
}

/* The join-key tuple of `row` (owned List of the key values), or NULL when a
 * join key is absent (such a row never matches). */
static Expr* join_tuple(const Expr* row, Expr** keys, size_t nk) {
    Expr** vals = malloc(sizeof(Expr*) * nk);
    for (size_t i = 0; i < nk; i++) {
        Expr* v = assoc_lookup_value(row, keys[i]);
        if (!v) { while (i--) expr_free(vals[i]); free(vals); return NULL; }
        vals[i] = expr_copy(v);
    }
    Expr* t = mk_fn(SYM_List, vals, nk);
    free(vals);
    return t;
}

/* True when `k` is a join key spelled the same on both sides. */
static bool js_shared_key(const JoinSpec* js, const Expr* k) {
    for (size_t i = 0; i < js->nk; i++)
        if (expr_eq(js->lk[i], js->rk[i]) && expr_eq(js->lk[i], (Expr*)k)) return true;
    return false;
}

/* Position of key `k` among the entries of `row` (a rule buffer), or SIZE_MAX. */
static size_t row_find(const ExprBuf* row, const Expr* k) {
    for (size_t i = 0; i < row->n; i++)
        if (expr_eq(RKEY(row->v[i]), (Expr*)k)) return i;
    return SIZE_MAX;
}

static void row_copy_entries(ExprBuf* row, const Expr* a) {
    for (size_t i = 0; i < ARGC(a); i++) eb_push(row, expr_copy(ARG(a, i)));
}

/* Keys of rows[0] present in every row (borrowed), and the first-appearance
 * union of the rows' keys (borrowed). */
static void side_keys(Expr** rows, size_t n, ExprBuf* common, ExprBuf* uni) {
    if (n == 0) return;
    for (size_t i = 0; i < ARGC(rows[0]); i++) {
        Expr* k = RKEY(ARG(rows[0], i));
        bool all = true;
        for (size_t j = 1; j < n && all; j++) all = assoc_lookup_value(rows[j], k) != NULL;
        if (all) eb_push(common, k);
    }
    size_t total = 0;
    for (size_t j = 0; j < n; j++) total += ARGC(rows[j]);
    ExprSet s;
    es_init(&s, total);
    for (size_t j = 0; j < n; j++)
        for (size_t i = 0; i < ARGC(rows[j]); i++) {
            bool added;
            es_add(&s, RKEY(ARG(rows[j], i)), &added);
            if (added) eb_push(uni, RKEY(ARG(rows[j], i)));
        }
    es_free(&s);
}

static bool valid_row_list(const Expr* l) {
    if (!is_list(l)) return false;
    for (size_t i = 0; i < ARGC(l); i++) if (!assoc_ok(ARG(l, i))) return false;
    return true;
}

enum { KCF_LEFT, KCF_RIGHT, KCF_FUNC };

/* Merge the entries of right row `r` into `row` (seeded with the left row). */
static bool join_merge(ExprBuf* row, const Expr* r, const JoinSpec* js, int kcf_mode,
                       const Expr* kcf) {
    for (size_t i = 0; i < ARGC(r); i++) {
        Expr* ent = ARG(r, i);
        Expr* k = RKEY(ent);
        if (js_shared_key(js, k)) continue;
        size_t p = row_find(row, k);
        if (p == SIZE_MAX) { eb_push(row, expr_copy(ent)); continue; }
        if (kcf_mode == KCF_LEFT) continue;
        if (kcf_mode == KCF_RIGHT) {
            expr_free(row->v[p]);
            row->v[p] = expr_copy(ent);
            continue;
        }
        Expr* pair = eval_call1(kcf, k);
        if (!is_list(pair) || ARGC(pair) != 2) {
            char* s = expr_to_string((Expr*)kcf);
            ops_msg("JoinAcross::jfun: The value %s of the option KeyCollisionFunction does not evaluate to a list of length 2.", s);
            free(s);
            expr_free(pair);
            return false;
        }
        Expr* lv = expr_copy(RVAL(row->v[p]));
        expr_free(row->v[p]);
        row->v[p] = mk_rule(expr_copy(ARG(pair, 0)), lv);
        eb_push(row, mk_rule(expr_copy(ARG(pair, 1)), expr_copy(RVAL(ent))));
        expr_free(pair);
    }
    return true;
}

static Expr* builtin_joinacross(Expr* res) {
    size_t argc = ARGC(res);
    /* Trailing KeyCollisionFunction -> f options. */
    int kcf_mode = KCF_LEFT;
    const Expr* kcf = NULL;
    while (argc > 3 && is_rule2(ARG(res, argc - 1))) {
        Expr* opt = ARG(res, argc - 1);
        if (!(RKEY(opt)->type == EXPR_SYMBOL &&
              RKEY(opt)->data.symbol.name == SYM_KeyCollisionFunction)) return NULL;
        Expr* v = RVAL(opt);
        if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_Left) kcf_mode = KCF_LEFT;
        else if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_Right) kcf_mode = KCF_RIGHT;
        else { kcf_mode = KCF_FUNC; kcf = v; }
        argc--;
    }
    if (argc != 3 && argc != 4) return NULL;
    Expr* L = ARG(res, 0);
    Expr* R = ARG(res, 1);
    for (int side = 0; side < 2; side++) {
        Expr* l = side ? R : L;
        if (!valid_row_list(l)) {
            char* s = expr_to_string(l);
            ops_msg("JoinAcross::invlc: The argument %s is not a list of associations.", s);
            free(s);
            return NULL;
        }
    }
    bool keep_left = false, keep_right = false;
    if (argc == 4) {
        Expr* t = ARG(res, 3);
        const char* ts = t->type == EXPR_STRING ? t->data.string : NULL;
        if (ts && strcmp(ts, "Inner") == 0) { /* default */ }
        else if (ts && strcmp(ts, "Left") == 0)  keep_left = true;
        else if (ts && strcmp(ts, "Right") == 0) keep_right = true;
        else if (ts && strcmp(ts, "Outer") == 0) keep_left = keep_right = true;
        else {
            char* s = ts ? NULL : expr_to_string(t);
            ops_msg("JoinAcross::jbspc: The join method %s is not \"Inner\", \"Outer\", \"Left\" or \"Right\".",
                    ts ? ts : s);
            free(s);
            return NULL;
        }
    }

    JoinSpec js = { NULL, NULL, 0 };
    Expr* spec = ARG(res, 2);
    if (is_list(spec)) {
        for (size_t i = 0; i < ARGC(spec); i++) js_add(&js, ARG(spec, i));
    } else {
        js_add(&js, spec);
    }
    if (js.nk == 0) { free(js.lk); free(js.rk); return NULL; }

    size_t nl = ARGC(L), nr = ARGC(R);
    Expr** lrows = L->data.function.args;
    Expr** rrows = R->data.function.args;

    /* Hash the right rows by join tuple; chain rows sharing a tuple in order. */
    Expr** rt = malloc(sizeof(Expr*) * (nr ? nr : 1));
    size_t* first = malloc(sizeof(size_t) * (nr ? nr : 1));
    size_t* last = malloc(sizeof(size_t) * (nr ? nr : 1));
    size_t* next = malloc(sizeof(size_t) * (nr ? nr : 1));
    bool* rmatched = calloc(nr ? nr : 1, sizeof(bool));
    ExprSet tix;
    es_init(&tix, nr);
    for (size_t j = 0; j < nr; j++) {
        rt[j] = join_tuple(rrows[j], js.rk, js.nk);
        next[j] = SIZE_MAX;
        if (!rt[j]) continue;
        bool added;
        size_t d = es_add(&tix, rt[j], &added);
        if (added) first[d] = last[d] = j;
        else { next[last[d]] = j; last[d] = j; }
    }

    ExprBuf lcommon = {0}, lunion = {0}, rcommon = {0}, runion = {0};
    side_keys(lrows, nl, &lcommon, &lunion);
    side_keys(rrows, nr, &rcommon, &runion);

    ExprBuf out = {0};
    bool* lmatched = calloc(nl ? nl : 1, sizeof(bool));
    bool failed = false;
    for (size_t i = 0; i < nl && !failed; i++) {
        Expr* t = join_tuple(lrows[i], js.lk, js.nk);
        if (!t) continue;
        size_t d = es_find(&tix, t);
        expr_free(t);
        if (d == SIZE_MAX) continue;
        lmatched[i] = true;
        for (size_t j = first[d]; j != SIZE_MAX; j = next[j]) {
            rmatched[j] = true;
            ExprBuf row = {0};
            row_copy_entries(&row, lrows[i]);
            if (!join_merge(&row, rrows[j], &js, kcf_mode, kcf)) {
                eb_free(&row); failed = true; break;
            }
            Expr* a = assoc_from_rules(row.v, row.n);
            eb_free(&row);
            eb_push(&out, a);
        }
    }
    if (!failed && keep_left) {
        for (size_t i = 0; i < nl; i++) {
            if (lmatched[i]) continue;
            ExprBuf row = {0};
            row_copy_entries(&row, lrows[i]);
            for (size_t c = 0; c < rcommon.n; c++)
                if (row_find(&row, rcommon.v[c]) == SIZE_MAX)
                    eb_push(&row, mk_rule(expr_copy(rcommon.v[c]), mk_missing("Unmatched")));
            eb_push(&out, assoc_from_rules(row.v, row.n));
            eb_free(&row);
        }
    }
    if (!failed && keep_right) {
        for (size_t j = 0; j < nr; j++) {
            if (rmatched[j]) continue;
            ExprBuf row = {0};
            for (size_t c = 0; c < lcommon.n; c++) {
                Expr* k = lcommon.v[c];
                Expr* v = js_shared_key(&js, k) ? assoc_lookup_value(rrows[j], k) : NULL;
                eb_push(&row, mk_rule(expr_copy(k), v ? expr_copy(v) : mk_missing("Unmatched")));
            }
            for (size_t i = 0; i < ARGC(rrows[j]); i++) {
                Expr* ent = ARG(rrows[j], i);
                if (row_find(&row, RKEY(ent)) == SIZE_MAX) eb_push(&row, expr_copy(ent));
            }
            eb_push(&out, assoc_from_rules(row.v, row.n));
            eb_free(&row);
        }
    }

    /* Pad to a common key set when the rows disagree. */
    if (!failed && out.n > 1) {
        bool uniform = true;
        for (size_t r = 1; r < out.n && uniform; r++) {
            if (ARGC(out.v[r]) != ARGC(out.v[0])) { uniform = false; break; }
            for (size_t i = 0; i < ARGC(out.v[0]) && uniform; i++)
                uniform = assoc_lookup_value(out.v[r], RKEY(ARG(out.v[0], i))) != NULL;
        }
        if (!uniform) {
            size_t total = lunion.n + runion.n;
            for (size_t r = 0; r < out.n; r++) total += ARGC(out.v[r]);
            ExprSet present, order;
            es_init(&present, total);
            es_init(&order, total);
            for (size_t r = 0; r < out.n; r++)
                for (size_t i = 0; i < ARGC(out.v[r]); i++)
                    es_add(&present, RKEY(ARG(out.v[r], i)), NULL);
            for (size_t i = 0; i < lunion.n; i++)
                if (es_find(&present, lunion.v[i]) != SIZE_MAX) es_add(&order, lunion.v[i], NULL);
            for (size_t i = 0; i < runion.n; i++)
                if (es_find(&present, runion.v[i]) != SIZE_MAX) es_add(&order, runion.v[i], NULL);
            for (size_t i = 0; i < present.n; i++) es_add(&order, present.items[i], NULL);
            /* `order` borrows keys from the old rows, so build every padded row
             * before any old row is released. */
            Expr** padded = malloc(sizeof(Expr*) * out.n);
            Expr** rules = malloc(sizeof(Expr*) * (order.n ? order.n : 1));
            for (size_t r = 0; r < out.n; r++) {
                for (size_t i = 0; i < order.n; i++) {
                    Expr* v = assoc_lookup_value(out.v[r], order.items[i]);
                    rules[i] = mk_rule(expr_copy(order.items[i]),
                                       v ? expr_copy(v) : mk_missing("NotAvailable"));
                }
                padded[r] = mk_fn(SYM_Association, rules, order.n);
            }
            free(rules);
            es_free(&present);
            es_free(&order);
            for (size_t r = 0; r < out.n; r++) { expr_free(out.v[r]); out.v[r] = padded[r]; }
            free(padded);
        }
    }

    free(lmatched);
    for (size_t j = 0; j < nr; j++) if (rt[j]) expr_free(rt[j]);
    free(rt); free(first); free(last); free(next); free(rmatched);
    es_free(&tix);
    free(lcommon.v); free(lunion.v); free(rcommon.v); free(runion.v);   /* borrowed */
    free(js.lk); free(js.rk);
    if (failed) { eb_free(&out); return NULL; }
    return eb_take(&out, SYM_List);
}

/* ======================================================================
 * Wrapped heads.  g_orig_* hold the builtin registered before assoc_ops_init.
 * ====================================================================== */
static BuiltinFunc g_orig_keys, g_orig_values, g_orig_keyunion, g_orig_lookup,
                   g_orig_merge, g_orig_associationmap, g_orig_positionindex,
                   g_orig_association, g_orig_normal, g_orig_transpose,
                   g_orig_deletemissing;

static Expr* call_orig(BuiltinFunc f, Expr* res) { return f ? f(res) : NULL; }

/* Keys[x] / Values[x] over any nesting of Lists of associations, rules and
 * lists of rules; with f, each key/value is wrapped as f[k].  Returns NULL
 * (after an invrl message when `report`) for an invalid element. */
static Expr* kv_extract(const Expr* x, bool keys, const Expr* f, const char* name, bool report) {
    if (is_rule2(x)) {
        Expr* v = expr_copy(keys ? RKEY(x) : RVAL(x));
        return f ? mk_call1(f, v) : v;
    }
    if (assoc_ok(x) || is_list(x)) {
        size_t n = ARGC(x);
        Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t i = 0; i < n; i++) {
            out[i] = is_association(x) ? kv_extract(ARG(x, i), keys, f, name, false)
                                       : kv_extract(ARG(x, i), keys, f, name, true);
            if (!out[i]) { while (i--) expr_free(out[i]); free(out); return NULL; }
        }
        Expr* l = mk_fn(SYM_List, out, n);
        free(out);
        return l;
    }
    if (report) {
        char* s = expr_to_string((Expr*)x);
        ops_msg("%s::invrl: The argument %s is not a valid Association or a list of rules.", name, s);
        free(s);
    }
    return NULL;
}

static Expr* keys_values(Expr* res, bool keys, BuiltinFunc orig) {
    size_t argc = ARGC(res);
    if (argc == 1) {
        Expr* r = call_orig(orig, res);
        if (r) return r;
        if (!is_list(ARG(res, 0))) return NULL;
        return kv_extract(ARG(res, 0), keys, NULL, keys ? "Keys" : "Values", false);
    }
    if (argc == 2)
        return kv_extract(ARG(res, 0), keys, ARG(res, 1), keys ? "Keys" : "Values", false);
    return NULL;
}

static Expr* ops_keys(Expr* res)   { return keys_values(res, true,  g_orig_keys); }
static Expr* ops_values(Expr* res) { return keys_values(res, false, g_orig_values); }

/* KeyUnion[{...}] accepting rules / rule lists, and KeyUnion[{...}, f] filling
 * each absent key k with f[k]. */
static Expr* ops_keyunion(Expr* res) {
    size_t argc = ARGC(res);
    if (argc != 1 && argc != 2) return NULL;
    Expr* lst = ARG(res, 0);
    if (argc == 1) {
        Expr* r = call_orig(g_orig_keyunion, res);
        if (r) return r;
    }
    size_t m;
    Expr** as = as_assoc_array(lst, &m);
    if (!as) return NULL;
    size_t total = 0;
    for (size_t j = 0; j < m; j++) total += ARGC(as[j]);
    ExprSet uk;
    es_init(&uk, total);
    for (size_t j = 0; j < m; j++)
        for (size_t i = 0; i < ARGC(as[j]); i++) es_add(&uk, RKEY(ARG(as[j], i)), NULL);
    ExprBuf out = {0};
    for (size_t j = 0; j < m; j++) {
        Expr** rules = malloc(sizeof(Expr*) * (uk.n ? uk.n : 1));
        for (size_t u = 0; u < uk.n; u++) {
            Expr* k = uk.items[u];
            Expr* v = assoc_lookup_value(as[j], k);
            Expr* val;
            if (v) val = expr_copy(v);
            else if (argc == 2) val = mk_call1(ARG(res, 1), expr_copy(k));
            else {
                Expr* ma[2] = { expr_new_string("KeyAbsent"), expr_copy(k) };
                val = mk_fn(SYM_Missing, ma, 2);
            }
            rules[u] = mk_rule(expr_copy(k), val);
        }
        eb_push(&out, mk_fn(SYM_Association, rules, uk.n));
        free(rules);
    }
    es_free(&uk);
    free_array(as, m);
    return eb_take(&out, SYM_List);
}

/* Lookup over a List mixing associations and lists of rules threads over the
 * elements: Lookup[{<|a -> 1|>, {a -> 3}, {}}, a, 0] -> {1, 3, 0}. */
static Expr* ops_lookup(Expr* res) {
    size_t argc = ARGC(res);
    if (argc >= 2 && argc <= 3 && is_list(ARG(res, 0)) && ARGC(ARG(res, 0)) > 0) {
        Expr* lst = ARG(res, 0);
        bool any_nonrule = false, any_assoc = false, bad = false;
        for (size_t i = 0; i < ARGC(lst); i++) {
            Expr* el = ARG(lst, i);
            if (is_rule2(el)) continue;
            any_nonrule = true;
            if (is_association(el) || is_list(el)) any_assoc = true;
            else bad = true;
        }
        if (any_nonrule && any_assoc && bad) {
            char* s = expr_to_string(lst);
            ops_msg("Lookup::invrl: The argument %s is not a valid Association or a list of rules.", s);
            free(s);
            return NULL;
        }
        bool mixed = false;           /* at least one rule-list element */
        for (size_t i = 0; i < ARGC(lst) && !mixed; i++) mixed = is_list(ARG(lst, i));
        if (any_nonrule && !bad && mixed) {
            size_t m = ARGC(lst);
            for (size_t i = 0; i < m; i++)
                if (is_rule2(ARG(lst, i))) return call_orig(g_orig_lookup, res);
            Expr** out = malloc(sizeof(Expr*) * m);
            for (size_t i = 0; i < m; i++) {
                Expr** la = malloc(sizeof(Expr*) * argc);
                la[0] = expr_copy(ARG(lst, i));
                for (size_t k = 1; k < argc; k++) la[k] = expr_copy(ARG(res, k));
                Expr* call = expr_new_function(expr_copy(res->data.function.head), la, argc);
                free(la);
                out[i] = evaluate(call);
                expr_free(call);
            }
            Expr* l = mk_fn(SYM_List, out, m);
            free(out);
            return l;
        }
    }
    return call_orig(g_orig_lookup, res);
}

/* Merge over rules and lists of rules: each rule contributes one value, so
 * Merge[{a -> 1, a -> 2}, Total] is <|a -> 3|>. */
static bool merge_flatten(const Expr* el, ExprBuf* out) {
    if (assoc_ok(el)) { eb_push(out, expr_copy((Expr*)el)); return true; }
    if (is_rule2(el)) { eb_push(out, assoc_from_rules((Expr**)&el, 1)); return true; }
    if (is_list(el)) {
        for (size_t i = 0; i < ARGC(el); i++)
            if (!merge_flatten(ARG(el, i), out)) return false;
        return true;
    }
    return false;
}

static Expr* ops_merge(Expr* res) {
    if (ARGC(res) == 2 && is_list(ARG(res, 0))) {
        Expr* lst = ARG(res, 0);
        bool all_assoc = true;
        for (size_t i = 0; i < ARGC(lst) && all_assoc; i++) all_assoc = is_association(ARG(lst, i));
        if (!all_assoc) {
            ExprBuf flat = {0};
            for (size_t i = 0; i < ARGC(lst); i++) {
                if (!merge_flatten(ARG(lst, i), &flat)) {
                    char* s = expr_to_string(ARG(lst, i));
                    ops_msg("Merge::list1: The argument %s is not a valid list of Associations or rules or lists of rules.", s);
                    free(s);
                    eb_free(&flat);
                    return NULL;
                }
            }
            Expr* margs[2] = { eb_take(&flat, SYM_List), expr_copy(ARG(res, 1)) };
            Expr* call = expr_new_function(expr_copy(res->data.function.head), margs, 2);
            Expr* r = call_orig(g_orig_merge, call);
            if (!r) r = expr_copy(call);
            expr_free(call);
            return r;
        }
    }
    return call_orig(g_orig_merge, res);
}

/* AssociationMap[f, assoc]: f is applied to each Rule; each result may be a
 * rule, a list of rules, an association, or {} / Nothing (dropped). */
static Expr* ops_associationmap(Expr* res) {
    if (ARGC(res) == 2 && assoc_ok(ARG(res, 1))) {
        Expr* f = ARG(res, 0);
        Expr* a = ARG(res, 1);
        size_t n = ARGC(a);
        Expr** raw = malloc(sizeof(Expr*) * (n ? n : 1));
        ExprBuf rules = {0};
        bool bad = false;
        for (size_t i = 0; i < n; i++) {
            Expr* ent = ARG(a, i);
            Expr* as_rule = mk_rule(expr_copy(RKEY(ent)), expr_copy(RVAL(ent)));
            raw[i] = eval_call1(f, as_rule);
            Expr* r = raw[i];
            if (is_rule2(r)) eb_push(&rules, expr_copy(r));
            else if (assoc_ok(r) || is_rule_list(r))
                for (size_t j = 0; j < ARGC(r); j++) eb_push(&rules, expr_copy(ARG(r, j)));
            else if (r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_Nothing) { }
            else {
                char* sf = expr_to_string(f);
                char* sr = expr_to_string(as_rule);
                char* sv = expr_to_string(r);
                ops_msg("AssociationMap::invrlf: Applying %s to %s yields %s, which is not a valid rule, list of rules or association.",
                        sf, sr, sv);
                free(sf); free(sr); free(sv);
                bad = true;
            }
            expr_free(as_rule);
        }
        Expr* out;
        if (bad) {
            /* Mathematica leaves the raw applications in an (invalid)
             * Association, which stays unevaluated. */
            out = mk_fn(SYM_Association, raw, n);
            free(raw);
        } else {
            out = assoc_from_rules(rules.v, rules.n);
            free_array(raw, n);
        }
        eb_free(&rules);
        return out;
    }
    return call_orig(g_orig_associationmap, res);
}

/* PositionIndex[assoc] -> <|value -> {keys where it occurs}|>: index the
 * values as a List with the original PositionIndex, then map positions to keys. */
static Expr* ops_positionindex(Expr* res) {
    if (ARGC(res) == 1 && assoc_ok(ARG(res, 0))) {
        Expr* a = ARG(res, 0);
        Expr* vals = assoc_values_list(a);
        Expr* call = mk_fn("PositionIndex", &vals, 1);
        Expr* pi = call_orig(g_orig_positionindex, call);
        expr_free(call);
        if (!pi || !assoc_ok(pi)) { if (pi) expr_free(pi); return NULL; }
        size_t n = ARGC(pi);
        Expr** rules = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t i = 0; i < n; i++) {
            Expr* ent = ARG(pi, i);
            Expr* pos = RVAL(ent);
            size_t m = ARGC(pos);
            Expr** ks = malloc(sizeof(Expr*) * (m ? m : 1));
            for (size_t j = 0; j < m; j++)
                ks[j] = expr_copy(RKEY(ARG(a, (size_t)ARG(pos, j)->data.integer - 1)));
            rules[i] = mk_rule(expr_copy(RKEY(ent)), mk_fn(SYM_List, ks, m));
            free(ks);
        }
        expr_free(pi);
        Expr* r = mk_fn(SYM_Association, rules, n);
        free(rules);
        return r;
    }
    return call_orig(g_orig_positionindex, res);
}

/* Association[...] with arbitrarily nested lists of rules and associations,
 * e.g. Association[{{a -> 1}, {b -> 2}}] -> <|a -> 1, b -> 2|>. */
static bool assoc_flatten(const Expr* a, ExprBuf* rules /* borrowed */) {
    if (is_rule2(a)) { eb_push(rules, (Expr*)a); return true; }
    if (is_list(a) || is_association(a)) {
        for (size_t i = 0; i < ARGC(a); i++)
            if (!assoc_flatten(ARG(a, i), rules)) return false;
        return true;
    }
    return false;
}

static Expr* ops_association(Expr* res) {
    bool nested = false;
    for (size_t i = 0; i < ARGC(res) && !nested; i++) {
        Expr* a = ARG(res, i);
        if (!is_list(a)) continue;
        for (size_t j = 0; j < ARGC(a) && !nested; j++)
            nested = is_list(ARG(a, j)) || is_association(ARG(a, j));
    }
    if (!nested) return call_orig(g_orig_association, res);
    ExprBuf rules = {0};
    bool ok = true;
    for (size_t i = 0; i < ARGC(res) && ok; i++) ok = assoc_flatten(ARG(res, i), &rules);
    Expr* r = ok ? assoc_from_rules(rules.v, rules.n) : NULL;
    free(rules.v);                                      /* borrowed entries */
    return r;
}

/* Normal reaches associations nested anywhere in an expression (outside held
 * heads), converting each to its rule list; like Mathematica it does not
 * descend into the values of a converted association. */
static bool is_held_head(const Expr* e) {
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL) return false;
    return (get_attributes(h->data.symbol.name) & (ATTR_HOLDALL | ATTR_HOLDALLCOMPLETE)) != 0;
}

static bool has_nested_assoc(const Expr* e) {
    if (!is_fn(e)) return false;
    if (assoc_ok(e)) return true;
    if (is_held_head(e)) return false;
    for (size_t i = 0; i < ARGC(e); i++)
        if (has_nested_assoc(ARG(e, i))) return true;
    return false;
}

static Expr* normal_assocs(const Expr* e) {
    if (!has_nested_assoc(e)) return expr_copy((Expr*)e);
    if (assoc_ok(e)) {
        Expr** rs = malloc(sizeof(Expr*) * (ARGC(e) ? ARGC(e) : 1));
        for (size_t i = 0; i < ARGC(e); i++) rs[i] = expr_copy(ARG(e, i));
        Expr* l = mk_fn(SYM_List, rs, ARGC(e));
        free(rs);
        return l;
    }
    Expr** as = malloc(sizeof(Expr*) * (ARGC(e) ? ARGC(e) : 1));
    for (size_t i = 0; i < ARGC(e); i++) as[i] = normal_assocs(ARG(e, i));
    Expr* r = expr_new_function(expr_copy(e->data.function.head), as, ARGC(e));
    free(as);
    return r;
}

static bool spec_has_association(const Expr* s) {
    if (s->type == EXPR_SYMBOL) return s->data.symbol.name == SYM_Association;
    if (is_list(s))
        for (size_t i = 0; i < ARGC(s); i++)
            if (spec_has_association(ARG(s, i))) return true;
    return false;
}

static Expr* ops_normal(Expr* res) {
    size_t argc = ARGC(res);
    if (argc == 2 && spec_has_association(ARG(res, 1)))
        return normal_assocs(ARG(res, 0));
    if (argc == 1 && !is_association(ARG(res, 0)) && has_nested_assoc(ARG(res, 0))) {
        Expr* conv = normal_assocs(ARG(res, 0));
        Expr* call = mk_fn("Normal", &conv, 1);
        Expr* r = call_orig(g_orig_normal, call);
        if (!r) r = expr_copy(ARG(call, 0));
        expr_free(call);
        return r;
    }
    return call_orig(g_orig_normal, res);
}

/* Transpose of an association of associations swaps the two key levels; the
 * inner associations must share their keys in the same order.  A List of
 * associations (no List elements) is a vector, which Transpose returns as is. */
static Expr* ops_transpose(Expr* res) {
    if (ARGC(res) == 1 && is_association(ARG(res, 0))) {
        Expr* a = ARG(res, 0);
        size_t n = ARGC(a);
        bool ok = n > 0 && assoc_ok(a);
        for (size_t i = 0; i < n && ok; i++) {
            Expr* v = RVAL(ARG(a, i));
            ok = assoc_ok(v) && same_keys(RVAL(ARG(a, 0)), v);
        }
        if (!ok) {
            char* s = expr_to_string(a);
            ops_msg("Transpose::nmtx: The first two levels of %s cannot be transposed.", s);
            free(s);
            return NULL;
        }
        Expr* inner0 = RVAL(ARG(a, 0));
        size_t m = ARGC(inner0);
        Expr** outer = malloc(sizeof(Expr*) * (m ? m : 1));
        Expr** col = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t j = 0; j < m; j++) {
            for (size_t i = 0; i < n; i++)
                col[i] = mk_rule(expr_copy(RKEY(ARG(a, i))),
                                 expr_copy(RVAL(ARG(RVAL(ARG(a, i)), j))));
            outer[j] = mk_rule(expr_copy(RKEY(ARG(inner0, j))), mk_fn(SYM_Association, col, n));
        }
        free(col);
        Expr* r = mk_fn(SYM_Association, outer, m);
        free(outer);
        return r;
    }
    if (ARGC(res) == 1 && is_list(ARG(res, 0))) {
        Expr* l = ARG(res, 0);
        bool any_assoc = false, any_list = false;
        for (size_t i = 0; i < ARGC(l); i++) {
            if (is_association(ARG(l, i))) any_assoc = true;
            if (is_list(ARG(l, i)) || is_packed_list(ARG(l, i))) any_list = true;
        }
        if (any_assoc && !any_list) return expr_copy(l);
    }
    return call_orig(g_orig_transpose, res);
}

/* DeleteMissing[expr, n] / DeleteMissing[expr, n, d]: delete, at levels 1..n
 * (n a positive integer or Infinity), every element that has a Missing[...]
 * at depth <= d inside it (d = 0: the element itself is Missing).  Levels are
 * processed deepest first, so an emptied sublist is kept.  Associations are
 * transparent: their values sit one level below them. */
static bool has_missing_within(const Expr* e, int64_t d) {
    if (head_is(e, SYM_Missing)) return true;
    if (d <= 0 || !is_fn(e)) return false;
    for (size_t i = 0; i < ARGC(e); i++)
        if (has_missing_within(elem_at(e, i), d - 1)) return true;
    return false;
}

static Expr* dm_rec(const Expr* e, int64_t level, int64_t maxlev, int64_t d) {
    bool assoc = assoc_ok(e);
    if (!assoc && !is_list(e)) return expr_copy((Expr*)e);
    if (level >= maxlev) return expr_copy((Expr*)e);
    size_t n = ARGC(e);
    Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
    size_t k = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* child = elem_at(e, i);
        Expr* nc = dm_rec(child, level + 1, maxlev, d);
        if (has_missing_within(nc, d)) { expr_free(nc); continue; }
        out[k++] = assoc ? assoc_entry_with_value(ARG(e, i), nc) : nc;
    }
    Expr* r = expr_new_function(expr_copy(e->data.function.head), out, k);
    free(out);
    return r;
}

static Expr* ops_deletemissing(Expr* res) {
    size_t argc = ARGC(res);
    Expr* x = argc ? ARG(res, 0) : NULL;
    if (argc >= 1 && argc <= 3 && !is_list(x) && !is_association(x) && !is_packed_list(x)) {
        char* s = expr_to_string(x);
        ops_msg("DeleteMissing::invrp: The argument %s is not a valid Association or a list.", s);
        free(s);
        return NULL;
    }
    if (argc != 2 && argc != 3) return call_orig(g_orig_deletemissing, res);
    Expr* ls = ARG(res, 1);
    int64_t maxlev, d = 0;
    if (ls->type == EXPR_INTEGER && ls->data.integer > 0) maxlev = ls->data.integer;
    else if (is_infinity_expr(ls)) maxlev = INT64_MAX;
    else {
        char* s = expr_to_string(ls);
        ops_msg("DeleteMissing::arg2: The second argument %s is expected to be a positive integer or Infinity.", s);
        free(s);
        return NULL;
    }
    if (argc == 3) {
        Expr* dd = ARG(res, 2);
        if (dd->type == EXPR_INTEGER && dd->data.integer >= 0) d = dd->data.integer;
        else if (is_infinity_expr(dd)) d = INT64_MAX;
        else return NULL;
    }
    if (is_packed_list(x)) return expr_copy(x);   /* a machine buffer holds no Missing */
    return dm_rec(x, 0, maxlev, d);
}

/* ======================================================================
 * Registration.
 * ====================================================================== */
static BuiltinFunc wrap(const char* name, BuiltinFunc fn) {
    SymbolDef* def = symtab_get_def(name);
    BuiltinFunc orig = def->builtin_func;
    def->builtin_func = fn;
    return orig;
}

void assoc_ops_init(void) {
    /* ---- New heads -------------------------------------------------- */
    symtab_add_builtin("KeyIntersection", builtin_keyintersection);
    symtab_get_def("KeyIntersection")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KeyIntersection",
        "KeyIntersection[{assoc1, assoc2, ...}]\n"
        "\tGives the list of associations restricted to the keys common to all of\n"
        "\tthem, each in the key order of assoc1. Elements may also be rules or\n"
        "\tlists of rules.");

    symtab_add_builtin("KeyComplement", builtin_keycomplement);
    symtab_get_def("KeyComplement")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KeyComplement",
        "KeyComplement[{assoc1, assoc2, ...}]\n"
        "\tGives the entries of assoc1 whose keys occur in none of the other\n"
        "\tassociations. Elements may also be rules or lists of rules.");

    symtab_add_builtin("MissingQ", builtin_missingq);
    symtab_get_def("MissingQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("MissingQ",
        "MissingQ[expr]\n\tGives True if expr has head Missing (Missing[], Missing[\"reason\"],\n"
        "\tMissing[\"KeyAbsent\", k], ...), and False otherwise.");

    symtab_add_builtin("Discard", builtin_discard);
    symtab_get_def("Discard")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Discard",
        "Discard[expr, crit]\n\tDrops the elements e of expr for which crit[e] is True (the\n"
        "\tcomplement of Select). On an association crit tests the values.\n"
        "Discard[expr, crit, n]\n\tDrops only the first n such elements (n may be Infinity).");

    symtab_add_builtin("CountDistinct", builtin_countdistinct);
    symtab_get_def("CountDistinct")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("CountDistinct",
        "CountDistinct[expr]\n\tGives the number of distinct elements of expr (of its values,\n"
        "\tfor an association). One hash pass, O(n).");

    symtab_add_builtin("CountDistinctBy", builtin_countdistinctby);
    symtab_get_def("CountDistinctBy")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("CountDistinctBy",
        "CountDistinctBy[expr, f]\n\tGives the number of distinct values of f[e] over the elements\n"
        "\te of expr (the values, for an association).");

    symtab_add_builtin("SubsetQ", builtin_subsetq);
    symtab_get_def("SubsetQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("SubsetQ",
        "SubsetQ[a, b]\n\tGives True if every element of b occurs in a (multiplicity ignored).\n"
        "\tLists and associations (compared by value) may be mixed; other\n"
        "\texpressions must share a head.");

    symtab_add_builtin("AssociationComap", builtin_associationcomap);
    symtab_get_def("AssociationComap")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("AssociationComap",
        "AssociationComap[{f1, f2, ...}, x]\n\tGives <|f1 -> f1[x], f2 -> f2[x], ...|>.");

    symtab_add_builtin("ApplyTo", builtin_applyto);
    symtab_get_def("ApplyTo")->attributes |= ATTR_HOLDFIRST | ATTR_PROTECTED;
    symtab_set_docstring("ApplyTo",
        "ApplyTo[x, f]\n\tSets x to f[x] and returns the new value. x may be a symbol with a\n"
        "\tvalue, a part s[[i]], or an association entry s[key].");

    symtab_add_builtin("JoinAcross", builtin_joinacross);
    symtab_get_def("JoinAcross")->attributes |= ATTR_PROTECTED;
    {
        Expr* opt = mk_rule(expr_new_symbol(SYM_KeyCollisionFunction), expr_new_symbol(SYM_Left));
        symtab_set_options("JoinAcross", mk_fn(SYM_List, &opt, 1));
    }
    symtab_set_docstring("JoinAcross",
        "JoinAcross[{a1, ...}, {b1, ...}, spec]\n"
        "\tJoins two lists of associations, merging each ai with every bj whose\n"
        "\tjoin keys agree (inner join). spec is a key, Key[k], k1 -> k2 (keys named\n"
        "\tdifferently on the two sides), or a list of these.\n"
        "JoinAcross[{a1, ...}, {b1, ...}, spec, type]\n"
        "\ttype is \"Inner\", \"Left\", \"Right\" or \"Outer\"; unmatched rows get\n"
        "\tMissing[\"Unmatched\"] for the other side's keys.\n"
        "\tOption KeyCollisionFunction -> Left | Right | f resolves a non-join key\n"
        "\tpresent on both sides (f[k] gives the pair of new keys).");

    symtab_get_def("Splice")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Splice",
        "Splice[{e1, e2, ...}]\n\tIs replaced by the sequence e1, e2, ... when it appears inside a\n"
        "\tList or an Association.\n"
        "Splice[{e1, e2, ...}, h]\n\tSplices into any head matching the pattern h.");

    /* Key and Missing are inert heads (Key[k][assoc] and Missing results are
     * handled by the evaluator / Part); they get attributes and docs here. */
    symtab_get_def("Key")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Key",
        "Key[k]\n\tRepresents the key k of an association. Key[k][assoc] gives the value\n"
        "\tat k (Missing[\"KeyAbsent\", k] if absent); assoc[[Key[k]]], Lookup and\n"
        "\tthe key-spec arguments of GroupBy, SortBy and JoinAcross accept it, and\n"
        "\tKey[k] names a key literally even when k is an integer or a list.");

    symtab_get_def("Missing")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Missing",
        "Missing[]  |  Missing[\"reason\"]  |  Missing[\"reason\", data]\n"
        "\tRepresents missing data. Lookup and key access give Missing[\"KeyAbsent\", k]\n"
        "\tfor an absent key; KeyUnion and JoinAcross fill gaps with Missing[...].\n"
        "\tTest with MissingQ, remove with DeleteMissing.");

    /* ---- Wrapped heads: extended forms --------------------------------- */
    g_orig_keys = wrap("Keys", ops_keys);
    symtab_set_docstring("Keys",
        "Keys[assoc]\n\tGives a list of the keys of an association (or of a rule or list of rules).\n"
        "Keys[{assoc1, assoc2, ...}]\n\tThreads over lists of associations and lists of rules.\n"
        "Keys[assoc, f]\n\tWraps each key: {f[k1], f[k2], ...}.");
    g_orig_values = wrap("Values", ops_values);
    symtab_set_docstring("Values",
        "Values[assoc]\n\tGives a list of the values of an association (or of a rule or list of rules).\n"
        "Values[{assoc1, assoc2, ...}]\n\tThreads over lists of associations and lists of rules.\n"
        "Values[assoc, f]\n\tWraps each value: {f[v1], f[v2], ...}.");
    g_orig_keyunion = wrap("KeyUnion", ops_keyunion);
    symtab_set_docstring("KeyUnion",
        "KeyUnion[{assoc1, assoc2, ...}]\n"
        "\tGives the list of associations padded to the union of all their keys;\n"
        "\ta key absent from an association is filled with Missing[\"KeyAbsent\", key].\n"
        "KeyUnion[{assoc1, assoc2, ...}, f]\n\tFills each absent key k with f[k] instead.");
    g_orig_lookup = wrap("Lookup", ops_lookup);
    symtab_set_docstring("Lookup",
        "Lookup[assoc, key]\n\tGives the value for key, or Missing[\"KeyAbsent\", key].\n"
        "Lookup[assoc, key, default]\n\tUses default when key is absent.\n"
        "Lookup[assoc, {k1, k2, ...}]\n\tLooks up several keys at once (O(n+m)).\n"
        "Lookup[{assoc1, assoc2, ...}, key, ...]\n\tThreads over a list of associations and/or\n"
        "\tlists of rules.");
    g_orig_merge = wrap("Merge", ops_merge);
    symtab_set_docstring("Merge",
        "Merge[{assoc1, assoc2, ...}, f]\n"
        "\tCombines associations, applying f to the list of values collected\n"
        "\tfor each key (e.g. Merge[{...}, Total]). The list may also hold rules\n"
        "\tand lists of rules, each rule contributing one value.");
    g_orig_associationmap = wrap("AssociationMap", ops_associationmap);
    symtab_set_docstring("AssociationMap",
        "AssociationMap[f, {k1, k2, ...}]\n\tGives <|k1 -> f[k1], k2 -> f[k2], ...|>.\n"
        "AssociationMap[f, assoc]\n\tApplies f to each rule k -> v of assoc; the results (rules,\n"
        "\tlists of rules or associations) form the new association.");
    g_orig_positionindex = wrap("PositionIndex", ops_positionindex);
    symtab_set_docstring("PositionIndex",
        "PositionIndex[list]\n\tGives <|value -> {positions}|> mapping each distinct\n"
        "\telement to the list of 1-based positions where it occurs. O(n).\n"
        "PositionIndex[assoc]\n\tGives <|value -> {keys}|>, the keys at which each value occurs.");
    g_orig_association = wrap("Association", ops_association);
    g_orig_normal = wrap("Normal", ops_normal);
    g_orig_transpose = wrap("Transpose", ops_transpose);
    g_orig_deletemissing = wrap("DeleteMissing", ops_deletemissing);
    symtab_set_docstring("DeleteMissing",
        "DeleteMissing[expr]\n\tRemoves all Missing[...] elements (equivalent to\n"
        "\tDeleteCases[expr, _Missing]). Over an association, drops entries whose\n"
        "\tvalue is Missing[...].\n"
        "DeleteMissing[expr, n]\n\tRemoves Missing[...] elements at levels 1 through n (n may be\n"
        "\tInfinity); association values count as one level down.\n"
        "DeleteMissing[expr, n, d]\n\tRemoves the elements at levels 1..n that contain a Missing[...]\n"
        "\tat depth d or less.");
}
