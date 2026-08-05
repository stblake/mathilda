/* ---------------------------------------------------------------------------
 * assoc.c — the Association (<| ... |>) data structure and its builtins.
 *
 * See assoc.h for the representation and performance rationale.  In short: an
 * association is Association[Rule[k,v], ...] with unique, insertion-ordered
 * keys, and every bulk operation is driven by a transient open-addressing
 * hash index (KeyIndex) so construction, de-duplication, Merge, KeyDrop,
 * Counts and GroupBy run in amortised O(n) rather than O(n^2).
 * -------------------------------------------------------------------------- */

#include "assoc.h"
#include "assoc_index.h"
#include "ndreduce.h"
#include "ndarray.h"
#include "sym_names.h"
#include "symtab.h"
#include "attr.h"
#include "common.h"
#include "eval.h"

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* ======================================================================
 * KeyIndex — a fixed-capacity open-addressing hash set over Expr* keys.
 *
 * The table stores 1-based indices into a caller-owned array of key
 * pointers (`keys`), so 0 marks an empty slot.  It is sized once, up
 * front, for the maximum number of distinct keys (never grows / rehashes),
 * which keeps every probe branch-predictable and the whole structure a
 * single malloc.  Keys are *borrowed* — the index never owns them.
 * ====================================================================== */
typedef struct {
    size_t* pos;   /* pos[slot] = (entry index into keys[]) + 1; 0 = empty */
    size_t  cap;   /* power of two */
    size_t  mask;  /* cap - 1 */
} KeyIndex;

/* Smallest power of two strictly greater than 2*n (load factor < 0.5),
 * with a floor of 8 so tiny associations still get a real table. */
static size_t ki_capacity_for(size_t n) {
    size_t want = (n < 4 ? 4 : n) * 2 + 1;
    size_t cap = 8;
    while (cap < want) cap <<= 1;
    return cap;
}

static bool ki_init(KeyIndex* ki, size_t n) {
    ki->cap = ki_capacity_for(n);
    ki->mask = ki->cap - 1;
    ki->pos = calloc(ki->cap, sizeof(size_t));
    return ki->pos != NULL;
}

static void ki_free(KeyIndex* ki) {
    free(ki->pos);
    ki->pos = NULL;
}

/* Look up `key`.  Returns the entry index if present, else SIZE_MAX and
 * writes the empty slot to *empty_slot (ready for ki_insert). */
static size_t ki_lookup(const KeyIndex* ki, Expr* const* keys,
                        const Expr* key, size_t* empty_slot) {
    size_t slot = (size_t)expr_hash(key) & ki->mask;
    while (ki->pos[slot] != 0) {
        size_t entry = ki->pos[slot] - 1;
        if (expr_eq(keys[entry], key)) return entry;
        slot = (slot + 1) & ki->mask;
    }
    *empty_slot = slot;
    return SIZE_MAX;
}

static void ki_insert(KeyIndex* ki, size_t slot, size_t entry_index) {
    ki->pos[slot] = entry_index + 1;
}

/* ======================================================================
 * Small expression constructors used throughout this module.
 * ====================================================================== */

/* Build List[...] adopting the `count` expressions in `elems` (the array
 * itself is not freed). */
static Expr* make_list(Expr** elems, size_t count) {
    return expr_new_function(expr_new_symbol(SYM_List), elems, count);
}

/* Build Rule[key, val] adopting both children. */
static Expr* make_rule(Expr* key, Expr* val) {
    Expr* args[2] = { key, val };
    return expr_new_function(expr_new_symbol(SYM_Rule), args, 2);
}

/* Build Missing["KeyAbsent", key] with a fresh copy of `key`. */
static Expr* make_missing(const Expr* key) {
    Expr* args[2] = { expr_new_string("KeyAbsent"), expr_copy((Expr*)key) };
    return expr_new_function(expr_new_symbol(SYM_Missing), args, 2);
}

/* ======================================================================
 * Rule / association predicates and accessors.
 * ====================================================================== */

bool is_association(const Expr* e) {
    return head_is(e, SYM_Association);
}

/* True for a two-argument Rule[k,v] or RuleDelayed[k,v]. */
bool is_rule2(const Expr* e) {
    return (head_is(e, SYM_Rule) || head_is(e, SYM_RuleDelayed)) &&
           e->data.function.arg_count == 2;
}

/* Rebuild an entry around a new value rather than assigning into the old one.
 * See the contract in assoc.h: entries are refcount-shared with the caller's
 * input, so an in-place write reaches every other holder of the association. */
Expr* assoc_entry_with_value(const Expr* entry, Expr* newval) {
    Expr* rargs[2] = { expr_copy(entry->data.function.args[0]), newval };  /* adopted */
    return expr_new_function(expr_copy(entry->data.function.head), rargs, 2);
}

static Expr* rule_key(const Expr* rule) { return rule->data.function.args[0]; }
static Expr* rule_val(const Expr* rule) { return rule->data.function.args[1]; }

/* ======================================================================
 * Canonicalisation: turn a flat list of Rule nodes into a deduplicated,
 * insertion-ordered Association.  First occurrence fixes position; last
 * occurrence fixes value.  Amortised O(count) via the KeyIndex.
 *
 * `rules` are borrowed; the result contains deep copies.
 * ====================================================================== */
Expr* assoc_from_rules(Expr** rules, size_t count) {
    KeyIndex ki;
    if (!ki_init(&ki, count)) return NULL;

    Expr** keys  = malloc(sizeof(Expr*) * (count ? count : 1)); /* borrowed */
    size_t* value_slot = malloc(sizeof(size_t) * (count ? count : 1));
    Expr** out   = malloc(sizeof(Expr*) * (count ? count : 1)); /* owned rules */
    size_t nout = 0;

    for (size_t i = 0; i < count; i++) {
        Expr* k = rule_key(rules[i]);
        Expr* v = rule_val(rules[i]);
        size_t slot;
        size_t idx = ki_lookup(&ki, keys, k, &slot);
        if (idx == SIZE_MAX) {
            keys[nout] = k;
            value_slot[nout] = nout;
            ki_insert(&ki, slot, nout);
            out[nout] = make_rule(expr_copy(k), expr_copy(v));
            nout++;
        } else {
            /* Key already present: overwrite the stored value (last wins). */
            size_t p = value_slot[idx];
            Expr* newrule = make_rule(expr_copy(k), expr_copy(v));
            expr_free(out[p]);
            out[p] = newrule;
        }
    }

    Expr* assoc = expr_new_function(expr_new_symbol(SYM_Association),out, nout);
    free(out);
    free(keys);
    free(value_slot);
    ki_free(&ki);
    return assoc;
}

/* ======================================================================
 * Association[...] — normalise arguments into a canonical association.
 *
 * Accepts, as arguments, any mix of Rule/RuleDelayed nodes, Lists of such
 * rules, and existing Associations (which are spliced).  Returns NULL when
 * the arguments are already a canonical association (all direct 2-arg rules
 * with unique keys) so the evaluator leaves the node untouched — the common
 * <|a->1, b->2|> literal path therefore costs nothing to re-evaluate.
 * ====================================================================== */

/* Append every rule reachable from `arg` to the (dynamic) collector.
 * Returns false if `arg` is not a rule / list-of-rules / association. */
static bool collect_entries(Expr* arg, Expr*** buf, size_t* n, size_t* cap,
                            bool* spliced) {
    if (is_rule2(arg)) {
        if (*n == *cap) { *cap *= 2; *buf = realloc(*buf, sizeof(Expr*) * *cap); }
        (*buf)[(*n)++] = arg;
        return true;
    }
    if (is_association(arg) || head_is(arg, SYM_List)) {
        *spliced = true;
        for (size_t i = 0; i < arg->data.function.arg_count; i++) {
            if (!is_rule2(arg->data.function.args[i])) return false;
            if (*n == *cap) { *cap *= 2; *buf = realloc(*buf, sizeof(Expr*) * *cap); }
            (*buf)[(*n)++] = arg->data.function.args[i];
        }
        return true;
    }
    return false;
}

Expr* builtin_association(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    size_t cap = argc ? argc : 1, n = 0;
    Expr** entries = malloc(sizeof(Expr*) * cap);
    bool spliced = false;

    for (size_t i = 0; i < argc; i++) {
        if (!collect_entries(res->data.function.args[i], &entries, &n, &cap, &spliced)) {
            free(entries);
            return NULL; /* invalid argument: leave unevaluated */
        }
    }

    /* Detect whether the input is already canonical: all direct rules (no
     * splicing) and no duplicate keys.  If so, no rebuild is needed. */
    bool changed = spliced;
    Expr* result = assoc_from_rules(entries, n);
    if (!changed && result && result->data.function.arg_count != n) {
        changed = true; /* de-duplication removed some keys */
    }
    free(entries);

    if (!changed) {
        expr_free(result);
        /* Already canonical: leave `res` as-is. We do NOT attach an index here —
         * evaluate()'s fixed-point step keeps the ORIGINAL node and discards this
         * rebuilt `res`, so an index built here would be thrown away. The index
         * is instead attached lazily by the first single-key read
         * (assoc_lookup_value), on whichever node actually survives. */
        return NULL;
    }
    return result;
}

/* ======================================================================
 * AssociationQ[expr] — True iff expr is an association.
 * ====================================================================== */
Expr* builtin_associationq(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    return expr_new_symbol(is_association(res->data.function.args[0])
                           ? SYM_True : SYM_False);
}

/* ======================================================================
 * Keys[assoc] / Values[assoc] — the keys / values as a List.
 * Also accept a single rule or a List of rules (Wolfram parity).
 * ====================================================================== */
static Expr* keys_or_values(Expr* res, bool want_keys) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* a = res->data.function.args[0];

    if (is_rule2(a))
        return expr_copy(want_keys ? rule_key(a) : rule_val(a));

    if (is_association(a) || head_is(a, SYM_List)) {
        size_t count = a->data.function.arg_count;
        Expr** out = malloc(sizeof(Expr*) * (count ? count : 1));
        for (size_t i = 0; i < count; i++) {
            Expr* el = a->data.function.args[i];
            if (!is_rule2(el)) { while (i--) expr_free(out[i]); free(out); return NULL; }
            out[i] = expr_copy(want_keys ? rule_key(el) : rule_val(el));
        }
        Expr* list = make_list(out, count);
        free(out);
        return list;
    }
    return NULL;
}

Expr* builtin_keys(Expr* res)   { return keys_or_values(res, true);  }
Expr* builtin_values(Expr* res) { return keys_or_values(res, false); }

/* ======================================================================
 * Lookup[assoc, key] / Lookup[assoc, key, default] / Lookup[assoc, {keys}].
 * Missing keys yield Missing["KeyAbsent", key] (or the supplied default).
 * A list of keys builds the index once → O(n + m).
 * ====================================================================== */

/* Linear-scan lookup for a single key (fast path; building an index for one
 * probe is not worth it).  Returns a *borrowed* value pointer or NULL. */
static Expr* assoc_scan(const Expr* assoc, const Expr* key) {
    for (size_t i = 0; i < assoc->data.function.arg_count; i++) {
        Expr* r = assoc->data.function.args[i];
        if (is_rule2(r) && expr_eq(rule_key(r), key)) return rule_val(r);
    }
    return NULL;
}

/* Borrowed value stored under `key`, or NULL if absent.  O(1) amortised when
 * the association carries a persistent index (built at canonicalisation), else
 * the O(n) assoc_scan.  Never mutates `assoc`: no lazy index build on a
 * borrowed/shared node, which keeps single-key reads pure and race-free under
 * the parallel compiled evaluator.  Also accepts a bare List of rules (index
 * absent -> scan), matching the callers that take is_assoc_or_rule_list. */
Expr* assoc_lookup_value(const Expr* assoc, const Expr* key) {
    AssocIndex* idx = assoc->data.function.index;
    if (!idx && is_association(assoc) && assoc->data.function.arg_count > 0) {
        /* Lazily build and cache the index on the FIRST single-key read. Eager
         * attachment at canonicalisation does not survive: evaluate()'s
         * fixed-point step keeps the original association node and frees the
         * rebuilt one an eager index would have been attached to. The reader,
         * by contrast, sees the surviving node — and evaluate()'s in-loop
         * timestamp short-circuit keeps that node stable across a Do/Table/Sum
         * loop, so this O(n) build happens once and every later probe is O(1).
         *
         * Mutating a `const` node is deliberate: the index is benign
         * acceleration metadata on an immutable association (the same discipline
         * as last_evaluated_at), and the interpreter is single-threaded. The
         * compiled evaluator pre-builds the index at its marshalling boundary,
         * so a shared association never reaches this lazy path from a worker
         * thread. assoc_index_build returns NULL for a non-canonical node
         * (an entry that is not a 2-arg rule), leaving the scan below correct. */
        idx = assoc_index_build(assoc->data.function.args, assoc->data.function.arg_count);
        ((Expr*)assoc)->data.function.index = idx;
    }
    if (idx) {
        int64_t e = assoc_index_lookup(idx, assoc->data.function.args, key);
        return e < 0 ? NULL : rule_val(assoc->data.function.args[e]);
    }
    return assoc_scan(assoc, key);
}

/* Build and cache the single-key index NOW.  The compiled evaluator calls this
 * at its marshalling boundary (cf_box), where evaluation is single-threaded, so
 * a subsequent parallel VM run never triggers assoc_lookup_value's lazy index
 * build — which mutates the shared node — from a worker thread.  Idempotent;
 * a no-op if the index already exists, the node is empty, or it is not a
 * canonical association (assoc_index_build returns NULL for a non-2-arg entry). */
void assoc_prebuild_index(const Expr* assoc) {
    if (!assoc || assoc->data.function.index) return;
    if (!is_association(assoc) || assoc->data.function.arg_count == 0) return;
    ((Expr*)assoc)->data.function.index =
        assoc_index_build(assoc->data.function.args, assoc->data.function.arg_count);
}

/* Machine-scalar key lookup for the compiled evaluator (B2 runtime keys).
 * Probes the index (O(1)) with a STACK Expr so a `Lookup[a, i]`-in-a-loop pays
 * no per-iteration malloc.  assoc_lookup_value reads only the key's type and
 * numeric payload (via expr_eq / expr_hash), never its refcount or metadata, so
 * a zeroed stack node is safe.  Returns a borrowed value pointer or NULL. */
Expr* assoc_lookup_value_i64(const Expr* assoc, int64_t key) {
    Expr k; memset(&k, 0, sizeof k);
    k.type = EXPR_INTEGER; k.refcount = 1; k.data.integer = key;
    return assoc_lookup_value(assoc, &k);
}
Expr* assoc_lookup_value_real(const Expr* assoc, double key) {
    Expr k; memset(&k, 0, sizeof k);
    k.type = EXPR_REAL; k.refcount = 1; k.data.real = key;
    return assoc_lookup_value(assoc, &k);
}

/* Lookup / KeyExistsQ accept an association or a bare list of rules (like
 * Keys/Values). A list with non-rule elements simply has no matching keys. */
static bool is_assoc_or_rule_list(const Expr* e) {
    return is_association(e) || head_is(e, SYM_List);
}

Expr* builtin_lookup(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;
    Expr* assoc = res->data.function.args[0];
    Expr* key   = res->data.function.args[1];
    Expr* deflt = argc == 3 ? res->data.function.args[2] : NULL;
    if (!is_assoc_or_rule_list(assoc)) return NULL;

    /* Lookup threads over a (non-empty) list of associations, extracting the
     * key from each: Lookup[{a1, a2, ...}, key] -> {Lookup[a1, key], ...}. This
     * is distinct from the rule-list form (whose elements are Rules, not
     * associations). Each element is delegated to Lookup, so a key-list and a
     * default thread through unchanged. */
    if (head_is(assoc, SYM_List) && assoc->data.function.arg_count > 0) {
        bool all_assoc = true;
        for (size_t i = 0; i < assoc->data.function.arg_count; i++)
            if (!is_association(assoc->data.function.args[i])) { all_assoc = false; break; }
        if (all_assoc) {
            size_t m = assoc->data.function.arg_count;
            Expr** out = malloc(sizeof(Expr*) * m);
            for (size_t i = 0; i < m; i++) {
                Expr** largs = malloc(sizeof(Expr*) * argc);
                largs[0] = expr_copy(assoc->data.function.args[i]);
                largs[1] = expr_copy(key);
                if (argc == 3) largs[2] = expr_copy(deflt);
                Expr* call = expr_new_function(expr_copy(res->data.function.head), largs, argc);
                free(largs);
                out[i] = evaluate(call);
                expr_free(call);
            }
            Expr* list = make_list(out, m);
            free(out);
            return list;
        }
    }

    /* Lookup over a list of keys: single index build, then O(1) per key. */
    if (head_is(key, SYM_List)) {
        size_t na = assoc->data.function.arg_count;
        size_t nk = key->data.function.arg_count;
        KeyIndex ki;
        if (!ki_init(&ki, na)) return NULL;
        Expr** akeys = malloc(sizeof(Expr*) * (na ? na : 1));
        for (size_t i = 0; i < na; i++) {
            Expr* r = assoc->data.function.args[i];
            if (!is_rule2(r)) continue;
            size_t slot, idx = ki_lookup(&ki, akeys, rule_key(r), &slot);
            if (idx == SIZE_MAX) { akeys[i] = rule_key(r); ki_insert(&ki, slot, i); }
        }
        Expr** out = malloc(sizeof(Expr*) * (nk ? nk : 1));
        for (size_t j = 0; j < nk; j++) {
            Expr* qk = key->data.function.args[j];
            size_t slot, idx = ki_lookup(&ki, akeys, qk, &slot);
            if (idx != SIZE_MAX)
                out[j] = expr_copy(rule_val(assoc->data.function.args[idx]));
            else
                out[j] = deflt ? expr_copy(deflt) : make_missing(qk);
        }
        Expr* list = make_list(out, nk);
        free(out); free(akeys); ki_free(&ki);
        return list;
    }

    Expr* v = assoc_lookup_value(assoc, key);
    if (v) return expr_copy(v);
    return deflt ? expr_copy(deflt) : make_missing(key);
}

/* ======================================================================
 * KeyExistsQ[assoc, key] — True iff the key is present.
 * ====================================================================== */
Expr* builtin_keyexistsq(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* assoc = res->data.function.args[0];
    Expr* key   = res->data.function.args[1];
    if (!is_assoc_or_rule_list(assoc)) return NULL;
    return expr_new_symbol(assoc_lookup_value(assoc, key) ? SYM_True : SYM_False);
}

/* KeyMemberQ[assoc, key] == KeyExistsQ; KeyFreeQ[assoc, key] is its complement.
 * All three accept an association or a bare list of rules (like Lookup). */
Expr* builtin_keymemberq(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* assoc = res->data.function.args[0];
    if (!is_assoc_or_rule_list(assoc)) return NULL;
    return expr_new_symbol(assoc_lookup_value(assoc, res->data.function.args[1]) ? SYM_True : SYM_False);
}

Expr* builtin_keyfreeq(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* assoc = res->data.function.args[0];
    if (!is_assoc_or_rule_list(assoc)) return NULL;
    return expr_new_symbol(assoc_lookup_value(assoc, res->data.function.args[1]) ? SYM_False : SYM_True);
}

/* ======================================================================
 * KeyDrop[assoc, key|{keys}] / KeyTake[assoc, key|{keys}].
 * Both preserve association order; a drop/keep set is indexed once.
 * ====================================================================== */

/* Native single-association core: a fresh Association with the keys in `karg`
 * (a key or a List of keys) dropped (take=false) or kept (take=true).  Borrows
 * `assoc` and `karg`; returns an owned node, or NULL if `assoc` is not an
 * association.  The compiled evaluator (B3) calls this directly — no evaluator,
 * no call-node round-trip — and the result carries no index yet (the caller
 * prebuilds one when the result feeds further O(1) reads). */
Expr* assoc_key_select(const Expr* assoc, const Expr* karg, bool take) {
    if (!is_association(assoc)) return NULL;

    /* Normalise the key argument to an array of key pointers. */
    Expr** wanted; size_t nwanted;
    if (head_is(karg, SYM_List)) {
        nwanted = karg->data.function.arg_count;
        wanted = karg->data.function.args;
    } else {
        nwanted = 1;
        wanted = (Expr**)&karg;
    }

    KeyIndex ki;
    if (!ki_init(&ki, nwanted)) return NULL;
    Expr** wk = malloc(sizeof(Expr*) * (nwanted ? nwanted : 1));
    for (size_t i = 0; i < nwanted; i++) {
        size_t slot, idx = ki_lookup(&ki, wk, wanted[i], &slot);
        if (idx == SIZE_MAX) { wk[i] = wanted[i]; ki_insert(&ki, slot, i); }
    }

    size_t na = assoc->data.function.arg_count;
    Expr** out = malloc(sizeof(Expr*) * (na ? na : 1));
    size_t nout = 0;
    for (size_t i = 0; i < na; i++) {
        Expr* r = assoc->data.function.args[i];
        size_t slot;
        bool present = ki_lookup(&ki, wk, rule_key(r), &slot) != SIZE_MAX;
        if (present == take) out[nout++] = expr_copy(r);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_Association), out, nout);
    free(out); free(wk); ki_free(&ki);
    return result;
}

static Expr* key_drop_take(Expr* res, bool take) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* assoc = res->data.function.args[0];
    Expr* karg  = res->data.function.args[1];

    /* Thread over a (non-empty) list of associations, dropping/keeping the keys
     * in each: KeyDrop[{a1, a2, ...}, keys] -> {KeyDrop[a1, keys], ...}. This is
     * the column-of-records form (matching Wolfram and Lookup's threading). Each
     * element is delegated to the same builtin, so a single key or a key-list
     * both thread. */
    if (head_is(assoc, SYM_List) && assoc->data.function.arg_count > 0) {
        bool all_assoc = true;
        for (size_t i = 0; i < assoc->data.function.arg_count; i++)
            if (!is_association(assoc->data.function.args[i])) { all_assoc = false; break; }
        if (all_assoc) {
            size_t m = assoc->data.function.arg_count;
            Expr** out = malloc(sizeof(Expr*) * m);
            for (size_t i = 0; i < m; i++) {
                Expr* a2[2] = { expr_copy(assoc->data.function.args[i]), expr_copy(karg) };
                Expr* call = expr_new_function(expr_copy(res->data.function.head), a2, 2);
                out[i] = evaluate(call);
                expr_free(call);
            }
            Expr* list = make_list(out, m);
            free(out);
            return list;
        }
    }

    return assoc_key_select(assoc, karg, take);
}

Expr* builtin_keydrop(Expr* res) { return key_drop_take(res, false); }
Expr* builtin_keytake(Expr* res) { return key_drop_take(res, true);  }

/* ======================================================================
 * KeyUnion[{assoc1, assoc2, ...}] — pad every association to the shared key
 * set (the union of all keys, in first-appearance order). A key absent from a
 * given association is filled with Missing["KeyAbsent", key]. Returns the list
 * of equalised associations, ready for row-wise / tabular processing.
 *
 * The union is collected in one O(sum of sizes) hash pass; each association is
 * then rebuilt against a hash index of its own keys, so the whole operation is
 * linear rather than the O(keys * entries) of repeated membership scans.
 * ====================================================================== */
Expr* builtin_keyunion(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;   /* list form only */
    Expr* lst = res->data.function.args[0];
    if (!head_is(lst, SYM_List)) return NULL;
    size_t m = lst->data.function.arg_count;
    for (size_t j = 0; j < m; j++)
        if (!is_association(lst->data.function.args[j])) return NULL;

    /* Upper bound on distinct keys: the total number of entries. */
    size_t total = 0;
    for (size_t j = 0; j < m; j++)
        total += lst->data.function.args[j]->data.function.arg_count;

    KeyIndex uki;
    if (!ki_init(&uki, total)) return NULL;
    Expr** ukeys = malloc(sizeof(Expr*) * (total ? total : 1)); /* borrowed */
    size_t nu = 0;
    for (size_t j = 0; j < m; j++) {
        Expr* a = lst->data.function.args[j];
        size_t na = a->data.function.arg_count;
        for (size_t i = 0; i < na; i++) {
            Expr* k = rule_key(a->data.function.args[i]);
            size_t slot;
            if (ki_lookup(&uki, ukeys, k, &slot) == SIZE_MAX) {
                ukeys[nu] = k; ki_insert(&uki, slot, nu); nu++;
            }
        }
    }

    Expr** outer = malloc(sizeof(Expr*) * (m ? m : 1));
    for (size_t j = 0; j < m; j++) {
        Expr* a = lst->data.function.args[j];
        size_t na = a->data.function.arg_count;
        /* Index this association's own keys so each union key is one lookup. */
        KeyIndex jki;
        ki_init(&jki, na);
        Expr** jkeys = malloc(sizeof(Expr*) * (na ? na : 1)); /* borrowed */
        for (size_t i = 0; i < na; i++) {
            Expr* k = rule_key(a->data.function.args[i]);
            size_t slot;
            if (ki_lookup(&jki, jkeys, k, &slot) == SIZE_MAX) {
                jkeys[i] = k; ki_insert(&jki, slot, i);
            }
        }
        Expr** entries = malloc(sizeof(Expr*) * (nu ? nu : 1));
        for (size_t u = 0; u < nu; u++) {
            Expr* k = ukeys[u];
            size_t slot, idx = ki_lookup(&jki, jkeys, k, &slot);
            Expr* val = (idx == SIZE_MAX)
                ? make_missing(k)
                : expr_copy(rule_val(a->data.function.args[idx]));
            Expr* rargs[2] = { expr_copy(k), val };
            entries[u] = expr_new_function(expr_new_symbol(SYM_Rule), rargs, 2);
        }
        outer[j] = expr_new_function(expr_new_symbol(SYM_Association),entries, nu);
        free(entries); free(jkeys); ki_free(&jki);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), outer, m);
    free(outer); free(ukeys); ki_free(&uki);
    return result;
}

/* ======================================================================
 * KeyValueMap[f, assoc] — List[f[k1,v1], f[k2,v2], ...].
 * The applications are left for the evaluator to reduce.
 * ====================================================================== */
Expr* builtin_keyvaluemap(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* f     = res->data.function.args[0];
    Expr* assoc = res->data.function.args[1];
    if (!is_association(assoc)) return NULL;

    size_t n = assoc->data.function.arg_count;
    Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr* r = assoc->data.function.args[i];
        Expr* fargs[2] = { expr_copy(rule_key(r)), expr_copy(rule_val(r)) };
        out[i] = expr_new_function(expr_copy(f), fargs, 2);
    }
    Expr* list = make_list(out, n);
    free(out);
    return list;
}

/* ======================================================================
 * AssociationThread[keys, values] / AssociationThread[keys -> values].
 * ====================================================================== */
Expr* builtin_associationthread(Expr* res) {
    Expr *keys = NULL, *vals = NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 1 && is_rule2(res->data.function.args[0])) {
        keys = rule_key(res->data.function.args[0]);
        vals = rule_val(res->data.function.args[0]);
    } else if (argc == 2) {
        keys = res->data.function.args[0];
        vals = res->data.function.args[1];
    } else {
        return NULL;
    }
    if (!head_is(keys, SYM_List) || !head_is(vals, SYM_List)) return NULL;
    size_t n = keys->data.function.arg_count;
    if (n != vals->data.function.arg_count) return NULL;

    Expr** rules = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        rules[i] = make_rule(expr_copy(keys->data.function.args[i]),
                             expr_copy(vals->data.function.args[i]));
    Expr* assoc = assoc_from_rules(rules, n);
    for (size_t i = 0; i < n; i++) expr_free(rules[i]);
    free(rules);
    return assoc;
}

/* ======================================================================
 * Counts[list] — <|element -> count, ...|> by first appearance.
 * Hash-indexed: one pass, O(n).
 * ====================================================================== */
/* Counts[buffer]: run Tally over the machine words and relabel its {key, count}
 * pairs as key -> count rules. Returns NULL if Tally did not answer with the
 * pair-list shape, in which case the caller carries on as before. */
static Expr* counts_from_ndarray(Expr* list) {
    Expr** targs = malloc(sizeof(Expr*) * 1);
    if (!targs) return NULL;
    targs[0] = expr_copy(list);
    Expr* tcall = expr_new_function(expr_new_symbol("Tally"), targs, 1);
    free(targs);
    Expr* tally = ndred_tally(tcall);
    expr_free(tcall);
    if (!tally) return NULL;

    /* Tally answers in EITHER of two shapes and both have to be read here.
     * Over an int64 buffer its {key, count} pairs are themselves machine
     * words, so it hands back a rank-2 ARRAY -- which is the whole point of
     * the direct-indexed path -- while the float64 route builds an ordinary
     * List of pairs. Reading only the List shape silently left Counts
     * unevaluated for exactly the integer data it exists for. */
    size_t nd;
    if (is_ndarray(tally)) {
        if (tally->data.ndarray.rank != 2 || tally->data.ndarray.dims[1] != 2) {
            expr_free(tally); return NULL;
        }
        nd = (size_t)tally->data.ndarray.dims[0];
    } else if (head_is(tally, SYM_List)) {
        nd = tally->data.function.arg_count;
    } else {
        expr_free(tally); return NULL;
    }

    Expr** rules = malloc(sizeof(Expr*) * (nd ? nd : 1));
    if (!rules) { expr_free(tally); return NULL; }
    for (size_t i = 0; i < nd; i++) {
        Expr *key, *cnt;
        if (is_ndarray(tally)) {
            /* ndarray_element_to_expr is the one place that decides an
             * element's head, so an int64 count comes back an exact Integer. */
            key = ndarray_element_to_expr(tally, 2 * i);
            cnt = ndarray_element_to_expr(tally, 2 * i + 1);
        } else {
            Expr* pair = tally->data.function.args[i];
            if (!head_is(pair, SYM_List) || pair->data.function.arg_count != 2) {
                for (size_t j = 0; j < i; j++) expr_free(rules[j]);
                free(rules); expr_free(tally);
                return NULL;
            }
            key = expr_copy(pair->data.function.args[0]);
            cnt = expr_copy(pair->data.function.args[1]);
        }
        rules[i] = make_rule(key, cnt);
    }
    Expr* assoc = expr_new_function(expr_new_symbol(SYM_Association),rules, nd);
    free(rules);
    expr_free(tally);
    return assoc;
}

/* Counts over a machine array -> <|element -> count|>, native (via ndred_tally,
 * no evaluator).  Borrows `arr`; owned result, or NULL if `arr` is not an
 * NDArray or the tally could not be formed.  Exposed for the compiled evaluator. */
Expr* assoc_counts_ndarray(const Expr* arr) {
    if (!arr || !is_ndarray(arr)) return NULL;
    return counts_from_ndarray((Expr*)arr);
}

Expr* builtin_counts(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* list = res->data.function.args[0];
    /* A buffer is counted by Tally's machine-word hash, then relabelled.
     *
     * Counts and Tally answer the same question in two presentations -- a
     * key -> count Association against a list of {key, count} pairs -- and
     * ndred_tally already does the counting on the raw int64/float64 words,
     * direct-indexed when the value range allows and hashed otherwise. Running
     * it and rewriting the pairs as rules costs one pass over the DISTINCT
     * values, which is the small number here; the alternative was ki_lookup
     * over 10^6 boxed Exprs. Tally declines exactly the dtypes and ranks it
     * cannot key faithfully (complex, rank > 1, non-finite floats) and returns
     * the List-path answer for them, which arrives here as an ordinary List of
     * pairs and is rewritten the same way. */
    if (is_ndarray(list)) { Expr* r = counts_from_ndarray(list); if (r) return r; }
    /* Counts over an association tallies its values (Counts[Values[assoc]]). */
    if (is_association(list)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }
    if (!head_is(list, SYM_List)) return NULL;
    size_t n = list->data.function.arg_count;

    KeyIndex ki;
    if (!ki_init(&ki, n)) return NULL;
    Expr** keys  = malloc(sizeof(Expr*) * (n ? n : 1)); /* borrowed */
    int64_t* cnt = malloc(sizeof(int64_t) * (n ? n : 1));
    size_t ndistinct = 0;

    for (size_t i = 0; i < n; i++) {
        Expr* x = list->data.function.args[i];
        size_t slot, idx = ki_lookup(&ki, keys, x, &slot);
        if (idx == SIZE_MAX) {
            keys[ndistinct] = x;
            cnt[ndistinct] = 1;
            ki_insert(&ki, slot, ndistinct);
            ndistinct++;
        } else {
            cnt[idx]++;
        }
    }

    Expr** rules = malloc(sizeof(Expr*) * (ndistinct ? ndistinct : 1));
    for (size_t i = 0; i < ndistinct; i++)
        rules[i] = make_rule(expr_copy(keys[i]), expr_new_integer(cnt[i]));
    Expr* assoc = expr_new_function(expr_new_symbol(SYM_Association),rules, ndistinct);
    free(rules); free(keys); free(cnt); ki_free(&ki);
    return assoc;
}

/* ======================================================================
 * GroupBy[list, f] — <|f[x] -> {x, ...}, ...|> preserving first-key order.
 * f is applied once per element (O(n) evaluations) and grouping is
 * hash-indexed, so the whole operation is O(n) plus f's cost.
 * ====================================================================== */
Expr* builtin_groupby(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 2 && argc != 3) return NULL;
    Expr* list = res->data.function.args[0];
    Expr* f    = res->data.function.args[1];
    Expr* reducer = (argc == 3) ? res->data.function.args[2] : NULL;
    /* GroupBy[assoc, f]: group the entries by f[value] into sub-associations
     * (keys preserved); GroupBy[assoc, f, red] reduces each sub-association. The
     * keyfn -> valfn transform form is list-only, so leave it unevaluated on an
     * association. */
    bool assoc_in = is_association(list);
    if (!head_is(list, SYM_List) && !assoc_in) return NULL;
    if (assoc_in && is_rule2(f)) return NULL;
    size_t n = list->data.function.arg_count;

    /* GroupBy[list, keyfn -> valfn, ...]: group by keyfn[x] but collect
     * valfn[x] in each group (Wolfram's key -> value transform form). */
    Expr* keyfn = f;
    Expr* valfn = NULL;
    if (is_rule2(f)) { keyfn = rule_key(f); valfn = rule_val(f); }

    KeyIndex ki;
    if (!ki_init(&ki, n)) return NULL;
    Expr** keys   = malloc(sizeof(Expr*) * (n ? n : 1)); /* owned group keys */
    Expr*** groups = malloc(sizeof(Expr**) * (n ? n : 1)); /* owned element copies */
    size_t* gcap  = malloc(sizeof(size_t) * (n ? n : 1));
    size_t* gcnt  = malloc(sizeof(size_t) * (n ? n : 1));
    size_t ngroups = 0;

    for (size_t i = 0; i < n; i++) {
        Expr* x = list->data.function.args[i];
        /* For an association, group by f applied to the value; the entry (the
         * whole Rule) is what gets collected so keys survive in each group. */
        Expr* key_input = assoc_in ? rule_val(x) : x;
        Expr* fx_args[1] = { expr_copy(key_input) };
        Expr* fx = expr_new_function(expr_copy(keyfn), fx_args, 1);
        Expr* key = evaluate(fx);      /* the group key */
        expr_free(fx);

        size_t slot, idx = ki_lookup(&ki, keys, key, &slot);
        if (idx == SIZE_MAX) {
            keys[ngroups] = key;       /* adopt */
            gcap[ngroups] = 4;
            gcnt[ngroups] = 0;
            groups[ngroups] = malloc(sizeof(Expr*) * gcap[ngroups]);
            ki_insert(&ki, slot, ngroups);
            idx = ngroups;
            ngroups++;
        } else {
            expr_free(key);            /* key already recorded */
        }
        if (gcnt[idx] == gcap[idx]) {
            gcap[idx] *= 2;
            groups[idx] = realloc(groups[idx], sizeof(Expr*) * gcap[idx]);
        }
        /* Collect valfn[x] when a value transform was given, else x itself. */
        if (valfn) {
            Expr* vx_args[1] = { expr_copy(x) };
            Expr* vx = expr_new_function(expr_copy(valfn), vx_args, 1);
            groups[idx][gcnt[idx]++] = evaluate(vx);
            expr_free(vx);
        } else {
            groups[idx][gcnt[idx]++] = expr_copy(x);
        }
    }

    Expr** rules = malloc(sizeof(Expr*) * (ngroups ? ngroups : 1));
    for (size_t i = 0; i < ngroups; i++) {
        /* Each group is a sub-association for an association input, else a list. */
        Expr* group_list = assoc_in
            ? expr_new_function(expr_new_symbol(SYM_Association),groups[i], gcnt[i])
            : make_list(groups[i], gcnt[i]);
        /* GroupBy[list, f, g] applies the reducer g to each group; the
         * g[{group}] application is left for the evaluator to reduce. */
        Expr* value = reducer
            ? expr_new_function(expr_copy(reducer), &group_list, 1)
            : group_list;
        rules[i] = make_rule(keys[i], value); /* adopts key + value */
        free(groups[i]);
    }
    Expr* assoc = expr_new_function(expr_new_symbol(SYM_Association),rules, ngroups);
    free(rules); free(keys); free(groups); free(gcap); free(gcnt); ki_free(&ki);
    return assoc;
}

/* ======================================================================
 * GatherBy[list, f] / Gather[list] — {group1, group2, ...}: gather elements
 * with equal f[x] into sublists, in first-appearance order (like GroupBy but
 * returning the groups as a plain list of lists). Hash-indexed, O(n) plus the
 * cost of f.
 *
 * assoc_gather_core is the single grouping engine behind both builtins. A NULL
 * `f` means "group by the element itself" — the identity key function — which
 * is exactly what Gather[list] needs; taking that path as a NULL check rather
 * than as f = Identity avoids n redundant Identity[x] evaluations and makes
 * Gather[l] === GatherBy[l, Identity] structural rather than coincidental.
 * ====================================================================== */
Expr* assoc_gather_core(Expr* list, Expr* f) {
    /* GatherBy[assoc, f] gathers the entries by f[value] into sub-associations
     * (keys preserved), returned as a list in first-appearance order -- like
     * GroupBy[assoc, f] but without the outer group keys. */
    bool assoc_in = is_association(list);
    if (!head_is(list, SYM_List) && !assoc_in) return NULL;
    size_t n = list->data.function.arg_count;

    KeyIndex ki;
    if (!ki_init(&ki, n)) return NULL;
    Expr**  keys   = malloc(sizeof(Expr*)  * (n ? n : 1)); /* owned group keys */
    Expr*** groups = malloc(sizeof(Expr**) * (n ? n : 1)); /* owned element copies */
    size_t* gcap   = malloc(sizeof(size_t) * (n ? n : 1));
    size_t* gcnt   = malloc(sizeof(size_t) * (n ? n : 1));
    size_t ng = 0;

    for (size_t i = 0; i < n; i++) {
        Expr* x = list->data.function.args[i];
        Expr* key_input = assoc_in ? rule_val(x) : x;   /* group by f[value] */
        /* Identity key (f == NULL): the element is its own group key, so no
         * function application is built or evaluated -- just a copy. */
        Expr* key;
        if (f) {
            Expr* fx_args[1] = { expr_copy(key_input) };
            Expr* fx = expr_new_function(expr_copy(f), fx_args, 1);
            key = evaluate(fx);
            expr_free(fx);
        } else {
            key = expr_copy(key_input);
        }
        size_t slot, idx = ki_lookup(&ki, keys, key, &slot);
        if (idx == SIZE_MAX) {
            keys[ng] = key; gcap[ng] = 4; gcnt[ng] = 0;
            groups[ng] = malloc(sizeof(Expr*) * gcap[ng]);
            ki_insert(&ki, slot, ng); idx = ng; ng++;
        } else {
            expr_free(key);
        }
        if (gcnt[idx] == gcap[idx]) { gcap[idx] *= 2; groups[idx] = realloc(groups[idx], sizeof(Expr*) * gcap[idx]); }
        groups[idx][gcnt[idx]++] = expr_copy(x);
    }

    Expr** outer = malloc(sizeof(Expr*) * (ng ? ng : 1));
    for (size_t i = 0; i < ng; i++) {
        outer[i] = assoc_in
            ? expr_new_function(expr_new_symbol(SYM_Association),groups[i], gcnt[i])
            : make_list(groups[i], gcnt[i]);
        expr_free(keys[i]);   /* the groups are returned bare; keys are dropped */
        free(groups[i]);
    }
    Expr* result = make_list(outer, ng);
    free(outer); free(keys); free(groups); free(gcap); free(gcnt); ki_free(&ki);
    return result;
}

Expr* builtin_gatherby(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    return assoc_gather_core(res->data.function.args[0],
                             res->data.function.args[1]);
}

/* ======================================================================
 * Merge[{assoc1, assoc2, ...}, f] — combine associations, applying f to the
 * List of values collected for each key (in first-seen key order).
 * ====================================================================== */
Expr* builtin_merge(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* col = res->data.function.args[0];
    Expr* f   = res->data.function.args[1];
    if (!head_is(col, SYM_List)) {
        /* Merge also accepts an association *of* associations — merge its
         * values. Delegate to the list form via Merge[Values[col], f]. */
        if (is_association(col)) {
            Expr* vals = assoc_values_list(col);
            Expr* margs[2] = { vals, expr_copy(f) };
            Expr* call = expr_new_function(expr_new_symbol(SYM_Merge), margs, 2);
            Expr* r = evaluate(call);
            expr_free(call);
            return r;
        }
        return NULL;
    }

    /* Upper bound on distinct keys = total rules across all associations. */
    size_t total = 0;
    for (size_t i = 0; i < col->data.function.arg_count; i++) {
        Expr* a = col->data.function.args[i];
        if (!is_association(a)) return NULL;
        total += a->data.function.arg_count;
    }

    KeyIndex ki;
    if (!ki_init(&ki, total)) return NULL;
    Expr**  keys  = malloc(sizeof(Expr*)  * (total ? total : 1)); /* borrowed */
    Expr*** vals  = malloc(sizeof(Expr**) * (total ? total : 1)); /* owned copies */
    size_t* vcap  = malloc(sizeof(size_t) * (total ? total : 1));
    size_t* vcnt  = malloc(sizeof(size_t) * (total ? total : 1));
    size_t ndistinct = 0;

    for (size_t i = 0; i < col->data.function.arg_count; i++) {
        Expr* a = col->data.function.args[i];
        for (size_t j = 0; j < a->data.function.arg_count; j++) {
            Expr* r = a->data.function.args[j];
            Expr* k = rule_key(r);
            size_t slot, idx = ki_lookup(&ki, keys, k, &slot);
            if (idx == SIZE_MAX) {
                keys[ndistinct] = k;
                vcap[ndistinct] = 4;
                vcnt[ndistinct] = 0;
                vals[ndistinct] = malloc(sizeof(Expr*) * vcap[ndistinct]);
                ki_insert(&ki, slot, ndistinct);
                idx = ndistinct;
                ndistinct++;
            }
            if (vcnt[idx] == vcap[idx]) {
                vcap[idx] *= 2;
                vals[idx] = realloc(vals[idx], sizeof(Expr*) * vcap[idx]);
            }
            vals[idx][vcnt[idx]++] = expr_copy(rule_val(r));
        }
    }

    Expr** rules = malloc(sizeof(Expr*) * (ndistinct ? ndistinct : 1));
    for (size_t i = 0; i < ndistinct; i++) {
        Expr* vlist = make_list(vals[i], vcnt[i]);
        Expr* fargs[1] = { vlist };
        Expr* fapp = expr_new_function(expr_copy(f), fargs, 1); /* f[{v...}] */
        rules[i] = make_rule(expr_copy(keys[i]), fapp);
        free(vals[i]);
    }
    Expr* assoc = expr_new_function(expr_new_symbol(SYM_Association),rules, ndistinct);
    free(rules); free(keys); free(vals); free(vcap); free(vcnt); ki_free(&ki);
    return assoc;
}

/* ======================================================================
 * AssociateTo[symbol, key->val | {rules}] — HoldFirst in-place update.
 * Mirrors AppendTo: read the symbol's current association, produce the
 * updated one, and assign it back.
 * ====================================================================== */
Expr* builtin_associate_to(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* sym = res->data.function.args[0];
    Expr* kv  = res->data.function.args[1];
    if (sym->type != EXPR_SYMBOL) return NULL;

    Expr* current = evaluate(sym);
    if (!current || !is_association(current)) { if (current) expr_free(current); return NULL; }

    /* Gather current rules plus the new one(s), then re-canonicalise. */
    size_t base = current->data.function.arg_count;
    size_t extra = head_is(kv, SYM_List) ? kv->data.function.arg_count : 1;
    Expr** all = malloc(sizeof(Expr*) * (base + extra));
    size_t n = 0;
    for (size_t i = 0; i < base; i++) all[n++] = current->data.function.args[i];
    if (head_is(kv, SYM_List)) {
        for (size_t i = 0; i < extra; i++) {
            if (!is_rule2(kv->data.function.args[i])) { free(all); expr_free(current); return NULL; }
            all[n++] = kv->data.function.args[i];
        }
    } else {
        if (!is_rule2(kv)) { free(all); expr_free(current); return NULL; }
        all[n++] = kv;
    }

    Expr* updated = assoc_from_rules(all, n);
    free(all);
    expr_free(current);

    /* Assign back to the symbol (HoldFirst guarantees `sym` is the symbol). */
    symtab_add_own_value(sym->data.symbol.name, sym, updated);
    return updated;
}

/* Apply f to a single argument and evaluate: evaluate(f[arg]). Owns result. */
static Expr* apply1(Expr* f, const Expr* arg) {
    Expr* a = expr_copy((Expr*)arg);
    Expr* call = expr_new_function(expr_copy(f), &a, 1);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* ======================================================================
 * Map / Select threading over association values (Wolfram semantics):
 * Map[f, <|k -> v|>]     -> <|k -> f[v]|>   (keys preserved)
 * Select[<|k -> v|>, p]  -> the entries where p[v] is True.
 * ====================================================================== */
Expr* assoc_map_values(Expr* f, const Expr* assoc) {
    size_t n = assoc->data.function.arg_count;
    Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr* r = assoc->data.function.args[i];
        /* is_association() only checks the head, and builtin_association leaves
         * a malformed Association[a, b] intact, so an entry is not guaranteed
         * to be a rule. Reading args[0]/args[1] of a symbol type-puns its
         * SymbolDef*; pass anything that is not a rule through untouched. */
        if (!is_rule2(r)) { out[i] = expr_copy(r); continue; }
        Expr* fv_arg = expr_copy(rule_val(r));
        Expr* fv = expr_new_function(expr_copy(f), &fv_arg, 1); /* f[v], evaluated later */
        out[i] = make_rule(expr_copy(rule_key(r)), fv);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_Association),out, n);
    free(out);
    return result;
}

Expr* assoc_select_values(Expr* pred, const Expr* assoc, int64_t max) {
    size_t n = assoc->data.function.arg_count;
    Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
    size_t nout = 0;
    for (size_t i = 0; i < n; i++) {
        if (max >= 0 && (int64_t)nout >= max) break;
        Expr* r = assoc->data.function.args[i];
        Expr* verdict = apply1(pred, rule_val(r));
        if (verdict->type == EXPR_SYMBOL && verdict->data.symbol.name == SYM_True)
            out[nout++] = expr_copy(r);
        expr_free(verdict);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_Association),out, nout);
    free(out);
    return result;
}

/* ======================================================================
 * Aggregation over values: Total / Min / Max / Mean / ... act on the values
 * of an association. assoc_apply_over_values rewrites head[assoc, rest] as
 * head[Values[assoc], rest] and re-evaluates, reusing the existing builtin.
 * ====================================================================== */
Expr* assoc_values_list(const Expr* assoc) {
    size_t n = assoc->data.function.arg_count;
    Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        out[i] = expr_copy(rule_val(assoc->data.function.args[i]));
    Expr* list = make_list(out, n);
    free(out);
    return list;
}

Expr* assoc_apply_over_values(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;
    if (!is_association(res->data.function.args[0])) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = malloc(sizeof(Expr*) * argc);
    args[0] = assoc_values_list(res->data.function.args[0]);
    for (size_t i = 1; i < argc; i++) args[i] = expr_copy(res->data.function.args[i]);
    Expr* call = expr_new_function(expr_copy(res->data.function.head), args, argc);
    free(args);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Re-pair the elements of a values List with the *trailing* keys of assoc,
 * producing an association: length n_result aligns with the last n_result keys,
 * so a windowed transform that dropped leading values drops the leading keys.
 * Returns NULL when `values` is not a List or is longer than assoc. */
Expr* assoc_rekey_from_list(const Expr* assoc, const Expr* values) {
    if (!head_is(values, SYM_List)) return NULL;
    size_t na   = assoc->data.function.arg_count;
    size_t nres = values->data.function.arg_count;
    if (nres > na) return NULL;
    size_t offset = na - nres;                     /* leading keys dropped, if any */
    Expr** entries = malloc(sizeof(Expr*) * (nres ? nres : 1));
    for (size_t i = 0; i < nres; i++) {
        Expr* key = rule_key(assoc->data.function.args[offset + i]);
        Expr* rargs[2] = { expr_copy(key), expr_copy(values->data.function.args[i]) };
        entries[i] = expr_new_function(expr_new_symbol(SYM_Rule), rargs, 2);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_Association),entries, nres);
    free(entries);
    return result;
}

/* Like assoc_apply_over_values, but for windowed value transforms whose result
 * is a List that stays parallel to the association's entries (Accumulate keeps
 * every entry; Differences drops the leading one). Returns NULL (leave to the
 * caller) when the first argument is not an association or the transform did
 * not yield a same-or-shorter List. */
Expr* assoc_rekey_over_values(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;
    Expr* assoc = res->data.function.args[0];
    if (!is_association(assoc)) return NULL;

    Expr* out = assoc_apply_over_values(res);      /* head[Values[assoc], ...] */
    if (!out) return NULL;
    Expr* result = assoc_rekey_from_list(assoc, out);
    expr_free(out);
    return result;                                 /* NULL if out was not a valid List */
}

/* ======================================================================
 * Sort[assoc] — reorder entries by their values in canonical order.
 * ====================================================================== */
/* A rule paired with its original position, so ties can be broken by position
 * — giving a stable sort by value (Sort[assoc] preserves input order for
 * equal values, matching Wolfram). */
typedef struct { Expr* rule; size_t idx; } RuleWithIndex;

static int rule_value_stable_cmp(const void* a, const void* b) {
    const RuleWithIndex* ra = (const RuleWithIndex*)a;
    const RuleWithIndex* rb = (const RuleWithIndex*)b;
    int c = expr_compare(rule_val(ra->rule), rule_val(rb->rule));
    if (c != 0) return c;
    return (ra->idx < rb->idx) ? -1 : (ra->idx > rb->idx ? 1 : 0);
}

Expr* assoc_sort_by_value(const Expr* assoc) {
    size_t n = assoc->data.function.arg_count;
    RuleWithIndex* pairs = malloc(sizeof(RuleWithIndex) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) { pairs[i].rule = assoc->data.function.args[i]; pairs[i].idx = i; }
    qsort(pairs, n, sizeof(RuleWithIndex), rule_value_stable_cmp);
    Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) out[i] = expr_copy(pairs[i].rule);
    Expr* result = expr_new_function(expr_new_symbol(SYM_Association),out, n);
    free(out); free(pairs);
    return result;
}

/* DeleteCases[assoc, patt] — drop entries whose VALUE matches patt, testing each
 * value with MatchQ[value, patt]. */
Expr* assoc_delete_cases(const Expr* assoc, Expr* pattern) {
    size_t n = assoc->data.function.arg_count;
    Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
    size_t nout = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* r = assoc->data.function.args[i];
        Expr* mq_args[2] = { expr_copy(rule_val(r)), expr_copy(pattern) };
        Expr* mq = expr_new_function(expr_new_symbol(SYM_MatchQ), mq_args, 2);
        Expr* verdict = evaluate(mq);
        expr_free(mq);
        bool matches = (verdict->type == EXPR_SYMBOL && verdict->data.symbol.name == SYM_True);
        expr_free(verdict);
        if (!matches) out[nout++] = expr_copy(r);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_Association),out, nout);
    free(out);
    return result;
}

/* DeleteDuplicates[assoc]: keep the first entry for each distinct value.
 * The seen values are indexed in a transient hash table so the whole scan is
 * O(n) rather than the O(n^2) of comparing every value against every kept one. */
Expr* assoc_delete_duplicate_values(const Expr* assoc) {
    size_t n = assoc->data.function.arg_count;
    KeyIndex ki;
    if (!ki_init(&ki, n)) return NULL;
    Expr** seen = malloc(sizeof(Expr*) * (n ? n : 1)); /* borrowed distinct values */
    size_t nseen = 0;
    Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));  /* kept rule copies */
    size_t nout = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* r = assoc->data.function.args[i];
        if (!is_rule2(r)) { out[nout++] = expr_copy(r); continue; }
        Expr* v = rule_val(r);
        size_t slot, idx = ki_lookup(&ki, seen, v, &slot);
        if (idx == SIZE_MAX) {          /* first time we see this value: keep it */
            seen[nseen] = v;
            ki_insert(&ki, slot, nseen);
            nseen++;
            out[nout++] = expr_copy(r);
        }
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_Association),out, nout);
    free(out);
    free(seen);
    ki_free(&ki);
    return result;
}

/* ======================================================================
 * KeySort[assoc] / KeySortBy[assoc, f] — order entries by key (or by f[key]).
 * ====================================================================== */
static int rule_key_cmp(const void* a, const void* b) {
    const Expr* ra = *(const Expr* const*)a;
    const Expr* rb = *(const Expr* const*)b;
    return expr_compare(rule_key(ra), rule_key(rb));
}

Expr* builtin_keysort(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* assoc = res->data.function.args[0];
    if (!is_association(assoc)) return NULL;
    size_t n = assoc->data.function.arg_count;
    Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) out[i] = expr_copy(assoc->data.function.args[i]);
    qsort(out, n, sizeof(Expr*), rule_key_cmp);   /* keys are distinct: total order */
    Expr* result = expr_new_function(expr_new_symbol(SYM_Association),out, n);
    free(out);
    return result;
}

Expr* builtin_keysortby(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* assoc = res->data.function.args[0];
    Expr* f     = res->data.function.args[1];
    if (!is_association(assoc)) return NULL;
    size_t n = assoc->data.function.arg_count;

    Expr** out  = malloc(sizeof(Expr*) * (n ? n : 1));
    Expr** fkey = malloc(sizeof(Expr*) * (n ? n : 1)); /* f[key], owned */
    for (size_t i = 0; i < n; i++) {
        out[i]  = expr_copy(assoc->data.function.args[i]);
        fkey[i] = apply1(f, rule_key(assoc->data.function.args[i]));
    }
    /* Stable insertion sort keyed by fkey (ties keep association order). */
    for (size_t i = 1; i < n; i++) {
        Expr* ro = out[i]; Expr* fo = fkey[i];
        size_t j = i;
        while (j > 0 && expr_compare(fkey[j - 1], fo) > 0) {
            out[j] = out[j - 1]; fkey[j] = fkey[j - 1]; j--;
        }
        out[j] = ro; fkey[j] = fo;
    }
    for (size_t i = 0; i < n; i++) expr_free(fkey[i]);
    free(fkey);
    Expr* result = expr_new_function(expr_new_symbol(SYM_Association),out, n);
    free(out);
    return result;
}

/* ======================================================================
 * KeyMap[f, assoc] — apply f to each key (values kept); re-canonicalise since
 * f may map distinct keys together.
 * ====================================================================== */
Expr* builtin_keymap(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* f     = res->data.function.args[0];
    Expr* assoc = res->data.function.args[1];
    if (!is_association(assoc)) return NULL;
    size_t n = assoc->data.function.arg_count;
    Expr** rules = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr* r = assoc->data.function.args[i];
        Expr* newkey = apply1(f, rule_key(r));
        rules[i] = make_rule(newkey, expr_copy(rule_val(r)));
    }
    Expr* result = assoc_from_rules(rules, n);
    for (size_t i = 0; i < n; i++) expr_free(rules[i]);
    free(rules);
    return result;
}

/* ======================================================================
 * KeySelect[assoc, pred] — keep entries whose key satisfies pred.
 * ====================================================================== */
Expr* builtin_keyselect(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* assoc = res->data.function.args[0];
    Expr* pred  = res->data.function.args[1];
    if (!is_association(assoc)) return NULL;
    size_t n = assoc->data.function.arg_count;
    Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
    size_t nout = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* r = assoc->data.function.args[i];
        Expr* verdict = apply1(pred, rule_key(r));
        if (verdict->type == EXPR_SYMBOL && verdict->data.symbol.name == SYM_True)
            out[nout++] = expr_copy(r);
        expr_free(verdict);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_Association),out, nout);
    free(out);
    return result;
}

/* ======================================================================
 * CountsBy[list, f] — <|f[x] -> count of elements with that f-value|>.
 * ====================================================================== */
Expr* builtin_countsby(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* list = res->data.function.args[0];
    Expr* f    = res->data.function.args[1];
    /* CountsBy over an association tallies its values by f. */
    if (is_association(list)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }
    if (!head_is(list, SYM_List)) return NULL;
    size_t n = list->data.function.arg_count;

    KeyIndex ki;
    if (!ki_init(&ki, n)) return NULL;
    Expr** keys  = malloc(sizeof(Expr*) * (n ? n : 1)); /* owned f-values */
    int64_t* cnt = malloc(sizeof(int64_t) * (n ? n : 1));
    size_t nd = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* key = apply1(f, list->data.function.args[i]);
        size_t slot, idx = ki_lookup(&ki, keys, key, &slot);
        if (idx == SIZE_MAX) { keys[nd] = key; cnt[nd] = 1; ki_insert(&ki, slot, nd); nd++; }
        else { cnt[idx]++; expr_free(key); }
    }
    Expr** rules = malloc(sizeof(Expr*) * (nd ? nd : 1));
    for (size_t i = 0; i < nd; i++) rules[i] = make_rule(keys[i], expr_new_integer(cnt[i]));
    Expr* assoc = expr_new_function(expr_new_symbol(SYM_Association),rules, nd);
    free(rules); free(keys); free(cnt); ki_free(&ki);
    return assoc;
}

/* ======================================================================
 * PositionIndex[list] — <|value -> {1-based positions}|>. Hash-indexed, O(n).
 * ====================================================================== */
Expr* builtin_positionindex(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* list = res->data.function.args[0];
    if (!head_is(list, SYM_List)) return NULL;
    size_t n = list->data.function.arg_count;

    KeyIndex ki;
    if (!ki_init(&ki, n)) return NULL;
    Expr**  keys = malloc(sizeof(Expr*)  * (n ? n : 1)); /* borrowed */
    Expr*** pos  = malloc(sizeof(Expr**) * (n ? n : 1)); /* owned integer positions */
    size_t* pcap = malloc(sizeof(size_t) * (n ? n : 1));
    size_t* pcnt = malloc(sizeof(size_t) * (n ? n : 1));
    size_t nd = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* x = list->data.function.args[i];
        size_t slot, idx = ki_lookup(&ki, keys, x, &slot);
        if (idx == SIZE_MAX) {
            keys[nd] = x; pcap[nd] = 4; pcnt[nd] = 0;
            pos[nd] = malloc(sizeof(Expr*) * pcap[nd]);
            ki_insert(&ki, slot, nd); idx = nd; nd++;
        }
        if (pcnt[idx] == pcap[idx]) { pcap[idx] *= 2; pos[idx] = realloc(pos[idx], sizeof(Expr*) * pcap[idx]); }
        pos[idx][pcnt[idx]++] = expr_new_integer((int64_t)i + 1);
    }
    Expr** rules = malloc(sizeof(Expr*) * (nd ? nd : 1));
    for (size_t i = 0; i < nd; i++) {
        Expr* plist = make_list(pos[i], pcnt[i]);
        rules[i] = make_rule(expr_copy(keys[i]), plist);
        free(pos[i]);
    }
    Expr* assoc = expr_new_function(expr_new_symbol(SYM_Association),rules, nd);
    free(rules); free(keys); free(pos); free(pcap); free(pcnt); ki_free(&ki);
    return assoc;
}

/* ======================================================================
 * AssociationMap[f, {k1, ...}] — <|k1 -> f[k1], ...|>.
 * ====================================================================== */
Expr* builtin_associationmap(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* f    = res->data.function.args[0];
    Expr* keys = res->data.function.args[1];
    if (!head_is(keys, SYM_List)) return NULL;
    size_t n = keys->data.function.arg_count;
    Expr** rules = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr* k = keys->data.function.args[i];
        Expr* karg = expr_copy(k);
        Expr* fv = expr_new_function(expr_copy(f), &karg, 1); /* f[k], evaluated later */
        rules[i] = make_rule(expr_copy(k), fv);
    }
    Expr* assoc = assoc_from_rules(rules, n);
    for (size_t i = 0; i < n; i++) expr_free(rules[i]);
    free(rules);
    return assoc;
}

/* ======================================================================
 * Registration.
 * ====================================================================== */
void assoc_init(void) {
    symtab_add_builtin("Association", builtin_association);
    symtab_get_def("Association")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Association",
        "Association[key1 -> val1, key2 -> val2, ...]  (also written <|...|>)\n"
        "\tRepresents an association mapping keys to values with unique,\n"
        "\tinsertion-ordered keys (last value wins on duplicates).\n"
        "\tArguments may be rules, lists of rules, or other associations.");

    symtab_add_builtin("AssociationQ", builtin_associationq);
    symtab_get_def("AssociationQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("AssociationQ",
        "AssociationQ[expr]\n\tGives True if expr is an Association, else False.");

    symtab_add_builtin("Keys", builtin_keys);
    symtab_get_def("Keys")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Keys",
        "Keys[assoc]\n\tGives a list of the keys of an association (or rules).");

    symtab_add_builtin("Values", builtin_values);
    symtab_get_def("Values")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Values",
        "Values[assoc]\n\tGives a list of the values of an association (or rules).");

    symtab_add_builtin("Lookup", builtin_lookup);
    symtab_get_def("Lookup")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Lookup",
        "Lookup[assoc, key]\n\tGives the value for key, or Missing[\"KeyAbsent\", key].\n"
        "Lookup[assoc, key, default]\n\tUses default when key is absent.\n"
        "Lookup[assoc, {k1, k2, ...}]\n\tLooks up several keys at once (O(n+m)).");

    symtab_add_builtin("KeyExistsQ", builtin_keyexistsq);
    symtab_get_def("KeyExistsQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KeyExistsQ",
        "KeyExistsQ[assoc, key]\n\tGives True if key is present in assoc, else False.");

    symtab_add_builtin("KeyMemberQ", builtin_keymemberq);
    symtab_get_def("KeyMemberQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KeyMemberQ",
        "KeyMemberQ[assoc, key]\n\tGives True if key is present in assoc (same as\n"
        "\tKeyExistsQ), else False.");

    symtab_add_builtin("KeyFreeQ", builtin_keyfreeq);
    symtab_get_def("KeyFreeQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KeyFreeQ",
        "KeyFreeQ[assoc, key]\n\tGives True if key is absent from assoc (the\n"
        "\tcomplement of KeyExistsQ), else False.");

    symtab_add_builtin("KeyDrop", builtin_keydrop);
    symtab_get_def("KeyDrop")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KeyDrop",
        "KeyDrop[assoc, key]  |  KeyDrop[assoc, {k1, ...}]\n"
        "\tGives assoc with the specified keys removed (order preserved).");

    symtab_add_builtin("KeyTake", builtin_keytake);
    symtab_get_def("KeyTake")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KeyTake",
        "KeyTake[assoc, {k1, ...}]\n"
        "\tGives the association of only the specified keys (order preserved).");

    symtab_add_builtin("KeyUnion", builtin_keyunion);
    symtab_get_def("KeyUnion")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KeyUnion",
        "KeyUnion[{assoc1, assoc2, ...}]\n"
        "\tGives the list of associations padded to the union of all their keys;\n"
        "\ta key absent from an association is filled with Missing[\"KeyAbsent\", key].");

    symtab_add_builtin("KeyValueMap", builtin_keyvaluemap);
    symtab_get_def("KeyValueMap")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KeyValueMap",
        "KeyValueMap[f, assoc]\n\tGives {f[k1, v1], f[k2, v2], ...}.");

    symtab_add_builtin("AssociationThread", builtin_associationthread);
    symtab_get_def("AssociationThread")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("AssociationThread",
        "AssociationThread[{k...}, {v...}]  |  AssociationThread[keys -> values]\n"
        "\tBuilds <|k1 -> v1, ...|> from parallel key and value lists.");

    symtab_add_builtin("Counts", builtin_counts);
    symtab_get_def("Counts")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Counts",
        "Counts[list]\n\tGives <|element -> count, ...|> tallying each distinct\n"
        "\telement. Hash-indexed: one O(n) pass.");

    symtab_add_builtin("GroupBy", builtin_groupby);
    symtab_get_def("GroupBy")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GroupBy",
        "GroupBy[list, f]\n\tGroups the elements of list by the value of f[element],\n"
        "\tgiving <|f[x] -> {matching elements}, ...|>.\n"
        "GroupBy[list, f, g]\n\tApplies the reducer g to each group, giving\n"
        "\t<|f[x] -> g[{matching elements}], ...|> (split-apply-combine).");

    symtab_add_builtin("GatherBy", builtin_gatherby);
    symtab_get_def("GatherBy")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GatherBy",
        "GatherBy[list, f]\n\tGathers elements with equal f[element] into sublists,\n"
        "\tin first-appearance order: {{group1}, {group2}, ...}.");

    symtab_add_builtin("Merge", builtin_merge);
    symtab_get_def("Merge")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Merge",
        "Merge[{assoc1, assoc2, ...}, f]\n"
        "\tCombines associations, applying f to the list of values collected\n"
        "\tfor each key (e.g. Merge[{...}, Total]).");

    symtab_add_builtin("AssociateTo", builtin_associate_to);
    symtab_get_def("AssociateTo")->attributes |= ATTR_HOLDFIRST | ATTR_PROTECTED;
    symtab_set_docstring("AssociateTo",
        "AssociateTo[s, key -> val]  |  AssociateTo[s, {rules}]\n"
        "\tAdds or updates key-value pairs in the association held by symbol s,\n"
        "\tmodifying s in place.");

    symtab_add_builtin("KeySort", builtin_keysort);
    symtab_get_def("KeySort")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KeySort",
        "KeySort[assoc]\n\tSorts an association into canonical key order.");

    symtab_add_builtin("KeySortBy", builtin_keysortby);
    symtab_get_def("KeySortBy")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KeySortBy",
        "KeySortBy[assoc, f]\n\tSorts an association by f applied to each key (stable).");

    symtab_add_builtin("KeyMap", builtin_keymap);
    symtab_get_def("KeyMap")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KeyMap",
        "KeyMap[f, assoc]\n\tApplies f to every key, keeping values (keys that\n"
        "\tcollide collapse with last-value-wins).");

    symtab_add_builtin("KeySelect", builtin_keyselect);
    symtab_get_def("KeySelect")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KeySelect",
        "KeySelect[assoc, pred]\n\tKeeps the entries whose key satisfies pred.");

    symtab_add_builtin("CountsBy", builtin_countsby);
    symtab_get_def("CountsBy")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("CountsBy",
        "CountsBy[list, f]\n\tGives <|f[x] -> count, ...|> tallying elements by f[x].");

    symtab_add_builtin("PositionIndex", builtin_positionindex);
    symtab_get_def("PositionIndex")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("PositionIndex",
        "PositionIndex[list]\n\tGives <|value -> {positions}|> mapping each distinct\n"
        "\telement to the list of 1-based positions where it occurs. O(n).");

    symtab_add_builtin("AssociationMap", builtin_associationmap);
    symtab_get_def("AssociationMap")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("AssociationMap",
        "AssociationMap[f, {k1, k2, ...}]\n\tGives <|k1 -> f[k1], k2 -> f[k2], ...|>.");
}
