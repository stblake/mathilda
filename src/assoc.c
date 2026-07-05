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
static bool is_rule2(const Expr* e) {
    return (head_is(e, SYM_Rule) || head_is(e, SYM_RuleDelayed)) &&
           e->data.function.arg_count == 2;
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

    Expr* assoc = expr_new_function(expr_new_symbol(SYM_Association), out, nout);
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
        if (expr_eq(rule_key(r), key)) return rule_val(r);
    }
    return NULL;
}

Expr* builtin_lookup(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;
    Expr* assoc = res->data.function.args[0];
    Expr* key   = res->data.function.args[1];
    Expr* deflt = argc == 3 ? res->data.function.args[2] : NULL;
    if (!is_association(assoc)) return NULL;

    /* Lookup over a list of keys: single index build, then O(1) per key. */
    if (head_is(key, SYM_List)) {
        size_t na = assoc->data.function.arg_count;
        size_t nk = key->data.function.arg_count;
        KeyIndex ki;
        if (!ki_init(&ki, na)) return NULL;
        Expr** akeys = malloc(sizeof(Expr*) * (na ? na : 1));
        for (size_t i = 0; i < na; i++) {
            Expr* r = assoc->data.function.args[i];
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

    Expr* v = assoc_scan(assoc, key);
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
    if (!is_association(assoc)) return NULL;
    return expr_new_symbol(assoc_scan(assoc, key) ? SYM_True : SYM_False);
}

/* ======================================================================
 * KeyDrop[assoc, key|{keys}] / KeyTake[assoc, key|{keys}].
 * Both preserve association order; a drop/keep set is indexed once.
 * ====================================================================== */
static Expr* key_drop_take(Expr* res, bool take) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* assoc = res->data.function.args[0];
    Expr* karg  = res->data.function.args[1];
    if (!is_association(assoc)) return NULL;

    /* Normalise the key argument to an array of key pointers. */
    Expr** wanted; size_t nwanted;
    if (head_is(karg, SYM_List)) {
        nwanted = karg->data.function.arg_count;
        wanted = karg->data.function.args;
    } else {
        nwanted = 1;
        wanted = &karg;
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

Expr* builtin_keydrop(Expr* res) { return key_drop_take(res, false); }
Expr* builtin_keytake(Expr* res) { return key_drop_take(res, true);  }

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
Expr* builtin_counts(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* list = res->data.function.args[0];
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
    Expr* assoc = expr_new_function(expr_new_symbol(SYM_Association), rules, ndistinct);
    free(rules); free(keys); free(cnt); ki_free(&ki);
    return assoc;
}

/* ======================================================================
 * GroupBy[list, f] — <|f[x] -> {x, ...}, ...|> preserving first-key order.
 * f is applied once per element (O(n) evaluations) and grouping is
 * hash-indexed, so the whole operation is O(n) plus f's cost.
 * ====================================================================== */
Expr* builtin_groupby(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* list = res->data.function.args[0];
    Expr* f    = res->data.function.args[1];
    if (!head_is(list, SYM_List)) return NULL;
    size_t n = list->data.function.arg_count;

    KeyIndex ki;
    if (!ki_init(&ki, n)) return NULL;
    Expr** keys   = malloc(sizeof(Expr*) * (n ? n : 1)); /* owned group keys */
    Expr*** groups = malloc(sizeof(Expr**) * (n ? n : 1)); /* owned element copies */
    size_t* gcap  = malloc(sizeof(size_t) * (n ? n : 1));
    size_t* gcnt  = malloc(sizeof(size_t) * (n ? n : 1));
    size_t ngroups = 0;

    for (size_t i = 0; i < n; i++) {
        Expr* x = list->data.function.args[i];
        Expr* fx_args[1] = { expr_copy(x) };
        Expr* fx = expr_new_function(expr_copy(f), fx_args, 1);
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
        groups[idx][gcnt[idx]++] = expr_copy(x);
    }

    Expr** rules = malloc(sizeof(Expr*) * (ngroups ? ngroups : 1));
    for (size_t i = 0; i < ngroups; i++) {
        Expr* group_list = make_list(groups[i], gcnt[i]);
        rules[i] = make_rule(keys[i], group_list); /* adopts key + list */
        free(groups[i]);
    }
    Expr* assoc = expr_new_function(expr_new_symbol(SYM_Association), rules, ngroups);
    free(rules); free(keys); free(groups); free(gcap); free(gcnt); ki_free(&ki);
    return assoc;
}

/* ======================================================================
 * Merge[{assoc1, assoc2, ...}, f] — combine associations, applying f to the
 * List of values collected for each key (in first-seen key order).
 * ====================================================================== */
Expr* builtin_merge(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* col = res->data.function.args[0];
    Expr* f   = res->data.function.args[1];
    if (!head_is(col, SYM_List)) return NULL;

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
    Expr* assoc = expr_new_function(expr_new_symbol(SYM_Association), rules, ndistinct);
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
    symtab_add_own_value(sym->data.symbol, sym, updated);
    return updated;
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
        "\tgiving <|f[x] -> {matching elements}, ...|>.");

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
}
