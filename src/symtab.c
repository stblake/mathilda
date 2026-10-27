#include "symtab.h"
#include "match.h"
#include "print.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "eval.h"
#include "arithmetic.h"
#include "options.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* Bucket count: a power of two, doubled whenever the table holds more than
 * twice as many nodes (symtab_grow).  Module mints a fresh symbol per local per
 * call and never removes it, so a long session interns 10^5-10^6 names; a fixed
 * 65535 buckets meant chains of dozens of cache-missing nodes on every new name
 * (issue #112).  Nodes never move in memory, so SymbolDef pointers stay valid. */
#define SYMTAB_INIT_SIZE 65536

/* Phase 2: one table backs both interning and definition lookup. Each bucket
 * chains SymbolDef nodes directly (the old SymEntry indirection is gone). A node
 * exists as soon as its name is interned; `has_def` marks whether it carries a
 * real definition. */
static SymbolDef** symtab = NULL;
static size_t symtab_size = 0;      /* buckets: a power of two */
static size_t symtab_nodes = 0;

static void symtab_alloc(void) {
    if (symtab) return;
    symtab_size = SYMTAB_INIT_SIZE;
    symtab = calloc(symtab_size, sizeof(SymbolDef*));
    if (!symtab) { fprintf(stderr, "Mathilda: out of memory (symbol table)\n"); abort(); }
}

/* Double the bucket array and relink every node by its stored hash. */
static void symtab_grow(void) {
    size_t ns = symtab_size * 2;
    SymbolDef** nb = calloc(ns, sizeof(SymbolDef*));
    if (!nb) return;                    /* keep the longer chains */
    for (size_t i = 0; i < symtab_size; i++) {
        SymbolDef* d = symtab[i];
        while (d) {
            SymbolDef* next = d->next;
            size_t j = d->name_hash & (ns - 1);
            d->next = nb[j];
            nb[j] = d;
            d = next;
        }
    }
    free(symtab);
    symtab = nb;
    symtab_size = ns;
}

static unsigned int name_hash_full(const char* s) {
    unsigned int h = 5381;
    while (*s) h = ((h << 5) + h) + *s++;
    return h;
}

static void free_rule_list(Rule* r) {
    while (r) {
        Rule* next = r->next;
        expr_free(r->pattern);
        expr_free(r->replacement);
        free(r);
        r = next;
    }
}

/* Running total of FormatValue rules (see symtab.h). */
size_t symtab_format_value_count = 0;
size_t symtab_n_value_count = 0;

/* Free an NValue list, keeping symtab_n_value_count in step. */
static void free_n_value_list(Rule* r) {
    while (r) {
        Rule* next = r->next;
        expr_free(r->pattern);
        expr_free(r->replacement);
        free(r);
        if (symtab_n_value_count) symtab_n_value_count--;
        r = next;
    }
}

/* Free a FormatValue list, keeping symtab_format_value_count in step. */
static void free_format_value_list(Rule* r) {
    while (r) {
        Rule* next = r->next;
        expr_free(r->pattern);
        expr_free(r->replacement);
        free(r);
        if (symtab_format_value_count) symtab_format_value_count--;
        r = next;
    }
}

#if UP_VALUES
/* Running total of UpValue rules (see symtab.h). */
size_t symtab_up_value_count = 0;

/* Free an UpValue list AND keep the global count in step. Every UpValue teardown
 * path (reset_node_payload, symtab_clear_symbol, symtab_set_up_values) routes
 * through here so the counter can never drift from the stored rules. */
static void free_up_value_list(Rule* r) {
    while (r) {
        Rule* next = r->next;
        expr_free(r->pattern);
        expr_free(r->replacement);
        free(r);
        if (symtab_up_value_count) symtab_up_value_count--;
        r = next;
    }
}
#endif

/* Reset a node to "interned only": free its definition payload but KEEP the node
 * and its canonical name string in the table (so SYM_* and any live Expr keep
 * valid pointers) and KEEP is_system. This is what symtab_clear / symtab_init /
 * symtab_remove_symbol leave behind -- observationally an undefined symbol. */
static void litidx_drop(SymbolDef* d);

static void reset_node_payload(SymbolDef* d) {
    litidx_drop(d);
    free_rule_list(d->own_values);   d->own_values = NULL;
    free_rule_list(d->down_values);  d->down_values = NULL;
    free_rule_list(d->sub_values);   d->sub_values = NULL;
#if UP_VALUES
    free_up_value_list(d->up_values); d->up_values = NULL;
#endif
    free_format_value_list(d->format_values); d->format_values = NULL;
    free_n_value_list(d->n_values); d->n_values = NULL;
    if (d->docstring) { free(d->docstring); d->docstring = NULL; }
    if (d->default_options) { expr_free(d->default_options); d->default_options = NULL; }
    d->builtin_func = NULL;
    d->attributes = 0;
    d->has_def = 0;
}

/* The single find-or-create used by both interning and definition lookup.
 * `materialize` marks the node as a real definition (has_def=1). An already
 * interned input is recognised by pointer identity first (the common hot case),
 * falling back to strcmp only for fresh string literals. */
/* Depth of intern_system_scope_push (see sym_intern.h). */
static int g_system_scope = 0;

void intern_system_scope_push(void) { g_system_scope++; }
void intern_system_scope_pop(void)  { if (g_system_scope > 0) g_system_scope--; }
int  intern_system_scope_depth(void) { return g_system_scope; }

/* A definition made inside a kernel-module load: flag the node System if its
 * bare name is spelled like a System name. */
static void note_materialized(SymbolDef* d) {
    const char* nm = d->symbol_name;
    if (g_system_scope > 0 && !d->has_def &&
        (nm[0] == '$' || (nm[0] >= 'A' && nm[0] <= 'Z')) && !strchr(nm, '`'))
        d->is_system = 1;
}

static SymbolDef* node_intern(const char* name, int materialize) {
    if (!name) return NULL;
    symtab_alloc();
    unsigned int h = name_hash_full(name);
    size_t idx = h & (symtab_size - 1);
    for (SymbolDef* d = symtab[idx]; d; d = d->next) {
        if (d->symbol_name == name ||
            (d->name_hash == h && strcmp(d->symbol_name, name) == 0)) {
            if (materialize) { note_materialized(d); d->has_def = 1; }
            return d;
        }
    }
    /* Not present: create a node owning a stable copy of the name. calloc zeroes
     * every field (own_values=NULL, attributes=0, is_system=0). */
    SymbolDef* d = calloc(1, sizeof(SymbolDef));
    size_t len = strlen(name);
    char* copy = malloc(len + 1);
    memcpy(copy, name, len + 1);
    d->symbol_name = copy;
    d->name_hash = h;
    if (materialize) note_materialized(d);
    d->has_def = (uint8_t)(materialize ? 1 : 0);
    d->next = symtab[idx];
    symtab[idx] = d;
    if (++symtab_nodes > 2 * symtab_size) symtab_grow();
    return d;
}

/* Non-creating find (any node, defined or interned-only). */
static SymbolDef* node_find(const char* name) {
    if (!name) return NULL;
    if (!symtab) return NULL;
    unsigned int h = name_hash_full(name);
    for (SymbolDef* d = symtab[h & (symtab_size - 1)]; d; d = d->next)
        if (d->symbol_name == name ||
            (d->name_hash == h && strcmp(d->symbol_name, name) == 0))
            return d;
    return NULL;
}

/* Bumped by every symtab_init.  A load-once memo (LoadModule's, the lazily
 * loaded Integrate packages') is valid only for the generation it was filled
 * in: a re-initialised table has lost the definitions the load created, and a
 * memo that survived it left those packages undefined -- their stage then
 * handed back its own unevaluated entry call as an answer. */
static unsigned g_symtab_generation = 0;
unsigned symtab_generation(void) { return g_symtab_generation; }

void symtab_init(void) {
    g_symtab_generation++;
    /* Reset every existing node to interned-only rather than memset-ing the
     * bucket array: the array now owns the interned name strings, and SYM_*
     * (populated by sym_names_init below) plus any live Expr point straight at
     * them, so a memset would strand those pointers. First call: table empty,
     * so this loop is a no-op. */
    for (size_t i = 0; i < symtab_size; i++)
        for (SymbolDef* d = symtab[i]; d; d = d->next)
            reset_node_payload(d);
    /* Populate the SYM_* canonical-pointer cache so that any code path that
     * compares e->data.symbol.name against SYM_* (match.c, eval.c, the arithmetic
     * heads, etc.) works as soon as the symtab is live -- not only after
     * core_init() has finished registering builtins. Tests that initialise just
     * the symtab depend on this. */
    sym_names_init();
}

SymbolDef* symtab_get_def(const char* symbol_name) {
    return node_intern(symbol_name, 1);
}

SymbolDef* symtab_lookup(const char* symbol_name) {
    /* Only materialized definitions "exist" for the context resolver and
     * Options[]; an interned-only node reads as absent, matching the pre-merge
     * behavior where such names had no SymbolDef at all. */
    SymbolDef* d = node_find(symbol_name);
    return (d && d->has_def) ? d : NULL;
}

/* ============================================================
 * String interning (Phase 2: merged into the symbol table).
 * intern_symbol canonicalizes a name to a stable pointer without materializing
 * a definition, so merely mentioning a symbol does not make it appear in
 * Names[]. These implement the sym_intern.h API over the unified table.
 * ============================================================ */
/* Direct-mapped cache of CANONICAL interned pointers.  Most calls pass a name
 * that is already the interned pointer (expr_new_symbol(SYM_Times), a copied
 * symbol's name), and hashing + strcmp-probing it again was ~10% of a symbolic
 * workload (issue #112).  Only pointers the table itself returned are stored,
 * so a hit means `name` IS the canonical pointer: an arbitrary caller buffer
 * that happens to be reused with other contents can never hit.  Interned names
 * live until intern_clear, which empties the cache. */
#define INTERN_PTR_CACHE 4096
static SymbolDef* g_intern_ptr_cache[INTERN_PTR_CACHE];   /* keyed by symbol_name */
static inline size_t intern_ptr_slot(const char* p) {
    return (size_t)(((uintptr_t)p >> 4) ^ ((uintptr_t)p >> 16)) & (INTERN_PTR_CACHE - 1);
}

/* The node for `name`, interning it if new (not materialised).  A pointer that
 * is already canonical hits the cache without hashing the string. */
SymbolDef* intern_symbol_node(const char* name) {
    if (!name) return NULL;
    SymbolDef* c = g_intern_ptr_cache[intern_ptr_slot(name)];
    if (c && c->symbol_name == name) return c;
    SymbolDef* d = node_intern(name, 0);
    if (!d) return NULL;
    g_intern_ptr_cache[intern_ptr_slot(d->symbol_name)] = d;
    return d;
}

const char* intern_symbol(const char* name) {
    SymbolDef* d = intern_symbol_node(name);
    return d ? d->symbol_name : NULL;
}

/* ============================================================
 * Module-local reclamation
 *
 * Module mints a fresh symbol per local per call (x$123, Temporary).  The
 * reference removes such a symbol once nothing refers to it; Mathilda used to
 * keep every one for the whole session -- +111,005 symbols and +44 MB after a
 * single issue #112 trace, and +548 bytes per Module call in general.
 *
 * Counting.  A local created by Module from a FRESH name is flagged
 * module_temp; from then on every live Expr symbol node naming it caches its
 * SymbolDef and is counted in temp_nodes (expr.c).  The references a symbol's
 * OWN definitions make (x$1 = {x$1}, the own-value pattern itself) must not
 * keep it alive, so a local is reclaimable when a walk of its own rule lists
 * finds ALL temp_nodes of its nodes, each reached only through nodes with
 * refcount 1.  A node under a shared parent may be reachable from outside, and
 * a rule list detached by Block or an iterator is invisible to the walk; both
 * leave the count unmatched and keep the symbol, which is the safe direction.
 *
 * Reclaiming frees the rule lists, unlinks the node and frees the SymbolDef.
 * The NAME string is kept (in a graveyard freed by intern_clear): a stale raw
 * name pointer can then never alias a later symbol's name.
 *
 * Timing.  Each local is tried at its own Module's exit -- safe, since no frame
 * outside that call can hold it unless it escaped, and then it is not
 * reclaimable.  One that escaped (returned, stored) is queued and retried when
 * the outermost evaluation ends, never mid-evaluation, where a frame might hold
 * a queued local's SymbolDef.
 * MATHILDA_NO_TEMP_RECLAIM=1 disables reclamation (A/B switch).
 * ============================================================ */

static int temp_reclaim_enabled(void) {
    static int on = -1;
    if (on < 0) {
        const char* v = getenv("MATHILDA_NO_TEMP_RECLAIM");
        on = (v && v[0] && v[0] != '0') ? 0 : 1;
    }
    return on;
}

static char** g_name_graveyard = NULL;
static size_t g_graveyard_n = 0, g_graveyard_cap = 0;

static SymbolDef** g_temp_pending = NULL;
static size_t g_temp_pending_n = 0, g_temp_pending_cap = 0;
static size_t g_temp_queued_since_sweep = 0;
unsigned long g_temp_dirty = 0;   /* nodes of queued locals freed since the last sweep */

bool symtab_name_interned(const char* name) {
    return node_find(name) != NULL;
}

bool symtab_temp_begin(SymbolDef* d) {
    if (!d || !temp_reclaim_enabled()) return false;
    d->module_temp = 1;
    d->temp_nodes = 0;
    d->temp_pending = 0;
    return true;
}

/* Nodes naming `name` in `e`, counted into *found when reached through
 * refcount-1 nodes only; false as soon as one sits under a shared node. */
static bool temp_walk(const Expr* e, const char* name, bool exclusive, uint32_t* found) {
    if (!e) return true;
    bool excl = exclusive && e->refcount == 1;
    if (e->type == EXPR_SYMBOL) {
        if (e->data.symbol.name != name) return true;
        if (!excl) return false;
        (*found)++;
        return true;
    }
    if (e->type != EXPR_FUNCTION) return true;
    if (!temp_walk(e->data.function.head, name, excl, found)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!temp_walk(e->data.function.args[i], name, excl, found)) return false;
    return true;
}

static bool temp_walk_rules(const Rule* r, const char* name, uint32_t* found) {
    for (; r; r = r->next) {
        if (!temp_walk(r->pattern, name, true, found)) return false;
        if (!temp_walk(r->replacement, name, true, found)) return false;
    }
    return true;
}

static void graveyard_keep(char* name) {
    if (g_graveyard_n == g_graveyard_cap) {
        size_t c = g_graveyard_cap ? 2 * g_graveyard_cap : 1024;
        char** a = realloc(g_name_graveyard, sizeof(char*) * c);
        if (!a) return;                 /* unreachable but harmless: name leaks */
        g_name_graveyard = a;
        g_graveyard_cap = c;
    }
    g_name_graveyard[g_graveyard_n++] = name;
}

/* Reclaim `d` if nothing outside its own definitions refers to it. */
static bool temp_try_reclaim(SymbolDef* d) {
    if (!d->module_temp) return false;
    uint32_t found = 0;
    const char* name = d->symbol_name;
    if (!temp_walk_rules(d->own_values, name, &found) ||
        !temp_walk_rules(d->down_values, name, &found) ||
        !temp_walk_rules(d->sub_values, name, &found) ||
#if UP_VALUES
        !temp_walk_rules(d->up_values, name, &found) ||
#endif
        !temp_walk_rules(d->n_values, name, &found) ||
        !temp_walk_rules(d->format_values, name, &found) ||
        !temp_walk(d->default_options, name, true, &found))
        return false;
    if (found != d->temp_nodes) return false;

    /* Nothing outside refers to it: drop the definitions (their nodes go to 0),
     * unlink and free the node; keep the name string. */
    reset_node_payload(d);
    size_t idx = d->name_hash & (symtab_size - 1);
    for (SymbolDef** pp = &symtab[idx]; *pp; pp = &(*pp)->next) {
        if (*pp == d) { *pp = d->next; break; }
    }
    if (symtab_nodes) symtab_nodes--;
    size_t slot = intern_ptr_slot(name);
    if (g_intern_ptr_cache[slot] == d) g_intern_ptr_cache[slot] = NULL;
    graveyard_keep(d->symbol_name);
    free(d);
    return true;
}

void symtab_temp_sweep(void) {
    size_t w = 0;
    for (size_t i = 0; i < g_temp_pending_n; i++) {
        SymbolDef* d = g_temp_pending[i];
        if (d->module_temp && (d->attributes & ATTR_TEMPORARY) && temp_try_reclaim(d))
            continue;
        if (!d->module_temp || !(d->attributes & ATTR_TEMPORARY)) {
            /* Temporary cleared by the user: keep it for good. */
            d->temp_pending = 0;
            continue;
        }
        g_temp_pending[w++] = d;
    }
    g_temp_pending_n = w;
    g_temp_queued_since_sweep = 0;
    g_temp_dirty = 0;
}

/* The top-level hook: sweep only when a local was queued, or a queued local
 * lost a reference, since the last sweep -- a long-lived escaped local is not
 * re-walked by statements that never touch it. */
void symtab_temp_sweep_if_queued(void) {
    if (g_temp_queued_since_sweep || g_temp_dirty) symtab_temp_sweep();
}

void symtab_temp_release(SymbolDef* d) {
    if (!d || !d->module_temp || !temp_reclaim_enabled()) return;
    if ((d->attributes & ATTR_TEMPORARY) && temp_try_reclaim(d)) return;
    if (!d->temp_pending) {
        if (g_temp_pending_n == g_temp_pending_cap) {
            size_t c = g_temp_pending_cap ? 2 * g_temp_pending_cap : 256;
            SymbolDef** a = realloc(g_temp_pending, sizeof(SymbolDef*) * c);
            if (!a) return;             /* not queued: kept, as before */
            g_temp_pending = a;
            g_temp_pending_cap = c;
        }
        g_temp_pending[g_temp_pending_n++] = d;
        d->temp_pending = 1;
        g_temp_queued_since_sweep++;
    }
    /* No sweep here: an inner frame up the stack may hold another pending
     * local's SymbolDef (a rule loop); sweeps run only at the top level. */
}

int intern_is_system(const char* name) {
    SymbolDef* d = node_find(name);
    return d ? d->is_system == 1 : 0;
}

int intern_is_mathilda(const char* name) {
    SymbolDef* d = node_find(name);
    return d ? d->is_system == 2 : 0;
}

void intern_mark_mathilda(const char* name) {
    SymbolDef* d = node_find(name);
    if (d) d->is_system = 2;
}

void intern_mark_all_system(void) {
    for (size_t i = 0; i < symtab_size; i++)
        for (SymbolDef* d = symtab[i]; d; d = d->next)
            d->is_system = 1;
}

void intern_clear(void) {
    memset(g_intern_ptr_cache, 0, sizeof g_intern_ptr_cache);
    free(g_temp_pending); g_temp_pending = NULL;
    g_temp_pending_n = g_temp_pending_cap = 0;
    for (size_t i = 0; i < g_graveyard_n; i++) free(g_name_graveyard[i]);
    free(g_name_graveyard); g_name_graveyard = NULL;
    g_graveyard_n = g_graveyard_cap = 0;
    /* Full teardown: free every node, its payload, and its name string. Any
     * previously-returned interned pointer (including SYM_*) becomes invalid. */
    for (size_t i = 0; i < symtab_size; i++) {
        SymbolDef* d = symtab[i];
        while (d) {
            SymbolDef* next = d->next;
            reset_node_payload(d);
            free(d->symbol_name);
            free(d);
            d = next;
        }
        symtab[i] = NULL;
    }
    symtab_nodes = 0;
}

void symtab_add_builtin(const char* symbol_name, BuiltinFunc func) {
    SymbolDef* def = symtab_get_def(symbol_name);
    def->builtin_func = func;
}

void symtab_set_docstring(const char* symbol_name, const char* doc) {
    SymbolDef* def = symtab_get_def(symbol_name);
    if (def->docstring) free(def->docstring);
    def->docstring = doc ? mathilda_strdup(doc) : NULL;
}

const char* symtab_get_docstring(const char* symbol_name) {
    SymbolDef* def = symtab_get_def(symbol_name);
    return def->docstring;
}

void symtab_set_ndarray_unary_kernel(const char* symbol_name, const void* kernel) {
    SymbolDef* def = symtab_get_def(symbol_name);
    def->ndarray_unary_kernel = kernel;
    /* A registered kernel reads the buffer directly and degrades to the List
     * path when it declines, so the head is packed-aware by construction. */
    def->packed_aware = 1;
}

void symtab_set_ndarray_binary_kernel(const char* symbol_name, const void* kernel) {
    SymbolDef* def = symtab_get_def(symbol_name);
    def->ndarray_binary_kernel = kernel;
    def->packed_aware = 1;
}

void symtab_set_packed_aware(const char* symbol_name) {
    symtab_get_def(symbol_name)->packed_aware = 1;
}

/* See SymbolDef.sparse_mode. */
void symtab_set_sparse_mode(const char* symbol_name, int mode) {
    symtab_get_def(symbol_name)->sparse_mode = (uint8_t)mode;
}

/* See SymbolDef.sparse_native. */
void symtab_set_sparse_native(const char* symbol_name, const void* fn) {
    symtab_get_def(symbol_name)->sparse_native = fn;
}

/* See SymbolDef.packed_broadcast_ok. Implies packed_aware: a head cannot claim
 * to handle a rank mismatch on a buffer without handling the buffer at all. */
void symtab_set_packed_broadcast_ok(const char* symbol_name) {
    SymbolDef* def = symtab_get_def(symbol_name);
    def->packed_broadcast_ok = 1;
    def->packed_aware = 1;
}

/* Undo packed_aware, for a head that acquired it implicitly from a kernel
 * registration but whose kernel answers with different element HEADS than the
 * ordinary List does (Floor of a Real is an Integer; Clip clamps to its exact
 * bounds). See the NOT_AWARE list in src/pack.c for why each one is there. */
void symtab_clear_packed_aware(const char* symbol_name) {
    SymbolDef* def = symtab_get_def(symbol_name);
    def->packed_aware = 0;
    def->packed_int64_ok = 0;
}

void symtab_set_packed_int64_ok(const char* symbol_name) {
    SymbolDef* def = symtab_get_def(symbol_name);
    def->packed_aware = 1;
    def->packed_int64_ok = 1;
}

void symtab_set_ndarray_nary_kernel(const char* symbol_name, const void* kernel) {
    SymbolDef* def = symtab_get_def(symbol_name);
    def->ndarray_nary_kernel = kernel;
    /* A Listable n-ary head maps over a buffer operand (ndarray_map_nary), so it
     * takes its packed arguments as buffers, like the unary / binary heads. */
    if (def->attributes & ATTR_LISTABLE) def->packed_aware = 1;
}

void symtab_set_options(const char* symbol_name, Expr* list) {
    SymbolDef* def = symtab_get_def(symbol_name);
    if (def->default_options) expr_free(def->default_options);
    def->default_options = list;
}

Expr* symtab_get_options(const char* symbol_name) {
    SymbolDef* def = symtab_lookup(symbol_name);
    return def ? def->default_options : NULL;
}

/* ============================================================
 * §3.5/§3.6: dispatch-key extraction and specificity scoring.
 *
 * These run once per rule at insertion time. The output is a cheap
 * pre-computed key (dispatch_arity, first_arg_head_canon) and an
 * integer specificity score that lets apply_down_values short-circuit
 * unmatchable rules and the rule list stay sorted "most specific first."
 * ============================================================ */

/* Strip transparent pattern wrappers down to the inner pattern.
 *
 * - HoldPattern / Verbatim are transparent for matching, so the
 *   dispatch key should look through them.
 * - Condition[p, c] and PatternTest[p, t] are guards that gate the
 *   match without altering its STRUCTURAL shape — the dispatch
 *   filter wants the arity and first-arg-head of the inner pattern
 *   `p`, not of the Condition/PatternTest wrapper itself.  The
 *   matcher (match.c:278-313 for Condition, 394-417 for PatternTest)
 *   evaluates the guard post-match; we just need dispatch not to
 *   pre-filter the rule away on a confused arity/head key.
 *
 * Stripping Condition here was the missing piece that prevented
 * every `f[...] /; cond := rhs` DownValue from firing (pre-fix the
 * rule's dispatch_arity was the Condition's arg count and its
 * first_arg_head_canon was the inner pattern's HEAD symbol — the
 * filter then routinely disagreed with an input's actual first-arg
 * head and skipped the rule entirely).  See tests/test_condition_
 * downvalue.c for the pinning suite. */
static const Expr* strip_pattern_wrappers(const Expr* p) {
    while (p && p->type == EXPR_FUNCTION && p->data.function.head &&
           p->data.function.head->type == EXPR_SYMBOL &&
           p->data.function.arg_count >= 1) {
        const char* h = p->data.function.head->data.symbol.name;
        if (h == SYM_HoldPattern || h == SYM_Verbatim) {
            p = p->data.function.args[0];
            continue;
        }
        if ((h == SYM_Condition || h == SYM_PatternTest) &&
            p->data.function.arg_count == 2) {
            p = p->data.function.args[0];
            continue;
        }
        break;
    }
    return p;
}

/* True when `p`, considered as a top-level argument slot in a DownValue
 * pattern, can absorb a variable number of input arguments. Patterns
 * that satisfy this force dispatch_arity = -1. */
static bool slot_is_variable_length(const Expr* p) {
    p = strip_pattern_wrappers(p);
    if (!p) return false;
    if (p->type != EXPR_FUNCTION || !p->data.function.head ||
        p->data.function.head->type != EXPR_SYMBOL) {
        return false;
    }
    const char* h = p->data.function.head->data.symbol.name;
    if (h == SYM_BlankSequence) return true;
    if (h == SYM_BlankNullSequence) return true;
    if (h == SYM_Optional) return true;
    if (h == SYM_Repeated) return true;
    if (h == SYM_RepeatedNull) return true;
    if (h == SYM_OptionsPattern) return true;
    /* A pattern sequence spans as many arguments as its elements do -- zero,
     * one or several -- so it cannot pin the arity (q[PatternSequence[a_, b_],
     * c_] takes three arguments, and was skipped for q[1, 2, 3]). */
    if (h == SYM_PatternSequence) return true;
    if (h == SYM_OrderlessPatternSequence) return true;
    if ((h == SYM_Longest || h == SYM_Shortest) && p->data.function.arg_count >= 1)
        return slot_is_variable_length(p->data.function.args[0]);
    /* Alternatives are variable iff some branch is (x__ | y_). */
    if (h == SYM_Alternatives) {
        for (size_t i = 0; i < p->data.function.arg_count; i++)
            if (slot_is_variable_length(p->data.function.args[i])) return true;
        return false;
    }
    /* Pattern[x, q] is variable iff q is. */
    if (h == SYM_Pattern && p->data.function.arg_count == 2) {
        return slot_is_variable_length(p->data.function.args[1]);
    }
    return false;
}

/* True when `p` is an Optional pattern (`x_.`, `Optional[p]`,
 * `Optional[p, def]`), possibly under a `Pattern[name, ...]` wrapper.
 * An Optional argument in a function-call pattern can be omitted, which
 * collapses the call to its base form (`Power[b_, n_.]` matches a bare
 * `b` as `b^1`; `Times[a_., x_]` matches `x`), so the concrete value
 * need not share the pattern's head. */
static bool pattern_arg_is_optional(const Expr* p) {
    if (!p || p->type != EXPR_FUNCTION || !p->data.function.head ||
        p->data.function.head->type != EXPR_SYMBOL) {
        return false;
    }
    const char* h = p->data.function.head->data.symbol.name;
    if (h == SYM_Optional) return true;
    if (h == SYM_Pattern && p->data.function.arg_count == 2) {
        return pattern_arg_is_optional(p->data.function.args[1]);
    }
    return false;
}

/* Extract the canonical head symbol that any concrete value matching
 * this pattern slot must have. Returns NULL when the slot is wildcard
 * with respect to head (Blank[], _, anything sequence-like, or anything
 * we are not willing to specialize on). */
static const char* pattern_arg_head_canon(const Expr* p) {
    p = strip_pattern_wrappers(p);
    if (!p) return NULL;

    /* Atoms in patterns match themselves -- their head equals the head
     * of the corresponding atomic input value. */
    if (p->type == EXPR_INTEGER || p->type == EXPR_BIGINT) {
        return intern_symbol("Integer");
    }
    if (p->type == EXPR_REAL) return intern_symbol("Real");
    if (p->type == EXPR_STRING) return intern_symbol("String");
    if (p->type == EXPR_SYMBOL) {
        /* Bare symbol used as a literal (e.g. f[True, x_]). Head is
         * "Symbol" -- the filter will only reject inputs whose first arg
         * is not a symbol (e.g. an integer or a function call). */
        return intern_symbol("Symbol");
    }
    if (p->type != EXPR_FUNCTION) return NULL;

    Expr* head = p->data.function.head;
    if (!head || head->type != EXPR_SYMBOL) return NULL;
    const char* h = head->data.symbol.name;

    /* Blank[] -> wildcard; Blank[h] -> h. */
    if (h == SYM_Blank) {
        if (p->data.function.arg_count == 0) return NULL;
        Expr* head_arg = p->data.function.args[0];
        if (head_arg && head_arg->type == EXPR_SYMBOL) {
            return intern_symbol(head_arg->data.symbol.name);
        }
        return NULL;
    }
    /* Pattern[x, q] -> use q's head. */
    if (h == SYM_Pattern && p->data.function.arg_count == 2) {
        return pattern_arg_head_canon(p->data.function.args[1]);
    }
    /* Sequence-like or optional patterns are wildcards w.r.t. head. */
    if (h == SYM_BlankSequence || h == SYM_BlankNullSequence ||
        h == SYM_Optional || h == SYM_Repeated ||
        h == SYM_RepeatedNull || h == SYM_OptionsPattern) {
        /* For typed BlankSequence[h] etc., picking up h would still be
         * unsafe because the matcher binds whole sequences. Stay safe. */
        return NULL;
    }
    /* Condition[p, c] / PatternTest[p, t] are transparent for head. */
    if ((h == SYM_Condition || h == SYM_PatternTest) &&
        p->data.function.arg_count >= 1) {
        return pattern_arg_head_canon(p->data.function.args[0]);
    }
    /* Alternatives[a, b, ...] -- only safe to specialize when every
     * branch agrees on a head. */
    if (h == SYM_Alternatives) {
        const char* acc = NULL;
        for (size_t i = 0; i < p->data.function.arg_count; i++) {
            const char* hi = pattern_arg_head_canon(p->data.function.args[i]);
            if (!hi) return NULL;          /* one branch is wildcard */
            if (acc && acc != hi) return NULL; /* branches disagree */
            acc = hi;
        }
        return acc;
    }
    /* Except[...] matches by *excluding* a pattern, and KeyValuePattern[...]
     * matches an association/list by its contents -- in neither case is the
     * head a literal that the input's first-arg head must equal. Treat them as
     * wildcards so the §3.5 dispatch filter never skips such rules. */
    if (h == SYM_Except || h == SYM_KeyValuePattern) {
        return NULL;
    }
    /* Plain function call in a pattern -- head must match, UNLESS one of
     * its arguments is Optional. An omitted Optional collapses the call to
     * its base (e.g. Power[b_, n_.] matches a bare `b` as b^1, Times[a_., x_]
     * matches `x`), so a concrete value need not share this head. Returning
     * the head here would make the §3.5 dispatch filter skip such rules for
     * any input that is not literally headed by `h`. */
    for (size_t i = 0; i < p->data.function.arg_count; i++) {
        if (pattern_arg_is_optional(p->data.function.args[i])) return NULL;
    }
    return intern_symbol(h);
}

/* §3.5 dispatch_arity: count of top-level arg slots, or -1 when any
 * slot is variable-length. Pattern is the WHOLE LHS, e.g. f[a, b]; we
 * inspect the top-level function's args. */
static int32_t pattern_dispatch_arity(const Expr* pattern) {
    pattern = strip_pattern_wrappers(pattern);
    if (!pattern || pattern->type != EXPR_FUNCTION) return -1;
    size_t n = pattern->data.function.arg_count;
    for (size_t i = 0; i < n; i++) {
        if (slot_is_variable_length(pattern->data.function.args[i])) return -1;
    }
    return (int32_t)n;
}

/* §3.5 first_arg_head_canon: NULL when arity is variable (a leading
 * sequence pattern would absorb the first arg) OR the first slot is
 * wildcard. */
static const char* pattern_first_arg_head(const Expr* pattern) {
    pattern = strip_pattern_wrappers(pattern);
    if (!pattern || pattern->type != EXPR_FUNCTION) return NULL;
    if (pattern->data.function.arg_count == 0) return NULL;
    /* If any slot at or before the first is variable-length, the "first
     * arg" of the input might end up bound to the sequence. Be safe. */
    if (slot_is_variable_length(pattern->data.function.args[0])) return NULL;
    return pattern_arg_head_canon(pattern->data.function.args[0]);
}

/* §3.6 specificity score. Computed recursively over the pattern tree.
 *
 *   literal atoms                  : +100 (most discriminating)
 *   typed Blank[h]                 :  +20
 *   plain Blank[]                  :  -10
 *   BlankSequence[]                : -100
 *   BlankNullSequence[]            : -200
 *   Optional[...]                  :  -50
 *   Repeated[p]                    : score(p) - 50
 *   RepeatedNull[p]                : score(p) - 100
 *   Condition / PatternTest        : score(inner) + 10 (extra constraint)
 *   HoldPattern / Verbatim         : score(inner) (transparent)
 *   Alternatives[a, b, ...]        : min(score(a_i)) (weakest branch wins)
 *   plain function call f[a, b...] :  +50 + sum(score(arg_i))
 *   Pattern[name, q]               : score(q) (named binding adds nothing)
 *
 * Scores are added across siblings, so a 3-arg literal call beats a
 * 2-arg literal call. Insertion order is used as the tie-break, so two
 * structurally-similar pattern rules keep the user's chosen ordering. */
static int32_t pattern_specificity(const Expr* p) {
    if (!p) return 0;
    if (p->type == EXPR_INTEGER || p->type == EXPR_BIGINT ||
        p->type == EXPR_REAL || p->type == EXPR_STRING) {
        return 100;
    }
    if (p->type == EXPR_SYMBOL) {
        return 100;
    }
    if (p->type != EXPR_FUNCTION) return 0;

    Expr* head = p->data.function.head;
    if (head && head->type == EXPR_SYMBOL) {
        const char* h = head->data.symbol.name;
        size_t ac = p->data.function.arg_count;

        if (h == SYM_HoldPattern || h == SYM_Verbatim) {
            return ac >= 1 ? pattern_specificity(p->data.function.args[0]) : 0;
        }
        if (h == SYM_Pattern) {
            /* Pattern[name, q] -- the name is irrelevant for matching power. */
            return ac >= 2 ? pattern_specificity(p->data.function.args[1]) : 0;
        }
        if (h == SYM_Blank) {
            if (ac == 0) return -10;
            return 20; /* typed blank */
        }
        if (h == SYM_BlankSequence) return -100;
        if (h == SYM_BlankNullSequence) return -200;
        if (h == SYM_Optional) return -50;
        if (h == SYM_Repeated) {
            return (ac >= 1 ? pattern_specificity(p->data.function.args[0]) : 0) - 50;
        }
        if (h == SYM_RepeatedNull) {
            return (ac >= 1 ? pattern_specificity(p->data.function.args[0]) : 0) - 100;
        }
        /* A condition does not make a rule more specific: the reference keeps
         * f[p[a_, a_]] := ... ahead of a later f[p[a1_, b1_]] := ... /; c, and a
         * rule written rhs /; c1 /; c2 is not promoted past ones with a single
         * condition (FeynCalc's fromGFAD rules depend on that order). Only
         * between two rules that differ in nothing but the condition does the
         * conditioned one go first -- see the insertion in add_rule. */
        if (h == SYM_Condition || h == SYM_PatternTest) {   /* x_?test likewise */
            return ac >= 1 ? pattern_specificity(p->data.function.args[0]) : 0;
        }
        if (h == SYM_Alternatives) {
            if (ac == 0) return 0;
            int32_t best = pattern_specificity(p->data.function.args[0]);
            for (size_t i = 1; i < ac; i++) {
                int32_t s = pattern_specificity(p->data.function.args[i]);
                if (s < best) best = s; /* weakest branch wins */
            }
            return best;
        }
        if (h == SYM_OptionsPattern) return -150;
    }

    /* Plain function call: literal head + sum of children. */
    int32_t score = 50;
    if (head) score += pattern_specificity(head);
    for (size_t i = 0; i < p->data.function.arg_count; i++) {
        score += pattern_specificity(p->data.function.args[i]);
    }
    /* A literal-headed call h[...] matches only expressions with head h, so it
     * is at least as specific as the typed blank _h (20) -- and so as any plain
     * blank -- however loose its arguments. Without this floor the sequence
     * penalties of {f___} (-50) sank it below a sibling f_ (-10), and
     * c[e_, f_, o:OptionsPattern[]] was tried before c[e_, {f___}, o:...],
     * binding f to the whole List. */
    if (head && head->type == EXPR_SYMBOL && score < 21) score = 21;
    return score;
}

/* ============================================================
 * §3.4 pattern canonicalization
 *
 * Stored DownValue / OwnValue patterns must share canonical form with
 * the runtime expressions they will eventually be matched against;
 * otherwise the §3.5 dispatch filter -- and in some cases the matcher
 * itself -- can silently reject a rule that the abstract pattern would
 * accept (i.e. one whose MatchQ form returns True against the input).
 *
 * The asymmetry exists because `Set` and `SetDelayed` deliberately hold
 * the LHS: every runtime input flowing into `apply_down_values` has
 * already been canonicalized by `evaluate_step` -- Flat heads flattened,
 * Orderless heads sorted, `Power[Times[a,...,z], n_Integer]` distributed
 * over the factors (src/power.c integer-exponent fast path), nested
 * `Power[Power[base, e1], e2]` collapsed -- but the stored pattern was
 * never put through that pass. The pre-2026-05 dispatch filter then
 * compared the parsed-but-uncanonicalized first-arg head of the rule
 * (e.g. `Power`) against the evaluated first-arg head of the input
 * (e.g. `Times`), found them unequal, and skipped the rule before the
 * matcher ran.
 *
 * The concrete bug that motivated this pass was the CRC integral table:
 * rules like
 *
 *     IntegrateTable[1/(x_ Sqrt[a_ + b_. x_ + c_. x_^2]), x_] /; ...
 *
 * parse as `Power[Times[x_, Sqrt[...]], -1]` at the first-arg slot, but
 * every runtime call arrives as `Times[Power[x,-1], Power[Sqrt[...],-1]]`
 * because Power-over-Times-integer distribution runs during normal
 * evaluation. The dispatch filter cached `Power` as the rule's first-arg
 * head, the input presented `Times`, and ~80 CRC formulas with similar
 * `1/(...)` LHSs went silently unmatched (see tasks/crc_corpus_2026-05-15.md
 * "Out-of-scope findings" -- which misattributed the failure to a 3-term
 * Plus matcher gap).
 *
 * What this pass structurally changes in the pattern matcher:
 *
 *   - Stored patterns now share canonical form with runtime inputs.
 *     `Power[Times[t1,...,tn], k]` and `Times[Power[t1,k],...,Power[tn,k]]`
 *     are no longer two distinguishable LHSs when `k` is an integer
 *     literal; the second form is the only one that survives storage.
 *   - Likewise, `Power[Power[b, e1], e2]` collapses to `Power[b, e1*e2]`
 *     when both exponents are numeric literals, matching the evaluator's
 *     nested-Power normalization.
 *   - Flat-headed pattern children (Plus, Times, ...) are flattened, and
 *     Orderless-headed pattern children are sorted with
 *     `eval_compare_expr_ptrs` -- the same comparator that orders input
 *     arguments -- so the §3.5 dispatch filter sees consistent shapes on
 *     both sides.
 *   - Pattern-marker wrappers (Pattern, Blank, BlankSequence,
 *     BlankNullSequence, Optional, Condition, PatternTest, HoldPattern,
 *     Verbatim, Alternatives, Repeated, RepeatedNull, OptionsPattern)
 *     are NOT structurally rewritten; the canonicalizer recurses into the
 *     inner pattern slot of each but never folds the wrapper itself. This
 *     preserves match-time semantics in full (a `/; cond` guard still
 *     evaluates with bindings in scope, an Optional default still binds
 *     when the slot is absent, HoldPattern still forces literal matching,
 *     etc.).
 *
 * The pass is idempotent: applying it twice to a pattern yields the same
 * tree as applying it once, so re-storing an already-canonicalized rule
 * is safe.
 *
 * Callers transfer ownership of `p` to `pattern_canonicalize`; the
 * returned tree may share or replace `p`. ============================================================ */
static Expr* pattern_canonicalize(Expr* p);

/* True for an exponent we are willing to fold across nested Powers.
 * Restricted to exact integer / rational literals so the collapse is
 * always sound -- floating-point or symbolic exponents are left alone. */
static bool is_foldable_numeric_exponent(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    int64_t n, d;
    return is_rational(e, &n, &d);
}

/* Returns Times[copies-of-a-and-b] reduced through `evaluate` so the
 * result is in canonical integer/rational form. Used to compose nested
 * Power exponents; both inputs must satisfy `is_foldable_numeric_exponent`. */
static Expr* compose_numeric_exponents(const Expr* a, const Expr* b) {
    Expr* args[2] = { expr_copy((Expr*)a), expr_copy((Expr*)b) };
    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times), args, 2);
    Expr* res = evaluate(prod);
    /* evaluate() makes its own copy of the input — release ours. */
    expr_free(prod);
    return res;
}

static Expr* pattern_canonicalize(Expr* p) {
    if (!p || p->type != EXPR_FUNCTION) return p;

    Expr* head = p->data.function.head;
    if (head && head->type == EXPR_SYMBOL) {
        const char* h = head->data.symbol.name;

        /* HoldPattern / Verbatim explicitly opt out of canonicalization;
         * their whole purpose is to match a literal subtree without
         * pattern-semantic interpretation. */
        if (h == SYM_HoldPattern || h == SYM_Verbatim) return p;

        /* Blank-family and OptionsPattern carry only type/head metadata
         * in their children -- no nested pattern that wants canonicalization. */
        if (h == SYM_Blank || h == SYM_BlankSequence ||
            h == SYM_BlankNullSequence || h == SYM_OptionsPattern) {
            return p;
        }

        /* Pattern[name, q]: name is held literally; canonicalize q only. */
        if (h == SYM_Pattern) {
            if (p->data.function.arg_count >= 2) {
                p = expr_unshare(p);
                p->data.function.args[1] =
                    pattern_canonicalize(p->data.function.args[1]);
            }
            return p;
        }

        /* Optional[p] / Optional[p, default]: canonicalize the inner
         * pattern. The default is a value, not a pattern; leave it alone
         * to avoid firing evaluator side effects before bindings exist. */
        if (h == SYM_Optional) {
            if (p->data.function.arg_count >= 1) {
                p = expr_unshare(p);
                p->data.function.args[0] =
                    pattern_canonicalize(p->data.function.args[0]);
            }
            return p;
        }

        /* Condition[p, guard] / PatternTest[p, test]: canonicalize the
         * inner pattern. The guard / test runs at match time with the
         * bindings substituted -- canonicalizing it here would force
         * evaluation against still-unbound variables. */
        if (h == SYM_Condition || h == SYM_PatternTest) {
            if (p->data.function.arg_count >= 1) {
                p = expr_unshare(p);
                p->data.function.args[0] =
                    pattern_canonicalize(p->data.function.args[0]);
            }
            return p;
        }

        /* Repeated[p] / RepeatedNull[p]: canonicalize the inner pattern. */
        if (h == SYM_Repeated || h == SYM_RepeatedNull) {
            if (p->data.function.arg_count >= 1) {
                p = expr_unshare(p);
                p->data.function.args[0] =
                    pattern_canonicalize(p->data.function.args[0]);
            }
            return p;
        }

        /* Alternatives[a, b, ...]: each branch is itself a pattern. */
        if (h == SYM_Alternatives) {
            p = expr_unshare(p);
            for (size_t i = 0; i < p->data.function.arg_count; i++) {
                p->data.function.args[i] =
                    pattern_canonicalize(p->data.function.args[i]);
            }
            return p;
        }
    }

    /* Plain function call: canonicalize all children, then apply the
     * structural rewrites that mirror evaluate_step. */
    p = expr_unshare(p);
    if (p->data.function.head) {
        p->data.function.head = pattern_canonicalize(p->data.function.head);
    }
    for (size_t i = 0; i < p->data.function.arg_count; i++) {
        p->data.function.args[i] =
            pattern_canonicalize(p->data.function.args[i]);
    }

    /* Look up structural attributes by head symbol. */
    head = p->data.function.head;
    uint32_t attrs = 0;
    if (head && head->type == EXPR_SYMBOL) {
        SymbolDef* def = symtab_lookup(head->data.symbol.name);
        if (def) attrs = def->attributes;
    }

    /* Flat: flatten nested same-head args (Plus[Plus[a,b],c] -> Plus[a,b,c]). */
    if (head && head->type == EXPR_SYMBOL && (attrs & ATTR_FLAT)) {
        eval_flatten_args(p, head->data.symbol.name);
    }

    /* Orderless: sort args using the same comparator the evaluator uses
     * for input expressions, so the §3.5 first-arg-head dispatch key and
     * the matcher both see identical orderings. */
    if (head && head->type == EXPR_SYMBOL && (attrs & ATTR_ORDERLESS)) {
        eval_sort_args(p);
    }

    /* Power[Times[t1,...,tn], k_Integer] -> Times[Power[t1,k],...,Power[tn,k]].
     * Mirrors src/power.c's integer-exponent fast path that fires
     * unconditionally during normal evaluation. */
    if (head && head->type == EXPR_SYMBOL && head->data.symbol.name == SYM_Power
        && p->data.function.arg_count == 2) {
        Expr* base = p->data.function.args[0];
        Expr* exp  = p->data.function.args[1];
        if (exp && (exp->type == EXPR_INTEGER || exp->type == EXPR_BIGINT)
            && base && base->type == EXPR_FUNCTION
            && base->data.function.head
            && base->data.function.head->type == EXPR_SYMBOL
            && base->data.function.head->data.symbol.name == SYM_Times) {
            size_t bc = base->data.function.arg_count;
            Expr** new_args = malloc(sizeof(Expr*) * bc);
            for (size_t i = 0; i < bc; i++) {
                Expr* pow_args[2] = {
                    expr_copy(base->data.function.args[i]),
                    expr_copy(exp)
                };
                new_args[i] = expr_new_function(
                    expr_new_symbol(SYM_Power), pow_args, 2);
            }
            Expr* times = expr_new_function(
                expr_new_symbol(SYM_Times), new_args, bc);
            free(new_args);
            expr_free(p);
            /* Recurse on the new Times so its children (each a Power
             * we just synthesized) themselves get folded if applicable
             * (e.g. Power[Power[plus, 1/2], -1] -> Power[plus, -1/2]). */
            return pattern_canonicalize(times);
        }
    }

    /* Power[Power[base, e1], e2] -> Power[base, e1*e2] when both
     * exponents are numeric literals. Matches the evaluator's nested-
     * Power collapse for sound bases. */
    if (head && head->type == EXPR_SYMBOL && head->data.symbol.name == SYM_Power
        && p->data.function.arg_count == 2) {
        Expr* base = p->data.function.args[0];
        Expr* outer_exp = p->data.function.args[1];
        if (base && base->type == EXPR_FUNCTION
            && base->data.function.head
            && base->data.function.head->type == EXPR_SYMBOL
            && base->data.function.head->data.symbol.name == SYM_Power
            && base->data.function.arg_count == 2) {
            Expr* inner_base = base->data.function.args[0];
            Expr* inner_exp  = base->data.function.args[1];
            if (is_foldable_numeric_exponent(inner_exp)
                && is_foldable_numeric_exponent(outer_exp)) {
                Expr* new_exp = compose_numeric_exponents(inner_exp, outer_exp);
                Expr* pow_args[2] = { expr_copy(inner_base), new_exp };
                Expr* folded = expr_new_function(
                    expr_new_symbol(SYM_Power), pow_args, 2);
                expr_free(p);
                /* The collapsed Power may now expose a further Power-of-
                 * Times opportunity if `inner_base` itself is a Times;
                 * recurse to let the previous branch fire. */
                return pattern_canonicalize(folded);
            }
        }
    }

    return p;
}

// Helper to add a rule to a list
static Expr* pattern_alpha_normalize(const Expr* p);
static uint64_t pattern_alpha_hash(const Expr* p);

/* ---- Conditions at the end of a scoping construct (`lhs := Module[{..}, v /; c]`) */

static bool head_sym_is(const Expr* e, const char* sym, size_t argc) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym && e->data.function.arg_count == argc;
}

/* The tail position of a scoping construct: Module/With/Block[vars, body] ->
 * body; CompoundExpression[..., last] -> last.  Returns the address of the
 * slot holding the tail of `*slot`, or NULL when `*slot` is neither. */
static Expr** scoped_tail_slot(Expr** slot) {
    Expr* e = *slot;
    if (head_sym_is(e, SYM_Module, 2) || head_sym_is(e, SYM_With, 2) || head_sym_is(e, SYM_Block, 2))
        return &e->data.function.args[1];
    if (e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_CompoundExpression && e->data.function.arg_count > 0)
        return &e->data.function.args[e->data.function.arg_count - 1];
    return NULL;
}

/* Does `rhs` -- a Module / With / Block, possibly nested and through a
 * CompoundExpression's last statement -- end in Condition[val, cond]? */
static bool rhs_has_scoped_condition(Expr* rhs) {
    if (!(head_sym_is(rhs, SYM_Module, 2) || head_sym_is(rhs, SYM_With, 2) || head_sym_is(rhs, SYM_Block, 2)))
        return false;
    Expr* cur = rhs;
    Expr** slot;
    while ((slot = scoped_tail_slot(&cur)) != NULL) cur = *slot;
    return head_sym_is(cur, SYM_Condition, 2);
}

/* An owned copy of the bound replacement with its tail Condition[val, cond]
 * turned into If[TrueQ[cond], val, $ScopedConditionFail$], so evaluating it
 * either yields the rule's value or the failure sentinel. */
static Expr* scoped_condition_rewrite(Expr* rhs) {
    Expr* root = expr_unshare(expr_copy(rhs));
    Expr** slot = &root;
    Expr** next;
    while ((next = scoped_tail_slot(slot)) != NULL) {
        *next = expr_unshare(*next);
        slot = next;
    }
    Expr* cond = *slot;
    if (!head_sym_is(cond, SYM_Condition, 2)) return root;
    Expr* tq = expr_new_function(expr_new_symbol("TrueQ"),
                                 (Expr*[]){ expr_copy(cond->data.function.args[1]) }, 1);
    Expr* ifargs[3] = { tq, expr_copy(cond->data.function.args[0]),
                        expr_new_symbol("$ScopedConditionFail$") };
    Expr* iff = expr_new_function(expr_new_symbol(SYM_If), ifargs, 3);
    expr_free(cond);
    *slot = iff;
    return root;
}

/* Try a conditional rule's bound replacement: its evaluated value, or NULL
 * when the scoped guard did not give True (the rule does not apply).  An
 * in-flight Throw/Abort is returned as the value. */
static Expr* try_scoped_condition(Expr* bound) {
    Expr* rw = scoped_condition_rewrite(bound);
    Expr* v = evaluate(rw);
    expr_free(rw);
    if (v && v->type == EXPR_SYMBOL && strcmp(v->data.symbol.name, "$ScopedConditionFail$") == 0) {
        expr_free(v);
        return NULL;
    }
    return v;
}

bool  symtab_rhs_has_scoped_condition(Expr* rhs) { return rhs_has_scoped_condition(rhs); }
Expr* symtab_try_scoped_condition(Expr* bound)   { return try_scoped_condition(bound); }

/* Does e mention the symbol Return anywhere (heads included)? */
static bool mentions_return(const Expr* e) {
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == SYM_Return;
    if (e->type != EXPR_FUNCTION) return false;
    if (mentions_return(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (mentions_return(e->data.function.args[i])) return true;
    return false;
}

/* The value of applying `rule` of `head` when its rewrite is `result`: for a
 * rule whose right-hand side mentions Return, the result is evaluated here and
 * a top-level Return[v] (or Return[v, head]) yields v, Return[] Null -- the
 * reference's return from a definition. Return[v, other] and a Return below
 * the top level (1 + Return[2]) are left as they are. Consumes result. */
static Expr* rule_result(const Rule* rule, Expr* result, const char* head) {
    if (!rule->has_return || !result) return result;
    Expr* ev = evaluate(result);
    expr_free(result);
    Expr* v = NULL;
    if (eval_classify_return(ev, head, &v) == EVAL_RETURN_CONSUME) {
        expr_free(ev);
        return v;
    }
    return ev;
}

/* ---- Rule order: pattern generality --------------------------------------
 *
 * The reference tries a new definition ahead of an existing one only when the
 * existing pattern is more general -- matches everything the new one matches
 * (f[x_] after f[2], f[x_Integer] or f[x_ /; x > 0]); two patterns neither of
 * which covers the other keep the order they were defined in, however
 * "specific" either looks (FeynCalc's h[b___G, di_G, c__G] /; .. stays ahead of
 * a later h[b___, G[5], c:G[_[__]].., d___]).  pattern_subsumes(g, s) answers
 * "does g match every expression s matches?" structurally and conservatively:
 * false whenever it cannot tell, which only ever keeps definition order. */

typedef enum { PK_ONE, PK_OPT, PK_SEQ0, PK_SEQ1, PK_OTHER } PatKind;

static bool pat_is(const Expr* e, const char* h) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == h;
}

/* Pattern names and HoldPattern change nothing about what matches. */
static const Expr* pat_unname(const Expr* p) {
    for (;;) {
        if (pat_is(p, SYM_Pattern) && p->data.function.arg_count == 2)
            p = p->data.function.args[1];
        else if (pat_is(p, SYM_HoldPattern) && p->data.function.arg_count == 1)
            p = p->data.function.args[0];
        else
            return p;
    }
}

/* The restrictions on an s-side element (a test or condition only narrows
 * what s matches, so g need not reproduce it). */
static const Expr* pat_untested(const Expr* p) {
    for (;;) {
        p = pat_unname(p);
        if ((pat_is(p, SYM_Condition) || pat_is(p, SYM_PatternTest)) &&
            p->data.function.arg_count == 2)
            p = p->data.function.args[0];
        else
            return p;
    }
}

static bool pat_is_construct(const Expr* e) {
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    return h == SYM_Blank || h == SYM_BlankSequence || h == SYM_BlankNullSequence ||
           h == SYM_Pattern || h == SYM_Optional || h == SYM_Alternatives ||
           h == SYM_Repeated || h == SYM_RepeatedNull || h == SYM_Condition ||
           h == SYM_PatternTest || h == SYM_HoldPattern || h == SYM_Verbatim ||
           h == SYM_Except || h == SYM_OptionsPattern || h == SYM_PatternSequence ||
           h == SYM_Longest || h == SYM_Shortest || h == SYM_OrderlessPatternSequence ||
           h == SYM_KeyValuePattern;
}

/* How many arguments an argument-position element takes, and (for a sequence)
 * the pattern each of them must match: *ep is NULL for OptionsPattern. */
static PatKind pat_kind(const Expr* core, const Expr** ep, const Expr** blank0) {
    size_t ac = core->type == EXPR_FUNCTION ? core->data.function.arg_count : 0;
    *ep = NULL;
    if (pat_is(core, SYM_BlankSequence) || pat_is(core, SYM_BlankNullSequence)) {
        *ep = core;               /* read as Blank[h] by pat_elem_covers */
        *blank0 = core;
        return pat_is(core, SYM_BlankSequence) ? PK_SEQ1 : PK_SEQ0;
    }
    if ((pat_is(core, SYM_Repeated) || pat_is(core, SYM_RepeatedNull)) && ac == 1) {
        *ep = core->data.function.args[0];
        return pat_is(core, SYM_Repeated) ? PK_SEQ1 : PK_SEQ0;
    }
    if (pat_is(core, SYM_OptionsPattern)) return PK_SEQ0;
    if (pat_is(core, SYM_Optional) && (ac == 1 || ac == 2)) {
        *ep = core->data.function.args[0];
        return PK_OPT;
    }
    if (pat_is(core, SYM_Repeated) || pat_is(core, SYM_RepeatedNull) ||
        pat_is(core, SYM_PatternSequence) || pat_is(core, SYM_Longest) ||
        pat_is(core, SYM_Shortest) || pat_is(core, SYM_OrderlessPatternSequence))
        return PK_OTHER;
    return PK_ONE;
}

static bool pat_sub1(const Expr* g, const Expr* s, int depth);

/* For Plus[..], Times[..] or Power[..] whose arguments are all Optional but
 * one, that one argument; NULL otherwise. */
static const Expr* pat_defaults_to_one_arg(const Expr* p) {
    if (!(pat_is(p, SYM_Plus) || pat_is(p, SYM_Times) || pat_is(p, SYM_Power))) return NULL;
    const Expr* req = NULL;
    for (size_t i = 0; i < p->data.function.arg_count; i++) {
        const Expr* a = pat_unname(p->data.function.args[i]);
        if (pat_is(a, SYM_Optional)) continue;
        if (req) return NULL;
        req = p->data.function.args[i];
    }
    return req;
}

/* A pattern-free left-hand side (f[1], f[g[a]]): it matches one expression,
 * and the reference tries every such definition before any pattern. */
static bool pat_is_literal(const Expr* p) {
    if (p->type != EXPR_FUNCTION) return true;
    if (pat_is_construct(p)) return false;
    if (!pat_is_literal(p->data.function.head)) return false;
    for (size_t i = 0; i < p->data.function.arg_count; i++)
        if (!pat_is_literal(p->data.function.args[i])) return false;
    return true;
}

/* The head an atom has to the matcher: _Integer matches 3. */
static const char* pat_atom_head(const Expr* e) {
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: return SYM_Integer;
        case EXPR_REAL: return SYM_Real;
        case EXPR_STRING: return SYM_String;
        case EXPR_SYMBOL: return SYM_Symbol;
        default: return NULL;
    }
}

/* Is a typed-blank element pattern Blank[h] / __h / ___h (h possibly absent)
 * general enough for s (one element)? */
static bool pat_blank_covers(const Expr* blank, const Expr* s) {
    if (blank->data.function.arg_count == 0) return true;
    const Expr* h = blank->data.function.args[0];
    if (h->type != EXPR_SYMBOL) return false;
    if (pat_is(s, SYM_Blank) || pat_is(s, SYM_BlankSequence) || pat_is(s, SYM_BlankNullSequence))
        return s->data.function.arg_count == 1 && expr_eq(s->data.function.args[0], h);
    if (s->type == EXPR_FUNCTION) {
        if (pat_is_construct(s)) return false;
        return s->data.function.head->type == EXPR_SYMBOL &&
               s->data.function.head->data.symbol.name == h->data.symbol.name;
    }
    return pat_atom_head(s) == h->data.symbol.name;
}

/* Does element pattern ep (from a g-side sequence; NULL = OptionsPattern)
 * cover the s-side argument element s? */
static bool pat_elem_covers(const Expr* ep, const Expr* s_elem, int depth) {
    const Expr* s = pat_untested(s_elem);
    const Expr* sep = NULL; const Expr* sb = NULL;
    PatKind sk = pat_kind(s, &sep, &sb);
    if (sk == PK_OTHER) return false;
    if (!ep) {                                         /* OptionsPattern[] */
        if (sk == PK_SEQ0 && !sep && !sb) return true;  /* s: OptionsPattern too */
        const Expr* one = sk == PK_OPT ? pat_untested(sep) : s;
        return (sk == PK_ONE || sk == PK_OPT) &&
               (pat_is(one, SYM_Rule) || pat_is(one, SYM_RuleDelayed)) && !pat_is_construct(one);
    }
    bool ep_blank = pat_is(ep, SYM_Blank) || pat_is(ep, SYM_BlankSequence) ||
                    pat_is(ep, SYM_BlankNullSequence);
    if (sk == PK_ONE || sk == PK_OPT) {
        const Expr* one = sk == PK_OPT ? pat_untested(sep) : s;
        return ep_blank ? pat_blank_covers(ep, one) : pat_sub1(ep, one, depth + 1);
    }
    /* s is itself a sequence: each of its elements must be covered. */
    /* __h, ___h: covered by a blank sequence; x:p.. counts as the narrower
     * of x:_.. and x__, as in the reference. */
    if (sb) return (pat_is(ep, SYM_BlankSequence) || pat_is(ep, SYM_BlankNullSequence)) &&
                   pat_blank_covers(ep, sb);
    if (!sep) return ep_blank && ep->data.function.arg_count == 0;   /* OptionsPattern */
    return ep_blank ? pat_blank_covers(ep, pat_untested(sep)) : pat_sub1(ep, sep, depth + 1);
}

/* Do g's argument elements G[i..] cover s's S[j..] -- every argument sequence
 * S matches, G matches?  memo[i * (ns + 1) + j]: 0 unknown, 1 yes, 2 no. */
static bool pat_cover(const Expr* const* G, size_t ng, size_t i,
                      const Expr* const* S, size_t ns, size_t j,
                      unsigned char* memo, int depth) {
    if (i == ng) return j == ns;
    unsigned char* m = &memo[i * (ns + 1) + j];
    if (*m) return *m == 1;
    bool res = false;
    const Expr* g = pat_unname(G[i]);
    const Expr* gep = NULL; const Expr* gb = NULL;
    PatKind gk = pat_kind(g, &gep, &gb);
    const Expr* sj = j < ns ? pat_untested(S[j]) : NULL;
    const Expr* sep = NULL; const Expr* sb = NULL;
    PatKind sk = sj ? pat_kind(sj, &sep, &sb) : PK_OTHER;
    switch (gk) {
        case PK_ONE:
            res = sj && sk == PK_ONE && pat_sub1(G[i], S[j], depth + 1) &&
                  pat_cover(G, ng, i + 1, S, ns, j + 1, memo, depth);
            break;
        case PK_OPT:
            res = pat_cover(G, ng, i + 1, S, ns, j, memo, depth) ||
                  (sj && (sk == PK_ONE || sk == PK_OPT) &&
                   pat_sub1(gep, sk == PK_OPT ? sep : S[j], depth + 1) &&
                   pat_cover(G, ng, i + 1, S, ns, j + 1, memo, depth));
            break;
        case PK_SEQ0: case PK_SEQ1: {
            bool nonempty = false;
            for (size_t k = j; ; k++) {
                if (k > j) {
                    const Expr* t = pat_untested(S[k - 1]);
                    const Expr* tep = NULL; const Expr* tb = NULL;
                    PatKind tk = pat_kind(t, &tep, &tb);
                    if (!pat_elem_covers(gep, S[k - 1], depth)) break;
                    if (tk == PK_ONE || tk == PK_SEQ1) nonempty = true;
                }
                if ((gk == PK_SEQ0 || nonempty) &&
                    pat_cover(G, ng, i + 1, S, ns, k, memo, depth)) { res = true; break; }
                if (k == ns) break;
            }
            break;
        }
        default:
            break;
    }
    *m = res ? 1 : 2;
    return res;
}

/* Does the single pattern g match every expression the single pattern s does? */
static bool pat_sub1(const Expr* g, const Expr* s, int depth) {
    if (depth > 64) return false;
    if (expr_eq(g, s)) return true;
    g = pat_unname(g);
    s = pat_unname(s);
    if (expr_eq(g, s)) return true;
    if (pat_is(g, SYM_Condition)) return false;                 /* g narrowed */
    if (pat_is(g, SYM_PatternTest) && g->data.function.arg_count == 2) {
        while (pat_is(s, SYM_Condition) && s->data.function.arg_count == 2)
            s = pat_unname(s->data.function.args[0]);
        return pat_is(s, SYM_PatternTest) && s->data.function.arg_count == 2 &&
               expr_eq(g->data.function.args[1], s->data.function.args[1]) &&
               pat_sub1(g->data.function.args[0], s->data.function.args[0], depth + 1);
    }
    if ((pat_is(s, SYM_Condition) || pat_is(s, SYM_PatternTest)) && s->data.function.arg_count == 2)
        return pat_sub1(g, s->data.function.args[0], depth + 1);
    if (pat_is(s, SYM_Alternatives)) {
        if (s->data.function.arg_count == 0) return false;
        for (size_t i = 0; i < s->data.function.arg_count; i++)
            if (!pat_sub1(g, s->data.function.args[i], depth + 1)) return false;
        return true;
    }
    if (pat_is(g, SYM_Alternatives)) {
        for (size_t i = 0; i < g->data.function.arg_count; i++)
            if (pat_sub1(g->data.function.args[i], s, depth + 1)) return true;
        return false;
    }
    const Expr* ep = NULL; const Expr* b0 = NULL;
    if (pat_kind(s, &ep, &b0) != PK_ONE) return false;      /* s is not one element */
    /* A sum, product or power with default arguments (x_^n_., a_. x_,
     * x_ + y_.) also matches what its one required argument matches, by the
     * defaults: it and a lone blank cover each other. */
    {
        const Expr* req = pat_defaults_to_one_arg(g);
        if (req) return pat_sub1(req, s, depth + 1);
    }
    if (pat_is(g, SYM_Blank)) return pat_blank_covers(g, s);
    if (pat_is(g, SYM_Optional) && g->data.function.arg_count >= 1)
        return pat_sub1(g->data.function.args[0], s, depth + 1);
    if (g->type != EXPR_FUNCTION || pat_is_construct(g)) return false;
    if (s->type != EXPR_FUNCTION || pat_is_construct(s)) return false;
    if (!pat_sub1(g->data.function.head, s->data.function.head, depth + 1)) return false;
    size_t ng = g->data.function.arg_count, ns = s->data.function.arg_count;
    if (ng > 64 || ns > 64) return false;
    unsigned char memo_buf[256];
    size_t need = (ng + 1) * (ns + 1);
    unsigned char* memo = need <= sizeof memo_buf ? memo_buf : calloc(need, 1);
    if (!memo) return false;
    if (memo == memo_buf) memset(memo_buf, 0, need);
    bool r = pat_cover((const Expr* const*)g->data.function.args, ng, 0,
                       (const Expr* const*)s->data.function.args, ns, 0, memo, depth);
    if (memo != memo_buf) free(memo);
    return r;
}

/* Collect pattern names; true if one occurs twice (a nonlinear pattern such as
 * f[x_, x_] matches less than its shape suggests). */
static bool pat_repeats_name(const Expr* p, const char** seen, size_t* n, size_t cap) {
    if (p->type != EXPR_FUNCTION) return false;
    if (pat_is(p, SYM_Pattern) && p->data.function.arg_count == 2 &&
        p->data.function.args[0]->type == EXPR_SYMBOL) {
        const char* nm = p->data.function.args[0]->data.symbol.name;
        for (size_t i = 0; i < *n; i++) if (seen[i] == nm) return true;
        if (*n == cap) return true;                     /* too many to track */
        seen[(*n)++] = nm;
    }
    if (pat_repeats_name(p->data.function.head, seen, n, cap)) return true;
    for (size_t i = 0; i < p->data.function.arg_count; i++)
        if (pat_repeats_name(p->data.function.args[i], seen, n, cap)) return true;
    return false;
}

static bool pattern_is_nonlinear(const Expr* p) {
    const char* seen[64]; size_t n = 0;
    return pat_repeats_name(p, seen, &n, 64);
}

/* Does pattern g match every expression pattern s matches?  g_nonlinear is
 * pattern_is_nonlinear(g), which the caller has cached. */
static bool pattern_subsumes(const Expr* g, bool g_nonlinear, const Expr* s) {
    if (expr_eq(g, s)) return true;
    if (g_nonlinear) return false;
    return pat_sub1(g, s, 0);
}

/* ---- Literal-rule index -------------------------------------------------
 *
 * Pattern-free definitions lead the DownValue list (see add_rule), and a
 * pattern-free rule matches exactly the expression equal to its left-hand side,
 * so they can be found by hash.  A memoised f[n_] := f[n] = ... over 10^4 values
 * otherwise scanned (and matched against) every stored value per call and per
 * new value: quadratic.  The index lists the rules of the literal prefix,
 * keyed by expr_hash of the pattern; `tail` is the prefix's last rule.  It is
 * trusted only while `head` is still the list head (Block and the optimisers
 * swap whole lists from outside), and every removal in this file drops it.
 * `multi`: two prefix rules share a left-hand side (conditional definitions
 * kept side by side) -- lookups then fall back to the scan. */
typedef struct LitIndex {
    const Rule* head;
    Rule* tail;
    Rule** slot;
    size_t n, cap;
    bool multi;
} LitIndex;

static void litidx_drop(SymbolDef* d) {
    if (!d || !d->lit_index) return;
    free(d->lit_index->slot);
    free(d->lit_index);
    d->lit_index = NULL;
}

static Rule* litidx_find(const LitIndex* ix, const Expr* e) {
    if (!ix->cap) return NULL;
    size_t mask = ix->cap - 1;
    for (size_t i = (size_t)expr_hash(e) & mask; ix->slot[i]; i = (i + 1) & mask)
        if (expr_eq(ix->slot[i]->pattern, e)) return ix->slot[i];
    return NULL;
}

static bool litidx_put(LitIndex* ix, Rule* r) {
    if (2 * (ix->n + 1) > ix->cap) {
        size_t nc = ix->cap ? 2 * ix->cap : 16;
        Rule** ns = calloc(nc, sizeof(Rule*));
        if (!ns) return false;
        for (size_t i = 0; i < ix->cap; i++) {
            Rule* q = ix->slot[i];
            if (!q) continue;
            size_t j = (size_t)expr_hash(q->pattern) & (nc - 1);
            while (ns[j]) j = (j + 1) & (nc - 1);
            ns[j] = q;
        }
        free(ix->slot);
        ix->slot = ns;
        ix->cap = nc;
    }
    if (litidx_find(ix, r->pattern)) ix->multi = true;
    size_t mask = ix->cap - 1;
    size_t i = (size_t)expr_hash(r->pattern) & mask;
    while (ix->slot[i]) i = (i + 1) & mask;
    ix->slot[i] = r;
    ix->n++;
    return true;
}

/* The valid index of d's DownValues, built when missing or stale; NULL on
 * allocation failure (callers then scan). */
static LitIndex* litidx_get(SymbolDef* d) {
    if (d->lit_index && d->lit_index->head == d->down_values) return d->lit_index;
    litidx_drop(d);
    LitIndex* ix = calloc(1, sizeof(LitIndex));
    if (!ix) return NULL;
    for (Rule* r = d->down_values; r && r->literal; r = r->next) {
        if (!litidx_put(ix, r)) { free(ix->slot); free(ix); return NULL; }
        ix->tail = r;
    }
    ix->head = d->down_values;
    d->lit_index = ix;
    return ix;
}

/* `owner`: the symbol whose DownValues `list` is, so its literal index is
 * used and kept; NULL for every other rule list.  `binding`: the symbol whose
 * OwnValues `list` is, NULL for every other rule list. */
static void add_rule_in(Rule** list, Expr* pattern, Expr* replacement, bool delayed,
                        SymbolDef* owner, SymbolDef* binding) {
    /* Any change to an OwnValue/DownValue list could change the meaning
     * of any expression that touches this symbol, so bump the global
     * eval clock. This invalidates every cached evaluation in one shot
     * (conservative variant of EVAL_IMPROVEMENTS_PLAN §3.3).  An OwnValue
     * binding bumps the clock too, but invalidates the INERT fixed points
     * only if one of them relied on this symbol being unbound (eval.c). */
    if (binding) eval_clock_bump_binding(binding);
    else         eval_clock_bump();

    /* §3.4 canonicalize the LHS before deriving any dispatch keys or
     * doing duplicate-rule detection. Held LHSs would otherwise sit in
     * parser-shape while every runtime input arrives in evaluated-shape;
     * see the long comment above pattern_canonicalize for the full
     * rationale. We retain ownership of `canon`. */
    Expr* canon = pattern_canonicalize(expr_copy(pattern));

    int32_t arity = pattern_dispatch_arity(canon);
    const char* first_head = (arity == -1) ? NULL : pattern_first_arg_head(canon);
    int32_t spec = pattern_specificity(canon);
    uint64_t ahash = 0;      /* computed on the first candidate, if any */

    /* A rule whose LHS is the same "up to renaming of patterns" (and the
     * same RHS condition, which canonicalisation folds into the pattern)
     * replaces the old one, pattern AND right-hand side -- the RHS refers to
     * the new names: f[x_] := x^2; f[y_] := y^3 leaves one rule, y^3. The
     * exact match is the cheap first test; the alpha-normalised comparison
     * (two deep copies) only runs on an alpha-hash hit, which every
     * alpha-equivalent pattern shares, so defining n distinct rules stays an
     * integer compare per stored rule. Both hashes are computed lazily, only
     * once a stored rule passes the cheap key test: an OwnValue update (the
     * same symbol, caught by expr_eq) never hashes at all. */
    Expr* key = NULL;
    uint8_t scoped = rhs_has_scoped_condition(replacement) ? 1 : 0;
    bool literal = pat_is_literal(canon);
    LitIndex* ix = owner ? litidx_get(owner) : NULL;
    if (ix && ix->multi) { litidx_drop(owner); ix = NULL; }
    Rule* scan_from = *list;
    if (ix) {
        /* A pattern-free left-hand side can only equal a rule of the literal
         * prefix, and a pattern one only a rule after it. */
        if (literal) {
            Rule* dup = litidx_find(ix, canon);
            scan_from = NULL;
            if (dup) {
                if ((dup->scoped_cond || scoped) && !expr_eq(dup->replacement, replacement)) {
                    litidx_drop(owner);       /* both kept: let the scan decide */
                    ix = NULL;
                    scan_from = *list;
                } else {
                    expr_free(dup->replacement);
                    dup->replacement = expr_copy(replacement);
                    dup->has_return = mentions_return(replacement);
                    dup->delayed = delayed;
                    expr_free(canon);
                    return;
                }
            }
        } else {
            scan_from = ix->tail ? ix->tail->next : *list;
        }
    }
    for (Rule* curr = scan_from; curr; curr = curr->next) {
        /* A conditional rule (scoped guard) is not replaced by a same-LHS
         * definition with a different RHS: both are kept, in order. */
        if ((curr->scoped_cond || scoped) && !expr_eq(curr->replacement, replacement))
            continue;
        if (expr_eq(curr->pattern, canon)) {
            // Replace replacement in place; sort key is unchanged.
            expr_free(curr->replacement);
            curr->replacement = expr_copy(replacement);
            curr->has_return = mentions_return(replacement);
            curr->delayed = delayed;
            expr_free(canon);
            if (key) expr_free(key);
            return;
        }
        if (curr->specificity != spec || curr->dispatch_arity != arity ||
            curr->first_arg_head_canon != first_head)
            continue;
        if (!ahash) ahash = pattern_alpha_hash(canon);
        if (!curr->alpha_hash) curr->alpha_hash = pattern_alpha_hash(curr->pattern);
        if (curr->alpha_hash != ahash)
            continue;
        if (!key) key = pattern_alpha_normalize(canon);
        Expr* curr_key = pattern_alpha_normalize(curr->pattern);
        bool same = expr_eq(curr_key, key);
        expr_free(curr_key);
        if (same) {
            expr_free(curr->pattern);
            curr->pattern = canon;          /* same alpha_hash: keep it */
            expr_free(curr->replacement);
            curr->replacement = expr_copy(replacement);
            curr->has_return = mentions_return(replacement);
            curr->delayed = delayed;
            expr_free(key);
            return;
        }
    }
    if (key) expr_free(key);

    Rule* new_rule = malloc(sizeof(Rule));
    new_rule->pattern = canon;
    new_rule->replacement = expr_copy(replacement);
    new_rule->dispatch_arity = arity;
    new_rule->first_arg_head_canon = first_head;
    new_rule->specificity = spec;
    new_rule->alpha_hash = ahash;
    new_rule->delayed = delayed;
    new_rule->scoped_cond = scoped;
    new_rule->has_return = mentions_return(replacement);
    new_rule->literal = literal ? 1 : 0;
    new_rule->nonlinear = pattern_is_nonlinear(canon) ? 1 : 0;
    new_rule->next = NULL;

    /* Insert before the first stored rule that is strictly more general (it
     * matches everything the new pattern does, and more): f[2] goes ahead of
     * f[x_], and so does f[y_] /; y > 0 (its guard narrows it).  Patterns
     * neither of which covers the other keep definition order, as in the
     * reference; see pattern_subsumes.  Pattern-free definitions come first
     * of all, in the order they were made. */
    Rule* prev = ix ? ix->tail : NULL;
    Rule* curr = prev ? prev->next : *list;
    while (curr && curr->literal) {                    /* literal rules lead */
        prev = curr;
        curr = curr->next;
    }
    Expr* new_key = NULL;
    while (!literal && curr) {
        /* Different fixed arities, or different heads of a fixed first
         * argument, match disjoint sets: no need to compare further. */
        bool disjoint =
            (arity >= 0 && curr->dispatch_arity >= 0 && arity != curr->dispatch_arity) ||
            (first_head && curr->first_arg_head_canon && first_head != curr->first_arg_head_canon);
        if (!disjoint && pattern_subsumes(curr->pattern, curr->nonlinear, canon)) {
            /* Strictly more general: the new rule goes first.  Of two that
             * cover each other (x_ and x_^n_., h[x___, y___] and h[z___]) the
             * newer goes first too, unless they are the same pattern up to
             * renaming (conditional definitions f[x_] := a /; c kept side by
             * side), which keep definition order. */
            if (!pattern_subsumes(canon, new_rule->nonlinear, curr->pattern)) break;
            if (!new_key) new_key = pattern_alpha_normalize(canon);
            Expr* ck = pattern_alpha_normalize(curr->pattern);
            bool same = expr_eq(ck, new_key);
            expr_free(ck);
            if (!same) break;
        }
        prev = curr;
        curr = curr->next;
    }
    if (new_key) expr_free(new_key);
    new_rule->next = curr;
    if (prev) {
        prev->next = new_rule;
    } else {
        *list = new_rule;
    }
    if (ix) {
        if (literal) {
            if (litidx_put(ix, new_rule)) ix->tail = new_rule;
            else { litidx_drop(owner); return; }
        }
        ix->head = *list;
    }
}

static void add_rule(Rule** list, Expr* pattern, Expr* replacement, bool delayed) {
    add_rule_in(list, pattern, replacement, delayed, NULL, NULL);
}
void symtab_add_own_value(const char* symbol_name, Expr* pattern, Expr* replacement) {
    symtab_add_own_value_as(symbol_name, pattern, replacement, false);
}

void symtab_add_down_value(const char* symbol_name, Expr* pattern, Expr* replacement) {
    symtab_add_down_value_as(symbol_name, pattern, replacement, true);
}

void symtab_add_own_value_as(const char* symbol_name, Expr* pattern, Expr* replacement,
                             bool delayed) {
    SymbolDef* def = symtab_get_def(symbol_name);
    add_rule_in(&def->own_values, pattern, replacement, delayed, NULL, def);
    /* An OwnValue binding is "soft": add_rule bumps only the eval clock, leaving
     * the rule epoch alone so a loop-invariant GROUND value survives an iterator
     * rebinding (the whole point of this path staying cheap). The one exception
     * is a Protected target: one of the pure constructor heads a GROUND node is
     * built from could in principle be bound by internal code, which would change
     * those nodes' meaning -- advance the rule epoch defensively so they
     * re-evaluate. (User Set on a Protected symbol never reaches here; it errors
     * first.) */
    if (get_attributes_def(def) & ATTR_PROTECTED) eval_rule_epoch_mark();
}

/* Unlink and return the DownValue of symbol_name whose left-hand side is
 * exactly lhs (after the canonicalisation every stored LHS gets), or NULL.
 * The caller owns the returned rule (free it with symtab_rule_free, or
 * re-install its pattern and replacement). Backs the iterators' indexed
 * variables (Table[..., {v[1], 3}]), which bind v[1] as a DownValue and put
 * back any v[1] rule that was there before. */
Rule* symtab_detach_down_value_exact(const char* symbol_name, Expr* lhs) {
    SymbolDef* def = symtab_get_def(symbol_name);
    if (!def->down_values) return NULL;
    Expr* canon = pattern_canonicalize(expr_copy(lhs));
    litidx_drop(def);
    Rule** pp = &def->down_values;
    Rule* found = NULL;
    for (; *pp; pp = &(*pp)->next) {
        if (expr_eq((*pp)->pattern, canon)) {
            found = *pp;
            *pp = found->next;
            found->next = NULL;
            eval_clock_bump();
            break;
        }
    }
    expr_free(canon);
    return found;
}

void symtab_rule_free(Rule* r) {
    if (!r) return;
    r->next = NULL;
    free_rule_list(r);
}

void symtab_add_down_value_as(const char* symbol_name, Expr* pattern, Expr* replacement,
                              bool delayed) {
    SymbolDef* def = symtab_get_def(symbol_name);
    /* `binding` = def: an INERT fixed point is invalidated only if it relied on
     * this symbol (eval_clock_bump_binding); the rule epoch below still
     * invalidates every GROUND node and every exact-clock stamp. */
    add_rule_in(&def->down_values, pattern, replacement, delayed, def, def);
    /* A DownValue changes how this head evaluates -- advance the rule epoch so
     * any GROUND node that uses it as a head is re-checked (add_rule already
     * bumped the eval clock). */
    eval_rule_epoch_mark();
}

void symtab_clear_symbol(const char* symbol_name) {
    /* Removing rules invalidates cached evaluations that might have relied on
     * them; bump the eval clock unconditionally. This also removes DownValues,
     * which changes how the head evaluates, so advance the rule epoch too. */
    eval_rule_epoch_bump();

    SymbolDef* def = symtab_get_def(symbol_name);
    Rule* curr = def->own_values;
    while (curr) {
        Rule* next = curr->next;
        expr_free(curr->pattern);
        expr_free(curr->replacement);
        free(curr);
        curr = next;
    }
    def->own_values = NULL;

    litidx_drop(def);
    curr = def->down_values;
    while (curr) {
        Rule* next = curr->next;
        expr_free(curr->pattern);
        expr_free(curr->replacement);
        free(curr);
        curr = next;
    }
    def->down_values = NULL;
    free_rule_list(def->sub_values);
    def->sub_values = NULL;

#if UP_VALUES
    free_up_value_list(def->up_values);
    def->up_values = NULL;
#endif
    free_format_value_list(def->format_values);
    def->format_values = NULL;
    free_n_value_list(def->n_values);
    def->n_values = NULL;
}

void symtab_remove_symbol(const char* symbol_name) {
    SymbolDef* d = node_find(symbol_name);
    if (!d || !d->has_def) return;   /* nothing materialized to remove */
    /* Removing the definition invalidates any cached evaluation that referenced
     * this symbol's rules, including its DownValues/attributes -- advance the
     * rule epoch as well. */
    eval_rule_epoch_bump();
    /* Reset to interned-only. The node and its name string remain in the table
     * so a later reference recreates a fresh, empty definition and, crucially,
     * any live Expr (and SYM_*) that already hold the interned name pointer stay
     * valid. is_system is preserved. */
    reset_node_payload(d);
}

/* ============================================================
 * Unset: remove a single rule by left-hand side
 *
 * Mathematica's `lhs =.` (Unset[lhs]) drops the one rule whose LHS is
 * identical to `lhs` "up to renaming of patterns" -- i.e. up to
 * alpha-equivalence of the bound pattern variables. `f[x_] =.` therefore
 * removes a rule defined as `f[y_] := ...`, but `fact[1] =.` leaves the
 * separate `fact[n_]` rule untouched.
 *
 * We compare by alpha-normalizing both the query LHS and each stored rule
 * LHS: every Pattern variable name is rewritten to a positional placeholder
 * ($Pat$0, $Pat$1, ...) in first-occurrence order, after which two
 * alpha-equivalent patterns become structurally equal under expr_eq.
 * ============================================================ */

#define ALPHA_MAX_VARS 64

/* Collect the distinct Pattern[name, ...] variable names of `p` into
 * `names`, in first-occurrence (pre-order) order. Names are interned, so
 * identity comparison suffices for the dedup check. */
static void alpha_collect_names(const Expr* p, const char** names, int* n) {
    if (!p || p->type != EXPR_FUNCTION) return;
    Expr* head = p->data.function.head;
    if (head && head->type == EXPR_SYMBOL && head->data.symbol.name == SYM_Pattern
        && p->data.function.arg_count >= 1
        && p->data.function.args[0]
        && p->data.function.args[0]->type == EXPR_SYMBOL) {
        const char* nm = p->data.function.args[0]->data.symbol.name;
        int seen = 0;
        for (int i = 0; i < *n; i++) if (names[i] == nm) { seen = 1; break; }
        if (!seen && *n < ALPHA_MAX_VARS) names[(*n)++] = nm;
    }
    alpha_collect_names(head, names, n);
    for (size_t i = 0; i < p->data.function.arg_count; i++)
        alpha_collect_names(p->data.function.args[i], names, n);
}

/* In-place rewrite of every symbol of `p` that names a bound pattern
 * variable (found in `names`) to its positional placeholder. Symbol names
 * are interned, so we only ever swap the interned pointer -- no allocation
 * is freed and the placeholder is itself interned for expr_eq identity. */
static void alpha_rename(Expr* p, const char** names, int n) {
    if (!p) return;
    /* This rewrites symbol names in place; every node on the path down to a
     * renamed leaf has its structural hash changed, so drop each visited node's
     * memoized hash (cold path: pattern setup, once per rule). */
    expr_invalidate_hash(p);
    if (p->type == EXPR_SYMBOL) {
        for (int i = 0; i < n; i++) {
            if (p->data.symbol.name == names[i]) {
                char buf[24];
                snprintf(buf, sizeof(buf), "$Pat$%d", i);
                if (p->data.symbol.def && p->data.symbol.def->module_temp &&
                    p->data.symbol.def->temp_nodes > 0)
                    p->data.symbol.def->temp_nodes--;   /* no longer names the local */
                p->data.symbol.name = (char*)intern_symbol(buf);
                p->data.symbol.def = NULL;   /* name changed: invalidate cached def */
                break;
            }
        }
        return;
    }
    if (p->type == EXPR_FUNCTION) {
        alpha_rename(p->data.function.head, names, n);
        for (size_t i = 0; i < p->data.function.arg_count; i++)
            alpha_rename(p->data.function.args[i], names, n);
    }
}

/* A fully private (deep) copy of `e`: the result shares no mutable node with `e`,
 * so it is safe to rewrite in place. Needed because expr_copy() is only a
 * refcount bump (COW) and alpha_rename() mutates leaves at every depth -- a bare
 * expr_copy would let alpha_rename corrupt the shared original (it did: a
 * non-matching `f[x_] =.` on `f[x_] := rhs /; c` renamed the STORED pattern's x
 * to $Pat$0 while the rhs still said x). expr_unshare privatises one level and
 * drops our extra ref, so recursing privatises the whole tree. Caller owns it. */
static Expr* deep_copy_private(Expr* e) {
    if (!e) return NULL;
    Expr* top = expr_unshare(expr_copy(e));   /* refcount-1 top; children shared */
    if (top->type == EXPR_FUNCTION) {
        Expr* oh = top->data.function.head;
        top->data.function.head = deep_copy_private(oh);
        expr_free(oh);
        for (size_t i = 0; i < top->data.function.arg_count; i++) {
            Expr* oa = top->data.function.args[i];
            top->data.function.args[i] = deep_copy_private(oa);
            expr_free(oa);
        }
    }
    return top;
}

/* Return a freshly allocated, alpha-normalized copy of `p`. Caller owns it, and
 * `p` is left untouched (deep copy -- see deep_copy_private). */
static Expr* pattern_alpha_normalize(const Expr* p) {
    Expr* copy = deep_copy_private((Expr*)p);
    const char* names[ALPHA_MAX_VARS];
    int n = 0;
    alpha_collect_names(copy, names, &n);
    alpha_rename(copy, names, n);
    return copy;
}

static uint64_t alpha_mix(uint64_t h, uint64_t v) {
    h ^= v + UINT64_C(0x9e3779b97f4a7c15) + (h << 6) + (h >> 2);
    return h;
}

/* Structural hash of `p` in which a bound pattern variable contributes its
 * first-occurrence index rather than its name -- equal for alpha-equivalent
 * patterns, like expr_eq after pattern_alpha_normalize, but with no copy. */
static uint64_t alpha_hash_rec(const Expr* p, const char** names, int n) {
    if (p->type == EXPR_SYMBOL) {
        for (int i = 0; i < n; i++)
            if (p->data.symbol.name == names[i]) return alpha_mix(UINT64_C(0xa1fa), (uint64_t)i);
        return expr_hash(p);
    }
    if (p->type != EXPR_FUNCTION) return expr_hash(p);
    uint64_t h = alpha_mix(UINT64_C(0xf00d), p->data.function.arg_count);
    h = alpha_mix(h, alpha_hash_rec(p->data.function.head, names, n));
    for (size_t i = 0; i < p->data.function.arg_count; i++)
        h = alpha_mix(h, alpha_hash_rec(p->data.function.args[i], names, n));
    return h;
}

static uint64_t pattern_alpha_hash(const Expr* p) {
    const char* names[ALPHA_MAX_VARS];
    int n = 0;
    alpha_collect_names(p, names, &n);
    return alpha_hash_rec(p, names, n) | 1;   /* never 0: 0 means "not computed" */
}

bool symtab_remove_matching_rule(const char* symbol_name, const Expr* lhs,
                                 bool own_value) {
    if (!symbol_name || !lhs) return false;

    /* Canonicalize the query LHS the same way add_rule() canonicalized the
     * stored patterns, then alpha-normalize so renamed pattern variables
     * still compare equal. */
    Expr* canon = pattern_canonicalize(expr_copy((Expr*)lhs));
    Expr* key = pattern_alpha_normalize(canon);
    expr_free(canon);

    SymbolDef* def = symtab_get_def(symbol_name);
    Rule** list = own_value ? &def->own_values : &def->down_values;
    if (!own_value) litidx_drop(def);

    Rule* prev = NULL;
    bool removed = false;
    for (Rule* curr = *list; curr; prev = curr, curr = curr->next) {
        Expr* curr_key = pattern_alpha_normalize(curr->pattern);
        int match = expr_eq(curr_key, key);
        expr_free(curr_key);
        if (!match) continue;

        if (prev) prev->next = curr->next;
        else *list = curr->next;
        expr_free(curr->pattern);
        expr_free(curr->replacement);
        free(curr);
        /* Removing a rule invalidates cached evaluations that relied on it. */
        eval_clock_bump();
        removed = true;
        break;
    }

    expr_free(key);
    return removed;
}

/* ---- SubValues ---------------------------------------------------------- */

const char* symtab_innermost_head(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    const Expr* h = e->data.function.head;
    if (!h || h->type != EXPR_FUNCTION) return NULL;      /* not compound */
    while (h->type == EXPR_FUNCTION && h->data.function.head) h = h->data.function.head;
    return h->type == EXPR_SYMBOL ? h->data.symbol.name : NULL;
}

void symtab_add_sub_value(const char* symbol_name, Expr* pattern, Expr* replacement,
                          bool delayed) {
    SymbolDef* def = symtab_get_def(symbol_name);
    add_rule(&def->sub_values, pattern, replacement, delayed);
}

Rule* symtab_get_sub_values(const char* symbol_name) {
    SymbolDef* d = symtab_lookup(symbol_name);
    return d ? d->sub_values : NULL;
}

bool symtab_remove_matching_sub_value(const char* symbol_name, const Expr* lhs) {
    if (!symbol_name || !lhs) return false;
    Expr* canon = pattern_canonicalize(expr_copy((Expr*)lhs));
    Expr* key = pattern_alpha_normalize(canon);
    expr_free(canon);
    SymbolDef* def = symtab_get_def(symbol_name);
    Rule* prev = NULL;
    bool removed = false;
    for (Rule* curr = def->sub_values; curr; prev = curr, curr = curr->next) {
        Expr* curr_key = pattern_alpha_normalize(curr->pattern);
        int same = expr_eq(curr_key, key);
        expr_free(curr_key);
        if (!same) continue;
        if (prev) prev->next = curr->next;
        else def->sub_values = curr->next;
        expr_free(curr->pattern);
        expr_free(curr->replacement);
        free(curr);
        eval_clock_bump();
        removed = true;
        break;
    }
    expr_free(key);
    return removed;
}

/* Try `def`'s SubValues on `expr` (a call with a compound head whose innermost
 * symbol is def's): the first matching rule's bound right-hand side, with
 * OptionValue resolved and a scoped guard tested, or NULL. */
Expr* apply_sub_values(SymbolDef* def, Expr* expr) {
    if (!def || !def->sub_values || !expr || expr->type != EXPR_FUNCTION) return NULL;
    for (Rule* rule = def->sub_values; rule; rule = rule->next) {
        MatchEnv* env = env_new();
        const char* prev_head = match_set_rule_head(def->symbol_name);
        Expr* result = match(expr, rule->pattern, env) ? replace_bindings(rule->replacement, env) : NULL;
        match_set_rule_head(prev_head);
        if (!result) { env_free(env); continue; }
        Expr* opts = env_get(env, "$OptionsPattern$");
        if (opts) {
            const char* hs = def->symbol_name;
            Expr* oh = env_get(env, "$OptionsPatternHead$");
            if (oh && oh->type == EXPR_SYMBOL) hs = oh->data.symbol.name;
            Expr* resolved = optionvalue_inject_context(result, hs, opts);
            expr_free(result);
            result = resolved;
        }
        env_free(env);
        if (rule->scoped_cond) {
            Expr* v = try_scoped_condition(result);
            expr_free(result);
            if (v) return v;
            continue;
        }
        return rule_result(rule, result, def->symbol_name);
    }
    return NULL;
}

void symtab_for_each(void (*visit)(const char* name, SymbolDef* def,
                                   void* user), void* user) {
    if (!visit) return;
    for (size_t i = 0; i < symtab_size; i++)
        for (SymbolDef* d = symtab[i]; d; d = d->next)
            if (d->has_def)   /* skip interned-only nodes: they are not "symbols" */
                visit(d->symbol_name, d, user);
}

void symtab_clear(void) {
    /* Reset every definition to interned-only. Nodes and their name strings are
     * retained so SYM_* and any Expr created before the clear keep valid
     * pointers -- callers routinely clear between subtests and keep evaluating.
     * Use intern_clear() for a full teardown that also frees the name strings. */
    for (size_t i = 0; i < symtab_size; i++)
        for (SymbolDef* d = symtab[i]; d; d = d->next)
            reset_node_payload(d);
}


Rule* symtab_get_own_values(const char* symbol_name) {
    SymbolDef* d = node_find(symbol_name);
    return d ? d->own_values : NULL;
}

Rule* symtab_get_down_values(const char* symbol_name) {
    SymbolDef* d = node_find(symbol_name);
    return d ? d->down_values : NULL;
}

/* §3.5 input-side dispatch key. Returns the canonical head of an
 * input-call's first argument, or NULL when we can't safely
 * specialize on it. The classification mirrors pattern_arg_head_canon
 * so a head set on a rule's pattern will line up with the head we
 * compute here for an input expression. */
static const char* input_arg_head_canon(const Expr* arg) {
    if (!arg) return NULL;
    if (arg->type == EXPR_INTEGER || arg->type == EXPR_BIGINT) {
        return intern_symbol("Integer");
    }
    if (arg->type == EXPR_REAL) return intern_symbol("Real");
    if (arg->type == EXPR_STRING) return intern_symbol("String");
    if (arg->type == EXPR_SYMBOL) return intern_symbol("Symbol");
    if (arg->type == EXPR_FUNCTION && arg->data.function.head &&
        arg->data.function.head->type == EXPR_SYMBOL) {
        return intern_symbol(arg->data.function.head->data.symbol.name);
    }
    return NULL;
}

/* Phase 3a: def-threaded entry point. The evaluator resolves the head's def
 * once per evaluate_step and passes it here, so DownValue dispatch does not
 * re-hash the head. `def` must be the definition of `expr`'s (symbol) head. */
Expr* apply_down_values_def(SymbolDef* def, Expr* expr) {
    if (!def || !expr || expr->type != EXPR_FUNCTION) return NULL;
    if (!def->down_values) return NULL;

    /* §3.5 dispatch filter: pre-compute the input's shape once, then
     * compare against each rule's cached key before invoking the
     * (much more expensive) matcher. A NULL key on either side means
     * "wildcard" -- we always scan in that case. */
    const int32_t input_arity = (int32_t)expr->data.function.arg_count;
    const char* input_first_head =
        (input_arity > 0)
            ? input_arg_head_canon(expr->data.function.args[0])
            : NULL;
    /* A Flat head rewrites a SUBSET of a longer operand sequence (Mathematica
     * semantics; see match_flat_leftover). Its dispatch key is then a LOWER
     * bound, not an exact arity, and the matched run need not start at arg 0,
     * so both pre-filters relax for a Flat head. */
    const bool down_flat = (def->attributes & ATTR_FLAT) != 0;
    const bool down_relaxed = (def->attributes & (ATTR_FLAT | ATTR_ORDERLESS)) != 0;

    /* Pattern-free rules lead the list and each matches only the expression
     * equal to it: try the one the index names, then go on after the prefix.
     * (Flat and Orderless heads match a literal f[a, b] against more than its
     * equal, so they scan.) */
    Rule* hit = NULL;
    Rule* start = def->down_values;
    if (start->literal && !down_relaxed) {
        LitIndex* ix = litidx_get(def);
        if (ix && !ix->multi) {
            hit = litidx_find(ix, expr);
            start = ix->tail ? ix->tail->next : start;
        }
    }

    for (Rule* rule = hit ? hit : start; rule;
         rule = (rule == hit) ? start : rule->next) {
        /* Arity filter: only skip when both sides committed to a
         * specific arity that disagrees. A rule with -1 arity (variable)
         * always passes; an input never has variable arity. For a Flat head
         * any input_arity >= the rule's arity can match (surplus -> leftover). */
        if (rule->dispatch_arity != -1 && rule->dispatch_arity != input_arity &&
            !(down_flat && input_arity >= rule->dispatch_arity)) {
            continue;
        }
        /* First-arg head filter: only skip when the rule expects a
         * specific head and we determined the input's first-arg head and
         * the two disagree. Pointer compare is sound because both are
         * interned. Skipped for a Flat head (the matched run can sit
         * anywhere) and for an Orderless one, whose matcher may bind the
         * pattern's first slot to any argument: Pair[n_ x_Momentum, y_] fires
         * on Pair[Momentum[p], lm Momentum[p]] through the second, as the
         * UpValue filter below already allowed. */
        if (!down_relaxed && rule->first_arg_head_canon && input_first_head &&
            rule->first_arg_head_canon != input_first_head) {
            continue;
        }

        MatchEnv* env = env_new();
        /* Phase 1 stays inline (identical to the pre-feature fast path); only a
         * failed match falls through to Flat-head leftover matching (Phase 2). */
        const char* prev_head = match_set_rule_head(def->symbol_name);
        Expr* result = match(expr, rule->pattern, env)
            ? replace_bindings(rule->replacement, env)
            : match_flat_leftover(expr, rule->pattern, rule->replacement, env);
        match_set_rule_head(prev_head);
        if (result && rule->scoped_cond) {
            /* Resolve OptionValue first (below), then test the scoped guard. */
            Expr* opts0 = env_get(env, "$OptionsPattern$");
            if (opts0) {
                const char* hs = def->symbol_name;
                Expr* oh = env_get(env, "$OptionsPatternHead$");
                if (oh && oh->type == EXPR_SYMBOL) hs = oh->data.symbol.name;
                Expr* resolved = optionvalue_inject_context(result, hs, opts0);
                expr_free(result);
                result = resolved;
            }
            Expr* v = try_scoped_condition(result);
            expr_free(result);
            env_free(env);
            if (v) return v;
            continue;                       /* guard failed: next rule */
        }
        if (result) {
            /* If the matched pattern carried an OptionsPattern, rewrite the
             * context-dependent OptionValue[...] nodes in the RHS into explicit
             * OptionValue[head, opts, name] form so they resolve against the
             * supplied options and the head's defaults. */
            Expr* opts = env_get(env, "$OptionsPattern$");
            if (opts) {
                /* OptionsPattern[s] resolves an unpassed option's default
                 * against s's Options (A9); a bare OptionsPattern[] falls back
                 * to the enclosing definition's own head. */
                const char* head_sym = def->symbol_name;
                Expr* optpat_head = env_get(env, "$OptionsPatternHead$");
                if (optpat_head && optpat_head->type == EXPR_SYMBOL)
                    head_sym = optpat_head->data.symbol.name;
                Expr* resolved =
                    optionvalue_inject_context(result, head_sym, opts);
                expr_free(result);
                result = resolved;
            }
            env_free(env);
            return rule_result(rule, result, def->symbol_name);
        }
        env_free(env);
    }
    return NULL;
}

/* Compatibility wrapper: resolve the head's def, then delegate. Used by callers
 * outside the evaluator hot path that only have the expression in hand. */
Expr* apply_down_values(Expr* expr) {
    if (!expr || expr->type != EXPR_FUNCTION) return NULL;
    Expr* head = expr->data.function.head;
    if (!head || head->type != EXPR_SYMBOL) return NULL;
    return apply_down_values_def(symtab_get_def(head->data.symbol.name), expr);
}

Expr* apply_own_values(Expr* expr) {
    if (!expr || expr->type != EXPR_SYMBOL) return NULL;

    /* Reuse the node's cached def pointer when present (same fast path the
     * evaluator uses for function heads), resolving and caching on first
     * touch. The def node is stable (never freed/reallocated), so the cache
     * stays valid even across redefinition; it is benign metadata and safe to
     * write on a shared node. Saves a symbol-table hash on every symbol read. */
    SymbolDef* def = expr->data.symbol.def;
    if (!def) {
        def = symtab_get_def(expr->data.symbol.name);
        expr->data.symbol.def = def;
    }
    Rule* rule = def->own_values;
    while (rule) {
        /* Fast path for the ubiquitous `x = value` OwnValue: the pattern is a
         * bare symbol keyed to this very symbol (own_values live on the
         * symbol's own def), so it matches unconditionally with no bindings.
         * Skip the matcher and its three-malloc MatchEnv entirely and clone
         * the replacement directly — this is the hottest read in tight loops
         * like Do[x = f[x], {n}]. The name compare is a cheap interned-pointer
         * check that also guards against any exotic caller. */
        if (rule->pattern->type == EXPR_SYMBOL &&
            rule->pattern->data.symbol.name == expr->data.symbol.name) {
            return rule_result(rule, expr_copy(rule->replacement), def->symbol_name);
        }
        MatchEnv* env = env_new();
        if (match(expr, rule->pattern, env)) {
            Expr* result = replace_bindings(rule->replacement, env);
            env_free(env);
            return rule_result(rule, result, def->symbol_name);
        }
        env_free(env);
        rule = rule->next;
    }
    return NULL;
}

#if UP_VALUES
/* ============================================================
 * UpValues (trial). Storage reuses the Rule list, add_rule, the dispatch
 * pre-filter and the matcher -- an UpValue differs from a DownValue only in
 * WHERE it is keyed (an argument symbol) and WHEN it fires (before the enclosing
 * head's DownValues). See symtab.h and apply_up_values below.
 * ============================================================ */

static size_t rule_list_len(const Rule* r) {
    size_t n = 0;
    while (r) { n++; r = r->next; }
    return n;
}

void symtab_add_up_value(const char* symbol_name, Expr* pattern, Expr* replacement) {
    symtab_add_up_value_as(symbol_name, pattern, replacement, true);
}

void symtab_add_up_value_as(const char* symbol_name, Expr* pattern, Expr* replacement,
                            bool delayed) {
    SymbolDef* def = symtab_get_def(symbol_name);
    /* add_rule either replaces a same-LHS rule in place (no count change) or
     * links a new node (+1). Measure the delta so symtab_up_value_count tracks
     * the real list length rather than assuming an insert. */
    size_t before = rule_list_len(def->up_values);
    add_rule(&def->up_values, pattern, replacement, delayed);
    symtab_up_value_count += rule_list_len(def->up_values) - before;
    /* An UpValue changes how an ENCLOSING head evaluates -- advance the rule
     * epoch so GROUND nodes re-check (add_rule already bumped the eval clock). */
    eval_rule_epoch_mark();
}

Rule* symtab_get_up_values(const char* symbol_name) {
    SymbolDef* d = node_find(symbol_name);
    return d ? d->up_values : NULL;
}

/* Scan one symbol's UpValue list against the enclosing call `expr`, with the same
 * arity/first-arg-head dispatch pre-filter and OptionsPattern injection as
 * apply_down_values_def. Returns a freshly-owned rewrite on first match, else
 * NULL. `expr` is borrowed. */
static Expr* apply_up_values_for_def(SymbolDef* def, Expr* expr, bool orderless) {
    if (!def || !def->up_values) return NULL;

    const int32_t input_arity = (int32_t)expr->data.function.arg_count;
    const char* input_first_head =
        (input_arity > 0) ? input_arg_head_canon(expr->data.function.args[0]) : NULL;
    /* Enclosing head Flat => the upvalue may fire on a longer operand sequence
     * than the pattern specifies (surplus -> leftover; see match_flat_leftover),
     * so the arity key is a lower bound and the matched run can sit anywhere. */
    const bool enclosing_flat =
        expr->data.function.head->type == EXPR_SYMBOL &&
        (get_attributes(expr->data.function.head->data.symbol.name) & ATTR_FLAT) != 0;

    for (Rule* rule = def->up_values; rule; rule = rule->next) {
        if (rule->dispatch_arity != -1 && rule->dispatch_arity != input_arity &&
            !(enclosing_flat && input_arity >= rule->dispatch_arity))
            continue;
        /* The first-arg-head filter assumes the pattern's first element aligns
         * positionally with the input's. That holds for a fixed head, but NOT for
         * an Orderless one: the element matching the pattern's (literal) first can
         * sit anywhere (e.g. `x + y_` keyed on x vs the input `1 + x`), so the
         * filter would wrongly reject. Skip it there -- the matcher still checks
         * fully; the arity filter remains sound since sorting preserves count.
         * A Flat (leftover) match likewise need not start at arg 0, so skip it. */
        if (!orderless && !enclosing_flat && rule->first_arg_head_canon && input_first_head &&
            rule->first_arg_head_canon != input_first_head)
            continue;

        MatchEnv* env = env_new();
        /* Phase 1 inline (unchanged fast path); failed match -> Flat leftover. */
        const char* prev_head = match_set_rule_head(
            (expr->data.function.head && expr->data.function.head->type == EXPR_SYMBOL)
                ? expr->data.function.head->data.symbol.name : def->symbol_name);
        Expr* result = match(expr, rule->pattern, env)
            ? replace_bindings(rule->replacement, env)
            : match_flat_leftover(expr, rule->pattern, rule->replacement, env);
        match_set_rule_head(prev_head);
        if (result && rule->scoped_cond) {
            Expr* v = try_scoped_condition(result);
            expr_free(result);
            env_free(env);
            if (v) return v;
            continue;
        }
        if (result) {
            Expr* opts = env_get(env, "$OptionsPattern$");
            if (opts) {
                /* For an upvalue, OptionsPattern[]'s enclosing head is `expr`'s
                 * head (the outer call), not the tag `def` we are keyed on. */
                const char* head_sym = def->symbol_name;
                Expr* eh = expr->data.function.head;
                if (eh && eh->type == EXPR_SYMBOL) head_sym = eh->data.symbol.name;
                Expr* optpat_head = env_get(env, "$OptionsPatternHead$");
                if (optpat_head && optpat_head->type == EXPR_SYMBOL)
                    head_sym = optpat_head->data.symbol.name;
                Expr* resolved = optionvalue_inject_context(result, head_sym, opts);
                expr_free(result);
                result = resolved;
            }
            env_free(env);
            return rule_result(rule, result, def->symbol_name);
        }
        env_free(env);
    }
    return NULL;
}

static Expr* apply_up_values_core(Expr* expr) {
    /* Cheap global early-out: with no UpValues defined anywhere, every function
     * call pays just this load+branch before DownValue dispatch. */
    if (symtab_up_value_count == 0) return NULL;
    if (!expr || expr->type != EXPR_FUNCTION) return NULL;

    /* Collect the distinct level-one candidate defs: for each argument, the arg
     * itself when it is a symbol, else its head when it is a symbol-headed call.
     * Dedupe on the (stable) def pointer so a symbol occurring several times is
     * consulted once, left-to-right. Lookups are NON-materializing: a symbol that
     * carries an UpValue already has a def (install went through symtab_get_def),
     * so node_find finds it, while a bare argument symbol is never promoted into
     * Names[] merely by appearing inside an evaluated call. */
    /* Fixed-size candidate buffer. If a single call somehow has more than 32
     * DISTINCT upvalue-bearing operand symbols at level one, the scan stops here
     * and later operands are silently ignored (their upvalues won't fire this
     * pass). This is a safety bound, not a semantic limit: 32 distinct
     * upvalue-carrying symbols in one call is implausible, and the alternative
     * (a heap allocation on every function-call evaluation) would tax the hot
     * path for a case that does not occur. Revisit only if it ever bites. */
    enum { UP_MAX_CANDIDATES = 32 };
    SymbolDef* cand[UP_MAX_CANDIDATES];
    int ncand = 0;

    size_t argc = expr->data.function.arg_count;
    for (size_t i = 0; i < argc && ncand < UP_MAX_CANDIDATES; i++) {
        Expr* a = expr->data.function.args[i];
        Expr* sym_node = NULL;
        if (a->type == EXPR_SYMBOL) {
            sym_node = a;
        } else if (a->type == EXPR_FUNCTION && a->data.function.head) {
            /* g[...] -> g; a compound head g[x][y] -> its innermost symbol g,
             * matching how the definition side keys the rule. */
            Expr* h = a->data.function.head;
            while (h->type == EXPR_FUNCTION && h->data.function.head) h = h->data.function.head;
            if (h->type == EXPR_SYMBOL) sym_node = h;
        }
        if (!sym_node) continue;

        SymbolDef* d = sym_node->data.symbol.def;          /* cached resolve */
        if (!d) d = node_find(sym_node->data.symbol.name); /* non-materializing */
        if (!d || !d->up_values) continue;

        int seen = 0;
        for (int k = 0; k < ncand; k++) if (cand[k] == d) { seen = 1; break; }
        if (!seen) cand[ncand++] = d;
    }
    if (ncand == 0) return NULL;

    /* WMA applies Orderless before upvalues. Plus/Times are excluded from the
     * evaluator's generic Orderless sort (they sort inside their builtin, which runs
     * AFTER this hook), so an upvalue on them would otherwise bind to the unsorted
     * argument order. Canonicalize now -- but only when an upvalue can actually fire
     * (ncand>0) on an Orderless head, so ordinary arithmetic is untouched. The
     * candidate SET is order-independent, so sorting after collection is safe.
     *
     * NOTE: this reorders `expr`'s argument array IN PLACE and does so even when
     * this function ultimately returns NULL (no upvalue fired). There is no
     * `*changed` out-parameter, so the caller is NOT told the order changed. This
     * is safe only because the single caller is the evaluator's pre-DownValue hook
     * and the Plus/Times builtin re-sorts its args immediately afterwards; do not
     * rely on `expr` keeping its original argument order across this call. */
    bool orderless = expr->data.function.head->type == EXPR_SYMBOL &&
        (get_attributes(expr->data.function.head->data.symbol.name) & ATTR_ORDERLESS) != 0;
    if (orderless && expr->data.function.arg_count > 1) {
        expr_orderless_sort(expr->data.function.args, expr->data.function.arg_count);
        expr_invalidate_hash(expr);
    }

    for (int k = 0; k < ncand; k++) {
        Expr* r = apply_up_values_for_def(cand[k], expr, orderless);
        if (r) return r;
    }
    return NULL;
}

/* The public entry: an N-value -- an UpValue matching N[x] or N[x, p], as
 * x /: N[x] = v or N[x] = v store it -- gives its value numericalized,
 * N[v, p], as the reference's NValues do (N[z4] = 5 makes N[z4] 5., not 5). */
Expr* apply_up_values(Expr* expr) {
    Expr* r = apply_up_values_core(expr);
    if (r && expr->data.function.head->type == EXPR_SYMBOL &&
        expr->data.function.head->data.symbol.name == SYM_N &&
        !(r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL &&
          r->data.function.head->data.symbol.name == SYM_N)) {
        size_t n = expr->data.function.arg_count;
        Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
        a[0] = r;
        for (size_t i = 1; i < n; i++) a[i] = expr_copy(expr->data.function.args[i]);
        r = expr_new_function(expr_new_symbol(SYM_N), a, n ? n : 1);
        free(a);
    }
    return r;
}

bool symtab_remove_matching_up_value(const char* symbol_name, const Expr* lhs) {
    if (!symbol_name || !lhs) return false;

    /* Same canonicalize + alpha-normalize as symtab_remove_matching_rule so a
     * renamed pattern variable still compares equal. */
    Expr* canon = pattern_canonicalize(expr_copy((Expr*)lhs));
    Expr* key = pattern_alpha_normalize(canon);
    expr_free(canon);

    SymbolDef* def = symtab_get_def(symbol_name);
    Rule* prev = NULL;
    bool removed = false;
    for (Rule* curr = def->up_values; curr; prev = curr, curr = curr->next) {
        Expr* curr_key = pattern_alpha_normalize(curr->pattern);
        int m = expr_eq(curr_key, key);
        expr_free(curr_key);
        if (!m) continue;

        if (prev) prev->next = curr->next;
        else def->up_values = curr->next;
        expr_free(curr->pattern);
        expr_free(curr->replacement);
        free(curr);
        if (symtab_up_value_count) symtab_up_value_count--;
        eval_clock_bump();
        removed = true;
        break;
    }

    expr_free(key);
    return removed;
}

/* Shared backend for OwnValues[f]=list / DownValues[f]=list / UpValues[f]=list.
 * which: 0=own, 1=down, 2=up, 3=format, 4=N. Frees the existing target list, then installs each
 * Rule/RuleDelayed[lhs,rhs] element (stripping an outer HoldPattern on lhs).
 * Malformed elements are skipped. Order follows add_rule's specificity sort
 * (insertion order the tie-break) -- identical to how a sequence of individual
 * assignments would land. */
static void symtab_set_values_common(const char* symbol_name, const Expr* list,
                                     int which) {
    SymbolDef* def = symtab_get_def(symbol_name);
    Rule** target = (which == 0) ? &def->own_values
                  : (which == 1) ? &def->down_values
                  : (which == 2) ? &def->up_values
                  : (which == 3) ? &def->format_values
                                 : &def->n_values;

    if (which == 1) litidx_drop(def);
    if (which == 2)      free_up_value_list(*target);
    else if (which == 3) free_format_value_list(*target);
    else if (which == 4) free_n_value_list(*target);
    else                 free_rule_list(*target);
    *target = NULL;

    if (list && list->type == EXPR_FUNCTION &&
        list->data.function.head->type == EXPR_SYMBOL &&
        list->data.function.head->data.symbol.name == SYM_List) {
        for (size_t i = 0; i < list->data.function.arg_count; i++) {
            Expr* el = list->data.function.args[i];
            if (!(el->type == EXPR_FUNCTION &&
                  el->data.function.head->type == EXPR_SYMBOL &&
                  (el->data.function.head->data.symbol.name == SYM_Rule ||
                   el->data.function.head->data.symbol.name == SYM_RuleDelayed) &&
                  el->data.function.arg_count == 2))
                continue;   /* skip a non-rule element */
            Expr* lhs = el->data.function.args[0];
            Expr* rhs = el->data.function.args[1];
            bool delayed = el->data.function.head->data.symbol.name == SYM_RuleDelayed;
            if (lhs->type == EXPR_FUNCTION &&
                lhs->data.function.head->type == EXPR_SYMBOL &&
                lhs->data.function.head->data.symbol.name == SYM_HoldPattern &&
                lhs->data.function.arg_count == 1)
                lhs = lhs->data.function.args[0];   /* unwrap HoldPattern */

            if (which >= 2) {
                size_t before = rule_list_len(*target);
                add_rule(target, lhs, rhs, delayed);
                size_t added = rule_list_len(*target) - before;
                if (which == 2)      symtab_up_value_count += added;
                else if (which == 3) symtab_format_value_count += added;
                else                 symtab_n_value_count += added;
            } else {
                add_rule_in(target, lhs, rhs, delayed, which == 1 ? def : NULL, NULL);
            }
        }
    }
    eval_clock_bump();
    eval_rule_epoch_mark();
}

void symtab_set_own_values(const char* symbol_name, const Expr* list) {
    symtab_set_values_common(symbol_name, list, 0);
}
void symtab_set_down_values(const char* symbol_name, const Expr* list) {
    symtab_set_values_common(symbol_name, list, 1);
}
void symtab_set_up_values(const char* symbol_name, const Expr* list) {
    symtab_set_values_common(symbol_name, list, 2);
}
void symtab_set_format_values(const char* symbol_name, const Expr* list) {
    symtab_set_values_common(symbol_name, list, 3);
}
void symtab_set_n_values(const char* symbol_name, const Expr* list) {
    symtab_set_values_common(symbol_name, list, 4);
}
#endif /* UP_VALUES */

void symtab_add_format_value(const char* symbol_name, Expr* pattern, Expr* replacement,
                             bool delayed) {
    SymbolDef* def = symtab_get_def(symbol_name);
    /* add_rule replaces a same-LHS rule in place or links a new one: count. */
    size_t before = 0, after = 0;
    for (const Rule* r = def->format_values; r; r = r->next) before++;
    add_rule(&def->format_values, pattern, replacement, delayed);
    for (const Rule* r = def->format_values; r; r = r->next) after++;
    symtab_format_value_count += after - before;
}

Rule* symtab_get_format_values(const char* symbol_name) {
    SymbolDef* d = node_find(symbol_name);
    return d ? d->format_values : NULL;
}

void symtab_add_n_value(const char* symbol_name, Expr* pattern, Expr* replacement,
                        bool delayed) {
    SymbolDef* def = symtab_get_def(symbol_name);
    size_t before = 0, after = 0;
    for (const Rule* r = def->n_values; r; r = r->next) before++;
    add_rule(&def->n_values, pattern, replacement, delayed);
    for (const Rule* r = def->n_values; r; r = r->next) after++;
    symtab_n_value_count += after - before;
}

Rule* symtab_get_n_values(const char* symbol_name) {
    SymbolDef* d = node_find(symbol_name);
    return d ? d->n_values : NULL;
}
