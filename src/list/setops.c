#include "list_common.h"
#include "setops.h"
#include "assoc.h"
#include "ndarray.h"    /* is_ndarray */
#include "ndarray_internal.h"  /* nd_sort_i64_asc — the packed set-op fast path */
#include "ndreduce.h"   /* ndred_tally — Tally's packed-buffer fast path */

/* ---------------------------------------------------------------------------
 *  Packed integer set operations.
 *
 *  Union / Intersection / Complement are written generically above: expr_copy
 *  every element, qsort through a function pointer calling expr_compare, then
 *  dedup with expr_eq. On a packed integer list that is one Expr allocation and
 *  two indirect calls per element -- Union of 10^6 integers cost 745 ms, against
 *  ~60 ms for np.unique on the same data -- and it is the shape that matters,
 *  because a set operation over integer LABELS is what a graph traversal and a
 *  k-mer count are made of. Breadth-first search is literally
 *      Complement[Union[Flatten[adj[[frontier]]]], visited]
 *  once per level.
 *
 *  Restricted to int64 on purpose. Two Reals that compare equal but print
 *  differently (0. and -0.) would let the buffer path keep a different
 *  representative than the List path keeps, and a set operation whose ANSWER
 *  depends on which of two equal elements survived is not a fast path, it is a
 *  second implementation. Integers have no such pair.
 * ------------------------------------------------------------------------- */

/* True when `e` is a rank-1 packed list of int64 -- the only shape below. */
static bool setop_i64(const Expr* e) {
    return is_packed_list(e) && e->data.ndarray.rank == 1 &&
           e->data.ndarray.dtype == NDT_INT64;
}

/* Sorted, deduplicated copy of `a`'s buffer. Caller owns `*out`. */
static bool setop_sorted_unique(const Expr* a, int64_t** out, size_t* n_out) {
    size_t n = (size_t)a->data.ndarray.dims[0];
    int64_t* v = malloc(sizeof(int64_t) * (n ? n : 1));
    if (!v) return false;
    memcpy(v, a->data.ndarray.data, sizeof(int64_t) * n);
    nd_sort_i64_asc(v, n);
    size_t m = 0;
    for (size_t i = 0; i < n; i++)
        if (m == 0 || v[i] != v[m - 1]) v[m++] = v[i];
    *out = v; *n_out = m;
    return true;
}

/* Wrap an owned int64 buffer as a packed list shaped like `src`. Frees `v` and
 * returns NULL on allocation failure, so callers can degrade. */
static Expr* setop_emit(const Expr* src, int64_t* v, size_t m) {
    int64_t dims[1] = { (int64_t)m };
    if (m == 0) { free(v); return expr_new_function(expr_new_symbol(SYM_List), NULL, 0); }
    int64_t* buf = malloc(sizeof(int64_t) * m);
    if (!buf) { free(v); return NULL; }
    memcpy(buf, v, sizeof(int64_t) * m);
    free(v);
    return expr_new_ndarray_like(src, 1, dims, buf, NDT_INT64);
}

/* Union / Intersection / Complement over `nl` packed int64 lists. `mode` is 0,
 * 1, 2 respectively. Returns NULL when the shape is outside the fast domain. */
static Expr* setop_packed(Expr** args, size_t nl, int mode) {
    for (size_t i = 0; i < nl; i++) if (!setop_i64(args[i])) return NULL;
    if (nl == 0) return NULL;

    int64_t* acc; size_t accn;
    if (!setop_sorted_unique(args[0], &acc, &accn)) return NULL;

    for (size_t i = 1; i < nl; i++) {
        int64_t* b; size_t bn;
        if (!setop_sorted_unique(args[i], &b, &bn)) { free(acc); return NULL; }
        size_t p = 0, q = 0, m = 0;
        if (mode == 0) {                              /* Union: merge */
            size_t cap = accn + bn;
            int64_t* merged = malloc(sizeof(int64_t) * (cap ? cap : 1));
            if (!merged) { free(acc); free(b); return NULL; }
            while (p < accn || q < bn) {
                int64_t v;
                if (q >= bn)                 v = acc[p++];
                else if (p >= accn)          v = b[q++];
                else if (acc[p] < b[q])      v = acc[p++];
                else if (b[q] < acc[p])      v = b[q++];
                else                       { v = acc[p++]; q++; }
                if (m == 0 || v != merged[m - 1]) merged[m++] = v;
            }
            free(acc); free(b); acc = merged; accn = m;
        } else if (mode == 1) {                       /* Intersection */
            while (p < accn && q < bn) {
                if (acc[p] < b[q]) p++;
                else if (b[q] < acc[p]) q++;
                else { acc[m++] = acc[p++]; q++; }
            }
            free(b); accn = m;
        } else {                                      /* Complement: acc \ b */
            while (p < accn) {
                while (q < bn && b[q] < acc[p]) q++;
                if (q < bn && b[q] == acc[p]) p++;
                else acc[m++] = acc[p++];
            }
            free(b); accn = m;
        }
    }
    return setop_emit(args[0], acc, accn);
}

/* True when any of args[0..n) is an ndarray -- the guard every generic set
 * operation needs before it runs, because the code below tests
 * `type != EXPR_FUNCTION` and an NDArray is not one, so it would answer with
 * the array itself or with {}. */
static bool setop_any_nd(Expr** args, size_t n) {
    for (size_t i = 0; i < n; i++) if (is_ndarray(args[i])) return true;
    return false;
}

typedef struct HashNode {
    Expr* key;
    size_t index; // Original index or position in unique_elems
    struct HashNode* next;
} HashNode;

typedef struct {
    HashNode** buckets;
    size_t size;
} HashTable;

static HashTable* ht_create(size_t size) {
    HashTable* ht = malloc(sizeof(HashTable));
    ht->size = size;
    ht->buckets = calloc(size, sizeof(HashNode*));
    return ht;
}

static void ht_free(HashTable* ht, bool free_keys) {
    for (size_t i = 0; i < ht->size; i++) {
        HashNode* node = ht->buckets[i];
        while (node) {
            HashNode* next = node->next;
            if (free_keys) expr_free(node->key);
            free(node);
            node = next;
        }
    }
    free(ht->buckets);
    free(ht);
}

static HashNode* ht_find(HashTable* ht, Expr* key) {
    uint64_t h = expr_hash(key);
    size_t bucket = (size_t)(h % ht->size);
    HashNode* node = ht->buckets[bucket];
    while (node) {
        if (expr_eq(node->key, key)) return node;
        node = node->next;
    }
    return NULL;
}

static void ht_insert(HashTable* ht, Expr* key, size_t index) {
    uint64_t h = expr_hash(key);
    size_t bucket = (size_t)(h % ht->size);
    HashNode* node = malloc(sizeof(HashNode));
    node->key = key;
    node->index = index;
    node->next = ht->buckets[bucket];
    ht->buckets[bucket] = node;
}

static int compare_expr_ptrs(const void* a, const void* b) {
    Expr* ea = *(Expr**)a;
    Expr* eb = *(Expr**)b;
    return expr_compare(ea, eb);
}

Expr* builtin_union(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;
    
    // Find options
    Expr* same_test = NULL;
    size_t last_arg = res->data.function.arg_count;
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        Expr* arg = res->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol.name == SYM_Rule &&
            arg->data.function.arg_count == 2 &&
            arg->data.function.args[0]->type == EXPR_SYMBOL &&
            arg->data.function.args[0]->data.symbol.name == SYM_SameTest) {
            same_test = arg->data.function.args[1];
            if (i < last_arg) last_arg = i;
        }
    }
    
    if (last_arg == 0) return NULL;

    /* Packed operands: take the buffer merge, or degrade. Never fall through --
     * the test just below is `type != EXPR_FUNCTION`, and an NDArray is not one,
     * so Union[packedList] would answer with the array itself unchanged. */
    if (setop_any_nd(res->data.function.args, last_arg)) {
        if (!same_test) {
            Expr* fast = setop_packed(res->data.function.args, last_arg, 0);
            if (fast) return fast;
        }
        return ndarray_delist_and_reeval(res);
    }

    // Check if first arg is a function
    Expr* first_list = res->data.function.args[0];
    if (first_list->type != EXPR_FUNCTION) return expr_copy(first_list);
    
    Expr* common_head = first_list->data.function.head;
    
    // Total count of elements
    size_t total_count = 0;
    for (size_t i = 0; i < last_arg; i++) {
        Expr* arg = res->data.function.args[i];
        if (arg->type != EXPR_FUNCTION || !expr_eq(arg->data.function.head, common_head)) {
            // Heads must match
            return NULL;
        }
        total_count += arg->data.function.arg_count;
    }
    
    if (total_count == 0) return expr_copy(first_list);
    
    Expr** all_args = malloc(sizeof(Expr*) * total_count);
    size_t idx = 0;
    for (size_t i = 0; i < last_arg; i++) {
        Expr* arg = res->data.function.args[i];
        for (size_t j = 0; j < arg->data.function.arg_count; j++) {
            all_args[idx++] = expr_copy(arg->data.function.args[j]);
        }
    }
    
    // Sort elements
    qsort(all_args, total_count, sizeof(Expr*), compare_expr_ptrs);
    
    // Remove duplicates
    Expr** unique_args = malloc(sizeof(Expr*) * total_count);
    size_t unique_count = 0;
    
    if (total_count > 0) {
        unique_args[unique_count++] = all_args[0];
        for (size_t i = 1; i < total_count; i++) {
            bool is_dup = false;
            if (same_test == NULL) {
                if (expr_eq(all_args[i], unique_args[unique_count - 1])) {
                    is_dup = true;
                }
            } else {
                Expr* call_args[2] = { expr_copy(all_args[i]), expr_copy(unique_args[unique_count - 1]) };
                Expr* call = expr_new_function(expr_copy(same_test), call_args, 2);
                Expr* eval_res = evaluate(call);
                if (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol.name == SYM_True) {
                    is_dup = true;
                }
                expr_free(eval_res);
                expr_free(call);
            }
            
            if (is_dup) {
                expr_free(all_args[i]);
            } else {
                unique_args[unique_count++] = all_args[i];
            }
        }
    }
    
    free(all_args);
    
    Expr* result = expr_new_function(expr_copy(common_head), unique_args, unique_count);
    if (unique_args) free(unique_args);

    return result;
}

/* Evaluate test[a, b] and report whether it yields the symbol True. Used by
 * Intersection's SameTest path (a and b are borrowed; the call is built from
 * fresh copies and freed here). */
static bool same_test_equal(Expr* test, Expr* a, Expr* b) {
    Expr* call_args[2] = { expr_copy(a), expr_copy(b) };
    Expr* call = expr_new_function(expr_copy(test), call_args, 2);
    Expr* r = evaluate(call);
    bool eq = (r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_True);
    expr_free(r);
    expr_free(call);
    return eq;
}

/* Intersection[l1, l2, ...] gives the sorted list of elements common to every
 * li, deduplicated, using the head of the first argument (need not be List).
 *
 * Default comparison is canonical structural equality, computed in O(total)
 * with the file-local hash set: deduplicate the first operand, then keep only
 * the candidates that survive a membership test against each later operand.
 *
 * SameTest -> f switches to an O(n^2) path using f[a,b]===True as the
 * equivalence relation. Wolfram keeps the canonically-greatest member of each
 * equivalence class as its representative, so we sort the first operand
 * ascending and, on a class collision, replace the stored representative with
 * the later (hence greater) element. */
Expr* builtin_intersection(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1)
        return builtin_arg_error("Intersection",
            (res->type == EXPR_FUNCTION) ? res->data.function.arg_count : 0,
            1, SIZE_MAX);

    /* Locate a trailing SameTest option, mirroring builtin_union. */
    Expr* same_test = NULL;
    size_t last_arg = res->data.function.arg_count;
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        Expr* arg = res->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol.name == SYM_Rule &&
            arg->data.function.arg_count == 2 &&
            arg->data.function.args[0]->type == EXPR_SYMBOL &&
            arg->data.function.args[0]->data.symbol.name == SYM_SameTest) {
            same_test = arg->data.function.args[1];
            if (i < last_arg) last_arg = i;
        }
    }

    if (last_arg == 0) return NULL;   /* only an option, no list operands */

    if (setop_any_nd(res->data.function.args, last_arg)) {
        if (!same_test) {
            Expr* fast = setop_packed(res->data.function.args, last_arg, 1);
            if (fast) return fast;
        }
        return ndarray_delist_and_reeval(res);
    }

    Expr* first = res->data.function.args[0];
    if (first->type != EXPR_FUNCTION) return expr_copy(first);
    Expr* common_head = first->data.function.head;

    /* Every operand must share the head of the first. */
    for (size_t i = 0; i < last_arg; i++) {
        Expr* arg = res->data.function.args[i];
        if (arg->type != EXPR_FUNCTION || !expr_eq(arg->data.function.head, common_head))
            return NULL;
    }

    size_t first_count = first->data.function.arg_count;
    Expr** cand = malloc(sizeof(Expr*) * (first_count ? first_count : 1));
    size_t ncand = 0;

    if (same_test == NULL) {
        /* --- Default path: hash-based, O(total). --- */
        HashTable* seen = ht_create(first_count * 2 + 1);
        for (size_t j = 0; j < first_count; j++) {
            Expr* e = first->data.function.args[j];
            if (!ht_find(seen, e)) {
                Expr* c = expr_copy(e);
                cand[ncand++] = c;
                ht_insert(seen, c, 0);
            }
        }
        ht_free(seen, false);

        for (size_t i = 1; i < last_arg && ncand > 0; i++) {
            Expr* arg = res->data.function.args[i];
            size_t m = arg->data.function.arg_count;
            HashTable* h = ht_create(m * 2 + 1);
            for (size_t j = 0; j < m; j++)
                if (!ht_find(h, arg->data.function.args[j]))
                    ht_insert(h, arg->data.function.args[j], 0);
            size_t w = 0;
            for (size_t k = 0; k < ncand; k++) {
                if (ht_find(h, cand[k])) cand[w++] = cand[k];
                else expr_free(cand[k]);
            }
            ncand = w;
            ht_free(h, false);
        }
    } else {
        /* --- SameTest path: O(n^2), greatest representative per class. --- */
        Expr** sorted0 = malloc(sizeof(Expr*) * (first_count ? first_count : 1));
        for (size_t j = 0; j < first_count; j++)
            sorted0[j] = first->data.function.args[j];   /* borrowed */
        qsort(sorted0, first_count, sizeof(Expr*), compare_expr_ptrs);
        for (size_t j = 0; j < first_count; j++) {
            Expr* e = sorted0[j];
            size_t match = SIZE_MAX;
            for (size_t k = 0; k < ncand; k++)
                if (same_test_equal(same_test, e, cand[k])) { match = k; break; }
            if (match == SIZE_MAX) {
                cand[ncand++] = expr_copy(e);
            } else {
                /* e is canonically >= cand[match]; keep the greater one. */
                expr_free(cand[match]);
                cand[match] = expr_copy(e);
            }
        }
        free(sorted0);

        for (size_t i = 1; i < last_arg && ncand > 0; i++) {
            Expr* arg = res->data.function.args[i];
            size_t m = arg->data.function.arg_count;
            size_t w = 0;
            for (size_t k = 0; k < ncand; k++) {
                bool present = false;
                for (size_t j = 0; j < m; j++)
                    if (same_test_equal(same_test, cand[k], arg->data.function.args[j])) {
                        present = true; break;
                    }
                if (present) cand[w++] = cand[k];
                else expr_free(cand[k]);
            }
            ncand = w;
        }
    }

    qsort(cand, ncand, sizeof(Expr*), compare_expr_ptrs);
    Expr* result = expr_new_function(expr_copy(common_head), cand, ncand);
    free(cand);
    return result;
}

/* Complement[eall, e1, e2, ...] gives the sorted list of the distinct elements
 * of eall that appear in none of the ei, using the head of the first argument
 * (need not be List). It is the structural twin of Intersection with the
 * membership test inverted: a candidate survives when it is ABSENT from every
 * later operand, rather than present in all of them. Unlike Intersection,
 * Complement is order-sensitive in its first argument, hence not Flat.
 *
 * Default comparison is canonical structural equality, computed in O(total)
 * with the file-local hash set: deduplicate the first operand, then drop the
 * candidates that hit any later operand.
 *
 * SameTest -> f switches to an O(n^2) path using f[a,b]===True as the
 * equivalence relation. Unlike Intersection (which keeps the canonically-
 * greatest member of each class), Wolfram's Complement keeps the canonically-
 * smallest member as the representative, so we sort the first operand ascending
 * and keep the first element seen for each class. SameTest -> Automatic is
 * treated as the default (hash) path. */
Expr* builtin_complement(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1)
        return builtin_arg_error("Complement",
            (res->type == EXPR_FUNCTION) ? res->data.function.arg_count : 0,
            1, SIZE_MAX);

    /* Locate a trailing SameTest option, mirroring builtin_intersection. */
    Expr* same_test = NULL;
    size_t last_arg = res->data.function.arg_count;
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        Expr* arg = res->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol.name == SYM_Rule &&
            arg->data.function.arg_count == 2 &&
            arg->data.function.args[0]->type == EXPR_SYMBOL &&
            arg->data.function.args[0]->data.symbol.name == SYM_SameTest) {
            same_test = arg->data.function.args[1];
            if (i < last_arg) last_arg = i;
        }
    }

    /* SameTest -> Automatic means "use the default comparison". */
    if (same_test != NULL && same_test->type == EXPR_SYMBOL &&
        same_test->data.symbol.name == SYM_Automatic)
        same_test = NULL;

    if (last_arg == 0) return NULL;   /* only an option, no list operands */

    if (setop_any_nd(res->data.function.args, last_arg)) {
        if (!same_test) {
            Expr* fast = setop_packed(res->data.function.args, last_arg, 2);
            if (fast) return fast;
        }
        return ndarray_delist_and_reeval(res);
    }

    Expr* first = res->data.function.args[0];
    if (first->type != EXPR_FUNCTION) return expr_copy(first);
    Expr* common_head = first->data.function.head;

    /* Every operand must share the head of the first. */
    for (size_t i = 0; i < last_arg; i++) {
        Expr* arg = res->data.function.args[i];
        if (arg->type != EXPR_FUNCTION || !expr_eq(arg->data.function.head, common_head))
            return NULL;
    }

    size_t first_count = first->data.function.arg_count;
    Expr** cand = malloc(sizeof(Expr*) * (first_count ? first_count : 1));
    size_t ncand = 0;

    if (same_test == NULL) {
        /* --- Default path: hash-based, O(total). --- */
        HashTable* seen = ht_create(first_count * 2 + 1);
        for (size_t j = 0; j < first_count; j++) {
            Expr* e = first->data.function.args[j];
            if (!ht_find(seen, e)) {
                Expr* c = expr_copy(e);
                cand[ncand++] = c;
                ht_insert(seen, c, 0);
            }
        }
        ht_free(seen, false);

        for (size_t i = 1; i < last_arg && ncand > 0; i++) {
            Expr* arg = res->data.function.args[i];
            size_t m = arg->data.function.arg_count;
            HashTable* h = ht_create(m * 2 + 1);
            for (size_t j = 0; j < m; j++)
                if (!ht_find(h, arg->data.function.args[j]))
                    ht_insert(h, arg->data.function.args[j], 0);
            size_t w = 0;
            for (size_t k = 0; k < ncand; k++) {
                if (!ht_find(h, cand[k])) cand[w++] = cand[k];
                else expr_free(cand[k]);
            }
            ncand = w;
            ht_free(h, false);
        }
    } else {
        /* --- SameTest path: O(n^2), smallest representative per class. --- */
        Expr** sorted0 = malloc(sizeof(Expr*) * (first_count ? first_count : 1));
        for (size_t j = 0; j < first_count; j++)
            sorted0[j] = first->data.function.args[j];   /* borrowed */
        qsort(sorted0, first_count, sizeof(Expr*), compare_expr_ptrs);
        for (size_t j = 0; j < first_count; j++) {
            Expr* e = sorted0[j];
            bool seen = false;
            for (size_t k = 0; k < ncand; k++)
                if (same_test_equal(same_test, e, cand[k])) { seen = true; break; }
            /* Sorted ascending: the first member seen for each class is its
             * canonically-smallest element, which we keep as representative. */
            if (!seen) cand[ncand++] = expr_copy(e);
        }
        free(sorted0);

        for (size_t i = 1; i < last_arg && ncand > 0; i++) {
            Expr* arg = res->data.function.args[i];
            size_t m = arg->data.function.arg_count;
            size_t w = 0;
            for (size_t k = 0; k < ncand; k++) {
                bool present = false;
                for (size_t j = 0; j < m; j++)
                    if (same_test_equal(same_test, cand[k], arg->data.function.args[j])) {
                        present = true; break;
                    }
                if (!present) cand[w++] = cand[k];
                else expr_free(cand[k]);
            }
            ncand = w;
        }
    }

    qsort(cand, ncand, sizeof(Expr*), compare_expr_ptrs);
    Expr* result = expr_new_function(expr_copy(common_head), cand, ncand);
    free(cand);
    return result;
}

Expr* builtin_tally(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    
    Expr* list = res->data.function.args[0];
    Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;

    /* Tally[assoc] (and Tally[assoc, test]) tallies the association's values. */
    if (is_association(list)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    /* A packed buffer hashes as machine words. A user-supplied test has to see
     * the elements as expressions, so that form is handed back as a List call --
     * NOT left to fall through, because everything below indexes `list` as an
     * EXPR_FUNCTION and an NDArray is not one: it would answer {}. */
    if (is_ndarray(list))
        return test ? ndarray_delist_and_reeval(res) : ndred_tally(res);

    if (list->type != EXPR_FUNCTION) return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);

    size_t count = list->data.function.arg_count;
    if (count == 0) return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);

    Expr** unique_elems = malloc(sizeof(Expr*) * count);
    int64_t* multiplicities = malloc(sizeof(int64_t) * count);
    size_t unique_count = 0;

    if (test == NULL) {
        HashTable* ht = ht_create(count * 2 + 1);
        for (size_t i = 0; i < count; i++) {
            Expr* elem = list->data.function.args[i];
            HashNode* node = ht_find(ht, elem);
            if (node) {
                multiplicities[node->index]++;
            } else {
                unique_elems[unique_count] = expr_copy(elem);
                multiplicities[unique_count] = 1;
                ht_insert(ht, unique_elems[unique_count], unique_count);
                unique_count++;
            }
        }
        ht_free(ht, false);
    } else {
        // Fallback to O(N^2) for custom test
        for (size_t i = 0; i < count; i++) {
            Expr* elem = list->data.function.args[i];
            int found_idx = -1;
            for (size_t j = 0; j < unique_count; j++) {
                Expr* call_args[2] = { expr_copy(elem), expr_copy(unique_elems[j]) };
                Expr* call = expr_new_function(expr_copy(test), call_args, 2);
                Expr* eval_res = evaluate(call);
                if (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol.name == SYM_True) {
                    found_idx = (int)j;
                    expr_free(eval_res);
                    expr_free(call);
                    break;
                }
                expr_free(eval_res);
                expr_free(call);
            }
            if (found_idx != -1) {
                multiplicities[found_idx]++;
            } else {
                unique_elems[unique_count] = expr_copy(elem);
                multiplicities[unique_count] = 1;
                unique_count++;
            }
        }
    }
    
    Expr** result_args = malloc(sizeof(Expr*) * unique_count);
    for (size_t i = 0; i < unique_count; i++) {
        Expr** pair_args = malloc(sizeof(Expr*) * 2);
        pair_args[0] = unique_elems[i];
        pair_args[1] = expr_new_integer(multiplicities[i]);
        result_args[i] = expr_new_function(expr_new_symbol(SYM_List), pair_args, 2);
        free(pair_args);
    }
    
    free(unique_elems);
    free(multiplicities);
    
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), result_args, unique_count);
    free(result_args);
    
    return result;
}

Expr* builtin_deleteduplicates(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    
    Expr* list = res->data.function.args[0];
    Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;

    /* DeleteDuplicates[assoc] keeps the first entry for each distinct value,
     * returning an association (Wolfram semantics). Only the default (no custom
     * test) case is hash-indexed here; a custom test falls through unhandled. */
    if (test == NULL && is_association(list)) return assoc_delete_duplicate_values(list);

    if (list->type != EXPR_FUNCTION) return expr_copy(list);
    
    size_t count = list->data.function.arg_count;
    if (count == 0) return expr_copy(list);
    
    Expr** unique_args = malloc(sizeof(Expr*) * count);
    size_t unique_count = 0;

    if (test == NULL) {
        HashTable* ht = ht_create(count * 2 + 1);
        for (size_t i = 0; i < count; i++) {
            Expr* elem = list->data.function.args[i];
            if (!ht_find(ht, elem)) {
                Expr* copy = expr_copy(elem);
                unique_args[unique_count++] = copy;
                ht_insert(ht, copy, 0);
            }
        }
        ht_free(ht, false);
    } else {
        // Fallback to O(N^2) for custom test
        for (size_t i = 0; i < count; i++) {
            Expr* elem = list->data.function.args[i];
            bool is_duplicate = false;
            for (size_t j = 0; j < unique_count; j++) {
                Expr* call_args[2] = { expr_copy(elem), expr_copy(unique_args[j]) };
                Expr* call = expr_new_function(expr_copy(test), call_args, 2);
                Expr* eval_res = evaluate(call);
                if (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol.name == SYM_True) {
                    is_duplicate = true;
                    expr_free(eval_res);
                    expr_free(call);
                    break;
                }
                expr_free(eval_res);
                expr_free(call);
            }
            if (!is_duplicate) {
                unique_args[unique_count++] = expr_copy(elem);
            }
        }
    }
    
    Expr* result = expr_new_function(expr_copy(list->data.function.head), unique_args, unique_count);
    if (unique_args) free(unique_args);

    return result;
}

/* DeleteDuplicatesBy[expr, f] keeps the first element for each distinct value of
 * f[element], preserving order. Over an association f is applied to each value
 * and the surviving entries are returned as an association (keys preserved).
 * The distinct f-values seen so far are compared directly (expr_eq); the count
 * of survivors is typically small, so this stays well within budget. */
Expr* builtin_deleteduplicatesby(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* coll = res->data.function.args[0];
    Expr* f    = res->data.function.args[1];
    bool assoc = is_association(coll);
    if (!assoc && coll->type != EXPR_FUNCTION) return NULL;

    size_t n = coll->data.function.arg_count;
    Expr** kept = malloc(sizeof(Expr*) * (n ? n : 1));   /* surviving elements/rules */
    Expr** keys = malloc(sizeof(Expr*) * (n ? n : 1));   /* their f-values (owned) */
    size_t nkept = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* elem = coll->data.function.args[i];
        Expr* val  = assoc ? elem->data.function.args[1] : elem;  /* value for assoc */
        Expr* fcall = expr_new_function(expr_copy(f), (Expr*[]){ expr_copy(val) }, 1);
        Expr* fk = evaluate(fcall);
        expr_free(fcall);
        bool dup = false;
        for (size_t j = 0; j < nkept; j++)
            if (expr_eq(fk, keys[j])) { dup = true; break; }
        if (dup) { expr_free(fk); continue; }
        keys[nkept] = fk;
        kept[nkept] = expr_copy(elem);
        nkept++;
    }
    Expr* head = assoc ? expr_new_symbol(SYM_Association)
                       : expr_copy(coll->data.function.head);
    Expr* out = expr_new_function(head, kept, nkept);
    for (size_t j = 0; j < nkept; j++) expr_free(keys[j]);
    free(keys); free(kept);
    return out;
}

typedef struct {
    Expr* element;
    int64_t count;
    size_t first_index;
} CommonestItem;

static int compare_commonest_items_desc(const void* a, const void* b) {
    const CommonestItem* item_a = (const CommonestItem*)a;
    const CommonestItem* item_b = (const CommonestItem*)b;
    if (item_a->count != item_b->count) {
        return (item_b->count > item_a->count) ? 1 : -1;
    }
    return (item_a->first_index > item_b->first_index) ? 1 : -1;
}

static int compare_commonest_items_index(const void* a, const void* b) {
    const CommonestItem* item_a = (const CommonestItem*)a;
    const CommonestItem* item_b = (const CommonestItem*)b;
    if (item_a->first_index == item_b->first_index) return 0;
    return (item_a->first_index > item_b->first_index) ? 1 : -1;
}

Expr* builtin_commonest(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    
    Expr* list = res->data.function.args[0];

    /* Commonest[assoc] (and Commonest[assoc, n]) uses the association's values. */
    if (is_association(list)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

    if (list->type != EXPR_FUNCTION) return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);

    size_t count = list->data.function.arg_count;
    if (count == 0) return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);

    Expr* n_arg = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    int64_t n = -1;
    bool n_upto = false;
    if (n_arg) {
        if (n_arg->type == EXPR_INTEGER) {
            n = n_arg->data.integer;
        } else if (n_arg->type == EXPR_FUNCTION && n_arg->data.function.head->type == EXPR_SYMBOL && 
                   n_arg->data.function.head->data.symbol.name == SYM_UpTo && n_arg->data.function.arg_count == 1) {
            if (n_arg->data.function.args[0]->type == EXPR_INTEGER) {
                n = n_arg->data.function.args[0]->data.integer;
                n_upto = true;
            } else return NULL;
        } else return NULL;
    }

    // Tally
    Expr** unique_elems = malloc(sizeof(Expr*) * count);
    int64_t* multiplicities = malloc(sizeof(int64_t) * count);
    size_t unique_count = 0;

    HashTable* ht = ht_create(count * 2 + 1);
    for (size_t i = 0; i < count; i++) {
        Expr* elem = list->data.function.args[i];
        HashNode* node = ht_find(ht, elem);
        if (node) {
            multiplicities[node->index]++;
        } else {
            unique_elems[unique_count] = expr_copy(elem);
            multiplicities[unique_count] = 1;
            ht_insert(ht, unique_elems[unique_count], unique_count);
            unique_count++;
        }
    }
    ht_free(ht, false);

    CommonestItem* items = malloc(sizeof(CommonestItem) * unique_count);
    for (size_t i = 0; i < unique_count; i++) {
        items[i].element = unique_elems[i];
        items[i].count = multiplicities[i];
        items[i].first_index = i;
    }
    free(multiplicities);
    free(unique_elems);

    // Sort by count DESC, first_index ASC
    qsort(items, unique_count, sizeof(CommonestItem), compare_commonest_items_desc);

    size_t target_n;
    if (n == -1) {
        // Just the most common ones (highest count)
        int64_t max_count = items[0].count;
        target_n = 0;
        while (target_n < unique_count && items[target_n].count == max_count) {
            target_n++;
        }
    } else {
        if (n < 0) n = 0;
        if ((size_t)n > unique_count) {
            if (!n_upto) {
                printf("Commonest::dstlms: The requested number of elements %" PRId64 " is greater than the number of distinct elements %zu. Only %zu elements will be returned.\n", n, unique_count, unique_count);
            }
            target_n = unique_count;
        } else {
            target_n = (size_t)n;
        }
    }

    // Sort target_n items by first_index ASC to preserve original order
    if (target_n > 0) {
        qsort(items, target_n, sizeof(CommonestItem), compare_commonest_items_index);
    }

    Expr** result_args = malloc(sizeof(Expr*) * target_n);
    for (size_t i = 0; i < target_n; i++) {
        result_args[i] = items[i].element;
    }
    // Free unused elements
    for (size_t i = target_n; i < unique_count; i++) {
        expr_free(items[i].element);
    }
    free(items);

    return expr_new_function(expr_new_symbol(SYM_List), result_args, target_n);
}
