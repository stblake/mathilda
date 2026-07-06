#include "list_common.h"
#include "setops.h"
#include "assoc.h"

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
            arg->data.function.head->data.symbol == SYM_Rule &&
            arg->data.function.arg_count == 2 &&
            arg->data.function.args[0]->type == EXPR_SYMBOL &&
            arg->data.function.args[0]->data.symbol == SYM_SameTest) {
            same_test = arg->data.function.args[1];
            if (i < last_arg) last_arg = i;
        }
    }
    
    if (last_arg == 0) return NULL;
    
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
                if (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol == SYM_True) {
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

Expr* builtin_tally(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    
    Expr* list = res->data.function.args[0];
    Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;

    /* Tally[assoc] (and Tally[assoc, test]) tallies the association's values. */
    if (is_association(list)) { Expr* r = assoc_apply_over_values(res); if (r) return r; }

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
                if (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol == SYM_True) {
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
                if (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol == SYM_True) {
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
                   n_arg->data.function.head->data.symbol == SYM_UpTo && n_arg->data.function.arg_count == 1) {
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
