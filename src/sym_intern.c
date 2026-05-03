/*
 * sym_intern.c
 *
 * Global string interner. A simple chained hash table keyed on the
 * symbol name. Each unique name is stored exactly once; intern_symbol()
 * returns the table-owned canonical pointer.
 *
 * Used by expr_new_symbol() so that every EXPR_SYMBOL Expr's data.symbol
 * is the same pointer for the same name. This reduces:
 *   - SYMBOL equality from strcmp() to pointer compare
 *   - expr_copy() of a symbol from strdup() to a pointer copy
 *   - expr_free() of a symbol from free() to nothing
 *
 * Also used by symtab.c so that SymbolDef::symbol_name shares the same
 * canonical pointer, enabling pointer-keyed lookups.
 */

#include "sym_intern.h"
#include <stdlib.h>
#include <string.h>

#define INTERN_BUCKETS 8191  /* prime; ~3.5x typical symbol-set size */

typedef struct InternEntry {
    char* name;                  /* malloc'd, owned by this table */
    struct InternEntry* next;
} InternEntry;

static InternEntry* g_buckets[INTERN_BUCKETS];
/* Static storage zero-initializes g_buckets to all-NULL on first use,
 * so no separate init step is required. */

/* djb2 hash, identical to the symtab hash so name -> bucket is uniform. */
static unsigned int hash_string(const char* s) {
    unsigned int h = 5381;
    while (*s) {
        h = ((h << 5) + h) + (unsigned char)*s++;
    }
    return h;
}

const char* intern_symbol(const char* name) {
    if (!name) return NULL;
    unsigned int idx = hash_string(name) % INTERN_BUCKETS;
    InternEntry* e = g_buckets[idx];
    while (e) {
        /* Pointer-compare first: a name already returned by this table
         * will hit instantly without walking the bytes. Only fall back
         * to strcmp for fresh string-literal lookups. */
        if (e->name == name || strcmp(e->name, name) == 0) {
            return e->name;
        }
        e = e->next;
    }
    /* Not present: allocate a stable copy and link it in. */
    size_t len = strlen(name);
    char* copy = malloc(len + 1);
    if (!copy) return NULL;
    memcpy(copy, name, len + 1);

    InternEntry* ne = malloc(sizeof(InternEntry));
    if (!ne) {
        free(copy);
        return NULL;
    }
    ne->name = copy;
    ne->next = g_buckets[idx];
    g_buckets[idx] = ne;
    return copy;
}

void intern_clear(void) {
    for (int i = 0; i < INTERN_BUCKETS; i++) {
        InternEntry* e = g_buckets[i];
        while (e) {
            InternEntry* n = e->next;
            free(e->name);
            free(e);
            e = n;
        }
        g_buckets[i] = NULL;
    }
}
