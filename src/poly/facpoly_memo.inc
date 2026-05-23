/* ===================================================================== */
/*  Per-call Factor memo                                                 */
/*                                                                       */
/*  Cache `Factor[poly]` results inside a single `Simplify` call.  See   */
/*  facpoly.h for the lifecycle contract.  The top-of-stack memo is      */
/*  consulted at entry to `builtin_factor`; on a hit we return a deep    */
/*  copy of the cached result.  On a miss we run the full pipeline,     */
/*  then store a deep copy of the input/output pair before returning.   */
/* ===================================================================== */

#define FACTOR_MEMO_BUCKETS 256

typedef struct FactorMemoEntry {
    Expr* key;
    Expr* value;
    struct FactorMemoEntry* next;
} FactorMemoEntry;

struct FactorMemo {
    FactorMemoEntry* buckets[FACTOR_MEMO_BUCKETS];
};

/* Stack of currently-active memos.  We use a stack to support
 * (in principle) nested Simplify calls, even though that's not the
 * common case.  In practice only the top entry is consulted. */
#define FACTOR_MEMO_STACK_DEPTH 8
static FactorMemo* g_memo_stack[FACTOR_MEMO_STACK_DEPTH] = {0};
static int g_memo_stack_top = -1;

FactorMemo* factor_memo_new(void) {
    FactorMemo* m = (FactorMemo*)calloc(1, sizeof(FactorMemo));
    return m;
}

void factor_memo_free(FactorMemo* m) {
    if (!m) return;
    for (int i = 0; i < FACTOR_MEMO_BUCKETS; i++) {
        FactorMemoEntry* e = m->buckets[i];
        while (e) {
            FactorMemoEntry* next = e->next;
            expr_free(e->key);
            expr_free(e->value);
            free(e);
            e = next;
        }
    }
    free(m);
}

void factor_memo_push(FactorMemo* m) {
    if (g_memo_stack_top + 1 >= FACTOR_MEMO_STACK_DEPTH) return;
    g_memo_stack[++g_memo_stack_top] = m;
}

void factor_memo_pop(void) {
    if (g_memo_stack_top < 0) return;
    g_memo_stack[g_memo_stack_top--] = NULL;
}

static FactorMemo* factor_memo_top(void) {
    if (g_memo_stack_top < 0) return NULL;
    return g_memo_stack[g_memo_stack_top];
}

FactorMemo* factor_memo_active(void) { return factor_memo_top(); }

const Expr* factor_memo_lookup(FactorMemo* m, Expr* key) {
    if (!m) return NULL;
    uint64_t h = expr_hash(key) % FACTOR_MEMO_BUCKETS;
    for (FactorMemoEntry* e = m->buckets[h]; e; e = e->next) {
        if (expr_eq(e->key, key)) return e->value;
    }
    return NULL;
}

void factor_memo_store(FactorMemo* m, Expr* key, Expr* value) {
    if (!m) return;
    uint64_t h = expr_hash(key) % FACTOR_MEMO_BUCKETS;
    FactorMemoEntry* e = (FactorMemoEntry*)malloc(sizeof(FactorMemoEntry));
    if (!e) return;
    e->key = expr_copy(key);
    e->value = expr_copy(value);
    e->next = m->buckets[h];
    m->buckets[h] = e;
}

/* Lookup; returns a borrowed pointer to the cached value (do NOT
 * free) or NULL on miss. */
static const Expr* factor_memo_get(FactorMemo* m, Expr* key) {
    if (!m) return NULL;
    uint64_t h = expr_hash(key) % FACTOR_MEMO_BUCKETS;
    for (FactorMemoEntry* e = m->buckets[h]; e; e = e->next) {
        if (expr_eq(e->key, key)) return e->value;
    }
    return NULL;
}

/* Store deep copies of key and value. */
static void factor_memo_put(FactorMemo* m, Expr* key, Expr* value) {
    if (!m) return;
    uint64_t h = expr_hash(key) % FACTOR_MEMO_BUCKETS;
    FactorMemoEntry* e = (FactorMemoEntry*)malloc(sizeof(FactorMemoEntry));
    if (!e) return;
    e->key = expr_copy(key);
    e->value = expr_copy(value);
    e->next = m->buckets[h];
    m->buckets[h] = e;
}
