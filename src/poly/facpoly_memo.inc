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
 * common case.  In practice only the top entry is consulted.
 *
 * A TimeConstrained timeout that siglongjmps out of a Simplify skips its
 * factor_memo_pop.  Before g_memo_stack_top was a registered unwind counter
 * (facpoly_init), that Simplify's memo stayed on top for the rest of the
 * session: every later Factor / Together / Cancel, in any computation and
 * under any $Assumptions, read and grew it.  It also held entries made while
 * the timeout was in flight, when the evaluator answers every call with the
 * Abort[] sentinel -- a Together whose inner evaluation aborted had stored a
 * value with Abort[] coefficients, a later PolynomialExtendedGCD read it back,
 * evaluated the Abort[], and aborted a computation no constraint owned
 * (EVAL-abort-1).  The jump now restores the depth.  A non-NULL slot above it
 * is the memo of a frame the jump abandoned (factor_memo_pop clears the slot
 * of a memo it hands back to its owner): the next push into the slot frees it. */
#define FACTOR_MEMO_STACK_DEPTH 8
static FactorMemo* g_memo_stack[FACTOR_MEMO_STACK_DEPTH] = {0};
static int g_memo_stack_top = -1;
int* factor_memo_depth_counter(void) { return &g_memo_stack_top; }

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
    ++g_memo_stack_top;
    FactorMemo* abandoned = g_memo_stack[g_memo_stack_top];   /* see above */
    if (abandoned && abandoned != m) factor_memo_free(abandoned);
    g_memo_stack[g_memo_stack_top] = m;
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
    /* A value computed while an abort is in flight may carry the Abort[]
     * sentinel inside it, where C code that ignored it built it in. */
    if (eval_abort_raised()) return;
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
    factor_memo_store(m, key, value);
}
