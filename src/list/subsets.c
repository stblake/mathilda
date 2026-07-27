/* Subsets — enumerate the sublists of an expression.
 *
 * Mathematica semantics:
 *
 *   Subsets[list]                 the power set, ordered by increasing length
 *                                 and lexicographically by original element
 *                                 position within each length:
 *                                 Subsets[{a,b,c}] ->
 *                                   {{}, {a}, {b}, {c}, {a,b}, {a,c}, {b,c},
 *                                    {a,b,c}}
 *   Subsets[list, n]              lengths 0 through n inclusive
 *   Subsets[list, {n}]            exactly length n
 *   Subsets[list, {nmin, nmax}]   lengths nmin..nmax inclusive
 *   Subsets[list, {nmin, nmax, d} lengths nmin, nmin+d, ... up to nmax
 *   Subsets[list, spec, s]        only the first s subsets the spec produces
 *
 * The head of the inner sublists is taken from the input expression, so
 * Subsets[f[a,b]] gives {f[], f[a], f[b], f[a,b]}. The outer wrapper is always
 * a List. Duplicate elements are treated as distinct by position: no dedup is
 * performed, hence Subsets[{a,a}] -> {{}, {a}, {a}, {a,a}}.
 *
 * PERFORMANCE
 *
 * The full result is exponential in Length[list], so the generator is lazy: it
 * walks index combinations with an odometer and stops the instant the `s`
 * budget is exhausted. Subsets[<30 elements>, All, 5] therefore costs five
 * sublist allocations, not 2^30. Output storage grows geometrically and is
 * never sized from the theoretical subset count (which would overflow for
 * moderate lengths anyway). */

#include "list_common.h"
#include "subsets.h"

/* ------------------------------------------------------------------------- */
/* Growable output buffer                                                    */
/* ------------------------------------------------------------------------- */

typedef struct {
    Expr** items;
    size_t count;
    size_t cap;
} SubsetBuf;

/* Append `e` (ownership transferred on success). Returns false on allocation
 * failure, in which case `e` is untouched and the caller must free it. */
static bool subsets_buf_push(SubsetBuf* b, Expr* e) {
    if (b->count == b->cap) {
        size_t ncap = b->cap ? b->cap * 2 : 8;
        Expr** grown = realloc(b->items, ncap * sizeof(Expr*));
        if (!grown) return false;
        b->items = grown;
        b->cap = ncap;
    }
    b->items[b->count++] = e;
    return true;
}

static void subsets_buf_free(SubsetBuf* b) {
    for (size_t i = 0; i < b->count; i++) expr_free(b->items[i]);
    free(b->items);
    b->items = NULL;
    b->count = b->cap = 0;
}

/* ------------------------------------------------------------------------- */
/* Length-spec decoding                                                      */
/* ------------------------------------------------------------------------- */

typedef enum {
    LEN_OK,     /* a usable subset length was read into *out */
    LEN_NONE,   /* a well-formed but negative length: selects no subsets */
    LEN_BAD     /* not a length at all: leave the expression unevaluated */
} LenStatus;

/* Read one component of a length spec as a subset length over a list of `n`
 * elements. Accepts a non-negative machine integer, a bigint (a positive one
 * necessarily exceeds any reachable length, so it saturates at `n`), and
 * Infinity. */
static LenStatus subsets_read_length(Expr* e, size_t n, size_t* out) {
    if (e->type == EXPR_INTEGER) {
        if (e->data.integer < 0) return LEN_NONE;
        *out = (size_t)e->data.integer;
        return LEN_OK;
    }
    if (e->type == EXPR_BIGINT) {
        if (mpz_sgn(e->data.bigint) < 0) return LEN_NONE;
        *out = n;               /* saturate: no list is this long */
        return LEN_OK;
    }
    if (is_infinity(e)) {
        *out = n;
        return LEN_OK;
    }
    return LEN_BAD;
}

/* Decode the optional second argument into a half-open-free inclusive length
 * range [*kmin, *kmax] walked in steps of *kstep.
 *
 * Returns false when the spec is not one Subsets understands, so the caller
 * can leave the whole expression unevaluated. Returns true with *empty set
 * when the spec is well-formed but selects no lengths at all (a negative
 * bound, nmin > nmax, or an exact length that exceeds Length[list]). */
static bool subsets_decode_lengths(Expr* spec, size_t n,
                                   size_t* kmin, size_t* kmax, size_t* kstep,
                                   bool* empty) {
    *kmin = 0;
    *kmax = n;
    *kstep = 1;
    *empty = false;

    /* Subsets[list] and Subsets[list, All] — the whole power set. */
    if (spec == NULL) return true;
    if (spec->type == EXPR_SYMBOL && spec->data.symbol.name == SYM_All) return true;

    /* Subsets[list, n] — lengths 0..n, clamped to Length[list]. */
    if (spec->type != EXPR_FUNCTION) {
        size_t k;
        switch (subsets_read_length(spec, n, &k)) {
            case LEN_NONE: *empty = true; return true;
            case LEN_BAD:  return false;
            case LEN_OK:   break;
        }
        *kmax = (k < n) ? k : n;
        return true;
    }

    if (!is_listq(spec)) return false;

    size_t parts = spec->data.function.arg_count;
    if (parts < 1 || parts > 3) return false;

    /* Subsets[list, {n}] — exactly length n. Unlike the bare-integer form this
     * does not clamp: an n past the end simply has no subsets. */
    if (parts == 1) {
        size_t k;
        switch (subsets_read_length(spec->data.function.args[0], n, &k)) {
            case LEN_NONE: *empty = true; return true;
            case LEN_BAD:  return false;
            case LEN_OK:   break;
        }
        if (k > n) { *empty = true; return true; }
        *kmin = *kmax = k;
        return true;
    }

    /* Subsets[list, {nmin, nmax}] and Subsets[list, {nmin, nmax, d}]. */
    size_t lo, hi;
    switch (subsets_read_length(spec->data.function.args[0], n, &lo)) {
        case LEN_NONE: *empty = true; return true;
        case LEN_BAD:  return false;
        case LEN_OK:   break;
    }
    switch (subsets_read_length(spec->data.function.args[1], n, &hi)) {
        case LEN_NONE: *empty = true; return true;
        case LEN_BAD:  return false;
        case LEN_OK:   break;
    }
    if (hi > n) hi = n;         /* an over-long upper bound clamps */
    if (parts == 3) {
        Expr* d = spec->data.function.args[2];
        if (d->type != EXPR_INTEGER) return false;
        if (d->data.integer <= 0) { *empty = true; return true; }
        *kstep = (size_t)d->data.integer;
    }
    if (lo > hi) { *empty = true; return true; }
    *kmin = lo;
    *kmax = hi;
    return true;
}

/* Decode the optional third argument, the number of subsets to keep. Accepts a
 * non-negative machine integer, All / Infinity / a positive bigint (all "no
 * limit"). Returns false for anything else, including a negative count. */
static bool subsets_decode_count(Expr* spec, size_t* limit) {
    if (spec->type == EXPR_INTEGER) {
        if (spec->data.integer < 0) return false;
        *limit = (size_t)spec->data.integer;
        return true;
    }
    if (spec->type == EXPR_BIGINT) {
        if (mpz_sgn(spec->data.bigint) < 0) return false;
        *limit = SIZE_MAX;
        return true;
    }
    if (is_infinity(spec)) { *limit = SIZE_MAX; return true; }
    if (spec->type == EXPR_SYMBOL && spec->data.symbol.name == SYM_All) {
        *limit = SIZE_MAX;
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------------- */
/* Combination odometer                                                      */
/* ------------------------------------------------------------------------- */

/* Advance `idx` (a strictly increasing k-tuple of indices drawn from
 * 0..n-1) to the next combination in lexicographic order. Returns false once
 * the final combination has been passed. Requires k <= n. */
static bool subsets_comb_next(size_t* idx, size_t k, size_t n) {
    if (k == 0) return false;   /* the empty combination is unique */
    size_t i = k;
    while (i > 0) {
        i--;
        /* idx[i] tops out at n - k + i, leaving room for the tail. */
        if (idx[i] < n - k + i) {
            idx[i]++;
            for (size_t j = i + 1; j < k; j++) idx[j] = idx[j - 1] + 1;
            return true;
        }
    }
    return false;
}

/* ------------------------------------------------------------------------- */
/* Builtin                                                                   */
/* ------------------------------------------------------------------------- */

Expr* builtin_subsets(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return NULL;

    Expr* list = res->data.function.args[0];
    if (list->type != EXPR_FUNCTION) return NULL;   /* atoms have no sublists */

    Expr* head = list->data.function.head;
    Expr** elems = list->data.function.args;
    size_t n = list->data.function.arg_count;

    size_t kmin, kmax, kstep;
    bool empty = false;
    if (!subsets_decode_lengths(argc >= 2 ? res->data.function.args[1] : NULL,
                                n, &kmin, &kmax, &kstep, &empty)) {
        return NULL;
    }

    size_t limit = SIZE_MAX;
    if (argc >= 3 && !subsets_decode_count(res->data.function.args[2], &limit)) return NULL;

    if (empty || limit == 0) return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);

    /* One scratch buffer sized for the largest length we will emit, reused
     * across every k so the generator allocates no per-combination state. */
    size_t* idx = NULL;
    if (kmax > 0) {
        idx = malloc(kmax * sizeof(size_t));
        if (!idx) return NULL;
    }

    SubsetBuf out = {NULL, 0, 0};
    bool ok = true;

    for (size_t k = kmin; k <= kmax && ok && out.count < limit; k += kstep) {
        /* Defensive: every spec path above clamps kmax to n, so this cannot
         * fire today. It stays because subsets_comb_next computes n - k + i in
         * size_t, which would wrap catastrophically rather than fail loudly if
         * a future spec form ever let k exceed n. */
        if (k > n) break;

        for (size_t i = 0; i < k; i++) idx[i] = i;
        do {
            Expr** sub = NULL;
            if (k > 0) {
                sub = malloc(k * sizeof(Expr*));
                if (!sub) { ok = false; break; }
                for (size_t i = 0; i < k; i++) sub[i] = expr_copy(elems[idx[i]]);
            }
            Expr* subset = expr_new_function(expr_copy(head), sub, k);
            free(sub);
            if (!subset || !subsets_buf_push(&out, subset)) {
                expr_free(subset);
                ok = false;
                break;
            }
            if (out.count == limit) break;
        } while (subsets_comb_next(idx, k, n));

        /* Stop before `k += kstep` can wrap past SIZE_MAX. The loop guard
         * (k <= kmax) would not catch a wrapped k. */
        if (kmax - k < kstep) break;
    }

    free(idx);
    if (!ok) {
        subsets_buf_free(&out);
        return NULL;
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), out.items, out.count);
    free(out.items);
    return result;
}
