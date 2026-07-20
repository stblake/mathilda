#include "expr.h"
#include "arithmetic.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "ndarray.h"   /* ndt_elem_size for dtype-aware copy/eq/hash */
#include <stdbool.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

/* Portable strdup (strdup is POSIX, not C99). See expr.h. */
char* mathilda_strdup(const char* s) {
    if (!s) return NULL;
    size_t n = strlen(s) + 1;
    char* r = (char*)malloc(n);
    if (r) memcpy(r, s, n);
    return r;
}

/* ------------------------------------------------------------------------
 *  Expr node pool (bounded free-list)
 *
 *  Profiling the logistic-map benchmark (Do[x = 3.5 x (1-x), {10^6}]) showed
 *  ~65% of wall time in the system allocator: the evaluator churns fixed-size
 *  Expr structs, allocating and freeing tens of nodes per iteration. Because
 *  every Expr is exactly sizeof(Expr) bytes, freed nodes are recycled through a
 *  singly-linked free-list instead of round-tripping through malloc/free. The
 *  link is stored in the dead node's payload (data.function.head, a genuine
 *  Expr* at union offset 0 — no aliasing games). The REPL is single-threaded,
 *  so no locking is needed.
 *
 *  The free-list is BOUNDED (EXPR_POOL_CAP nodes). This is the crucial detail:
 *  an unbounded pool accumulates every freed node and, after mixed churn,
 *  threads across memory in scrambled order — so a large batch allocation (e.g.
 *  KeySort building a 40k-key association, then sorting it) pulls scattered
 *  nodes and the sort thrashes the cache (its O(n log n) doubling ratio
 *  regressed from ~2.2 to ~3.8). Capping keeps the recycled set small and
 *  cache-hot: the evaluator's working set (tens of nodes) always hits it, while
 *  a large batch drains the cap and then falls to fresh malloc, whose natural
 *  burst locality restores the sort's scaling. Frees past the cap go straight
 *  back to the system allocator.
 *
 *  expr_pool_free_all() drains the (bounded) free-list back to the OS; it is
 *  registered with atexit() on first use so every binary leaves a clean heap
 *  for valgrind. Nodes still in use at exit are freed by their owners exactly
 *  as before — the pool only holds nodes that have already been expr_free'd.
 * ---------------------------------------------------------------------- */
#define EXPR_POOL_CAP 8192                /* max recycled nodes held (~0.5 MB) */
static Expr*  g_expr_pool = NULL;         /* free-list head */
static size_t g_expr_pool_size = 0;       /* current free-list length */
static bool   g_expr_pool_atexit_registered = false;

/* Drain the free-list back to the system allocator. Idempotent; registered
 * with atexit(). Only touches already-freed (pooled) nodes. */
void expr_pool_free_all(void) {
    Expr* e = g_expr_pool;
    while (e) {
        Expr* next = e->data.function.head;
        free(e);
        e = next;
    }
    g_expr_pool = NULL;
    g_expr_pool_size = 0;
}

static Expr* expr_alloc_node(void) {
    Expr* e = g_expr_pool;
    if (e) {
        g_expr_pool = e->data.function.head;   /* pop */
        g_expr_pool_size--;
        return e;
    }
    if (!g_expr_pool_atexit_registered) {
        g_expr_pool_atexit_registered = true;
        atexit(expr_pool_free_all);
    }
    return (Expr*)malloc(sizeof(Expr));
}

/* Return a physically-dead node to the pool, or to the OS once the pool is at
 * capacity. Caller must have already released any owned payload. */
static void expr_release_node(Expr* e) {
    if (g_expr_pool_size >= EXPR_POOL_CAP) {
        free(e);
        return;
    }
    e->data.function.head = g_expr_pool;       /* push */
    g_expr_pool = e;
    g_expr_pool_size++;
}

// Create/allocate a new integer expression.
Expr* expr_new_integer(int64_t value) {
    Expr* e = expr_alloc_node();
    if (!e) return NULL;
    e->type = EXPR_INTEGER;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    e->data.integer = value;
    return e;
}

// Create/allocate a new real (double) expression.
Expr* expr_new_real(double value) {
    Expr* e = expr_alloc_node();
    if (!e) return NULL;
    e->type = EXPR_REAL;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    e->data.real = value;
    return e;
}

// Create/allocate a new symbol expression.
//
// The symbol name is funneled through the global interner so that two
// symbols with the same name share the same `const char*`. This makes:
//   - expr_eq() on symbols a pointer compare,
//   - expr_copy() of a symbol a pointer copy,
//   - expr_free() of a symbol a no-op for the name field.
//
// The cast to `char*` is intentional: `data.symbol.name` retains its existing
// type for ABI stability, but the memory is owned by the interner and
// must never be freed by Expr code.
Expr* expr_new_symbol(const char* name) {
    Expr* e = expr_alloc_node();
    if (!e) return NULL;
    e->type = EXPR_SYMBOL;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    e->data.symbol.name = (char*)intern_symbol(name);
    e->data.symbol.def = NULL;   /* Phase 3b: resolved + cached lazily on first eval use */
    if (!e->data.symbol.name) {
        free(e);
        return NULL;
    }
    return e;
}

// Create/allocate a new string expression.
Expr* expr_new_string(const char* str) {
    Expr* e = expr_alloc_node();
    if (!e) return NULL;

    e->type = EXPR_STRING;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    e->data.string = mathilda_strdup(str);
    if (!e->data.string) {
        free(e);
        return NULL;
    }
    return e;
}

// Create/allocate an expression: h[arg1, arg2, ...]
Expr* expr_new_function(Expr* head, Expr** args, size_t arg_count) {
    Expr* e = expr_alloc_node();
    if (!e) return NULL;

    e->type = EXPR_FUNCTION;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    e->data.function.head = head;
    if (arg_count > 0) {
        e->data.function.args = calloc(arg_count, sizeof(Expr*));
        if (!e->data.function.args) {
            free(e);
            return NULL;
        }
        if (args) memcpy(e->data.function.args, args, sizeof(Expr*) * arg_count);
    } else {
        e->data.function.args = NULL;
    }
    e->data.function.arg_count = arg_count;
    return e;
}

// BigInt constructors
Expr* expr_new_bigint_from_mpz(const mpz_t val) {
    Expr* e = expr_alloc_node();
    if (!e) return NULL;
    e->type = EXPR_BIGINT;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    mpz_init_set(e->data.bigint, val);
    return e;
}

Expr* expr_new_bigint_from_int64(int64_t val) {
    Expr* e = expr_alloc_node();
    if (!e) return NULL;
    e->type = EXPR_BIGINT;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    mpz_init_set_si(e->data.bigint, val);
    return e;
}

Expr* expr_new_bigint_from_str(const char* str) {
    Expr* e = expr_alloc_node();
    if (!e) return NULL;
    e->type = EXPR_BIGINT;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    if (mpz_init_set_str(e->data.bigint, str, 10) == -1) {
        mpz_clear(e->data.bigint);
        free(e);
        return NULL;
    }
    return e;
}

Expr* expr_new_ndarray(int rank, const int64_t* dims, void* data, NDType dtype) {
    Expr* e = expr_alloc_node();
    if (!e) return NULL;
    e->type = EXPR_NDARRAY;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    e->data.ndarray.rank = rank;
    e->data.ndarray.dims = malloc(sizeof(int64_t) * (size_t)rank);
    if (!e->data.ndarray.dims) { free(e); return NULL; }
    memcpy(e->data.ndarray.dims, dims, sizeof(int64_t) * (size_t)rank);
    e->data.ndarray.data = data;
    e->data.ndarray.dtype = dtype;
    return e;
}

#ifdef USE_MPFR
/* MPFR constructors. All allocate an Expr and initialize the payload
 * `mpfr_t` at the requested precision; the caller owns the result and
 * should free it with `expr_free`, which calls `mpfr_clear`. */
Expr* expr_new_mpfr_bits(mpfr_prec_t bits) {
    Expr* e = expr_alloc_node();
    if (!e) return NULL;
    e->type = EXPR_MPFR;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    mpfr_init2(e->data.mpfr, bits);
    mpfr_set_zero(e->data.mpfr, +1);
    return e;
}
Expr* expr_new_mpfr_from_d(double v, mpfr_prec_t bits) {
    Expr* e = expr_new_mpfr_bits(bits);
    if (e) mpfr_set_d(e->data.mpfr, v, MPFR_RNDN);
    return e;
}
Expr* expr_new_mpfr_from_si(long v, mpfr_prec_t bits) {
    Expr* e = expr_new_mpfr_bits(bits);
    if (e) mpfr_set_si(e->data.mpfr, v, MPFR_RNDN);
    return e;
}
Expr* expr_new_mpfr_from_mpz(const mpz_t z, mpfr_prec_t bits) {
    Expr* e = expr_new_mpfr_bits(bits);
    if (e) mpfr_set_z(e->data.mpfr, z, MPFR_RNDN);
    return e;
}
Expr* expr_new_mpfr_from_str(const char* str, mpfr_prec_t bits) {
    Expr* e = expr_new_mpfr_bits(bits);
    if (!e) return NULL;
    if (mpfr_set_str(e->data.mpfr, str, 10, MPFR_RNDN) != 0) {
        mpfr_clear(e->data.mpfr);
        free(e);
        return NULL;
    }
    return e;
}
Expr* expr_new_mpfr_move(mpfr_t src) {
    Expr* e = expr_alloc_node();
    if (!e) { mpfr_clear(src); return NULL; }
    e->type = EXPR_MPFR;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    /* mpfr_t is an array-type alias for __mpfr_struct[1]; `memcpy` moves
     * the header — MPFR's documentation (mpfr_swap, GMP's mpz semantics)
     * sanctions this as long as the source is then treated as uninit. */
    memcpy(e->data.mpfr, src, sizeof(e->data.mpfr));
    return e;
}
Expr* expr_new_mpfr_copy(const mpfr_t src) {
    Expr* e = expr_new_mpfr_bits(mpfr_get_prec(src));
    if (e) mpfr_set(e->data.mpfr, src, MPFR_RNDN);
    return e;
}
#endif

void expr_to_mpz(const Expr* e, mpz_t out) {
    if (e->type == EXPR_INTEGER) {
        mpz_init_set_si(out, e->data.integer);
    } else { // EXPR_BIGINT
        mpz_init_set(out, e->data.bigint);
    }
}

bool expr_is_integer_like(const Expr* e) {
    return e && (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT);
}

bool expr_is_numeric_like(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:
        case EXPR_REAL:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        case EXPR_FUNCTION:
            /* Rational[n,d] and Complex[re,im] with numeric components.
             * Use the bigint-aware rational check so rationals whose
             * components have overflowed int64 still fold via the GMP
             * fallback in multiply_numbers / add_numbers (otherwise
             * Times[BigInt, Rational[1, BigInt]] survives unsimplified
             * and breaks exact polynomial division). */
            if (is_rational_like((Expr*)e)) return true;
            {
                Expr *re, *im;
                if (is_complex((Expr*)e, &re, &im)) {
                    return expr_is_numeric_like(re) && expr_is_numeric_like(im);
                }
            }
            return false;
        default:
            return false;
    }
}

Expr* expr_bigint_normalize(Expr* e) {
    if (e->type == EXPR_BIGINT && mpz_fits_slong_p(e->data.bigint)) {
        long val = mpz_get_si(e->data.bigint);
        Expr* result = expr_new_integer((int64_t)val);
        /* Refcount-safe: dec-ref the BigInt; if shared, the original
         * stays alive for other holders (each will normalize on their
         * own); if unique (refcount==1, the common case), expr_free
         * physically clears the mpz_t and frees the node. */
        expr_free(e);
        return result;
    }
    return e;
}

/* Inc-ref. NULL passes through. */
Expr* expr_ref(Expr* e) {
    if (e) e->refcount++;
    return e;
}

/* M3 phase-2 copy-on-write helper. See header doc. */
Expr* expr_unshare(Expr* e) {
    if (!e) return NULL;
    if (e->refcount == 1) return e;

    /* Build a one-level private copy. Children stay shared (inc-ref). */
    Expr* fresh = expr_alloc_node();
    if (!fresh) return e;  /* OOM: fall back to (still-shared) original */
    fresh->type = e->type;
    fresh->refcount = 1;
    /* Reset the eval timestamp on unshare. The caller is preparing to
     * mutate this node; even if they do not, treating the node as
     * "never evaluated" only costs one extra evaluation, whereas
     * inheriting the timestamp risks a false cache hit on a node that
     * the caller subsequently rewrites. */
    fresh->last_evaluated_at = 0;

    switch (e->type) {
        case EXPR_INTEGER:
            fresh->data.integer = e->data.integer;
            break;
        case EXPR_REAL:
            fresh->data.real = e->data.real;
            break;
        case EXPR_SYMBOL:
            fresh->data.symbol.name = e->data.symbol.name;  /* interned */
            fresh->data.symbol.def = e->data.symbol.def;    /* same symbol: cache carries over */
            break;
        case EXPR_STRING:
            fresh->data.string = e->data.string ? mathilda_strdup(e->data.string) : NULL;
            break;
        case EXPR_BIGINT:
            mpz_init(fresh->data.bigint);
            mpz_set(fresh->data.bigint, e->data.bigint);
            break;
        case EXPR_NDARRAY: {
            int rank = e->data.ndarray.rank;
            NDType dt = e->data.ndarray.dtype;
            size_t n = 1;
            for (int i = 0; i < rank; i++) n *= (size_t)e->data.ndarray.dims[i];
            size_t bytes = ndt_elem_size(dt) * n;
            fresh->data.ndarray.rank = rank;
            fresh->data.ndarray.dtype = dt;
            fresh->data.ndarray.dims = malloc(sizeof(int64_t) * (size_t)rank);
            memcpy(fresh->data.ndarray.dims, e->data.ndarray.dims, sizeof(int64_t) * (size_t)rank);
            fresh->data.ndarray.data = malloc(bytes);
            memcpy(fresh->data.ndarray.data, e->data.ndarray.data, bytes);
            break;
        }
#ifdef USE_MPFR
        case EXPR_MPFR:
            mpfr_init2(fresh->data.mpfr, mpfr_get_prec(e->data.mpfr));
            mpfr_set(fresh->data.mpfr, e->data.mpfr, MPFR_RNDN);
            break;
#endif
        case EXPR_FUNCTION:
            fresh->data.function.head = expr_copy(e->data.function.head);
            fresh->data.function.arg_count = e->data.function.arg_count;
            if (e->data.function.arg_count > 0) {
                fresh->data.function.args = malloc(sizeof(Expr*) * e->data.function.arg_count);
                if (!fresh->data.function.args) {
                    if (fresh->data.function.head) expr_free(fresh->data.function.head);
                    free(fresh);
                    return e;
                }
                for (size_t i = 0; i < e->data.function.arg_count; i++) {
                    fresh->data.function.args[i] = expr_copy(e->data.function.args[i]);
                }
            } else {
                fresh->data.function.args = NULL;
            }
            break;
        default:
            break;
    }

    /* Drop the caller's ref on the original. */
    expr_free(e);
    return fresh;
}

// Deallocate an expression. Decrements the refcount; only physically
// destroys the node (and its children, for FUNCTION) when refcount
// drops to 0. Atoms are eligible to have refcount > 1 starting in M3
// (see expr_copy). FUNCTION nodes always have refcount == 1 today
// because expr_copy still deep-copies them.
void expr_free(Expr* e) {
    if (!e) return;
    if (e->refcount > 1) {
        e->refcount--;
        return;
    }
    /* refcount == 1 (or, defensively, 0 — should not happen, but treat
     * as last reference). Physically free. */

    switch (e->type) {
        case EXPR_SYMBOL:
            /* data.symbol.name is interned (owned by sym_intern); never free. */
            break;
        case EXPR_STRING:
            if (e->data.string) free(e->data.string);
            break;
        case EXPR_FUNCTION:
            if (e->data.function.head) expr_free(e->data.function.head);
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (e->data.function.args && e->data.function.args[i]) {
                    expr_free(e->data.function.args[i]);
                }
            }
            if (e->data.function.args) free(e->data.function.args);
            break;
        case EXPR_BIGINT:
            mpz_clear(e->data.bigint);
            break;
        case EXPR_NDARRAY:
            free(e->data.ndarray.dims);
            free(e->data.ndarray.data);
            break;
#ifdef USE_MPFR
        case EXPR_MPFR:
            mpfr_clear(e->data.mpfr);
            break;
#endif
        default:
            break;
    }
    /* Recycle the fixed-size node through the pool instead of free() —
     * see the free-list note above expr_new_integer. */
    expr_release_node(e);
}



// Logically copy an expression. M3 phase-2: every node type, INCLUDING
// FUNCTION, is shared via inc-ref. Mutating helpers must call
// expr_unshare() to obtain a private (refcount==1) version before
// rewriting fields in place. The audit (commit notes) verified that
// the only sites that mutated a possibly-shared FUNCTION node lived
// in print.c; those now wrap their expr_copy() in expr_unshare().
// Other mutators (eval_flatten_args, flatten_sequences, builtin_plus
// numeric-contagion args[i] writes, parse.c arg appending, match.c
// sequence-binding writes, core.c QuotientRemainder zeroing, etc.)
// all act on freshly-constructed FUNCTION nodes whose refcount is
// guaranteed to be 1 at the point of mutation.
Expr* expr_copy(Expr* e) {
    if (!e) return NULL;
    e->refcount++;
    return e;
}

bool expr_eq(const Expr* a, const Expr* b) {
    if (a == b) return true;
    if (!a || !b) return false;
    if (a->type != b->type) {
        // Cross-type equality: EXPR_INTEGER vs EXPR_BIGINT
        if ((a->type == EXPR_INTEGER && b->type == EXPR_BIGINT) ||
            (a->type == EXPR_BIGINT && b->type == EXPR_INTEGER)) {
            mpz_t va, vb;
            expr_to_mpz(a, va);
            expr_to_mpz(b, vb);
            bool result = (mpz_cmp(va, vb) == 0);
            mpz_clear(va);
            mpz_clear(vb);
            return result;
        }
        return false;
    }

    switch (a->type) {
        case EXPR_INTEGER:
            return a->data.integer == b->data.integer;
        case EXPR_REAL:
            if (isnan(a->data.real) && isnan(b->data.real)) return true;
            return a->data.real == b->data.real;
        case EXPR_SYMBOL:
            /* Interned: same name guarantees same pointer. */
            return a->data.symbol.name == b->data.symbol.name;
        case EXPR_STRING:
            return strcmp(a->data.string, b->data.string) == 0;
        case EXPR_FUNCTION: {
            if (!expr_eq(a->data.function.head, b->data.function.head)) return false;
            if (a->data.function.arg_count != b->data.function.arg_count) return false;
            for (size_t i = 0; i < a->data.function.arg_count; i++) {
                if (!expr_eq(a->data.function.args[i], b->data.function.args[i])) return false;
            }
            return true;
        }
        case EXPR_BIGINT:
            return mpz_cmp(a->data.bigint, b->data.bigint) == 0;
        case EXPR_NDARRAY: {
            /* dtype is part of identity (a float32 array is not SameQ to a
             * float64 one with equal values, matching the MPFR-precision rule
             * below). */
            if (a->data.ndarray.dtype != b->data.ndarray.dtype) return false;
            if (a->data.ndarray.rank != b->data.ndarray.rank) return false;
            size_t n = 1;
            for (int i = 0; i < a->data.ndarray.rank; i++) {
                if (a->data.ndarray.dims[i] != b->data.ndarray.dims[i]) return false;
                n *= (size_t)a->data.ndarray.dims[i];
            }
            size_t bytes = ndt_elem_size(a->data.ndarray.dtype) * n;
            return memcmp(a->data.ndarray.data, b->data.ndarray.data, bytes) == 0;
        }
#ifdef USE_MPFR
        case EXPR_MPFR:
            /* Equal iff same precision AND same value (matches SameQ
             * semantics: 1.`20 and 1.`30 are not SameQ even though
             * their values agree). */
            if (mpfr_get_prec(a->data.mpfr) != mpfr_get_prec(b->data.mpfr)) return false;
            return mpfr_equal_p(a->data.mpfr, b->data.mpfr) != 0;
#endif
    }
    return false;
}

/* expr_compare lives in sort.c so the canonical-order policy stays in
 * one place alongside the user-facing Sort/OrderedQ builtins.  The
 * helpers it needs (is_atomic_numeric, expr_poly_degree, etc.) are
 * static there. */

uint64_t expr_hash(const Expr* e) {
    if (!e) return 0;
    uint64_t h = 14695981039346656037ULL;
    const uint64_t prime = 1099511628211ULL;

    h ^= (uint64_t)e->type;
    h *= prime;

    switch (e->type) {
        case EXPR_INTEGER: {
            uint64_t v;
            memcpy(&v, &e->data.integer, 8);
            h ^= v; h *= prime;
            break;
        }
        case EXPR_REAL: {
            uint64_t v;
            memcpy(&v, &e->data.real, 8);
            h ^= v; h *= prime;
            break;
        }
        case EXPR_SYMBOL:
        case EXPR_STRING: {
            const char* s = e->data.symbol.name;
            if (s) {
                while (*s) {
                    h ^= (uint8_t)(*s++);
                    h *= prime;
                }
            }
            break;
        }
        case EXPR_FUNCTION: {
            h ^= expr_hash(e->data.function.head);
            h *= prime;
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                h ^= expr_hash(e->data.function.args[i]);
                h *= prime;
            }
            break;
        }
        case EXPR_BIGINT: {
            h ^= (uint64_t)(mpz_sgn(e->data.bigint) + 2);
            h *= prime;
            size_t nlimbs = mpz_size(e->data.bigint);
            for (size_t i = 0; i < nlimbs; i++) {
                h ^= mpz_getlimbn(e->data.bigint, i);
                h *= prime;
            }
            break;
        }
        case EXPR_NDARRAY: {
            h ^= (uint64_t)e->data.ndarray.rank;
            h *= prime;
            h ^= (uint64_t)e->data.ndarray.dtype;
            h *= prime;
            size_t n = 1;
            for (int i = 0; i < e->data.ndarray.rank; i++) {
                h ^= (uint64_t)e->data.ndarray.dims[i];
                h *= prime;
                n *= (size_t)e->data.ndarray.dims[i];
            }
            /* Hash the raw buffer as bytes so all dtypes (incl. 4-byte float /
             * float32-complex) hash correctly. */
            const unsigned char* bytes = (const unsigned char*)e->data.ndarray.data;
            size_t nbytes = ndt_elem_size(e->data.ndarray.dtype) * n;
            for (size_t i = 0; i < nbytes; i++) {
                h ^= (uint64_t)bytes[i]; h *= prime;
            }
            break;
        }
#ifdef USE_MPFR
        case EXPR_MPFR: {
            /* Hash precision + IEEE double approximation. Not perfectly
             * collision-free across precisions with equal double value,
             * but it's a hash, not an equality. */
            h ^= (uint64_t)mpfr_get_prec(e->data.mpfr);
            h *= prime;
            double d = mpfr_get_d(e->data.mpfr, MPFR_RNDN);
            uint64_t v; memcpy(&v, &d, 8);
            h ^= v; h *= prime;
            break;
        }
#endif
    }
    return h;
}

