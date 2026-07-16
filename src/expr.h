
// Mathilda

#ifndef EXPR_H
#define EXPR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdarg.h>
#include <math.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

typedef enum {
    EXPR_INTEGER,
    EXPR_REAL,
    EXPR_SYMBOL,
    EXPR_STRING,
    EXPR_FUNCTION,
    EXPR_BIGINT,
    EXPR_NDARRAY               /* dense machine-precision ndarray, see ndarray.h */
#ifdef USE_MPFR
    , EXPR_MPFR                /* arbitrary-precision real (MPFR) */
#endif
} ExprType;

/* Element data type of an EXPR_NDARRAY, analogous to numpy's dtype. The four
 * types map one-to-one onto BLAS's four precisions (s/d/c/z), where the numeric
 * suffix is the *per-component* bit width:
 *   NDT_FLOAT64   double            (BLAS d) — DEFAULT, value 0 so a
 *                                   zero-initialized/legacy struct behaves as
 *                                   today's double-only NDArray.
 *   NDT_FLOAT32   float             (BLAS s)
 *   NDT_COMPLEX64 2x double (re,im) (BLAS z)
 *   NDT_COMPLEX32 2x float  (re,im) (BLAS c)
 * Complex is stored as interleaved (re, im) plain floats/doubles — no
 * <complex.h>, matching the existing BLAS interleaved-buffer convention and
 * keeping the code strict-C99 portable. See ndarray.h for the ndt_* helpers. */
typedef enum {
    NDT_FLOAT64 = 0,
    NDT_FLOAT32,
    NDT_COMPLEX64,
    NDT_COMPLEX32
} NDType;

/* Dense machine-precision ndarray payload for EXPR_NDARRAY. Row-major flat
 * storage: `data` holds dims[0]*dims[1]*...*dims[rank-1] elements of type
 * `dtype` (each element is ndt_elem_size(dtype) bytes). Always privately owned
 * (dims/data are never aliased between Expr nodes) so C fast paths can
 * read/write the buffer directly without copy-on-write bookkeeping beyond the
 * usual expr_unshare at the Expr level. */
typedef struct {
    int rank;        /* >= 1 */
    int64_t* dims;   /* rank entries, owned */
    void* data;      /* owned, row-major; element type given by dtype */
    NDType dtype;    /* element data type */
} NDArrayData;

typedef struct Expr {
    ExprType type;
    /* Reference count for shared-ownership / copy-on-write semantics
     * (M3 milestone). Every successful expr_new_* / expr_copy returns a
     * node with refcount == 1; expr_ref bumps it; expr_free decrements
     * and only physically frees on transition to 0. Atoms are eligible
     * to be shared today (see expr_copy); compound (FUNCTION) nodes
     * remain deep-copied for now and therefore always have refcount==1
     * in current code paths. */
    unsigned refcount;
    /* M3 phase-3 evaluation timestamp. The value of `eval_clock_get()` at
     * the moment this node was last evaluated to a fixed point. Fresh
     * (parser- or builtin-constructed) nodes carry 0, which never matches
     * the live clock (clock starts at 1) so they always evaluate on the
     * first call. Once `evaluate(e)` reaches fixed point it stamps the
     * result with the current clock; a subsequent `evaluate(e)` with the
     * same clock returns an inc-ref'd view immediately, skipping the
     * outer fixed-point loop and all of `evaluate_step`. The clock is
     * bumped by every symbol-table mutation (Set, SetDelayed, Clear,
     * SetAttributes, ClearAttributes), so any user definition change
     * cleanly invalidates every cached evaluation in one shot. The
     * field is benign metadata and is intentionally NOT considered by
     * `expr_eq` / `expr_hash`; updating it on a shared node is safe. */
    uint64_t last_evaluated_at;
    union {
        int64_t integer;
        double real;
        /* Phase 3b (EVAL_SYMTAB_IMPROVEMENTS): an EXPR_SYMBOL carries the
         * interned `name` AND a lazily-resolved, cached pointer to its
         * definition cell, so the evaluator reaches attributes / DownValues /
         * builtin by pointer instead of a symbol-table lookup. `name` stays the
         * FIRST member (union offset 0), so the historical EXPR_STRING type-pun
         * (reading `data.symbol.name.name` for a string yields `data.string`) still
         * holds. `def` is benign, lazily-filled metadata (NULL until first
         * resolve), NOT considered by expr_eq/expr_hash, and MUST be reset to
         * NULL wherever `name` is reassigned in place. */
        struct {
            char* name;
            struct SymbolDef* def;
        } symbol;
        char* string;
        struct {
            struct Expr* head;
            struct Expr** args;
            size_t arg_count;
        } function;
        mpz_t bigint;
        NDArrayData ndarray;
#ifdef USE_MPFR
        mpfr_t mpfr;          /* carries its own precision in bits */
#endif
    } data;
} Expr;

Expr* expr_new_integer(int64_t value);
Expr* expr_new_real(double value);
Expr* expr_new_symbol(const char* name);
Expr* expr_new_string(const char* str);
Expr* expr_new_function(Expr* head, Expr** args, size_t arg_count);
void expr_free(Expr* e);
Expr* expr_copy(Expr* e);
/* Inc-ref `e` (no copy) and return the same pointer. Pairs with
 * expr_free, which dec-refs. NULL is passed through. Use this in places
 * that previously called expr_copy on an immutable atom. */
Expr* expr_ref(Expr* e);

/* M3 phase-2 copy-on-write helper. Returns a refcount==1 Expr that is
 * logically equal to `e`, consuming the input ref. If `e->refcount` is
 * already 1, returns `e` unchanged (zero work). Otherwise builds a
 * one-level private copy: a fresh top-level node with a private args[]
 * array (FUNCTION) or owned mpz_t/mpfr_t/string payload (others), with
 * children inc-ref'd via expr_copy. The caller may then mutate the
 * returned node's *direct* fields (args[i], arg_count, head, etc.).
 * To mutate deeper than one level, unshare each interior node along
 * the path. NULL passes through. */
Expr* expr_unshare(Expr* e);

#include <stdbool.h>
bool expr_eq(const Expr* a, const Expr* b);
int expr_compare(const Expr* a, const Expr* b);
uint64_t expr_hash(const Expr* e);

/* BigInt constructors */
Expr* expr_new_bigint_from_mpz(const mpz_t val);
Expr* expr_new_bigint_from_int64(int64_t val);
Expr* expr_new_bigint_from_str(const char* str);

/* NDArray constructor. Takes ownership of `dims` (copied internally, caller
 * keeps its own copy) and `data` (moved, not copied — caller must not free
 * or reuse it afterwards). `rank` >= 1, `dims[i]` >= 1, `data` must have
 * exactly product(dims) elements of `dtype` (ndt_elem_size(dtype) bytes each). */
Expr* expr_new_ndarray(int rank, const int64_t* dims, void* data, NDType dtype);

#ifdef USE_MPFR
/* MPFR constructors. Each allocates an Expr whose payload `mpfr_t` is
 * initialized to the requested precision. Caller must `expr_free` when
 * done, which calls `mpfr_clear`. */
Expr* expr_new_mpfr_bits(mpfr_prec_t bits);                       /* zero */
Expr* expr_new_mpfr_from_d(double v, mpfr_prec_t bits);
Expr* expr_new_mpfr_from_si(long v, mpfr_prec_t bits);
Expr* expr_new_mpfr_from_mpz(const mpz_t z, mpfr_prec_t bits);
Expr* expr_new_mpfr_from_str(const char* str, mpfr_prec_t bits);
/* Build an Expr taking ownership of `src`. The mpfr_t is moved, not copied;
 * afterwards the caller must not touch `src`. Precision is inherited. */
Expr* expr_new_mpfr_move(mpfr_t src);
/* Copy constructor: new Expr with an independent mpfr_t at the same
 * precision as `src`. */
Expr* expr_new_mpfr_copy(const mpfr_t src);
#endif

/* Portable strdup: `strdup` is POSIX, not C99, and glibc hides it under
 * -std=c99.  Returns a malloc'd copy of `s` (NULL if `s` is NULL). Caller
 * frees with free(). Use this instead of strdup throughout the codebase. */
char* mathilda_strdup(const char* s);

/* Portable "this static function/variable may be unused" marker. Wraps the
 * GNU/Clang attribute behind a guard so strict C99 elsewhere still builds
 * (SPEC §10 forbids an unguarded __attribute__). */
#if defined(__GNUC__) || defined(__clang__)
#define MATHILDA_MAYBE_UNUSED __attribute__((unused))
#else
#define MATHILDA_MAYBE_UNUSED
#endif

/* Helpers used by arithmetic modules */
void  expr_to_mpz(const Expr* e, mpz_t out);
bool  expr_is_integer_like(const Expr* e);
/* True if `e` represents any concrete number: Integer, BigInt, Real,
 * Rational[n,d], Complex with numeric parts, or (with USE_MPFR) MPFR. */
bool  expr_is_numeric_like(const Expr* e);
Expr* expr_bigint_normalize(Expr* e);

#endif // EXPR_H
