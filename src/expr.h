
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
    EXPR_NDARRAY,              /* dense machine-precision ndarray, see ndarray.h */
    EXPR_COMPILED              /* compiled numeric function, see compile/compiled_function.h */
#ifdef USE_MPFR
    , EXPR_MPFR                /* arbitrary-precision real (MPFR) */
#endif
} ExprType;

/* Element data type of an EXPR_NDARRAY, analogous to numpy's dtype. The four
 * floating types map one-to-one onto BLAS's four precisions (s/d/c/z), where the
 * numeric suffix is the *per-component* bit width:
 *   NDT_FLOAT64   double            (BLAS d) — DEFAULT, value 0 so a
 *                                   zero-initialized/legacy struct behaves as
 *                                   today's double-only NDArray.
 *   NDT_FLOAT32   float             (BLAS s)
 *   NDT_COMPLEX64 2x double (re,im) (BLAS z)
 *   NDT_COMPLEX32 2x float  (re,im) (BLAS c)
 * Complex is stored as interleaved (re, im) plain floats/doubles — no
 * <complex.h>, matching the existing BLAS interleaved-buffer convention and
 * keeping the code strict-C99 portable. See ndarray.h for the ndt_* helpers.
 *
 *   NDT_INT64     int64_t           — COMPILER-INTERNAL, no BLAS counterpart.
 *
 * NDT_INT64 exists so Compile[] can hold a packed buffer of machine INTEGERS:
 * `Range[n]` and `Table[i, {i, 1, n}]` are exact Integers in the interpreter, and
 * before this there was no integer dtype to put them in, so compiling them would
 * have answered with a different element type rather than merely faster.
 *
 * It is deliberately NOT reachable from user syntax — `ndt_from_string` does not
 * produce it and `NDArray[...]` never infers it — and the Compile[] boundary
 * always unpacks an int64 buffer to a List of Integers, so one never escapes
 * into the interpreter. That containment is what lets the generic ndt_get/ndt_set
 * pair keep its `double` signature: those two are exact only to 2^53, and the
 * exact accessors (ndt_get_i/ndt_set_i) are used everywhere an int64 array is
 * actually read or written. */
typedef enum {
    NDT_FLOAT64 = 0,
    NDT_FLOAT32,
    NDT_COMPLEX64,
    NDT_COMPLEX32,
    NDT_INT64
} NDType;

/* How an EXPR_NDARRAY presents itself to the user. ONE node type, TWO surfaces.
 *
 *   NDA_HEAD_NDARRAY  the explicit `NDArray[...]` a user typed. Head is
 *                     NDArray, prints as NDArray[{...}], AtomQ is True, ListQ
 *                     is False, and a Part assignment COERCES the rhs into the
 *                     machine buffer (the user asked for a machine buffer).
 *
 *   NDA_HEAD_LIST     a PACKED LIST: an ordinary List that the system chose to
 *                     store as a dense buffer, or that ToNDArray[] packed on
 *                     request. Indistinguishable from the equivalent nested
 *                     List -- Head is List, prints as {1., 2.}, ListQ/VectorQ/
 *                     MatrixQ are True, AtomQ is False, patterns match, and
 *                     expr_eq/expr_hash/expr_compare agree with the plain List
 *                     element for element. Only NDArrayQ[] tells them apart.
 *                     A Part assignment that is not exactly representable at
 *                     the buffer's dtype (INCLUDING its head: writing an exact
 *                     Integer into a float64 buffer counts) UNPACKS first, so
 *                     packing can never change a value or a head.
 *
 * Zero is NDA_HEAD_NDARRAY so every pre-existing construction site keeps its
 * historical meaning. Deriving a new array FROM an existing one must inherit
 * the field -- use expr_new_ndarray_like, not expr_new_ndarray_raw, or a packed
 * list turns visible at the first operation applied to it.
 *
 * This IS part of a value's identity, unlike `refcount` or `last_evaluated_at`:
 * it decides the value's Head, and two values with different Heads are never
 * SameQ. So `NDArray[{1.,2.}] === {1.,2.}` is False (NDArray vs List) while
 * `ToNDArray[{1.,2.}] === {1.,2.}` is True. expr_eq, expr_hash and expr_compare
 * each fork on it:
 *   - packed vs packed     elementwise; dtype is NOT observable (a packed list
 *                          of float32 equals the float64 one with equal values,
 *                          because both materialise to the same Reals)
 *   - packed vs plain List elementwise against the List's leaves
 *   - visible vs anything  today's rule: dtype and raw bytes are identity
 * and expr_hash of a packed list must equal expr_hash of the List it
 * materialises to, or Association keys and Union stop working.
 * See docs/design/packed_arrays.md. */
typedef enum {
    NDA_HEAD_NDARRAY = 0,
    NDA_HEAD_LIST
} NDPresentation;

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
    NDPresentation present_as;  /* NDArray[...] or a packed List; see above */
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
        /* EXPR_COMPILED: opaque, reference-counted, immutable-after-build
         * compiled-function payload.  Defined in compile/compiled_function.h;
         * only a forward declaration is visible here to keep the layering clean
         * (a bare pointer never grows the union). */
        struct CompiledFunction* compiled;
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

/* Drain the internal Expr free-list, returning every recycled node to the
 * system allocator. Registered with atexit() automatically on first
 * allocation, so callers rarely need it; exposed for tests and for an
 * explicit end-of-run cleanup. Safe to call repeatedly. Nodes still in use
 * are unaffected — only already-freed (pooled) nodes are released. */
void expr_pool_free_all(void);
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
 * exactly product(dims) elements of `dtype` (ndt_elem_size(dtype) bytes each).
 *
 * The result presents as NDA_HEAD_NDARRAY — a VISIBLE `NDArray[...]`. That is
 * correct only for a genuine constructor, i.e. one building an array from
 * something that was not already an array: `NDArray[...]` itself, a producer
 * like Range/Table, an Import, the Compile[] boundary.
 *
 * If you are computing a NEW array FROM an existing one — elementwise map,
 * reduction, reshape, sub-array, transpose, anything — use
 * expr_new_ndarray_like instead. Using this function there silently turns a
 * packed List back into a visible NDArray at the first operation applied to it,
 * which is a user-visible wrong answer with no crash to find it by. The name
 * carries `_raw` so the choice has to be made deliberately at every site. */
Expr* expr_new_ndarray_raw(int rank, const int64_t* dims, void* data, NDType dtype);

/* Derive a new NDArray FROM `src`, inheriting its presentation (visible
 * NDArray[...] vs packed List). Identical to expr_new_ndarray_raw otherwise;
 * `src` is only read for its presentation and may be NULL (treated as raw). */
Expr* expr_new_ndarray_like(const Expr* src, int rank, const int64_t* dims,
                            void* data, NDType dtype);

/* EXPR_COMPILED constructor.  Takes ownership of one reference to `cf` (the
 * caller must not free it afterwards; expr_free will dec-ref).  See
 * compile/compiled_function.h. */
Expr* expr_new_compiled(struct CompiledFunction* cf);
/* A stable per-object identity (used by expr_eq/hash/compare for EXPR_COMPILED).
 * Implemented in compile/compiled_function.c. */
uint64_t compiled_function_identity(const struct CompiledFunction* cf);

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
