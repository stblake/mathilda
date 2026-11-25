
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
    EXPR_COMPILED,             /* compiled numeric function, see compile/compiled_function.h */
    EXPR_SPARSE,               /* SparseArray raw object, see sparse/sparse.h */
    EXPR_LINKEDLIST,           /* mutable LinkedList object, see linkedlist.h */
    EXPR_PVEC                  /* persistent vector value, see pvec.h */
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
 * actually read or written.
 *
 *   NDT_BOOL      uint8_t (0/1)     — no BLAS counterpart.
 *
 * NDT_BOOL is one byte per element and materialises as the symbols True/False
 * (ndarray_buffer_element_to_expr), so a bool array is observably the List of
 * True/False it packs. Unlike NDT_INT64 it IS reachable from user syntax
 * (`DataType -> "bool"`/`"Boolean"`), is inferred by the auto-packer from an
 * all-True/False list, and is produced by the elementwise sign predicates. It is
 * deliberately NOT numeric: arithmetic over a bool array delists to the symbolic
 * List (Mathematica leaves `True + True` unevaluated), and the compiler mirrors
 * this with its "Bool is not numeric" type-lattice invariant. */
typedef enum {
    NDT_FLOAT64 = 0,
    NDT_FLOAT32,
    NDT_COMPLEX64,
    NDT_COMPLEX32,
    NDT_INT64,
    NDT_BOOL
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
    /* Layout: this struct is the largest member of Expr's union, so its size
     * IS sizeof(Expr) - 24 (tests/test_symsum.c pins 56). dtype and present_as
     * are stored in one byte each, beside rank, so cap_rows fits in 32 bytes;
     * they hold NDType / NDPresentation values (read and assign as usual). */
    int rank;        /* >= 1 */
    uint8_t dtype;      /* NDType: element data type */
    uint8_t present_as; /* NDPresentation: NDArray[...] or a packed List; see above */
    int64_t* dims;   /* rank entries, owned */
    void* data;      /* owned, row-major; element type given by dtype */
    /* Allocated length of `data` in ROWS (dims[0] units), 0 = exactly dims[0].
     * Only an in-place rank-1 AppendTo (ndstruct_append_inplace) ever sets it,
     * to give amortised O(1) growth; every other site leaves it 0 and treats
     * dims[0] as the buffer length, which stays true (cap >= dims[0] always).
     * A copy of the node gets an exact buffer and cap 0 (expr_unshare). */
    size_t cap_rows;
} NDArrayData;

typedef struct Expr {
    ExprType type;
    /* Reference count for shared-ownership / copy-on-write semantics
     * (M3 milestone). Every successful expr_new_* / expr_copy returns a
     * node with refcount == 1; expr_ref bumps it; expr_free decrements
     * and only physically frees on transition to 0. Atoms AND compound
     * (FUNCTION) nodes alike are shared by inc-ref in expr_copy today; a
     * node with refcount > 1 is immutable, so any mutator calls
     * expr_unshare() to obtain a private (refcount==1) copy before
     * rewriting fields in place. That immutable-while-shared invariant is
     * what makes the cached hash below safe. */
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
     * `expr_eq` / `expr_hash`; updating it on a shared node is safe.
     *
     * The TOP bit (1<<63) is a benign GROUND flag (see eval.c): set on a
     * fixed point built solely from literals under the six pure structural
     * constructors, it lets the node be re-validated as a fixed point after
     * the eval clock churns (loop-invariant O(1) re-check). The clock is a
     * monotone counter that never reaches 2^63, so the low 63 bits are the
     * stamp; every eval-clock comparison masks the flag off. Read the field
     * through eval_node_stamp() / eval_node_is_ground(), not directly. */
    uint64_t last_evaluated_at;
    /* Lazily-memoized structural hash (== expr_hash(this) recomputed fresh).
     * 0 means "not yet computed" (a genuine hash of 0 simply never cache-hits,
     * which is harmless — the value returned is always the true hash). expr_hash
     * fills it on first call and reuses it thereafter, turning the per-element
     * subtree walks in Plus/Times grouping, Association indexing and set ops into
     * O(1) reads. Like last_evaluated_at / symbol.def it is benign metadata NOT
     * considered by expr_eq/expr_compare, and — because the cache is a function
     * of the node's STRUCTURE — it MUST be reset to 0 (expr_invalidate_hash)
     * anywhere a live node's head/args/arg_count/payload is rewritten in place.
     * Fresh nodes get 0 from expr_alloc_node; expr_unshare's private copy starts
     * at 0. The MATHILDA_HASH_VERIFY build recomputes and asserts every hit to
     * catch a missed invalidation deterministically. */
    uint64_t hash_cache;
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
            /* Read through expr_fn_assoc_index(): this slot also carries the
             * scoping walkers' binder summary, as an odd word (symsum.h).
             * Persistent key->position hash index, non-NULL ONLY for canonical
             * Association nodes (built at canonicalisation; see assoc_index.h).
             * Benign acceleration metadata: NULL for every other function node,
             * ignored by expr_eq/expr_hash/expr_compare, freed with the node,
             * and reset to NULL on any physical copy (expr_unshare). Occupies
             * the union's existing slack (NDArrayData/mpfr_t are larger), so
             * sizeof(Expr) is unchanged. */
            struct AssocIndex* index;
        } function;
        mpz_t bigint;
        NDArrayData ndarray;
        /* EXPR_COMPILED: opaque, reference-counted, immutable-after-build
         * compiled-function payload.  Defined in compile/compiled_function.h;
         * only a forward declaration is visible here to keep the layering clean
         * (a bare pointer never grows the union). */
        struct CompiledFunction* compiled;
        /* EXPR_SPARSE: opaque, reference-counted, immutable SparseArray
         * payload (sparse/sparse.h); a bare pointer never grows the union. */
        struct SparseArray* sparse;
        /* EXPR_LINKEDLIST: the ONE mutable payload. Reference-counted and
         * shared by every node that refers to it (expr_unshare shares it too):
         * mutators write through it, so all aliases see the change. Identity
         * is the payload (expr_eq/hash/compare use its stable id), which is
         * why a mutation never invalidates a cached hash. linkedlist.h. */
        struct LinkedList* llist;
        /* EXPR_PVEC: a persistent (immutable, structurally shared) vector.
         * Reference-counted; copying a node retains the payload, and every
         * update builds a new payload that shares untouched nodes. pvec.h. */
        struct PVec* pvec;
#ifdef USE_MPFR
        mpfr_t mpfr;          /* carries its own precision in bits */
#endif
    } data;
} Expr;

Expr* expr_new_integer(int64_t value);
Expr* expr_new_real(double value);
Expr* expr_new_symbol(const char* name);
/* A symbol node for a SymbolDef already in hand (no interning lookup). */
struct SymbolDef;
Expr* expr_new_symbol_def(struct SymbolDef* d);
Expr* expr_new_string(const char* str);
Expr* expr_new_string_adopt(char* buf);              /* takes ownership of buf */
Expr* expr_new_string_n(const char* s, size_t n);    /* copy of s[0, n) */
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
/* Does `e` (heads included) mention a symbol named names[i] for some i < n?
 * Names are interned, so this compares pointers.  Iterative: constant C stack
 * for any depth of nesting. */
bool expr_mentions_any(const Expr* e, const char* const* names, size_t n);
int expr_compare(const Expr* a, const Expr* b);
uint64_t expr_hash(const Expr* e);

/* Drop a node's memoized structural hash (a pure function of the node's
 * structure). MUST be
 * called by any code that rewrites a live node's head/args/arg_count/scalar
 * payload IN PLACE (rather than building a fresh node), so a later expr_hash /
 * Orderless comparison recomputes instead of using the pre-mutation value. A
 * no-op on NULL and on nodes never hashed/compared. */
static inline void expr_invalidate_hash(Expr* e) {
    if (!e) return;
    e->hash_cache = 0;
    /* The binder summary (symsum.h) describes the same structure: drop it too
     * (an odd word in the index slot; an AssocIndex pointer is left alone). */
    if (e->type == EXPR_FUNCTION && ((uintptr_t)e->data.function.index & 1))
        e->data.function.index = NULL;
}

/* The `data.function.index` slot holds one of three things: NULL; an
 * AssocIndex pointer (canonical Associations only; malloc-aligned, so even);
 * or a binder summary (symsum.h; odd).  Read and write the AssocIndex only
 * through these. */
static inline struct AssocIndex* expr_fn_assoc_index(const Expr* e) {
    uintptr_t w = (uintptr_t)e->data.function.index;
    return (w & 1) ? NULL : (struct AssocIndex*)w;
}
static inline void expr_fn_set_assoc_index(Expr* e, struct AssocIndex* ix) {
    e->data.function.index = ix;
}

/* BigInt constructors */
Expr* expr_new_bigint_from_mpz(const mpz_t val);
Expr* expr_new_bigint_from_int64(int64_t val);
Expr* expr_new_bigint_from_str(const char* str);
/* Move constructor: the node takes over `val`'s limbs (no copy). `val` is
 * consumed -- the caller must neither mpz_clear nor reuse it. The result is
 * NOT normalised; wrap in expr_bigint_normalize for a value that may fit. */
Expr* expr_new_bigint_take(mpz_t val);

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
/* Bytes a NDArray node accounts for under MemoryConstrained (buffer + dims). */
int64_t expr_ndarray_payload_bytes(const Expr* e);
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
/* The structural hash of a CompiledFunction (used by expr_hash/compare for
 * EXPR_COMPILED), and its value equality (expr_eq): the same signature, body
 * and options.  Implemented in compile/compiled_function.c. */
uint64_t compiled_function_identity(const struct CompiledFunction* cf);
bool compiled_function_same(const struct CompiledFunction* a, const struct CompiledFunction* b);

/* EXPR_SPARSE constructor. Takes ownership of one reference to `sa` (expr_free
 * releases it). See sparse/sparse.h. */
Expr* expr_new_sparse(struct SparseArray* sa);

/* EXPR_LINKEDLIST constructor. Takes ownership of one reference to `ll`
 * (expr_free releases it). See linkedlist.h. */
Expr* expr_new_linkedlist(struct LinkedList* ll);

/* EXPR_PVEC constructor. Takes ownership of one reference to `pv`. See pvec.h. */
Expr* expr_new_pvec(struct PVec* pv);
/* Has any PVector value been created yet?  (Lets the evaluator skip its
 * PVector gate with one load until the first one exists.) */
extern bool g_pvec_created;

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
/* Install the $MinPrecision / $MaxPrecision width clamp applied by every
 * arbitrary-precision constructor (numeric.c, numeric_mpfr_prec_clamp). */
void  expr_set_mpfr_prec_clamp(mpfr_prec_t (*clamp)(mpfr_prec_t));
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
/* Read-only mpz view of an integer-like `e` WITHOUT copying a BigInt's limbs:
 * returns e's own mpz for a BigInt, or `scratch` initialised from an Integer.
 * Pair every call with expr_mpz_release(e, scratch). expr_to_mpz copies the
 * whole number, which for a several-thousand-digit operand of a + b costs
 * more than the addition itself. */
mpz_srcptr expr_mpz_borrow(const Expr* e, mpz_t scratch);
void  expr_mpz_release(const Expr* e, mpz_t scratch);
bool  expr_is_integer_like(const Expr* e);
/* True if `e` represents any concrete number: Integer, BigInt, Real,
 * Rational[n,d], Complex with numeric parts, or (with USE_MPFR) MPFR. */
bool  expr_is_numeric_like(const Expr* e);
/* True and *out = z when z fits int64 (portable: not via `long`). */
bool  expr_mpz_get_int64(const mpz_t z, int64_t* out);
/* True for an Integer or BigInt e, with *out its value saturated to
 * [-INT64_MAX, INT64_MAX]: a BigInt beyond int64 reads as the bound of its
 * sign.  For a count, length or index argument whose every value past the
 * int64 range means the same thing as the bound ("at most 10^20 edges",
 * "the 2^70-th diagonal"), so a huge spec answers instead of staying
 * unevaluated.  -INT64_MAX, not INT64_MIN, so negating *out never overflows.
 * False (and *out untouched) for anything else. */
bool  expr_get_int64_saturated(const Expr* e, int64_t* out);
/* z = v for any int64 v.  mpz_set_si/mpz_init_set_si take a `long`, which is
 * only 32 bits on LLP64 Windows, so an int64 above 2^31 would be truncated;
 * these are exact on every data model. */
void  mth_mpz_set_i64(mpz_t z, int64_t v);
void  mth_mpz_init_set_i64(mpz_t z, int64_t v);
/* The rest of the int64/uint64 siblings of GMP's `_si`/`_ui` calls, for the
 * same reason: on LLP64 a `long` operand is 32 bits.  Each is the plain
 * `_si`/`_ui` call where `long` is 64 bits (and, on LLP64, whenever the value
 * fits a `long`), otherwise exact through an mpz.  mth_mpq_set_i64 does not
 * canonicalise, exactly like mpq_set_si. */
bool  mth_u64_fits_ulong(uint64_t v);      /* v fits an `unsigned long` */
void  mth_mpz_set_u64(mpz_t z, uint64_t v);
void  mth_mpz_init_set_u64(mpz_t z, uint64_t v);
void  mth_mpz_mul_i64(mpz_t r, const mpz_t a, int64_t v);
void  mth_mpz_mul_u64(mpz_t r, const mpz_t a, uint64_t v);
void  mth_mpz_add_u64(mpz_t r, const mpz_t a, uint64_t v);
void  mth_mpz_sub_u64(mpz_t r, const mpz_t a, uint64_t v);
int   mth_mpz_cmp_i64(const mpz_t a, int64_t v);
int   mth_mpz_cmp_u64(const mpz_t a, uint64_t v);
int   mth_mpz_divisible_u64_p(const mpz_t a, uint64_t d);
void  mth_mpz_divexact_u64(mpz_t r, const mpz_t a, uint64_t d);
uint64_t mth_mpz_fdiv_u64(const mpz_t a, uint64_t d);   /* a mod d in [0, d) */
uint64_t mth_mpz_get_u64(const mpz_t z);   /* low 64 bits of |z|, as mpz_get_ui */
void  mth_mpq_set_i64(mpq_t q, int64_t num, uint64_t den);
void  mth_mpq_set_u64(mpq_t q, uint64_t num, uint64_t den);
#ifdef USE_MPFR
int   mth_mpfr_set_i64(mpfr_t x, int64_t v, mpfr_rnd_t rnd);
int   mth_mpfr_set_u64(mpfr_t x, uint64_t v, mpfr_rnd_t rnd);
int   mth_mpfr_mul_i64(mpfr_t r, const mpfr_t a, int64_t v, mpfr_rnd_t rnd);
int   mth_mpfr_div_i64(mpfr_t r, const mpfr_t a, int64_t v, mpfr_rnd_t rnd);
int   mth_mpfr_add_i64(mpfr_t r, const mpfr_t a, int64_t v, mpfr_rnd_t rnd);
int   mth_mpfr_pow_i64(mpfr_t r, const mpfr_t a, int64_t v, mpfr_rnd_t rnd);
#endif
Expr* expr_bigint_normalize(Expr* e);

#endif // EXPR_H
