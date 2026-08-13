/* Automatic packed arrays — the packing decision, the producer API, and the
 * ToNDArray/FromNDArray builtins. See pack.h for the model and the contract. */

#include "pack.h"
#include "ndarray.h"
#include "symtab.h"
#include "sym_names.h"
#include "attr.h"
#include "eval.h"
#include "common.h"   /* head_is */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

/* --------------------------------------------------------------------------
 * Switches
 * -------------------------------------------------------------------------- */

/* -1 = not yet read from the environment. Matches src/numloop.c's shape. */
static int    g_pack_enabled = -1;
static size_t g_pack_min = 0;          /* 0 => PACK_MIN_ELEMENTS */

bool pack_g_any_created = false;

bool pack_enabled(void) {
    if (g_pack_enabled < 0)
        g_pack_enabled = (getenv("MATHILDA_NO_PACK") == NULL) ? 1 : 0;
    return g_pack_enabled != 0;
}

void pack_set_enabled(bool on) { g_pack_enabled = on ? 1 : 0; }

void   pack_set_min_elements(size_t n) { g_pack_min = n; }

/* MATHILDA_PACK_MIN overrides the threshold for a whole session, the same
 * family of knob as MATHILDA_NO_PACK and MATHILDA_PACK_DIAG. It exists because
 * the threshold is the ONE packing parameter whose right value is empirical:
 * it trades the buffer path on small values against the cost of routing
 * symbolic and pattern-matching code through a second representation, and that
 * trade can only be swept, not derived. Read once, lazily; an unparseable or
 * zero value leaves the compiled-in default alone. The C setter still wins, so
 * the differential suite is unaffected. */
size_t pack_min_elements(void) {
    if (g_pack_min) return g_pack_min;
    static int env_min = -1;                       /* -1 = not yet read */
    if (env_min < 0) {
        const char* s = getenv("MATHILDA_PACK_MIN");
        long v = s ? strtol(s, NULL, 10) : 0;
        env_min = (v > 0 && v < 1000000000L) ? (int)v : 0;
    }
    return env_min ? (size_t)env_min : PACK_MIN_ELEMENTS;
}

/* --------------------------------------------------------------------------
 * Sniff: shape, rectangularity and element class in ONE pass, allocating
 * nothing and bailing at the first element that rules packing out.
 *
 * Deliberately not ndarray_from_nested_list's probe_dims: that takes the dtype
 * as a parameter (the caller must already know it), so inferring the dtype
 * around it would cost a third traversal. It is also more permissive than
 * packing can afford -- it accepts an Integer into a float dtype, which is
 * exactly the silent head change ({1, 2, 3} becoming {1., 2., 3.}) that a
 * representation choice may never make.
 * -------------------------------------------------------------------------- */

typedef enum { PK_NONE = 0, PK_INT, PK_REAL, PK_BOOL } PackClass;

static PackClass leaf_class(const Expr* e) {
    if (!e) return PK_NONE;
    if (e->type == EXPR_INTEGER) return PK_INT;
    if (e->type == EXPR_REAL)    return PK_REAL;
    /* The two boolean literals pack into a one-byte bool buffer; every other
     * symbol still declines. A bool list round-trips exactly (the buffer
     * materialises back to True/False), so unlike Complex it is safe to pack
     * automatically. */
    if (e->type == EXPR_SYMBOL &&
        (e->data.symbol.name == SYM_True || e->data.symbol.name == SYM_False))
        return PK_BOOL;
    /* Everything else declines, and each for a reason:
     *   BigInt / Rational / MPFR  no machine representation
     *   Complex                   materialises as Complex[re, im] even when
     *                             im is 0, which the evaluator never produces
     *                             (it folds to a Real) -- so a packed complex
     *                             list would not round-trip. Explicit
     *                             NDArray[..., DataType -> "complex64"] still
     *                             works; packing it is a follow-up that needs
     *                             a presentation-aware zero-imaginary fold.
     *   other symbols / expressions  not machine values at all */
    return PK_NONE;
}

/* Returns rank (>= 1) on success, -1 when not packable. `*cls` accumulates the
 * element class and must start at PK_NONE. Mixed Integer/Real declines: a
 * uniform buffer cannot hold an exact head on one element and an inexact head
 * on another (src/numloop.h states the same rule for compiled loops). */
static int pack_sniff(const Expr* e, int64_t* dims, int maxrank, PackClass* cls) {
    if (!e) return -1;
    /* An already-packed row. Reached when a producer assembles a List of results
     * that each packed on their own -- Table[i j, {i,300}, {j,300}] evaluates the
     * inner Table per outer index and gets 300 packed rows back. Absorbing them
     * gives the rank-2 array that shape deserves; refusing would leave a List of
     * buffers, which the gate then has to materialise row by row. */
    if (is_packed_list(e)) {
        int rank = e->data.ndarray.rank;
        if (rank > maxrank) return -1;
        PackClass c = (e->data.ndarray.dtype == NDT_BOOL)  ? PK_BOOL
                    : (e->data.ndarray.dtype == NDT_INT64) ? PK_INT : PK_REAL;
        if (*cls == PK_NONE) *cls = c;
        else if (*cls != c) return -1;          /* mixed exact/inexact */
        for (int i = 0; i < rank; i++) dims[i] = e->data.ndarray.dims[i];
        return rank;
    }
    if (e->type != EXPR_FUNCTION ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_List) {
        PackClass c = leaf_class(e);
        if (c == PK_NONE) return -1;
        if (*cls == PK_NONE) *cls = c;
        else if (*cls != c) return -1;      /* mixed exact/inexact */
        return 0;                            /* a leaf: rank 0 below this axis */
    }
    size_t len = e->data.function.arg_count;
    if (len == 0 || maxrank <= 0) return -1;         /* {} never packs */
    int64_t sub[NDARRAY_MAX_RANK];
    int sub_rank = pack_sniff(e->data.function.args[0], sub, maxrank - 1, cls);
    if (sub_rank < 0) return -1;
    for (size_t i = 1; i < len; i++) {
        int64_t other[NDARRAY_MAX_RANK];
        int r = pack_sniff(e->data.function.args[i], other, maxrank - 1, cls);
        if (r != sub_rank) return -1;                /* ragged */
        for (int j = 0; j < r; j++)
            if (other[j] != sub[j]) return -1;       /* ragged */
    }
    dims[0] = (int64_t)len;
    for (int j = 0; j < sub_rank; j++) dims[j + 1] = sub[j];
    return sub_rank + 1;
}

/* Second pass: write the leaves into the buffer in row-major order. The shape
 * was already validated, so this cannot fail. */
static void pack_flatten(const Expr* e, void* buf, NDType dt, size_t* k) {
    if (is_packed_list(e)) {
        /* Row already in a buffer. Same dtype (pack_sniff made the classes
         * agree) copies straight across; a widening int64 -> float64 goes
         * element by element through the exact accessor. */
        size_t n = 1;
        for (int i = 0; i < e->data.ndarray.rank; i++) n *= (size_t)e->data.ndarray.dims[i];
        const void* src = e->data.ndarray.data;
        NDType sdt = e->data.ndarray.dtype;
        if (sdt == dt) {
            memcpy((char*)buf + ndt_elem_size(dt) * (*k), src, ndt_elem_size(dt) * n);
            *k += n;
        } else if (sdt == NDT_INT64) {
            for (size_t i = 0; i < n; i++)
                ndt_set(buf, (*k)++, dt, (double)ndt_get_i(src, i, sdt), 0.0);
        } else {
            for (size_t i = 0; i < n; i++) {
                double re, im;
                ndt_get(src, i, sdt, &re, &im);
                ndt_set(buf, (*k)++, dt, re, im);
            }
        }
        return;
    }
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_List) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            pack_flatten(e->data.function.args[i], buf, dt, k);
        return;
    }
    if (dt == NDT_BOOL)
        ndt_set_i(buf, (*k)++, dt, e->data.symbol.name == SYM_True ? 1 : 0);
    else if (dt == NDT_INT64) ndt_set_i(buf, (*k)++, dt, e->data.integer);
    else ndt_set(buf, (*k)++, dt,
                 (e->type == EXPR_INTEGER) ? (double)e->data.integer : e->data.real,
                 0.0);
}

/* The shared core. `min_elems` of 0 means "no threshold" (ToNDArray).
 * `want` overrides the inferred dtype when `have_want` is set. Borrows `list`;
 * returns a new packed list or NULL. */
static Expr* pack_build(const Expr* list, size_t min_elems,
                        bool have_want, NDType want) {
    int64_t dims[NDARRAY_MAX_RANK];
    PackClass cls = PK_NONE;
    int rank = pack_sniff(list, dims, NDARRAY_MAX_RANK, &cls);
    if (rank < 1 || cls == PK_NONE) return NULL;

    size_t n = 1;
    for (int i = 0; i < rank; i++) n *= (size_t)dims[i];
    if (n < min_elems) return NULL;      /* declined before any allocation */

    NDType dt = (cls == PK_BOOL) ? NDT_BOOL
              : (cls == PK_INT)  ? NDT_INT64 : NDT_FLOAT64;
    if (have_want) {
        /* An explicit DataType may WIDEN an exact list to a float buffer -- the
         * user asked for that -- but never the reverse: rounding Reals into an
         * int64 buffer would change values, not just storage. A bool buffer is
         * neither a widening nor a narrowing of a number, so it is only legal
         * over an all-True/False list (and, conversely, a True/False list only
         * packs to bool). */
        if (want == NDT_INT64 && cls != PK_INT) return NULL;
        if (ndt_is_complex(want)) return NULL;   /* see leaf_class */
        if ((want == NDT_BOOL) != (cls == PK_BOOL)) return NULL;
        dt = want;
    }

    void* buf = malloc(ndt_elem_size(dt) * n);
    if (!buf) return NULL;
    size_t k = 0;
    pack_flatten(list, buf, dt, &k);

    Expr* out = expr_new_ndarray_raw(rank, dims, buf, dt);
    if (!out) { free(buf); return NULL; }
    out->data.ndarray.present_as = NDA_HEAD_LIST;
    pack_g_any_created = true;
    return out;
}

Expr* pack_offer(Expr* list) {
    if (!list || !pack_enabled()) return list;
    /* O(1) pre-test so a short list is refused without walking anything: either
     * it is long enough on its own, or it is nested and its rows might be. */
    if (list->type != EXPR_FUNCTION) return list;
    size_t len = list->data.function.arg_count;
    size_t min = pack_min_elements();
    if (len < min) {
        if (len < 2) return list;
        /* Short on its own, but a nested one may still clear the threshold on
         * total elements -- 3 x 100000 is the best-value shape in the space. A
         * packed first element counts as nested too (a List of packed rows). */
        const Expr* first = list->data.function.args[0];
        if (!first || (first->type != EXPR_FUNCTION && !is_packed_list(first)))
            return list;
    }
    Expr* packed = pack_build(list, min, false, NDT_FLOAT64);
    if (!packed) return list;
    expr_free(list);
    return packed;
}

Expr* pack_force(Expr* list, bool have_dtype, NDType dtype) {
    if (!list) return list;
    Expr* packed = pack_build(list, 0, have_dtype, dtype);
    if (!packed) return list;
    expr_free(list);
    return packed;
}

Expr* pack_repack_like(const Expr* src, Expr* list) {
    if (!list || !is_packed_list(src)) return list;
    /* No threshold: the input was already a buffer, so staying one costs nothing
     * and materialising it would be a silent downgrade. */
    Expr* packed = pack_build(list, 0, false, NDT_FLOAT64);
    if (!packed) return list;
    expr_free(list);
    return packed;
}

/* Pack every plain numeric List argument of a Listable call up to a buffer, so
 * the head's own NDArray kernel handles the call instead of apply_listable
 * materialising the buffer that is already there. See the call site in
 * evaluate_step for why this direction is the right one.
 *
 * All or nothing: a partial lift would leave a plain List behind, `listable_mixed`
 * would stay true, and the gate would then materialise the arrays we had just
 * built. Declines outright -- leaving the call exactly as it was, to thread --
 * when any of these hold:
 *
 *   - an operand is symbolic or otherwise not a machine scalar. Threading is
 *     what makes `packedA + plainB + x` work at all: hand the same call to
 *     ndarray_elementwise and it returns NULL on the symbolic term, and Plus
 *     warns and leaves the sum unevaluated.
 *   - a List does not pack (Rationals, symbols, mixed exact/inexact).
 *   - a List packs to int64 and the head is not exact on int64 buffers.
 *   - the packed shape is not one the head's kernel covers: equal to the widest
 *     array operand's, or -- with packed_broadcast_ok -- a prefix of it. */
void pack_lift_listable_args(Expr* call, bool broadcast_ok, bool int64_ok,
                             bool* changed) {
    if (!call || call->type != EXPR_FUNCTION) return;
    size_t n = call->data.function.arg_count;

    const Expr* ref = NULL;                      /* widest array operand */
    size_t nlists = 0;
    for (size_t i = 0; i < n; i++) {
        const Expr* a = call->data.function.args[i];
        if (!a) return;
        if (a->type == EXPR_NDARRAY) {
            if (!ref || a->data.ndarray.rank > ref->data.ndarray.rank) ref = a;
        } else if (head_is(a, SYM_List)) {
            nlists++;
        } else if (a->type != EXPR_INTEGER && a->type != EXPR_REAL) {
            return;
        }
    }
    if (!ref || nlists == 0) return;             /* nothing to lift toward */

    Expr** lifted = calloc(n, sizeof(Expr*));
    if (!lifted) return;
    bool ok = true;
    for (size_t i = 0; i < n && ok; i++) {
        Expr* a = call->data.function.args[i];
        if (!head_is(a, SYM_List)) continue;
        Expr* p = pack_force(expr_copy(a), false, NDT_FLOAT64);
        if (!p || p->type != EXPR_NDARRAY ||
            (p->data.ndarray.dtype == NDT_INT64 && !int64_ok)) {
            expr_free(p);
            ok = false;
            break;
        }
        int rp = p->data.ndarray.rank, rr = ref->data.ndarray.rank;
        bool shape_ok = (rp == rr) || (broadcast_ok && rp < rr);
        for (int k = 0; shape_ok && k < rp; k++)
            if (p->data.ndarray.dims[k] != ref->data.ndarray.dims[k]) shape_ok = false;
        if (!shape_ok) { expr_free(p); ok = false; break; }
        lifted[i] = p;
    }
    if (ok) {
        for (size_t i = 0; i < n; i++) {
            if (!lifted[i]) continue;
            expr_free(call->data.function.args[i]);
            call->data.function.args[i] = lifted[i];
            lifted[i] = NULL;
            *changed = true;
        }
    }
    for (size_t i = 0; i < n; i++) expr_free(lifted[i]);
    free(lifted);
}

Expr* pack_unpack(const Expr* e) {
    if (!is_ndarray(e)) return NULL;
    return ndarray_to_nested_list(e);
}

Expr* pack_eval_plain(Expr* e) {
    Expr* r = evaluate(e);
    if (!is_packed_list(r)) return r;
    Expr* plain = ndarray_to_nested_list(r);
    expr_free(r);
    return plain;
}

/* --------------------------------------------------------------------------
 * Direct construction
 * -------------------------------------------------------------------------- */

Expr* ndbuild_open(int rank, const int64_t* dims, NDType dtype, void** buf_out) {
    *buf_out = NULL;
    if (!pack_enabled() || rank < 1 || rank >= NDARRAY_MAX_RANK) return NULL;
    size_t n = 1;
    for (int i = 0; i < rank; i++) {
        if (dims[i] < 1) return NULL;
        n *= (size_t)dims[i];
    }
    if (n < pack_min_elements()) return NULL;
    void* buf = malloc(ndt_elem_size(dtype) * n);
    if (!buf) return NULL;
    Expr* e = expr_new_ndarray_raw(rank, dims, buf, dtype);
    if (!e) { free(buf); return NULL; }
    e->data.ndarray.present_as = NDA_HEAD_LIST;
    pack_g_any_created = true;
    *buf_out = buf;
    return e;
}

Expr* ndbuild_open_like(const Expr* src, int rank, const int64_t* dims,
                        NDType dtype, void** buf_out) {
    *buf_out = NULL;
    if (rank < 1 || rank >= NDARRAY_MAX_RANK) return NULL;
    size_t n = 1;
    for (int i = 0; i < rank; i++) {
        if (dims[i] < 1) return NULL;
        n *= (size_t)dims[i];
    }
    void* buf = malloc(ndt_elem_size(dtype) * n);
    if (!buf) return NULL;
    Expr* e = expr_new_ndarray_like(src, rank, dims, buf, dtype);
    if (!e) { free(buf); return NULL; }
    if (e->data.ndarray.present_as == NDA_HEAD_LIST) pack_g_any_created = true;
    *buf_out = buf;
    return e;
}

Expr* ndbuild_open_f64(int64_t n, double** buf_out) {
    void* b = NULL;
    Expr* e = ndbuild_open(1, &n, NDT_FLOAT64, &b);
    *buf_out = (double*)b;
    return e;
}

Expr* ndbuild_open_i64(int64_t n, int64_t** buf_out) {
    void* b = NULL;
    Expr* e = ndbuild_open(1, &n, NDT_INT64, &b);
    *buf_out = (int64_t*)b;
    return e;
}

void ndbuild_abandon(Expr* node, size_t written, Expr** args_out) {
    if (!node) return;
    NDType dt = node->data.ndarray.dtype;
    const void* buf = node->data.ndarray.data;
    for (size_t i = 0; i < written; i++)
        args_out[i] = ndarray_buffer_element_to_expr(buf, i, dt);
    expr_free(node);
}

/* --------------------------------------------------------------------------
 * Builtins
 * -------------------------------------------------------------------------- */

/* Strip a trailing DataType -> "..." rule, mirroring NDArray[]'s option
 * handling. Returns the number of leading positional args. */
static size_t pack_take_dtype(Expr* res, bool* have, NDType* dt) {
    size_t argc = res->data.function.arg_count;
    *have = false;
    while (argc > 1) {
        Expr* last = res->data.function.args[argc - 1];
        if (last->type != EXPR_FUNCTION || last->data.function.arg_count != 2) break;
        const Expr* h = last->data.function.head;
        if (h->type != EXPR_SYMBOL ||
            (h->data.symbol.name != SYM_Rule && h->data.symbol.name != SYM_RuleDelayed))
            break;
        const Expr* key = last->data.function.args[0];
        if (key->type != EXPR_SYMBOL || key->data.symbol.name != SYM_DataType) break;
        const Expr* val = last->data.function.args[1];
        if (val->type == EXPR_STRING && !*have)      /* rightmost wins */
            *have = ndt_from_string(val->data.string, dt);
        argc--;
    }
    return argc;
}

static Expr* builtin_tondarray(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;
    bool have = false;
    NDType dt = NDT_FLOAT64;
    size_t argc = pack_take_dtype(res, &have, &dt);
    if (argc != 1) return NULL;

    Expr* arg = res->data.function.args[0];
    /* Already packed, and no conflicting DataType asked for: nothing to do. */
    if (is_packed_list(arg) && (!have || arg->data.ndarray.dtype == dt))
        return expr_copy(arg);
    /* A visible NDArray[...] repacks through its List form, so ToNDArray is a
     * way to say "same values, but as a List". */
    Expr* list = is_ndarray(arg) ? ndarray_to_nested_list(arg) : expr_copy(arg);
    return pack_force(list, have, dt);
}

static Expr* builtin_fromndarray(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    if (!is_ndarray(arg)) return expr_copy(arg);   /* already an ordinary List */
    return ndarray_to_nested_list(arg);
}

/* --------------------------------------------------------------------------
 * The packed-aware head list.
 *
 * Every head named here is a claim that it has been CHECKED to handle a packed
 * list -- because the evaluator will then hand it one instead of materialising
 * (see the gate in evaluate_step). Getting an entry wrong is a silent wrong
 * answer, so the list lives here as one reviewable block rather than as a
 * symtab_set_packed_aware call buried in each module.
 *
 * The ~80 heads with an element-wise kernel (src/ndkernels.c) opt in
 * automatically via symtab_set_ndarray_*_kernel, and are not repeated here.
 *
 * DELIBERATELY ABSENT, and each for a reason:
 *   List, Rule, RuleDelayed, Association, Hold, HoldForm
 *       These enforce the NO-NESTING INVARIANT. Because List is not aware,
 *       {a, packed} materialises its packed element as the enclosing List is
 *       built, so a packed node can never sit inside a plain EXPR_FUNCTION
 *       tree -- which is exactly what makes the gate's O(argc) top-level scan
 *       complete rather than needing to be O(tree) on every unaware head.
 *   ListQ, VectorQ, MatrixQ, ArrayQ, AtomQ, Insert, Delete, Position,
 *   ReplacePart, Append, Prepend, AppendTo, PrependTo, LeafCount, Simplify
 *       Correct by omission: they get a materialised List and answer normally.
 *       Several would be actively wrong if marked aware -- is_listq is a bare
 *       head test that reports False, Append returns unevaluated on an
 *       ndarray, and LeafCount would count a 10^6-element list as one node and
 *       skew Simplify's complexity metric. Marking them aware is an
 *       optimisation with a correctness precondition, not a free win.
 * -------------------------------------------------------------------------- */
static void pack_mark_aware_heads(void) {
    static const char* const AWARE[] = {
        /* Arithmetic and algebra: scan for an ndarray operand themselves
         * (src/plus.c, times.c, power.c, linalg/dot.c). */
        "Plus", "Times", "Power", "Dot",
        /* Subtract and Divide joined on 2026-08-02, found by the gate pass of
         * tools/nd_fastpath_sweep.py. Neither reads an element: builtin_subtract
         * ALWAYS rewrites to Plus[a, Times[-1, b]], and builtin_divide falls
         * through to Times[a, Power[b, -1]] for anything that is not a scalar
         * Real or Rational -- which an NDArray never is. So the fast path was
         * already there, one rewrite away, and the gate was firing at the head
         * that does the rewriting rather than at the head that can use the
         * buffer. Exactly the MinMax / Nest / Fourier defect, and `a - 1` and
         * `a/2` are about as common as an array expression gets: the gate boxed
         * every element, and Listable then threaded one Subtract call PER
         * ELEMENT over the boxed copy.
         *
         * Subtract is int64_ok below (integer minus integer is an integer);
         * Divide deliberately is NOT, because Range[10]/2 is a list of exact
         * Rationals and no buffer holds one. */
        "Subtract", "Divide",
        /* Structural ops over the buffer (src/ndstruct.c). ConjugateTranspose
         * joined 2026-08-01: it delegates to Transpose and Conjugate, both of
         * which have had buffer paths all along, but its rank-1 one-argument
         * form handed a vector to Transpose (which declines one) and came back
         * unevaluated -- so it was EXEMPT rather than aware, and a packed
         * 1000x1000 cost 430 ms against NumPy's 2.20 ms. */
        "Reverse", "Sort", "Flatten", "Transpose", "ConjugateTranspose",
        "Take", "Drop", "Diagonal",
        /* Ordering (src/ndstruct.c) argsorts the buffer and returns the int64
         * permutation of positions -- so, unlike Sort, its result dtype is int64
         * regardless of the input dtype. Still a single-headed (Integer) result,
         * so it stays exact and packable. */
        "Ordering",
        "RotateLeft", "RotateRight", "Partition", "Riffle", "Join",
        "Differences", "First", "Last", "Most", "Rest",
        "PadLeft", "PadRight",
        /* ListConvolve / ListCorrelate (src/convolutions.c). Both engines --
         * the direct O(L*m) loop and the FFT -- already worked on flat double
         * arrays; what they lacked was a way to GET one without boxing every
         * pixel first. A 1024^2 float64 image through a 5x5 kernel spent its
         * time building 2 million Expr nodes, not convolving. */
        "ListConvolve", "ListCorrelate",
        /* Reductions and order statistics (src/ndreduce.c, src/stats.c).
         * MinMax joined on 2026-07-31 from the coverage sweep (experiment 20).
         * It delegates to Min and Max, each of which has had a buffer path all
         * along -- but the head was not here, so the gate materialised 10^6
         * Expr nodes first and both delegates then ran the generic List code on
         * the boxed copy: 307 ms against 426 us for NumPy's
         * [v.min(), v.max()], while Min[v] and Max[v] called directly cost
         * 255 us and 232 us. The same defect as Nest, Fourier and Outer before
         * it, found this time by enumeration rather than by a workload. */
        "Total", "Mean", "Min", "Max", "MinMax", "Median", "Variance",
        "StandardDeviation", "RootMeanSquare", "Quartiles", "Accumulate",
        "MovingAverage", "MovingMedian", "ExponentialMovingAverage",
        /* Moment[v, r] — the r-th RAW (power) moment on the buffer (ndred_moment),
         * (1/n)Sum[x^r]. CentralMoment[v, r] — the r-th moment about the mean
         * (ndred_central_moment), like Variance but /n and power r. NEITHER is on
         * INT64_OK: an integer vector's exact moment is a Rational no machine slot
         * holds, so both degrade to the exact List path. Skewness / Kurtosis are
         * standardized central moments on the buffer (ndred_skewness /
         * ndred_kurtosis); an integer sample's value is a radical, so they degrade
         * the same way and are likewise not INT64_OK. */
        "Moment", "CentralMoment", "Skewness", "Kurtosis",
        /* Covariance / Correlation (src/linalg/ndcorrcov.c). Covariance[v, w] is a
         * threaded centered inner product off the buffer; the matrix forms are a
         * BLAS gram A_c^T B_c. Off this list the gate would materialise both 10^7
         * vectors just to walk them elementwise. NOT int64_ok: an integer sample's
         * covariance is a Rational no float64 slot holds, so an integer buffer
         * degrades to the exact List path, exactly like Variance. */
        "Covariance", "Correlation",
        /* RankedMin[v, n] / RankedMax[v, n] — the n-th order statistic, selected
         * straight off the buffer (ndred_ranked_*): int64 exactly, real by an
         * O(n) quickselect. Off this list the gate would materialise the whole
         * vector just to pick one element, the same defect as MinMax above. */
        "RankedMin", "RankedMax",
        /* Tally is a keyed (irregular) reduction, not a sweep, but the reason it
         * belongs here is the same: ndred_tally hashes the machine words, where
         * materialising boxed one Expr per element just to hash it. */
        "Tally",
        /* Commonest joined it on 2026-08-02: it IS a tally plus a top-n
         * selection, and shares ndred_tally's counting routine. Off this list it
         * paid 1.19 s for 10^7 int64 where Tally of the same buffer cost 22.6 ms
         * -- the entire difference being the 10^7 Expr the gate built first. */
        "Commonest",
        /* The set operations (src/list/setops.c), for the same reason as Tally:
         * a sorted merge over int64 words, where the generic path allocates one
         * Expr per element and sorts through expr_compare. Each guards itself
         * with setop_any_nd and degrades for every shape outside int64. */
        "Union", "Intersection", "Complement",
        /* DeleteDuplicates joined them on 2026-08-01. It was the last head the
         * gate diagnostic still named on a set-op workload: with no NDArray path
         * it materialised the whole buffer just to hash it, 1.31 s over 10^7
         * int64. It now dedups on the machine words -- direct-indexed when the
         * value range is bounded, an inline-key hash set otherwise. */
        "DeleteDuplicates",
        /* Clip came BACK on 2026-07-31. It was pulled because an exact bound
         * makes Clip[reals, {1, 2}] answer with exact Integers where it clipped,
         * which no uniform buffer holds -- but clearing the whole head to fix the
         * bounds cost 356 ms on a 10^6 vector (np.clip: 4.7 ms). ndstruct_clip
         * now gates on the BOUNDS instead: Real bounds stay on the buffer, and an
         * exact bound is taken only when an in-range scan proves nothing is
         * clipped, where the answer is the input. Everything else degrades. */
        "Clip",
        /* Rescale rewrites its call into (x - min)/(max - min) over the WHOLE
         * argument and lets Listable Plus/Times thread it, so it never reads an
         * element and the buffer flows straight into the elementwise kernels.
         * Deliberately not int64_ok: Rescale[Range[10]] is a list of exact
         * Rationals, which no buffer holds. */
        "Rescale",
        /* Chop joined on 2026-08-12. Like Clip's exact bounds, a chopped element
         * is the exact Integer 0, which no uniform Real buffer holds -- so
         * chop_ndarray (src/core.c) scans and keeps the buffer only when NOTHING
         * chops, degrading to the exact List path otherwise. Before it opted in,
         * Chop[largeList] materialised at ~40x an aware head. */
        "Chop",
        /* Note the absentees just below this list -- Precision and Accuracy
         * USED to be here and had to come out. */
        /* Functional heads with packed paths (src/funcprog.c). */
        "Map", "MapIndexed", "MapAt", "MapAll", "Select", "TakeWhile",
        "LengthWhile", "SelectFirst", "Scan", "MapThread", "Fold", "FoldList",
        /* Outer, for nd_outer2 -- the all-pairs buffer path. Same story as
         * Fourier below: the fast path was written first and was unreachable
         * until the head opted in, because the gate had already turned both
         * operands into Lists. builtin_outer materialises again for every case
         * nd_outer2 declines (symbolic f, integer dtypes, three or more
         * tensors), so only the float64 arithmetic shapes see a buffer. */
        "Outer",
        /* The iteration family. These NEVER inspect the value's structure --
         * nest_impl expr_copy's it and hands it to the function, so the next head
         * is what the gate actually has to protect, exactly as for Fold.
         *
         * Leaving them off cost more than any other omission. Nest is THE idiom
         * for iterating an array-valued step, and it was materialising the buffer
         * once per iteration: a 512x512 Jacobi sweep body evaluates in 1.10 ms
         * standalone, but Nest[f, u0, 10] took 2.19 s -- 200x -- against 0.017 s
         * for the identical Do[u = f[u], {10}] loop, and answered with an
         * unpacked List. Fold was aware; Nest was not, for no reason. */
        "Nest", "NestList", "NestWhile", "NestWhileList",
        "FixedPoint", "FixedPointList",
        /* Metadata: read rank/dims/dtype straight off the node, so
         * materialising for them would be pure waste. */
        "Length", "Dimensions", "Depth", "Part", "Head", "ByteCount",
        "NDArrayQ", "DataType", "Normal", "ToNDArray", "FromNDArray",
        "ToPackedArray", "FromPackedArray",
        /* Spectral transforms. fourier.c has had a complete NDArray fast path
         * all along (`data->type == EXPR_NDARRAY` -> machine_path_ndarray,
         * FFTW under USE_FFTW), but the head was never marked aware, so the
         * gate materialised the buffer into one Expr per element before the
         * builtin could see it -- and the fast path was unreachable from
         * ordinary list code. Measured on 2^20 reals: Fourier[packedList]
         * 986 ms against Fourier[NDArray[...]] 62.8 ms on identical data, the
         * two differing only in `present_as`. */
        "Fourier", "InverseFourier", "FourierDCT", "FourierDST",
        /* Heads with an NDArray path that tools/check_packed_aware.py surfaced
         * and a differential sweep then proved safe: each answers identically
         * for a buffer and for the same value as a plain List. The six that did
         * NOT are listed in that script's EXEMPT table with the observed
         * difference. ReverseSort is aware because the sweep found and got fixed
         * a live wrong answer in it (see reverse_top_level in src/sort.c). */
        "AtomQ", "MatrixQ", "VectorQ", "N", "TeXForm", "Extract",
        "QuotientRemainder", "AllTrue", "AnyTrue", "NoneTrue", "ReverseSort",
        /* Pass-through: the value is moved, not inspected, so a packed list
         * survives assignment and control flow. The one place assignment reads
         * the RHS *structure* is list destructuring ({a, b} = {...}); leaving
         * Set aware would strand a packed RHS there, so apply_assignment
         * (src/eval.c) normalises it to a List of its (still-packed) top-level
         * slices before destructuring. Without that, {xc, yc} = {Range[m],
         * Range[n]} silently bound nothing once the RHS packed to a rank-2
         * buffer. */
        "Set", "SetDelayed", "CompoundExpression", "If", "Which", "Module",
        "Block", "With", "Return", "Print", "Echo",
        /* ------------------------------------------------------------------
         * Linear algebra. Every one of these opens with
         *     if (linalg_call_has_ndarray(res)) return ndla_<op>(res);
         * or the same test routing to linalg_delist_and_reeval, so each already
         * handles a buffer -- that test is tag-based and so cannot tell the
         * visible NDArray[...] from a packed list.
         *
         * Leaving them off this list did not merely cost the fast path, it cost
         * CORRECTNESS: the gate materialised the buffer, linalg_call_has_ndarray
         * then read false, and a machine-real matrix went down the exact
         * fraction-free (Bareiss) path, which calls polynomial GCD per pivot and
         * recurses. LinearSolve / Det / Inverse of a 90x90 Real matrix
         * STACK-OVERFLOWED. Reaching LAPACK is what keeps a dense solve linear
         * in memory as well as fast.
         *
         * Marked aware but deliberately NOT int64_ok (see the list below): the
         * exact determinant of an integer matrix is an Integer and its inverse
         * is a matrix of Rationals, neither of which a float64 kernel can
         * express, so an integer buffer still materialises and takes the exact
         * path. That split is the whole reason the two flags are separate. */
        "Det", "Inverse", "PseudoInverse", "LinearSolve", "LeastSquares",
        "RowReduce", "NullSpace", "MatrixRank", "MatrixPower",
        "LUDecomposition", "QRDecomposition", "SingularValueDecomposition",
        "Eigenvalues", "Eigenvectors", "Tr", "Norm", "Normalize", "Cross",
        "DiagonalMatrix", "HankelMatrix", "ToeplitzMatrix", "VandermondeMatrix",
        "PositiveDefiniteMatrixQ", "NegativeDefiniteMatrixQ",
        "LatticeReduce", "FindIntegerNullVector",
        /* Fit reads its data straight off the buffer (fit_read_numeric in
         * src/fit.c) and builds a packed float64 design matrix column-wise, so
         * the gate must not materialise a packed data matrix into 2*npts boxed
         * Exprs first. A fit is inherently a machine-real computation, so an
         * int64 data buffer gives the same float coefficients as materialising
         * would (see INT64_OK); WorkingPrecision -> Infinity bails from the fast
         * path and re-delists for the exact rational solve.  DesignMatrix accepts
         * the NDArray surface too (fit_normalize_data delists it) and so is aware
         * for the same reason -- it keeps its exact List output, it just no
         * longer errors on a packed/NDArray argument. */
        "Fit", "DesignMatrix",
        /* Interpolation reads a packed/NDArray data table straight through
         * ndarray_to_nested_list (src/interp.c) instead of erroring, so the gate
         * must not pre-materialise it.  Separately, applying the resulting
         * InterpolatingFunction to a packed array of query points is kept off
         * this gate by the interp_head branch in eval.c (the head is the object,
         * not a symbol), where the vectorised evaluator reads the point buffer
         * directly. */
        "Interpolation",
        /* InterpolatingPolynomial reads a packed/NDArray data table straight
         * through ndarray_to_nested_list (src/interp.c) as well, so the gate
         * must not pre-materialise it. */
        "InterpolatingPolynomial",
        /* ------------------------------------------------------------------
         * The integer domain (src/ndinteger.c). Every one of these answered a
         * packed List with one Expr per element and a visible NDArray with the
         * unevaluated call, because neither surface had a path at all -- the
         * class tools/check_packed_aware.py is blind to by construction, since
         * there was no dispatch site for it to read.
         *
         * Prime, PowerMod and IntegerDigits are the ones that needed care: each
         * has a head-level path rather than an element kernel, and each tests
         * for EXPR_INTEGER below it. Marking them aware stops the gate
         * materialising, so their buffer now arrives where an integer test
         * fails -- they degrade through ndarray_delist_and_reeval explicitly
         * instead of falling out unevaluated. */
        "GCD", "LCM", "DivisorSigma", "EulerPhi", "MoebiusMu",
        "IntegerLength", "PowerMod", "Prime", "IntegerDigits",
        /* The elementwise sign predicates. Now buffer on BOTH sides (§13 gap C.1
         * closed): the input is read straight off the buffer and the answer is a
         * packed NDT_BOOL array of True/False, never 10^6 Real comparisons nor
         * 10^6 True/False Expr nodes (see ndint_sign_predicate). */
        "Positive", "Negative", "NonNegative", "NonPositive",
        /* Boole maps a bool buffer straight to an int64 one (True -> 1); its
         * buffer path answers with a single head (Integer), matching Boole of the
         * List element for element. */
        "Boole",
        /* ------------------------------------------------------------------
         * The structural batch (src/ndstruct.c, src/random.c, src/assoc.c).
         *
         * Append and Prepend were named in the "correct by omission" note just
         * above this list -- the gate materialised, the generic walk answered,
         * and adding one element to a 10^6-element vector cost 10^6 Expr
         * allocations. They now have buffer paths, so the note no longer
         * applies to them and they have moved here.
         *
         * Catenate reads a rank-2 buffer as a reshape; Ratios divides in
         * place; TakeLargest/TakeSmallest keep a bounded heap; RandomSample and
         * RandomChoice gather from the buffer through the SAME draw sequence
         * the List path uses, so a seeded run agrees across surfaces; Counts
         * relabels Tally's machine-word count. */
        "Append", "Prepend", "Catenate", "Ratios",
        "TakeLargest", "TakeSmallest",
        "RandomSample", "RandomChoice", "Counts", "Inner",
        /* The hypergeometric family. The three convenience heads rewrite to
         * HypergeometricPFQ and pass the argument straight through, so ALL FOUR
         * have to be aware: the gate fires at whichever head is holding the
         * buffer, and marking only PFQ would have materialised at the rewrite,
         * one step before the path that can use it. Deliberately not int64_ok:
         * the series runs in double, and an exact-integer argument is evaluated
         * in its own right by the List path. */
        "HypergeometricPFQ", "Hypergeometric0F1",
        "Hypergeometric1F1", "Hypergeometric2F1",
        /* ArrayPlot reads its data table straight off the buffer via
         * na_load_matrix (src/graphics/arrayplot.c, the same loader Fit and
         * DesignMatrix use above) when given an EXPR_NDARRAY, so the gate
         * must not pre-materialise a packed matrix into rows*cols boxed
         * Exprs first just to render a heatmap. A nested-List argument takes
         * a separate per-cell path (arrayplot_load) that is unaffected
         * either way. */
        "ArrayPlot",
    };
    for (size_t i = 0; i < sizeof(AWARE) / sizeof(AWARE[0]); i++)
        symtab_set_packed_aware(AWARE[i]);

    /* ------------------------------------------------------------------
     * Heads whose buffer fast path answers with different ELEMENT HEADS than
     * the ordinary list does, and which therefore must NOT keep a buffer.
     *
     * These were found by sweeping every aware head and every registered
     * kernel packed-against-plain -- not by reading the code, which is why the
     * list is this specific. Each was a live wrong answer:
     *
     *   Floor, Ceiling, Round, IntegerPart, Sign, Im
     *       Real in, INTEGER out. The kernel is "real closed" so it keeps the
     *       float64 dtype and writes 1.0 where the list gives the exact 1:
     *       Floor[{0.25, 0.5, 1.0}] is {0, 0, 1}, not {0., 0., 1.}.
     *   Clip
     *       Clip[x] clamps to [-1, 1] using the EXACT bounds, so a clipped
     *       element comes back as the Integer 1 -- Clip[{1., 2., 3.}] is
     *       {1., 1, 1}. The buffer makes all three Real.
     *   Precision, Accuracy
     *       Listable in Mathilda, so the list answers with one value PER
     *       ELEMENT while the buffer path answers with a single scalar.
     *
     * The kernel heads acquire packed_aware implicitly from
     * symtab_set_ndarray_*_kernel, so they need clearing rather than omitting.
     * All of them keep their fast path for a VISIBLE NDArray[...], where
     * machine-buffer coercion is the documented behaviour of naming that head;
     * only the invisible surface materialises. Giving them exact paths (an
     * int64 output buffer for Floor and friends) would let them back on this
     * list -- an optimisation, not a correctness question. */
    static const char* const NOT_AWARE[] = {
        /* Floor, Ceiling, Round, IntegerPart and Sign came OFF this list on
         * 2026-07-30: they now have narrowing kernels (real in, NDT_INT64 out --
         * see NDUnaryKernel.to_int in ndarray.h), so the buffer path produces the
         * exact Integers the List does instead of 1.0 where the List gives 1.
         * That was the "optimisation, not a correctness question" the note below
         * anticipated. `Im` joined them on 2026-08-12: Im of any real (or exact
         * integer) is the Integer 0, so its narrowing pair (NDKU_ImInt in
         * ndkernels.c) writes an NDT_INT64 zero buffer that matches the List path
         * -- previously it materialised, at ~22x the cost of an aware head.
         *
         * Clip left this list on 2026-07-31 for the same reason -- the exactness
         * problem was real but belonged to the BOUNDS, not to the head. See the
         * note beside it in AWARE. */
        "Precision", "Accuracy",
        /* Quotient joined this list on 2026-07-30 and came OFF it on
         * 2026-08-01, the same way Floor and friends did. It was here because
         * real in / exact Integer out could not be expressed -- and it was
         * OBSERVABLE from ordinary code, because packing keys on length:
         *     Quotient[Range[1., 300.], 2]  ->  {0., 1., 1., 2.}   (packed)
         *     Quotient[Range[1., 200.], 2]  ->  {0, 1, 1, 2}       (not packed)
         * the same expression answering with different element heads either
         * side of the threshold, where Mathematica gives the exact Integers
         * throughout. NDBinaryKernel now carries the same to_int trio as
         * NDUnaryKernel (ndarray.h) and NDKB_Quotient supplies both arms, so
         * the buffer path produces those exact Integers itself instead of
         * having to be kept away from the data. */
    };
    for (size_t i = 0; i < sizeof(NOT_AWARE) / sizeof(NOT_AWARE[0]); i++)
        symtab_clear_packed_aware(NOT_AWARE[i]);

    /* Exact on an int64 buffer as well -- see SymbolDef.packed_int64_ok. This
     * list is deliberately short and grows only as sites are given exact int64
     * paths. Everything absent from it materialises an integer argument, which
     * costs the computation its fast path but keeps the answer identical to the
     * ordinary list's, heads included.
     *
     * Verified for each entry below: it touches the buffer only through
     * ndt_get_i or memcpy, or does not read elements at all.
     *
     * Note what is NOT here and why, since these are the ones that bite:
     *   the kernel heads       a real_closed kernel keeps the input dtype and
     *                          truncates -- and Sin[{1,2,3}] is SYMBOLIC in the
     *                          interpreter, so the right answer is to degrade,
     *                          not to compute in float
     *   Variance/Std/RMS       the exact answer is a Rational or a radical
     *                          (Variance[Range[10]] is 55/6), so there is
     *                          nothing a buffer could hold
     *   Quartiles              interpolates between order statistics; the exact
     *                          answer is {151/2, 301/2, 451/2}
     *   MovingAverage/Median/  same: exact integer inputs give Rational windows
     *   ExponentialMovingAverage
     *   Clip                   the clamp runs through doubles, and the bounds
     *                          may themselves be exact
     *   Precision/Accuracy     would answer MachinePrecision where an exact
     *                          integer list is Infinity
     *   Flatten/Take/Drop/     safe (they memcpy) but not yet re-verified with
     *   Reverse/Join/...       the int64 write paths; they materialise for now
     * Run the suite with MATHILDA_PACK_DIAG=1 to see the remaining lossy sites;
     * `nd_int64_lossy_hit` in ndarray.c is a breakpoint target. */
    static const char* const INT64_OK[] = {
        /* Read rank/dims/dtype only, never an element. */
        "Length", "Dimensions", "Depth", "Head", "ByteCount",
        "NDArrayQ", "DataType",
        /* N reads every element, but its answer on an int64 buffer IS a real
         * array -- the machine reals the interpreter's N gives element by
         * element -- so nothing truncates and no element's head changes. Its
         * EXPR_NDARRAY case in numericalize_rec (src/numeric.c) widens the
         * whole buffer to float64 in one pass, inheriting presentation, so
         * N[Range[10^6]] stays a packed array instead of materialising to a
         * list of boxed reals. The (double) cast rounds int64 past 2^53 exactly
         * as N does to the same integer scalar, so the surfaces agree. */
        "N",
        /* Materialise or index through the exact accessors. Extract joined
         * 2026-08-01 with its buffer gather: it is Part's own selector, so the
         * justification is Part's -- ndarray_element_to_expr yields an exact
         * Integer from an int64 buffer, so no element's head changes. */
        "Part", "Extract", "Normal", "ToNDArray", "FromNDArray",
        "ToPackedArray", "FromPackedArray",
        /* Fit reads an int64 data buffer through ndt_get exactly as it would the
         * materialised integers, and a least-squares fit is a machine-real
         * computation either way, so the float coefficients are identical -- the
         * flag only spares an int64 data matrix from materialising. */
        "Fit",
        /* InterpolatingPolynomial reads int64 data exactly and, on exact data,
         * returns an exact polynomial -- so an int64 buffer needs no
         * materialisation and the result is identical. */
        "InterpolatingPolynomial",
        /* ArrayPlot reads an int64 data buffer through na_load_matrix exactly
         * as it would the materialised integers (arrayplot_load's numeric arm
         * is expr_to_real_double either way), and cell colouring is a machine-
         * real computation regardless of the input's exactness -- ArrayPlot
         * never produces an exact result, so nothing is lost that the
         * materialised path would have kept. The flag only spares an int64
         * array from materialising for the common RandomInteger[...]/
         * Table[Mod[...], ...] usage shown in its own docstring. */
        "ArrayPlot",
        /* Move the value without inspecting it. (List destructuring is the one
         * exception -- apply_assignment reads an int64 slice back as an exact
         * Integer via ndarray_part, so {a, b, c, d} = Range[4] binds exact
         * Integers, not coerced Reals; see src/eval.c.) */
        "Set", "SetDelayed", "CompoundExpression", "If", "Which", "Module",
        "Block", "With", "Return", "Print", "Echo",
        /* Exact int64 arithmetic (src/ndarray.c): the accumulate uses ci_*_i64
         * and abandons the whole result on the first overflow, so
         * Total[Range[10^6]^3] still reaches GMP and answers exactly. A Real
         * scalar widens the buffer to float64; a Rational, BigInt, MPFR or
         * Complex scalar bails, because Range[10]/2 is a list of Rationals. */
        "Plus", "Times", "Power", "Dot",
        /* Subtract rewrites to Plus[a, Times[-1, b]] and inherits both arms'
         * exactness, so an integer buffer stays one. Divide is NOT here on
         * purpose: it rewrites to Times[a, Power[b, -1]], and Range[10]/2 is a
         * list of exact Rationals -- the integer buffer has to materialise for
         * the List path to give them. */
        "Subtract",
        /* Exact int64 reductions (src/ndreduce.c). Total/Max/Min/Accumulate
         * answer with Integers; Mean and Median build the exact reduced
         * Rational (Mean[Range[10]] is 11/2, Median[Range[300]] is 301/2). */
        "Total", "Mean", "Median", "Max", "Min", "MinMax", "Accumulate",
        /* RankedMin/RankedMax SELECT an element (like Max/Min), so the r-th order
         * statistic of an int64 vector is an Integer taken straight off the
         * buffer -- exact past 2^53, no rounded Real. */
        "RankedMin", "RankedMax",
        /* Tally keys on the raw int64 word, so two integers past 2^53 stay
         * distinct -- the failure a float64 gather would have introduced -- and
         * each distinct value is rebuilt with expr_new_integer, so no element's
         * head changes. Integer data is also the usual case: an int64 buffer
         * that could not keep its dtype here would leave the common Tally
         * exactly where it started. */
        "Tally",
        /* Commonest keys on the same raw word, and every element it returns was
         * COPIED out of the int64 buffer -- there is no arithmetic to be exact
         * about, and no element's head can change. */
        "Commonest",
        /* Exact int64 structure (src/ndstruct.c): Sort orders the buffer in
         * int64 -- through doubles, two integers past 2^53 would compare equal
         * and be reordered -- and Transpose moves elements by memcpy. Ordering
         * argsorts in int64 for the same reason, so the permutation is correct
         * past 2^53 too. */
        "Sort", "Transpose", "Ordering",
        /* Pure element MOVES: whole rows by memcpy, never through ndt_get, so
         * exactness is not in question -- there is no arithmetic to be exact
         * about. Verified individually in src/ndstruct.c; Reverse also keeps
         * Sort[Reverse[intList]] on the buffer end to end, which was the last
         * measured place where packing was slower than a plain List. */
        "Reverse", "Flatten", "Take", "Drop", "Diagonal",
        /* The leading-axis slices, added 2026-07-31 with the fourth sweep's
         * ndstruct_head_tail. A rank-1 First/Last is one element through
         * ndarray_element_to_expr, which yields an exact Integer from an int64
         * buffer; rank >= 2, and Most/Rest at any rank, are whole-row memcpys.
         * There is no arithmetic to be exact about, and leaving them off cost
         * the sequence-alignment kernel its whole inner row: Most[prev] and
         * Rest[prev] on a 2000-element int64 buffer materialised twice per row,
         * 2000 rows, and then poisoned the MapThread that consumed them. */
        "First", "Last", "Most", "Rest",
        /* Abs, with the integer arm in ndkernels.c (NDKU_AbsInt). */
        "Abs",
        /* Also pure moves, added with the native buffer paths in ndstruct.c:
         * Rotate/Join/Partition/Riffle/Pad move whole elements by memcpy and
         * refuse any fill or separator that is not an exact Integer at int64
         * (nd_fill_value), so no element's head can change. Differences is the
         * one that computes: it subtracts in int64 via ci_sub_i64 and abandons
         * the whole result on the first overflow, so the List path re-runs it and
         * GMP answers exactly -- never a wrapped difference. */
        "RotateLeft", "RotateRight", "Join", "Partition", "Riffle",
        "PadLeft", "PadRight", "Differences",
        /* The iteration family reads no element at all -- it copies the value and
         * applies a function, so exactness is entirely the function's business. */
        "Nest", "NestList", "NestWhile", "NestWhileList",
        "FixedPoint", "FixedPointList",
        /* The scan and the elementwise thread, added 2026-07-31 with the integer
         * arms in ndreduce.c and funcprog.c. Each was verified individually:
         *
         *   Fold / FoldList  ndred_scan's want_exact branch runs the whole scan
         *                    in int64 and requires an exact SEED, so a Real seed
         *                    over an integer buffer declines rather than
         *                    silently returning Reals where Max picked the seed.
         *                    Plus and Times go through ci_add_i64 / ci_mul_i64
         *                    and abandon the array on overflow. The compiled-VM
         *                    path beside it already refuses an int64 buffer
         *                    outright (elems_inexact is false), so there is no
         *                    second route that could compute in double.
         *   MapThread       nd_mapthread2's int64 arm is Min/Max by comparison
         *                    (exact by construction -- the answer is one of the
         *                    inputs) and Plus/Times through the same checked
         *                    helpers.
         *
         * Needleman-Wunsch is the workload that found these: its inner row is
         * MapThread[Max, {diag, up}] followed by FoldList[Max, ...], both int64,
         * and materialising for each cost 853 ms and 302 ms per 10^6 against
         * 12 ms and 1.0 ms for the identical Real pair. */
        "Fold", "FoldList", "MapThread",
        /* The set operations: int64 IS their fast domain (setop_packed handles
         * nothing else), and each element is rebuilt with expr_new_integer, so
         * no head changes. DeleteDuplicates is int64-only for the same reason
         * spelled out at the top of setops.c: 0. and -0. compare equal but print
         * differently, so over Reals WHICH of two equal elements survives would
         * differ between the buffer and the List path. */
        "Union", "Intersection", "Complement", "DeleteDuplicates",
        /* The narrowing kernels. Exact on an int64 buffer because that arm is
         * written in int64 throughout: Floor/Ceiling/Round/IntegerPart are the
         * identity on an integer, Sign and UnitStep are trivial comparisons, and
         * a real input that will not fit an int64 abandons the whole array so
         * the List path answers with a bignum. */
        "Floor", "Ceiling", "Round", "IntegerPart", "Sign", "UnitStep",
        /* Mod and Quotient, once ndkernels.c gave them exact int64 arms
         * (NDKB_Mod.to_int_i, NDKB_Quotient.to_int_i / .to_int_r). Until then
         * their kernel computed m - n*floor(m/n) in double and would have
         * written {0., 1.} where the List gives the exact {0, 1}, so they had
         * to stay off this list -- and the cost of that landed on their
         * CONSUMERS, not on them: Mod[Range[10^6] 7919, 1000] answered with
         * 10^6 boxed Integers and the Union / Tally / DeleteDuplicates that
         * read it ran 182x / 45.8x / 34.0x slower than the same data packed.
         *
         * An int64 buffer with a REAL scalar (Mod[{1,2,3}, 2.] = {1., 0., 1.})
         * has no int64 answer; ndarray_map_binary declines it and the List path
         * answers with the Real elements, exactly as before. */
        "Mod", "Quotient",
        /* The number-theoretic heads, once src/ndinteger.c gave them int64
         * kernels. int64 IS their domain -- EulerPhi[3.] is not EulerPhi[3] --
         * so unlike every other entry on this list the int64 claim is not an
         * extension of a float64 one, it is the only claim there is. Each
         * kernel is written in int64 end to end through the ci_*_i64 helpers
         * and abandons the whole array on overflow, on a factor past the
         * trial-division ceiling, or on an argument the scalar builtin leaves
         * unevaluated (MoebiusMu[0], DivisorSigma[k, 0]) -- so the List path
         * and GMP still answer every case the buffer cannot. */
        "GCD", "LCM", "DivisorSigma", "EulerPhi", "MoebiusMu",
        "IntegerLength", "PowerMod", "Prime", "IntegerDigits",
        /* The sign predicates answer with True/False, which no dtype holds --
         * but they READ an int64 buffer exactly, by comparison, and comparison
         * is the one operation on integers that needs no arithmetic to be
         * exact about. */
        "Positive", "Negative", "NonNegative", "NonPositive",
        /* Pure element MOVES over an int64 buffer, exact for the same reason
         * Join and Take are: whole elements travel by memcpy and there is no
         * arithmetic to be exact about. Append/Prepend additionally refuse any
         * appended value that is not an exact Integer at int64 (nd_fill_value),
         * so Append[Range[10], 1.] declines rather than coercing the 1. into an
         * integer slot. Counts keys on the raw word exactly as Tally does.
         *
         * Deliberate absentees: RATIOS, whose int64 answer is a list of exact
         * Rationals -- Ratios[{1,2,3}] is {2, 3/2} -- and TAKELARGEST /
         * TAKESMALLEST, which order through `double` and would compare two
         * integers past 2^53 equal. Both decline an int64 buffer in the kernel
         * as well; leaving them off here is what routes the whole call to the
         * List path instead of round-tripping through a decline. */
        "Append", "Prepend", "Catenate",
        "RandomSample", "RandomChoice", "Counts",
    };
    for (size_t i = 0; i < sizeof(INT64_OK) / sizeof(INT64_OK[0]); i++)
        symtab_set_packed_int64_ok(INT64_OK[i]);

    /* ------------------------------------------------------------------
     * Heads whose buffer path also handles operands of DIFFERENT RANK --
     * see SymbolDef.packed_broadcast_ok. Only these two: ndarray_elementwise
     * grew a leading-axis broadcast pre-pass, ndarray_elementwise_power did
     * not and still requires same_shape, so Power stays off this list and a
     * mixed-rank Power keeps threading through the interpreter.
     * ------------------------------------------------------------------ */
    static const char* const BROADCAST_OK[] = { "Plus", "Times" };
    for (size_t i = 0; i < sizeof(BROADCAST_OK) / sizeof(BROADCAST_OK[0]); i++)
        symtab_set_packed_broadcast_ok(BROADCAST_OK[i]);
}

/* --------------------------------------------------------------------------
 * MATHILDA_PACK_DIAG=gate -- which heads are making the gate materialise?
 *
 * The transparency gate's rule is "materialise for any head that has not opted
 * in", which makes a missing opt-in fail SAFE and therefore SILENT: right
 * answer, slow path, nothing to see. Four separate subsystems sat like that at
 * 30x-658x, and none was found by reading code.
 *
 * tools/check_packed_aware.py catches the case where a head HAS a fast path and
 * has not opted in. This catches the other half, which no static check can: a
 * head with no fast path at all, that a real workload keeps handing buffers to.
 * `Fourier` is the standing example -- 40x-80x behind Mathematica, and it does
 * not appear in the static audit because there is no dispatch site to find.
 *
 * Run a corpus under MATHILDA_PACK_DIAG=gate and the exit summary is a work
 * list, ordered by how many elements each head actually cost.
 * -------------------------------------------------------------------------- */

#define PACK_DIAG_MAX_HEADS 256

typedef struct {
    const char* name;      /* interned, so a pointer compare identifies it */
    uint64_t    calls;
    uint64_t    elements;
} PackDiagEntry;

static PackDiagEntry g_diag[PACK_DIAG_MAX_HEADS];
static size_t        g_diag_n = 0;
static int           g_diag_on = -1;      /* -1 = not yet read from the env */
static bool          g_diag_atexit = false;

static int pack_diag_cmp(const void* a, const void* b) {
    const PackDiagEntry* x = (const PackDiagEntry*)a;
    const PackDiagEntry* y = (const PackDiagEntry*)b;
    if (x->elements != y->elements) return (x->elements < y->elements) ? 1 : -1;
    return 0;
}

static void pack_diag_dump(void) {
    if (!g_diag_n) {
        fprintf(stderr, "\n[pack] gate: no packed argument was ever materialised.\n");
        return;
    }
    qsort(g_diag, g_diag_n, sizeof(g_diag[0]), pack_diag_cmp);
    fprintf(stderr,
            "\n[pack] transparency gate materialised a packed argument for these "
            "heads\n       (ordered by elements -- each one is a fast path not "
            "taken):\n\n");
    fprintf(stderr, "       %-28s %10s %14s\n", "head", "calls", "elements");
    for (size_t i = 0; i < g_diag_n; i++)
        fprintf(stderr, "       %-28s %10llu %14llu\n", g_diag[i].name,
                (unsigned long long)g_diag[i].calls,
                (unsigned long long)g_diag[i].elements);
    fprintf(stderr,
            "\n       A head here either needs an NDArray fast path, or needs "
            "adding to\n       AWARE / INT64_OK in src/pack.c once it has one. "
            "See\n       docs/design/packed_arrays.md and "
            "plans/HPC_IMPROVEMENT_PLAN.md.\n");
}

void pack_gate_report(const Expr* head, const Expr* arg) {
    if (g_diag_on < 0) {
        const char* v = getenv("MATHILDA_PACK_DIAG");
        g_diag_on = (v && (strstr(v, "gate") || strstr(v, "all"))) ? 1 : 0;
    }
    if (!g_diag_on || !arg || !is_ndarray(arg)) return;

    const char* name = "(non-symbol head)";
    if (head && head->type == EXPR_SYMBOL) name = head->data.symbol.name;

    size_t n = 1;
    for (int i = 0; i < arg->data.ndarray.rank; i++)
        n *= (size_t)arg->data.ndarray.dims[i];

    for (size_t i = 0; i < g_diag_n; i++) {
        if (g_diag[i].name == name || strcmp(g_diag[i].name, name) == 0) {
            g_diag[i].calls++;
            g_diag[i].elements += (uint64_t)n;
            return;
        }
    }
    if (g_diag_n < PACK_DIAG_MAX_HEADS) {
        g_diag[g_diag_n].name = name;
        g_diag[g_diag_n].calls = 1;
        g_diag[g_diag_n].elements = (uint64_t)n;
        g_diag_n++;
    }
    if (!g_diag_atexit) { g_diag_atexit = true; atexit(pack_diag_dump); }
}

void pack_init(void) {
    pack_mark_aware_heads();

    symtab_add_builtin("ToNDArray", builtin_tondarray);
    set_attributes("ToNDArray", ATTR_PROTECTED);
    symtab_set_docstring("ToNDArray",
        "ToNDArray[list] returns list stored as a dense machine-precision "
        "buffer. The result is still a List -- same Head, same printed form, "
        "same elements -- but NDArrayQ gives True for it. "
        "ToNDArray[list, DataType -> \"float64\"] forces the element type. "
        "Returns list unchanged when it is not rectangular, is empty, or holds "
        "anything other than uniformly Integer or uniformly Real machine "
        "values. Unlike automatic packing it ignores the size threshold.");

    /* Mathematica's spellings for the same two operations, so code written
     * against Developer`ToPackedArray / Developer`FromPackedArray reads across.
     * Registered as the SAME builtin rather than as a rule that rewrites to
     * ToNDArray: an alias installed as a DownValue would show up in traces, cost
     * an evaluation pass, and could be shadowed by a user definition on the
     * target. These are genuinely the same function under two names. */
    symtab_add_builtin("ToPackedArray", builtin_tondarray);
    set_attributes("ToPackedArray", ATTR_PROTECTED);
    symtab_set_docstring("ToPackedArray",
        "ToPackedArray[list] is ToNDArray[list]: it returns list stored as a "
        "dense machine-precision buffer. The result is still a List -- same "
        "Head, same printed form, same elements -- but NDArrayQ gives True for "
        "it. Provided under Mathematica's name for the same operation; see "
        "ToNDArray.");

    symtab_add_builtin("FromNDArray", builtin_fromndarray);
    set_attributes("FromNDArray", ATTR_PROTECTED);
    symtab_set_docstring("FromNDArray",
        "FromNDArray[expr] returns expr with any dense buffer storage undone: "
        "a packed List becomes an ordinary List of separate elements, and an "
        "NDArray[...] becomes the nested List of its entries. Anything else is "
        "returned unchanged. Inverse of ToNDArray.");

    symtab_add_builtin("FromPackedArray", builtin_fromndarray);
    set_attributes("FromPackedArray", ATTR_PROTECTED);
    symtab_set_docstring("FromPackedArray",
        "FromPackedArray[expr] is FromNDArray[expr]: it returns expr with any "
        "dense buffer storage undone, so a packed List becomes an ordinary List "
        "of separate elements and an NDArray[...] becomes the nested List of its "
        "entries. Provided under Mathematica's name for the same operation; see "
        "FromNDArray.");
}
