/* NDArray structural fast paths — see ndstruct.h. Each op copies from the input
 * buffer into a fresh output buffer (inputs are immutable during evaluation)
 * and returns a new EXPR_NDARRAY; anything outside the fast domain degrades to
 * ndarray_delist_and_reeval so the result matches the List path exactly. */

#include "ndstruct.h"
#include "pack.h"    /* pack_repack_like — a packed source re-sniffs its dtype */
#include "checked_int.h"  /* ci_sub_i64 — Differences abandons on int64 overflow */
#include "ndarray.h"
#include "ndarray_internal.h"
#include "take_drop.h"  /* get_seq_spec_indices — Ordering[a, seq] slices the permutation */
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>   /* HUGE_VAL -- the infinite Clip bound and Ramp's upper limit */

/* Product of dims[lo..hi). */
static size_t nd_prod(const int64_t* dims, int lo, int hi) {
    size_t p = 1;
    for (int i = lo; i < hi; i++) p *= (size_t)dims[i];
    return p;
}

/* An Integer/Real Expr as a double (Clip bounds); false for anything else.
 *
 * Infinity and -Infinity are accepted and map to +/-HUGE_VAL. A one-sided Clip
 * is how the positive part is spelled -- Clip[x, {0., Infinity}] -- and without
 * this the bound simply failed to parse, so `Clip` returned UNEVALUATED on every
 * ReLU in the system. `*infinite` reports which bounds were infinite, because an
 * infinite bound is never ATTAINED and so cannot put its own head into the
 * answer: it is exempt from the exactness gate below. */
static bool nd_real_value(const Expr* e, double* out, bool* infinite) {
    if (infinite) *infinite = false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real; return true; }
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Infinity) {
        *out = HUGE_VAL; if (infinite) *infinite = true; return true;
    }
    /* -Infinity is Times[-1, Infinity] after canonicalisation. */
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Times &&
        e->data.function.args[0]->type == EXPR_INTEGER &&
        e->data.function.args[0]->data.integer == -1 &&
        e->data.function.args[1]->type == EXPR_SYMBOL &&
        e->data.function.args[1]->data.symbol.name == SYM_Infinity) {
        *out = -HUGE_VAL; if (infinite) *infinite = true; return true;
    }
    return false;
}

bool ndstruct_call_has_ndarray(const Expr* res) {
    return res && res->type == EXPR_FUNCTION &&
           res->data.function.arg_count >= 1 &&
           is_ndarray(res->data.function.args[0]);
}

/* ------------------------------------------------------------------- Sort */

/* Sort[a]: real rank-1 only — sort the values ascending and rebuild an NDArray
 * of the same dtype (float32 values round-trip through double exactly). A custom
 * ordering function (2-arg), a complex dtype (canonical complex order differs
 * from numeric order), or rank >= 2 (Sort orders whole rows by canonical Expr
 * order) all degrade to the List path. */
Expr* ndstruct_sort(Expr* res) {
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    NDType dt = a->data.ndarray.dtype;
    if (a->data.ndarray.rank != 1 || ndt_is_complex(dt))
        return ndarray_delist_and_reeval(res);

    size_t n = (size_t)a->data.ndarray.dims[0];
    if (dt == NDT_INT64) {
        /* Sort an integer buffer in place, in int64: the double gather is exact
         * only to 2^53, so two large integers would compare equal and the sort
         * would silently reorder them -- and every element would come back
         * rounded. */
        int64_t* o = malloc(sizeof(int64_t) * (n ? n : 1));
        if (!o) return ndarray_delist_and_reeval(res);
        memcpy(o, a->data.ndarray.data, sizeof(int64_t) * n);
        nd_sort_i64_asc(o, n);
        int64_t idims[1] = { (int64_t)n };
        return expr_new_ndarray_like(a, 1, idims, o, dt);
    }
    double* s = malloc(sizeof(double) * n);
    if (!s) return ndarray_delist_and_reeval(res);
    nd_gather_real(a->data.ndarray.data, dt, 0, 1, n, s);
    nd_sort_ascending(s, n);

    void* out = malloc(ndt_elem_size(dt) * n);
    if (!out) { free(s); return ndarray_delist_and_reeval(res); }
    for (size_t i = 0; i < n; i++) ndt_set(out, i, dt, s[i], 0.0);
    free(s);
    int64_t dims[1] = { (int64_t)n };
    return expr_new_ndarray_like(a, 1, dims, out, dt);
}

/* --------------------------------------------------------------- Ordering */

/* Ordering[a] / Ordering[a, seq]: the argsort permutation of 1-based positions,
 * straight off the buffer. The result is ALWAYS an NDT_INT64 array regardless of
 * the input dtype (a permutation is integer), which is why it needs its own entry
 * point rather than riding ndstruct_sort's same-dtype rebuild.
 *
 * Real rank-1 only: rank >= 2 (Ordering compares whole rows by canonical Expr
 * order), a complex dtype (canonical complex order differs from numeric order),
 * or a custom ordering function (the 3-arg form) all degrade to the List path.
 * The stable argsort (ties by original position) matches the interpreter's
 * Ordering exactly, which the Compile[] delegation depends on. */
Expr* ndstruct_ordering(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc > 2) return ndarray_delist_and_reeval(res);   /* custom comparator */
    Expr* a = res->data.function.args[0];
    NDType dt = a->data.ndarray.dtype;
    if (a->data.ndarray.rank != 1 || ndt_is_complex(dt))
        return ndarray_delist_and_reeval(res);

    size_t n = (size_t)a->data.ndarray.dims[0];
    int64_t* idx = malloc(sizeof(int64_t) * (n ? n : 1));
    if (!idx) return ndarray_delist_and_reeval(res);

    bool ok;
    if (dt == NDT_INT64) {
        /* Argsort in int64 -- the double gather is exact only to 2^53, so two
         * large integers would tie and the permutation would order them wrong. */
        ok = nd_argsort_i64((const int64_t*)a->data.ndarray.data, n, idx);
    } else {
        double* vals = malloc(sizeof(double) * (n ? n : 1));
        if (!vals) { free(idx); return ndarray_delist_and_reeval(res); }
        nd_gather_real(a->data.ndarray.data, dt, 0, 1, n, vals);
        ok = nd_argsort_real(vals, n, idx);
        free(vals);
    }
    if (!ok) { free(idx); return ndarray_delist_and_reeval(res); }

    for (size_t i = 0; i < n; i++) idx[i] += 1;            /* 1-based positions */

    if (argc == 2) {                                       /* Ordering[a, seq] */
        int64_t* sel = NULL;
        size_t sel_count = 0;
        if (!get_seq_spec_indices(res->data.function.args[1], (int64_t)n, &sel, &sel_count)) {
            free(idx);
            return ndarray_delist_and_reeval(res);
        }
        /* An empty selection has no buffer shape (dims must be >= 1); let the
         * List path answer {} for it. */
        if (sel_count == 0) { free(sel); free(idx); return ndarray_delist_and_reeval(res); }
        int64_t* out = malloc(sizeof(int64_t) * sel_count);
        if (!out) { free(sel); free(idx); return ndarray_delist_and_reeval(res); }
        for (size_t i = 0; i < sel_count; i++) out[i] = idx[sel[i] - 1];
        free(sel);
        free(idx);
        int64_t dims[1] = { (int64_t)sel_count };
        return expr_new_ndarray_like(a, 1, dims, out, NDT_INT64);
    }

    int64_t dims[1] = { (int64_t)n };
    return expr_new_ndarray_like(a, 1, dims, idx, NDT_INT64);
}

/* ---------------------------------------------------------------- Reverse */

/* Reverse[a]: reverse the order of the leading-axis rows (each row a contiguous
 * trailing block). A level/axis spec (2-arg) degrades. */
Expr* ndstruct_reverse(Expr* res) {
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;
    NDType dt = a->data.ndarray.dtype;
    const char* buf = (const char*)a->data.ndarray.data;

    size_t blocks = (size_t)dims[0];
    size_t rowbytes = nd_prod(dims, 1, rank) * ndt_elem_size(dt);
    char* out = malloc(rowbytes * blocks);
    if (!out) return ndarray_delist_and_reeval(res);
    for (size_t b = 0; b < blocks; b++)
        memcpy(out + b * rowbytes, buf + (blocks - 1 - b) * rowbytes, rowbytes);
    return expr_new_ndarray_like(a, rank, dims, out, dt);
}

/* ----------------------------------------------------------------- Rotate */

/*
 * RotateLeft[a, n] / RotateLeft[a, {n1, n2, ...}] straight on the buffer.
 *
 * WHY THIS EXISTS. RotateLeft had no buffer walk, so it went through
 * ndstruct_delist_repack: materialise the array into one Expr per element, run
 * the generic recursive List rotate, re-sniff, re-pack. On a 512x512 float64
 * matrix that is 262144 allocations and it measured **42.6 ms** -- against
 * 0.50 ms for `u + u` and 0.66 ms for Transpose[u] on the very same array, i.e.
 * ~85x the cost of touching every element once.
 *
 * That single number is what made the classical 5-point Jacobi stencil
 *     Nest[(RotateLeft[u,{1,0}] + RotateRight[u,{1,0}]
 *         + RotateLeft[u,{0,1}] + RotateRight[u,{0,1}])/4. &, u0, 100]
 * cost 23 s: four rotates per sweep x 100 sweeps is 17 s of pure marshalling.
 * Rotation is the idiom for a shifted neighbour in array-style numerical code,
 * so it sits in the inner loop of every stencil, convolution and cyclic-shift
 * algorithm anyone writes.
 *
 * HOW. A rotation is a permutation of whole contiguous blocks, never an
 * element-by-element computation. Viewing the array along axis `ax` as
 * [outer][dims[ax]][inner], output slice i is input slice (i + r) mod dims[ax],
 * copied as inner*elemsize contiguous bytes. Axes are rotated one at a time
 * through a scratch buffer, so a rank-N spec costs one pass per NONZERO axis
 * (all-zero axes are skipped, which is the common {1,0} / {0,1} stencil case).
 * For axis 0 the inner extent is the whole row block, so the pass degenerates
 * to two memcpys.
 *
 * Degrades (returns NULL, caller falls back) on a non-Integer spec entry or a
 * spec longer than the rank -- the List path defines those and this must not
 * guess.
 */
static Expr* nd_rotate_axes(Expr* a, const int64_t* rots, int nrots) {
    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;
    NDType dt = a->data.ndarray.dtype;
    size_t esz = ndt_elem_size(dt);
    size_t total = nd_prod(dims, 0, rank);
    if (total == 0) return NULL;

    char* cur = malloc(total * esz);
    if (!cur) return NULL;
    memcpy(cur, a->data.ndarray.data, total * esz);

    char* tmp = NULL;
    for (int ax = 0; ax < nrots && ax < rank; ax++) {
        int64_t d = dims[ax];
        if (d <= 1) continue;
        int64_t r = rots[ax] % d;
        if (r < 0) r += d;
        if (r == 0) continue;                       /* identity on this axis */

        if (!tmp) {
            tmp = malloc(total * esz);
            if (!tmp) { free(cur); return NULL; }
        }
        size_t outer = nd_prod(dims, 0, ax);
        size_t inner = nd_prod(dims, ax + 1, rank) * esz;
        size_t slice = (size_t)d * inner;           /* bytes per outer block */
        for (size_t o = 0; o < outer; o++) {
            const char* src = cur + o * slice;
            char* dst = tmp + o * slice;
            /* Rotating LEFT by r: out[i] = in[i + r]. Two contiguous runs. */
            size_t head = (size_t)(d - r) * inner;  /* in[r .. d) -> out[0 ..) */
            memcpy(dst, src + (size_t)r * inner, head);
            memcpy(dst + head, src, (size_t)r * inner);
        }
        char* swap = cur; cur = tmp; tmp = swap;
    }
    free(tmp);
    return expr_new_ndarray_like(a, rank, dims, cur, dt);
}

/* Parse a rotation spec into `rots` (caller-provided, capacity >= rank).
 * An omitted spec means 1 on axis 0. `sign` is +1 for RotateLeft, -1 for
 * RotateRight. Returns the number of axes covered, or -1 to degrade. */
static int nd_rotate_spec(const Expr* n_spec, int rank, int64_t sign,
                          int64_t* rots) {
    for (int i = 0; i < rank; i++) rots[i] = 0;
    if (!n_spec) { rots[0] = sign; return 1; }
    if (n_spec->type == EXPR_INTEGER) {
        rots[0] = sign * n_spec->data.integer;
        return 1;
    }
    if (n_spec->type == EXPR_FUNCTION && n_spec->data.function.head->type == EXPR_SYMBOL
        && n_spec->data.function.head->data.symbol.name == SYM_List) {
        size_t k = n_spec->data.function.arg_count;
        if (k > (size_t)rank) return -1;            /* spec deeper than the array */
        for (size_t i = 0; i < k; i++) {
            const Expr* e = n_spec->data.function.args[i];
            if (e->type != EXPR_INTEGER) return -1;
            rots[i] = sign * e->data.integer;
        }
        return (int)k;
    }
    return -1;
}

Expr* ndstruct_rotate(Expr* res, bool left) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    Expr* a = res->data.function.args[0];
    if (!is_ndarray(a)) return NULL;
    int rank = a->data.ndarray.rank;
    if (rank < 1 || rank > NDARRAY_MAX_RANK) return NULL;

    int64_t rots[NDARRAY_MAX_RANK];
    int n = nd_rotate_spec(argc == 2 ? res->data.function.args[1] : NULL,
                           rank, left ? 1 : -1, rots);
    if (n < 0) return NULL;
    return nd_rotate_axes(a, rots, n);
}

/* ------------------------------------- Join / Partition / Differences /    */
/* ------------------------------------- Riffle / PadLeft / PadRight        */

/*
 * Five more structural heads straight on the buffer.
 *
 * All five previously went through ndstruct_delist_repack -- one Expr per
 * element, the generic List implementation, then a re-sniff and re-pack. Against
 * Mathematica 14.0 on the same machine that cost (Mathilda / Mathematica):
 *
 *     Differences, 10^6 reals      834 ms / 0.52 ms     1600x
 *     Partition[v, 2], 10^6        257 ms / 0.61 ms      420x
 *     Join, two 10^6 vectors       352 ms / 2.5 ms       140x
 *     Riffle[v, 0.], 10^6          206 ms / 17.6 ms       12x
 *     PadRight, 1000x1000           99 ms / 51 ms          2x
 *
 * None of that was algorithmic -- every one of these is a block move plus, for
 * Differences, one subtraction per element. The marshalling WAS the cost.
 *
 * THE EXACTNESS RULE, which is what most of the care below is about. A head that
 * introduces a NEW element (Riffle's separator, Pad's fill) may only use the
 * buffer when that element is exactly representable at the buffer's dtype AND
 * has the matching exactness. PadLeft[{1.,2.,3.}, 5] is {0, 0, 1., 2., 3.} --
 * exact Integer zeros beside Reals -- in Mathilda AND in Mathematica, and no
 * uniform buffer can hold that, so it must decline and stay a plain List.
 * (Mathematica's own PadRight is unpacked and slow for exactly this reason.)
 * Riffle[v, 0.] with a Real separator into a float64 buffer is fine; Riffle[v, 0]
 * is not. Same rule as pack_repack_like's re-sniff, applied before the work
 * instead of after.
 *
 * Each returns NULL to decline, leaving the caller's existing
 * ndstruct_delist_repack fallback to define every case not handled here.
 */

/* The scalar `e` as a raw buffer element of dtype `dt`, but only when it is
 * exactly representable AND its head matches what `dt` yields on readback:
 * EXPR_REAL for a float dtype, EXPR_INTEGER for NDT_INT64. Anything else --
 * an Integer into float64, a Real into int64, Rational, BigInt, MPFR, Complex,
 * symbol -- is refused, because it would change an element's head. */
static bool nd_fill_value(const Expr* e, NDType dt, double* re, int64_t* iv) {
    if (!e) return false;
    if (dt == NDT_INT64) {
        if (e->type != EXPR_INTEGER) return false;
        *iv = e->data.integer; return true;
    }
    if (ndt_is_complex(dt)) return false;      /* complex fill: not yet */
    if (e->type != EXPR_REAL) return false;
    *re = e->data.real; return true;
}

/* Write `count` copies of the fill value at element offset `at`. */
static void nd_fill_run(void* buf, NDType dt, size_t at, size_t count,
                        double re, int64_t iv) {
    if (dt == NDT_INT64) {
        int64_t* d = (int64_t*)buf;
        for (size_t i = 0; i < count; i++) d[at + i] = iv;
    } else {
        for (size_t i = 0; i < count; i++) ndt_set(buf, at + i, dt, re, 0.0);
    }
}

/* Join[a1, a2, ...] along the leading axis.
 *
 * Every argument must be an ndarray of the SAME dtype and the same trailing
 * dims. Same dtype rather than a widening promotion on purpose: Join of an
 * int64 and a float64 buffer is a list of Integers followed by Reals in the
 * interpreter, and widening would silently turn the Integers into Reals. */
/* The concatenation itself, over an argument VECTOR rather than a call node, so
 * Catenate can reach it: `Catenate[{v, w}]` is `Join[v, w]` with the operands
 * one level deeper, and reproducing 60 lines of pack-lift-and-memcpy to say that
 * would be two implementations of one thing. */
static Expr* nd_concat_leading(Expr** in, size_t argc);

Expr* ndstruct_join(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;
    return nd_concat_leading(res->data.function.args, argc);
}

static Expr* nd_concat_leading(Expr** in, size_t argc) {

    /* PACK THE SMALL OPERANDS UP. Join's cost is set by its LARGEST operand,
     * but the packing decision was made about each one in isolation -- so a
     * two-element boundary list beside a 100000-element buffer used to send the
     * whole call down the List path and materialise the buffer. That is exactly
     * how an explicit finite-difference sweep is written:
     *     v = Join[{lo}, interior, {hi}]
     * once per time step, and it was the dominant cost of a 25000-step American
     * option pricer. Same rule as Dot, Plus/Times, Outer and the Listable gate:
     * pack the small operand up, never materialise the large one down.
     *
     * The lifted list must SNIFF to the same dtype -- never coerce. Join[{1},
     * realBuffer] is a mixed exact/inexact answer in Mathematica, and forcing
     * the exact 1 into a float64 slot would silently make it 1.. An exact list
     * therefore sniffs to int64, fails the dtype test, and the whole call
     * declines to the List path, which gives the mixed answer. */
    Expr** args = malloc(sizeof(Expr*) * argc);
    bool* owned = calloc(argc, sizeof(bool));
    if (!args || !owned) { free(args); free(owned); return NULL; }
    for (size_t i = 0; i < argc; i++) args[i] = in[i];

    Expr* a0 = NULL;
    for (size_t i = 0; i < argc && !a0; i++)
        if (is_ndarray(args[i])) a0 = args[i];
    if (!a0) { free(args); free(owned); return NULL; }

    int rank = a0->data.ndarray.rank;
    NDType dt = a0->data.ndarray.dtype;
    size_t esz = ndt_elem_size(dt);
    size_t rowelts = nd_prod(a0->data.ndarray.dims, 1, rank);

    if (pack_enabled()) {
        for (size_t i = 0; i < argc; i++) {
            if (is_ndarray(args[i])) continue;
            Expr* p = pack_force(expr_copy(args[i]), false, dt);
            if (is_ndarray(p) && p->data.ndarray.dtype == dt) {
                args[i] = p; owned[i] = true;
            } else {
                expr_free(p);
                break;                      /* declines below on the type test */
            }
        }
    }

    int64_t lead = 0;
    bool ok = true;
    for (size_t i = 0; i < argc && ok; i++) {
        Expr* a = args[i];
        if (!is_ndarray(a) || a->data.ndarray.rank != rank
            || a->data.ndarray.dtype != dt) { ok = false; break; }
        for (int d = 1; d < rank; d++)
            if (a->data.ndarray.dims[d] != a0->data.ndarray.dims[d]) { ok = false; break; }
        if (ok) lead += a->data.ndarray.dims[0];
    }

    char* out = NULL;
    if (ok) {
        size_t total_bytes = (size_t)lead * rowelts * esz;
        out = malloc(total_bytes ? total_bytes : 1);  /* malloc(0) is impl-defined */
        if (!out) ok = false;
    }
    if (!ok) {
        for (size_t i = 0; i < argc; i++) if (owned[i]) expr_free(args[i]);
        free(args); free(owned);
        return NULL;
    }

    size_t off = 0;
    for (size_t i = 0; i < argc; i++) {
        Expr* a = args[i];
        size_t bytes = (size_t)a->data.ndarray.dims[0] * rowelts * esz;
        memcpy(out + off, a->data.ndarray.data, bytes);
        off += bytes;
    }
    int64_t odims[NDARRAY_MAX_RANK];
    for (int d = 0; d < rank; d++) odims[d] = a0->data.ndarray.dims[d];
    odims[0] = lead;

    for (size_t i = 0; i < argc; i++) if (owned[i]) expr_free(args[i]);
    free(args); free(owned);
    return expr_new_ndarray_like(a0, rank, odims, out, dt);
}

/* Partition[a, k] on a rank-1 buffer: rank 2, dims {n/k, k}, one memcpy of the
 * whole usable prefix. The offset form Partition[a, k, d] and rank >= 2 decline. */
/* Partition[a, k] and Partition[a, k, d] on a rank-1 buffer: rank 2, dims
 * {rows, k}. With d == k (the default) the rows tile the input and the whole
 * thing is one memcpy; with any other d they OVERLAP and it is `rows` strided
 * copies of k elements each.
 *
 * The offset form is the sliding window -- rolling statistics, n-grams, time
 * series embedding, the setup for a correlation -- and it had no buffer path at
 * all, so `Partition[seq, 12, 1]` over 500000 integers materialised 6 million
 * elements and cost 439 ms. That single call was the whole of a k-mer counting
 * kernel; the Dot that consumed it costs 5.3 ms and the Union after that 11.7.
 *
 * Only complete rows are produced, matching the 2- and 3-argument List path.
 * The padded forms (Partition[a, k, d, spec] and a cyclic overhang) take four
 * or more arguments and are declined here. */
Expr* ndstruct_partition(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 2 && argc != 3) return NULL;
    Expr* a = res->data.function.args[0];
    Expr* ks = res->data.function.args[1];
    if (!is_ndarray(a) || a->data.ndarray.rank != 1) return NULL;
    if (ks->type != EXPR_INTEGER || ks->data.integer <= 0) return NULL;
    int64_t k = ks->data.integer;
    int64_t d = k;
    if (argc == 3) {
        Expr* ds = res->data.function.args[2];
        if (ds->type != EXPR_INTEGER || ds->data.integer <= 0) return NULL;
        d = ds->data.integer;
    }
    int64_t n = a->data.ndarray.dims[0];
    if (n < k) return NULL;                    /* {} -- let the List path answer */
    int64_t rows = (n - k) / d + 1;            /* complete rows only */
    if (rows <= 0) return NULL;
    NDType dt = a->data.ndarray.dtype;
    size_t esz = ndt_elem_size(dt);
    size_t rowbytes = (size_t)k * esz;
    char* out = malloc((size_t)rows * rowbytes);
    if (!out) return NULL;
    const char* in = (const char*)a->data.ndarray.data;
    if (d == k) {
        memcpy(out, in, (size_t)rows * rowbytes);      /* contiguous tiling */
    } else {
        for (int64_t r = 0; r < rows; r++)
            memcpy(out + (size_t)r * rowbytes, in + (size_t)(r * d) * esz, rowbytes);
    }
    int64_t odims[2] = { rows, k };
    return expr_new_ndarray_like(a, 2, odims, out, dt);
}

/* Differences[a] on a rank-1 buffer: n-1 successive differences.
 *
 * int64 subtracts in int64 and abandons the whole result on the first overflow
 * (ci_sub_i64), so the List path re-runs it and GMP gives the exact answer --
 * never a wrapped one. Higher-order and level-spec forms decline. */
Expr* ndstruct_differences(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* a = res->data.function.args[0];
    if (!is_ndarray(a) || a->data.ndarray.rank != 1) return NULL;
    NDType dt = a->data.ndarray.dtype;
    if (ndt_is_complex(dt)) return NULL;
    int64_t n = a->data.ndarray.dims[0];
    if (n < 2) return NULL;                    /* {} -- let the List path answer */
    int64_t m = n - 1;
    int64_t odims[1] = { m };

    if (dt == NDT_INT64) {
        const int64_t* s = (const int64_t*)a->data.ndarray.data;
        int64_t* o = malloc(sizeof(int64_t) * (size_t)m);
        if (!o) return NULL;
        for (int64_t i = 0; i < m; i++) {
            /* ci_sub_i64 is __builtin_sub_overflow: it returns TRUE on
             * OVERFLOW, not on success. Abandon the whole result on the first
             * one so the List path re-runs it and GMP answers exactly. */
            if (ci_sub_i64(s[i + 1], s[i], &o[i])) { free(o); return NULL; }
        }
        return expr_new_ndarray_like(a, 1, odims, o, dt);
    }
    if (dt == NDT_FLOAT64) {
        /* The ndt_get/ndt_set choke point is the right default -- it is the one
         * place that knows every dtype -- but on float64 it is two indirect
         * calls where the loop is one subtract, and the compiler can neither
         * inline through it nor vectorise around it. float64 is also the dtype
         * essentially all of this data is in. Same values either way. */
        const double* sd = (const double*)a->data.ndarray.data;
        double* od = malloc(sizeof(double) * (size_t)m);
        if (!od) return NULL;
        for (int64_t i = 0; i < m; i++) od[i] = sd[i + 1] - sd[i];
        return expr_new_ndarray_like(a, 1, odims, od, dt);
    }
    void* o = malloc(ndt_elem_size(dt) * (size_t)m);
    if (!o) return NULL;
    for (int64_t i = 0; i < m; i++) {
        double hi, lo, dummy;
        ndt_get(a->data.ndarray.data, (size_t)(i + 1), dt, &hi, &dummy);
        ndt_get(a->data.ndarray.data, (size_t)i, dt, &lo, &dummy);
        ndt_set(o, (size_t)i, dt, hi - lo, 0.0);
    }
    return expr_new_ndarray_like(a, 1, odims, o, dt);
}

/* Riffle[a, x] on a rank-1 buffer: a[1], x, a[2], x, ..., a[n] -- 2n-1 elements,
 * no trailing separator. `x` must pass nd_fill_value; a List separator (which
 * cycles) declines. */
Expr* ndstruct_riffle(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* a = res->data.function.args[0];
    Expr* sep = res->data.function.args[1];
    if (!is_ndarray(a) || a->data.ndarray.rank != 1) return NULL;
    NDType dt = a->data.ndarray.dtype;
    double re = 0.0; int64_t iv = 0;
    if (!nd_fill_value(sep, dt, &re, &iv)) return NULL;
    int64_t n = a->data.ndarray.dims[0];
    if (n < 1) return NULL;
    size_t esz = ndt_elem_size(dt);
    int64_t m = 2 * n - 1;
    char* out = malloc((size_t)m * esz);
    if (!out) return NULL;
    const char* s = (const char*)a->data.ndarray.data;
    for (int64_t i = 0; i < n; i++) {
        memcpy(out + (size_t)(2 * i) * esz, s + (size_t)i * esz, esz);
        if (i + 1 < n) nd_fill_run(out, dt, (size_t)(2 * i + 1), 1, re, iv);
    }
    int64_t odims[1] = { m };
    return expr_new_ndarray_like(a, 1, odims, out, dt);
}

/* PadLeft[a, n] / PadRight[a, n] / with an explicit fill, on a rank-1 buffer.
 *
 * The fill is only consulted when padding actually happens: PadLeft[a, n] with
 * n <= Length[a] is pure truncation, which needs no fill and so stays on the
 * buffer even with the default exact 0 that would otherwise force a decline.
 * A List size spec (rank >= 2 padding) declines. */
Expr* ndstruct_pad(Expr* res, bool left) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;
    Expr* a = res->data.function.args[0];
    Expr* ns = res->data.function.args[1];
    if (!is_ndarray(a) || a->data.ndarray.rank != 1) return NULL;
    if (ns->type != EXPR_INTEGER || ns->data.integer < 0) return NULL;
    int64_t want = ns->data.integer;
    int64_t n = a->data.ndarray.dims[0];
    NDType dt = a->data.ndarray.dtype;
    size_t esz = ndt_elem_size(dt);
    const char* s = (const char*)a->data.ndarray.data;
    int64_t odims[1] = { want };

    if (want <= n) {
        /* Truncation: PadLeft keeps the LAST `want`, PadRight the first. */
        if (want == 0) return NULL;            /* {} -- let the List path answer */
        char* out = malloc((size_t)want * esz);
        if (!out) return NULL;
        memcpy(out, s + (left ? (size_t)(n - want) * esz : 0), (size_t)want * esz);
        return expr_new_ndarray_like(a, 1, odims, out, dt);
    }

    double re = 0.0; int64_t iv = 0;
    const Expr* fill = (argc == 3) ? res->data.function.args[2] : NULL;
    if (argc == 3) {
        if (!nd_fill_value(fill, dt, &re, &iv)) return NULL;
    } else if (dt != NDT_INT64) {
        /* Default fill is the exact Integer 0. In a float64 buffer that would
         * become 0.0, changing {0, 0, 1., 2., 3.} into {0., 0., 1., 2., 3.} --
         * a different value. Decline; the List path produces the mixed answer,
         * which is what Mathematica gives too. */
        return NULL;
    }
    char* out = malloc((size_t)want * esz);
    if (!out) return NULL;
    int64_t pad = want - n;
    if (left) {
        nd_fill_run(out, dt, 0, (size_t)pad, re, iv);
        memcpy(out + (size_t)pad * esz, s, (size_t)n * esz);
    } else {
        memcpy(out, s, (size_t)n * esz);
        nd_fill_run(out, dt, (size_t)n, (size_t)pad, re, iv);
    }
    return expr_new_ndarray_like(a, 1, odims, out, dt);
}

/* -------------------------------------------------------------- Transpose */

/* Transpose[a] for a rank-2 array: out[c, r] = in[r, c]. Higher rank or an
 * explicit permutation spec degrades. */
/* Cache-blocked transpose tile size (elements). A 32x32 f64 tile is 8 KB per
 * source/dest window — comfortably inside L1, so each tile is read and written
 * with unit stride on one side and short strided hops on the other, instead of
 * the full-matrix stride that thrashes cache on a naive i/j double loop. */
#define ND_TRANSPOSE_TILE 32

Expr* ndstruct_transpose(Expr* res) {
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    if (a->data.ndarray.rank != 2) return ndarray_delist_and_reeval(res);
    const int64_t* dims = a->data.ndarray.dims;
    NDType dt = a->data.ndarray.dtype;
    const void* buf = a->data.ndarray.data;
    size_t R = (size_t)dims[0], C = (size_t)dims[1];   /* in[r,c] -> out[c,r] */

    void* out = malloc(ndt_elem_size(dt) * R * C);
    if (!out) return ndarray_delist_and_reeval(res);
    const size_t T = ND_TRANSPOSE_TILE;
    if (dt == NDT_FLOAT64) {                            /* raw-double blocked path */
        const double* p = (const double*)buf;
        double* o = (double*)out;
        for (size_t rr = 0; rr < R; rr += T)
            for (size_t cc = 0; cc < C; cc += T) {
                size_t rmax = rr + T < R ? rr + T : R;
                size_t cmax = cc + T < C ? cc + T : C;
                for (size_t r = rr; r < rmax; r++)
                    for (size_t c = cc; c < cmax; c++)
                        o[c * R + r] = p[r * C + c];
            }
    } else {                                            /* dtype-generic blocked */
        for (size_t rr = 0; rr < R; rr += T)
            for (size_t cc = 0; cc < C; cc += T) {
                size_t rmax = rr + T < R ? rr + T : R;
                size_t cmax = cc + T < C ? cc + T : C;
                /* memcpy per element, not a ndt_get/ndt_set round trip: a
                 * transpose only MOVES values, and routing them through a double
                 * would round an int64 buffer past 2^53. */
                size_t esz = ndt_elem_size(dt);
                for (size_t r = rr; r < rmax; r++)
                    for (size_t c = cc; c < cmax; c++)
                        memcpy((char*)out + (c * R + r) * esz,
                               (const char*)buf + (r * C + c) * esz, esz);
            }
    }
    int64_t odims[2] = { (int64_t)C, (int64_t)R };
    return expr_new_ndarray_like(a, 2, odims, out, dt);
}

/* --------------------------------------------------------------- Diagonal */

/* Diagonal[a] / Diagonal[a, k]: the k-th diagonal of a rank-2 buffer, returned
 * as a fresh rank-1 NDArray. k defaults to 0 (leading diagonal); k > 0 selects a
 * superdiagonal, k < 0 a subdiagonal. Elements are moved by memcpy (never through
 * ndt_get), so an int64 buffer stays exact past 2^53 -- a diagonal only MOVES
 * values, exactly like Transpose.
 *
 * Degrades (ndarray_delist_and_reeval) for a non-2-D operand, a non-integer k, an
 * arity outside {1, 2}, or an EMPTY diagonal (|k| beyond the matrix): the List
 * path answers {} there, which no rank-1 buffer shape can express. The empty
 * guard `k <= -R || k >= C` also keeps the -k below from overflowing at k=INT64_MIN. */
Expr* ndstruct_diagonal(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    if (a->data.ndarray.rank != 2) return ndarray_delist_and_reeval(res);

    int64_t k = 0;
    if (argc == 2) {
        Expr* ke = res->data.function.args[1];
        if (ke->type != EXPR_INTEGER) return ndarray_delist_and_reeval(res);
        k = ke->data.integer;
    }

    const int64_t* dims = a->data.ndarray.dims;
    int64_t R = dims[0], C = dims[1];
    if (k <= -R || k >= C) return ndarray_delist_and_reeval(res);   /* empty diagonal */

    int64_t start_row = (k < 0) ? -k : 0;
    int64_t start_col = (k > 0) ?  k : 0;
    int64_t len_r = R - start_row, len_c = C - start_col;
    int64_t len = (len_r < len_c) ? len_r : len_c;                  /* guaranteed >= 1 */

    NDType dt = a->data.ndarray.dtype;
    size_t esz = ndt_elem_size(dt);
    void* out = malloc(esz * (size_t)len);
    if (!out) return ndarray_delist_and_reeval(res);
    const char* buf = (const char*)a->data.ndarray.data;
    for (int64_t t = 0; t < len; t++) {
        int64_t i = start_row + t, j = start_col + t;
        memcpy((char*)out + (size_t)t * esz,
               buf + ((size_t)i * (size_t)C + (size_t)j) * esz, esz);
    }
    int64_t odims[1] = { len };
    return expr_new_ndarray_like(a, 1, odims, out, dt);
}

/* ---------------------------------------------------------------- Flatten */

/* Flatten[a]: collapse every axis into a single rank-1 array. The buffer is
 * already contiguous row-major, so this is a straight copy + reshape. A partial
 * Flatten[a, n] (or a head arg) degrades. */
Expr* ndstruct_flatten(Expr* res) {
    if (res->data.function.arg_count != 1) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    NDType dt = a->data.ndarray.dtype;
    size_t sz = ndarray_size(a);
    size_t bytes = ndt_elem_size(dt) * sz;
    void* out = malloc(bytes);
    if (!out) return ndarray_delist_and_reeval(res);
    memcpy(out, a->data.ndarray.data, bytes);
    int64_t odims[1] = { (int64_t)sz };
    return expr_new_ndarray_like(a, 1, odims, out, dt);
}

/* ------------------------------------------------------------- Take / Drop */

/* Build a new NDArray of `count` leading rows starting at row `start`, copied
 * from `a` (rows are contiguous trailing blocks). */
static Expr* nd_rows(const Expr* a, size_t start, size_t count) {
    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;
    NDType dt = a->data.ndarray.dtype;
    size_t rowbytes = nd_prod(dims, 1, rank) * ndt_elem_size(dt);
    char* out = malloc(rowbytes * count);
    if (!out) return NULL;
    memcpy(out, (const char*)a->data.ndarray.data + start * rowbytes, rowbytes * count);
    int64_t odims[NDARRAY_MAX_RANK];
    for (int i = 0; i < rank; i++) odims[i] = dims[i];
    odims[0] = (int64_t)count;
    return expr_new_ndarray_like(a, rank, odims, out, dt);
}

/* ------------------------------------------------- First / Last / Most / Rest
 *
 * All four are leading-axis slices of exactly the shape nd_rows already cuts —
 * one row for First/Last, all-but-one for Most/Rest — so each is a pointer read
 * or a single memcpy. They reached the buffer through ndstruct_delist_repack
 * instead, which materialises every element to build a List and then repacks the
 * answer, and that made First and Last **asymptotically wrong**: reading one
 * element of a 10^6 float64 vector cost 123 ms, against 0.88 ms for the
 * identical Drop[v, 250] beside it and an unmeasurable v[0] in NumPy.
 *
 * A rank-1 First/Last yields a SCALAR, so it leaves the buffer entirely and goes
 * through ndarray_element_to_expr — the one place that decides an element's head
 * (ndarray.h:87), so an int64 buffer still yields an exact Integer.
 *
 * Empty results degrade: Rest and Most of a one-row array are {}, and an empty
 * NDArray is not a shape this layer builds. */
Expr* ndstruct_head_tail(Expr* res, NDHeadTail which) {
    Expr* a = res->data.function.args[0];
    /* First[a, default] / Last[a, default] never reach the default on a packed
     * array (it is non-empty by construction), so the 2-arg form is fine here;
     * every other arity is not ours. */
    size_t argc = res->data.function.arg_count;
    bool is_end = (which == ND_FIRST || which == ND_LAST);
    if (argc != 1 && !(is_end && argc == 2)) return NULL;

    int rank = a->data.ndarray.rank;
    size_t blocks = (size_t)a->data.ndarray.dims[0];
    if (blocks == 0) return NULL;

    switch (which) {
    case ND_FIRST:
    case ND_LAST: {
        size_t row = (which == ND_FIRST) ? 0 : blocks - 1;
        if (rank == 1) return ndarray_element_to_expr(a, row);
        /* Rank >= 2: one row, with the (now length-1) leading axis DROPPED —
         * First of a 3x4 is a length-4 vector, not a 1x4 matrix. */
        const int64_t* dims = a->data.ndarray.dims;
        NDType dt = a->data.ndarray.dtype;
        size_t rowbytes = nd_prod(dims, 1, rank) * ndt_elem_size(dt);
        char* out = malloc(rowbytes);
        if (!out) return NULL;
        memcpy(out, (const char*)a->data.ndarray.data + row * rowbytes, rowbytes);
        int64_t odims[NDARRAY_MAX_RANK];
        for (int i = 1; i < rank; i++) odims[i - 1] = dims[i];
        return expr_new_ndarray_like(a, rank - 1, odims, out, dt);
    }
    case ND_REST:
    case ND_MOST: {
        if (blocks < 2) return NULL;   /* -> {} */
        return nd_rows(a, which == ND_REST ? 1 : 0, blocks - 1);
    }
    }
    return NULL;
}

/* Resolve a leading-axis spec into a 1-based inclusive [lo, hi] range over
 * `blocks` rows. Handles an integer n (n>0: first n; n<0: last |n|) and a list
 * {i} / {i, j} with 1-based, possibly-negative endpoints. Returns false (caller
 * degrades) for anything else or an out-of-range / empty selection. */
static bool nd_span(const Expr* spec, size_t blocks, int64_t* lo, int64_t* hi) {
    if (spec->type == EXPR_INTEGER) {
        int64_t n = spec->data.integer;
        if (n > 0)      { *lo = 1; *hi = n; }
        else if (n < 0) { *lo = (int64_t)blocks + n + 1; *hi = (int64_t)blocks; }
        else return false;                            /* n == 0 */
    } else if (spec->type == EXPR_FUNCTION &&
               spec->data.function.head->type == EXPR_SYMBOL &&
               spec->data.function.head->data.symbol.name == SYM_List &&
               (spec->data.function.arg_count == 1 || spec->data.function.arg_count == 2)) {
        Expr* e0 = spec->data.function.args[0];
        Expr* e1 = spec->data.function.args[spec->data.function.arg_count == 2 ? 1 : 0];
        if (e0->type != EXPR_INTEGER || e1->type != EXPR_INTEGER) return false;
        int64_t i = e0->data.integer, j = e1->data.integer;
        *lo = (i > 0) ? i : (int64_t)blocks + i + 1;
        *hi = (j > 0) ? j : (int64_t)blocks + j + 1;
    } else {
        return false;
    }
    if (*lo < 1 || *hi > (int64_t)blocks || *lo > *hi) return false;
    return true;
}

Expr* ndstruct_take(Expr* res) {
    if (res->data.function.arg_count != 2) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    size_t blocks = (size_t)a->data.ndarray.dims[0];
    int64_t lo, hi;
    if (!nd_span(res->data.function.args[1], blocks, &lo, &hi))
        return ndarray_delist_and_reeval(res);
    Expr* r = nd_rows(a, (size_t)(lo - 1), (size_t)(hi - lo + 1));
    return r ? r : ndarray_delist_and_reeval(res);
}

Expr* ndstruct_drop(Expr* res) {
    if (res->data.function.arg_count != 2) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;
    NDType dt = a->data.ndarray.dtype;
    size_t blocks = (size_t)dims[0];
    int64_t lo, hi;
    if (!nd_span(res->data.function.args[1], blocks, &lo, &hi))
        return ndarray_delist_and_reeval(res);

    /* Keep rows [1, lo) and (hi, blocks] (1-based); drop the middle span. */
    size_t before = (size_t)(lo - 1);
    size_t after = blocks - (size_t)hi;
    size_t kept = before + after;
    if (kept == 0) return ndarray_delist_and_reeval(res);   /* empty -> List path */
    size_t rowbytes = nd_prod(dims, 1, rank) * ndt_elem_size(dt);
    const char* buf = (const char*)a->data.ndarray.data;
    char* out = malloc(rowbytes * kept);
    if (!out) return ndarray_delist_and_reeval(res);
    memcpy(out, buf, rowbytes * before);
    memcpy(out + rowbytes * before, buf + rowbytes * (size_t)hi, rowbytes * after);
    int64_t odims[NDARRAY_MAX_RANK];
    for (int i = 0; i < rank; i++) odims[i] = dims[i];
    odims[0] = (int64_t)kept;
    return expr_new_ndarray_like(a, rank, odims, out, dt);
}

/* ------------------------------------------------------------------- Clip */

/* Elementwise clamp chunk. float64 clamps a raw double buffer (vectorizes);
 * other real dtypes go through the ndt_get/ndt_set choke point. */
typedef struct { const void* buf; void* out; NDType dt; double lo, hi; } nd_clip_ctx;
static bool nd_clip_chunk(void* c, size_t lo_k, size_t hi_k) {
    const nd_clip_ctx* x = (const nd_clip_ctx*)c;
    double lo = x->lo, hi = x->hi;
    if (x->dt == NDT_FLOAT64) {
        const double* p = (const double*)x->buf;
        double* o = (double*)x->out;
        for (size_t k = lo_k; k < hi_k; k++) {
            double r = p[k];
            r = r < lo ? lo : (r > hi ? hi : r);
            o[k] = r;
        }
    } else {
        for (size_t k = lo_k; k < hi_k; k++) {
            double r, im;
            ndt_get(x->buf, k, x->dt, &r, &im);
            r = r < lo ? lo : (r > hi ? hi : r);
            ndt_set(x->out, k, x->dt, r, 0.0);
        }
    }
    return true;
}

/* Clip[a] clamps to [-1, 1]; Clip[a, {min, max}] to [min, max]. Elementwise,
 * real dtype only, threaded. The 3-arg replacement form and complex dtypes degrade.
 *
 * THE EXACTNESS GATE. Clip returns the BOUND at every clipped position, with the
 * bound's own head — so on a Real array an exact bound produces exact Integers:
 *
 *     Clip[{0.5, 1.5, -2.5}, {1, 2}]   ->  {1, 1.5, 1}     (mixed, unpackable)
 *     Clip[{0.5, 1.5, -2.5}, {1., 2.}] ->  {1., 1.5, 1.}   (uniform float64)
 *     Clip[{0.5, 1.5, -2.5}]           ->  {0.5, 1, -1}    (default bounds are exact)
 *
 * Clip was removed from pack.c's AWARE list over exactly this, which cost the
 * head its buffer path entirely: 356 ms on a 10^6 float64 vector against 4.7 ms
 * for np.clip. The correct gate is narrower than the head — Real bounds are
 * uniform and safe, and an EXACT bound is still safe when NOTHING is clipped,
 * because then the answer just IS the input. One in-range scan decides that, and
 * it is the same pass the clamp would have run. */
Expr* ndstruct_clip(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return ndarray_delist_and_reeval(res);
    Expr* a = res->data.function.args[0];
    NDType dt = a->data.ndarray.dtype;
    /* An integer Clip is exact and its bounds may themselves be exact
     * (Clip[Range[10], {3, 7}] is a list of Integers) -- but the clamp runs
     * through doubles below, so hand integer buffers to the List path rather than
     * round them. */
    if (ndt_is_complex(dt) || dt == NDT_INT64) return ndarray_delist_and_reeval(res);

    double lo = -1.0, hi = 1.0;
    /* The 1-arg default bounds are the exact Integers -1 and 1. */
    bool bounds_real = false;
    if (argc == 2) {
        Expr* iv = res->data.function.args[1];
        bool lo_inf = false, hi_inf = false;
        if (iv->type != EXPR_FUNCTION ||
            iv->data.function.head->type != EXPR_SYMBOL ||
            iv->data.function.head->data.symbol.name != SYM_List ||
            iv->data.function.arg_count != 2 ||
            !nd_real_value(iv->data.function.args[0], &lo, &lo_inf) ||
            !nd_real_value(iv->data.function.args[1], &hi, &hi_inf))
            return ndarray_delist_and_reeval(res);
        /* Per side: a Real bound is safe because the clipped positions get a
         * Real; an INFINITE bound is safe because no position is ever clipped to
         * it. Anything else (an exact Integer or Rational bound) can put its own
         * exact head into a Real answer, and falls to the in-range scan below. */
        bounds_real = (iv->data.function.args[0]->type == EXPR_REAL || lo_inf) &&
                      (iv->data.function.args[1]->type == EXPR_REAL || hi_inf);
    }

    size_t sz = ndarray_size(a);

    if (!bounds_real) {
        /* An exact bound would put an Integer into the result at every clipped
         * position. Safe only if there is no such position — in which case the
         * result is the input, unchanged and still uniformly Real. */
        for (size_t k = 0; k < sz; k++) {
            double r, im;
            ndt_get(a->data.ndarray.data, k, dt, &r, &im);
            if (r < lo || r > hi) return ndarray_delist_and_reeval(res);
        }
        void* same = malloc(ndt_elem_size(dt) * (sz ? sz : 1));
        if (!same) return ndarray_delist_and_reeval(res);
        memcpy(same, a->data.ndarray.data, ndt_elem_size(dt) * sz);
        return expr_new_ndarray_like(a, a->data.ndarray.rank, a->data.ndarray.dims,
                                     same, dt);
    }

    void* out = malloc(ndt_elem_size(dt) * sz);
    if (!out) return ndarray_delist_and_reeval(res);
    nd_clip_ctx c = { a->data.ndarray.data, out, dt, lo, hi };
    nd_parallel_for(sz, nd_clip_chunk, &c);
    return expr_new_ndarray_like(a, a->data.ndarray.rank, a->data.ndarray.dims, out, dt);
}

/* ------------------------------------------------------------------ *
 *  Delist-and-repack — see ndstruct.h                                 *
 * ------------------------------------------------------------------ */
Expr* ndstruct_delist_repack(const Expr* call, const Expr* src) {
    Expr* out = ndarray_delist_and_reeval(call);
    if (!out || !src || !is_ndarray(src)) return out;
    /* A scalar, a Missing, a symbol — nothing to pack, and the List
     * implementation's answer is already the right one. */
    if (out->type != EXPR_FUNCTION) return out;
    /* A PACKED source re-sniffs the dtype: the operation may have introduced an
     * element of a different exactness, and repacking at the source's dtype would
     * coerce it. Join[Range[1., 300.], {1}] gave 1. for the appended exact 1, and
     * Riffle[Range[1., 300.], 1] turned every interleaved 1 into 1.; both now
     * decline to pack and keep the ordinary List, which is the right value.
     *
     * A VISIBLE NDArray[...] keeps repacking at its own dtype -- coercing what is
     * put into a machine buffer is what naming that head asks for. */
    if (is_packed_list(src)) return pack_repack_like(src, out);
    Expr* packed = ndarray_from_nested_list_like(src, out, src->data.ndarray.dtype);
    if (!packed) return out;            /* symbolic, ragged, or empty */
    expr_free(out);
    return packed;
}

/* ==================================================================== *
 *  The ninth round: heads that had no buffer path on EITHER surface.
 *
 *  Each of these answered a visible NDArray with the unevaluated call and a
 *  packed List with one Expr per element, because there was no fast path to
 *  reach -- the class both static audits are blind to (see ndstruct.h).
 * ==================================================================== */

/* ---------------------------------------------------------------- Ratios */

/* Ratios[a]: the multiplicative Differences, out[i] = a[i+1]/a[i].
 *
 * INEXACT dtypes only, and that is the whole subtlety. Ratios of an int64
 * buffer is a list of exact Rationals -- Ratios[{1, 2, 3}] is {2, 3/2} -- which
 * no machine buffer holds, so an integer argument declines and the List path
 * answers exactly. A zero element declines for the same reason: the answer is
 * ComplexInfinity, not a double. */
Expr* ndstruct_ratios(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* a = res->data.function.args[0];
    if (!is_ndarray(a) || a->data.ndarray.rank != 1) return NULL;
    NDType dt = a->data.ndarray.dtype;
    if (dt == NDT_INT64 || ndt_is_complex(dt)) return NULL;
    int64_t n = a->data.ndarray.dims[0];
    if (n < 2) return NULL;                    /* {} -- let the List path answer */
    int64_t m = n - 1, odims[1];
    odims[0] = m;

    if (dt == NDT_FLOAT64) {
        const double* s = (const double*)a->data.ndarray.data;
        for (int64_t i = 0; i < m; i++) if (s[i] == 0.0) return NULL;
        double* o = malloc(sizeof(double) * (size_t)m);
        if (!o) return NULL;
        for (int64_t i = 0; i < m; i++) o[i] = s[i + 1] / s[i];
        return expr_new_ndarray_like(a, 1, odims, o, dt);
    }
    void* o = malloc(ndt_elem_size(dt) * (size_t)m);
    if (!o) return NULL;
    for (int64_t i = 0; i < m; i++) {
        double hi, lo, dummy;
        ndt_get(a->data.ndarray.data, (size_t)(i + 1), dt, &hi, &dummy);
        ndt_get(a->data.ndarray.data, (size_t)i, dt, &lo, &dummy);
        if (lo == 0.0) { free(o); return NULL; }
        ndt_set(o, (size_t)i, dt, hi / lo, 0.0);
    }
    return expr_new_ndarray_like(a, 1, odims, o, dt);
}

/* ------------------------------------------------------- Append / Prepend */

/* Append[a, x] / Prepend[a, x] on the buffer: one allocation and two memcpys.
 *
 * These were on pack.c's "correct by omission" list -- the gate materialised
 * for them and the List code answered -- which is right for correctness and
 * costs a full 10^6-element boxing to add one element. numpy calls the same
 * thing np.append, and it is how a growing series is written.
 *
 * Two element shapes:
 *   rank 1, x a scalar   nd_fill_value's exactness rule applies verbatim, so
 *                        an exact Integer into a float64 buffer DECLINES.
 *                        Append[Range[1., 10.], 0] is a mixed exact/inexact
 *                        list in Mathematica and coercing the 0 to 0. would be
 *                        a different answer.
 *   rank N, x an array   Append[matrix, row]: x's shape must equal a's
 *                        trailing dims and its dtype must match exactly, the
 *                        same never-widen rule as Join. */
Expr* ndstruct_append(Expr* res, bool prepend) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* a = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (!is_ndarray(a)) return NULL;
    NDType dt = a->data.ndarray.dtype;
    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;
    size_t esz = ndt_elem_size(dt);
    size_t rowelts = nd_prod(dims, 1, rank);
    size_t rowbytes = rowelts * esz;
    size_t oldbytes = (size_t)dims[0] * rowbytes;

    char* out = NULL;
    if (rank == 1) {
        double re = 0.0; int64_t iv = 0;
        if (!nd_fill_value(x, dt, &re, &iv)) return NULL;
        out = malloc(oldbytes + esz);
        if (!out) return NULL;
        size_t at = prepend ? 0 : (size_t)dims[0];
        memcpy(out + (prepend ? esz : 0), a->data.ndarray.data, oldbytes);
        nd_fill_run(out, dt, at, 1, re, iv);
    } else {
        if (!is_ndarray(x) || x->data.ndarray.dtype != dt) return NULL;
        if (x->data.ndarray.rank != rank - 1) return NULL;
        for (int d = 1; d < rank; d++)
            if (x->data.ndarray.dims[d - 1] != dims[d]) return NULL;
        out = malloc(oldbytes + rowbytes);
        if (!out) return NULL;
        if (prepend) {
            memcpy(out, x->data.ndarray.data, rowbytes);
            memcpy(out + rowbytes, a->data.ndarray.data, oldbytes);
        } else {
            memcpy(out, a->data.ndarray.data, oldbytes);
            memcpy(out + oldbytes, x->data.ndarray.data, rowbytes);
        }
    }
    int64_t odims[NDARRAY_MAX_RANK];
    for (int d = 0; d < rank; d++) odims[d] = dims[d];
    odims[0] = dims[0] + 1;
    return expr_new_ndarray_like(a, rank, odims, out, dt);
}

/* -------------------------------------------------------------- Catenate */

/* Catenate[{l1, l2, ...}] = Join[l1, l2, ...], one level down.
 *
 * Two argument shapes reach here and they are the same operation seen twice:
 *
 *   an ARRAY of rank >= 2   the sublists are its rows, and a row-major buffer
 *                           already stores them consecutively -- so catenating
 *                           is a RESHAPE, dims {r, c, ...} -> {r*c, ...}, and
 *                           the copy exists only because inputs are immutable.
 *                           This is the shape a packed argument arrives in:
 *                           pack_sniff absorbs a List of packed vectors into
 *                           one rank-2 array before Catenate is ever called.
 *   a LIST of arrays        which is what a VISIBLE NDArray argument gives,
 *                           since the packing gate never absorbs those. Exactly
 *                           Join over the elements. */
Expr* ndstruct_catenate(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* a = res->data.function.args[0];

    if (is_ndarray(a)) {
        int rank = a->data.ndarray.rank;
        if (rank < 2) return NULL;             /* Catenate of a vector is an error */
        NDType dt = a->data.ndarray.dtype;
        size_t total = nd_prod(a->data.ndarray.dims, 0, rank);
        size_t bytes = total * ndt_elem_size(dt);
        void* out = malloc(bytes ? bytes : 1);
        if (!out) return NULL;
        memcpy(out, a->data.ndarray.data, bytes);
        int64_t odims[NDARRAY_MAX_RANK];
        odims[0] = a->data.ndarray.dims[0] * a->data.ndarray.dims[1];
        for (int d = 2; d < rank; d++) odims[d - 1] = a->data.ndarray.dims[d];
        return expr_new_ndarray_like(a, rank - 1, odims, out, dt);
    }

    if (a->type != EXPR_FUNCTION || a->data.function.head->type != EXPR_SYMBOL ||
        a->data.function.head->data.symbol.name != SYM_List) return NULL;
    size_t argc = a->data.function.arg_count;
    if (argc < 1) return NULL;
    return nd_concat_leading(a->data.function.args, argc);
}

/* ------------------------------------------ TakeLargest / TakeSmallest */

/* A bounded binary heap over (value) doubles, used to keep the running best k.
 * `want_max` selects a MIN-heap (so the root is the weakest of the current best
 * and is what a new candidate displaces) and vice versa. */
static void nd_heap_sift_down(double* h, size_t n, size_t i, bool minheap) {
    for (;;) {
        size_t l = 2 * i + 1, r = l + 1, best = i;
        if (l < n && (minheap ? h[l] < h[best] : h[l] > h[best])) best = l;
        if (r < n && (minheap ? h[r] < h[best] : h[r] > h[best])) best = r;
        if (best == i) return;
        double t = h[i]; h[i] = h[best]; h[best] = t;
        i = best;
    }
}

/* TakeLargest[a, k] / TakeSmallest[a, k] over a real rank-1 buffer.
 *
 * O(n log k), not O(n log n): the whole point of asking for the top ten of a
 * million is that sorting the million is unnecessary. The List path builds n
 * Exprs and sorts them through expr_compare; NumPy's own idiom for this row is
 * a full np.sort, so the heap is ahead of the reference as well as of the
 * interpreter.
 *
 * Declines a NaN element -- an unordered value has no k-th largest, and the
 * List path's answer for it is the one to keep. Ties are by value, and since
 * every element is a machine number two equal elements are indistinguishable,
 * so which one is reported cannot be observed. */
Expr* ndstruct_take_extreme(Expr* res, bool largest) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* a = res->data.function.args[0];
    Expr* kx = res->data.function.args[1];
    if (!is_ndarray(a) || a->data.ndarray.rank != 1) return NULL;
    NDType dt = a->data.ndarray.dtype;
    if (ndt_is_complex(dt)) return NULL;
    /* int64 is exact only to 2^53 through `double`, so ordering it here could
     * reorder two large integers that compare equal. The List path sorts them
     * exactly; leave it to do so. */
    if (dt == NDT_INT64) return NULL;
    if (kx->type != EXPR_INTEGER || kx->data.integer < 0) return NULL;
    size_t n = (size_t)a->data.ndarray.dims[0];
    size_t k = (size_t)kx->data.integer;
    if (k == 0 || k > n) return NULL;     /* Take::take diagnostics: List path */

    double* h = malloc(sizeof(double) * k);
    if (!h) return NULL;
    for (size_t i = 0; i < n; i++) {
        double v, im;
        ndt_get(a->data.ndarray.data, i, dt, &v, &im);
        if (!(v == v)) { free(h); return NULL; }          /* NaN: degrade */
        if (i < k) {
            h[i] = v;
            if (i + 1 == k) for (size_t j = k / 2; j-- > 0; )
                nd_heap_sift_down(h, k, j, largest);
        } else if (largest ? (v > h[0]) : (v < h[0])) {
            h[0] = v;
            nd_heap_sift_down(h, k, 0, largest);
        }
    }
    /* Heap-sort the k survivors in place, and the two cases need no separate
     * code because the heap already differs. `largest` used a MIN-heap, so
     * moving the root to the end of the shrinking region deposits the smallest
     * survivor last and leaves the array DESCENDING -- what TakeLargest wants.
     * `!largest` used a max-heap and leaves it ASCENDING, what TakeSmallest
     * wants. Both are then read front to back. */
    for (size_t m = k; m > 1; m--) {
        double t = h[0]; h[0] = h[m - 1]; h[m - 1] = t;
        nd_heap_sift_down(h, m - 1, 0, largest);
    }
    void* out = malloc(ndt_elem_size(dt) * k);
    if (!out) { free(h); return NULL; }
    for (size_t i = 0; i < k; i++) ndt_set(out, i, dt, h[i], 0.0);
    free(h);
    int64_t odims[1];
    odims[0] = (int64_t)k;
    return expr_new_ndarray_like(a, 1, odims, out, dt);
}
