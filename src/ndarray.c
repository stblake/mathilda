/* NDArray — see ndarray.h for the design rationale. */

#include "ndarray.h"
#include "checked_int.h"   /* ci_*_i64 — exact int64 arithmetic, bail on overflow */
#include "ndarray_internal.h"
#include "sym_names.h"
#include "pack.h"     /* pack_force — pack a small List operand meeting a buffer */
#include "attr.h"
#include "symtab.h"
#include "common.h"
#include "arithmetic.h"   /* arith_warnings_muted() */
#include "eval.h"         /* eval_clock_get() — dedup warnings within one eval */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

bool is_ndarray(const Expr* e) {
    return e && e->type == EXPR_NDARRAY;
}

bool is_packed_list(const Expr* e) {
    return e && e->type == EXPR_NDARRAY &&
           e->data.ndarray.present_as == NDA_HEAD_LIST;
}

const Expr* nd_present_src2(const Expr* a, const Expr* b) {
    /* Visible dominates packed, so scan for a visible array first. */
    if (is_ndarray(a) && !is_packed_list(a)) return a;
    if (is_ndarray(b) && !is_packed_list(b)) return b;
    if (is_ndarray(a)) return a;
    if (is_ndarray(b)) return b;
    return NULL;
}

size_t ndarray_size(const Expr* a) {
    size_t n = 1;
    for (int i = 0; i < a->data.ndarray.rank; i++) {
        n *= (size_t)a->data.ndarray.dims[i];
    }
    return n;
}

/* ------------------------------- dtype helpers --------------------------- */

bool ndt_is_complex(NDType dt) {
    return dt == NDT_COMPLEX64 || dt == NDT_COMPLEX32;
}

int ndt_components(NDType dt) {
    return ndt_is_complex(dt) ? 2 : 1;
}

size_t ndt_comp_size(NDType dt) {
    if (dt == NDT_BOOL) return 1;               /* one uint8_t per element */
    if (dt == NDT_INT64) return sizeof(int64_t);
    return (dt == NDT_FLOAT64 || dt == NDT_COMPLEX64) ? sizeof(double)
                                                      : sizeof(float);
}

size_t ndt_elem_size(NDType dt) {
    return ndt_comp_size(dt) * (size_t)ndt_components(dt);
}

/* ---------------------------------------------------------------------------
 * int64 audit tripwire.
 *
 * ndt_get/ndt_set route through `double` and are therefore exact only to 2^53.
 * Historically that was safe because NDT_INT64 was reachable only from
 * Compile[], which uses the exact ndt_get_i/ndt_set_i pair. Packed lists infer
 * int64 from an all-Integer list (see pack.c), so integer buffers can now reach
 * code that was written when they could not, and every lossy site is a wrong
 * answer waiting for one -- silently, since the values below 2^53 look right.
 *
 * This makes the audit executable instead of by-inspection. Run the suite with
 * MATHILDA_PACK_DIAG=1: every lossy touch of an int64 buffer prints and, with
 * MATHILDA_PACK_DIAG=abort, aborts so a debugger shows the caller. Either way
 * `nd_int64_lossy_hit` is a stable breakpoint target -- `break
 * nd_int64_lossy_hit` catches them without any environment variable at all.
 *
 * Costs nothing for other dtypes: the check lives INSIDE the existing
 * `case NDT_INT64` arm, which the switch was already going to select.
 * -------------------------------------------------------------------------- */
unsigned long nd_int64_lossy_count = 0;

void nd_int64_lossy_hit(const char* what) {
    static int mode = -1;      /* -1 unread, 0 off, 1 warn, 2 abort */
    nd_int64_lossy_count++;
    if (mode < 0) {
        const char* v = getenv("MATHILDA_PACK_DIAG");
        mode = !v ? 0 : (strcmp(v, "abort") == 0 ? 2 : 1);
    }
    if (mode == 0) return;
    fprintf(stderr, "PACKDIAG: lossy %s on an NDT_INT64 buffer "
                    "(exact only to 2^53) -- hit #%lu\n", what, nd_int64_lossy_count);
    if (mode == 2) abort();
}

void ndt_get(const void* buf, size_t k, NDType dt, double* re, double* im) {
    switch (dt) {
        case NDT_FLOAT64:
            *re = ((const double*)buf)[k]; *im = 0.0; break;
        case NDT_FLOAT32:
            *re = (double)((const float*)buf)[k]; *im = 0.0; break;
        case NDT_COMPLEX64:
            *re = ((const double*)buf)[2 * k];
            *im = ((const double*)buf)[2 * k + 1]; break;
        case NDT_COMPLEX32:
            *re = (double)((const float*)buf)[2 * k];
            *im = (double)((const float*)buf)[2 * k + 1]; break;
        case NDT_INT64:
            /* Exact only to 2^53. Every read of an int64 array that must be
             * exact goes through ndt_get_i instead; this case exists so the
             * generic dtype-agnostic consumers stay total. */
            nd_int64_lossy_hit("read");
            *re = (double)((const int64_t*)buf)[k]; *im = 0.0; break;
        case NDT_BOOL:
            *re = ((const uint8_t*)buf)[k] ? 1.0 : 0.0; *im = 0.0; break;
        default:
            /* Unreachable for a valid dtype, but makes the switch total so the
             * compiler can prove re and im are always written — otherwise every
             * caller reading re and im warns under -O3 -Wmaybe-uninitialized. */
            *re = 0.0; *im = 0.0; break;
    }
}

void ndt_set(void* buf, size_t k, NDType dt, double re, double im) {
    switch (dt) {
        case NDT_FLOAT64:
            ((double*)buf)[k] = re; break;
        case NDT_FLOAT32:
            ((float*)buf)[k] = (float)re; break;
        case NDT_COMPLEX64:
            ((double*)buf)[2 * k] = re;
            ((double*)buf)[2 * k + 1] = im; break;
        case NDT_COMPLEX32:
            ((float*)buf)[2 * k] = (float)re;
            ((float*)buf)[2 * k + 1] = (float)im; break;
        case NDT_INT64:
            nd_int64_lossy_hit("write");
            ((int64_t*)buf)[k] = (int64_t)re; break;      /* see ndt_set_i */
        case NDT_BOOL:
            ((uint8_t*)buf)[k] = (re != 0.0) ? 1 : 0; break;
    }
}

/* EXACT int64 element access, the pair ndt_get/ndt_set cannot provide because
 * they route through `double` and so lose the top 11 bits of an int64. Anything
 * that reads or writes an NDT_INT64 buffer for real — the Compile[] boundary,
 * the VM's array opcodes — uses these. */
int64_t ndt_get_i(const void* buf, size_t k, NDType dt) {
    if (dt == NDT_INT64) return ((const int64_t*)buf)[k];
    if (dt == NDT_BOOL)  return ((const uint8_t*)buf)[k] ? 1 : 0;
    double re, im; ndt_get(buf, k, dt, &re, &im); (void)im;
    return (int64_t)re;
}

void ndt_set_i(void* buf, size_t k, NDType dt, int64_t v) {
    if (dt == NDT_INT64) { ((int64_t*)buf)[k] = v; return; }
    if (dt == NDT_BOOL)  { ((uint8_t*)buf)[k] = v ? 1 : 0; return; }
    ndt_set(buf, k, dt, (double)v, 0.0);
}

bool ndt_from_string(const char* s, NDType* out) {
    if (!s) return false;
    if (strcmp(s, "float64") == 0) { *out = NDT_FLOAT64; return true; }
    if (strcmp(s, "float32") == 0) { *out = NDT_FLOAT32; return true; }
    if (strcmp(s, "complex64") == 0) { *out = NDT_COMPLEX64; return true; }
    if (strcmp(s, "complex32") == 0) { *out = NDT_COMPLEX32; return true; }
    /* "int64" is accepted because ndt_to_string already EMITS it, so
     * DataType[x] must round-trip back through DataType -> ... . Packed lists
     * infer this dtype from an all-Integer list anyway (see pack.c), so the
     * option adds no reach the system did not already have. */
    if (strcmp(s, "int64") == 0) { *out = NDT_INT64; return true; }
    /* One byte per element, materialising as True/False. Both the lowercase
     * "bool" (which ndt_to_string emits, so DataType[x] round-trips) and the
     * WL-style "Boolean" are accepted. */
    if (strcmp(s, "bool") == 0 || strcmp(s, "Boolean") == 0) {
        *out = NDT_BOOL; return true;
    }
    return false;
}

const char* ndt_to_string(NDType dt) {
    switch (dt) {
        case NDT_FLOAT64:   return "float64";
        case NDT_FLOAT32:   return "float32";
        case NDT_COMPLEX64: return "complex64";
        case NDT_COMPLEX32: return "complex32";
        case NDT_INT64:     return "int64";
        case NDT_BOOL:      return "bool";
    }
    return "float64";
}

NDType ndt_promote(NDType a, NDType b) {
    /* bool is weaker than int64 and, like numpy, promotes to int64 the moment it
     * meets arithmetic. This is defensive only: a bool operand makes an arithmetic
     * head decline its fast path and delist to the symbolic List (True + True is
     * unevaluated in Mathematica), so bool rarely reaches here. Folding it to int64
     * keeps ndt_promote total and matches numpy for the paths that do (e.g. a bool
     * array meeting a float array through a listable numeric kernel). */
    if (a == NDT_BOOL) a = NDT_INT64;
    if (b == NDT_BOOL) b = NDT_INT64;
    /* int64 is the weakest numeric dtype: meeting any float type yields that float
     * type, and int64 with itself stays exact. */
    if (a == NDT_INT64 && b == NDT_INT64) return NDT_INT64;
    if (a == NDT_INT64) return b;
    if (b == NDT_INT64) return a;
    bool cplx = ndt_is_complex(a) || ndt_is_complex(b);
    bool wide = ndt_comp_size(a) == sizeof(double) ||
                ndt_comp_size(b) == sizeof(double);
    if (cplx) return wide ? NDT_COMPLEX64 : NDT_COMPLEX32;
    return wide ? NDT_FLOAT64 : NDT_FLOAT32;
}

/* Move a dtype onto the complex axis, preserving its component width: float32 ->
 * complex32, float64 -> complex64 (complex dtypes unchanged). Used when a
 * complex scalar or a real->complex escape must not widen the array's precision
 * (numpy value-based casting: a float32 array meeting a complex value stays at
 * float32 components, i.e. our complex32). */
NDType ndt_as_complex(NDType dt) {
    if (dt == NDT_INT64 || dt == NDT_BOOL) return NDT_COMPLEX64;
    return (ndt_comp_size(dt) == sizeof(double)) ? NDT_COMPLEX64 : NDT_COMPLEX32;
}

Expr* ndarray_buffer_element_to_expr(const void* buf, size_t k, NDType dt) {
    /* THE single decision of what head an element of dtype `dt` materialises
     * as. Every unpack path in the tree routes through here so none of them can
     * drift: an NDT_INT64 element is an exact Integer read through the exact
     * accessor (ndt_get/ndt_set are exact only to 2^53), everything else is a
     * Real or a Complex[Real, Real]. */
    if (dt == NDT_INT64)
        return expr_new_integer(ndt_get_i(buf, k, dt));
    if (dt == NDT_BOOL)
        return expr_new_symbol(((const uint8_t*)buf)[k] ? SYM_True : SYM_False);
    double re, im;
    ndt_get(buf, k, dt, &re, &im);
    if (ndt_is_complex(dt)) {
        Expr* args[2] = { expr_new_real(re), expr_new_real(im) };
        return expr_new_function(expr_new_symbol(SYM_Complex), args, 2);
    }
    return expr_new_real(re);
}

Expr* ndarray_element_to_expr(const Expr* a, size_t k) {
    return ndarray_buffer_element_to_expr(a->data.ndarray.data, k,
                                          a->data.ndarray.dtype);
}

/* Integer or Real scalar → double. Used for the components of a Complex[] leaf
 * as well as bare real leaves. */
static bool leaf_real(const Expr* e, double* out) {
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL) { *out = e->data.real; return true; }
    return false;
}

/* Machine-precision-packable leaf for dtype `dt`, decoded into a (re, im) pair.
 * Real dtypes accept Integer/Real only. Complex dtypes additionally accept a
 * Complex[re, im] literal (and a bare real, im=0). BigInt/Rational/MPFR/symbolic
 * all fail packing so the caller falls back to a nested List. */
/* An element of an NDT_INT64 buffer: an exact machine Integer and nothing else.
 *
 * Not `leaf_to_component`, which routes through a double — that would round
 * 9007199254740993 on the way in, and would silently accept 1.5 as an element of
 * an integer array.  A bigint leaf does not fit and is refused here, so the call
 * falls back to the interpreter, which has GMP. */
static bool leaf_to_int64(const Expr* e, int64_t* out) {
    if (!e || e->type != EXPR_INTEGER) return false;
    *out = e->data.integer;
    return true;
}

static bool leaf_to_component(const Expr* e, NDType dt, double* re, double* im) {
    if (dt == NDT_BOOL) {
        /* A bool buffer accepts the two symbols True/False and nothing else, so
         * `NDArray[{1, 0}, DataType -> "bool"]` declines to the interpreter rather
         * than silently reinterpreting integers as truth values. */
        if (e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_True)  { *re = 1.0; *im = 0.0; return true; }
        if (e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_False) { *re = 0.0; *im = 0.0; return true; }
        return false;
    }
    if (dt == NDT_INT64) {
        int64_t v;
        if (!leaf_to_int64(e, &v)) return false;
        *re = (double)v; *im = 0.0;      /* probe_dims only asks "is this legal" */
        return true;
    }
    *im = 0.0;
    if (leaf_real(e, re)) return true;
    if (ndt_is_complex(dt) && head_is(e, SYM_Complex) &&
        e->data.function.arg_count == 2) {
        return leaf_real(e->data.function.args[0], re) &&
               leaf_real(e->data.function.args[1], im);
    }
    return false;
}

/* Recursive rank/shape probe, mirroring linalg/util.c's get_tensor_dims but
 * over `int64_t dims[NDARRAY_MAX_RANK]` and rejecting leaves not packable at
 * dtype `dt` up front. Returns 0 for a non-List leaf, -1 for jagged/non-numeric,
 * else the rank. */
static int probe_dims(const Expr* e, int64_t* dims, NDType dt) {
    if (!head_is(e, SYM_List)) return 0;
    int64_t len = (int64_t)e->data.function.arg_count;
    dims[0] = len;
    if (len == 0) return -1; /* empty: nothing to pack, degrade to List */

    int sub_rank = probe_dims(e->data.function.args[0], dims + 1, dt);
    if (sub_rank == 0) {
        double r, im;
        if (!leaf_to_component(e->data.function.args[0], dt, &r, &im)) return -1;
    } else if (sub_rank < 0) {
        return -1;
    }
    for (int64_t i = 1; i < len; i++) {
        int64_t cur[NDARRAY_MAX_RANK];
        int cur_rank = probe_dims(e->data.function.args[i], cur, dt);
        if (cur_rank != sub_rank) return -1;
        if (cur_rank == 0) {
            double r, im;
            if (!leaf_to_component(e->data.function.args[i], dt, &r, &im)) return -1;
        }
        for (int j = 0; j < sub_rank; j++) {
            if (cur[j] != dims[j + 1]) return -1;
        }
    }
    return sub_rank + 1;
}

static void flatten_into(const Expr* e, void* flat, size_t* idx, NDType dt) {
    if (head_is(e, SYM_List)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            flatten_into(e->data.function.args[i], flat, idx, dt);
        }
    } else if (dt == NDT_INT64) {
        int64_t v = 0;
        leaf_to_int64(e, &v);               /* already validated by probe_dims */
        ndt_set_i(flat, (*idx)++, dt, v);   /* EXACT: ndt_set would round */
    } else {
        double re = 0.0, im = 0.0;
        leaf_to_component(e, dt, &re, &im); /* already validated by probe_dims */
        ndt_set(flat, (*idx)++, dt, re, im);
    }
}

Expr* ndarray_from_nested_list(const Expr* list, NDType dtype) {
    if (!head_is(list, SYM_List)) return NULL;
    int64_t dims[NDARRAY_MAX_RANK];
    int rank = probe_dims(list, dims, dtype);
    if (rank <= 0 || rank >= NDARRAY_MAX_RANK) return NULL;

    size_t n = 1;
    for (int i = 0; i < rank; i++) n *= (size_t)dims[i];
    void* data = malloc(ndt_elem_size(dtype) * n);
    if (!data) return NULL;
    size_t idx = 0;
    flatten_into(list, data, &idx, dtype);

    return expr_new_ndarray_raw(rank, dims, data, dtype); /* takes ownership of data */
}

/* Recursive rebuild of one axis-level of the nested List tree, mirroring
 * linalg/dot.c's build_tensor but reading from a flat double* rather than
 * an Expr** buffer. `idx` is the running cursor into `data`. */
static Expr* rebuild_level(const int64_t* dims, int rank, int level,
                            const void* data, NDType dt, size_t* idx) {
    if (level == rank) {
        /* The element HEAD is decided in exactly one place so no rebuild path can
         * drift: an int64 buffer rebuilds as exact Integers (through the EXACT
         * accessor, since ndt_get loses anything above 2^53), a bool buffer as
         * True/False, everything else as a Real or Complex[re, im]. This is the
         * path a Compile[] array result comes back on. */
        return ndarray_buffer_element_to_expr(data, (*idx)++, dt);
    }
    int64_t len = dims[level];
    Expr** args = malloc(sizeof(Expr*) * (size_t)len);
    for (int64_t i = 0; i < len; i++) {
        args[i] = rebuild_level(dims, rank, level + 1, data, dt, idx);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), args, (size_t)len);
    free(args);
    return result;
}

Expr* ndarray_from_nested_list_like(const Expr* src, const Expr* list, NDType dtype) {
    Expr* out = ndarray_from_nested_list(list, dtype);
    if (out && src && is_ndarray(src))
        out->data.ndarray.present_as = src->data.ndarray.present_as;
    return out;
}

Expr* ndarray_to_nested_list(const Expr* a) {
    size_t idx = 0;
    return rebuild_level(a->data.ndarray.dims, a->data.ndarray.rank, 0,
                          a->data.ndarray.data, a->data.ndarray.dtype, &idx);
}

Expr* ndarray_unpack_top_level(const Expr* a) {
    if (!is_ndarray(a) || a->data.ndarray.rank < 1) return NULL;
    int64_t n = a->data.ndarray.dims[0];
    Expr** rows = malloc(sizeof(Expr*) * (n > 0 ? (size_t)n : 1));
    if (!rows) return NULL;
    for (int64_t i = 0; i < n; i++) {
        /* a[[i+1]]: 1-based; a scalar leaf for rank 1, else a packed sub-array
         * inheriting `a`'s presentation (ndarray_part builds through
         * expr_new_ndarray_like). No non-integer subscript, so *degrade is only
         * ever set on an internal failure -- treated the same as NULL below. */
        Expr* idx = expr_new_integer(i + 1);
        bool degrade = false;
        Expr* row = ndarray_part(a, &idx, 1, &degrade);
        expr_free(idx);
        if (!row) {
            for (int64_t j = 0; j < i; j++) expr_free(rows[j]);
            free(rows);
            return NULL;
        }
        rows[i] = row;
    }
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)(n > 0 ? n : 0));
    free(rows);
    return list;
}

/* Per-axis index selection for ndarray_part's general (sliced) path. `keep` is
 * false for an integer subscript (the axis is dropped from the result) and true
 * for All / Span / a List of positions (the axis survives with `n` selected
 * source positions listed 0-based in `pos`). `pos` is always heap-allocated
 * (length 1 for a dropped axis) so cleanup is uniform. */
typedef struct { int64_t* pos; int64_t n; bool keep; } NDAxisSel;

/* Status codes for build_axis_selector: OK, or leave the whole Part alone. */
#define NDPART_OK       0
#define NDPART_DEGRADE  1   /* spec is valid but not natively handled -> delist */
#define NDPART_INVALID  2   /* out-of-range subscript -> leave Part unevaluated */

/* Resolve a 1-based (possibly negative) subscript against axis length `len` to
 * a 0-based position, or return false if out of range. */
static bool nd_resolve_index(int64_t k, int64_t len, int64_t* pos0) {
    if (k < 0) k = len + k + 1;
    if (k < 1 || k > len) return false;
    *pos0 = k - 1;
    return true;
}

/* Fill `sel` for one axis of length `len` from `spec` (NULL means an implicit
 * trailing All). Mirrors the List-path Span/All/List semantics in part.c. */
static int build_axis_selector(const Expr* spec, int64_t len, NDAxisSel* sel) {
    sel->pos = NULL; sel->n = 0; sel->keep = false;

    /* Implicit trailing All, or an explicit All symbol: keep every position. */
    if (spec == NULL ||
        (spec->type == EXPR_SYMBOL && spec->data.symbol.name == SYM_All)) {
        sel->pos = malloc(sizeof(int64_t) * (size_t)(len > 0 ? len : 1));
        if (!sel->pos) return NDPART_DEGRADE;
        for (int64_t i = 0; i < len; i++) sel->pos[i] = i;
        sel->n = len; sel->keep = true;
        return NDPART_OK;
    }

    /* Integer subscript: drop the axis, fix a single position. */
    if (spec->type == EXPR_INTEGER) {
        int64_t p;
        if (!nd_resolve_index(spec->data.integer, len, &p)) return NDPART_INVALID;
        sel->pos = malloc(sizeof(int64_t));
        if (!sel->pos) return NDPART_DEGRADE;
        sel->pos[0] = p; sel->n = 1; sel->keep = false;
        return NDPART_OK;
    }

    if (spec->type == EXPR_FUNCTION && spec->data.function.head->type == EXPR_SYMBOL) {
        const char* h = spec->data.function.head->data.symbol.name;

        /* Span[start, end, step] with the same start/end/step resolution the
         * List path uses (Integer / All / UpTo per component). */
        if (h == SYM_Span) {
            int64_t start = 1, end = len, step = 1;
            size_t argc = spec->data.function.arg_count;
            const Expr* const* sa = (const Expr* const*)spec->data.function.args;
            for (size_t c = 0; c < argc && c < 3; c++) {
                const Expr* e = sa[c];
                bool is_all = (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_All);
                bool is_upto = (e->type == EXPR_FUNCTION &&
                                e->data.function.head->type == EXPR_SYMBOL &&
                                e->data.function.head->data.symbol.name == SYM_UpTo &&
                                e->data.function.arg_count == 1 &&
                                e->data.function.args[0]->type == EXPR_INTEGER);
                int64_t v = 0;
                if (e->type == EXPR_INTEGER) v = e->data.integer;
                else if (is_upto) v = e->data.function.args[0]->data.integer;
                else if (!is_all) return NDPART_DEGRADE;

                if (c == 0) {          /* start */
                    if (is_all) start = 1;
                    else { if (is_upto && v > len) v = len; start = v < 0 ? len + v + 1 : v; }
                } else if (c == 1) {   /* end */
                    if (is_all) end = len;
                    else { if (is_upto && v > len) v = len; end = v < 0 ? len + v + 1 : v; }
                } else {               /* step */
                    if (is_all) step = 1;
                    else if (is_upto) return NDPART_DEGRADE;
                    else step = v;
                    if (step == 0) return NDPART_INVALID;
                }
            }
            int64_t count = 0;
            if (step > 0) { if (start <= end && start >= 1 && end <= len) count = (end - start) / step + 1; }
            else          { if (start >= end && start <= len && end >= 1) count = (start - end) / (-step) + 1; }
            if (count <= 0) return NDPART_DEGRADE;   /* empty span -> exact List {} */
            sel->pos = malloc(sizeof(int64_t) * (size_t)count);
            if (!sel->pos) return NDPART_DEGRADE;
            int64_t cur = start;
            for (int64_t i = 0; i < count; i++) { sel->pos[i] = cur - 1; cur += step; }
            sel->n = count; sel->keep = true;
            return NDPART_OK;
        }

        /* A List of plain integer positions (fancy gather). Any non-integer
         * entry, or an empty list, degrades to the exact List-path result. */
        if (h == SYM_List) {
            size_t m = spec->data.function.arg_count;
            if (m == 0) return NDPART_DEGRADE;
            sel->pos = malloc(sizeof(int64_t) * m);
            if (!sel->pos) return NDPART_DEGRADE;
            for (size_t i = 0; i < m; i++) {
                const Expr* e = spec->data.function.args[i];
                int64_t p;
                if (e->type != EXPR_INTEGER || !nd_resolve_index(e->data.integer, len, &p)) {
                    free(sel->pos); sel->pos = NULL;
                    return (e->type == EXPR_INTEGER) ? NDPART_INVALID : NDPART_DEGRADE;
                }
                sel->pos[i] = p;
            }
            sel->n = (int64_t)m; sel->keep = true;
            return NDPART_OK;
        }
    }

    /* A PACKED list of positions -- the same fancy gather as the List arm
     * above, arriving in the representation the rest of the system now hands
     * back. Flatten, Range, RandomInteger, Ordering and every int64 kernel
     * produce a packed list, so without this arm `x[[idx]]` -- the operation
     * every sparse-matrix and graph kernel is built out of -- degraded the
     * whole Part and materialised BOTH arrays: a 1.6e6-index gather out of a
     * 100000-element vector cost 389 ms against 15.7 ms once the index buffer
     * is read directly, and a PageRank built on it 14.1 s against 486 ms.
     *
     * is_packed_list, not is_ndarray: a VISIBLE NDArray[...] is not a List and
     * the List path would not have accepted one, so the two surfaces must not
     * start disagreeing here. Rank 1 only, and int64 only -- a Real position is
     * Part::pkspec1 in Mathematica, and the List arm above rejects it too. */
    if (is_packed_list(spec) &&
        spec->data.ndarray.rank == 1 &&
        spec->data.ndarray.dtype == NDT_INT64) {
        int64_t m = spec->data.ndarray.dims[0];
        if (m <= 0) return NDPART_DEGRADE;          /* {} -> exact List path */
        sel->pos = malloc(sizeof(int64_t) * (size_t)m);
        if (!sel->pos) return NDPART_DEGRADE;
        const int64_t* src = (const int64_t*)spec->data.ndarray.data;
        for (int64_t i = 0; i < m; i++) {
            if (!nd_resolve_index(src[i], len, &sel->pos[i])) {
                free(sel->pos); sel->pos = NULL;
                return NDPART_INVALID;
            }
        }
        sel->n = m; sel->keep = true;
        return NDPART_OK;
    }

    return NDPART_DEGRADE;   /* unrecognised spec (Key, pattern, ...) */
}

Expr* ndarray_part(const Expr* a, Expr** indices, size_t nindices, bool* degrade) {
    *degrade = false;
    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;
    NDType dt = a->data.ndarray.dtype;

    if (nindices == 0) return NULL;               /* caller supplies >= 1 index */
    if (nindices > (size_t)rank) return NULL;     /* too many subscripts: partw */

    /* Fast path: every subscript is a plain integer. Fixing the leading axes of
     * a row-major buffer selects a contiguous block, so this is O(1) for a
     * scalar leaf and a single memcpy for a sub-array. */
    bool all_int = true;
    for (size_t i = 0; i < nindices; i++)
        if (indices[i]->type != EXPR_INTEGER) { all_int = false; break; }

    if (all_int) {
        size_t offset = 0;
        for (size_t i = 0; i < nindices; i++) {
            int64_t p;
            if (!nd_resolve_index(indices[i]->data.integer, dims[i], &p)) return NULL;
            int64_t stride = 1;
            for (int j = (int)i + 1; j < rank; j++) stride *= dims[j];
            offset += (size_t)p * (size_t)stride;
        }
        if ((size_t)rank == nindices)             /* fully indexed -> scalar leaf */
            return ndarray_buffer_element_to_expr(a->data.ndarray.data, offset, dt);
        int subrank = rank - (int)nindices;
        int64_t subdims[NDARRAY_MAX_RANK];
        size_t subsize = 1;
        for (int j = 0; j < subrank; j++) {
            subdims[j] = dims[nindices + (size_t)j];
            subsize *= (size_t)subdims[j];
        }
        size_t esz = ndt_elem_size(dt);
        void* out = malloc(esz * subsize);
        if (!out) return NULL;
        memcpy(out, (const char*)a->data.ndarray.data + offset * esz, esz * subsize);
        return expr_new_ndarray_like(a, subrank, subdims, out, dt); /* takes ownership */
    }

    /* General path: a mix of Integer / All / Span / List subscripts (plus
     * implicit trailing All axes). Build a per-axis list of source positions,
     * then gather into a packed result via a mixed-radix walk over all axes
     * (dropped integer axes have radix 1, so they contribute a fixed digit and
     * do not appear in the output shape). */
    NDAxisSel sel[NDARRAY_MAX_RANK];
    for (int i = 0; i < rank; i++) { sel[i].pos = NULL; sel[i].n = 0; sel[i].keep = false; }

    int status = NDPART_OK;
    for (int i = 0; i < rank; i++) {
        const Expr* spec = ((size_t)i < nindices) ? indices[i] : NULL;
        status = build_axis_selector(spec, dims[i], &sel[i]);
        if (status != NDPART_OK) break;
    }
    if (status != NDPART_OK) {
        for (int i = 0; i < rank; i++) free(sel[i].pos);
        if (status == NDPART_DEGRADE) *degrade = true;
        return NULL;   /* INVALID -> leave unevaluated; DEGRADE -> caller delists */
    }

    int subrank = 0;
    int64_t subdims[NDARRAY_MAX_RANK];
    size_t total = 1;
    for (int i = 0; i < rank; i++) {
        if (sel[i].keep) subdims[subrank++] = sel[i].n;
        total *= (size_t)sel[i].n;
    }

    int64_t stride[NDARRAY_MAX_RANK];
    { int64_t s = 1; for (int i = rank - 1; i >= 0; i--) { stride[i] = s; s *= dims[i]; } }

    size_t esz = ndt_elem_size(dt);
    void* out = malloc(esz * (total ? total : 1));
    if (!out) { for (int i = 0; i < rank; i++) free(sel[i].pos); return NULL; }

    /* A selector that is a CONTIGUOUS unit-step run on every axis is a block
     * copy, not a gather. The general walk below costs a modulo and a division
     * PER AXIS PER ELEMENT plus an indirect position lookup, which is what made
     * v[[2 ;; -1]] cost 14 ms where the identical Drop[v, 1] cost 0.9 --
     * a plain memcpy of the same bytes. Spans are how array code names a
     * sub-range, so this is the common shape, not a corner. */
    bool contig = (total > 0);
    for (int i = 0; i < rank && contig; i++)
        for (int64_t k = 1; k < sel[i].n; k++)
            if (sel[i].pos[k] != sel[i].pos[k - 1] + 1) { contig = false; break; }
    if (contig) {
        /* Innermost axis supplies the contiguous run; the outer axes index which
         * run, so there are total/run copies, not total. */
        size_t run = (size_t)sel[rank - 1].n;
        size_t nrows = total / run;
        for (size_t ri = 0; ri < nrows; ri++) {
            size_t rem = ri, src = (size_t)sel[rank - 1].pos[0] * (size_t)stride[rank - 1];
            for (int i = rank - 2; i >= 0; i--) {
                int64_t d = sel[i].n;
                src += (size_t)sel[i].pos[rem % (size_t)d] * (size_t)stride[i];
                rem /= (size_t)d;
            }
            memcpy((char*)out + ri * run * esz,
                   (const char*)a->data.ndarray.data + src * esz, run * esz);
        }
        for (int i = 0; i < rank; i++) free(sel[i].pos);
        return expr_new_ndarray_like(a, subrank, subdims, out, dt);
    }

    for (size_t oi = 0; oi < total; oi++) {
        size_t rem = oi, src = 0;
        for (int i = rank - 1; i >= 0; i--) {
            int64_t d = sel[i].n;                 /* radix (>=1) */
            int64_t digit = (int64_t)(rem % (size_t)d);
            rem /= (size_t)d;
            src += (size_t)sel[i].pos[digit] * (size_t)stride[i];
        }
        /* Source and destination share `dt`, so this is a copy, not a
         * conversion -- move the bytes rather than round-tripping through the
         * (double, double) pair, which is both slower and inexact above 2^53
         * for an int64 buffer. */
        memcpy((char*)out + oi * esz,
               (const char*)a->data.ndarray.data + src * esz, esz);
    }

    for (int i = 0; i < rank; i++) free(sel[i].pos);

    if (subrank == 0) {                           /* all axes dropped -> scalar */
        Expr* leaf = ndarray_buffer_element_to_expr(out, 0, dt);
        free(out);
        return leaf;
    }
    return expr_new_ndarray_like(a, subrank, subdims, out, dt); /* takes ownership */
}

/* ---------------------------------------------------------------------------
 * Part assignment, IN PLACE.
 *
 * `a[[spec...]] = rhs` writes the selected positions of `a`'s buffer directly
 * instead of rebuilding the array. Mutating in place is only sound because the
 * one caller — compiled bytecode — emits this exclusively for arrays the
 * running program OWNS (a Module local or a temporary). An argument array is
 * borrowed and is rejected at emit time.
 *
 * The positions come from the SAME build_axis_selector and mixed-radix walk
 * that ndarray_part gathers through, so a read and a write of one spec can
 * never disagree about which elements that spec names.
 *
 * `rhs` is either a single number, broadcast to every selected position, or an
 * NDArray whose element count matches the selection (WL assigns a matching
 * sub-shape element-wise). Returns false for anything else — an unsupported
 * spec, an out-of-range subscript, a count mismatch, or a complex value into a
 * real buffer — having written nothing, so the compiled call fails cleanly and
 * the interpreter re-runs the whole body.
 * ------------------------------------------------------------------------ */
bool ndarray_part_set(Expr* a, Expr* const* indices, size_t nindices,
                      const Expr* rhs) {
    if (!a || a->type != EXPR_NDARRAY || !rhs) return false;
    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;
    NDType dt = a->data.ndarray.dtype;
    if (nindices == 0 || nindices > (size_t)rank) return false;

    NDAxisSel sel[NDARRAY_MAX_RANK];
    for (int i = 0; i < rank; i++) { sel[i].pos = NULL; sel[i].n = 0; sel[i].keep = false; }

    int status = NDPART_OK;
    for (int i = 0; i < rank; i++) {
        const Expr* spec = ((size_t)i < nindices) ? indices[i] : NULL;
        status = build_axis_selector(spec, dims[i], &sel[i]);
        if (status != NDPART_OK) break;
    }
    if (status != NDPART_OK) {
        for (int i = 0; i < rank; i++) free(sel[i].pos);
        return false;
    }

    size_t total = 1;
    for (int i = 0; i < rank; i++) total *= (size_t)sel[i].n;

    /* Validate the whole right-hand side BEFORE the first store, so a rejected
     * assignment leaves the buffer exactly as it was. */
    bool rhs_is_arr = (rhs->type == EXPR_NDARRAY);
    double sre = 0.0, sim = 0.0;
    if (rhs_is_arr) {
        size_t rn = ndarray_size(rhs);
        if (rn != total) { for (int i = 0; i < rank; i++) free(sel[i].pos); return false; }
        if (ndt_is_complex(rhs->data.ndarray.dtype) && !ndt_is_complex(dt)) {
            for (int i = 0; i < rank; i++) free(sel[i].pos);
            return false;
        }
    } else if (!leaf_to_component(rhs, dt, &sre, &sim)) {
        /* Same packability rule a List uses to become a buffer of this dtype,
         * so it also rejects a complex value into a real array. */
        for (int i = 0; i < rank; i++) free(sel[i].pos);
        return false;
    }

    int64_t stride[NDARRAY_MAX_RANK];
    { int64_t s = 1; for (int i = rank - 1; i >= 0; i--) { stride[i] = s; s *= dims[i]; } }

    for (size_t oi = 0; oi < total; oi++) {
        size_t rem = oi, dstpos = 0;
        for (int i = rank - 1; i >= 0; i--) {
            int64_t d = sel[i].n;                 /* radix (>=1) */
            int64_t digit = (int64_t)(rem % (size_t)d);
            rem /= (size_t)d;
            dstpos += (size_t)sel[i].pos[digit] * (size_t)stride[i];
        }
        double re = sre, im = sim;
        if (rhs_is_arr) ndt_get(rhs->data.ndarray.data, oi, rhs->data.ndarray.dtype, &re, &im);
        ndt_set(a->data.ndarray.data, dstpos, dt, re, im);
    }

    for (int i = 0; i < rank; i++) free(sel[i].pos);
    return true;
}

Expr* ndarray_dot2(const Expr* a, const Expr* b, bool* shape_error) {
    int rankA = a->data.ndarray.rank;
    int rankB = b->data.ndarray.rank;
    if (rankA < 1 || rankA > 2 || rankB < 1 || rankB > 2) return NULL;

    const int64_t* dimsA = a->data.ndarray.dims;
    const int64_t* dimsB = b->data.ndarray.dims;
    int64_t K = dimsA[rankA - 1];
    if (K != dimsB[0]) { *shape_error = true; return NULL; }

    int64_t Ra = (rankA == 2) ? dimsA[0] : 1;
    int64_t Sb = (rankB == 2) ? dimsB[1] : 1;

    NDType dta = a->data.ndarray.dtype, dtb = b->data.ndarray.dtype;
    if (dta == NDT_BOOL || dtb == NDT_BOOL) return NULL;  /* not numeric; see ndarray_map_unary */
    NDType dtc = ndt_promote(dta, dtb);
    const void* pa = a->data.ndarray.data;
    const void* pb = b->data.ndarray.data;
    size_t count = (size_t)(Ra * Sb);
    void* out = malloc(ndt_elem_size(dtc) * count);

    if (dta == NDT_INT64 || dtb == NDT_INT64) {
        /* Range[10^6] . Range[10^6] is about 3.3e17 -- well past 2^53, so the
         * double path would answer with a rounded Real. Accumulate exactly, and
         * abandon on the first overflow so the List path gives the bigint. Mixed
         * int/float is left to the promoted double path below, where inexactness
         * is already the right answer. */
        if (dta == NDT_INT64 && dtb == NDT_INT64) {
            int64_t* od = (int64_t*)out;
            for (int64_t r = 0; r < Ra; r++)
                for (int64_t s = 0; s < Sb; s++) {
                    int64_t sum = 0;
                    for (int64_t k = 0; k < K; k++) {
                        int64_t t;
                        if (ci_mul_i64(((const int64_t*)pa)[r * K + k],
                                       ((const int64_t*)pb)[k * Sb + s], &t) ||
                            ci_add_i64(sum, t, &sum)) { free(out); return NULL; }
                    }
                    od[r * Sb + s] = sum;
                }
            goto dot_done;
        }
    }
    if (dta == NDT_FLOAT64 && dtb == NDT_FLOAT64) {
        /* Preserved hot path: real double triple loop, no per-element widening. */
        double* od = (double*)out;
        for (int64_t r = 0; r < Ra; r++)
            for (int64_t s = 0; s < Sb; s++) {
                double sum = 0.0;
                for (int64_t k = 0; k < K; k++)
                    sum += ((const double*)pa)[r * K + k] *
                           ((const double*)pb)[k * Sb + s];
                od[r * Sb + s] = sum;
            }
    } else {
        for (int64_t r = 0; r < Ra; r++)
            for (int64_t s = 0; s < Sb; s++) {
                double sre = 0.0, sim = 0.0;
                for (int64_t k = 0; k < K; k++) {
                    double ar, ai, br, bi;
                    ndt_get(pa, (size_t)(r * K + k), dta, &ar, &ai);
                    ndt_get(pb, (size_t)(k * Sb + s), dtb, &br, &bi);
                    sre += ar * br - ai * bi;
                    sim += ar * bi + ai * br;
                }
                ndt_set(out, (size_t)(r * Sb + s), dtc, sre, sim);
            }
    }

dot_done:;
    int rankC = rankA + rankB - 2;
    if (rankC == 0) {
        Expr* leaf = ndarray_buffer_element_to_expr(out, 0, dtc);
        free(out);
        return leaf;
    }
    int64_t dimsC[2];
    if (rankC == 1) dimsC[0] = (rankA == 2) ? Ra : Sb;
    else { dimsC[0] = Ra; dimsC[1] = Sb; }
    return expr_new_ndarray_like(nd_present_src2(a, b), rankC, dimsC, out, dtc); /* takes ownership of out */
}

static bool same_shape(const Expr* a, const Expr* b) {
    if (a->data.ndarray.rank != b->data.ndarray.rank) return false;
    for (int i = 0; i < a->data.ndarray.rank; i++) {
        if (a->data.ndarray.dims[i] != b->data.ndarray.dims[i]) return false;
    }
    return true;
}

/* Decode a numeric scalar operand (Integer / Real / Complex[re,im]) into a
 * (re, im) pair. Returns false for a symbolic or non-machine scalar, so a Plus /
 * Times containing such a term falls through to the generic symbolic path. */
static bool scalar_value(const Expr* e, double* re, double* im, bool* is_cplx) {
    int64_t rn, rd;
    *im = 0.0; *is_cplx = false;
    if (leaf_real(e, re)) return true;
    /* An NDArray is a machine-float buffer, so a Rational / BigInt / MPFR scalar
     * broadcasts at its machine value (no exactness to preserve). This is what
     * folds 1/2 NDArray[...] (= Times[1/2, NDArray]) and Mod[NDArray, 1/2]
     * instead of leaving them symbolic. Note leaf_real (used for List->buffer
     * packing) deliberately stays Integer/Real-only. */
    if (e->type == EXPR_BIGINT) { *re = mpz_get_d(e->data.bigint); return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)   { *re = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true; }
#endif
    if (is_rational(e, &rn, &rd)) { *re = (double)rn / (double)rd; return true; }
    if (head_is(e, SYM_Complex) && e->data.function.arg_count == 2 &&
        leaf_real(e->data.function.args[0], re) &&
        leaf_real(e->data.function.args[1], im)) {
        *is_cplx = true;
        return true;
    }
    return false;
}

/* Elementwise Plus (is_plus) / Times over a mix of NDArray and numeric-scalar
 * operands. numpy-style scalar broadcasting: a scalar operand (Integer / Real /
 * Complex) applies to every element; at least one operand must be an NDArray and
 * all NDArray operands must share a shape. A symbolic operand, or no NDArray at
 * all, returns NULL so the caller uses the generic symbolic path. Scalars are
 * "weak": they set the complex axis (a Complex scalar promotes the result to
 * complex) but do not widen the float precision, which is taken from the NDArray
 * operands (float32 array + 1 stays float32; float32 array + I -> complex32). */
/* Chunk context: the NDArray operands (ndops[0..mnd)), the once-folded scalar
 * constant (sc_re, sc_im — the sum of the scalar operands for Plus, their
 * product for Times), and the output. Pre-folding the scalars means the hot
 * loop only touches array buffers, so `c * arr` / `arr + c` vectorize. */
typedef struct {
    Expr** ndops; size_t mnd; double sc_re, sc_im;
    bool is_plus, all_f64; void* out; NDType dtc;
} nd_ew_ctx;

static bool nd_ew_chunk(void* c, size_t lo, size_t hi) {
    const nd_ew_ctx* x = (const nd_ew_ctx*)c;
    if (x->all_f64) {
        double* od = (double*)x->out;
        for (size_t k = lo; k < hi; k++) {
            double acc = x->sc_re;
            for (size_t j = 0; j < x->mnd; j++) {
                double v = ((const double*)x->ndops[j]->data.ndarray.data)[k];
                acc = x->is_plus ? acc + v : acc * v;
            }
            od[k] = acc;
        }
    } else {
        for (size_t k = lo; k < hi; k++) {
            double are = x->sc_re, aim = x->sc_im;
            for (size_t j = 0; j < x->mnd; j++) {
                double vre, vim;
                ndt_get(x->ndops[j]->data.ndarray.data, k,
                        x->ndops[j]->data.ndarray.dtype, &vre, &vim);
                if (x->is_plus) { are += vre; aim += vim; }
                else {
                    double nre = are * vre - aim * vim;
                    double nim = are * vim + aim * vre;
                    are = nre; aim = nim;
                }
            }
            ndt_set(x->out, k, x->dtc, are, aim);
        }
    }
    return true;
}

/* Repeat each element of `small` across the trailing axes to fill `dims`.
 * `small`'s dims must already be a prefix of `dims` (the caller checks). */
static Expr* nd_broadcast_to(const Expr* small, const int64_t* dims, int rank) {
    int r = small->data.ndarray.rank;
    NDType dt = small->data.ndarray.dtype;
    size_t esz = ndt_elem_size(dt);
    size_t nsmall = ndarray_size(small);
    size_t inner = 1;
    for (int i = r; i < rank; i++) inner *= (size_t)dims[i];
    size_t total = nsmall * inner;
    void* out = malloc(esz * (total ? total : 1));
    if (!out) return NULL;
    if (esz == sizeof(uint64_t)) {           /* float64 and int64: the common case */
        const uint64_t* s = (const uint64_t*)small->data.ndarray.data;
        uint64_t* d = (uint64_t*)out;
        for (size_t i = 0; i < nsmall; i++) {
            uint64_t v = s[i];
            uint64_t* row = d + i * inner;
            for (size_t j = 0; j < inner; j++) row[j] = v;
        }
    } else {
        const char* s = (const char*)small->data.ndarray.data;
        char* d = (char*)out;
        for (size_t i = 0; i < nsmall; i++)
            for (size_t j = 0; j < inner; j++)
                memcpy(d + (i * inner + j) * esz, s + i * esz, esz);
    }
    return expr_new_ndarray_like(small, rank, dims, out, dt);
}

Expr* ndarray_elementwise(Expr** args, size_t n, bool is_plus) {
    if (n == 0) return NULL;

    /* A bool array is not a numeric operand: {True, False} + {True, True} threads
     * to a List of unevaluated Plus[...] in Mathematica, it does not add as 1/0.
     * Decline so the arithmetic head's ndarray_int64_delist_retry (which also
     * covers bool) delists the operands and the List path reproduces that. */
    for (size_t i = 0; i < n; i++)
        if (args[i]->type == EXPR_NDARRAY &&
            args[i]->data.ndarray.dtype == NDT_BOOL) return NULL;

    /* A plain numeric List meeting a buffer: pack it first.
     *
     * The same trap Dot fell into (see dot_pack_operand in src/linalg/dot.c).
     * PACK_MIN_ELEMENTS judges each value on its own size, which is right for
     * blast radius and wrong for a binary operation, because the work is set by
     * the LARGEST operand. `RandomReal[{0,1},{8,100000}] - ConstantArray[0.3,8]`
     * has an 8-element right-hand side that is entirely correct to leave
     * unpacked, and that one decision put 800,000 elements back on the symbolic
     * path. Packing 8 doubles to avoid it is never a bad trade. */
    if (pack_any_created()) {
        bool have_nd = false, have_list = false;
        for (size_t i = 0; i < n; i++) {
            if (args[i]->type == EXPR_NDARRAY) have_nd = true;
            else if (head_is(args[i], SYM_List)) have_list = true;
        }
        if (have_nd && have_list) {
            Expr** packed = calloc(n, sizeof(Expr*));
            if (!packed) return NULL;
            bool changed = false;
            for (size_t i = 0; i < n; i++) {
                if (head_is(args[i], SYM_List)) {
                    packed[i] = pack_force(expr_copy(args[i]), false, NDT_FLOAT64);
                    if (packed[i] && packed[i]->type == EXPR_NDARRAY) changed = true;
                } else {
                    packed[i] = expr_copy(args[i]);
                }
            }
            /* `changed` also guards the recursion: without it a List that cannot
             * pack (symbolic elements) would re-enter here forever. */
            Expr* r = changed ? ndarray_elementwise(packed, n, is_plus) : NULL;
            for (size_t i = 0; i < n; i++) expr_free(packed[i]);
            free(packed);
            if (changed) return r;      /* NULL here means genuinely not ours */
        }
    }

    /* Leading-axis broadcast, i.e. Listable threading across levels: a rank-2
     * operand and a rank-1 one pair element i of the vector with ROW i of the
     * matrix, so {{1,2,3},{4,5,6}} + {10,20} is {{11,12,13},{24,25,26}}. This is
     * NOT numpy broadcasting -- the smaller operand's dims must be a PREFIX of
     * the larger's, and a length mismatch is an error in both this language and
     * Wolfram's -- but it is the same buffer job, and it had no buffer path at
     * all: `RandomReal[{0,1},{8,100000}] - ConstantArray[0.3,8]` measured 322 ms
     * for 800,000 elements, ~357 ns each, because it fell through to generic
     * symbolic threading. Subtracting a row vector from a matrix is about as
     * common as array code gets (centroid distances, mean removal, per-column
     * scaling).
     *
     * Done by expanding the smaller operand and reusing the same-shape kernel
     * below, rather than by giving that kernel a per-operand stride: Plus and
     * Times are on nearly every evaluation in the system and their inner loop is
     * not worth complicating for this. The cost is one extra pass over the
     * result -- roughly 2x the ideal -- against the ~80x this removes. */
    {
        const Expr* big = NULL;
        bool ranks_differ = false;
        for (size_t i = 0; i < n; i++) {
            if (args[i]->type != EXPR_NDARRAY) continue;
            if (!big) { big = args[i]; continue; }
            int rb = big->data.ndarray.rank, ri = args[i]->data.ndarray.rank;
            if (ri != rb) ranks_differ = true;
            if (ri > rb) big = args[i];
        }
        if (ranks_differ) {
            int rank = big->data.ndarray.rank;
            const int64_t* dims = big->data.ndarray.dims;
            /* Every array operand's shape must be a prefix of the widest one.
             * Anything else is a genuine shape error, which the generic path
             * reports; silently broadcasting it would invent a value. */
            for (size_t i = 0; i < n; i++) {
                if (args[i]->type != EXPR_NDARRAY) continue;
                int ri = args[i]->data.ndarray.rank;
                if (ri > rank) return NULL;
                for (int k = 0; k < ri; k++)
                    if (args[i]->data.ndarray.dims[k] != dims[k]) return NULL;
            }
            Expr** wide = calloc(n, sizeof(Expr*));   /* calloc: the cleanup below
                                                       * frees every non-NULL slot,
                                                       * including on a partial fill */
            if (!wide) return NULL;
            bool ok = true;
            for (size_t i = 0; i < n; i++) {
                if (args[i]->type == EXPR_NDARRAY &&
                    args[i]->data.ndarray.rank < rank) {
                    wide[i] = nd_broadcast_to(args[i], dims, rank);
                    if (!wide[i]) { ok = false; break; }
                } else {
                    wide[i] = expr_copy(args[i]);
                }
            }
            Expr* r = ok ? ndarray_elementwise(wide, n, is_plus) : NULL;
            for (size_t i = 0; i < n; i++) expr_free(wide[i]);
            free(wide);
            return r;
        }
    }

    const Expr* first_nd = NULL;
    /* Presentation of the result, joined over every array operand: a visible
     * NDArray on any side dominates a packed List (see nd_present_src2).
     * Distinct from first_nd, which is the SHAPE reference. */
    const Expr* present_src = NULL;
    NDType dtc = NDT_FLOAT64;
    bool have_dt = false, scalar_complex = false, all_nd_f64 = true;
    /* Exactness bookkeeping for an int64 result. `all_nd_int` == every array
     * operand is an integer buffer; the two scalar flags say whether the folded
     * scalar keeps the result exact (all operands exact Integers) or at least
     * machine-representable (a plain Real, which widens the result to float64).
     * Anything else -- Rational, BigInt, MPFR, Complex -- has no place in a
     * machine buffer alongside exact integers, because Range[10]/2 is a list of
     * Rationals and Range[10] + I a list of exact Complex; those go to the List
     * path where the interpreter builds them. */
    bool all_nd_int = true, sc_all_exact_int = true, sc_all_machine = true;
    int64_t sc_i = is_plus ? 0 : 1;
    /* Separate array operands from scalars; fold the scalars into one constant. */
    Expr** ndops = malloc(n * sizeof(Expr*));
    size_t mnd = 0;
    double sc_re = is_plus ? 0.0 : 1.0, sc_im = 0.0;
    for (size_t i = 0; i < n; i++) {
        if (args[i]->type == EXPR_NDARRAY) {
            if (!first_nd) first_nd = args[i];
            else if (!same_shape(first_nd, args[i])) { free(ndops); return NULL; }
            present_src = nd_present_src2(present_src, args[i]);
            NDType d = args[i]->data.ndarray.dtype;
            if (d != NDT_FLOAT64) all_nd_f64 = false;
            if (d != NDT_INT64)   all_nd_int = false;
            dtc = have_dt ? ndt_promote(dtc, d) : d;
            have_dt = true;
            ndops[mnd++] = args[i];
        } else {
            double re, im; bool isc;
            if (args[i]->type == EXPR_INTEGER) {
                int64_t v = args[i]->data.integer;
                if (is_plus ? ci_add_i64(sc_i, v, &sc_i) : ci_mul_i64(sc_i, v, &sc_i))
                    sc_all_exact_int = false;      /* folded past int64 */
            } else {
                sc_all_exact_int = false;
                if (args[i]->type != EXPR_REAL) sc_all_machine = false;
            }
            if (!scalar_value(args[i], &re, &im, &isc)) { free(ndops); return NULL; }
            if (isc) scalar_complex = true;
            if (is_plus) { sc_re += re; sc_im += im; }
            else {
                double nr = sc_re * re - sc_im * im;
                double ni = sc_re * im + sc_im * re;
                sc_re = nr; sc_im = ni;
            }
        }
    }
    if (!first_nd) { free(ndops); return NULL; }  /* no array operand — not ours */

    if (all_nd_int) {
        if (sc_all_exact_int) {
            /* Exact: Total[Range[10^6]^3] must promote to a BigInt rather than
             * wrap, so the first overflow abandons the whole result and the List
             * path recomputes it under GMP. */
            size_t szi = ndarray_size(first_nd);
            int64_t* o = malloc(sizeof(int64_t) * (szi ? szi : 1));
            if (!o) { free(ndops); return NULL; }
            for (size_t k = 0; k < szi; k++) {
                int64_t acc = sc_i;
                for (size_t j = 0; j < mnd; j++) {
                    int64_t v = ((const int64_t*)ndops[j]->data.ndarray.data)[k];
                    if (is_plus ? ci_add_i64(acc, v, &acc) : ci_mul_i64(acc, v, &acc)) {
                        free(o); free(ndops); return NULL;
                    }
                }
                o[k] = acc;
            }
            free(ndops);
            return expr_new_ndarray_like(present_src, first_nd->data.ndarray.rank,
                                    first_nd->data.ndarray.dims, o, NDT_INT64);
        }
        if (!sc_all_machine || scalar_complex) { free(ndops); return NULL; }
        /* A Real scalar makes every element inexact, so the buffer widens. Without
         * this the result kept dtype int64 and ndt_set TRUNCATED: Range[10] * 2.5
         * came out as {2, 5, 7, ...}. */
        dtc = NDT_FLOAT64;
        all_nd_f64 = false;
    }
    if (scalar_complex && !ndt_is_complex(dtc)) dtc = ndt_as_complex(dtc);

    size_t sz = ndarray_size(first_nd);
    void* out = malloc(ndt_elem_size(dtc) * sz);
    nd_ew_ctx ctx = { ndops, mnd, sc_re, sc_im, is_plus,
                      (all_nd_f64 && dtc == NDT_FLOAT64), out, dtc };
    nd_parallel_for(sz, nd_ew_chunk, &ctx);
    free(ndops);
    /* takes ownership of out */
    return expr_new_ndarray_like(present_src, first_nd->data.ndarray.rank,
                            first_nd->data.ndarray.dims, out, dtc);
}

/* --------------------------- elementwise Power --------------------------- */

/* Complex power (br+bi i) = (ar+ai i)^(er+ei i) via polar form, avoiding
 * <complex.h>. Real, non-negative real base with integer/real exponent stays
 * on the real axis; a negative real base with a non-integer exponent yields a
 * genuinely complex value (caller promotes the result dtype accordingly). */
static void cpow_pair(double ar, double ai, double er, double ei,
                      double* outr, double* outi) {
    if (ai == 0.0 && ei == 0.0 && (ar >= 0.0 || er == floor(er))) {
        /* Stays real. Small integer powers (the common x^2, x^3, 1/x, ...) are
         * a few multiplies — far cheaper than a libm pow() call per element. */
        if (er == 2.0)       { *outr = ar * ar;           *outi = 0.0; return; }
        if (er == 3.0)       { *outr = ar * ar * ar;      *outi = 0.0; return; }
        if (er == 0.5 && ar >= 0.0) { *outr = sqrt(ar);   *outi = 0.0; return; }
        if (er == 1.0)       { *outr = ar;                *outi = 0.0; return; }
        if (er == -1.0)      { *outr = 1.0 / ar;          *outi = 0.0; return; }
        if (er == 0.0)       { *outr = 1.0;               *outi = 0.0; return; }
        *outr = pow(ar, er); *outi = 0.0; return;
    }
    /* Complex base with an integer exponent: binary exponentiation gives an
     * exact result (e.g. I^2 == -1 + 0 I) without polar-form round-off. */
    if (ei == 0.0 && er == floor(er) && fabs(er) <= 1024.0) {
        long p = (long)er;
        unsigned long m = (p < 0) ? (unsigned long)(-p) : (unsigned long)p;
        double rr = 1.0, ri = 0.0, br = ar, bi = ai;
        while (m) {
            if (m & 1ul) {
                double nr = rr * br - ri * bi;
                double ni = rr * bi + ri * br;
                rr = nr; ri = ni;
            }
            double sr = br * br - bi * bi;
            double si = 2.0 * br * bi;
            br = sr; bi = si;
            m >>= 1;
        }
        if (p < 0) {
            double d = rr * rr + ri * ri;
            *outr = rr / d; *outi = -ri / d;
        } else { *outr = rr; *outi = ri; }
        return;
    }
    if (ar == 0.0 && ai == 0.0) { *outr = 0.0; *outi = 0.0; return; }
    double logr = 0.5 * log(ar * ar + ai * ai);   /* log|a|            */
    double th = atan2(ai, ar);                     /* arg(a)            */
    /* log(a) = logr + i th ; (er+ei i)*log(a) = p + i q */
    double p = er * logr - ei * th;
    double q = er * th + ei * logr;
    double mag = exp(p);
    *outr = mag * cos(q);
    *outi = mag * sin(q);
}

/* True when raising a real-dtype array (base) to `(er,ei)` can leave the real
 * axis: a complex exponent, or a non-integer exponent applied to a possibly
 * negative base. Used to decide whether the result must promote to complex. */
static bool real_power_escapes(const Expr* base, double er, double ei) {
    if (ei != 0.0) return true;
    if (er == floor(er)) return false;             /* integer power stays real */
    size_t sz = ndarray_size(base);
    NDType dt = base->data.ndarray.dtype;
    for (size_t k = 0; k < sz; k++) {
        double re, im;
        ndt_get(base->data.ndarray.data, k, dt, &re, &im);
        if (re < 0.0) return true;
    }
    return false;
}

/* Threaded power over one output element per index. Either operand may be an
 * NDArray (buf non-NULL, dtype dt) or a broadcast scalar (buf NULL, use the
 * s* value). cpow_pair is the pure per-element kernel. */
typedef struct {
    const void* abuf; NDType adt; double sar, sai;   /* base:  buf or scalar */
    const void* bbuf; NDType bdt; double sbr, sbi;   /* exp:   buf or scalar */
    void* out; NDType dtc;
} nd_pow_ctx;
static bool nd_pow_chunk(void* c, size_t lo, size_t hi) {
    const nd_pow_ctx* x = (const nd_pow_ctx*)c;
    for (size_t k = lo; k < hi; k++) {
        double ar, ai, er, ei, rr, ri;
        if (x->abuf) ndt_get(x->abuf, k, x->adt, &ar, &ai); else { ar = x->sar; ai = x->sai; }
        if (x->bbuf) ndt_get(x->bbuf, k, x->bdt, &er, &ei); else { er = x->sbr; ei = x->sbi; }
        cpow_pair(ar, ai, er, ei, &rr, &ri);
        ndt_set(x->out, k, x->dtc, rr, ri);
    }
    return true;
}

Expr* ndarray_elementwise_power(const Expr* a, const Expr* b) {
    if (!is_ndarray(a) || !is_ndarray(b) || !same_shape(a, b)) return NULL;
    NDType dta = a->data.ndarray.dtype, dtb = b->data.ndarray.dtype;
    if (dta == NDT_BOOL || dtb == NDT_BOOL) return NULL;  /* not numeric; see ndarray_map_unary */
    size_t sz = ndarray_size(a);

    if (dta == NDT_INT64 || dtb == NDT_INT64) {
        /* Exact only when BOTH sides are integer buffers and every exponent is
         * non-negative; a mixed int/float pair or any negative exponent has an
         * exact answer the buffer cannot express. */
        if (dta != NDT_INT64 || dtb != NDT_INT64) return NULL;
        const int64_t* pa = (const int64_t*)a->data.ndarray.data;
        const int64_t* pb = (const int64_t*)b->data.ndarray.data;
        for (size_t k = 0; k < sz; k++) if (pb[k] < 0) return NULL;
        int64_t* o = malloc(sizeof(int64_t) * (sz ? sz : 1));
        if (!o) return NULL;
        for (size_t k = 0; k < sz; k++)
            if (ci_powi_i64(pa[k], pb[k], &o[k])) { free(o); return NULL; }
        return expr_new_ndarray_like(a, a->data.ndarray.rank, a->data.ndarray.dims,
                                o, NDT_INT64);
    }

    /* Result promotes to complex if either operand is complex, or a real base
     * escapes the real axis at any element. */
    NDType dtc = ndt_promote(dta, dtb);
    if (!ndt_is_complex(dtc)) {
        for (size_t k = 0; k < sz; k++) {
            double ar, ai, er, ei;
            ndt_get(a->data.ndarray.data, k, dta, &ar, &ai);
            ndt_get(b->data.ndarray.data, k, dtb, &er, &ei);
            if (er != floor(er) && ar < 0.0) { dtc = ndt_as_complex(dtc); break; }
        }
    }

    void* out = malloc(ndt_elem_size(dtc) * sz);
    nd_pow_ctx c = { a->data.ndarray.data, dta, 0, 0,
                     b->data.ndarray.data, dtb, 0, 0, out, dtc };
    nd_parallel_for(sz, nd_pow_chunk, &c);
    return expr_new_ndarray_like(a, a->data.ndarray.rank, a->data.ndarray.dims, out, dtc);
}

Expr* ndarray_base_scalar_power(double br, double bi, const Expr* exp_arr) {
    if (!is_ndarray(exp_arr)) return NULL;
    NDType dte = exp_arr->data.ndarray.dtype;
    if (dte == NDT_BOOL) return NULL;      /* not numeric; see ndarray_map_unary */
    size_t sz = ndarray_size(exp_arr);

    /* 2^Range[300] is a list of exact Integers, but the base arrives here already
     * flattened to a double, so this routine cannot tell an Integer base from a
     * Real one. Rather than guess, hand an integer exponent buffer to the List
     * path -- rare shape, exact answer. */
    if (dte == NDT_INT64) return NULL;

    /* Result is complex if the base is complex, the exponent array is complex,
     * or a negative real base is raised to a non-integer exponent anywhere. */
    NDType dtc = dte;
    if (bi != 0.0) dtc = ndt_as_complex(dtc);
    else if (!ndt_is_complex(dte) && br < 0.0) {
        for (size_t k = 0; k < sz; k++) {
            double er, ei;
            ndt_get(exp_arr->data.ndarray.data, k, dte, &er, &ei);
            if (er != floor(er)) { dtc = ndt_as_complex(dtc); break; }
        }
    }

    void* out = malloc(ndt_elem_size(dtc) * sz);
    nd_pow_ctx c = { NULL, dte, br, bi,             /* scalar base */
                     exp_arr->data.ndarray.data, dte, 0, 0, out, dtc };
    nd_parallel_for(sz, nd_pow_chunk, &c);
    return expr_new_ndarray_like(exp_arr, exp_arr->data.ndarray.rank,
                            exp_arr->data.ndarray.dims, out, dtc);
}

Expr* ndarray_scalar_power(const Expr* a, double er, double ei) {
    if (!is_ndarray(a)) return NULL;
    NDType dta = a->data.ndarray.dtype;
    if (dta == NDT_BOOL) return NULL;      /* not numeric; see ndarray_map_unary */
    size_t sz = ndarray_size(a);

    if (dta == NDT_INT64) {
        /* Range[300]^2 must stay a list of Integers. Only a NON-NEGATIVE integer
         * exponent keeps an integer result: a negative one gives Rationals and a
         * fractional one gives radicals, neither of which a buffer can hold, so
         * both go to the List path. ci_powi_i64 abandons on overflow the same
         * way, so Total[Range[10^6]^3] still reaches GMP and answers exactly. */
        if (ei != 0.0 || er != floor(er) || er < 0.0) return NULL;
        int64_t* o = malloc(sizeof(int64_t) * (sz ? sz : 1));
        if (!o) return NULL;
        const int64_t* p = (const int64_t*)a->data.ndarray.data;
        for (size_t k = 0; k < sz; k++)
            if (ci_powi_i64(p[k], (int64_t)er, &o[k])) { free(o); return NULL; }
        return expr_new_ndarray_like(a, a->data.ndarray.rank, a->data.ndarray.dims,
                                o, NDT_INT64);
    }

    NDType dtc = dta;
    if (ei != 0.0) dtc = ndt_as_complex(dtc);
    else if (!ndt_is_complex(dta) && real_power_escapes(a, er, ei))
        dtc = ndt_as_complex(dtc);

    void* out = malloc(ndt_elem_size(dtc) * sz);
    nd_pow_ctx c = { a->data.ndarray.data, dta, 0, 0,
                     NULL, dta, er, ei, out, dtc };   /* scalar exponent */
    nd_parallel_for(sz, nd_pow_chunk, &c);
    return expr_new_ndarray_like(a, a->data.ndarray.rank, a->data.ndarray.dims, out, dtc);
}

/* -------------------- parallel element-wise map ------------------------- */
/* The element-wise kernel loops below are embarrassingly parallel: each output
 * element depends only on the same-index input element, and every kernel is a
 * pure function (libm / arithmetic, thread-safe). For large arrays we split the
 * flat index range across worker threads, which is where the wall-clock cost of
 * e.g. Erf[NDArray[...]] over millions of elements actually lives. Guarded by
 * MATHILDA_THREADS (set by the makefile with -pthread); without it every map
 * runs serially, so the CMake test build and any thread-less platform are
 * unaffected. */

/* A chunk callback processes the half-open flat range [lo, hi). It returns false
 * to signal that the kernel declined an element (out-of-fast-domain); the caller
 * then frees its buffer and degrades to the exact per-element List path, exactly
 * as the serial loops did on a false kernel return. (nd_chunk_fn and
 * NDARRAY_MAX_THREADS are declared in ndarray_internal.h.) */

#ifdef MATHILDA_THREADS
#include <pthread.h>
#include <unistd.h>

/* Arrays smaller than this run serially: below it the thread spawn/join cost
 * dominates the actual kernel work. */
typedef struct { nd_chunk_fn fn; void* ctx; size_t lo, hi; bool ok; } nd_thread_job;

static void* nd_thread_run(void* p) {
    nd_thread_job* j = (nd_thread_job*)p;
    j->ok = j->fn(j->ctx, j->lo, j->hi);
    return NULL;
}

/* Threads to use for `n` elements: 1 for small arrays, else one per online CPU
 * (capped), and never so many that a chunk falls below the threshold. */
static int nd_thread_count(size_t n) {
    if (n < NDARRAY_THREAD_THRESHOLD) return 1;
    long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
    if (ncpu < 1) ncpu = 1;
    if (ncpu > NDARRAY_MAX_THREADS) ncpu = NDARRAY_MAX_THREADS;
    size_t by_work = n / NDARRAY_THREAD_THRESHOLD;
    if (by_work < 1) by_work = 1;
    if ((size_t)ncpu > by_work) ncpu = (long)by_work;
    return (int)ncpu;
}
#endif /* MATHILDA_THREADS */

/* Run `fn` over [0, n), in parallel for large arrays when threads are available,
 * serially otherwise. Returns true iff every chunk returned true. */
bool nd_parallel_for(size_t n, nd_chunk_fn fn, void* ctx) {
#ifdef MATHILDA_THREADS
    int nt = nd_thread_count(n);
    if (nt > 1) {
        pthread_t th[NDARRAY_MAX_THREADS];
        nd_thread_job jb[NDARRAY_MAX_THREADS];
        bool threaded[NDARRAY_MAX_THREADS];
        size_t chunk = (n + (size_t)nt - 1) / (size_t)nt;
        int nc = 0;
        for (int t = 0; t < nt; t++) {
            size_t lo = (size_t)t * chunk;
            if (lo >= n) break;
            size_t hi = lo + chunk;
            if (hi > n) hi = n;
            jb[nc].fn = fn; jb[nc].ctx = ctx; jb[nc].lo = lo; jb[nc].hi = hi;
            jb[nc].ok = true;
            threaded[nc] = (pthread_create(&th[nc], NULL, nd_thread_run, &jb[nc]) == 0);
            if (!threaded[nc]) jb[nc].ok = fn(ctx, lo, hi); /* spawn failed: inline */
            nc++;
        }
        bool all = true;
        for (int t = 0; t < nc; t++) {
            if (threaded[t]) pthread_join(th[t], NULL);
            all = all && jb[t].ok;
        }
        return all;
    }
#endif
    return fn(ctx, 0, n);
}

/* Parallel reduction — each chunk folds its range into a private slot; the
 * caller combines the returned slots. See ndarray_internal.h. */
#ifdef MATHILDA_THREADS
typedef struct {
    nd_reduce_chunk_fn fn; void* ctx; size_t lo, hi; double* slot;
} nd_reduce_job;
static void* nd_reduce_run(void* p) {
    nd_reduce_job* j = (nd_reduce_job*)p;
    j->fn(j->ctx, j->lo, j->hi, j->slot);
    return NULL;
}
#endif

int nd_parallel_reduce(size_t n, nd_reduce_chunk_fn fn, void* ctx,
                       int ncomp, double* slots) {
    (void)ncomp;
#ifdef MATHILDA_THREADS
    int nt = nd_thread_count(n);
    if (nt > 1) {
        pthread_t th[NDARRAY_MAX_THREADS];
        nd_reduce_job jb[NDARRAY_MAX_THREADS];
        bool threaded[NDARRAY_MAX_THREADS];
        size_t chunk = (n + (size_t)nt - 1) / (size_t)nt;
        int nc = 0;
        for (int t = 0; t < nt; t++) {
            size_t lo = (size_t)t * chunk;
            if (lo >= n) break;
            size_t hi = lo + chunk;
            if (hi > n) hi = n;
            jb[nc].fn = fn; jb[nc].ctx = ctx; jb[nc].lo = lo; jb[nc].hi = hi;
            jb[nc].slot = slots + (size_t)nc * (size_t)ncomp;
            threaded[nc] = (pthread_create(&th[nc], NULL, nd_reduce_run, &jb[nc]) == 0);
            if (!threaded[nc]) fn(ctx, lo, hi, jb[nc].slot); /* spawn failed: inline */
            nc++;
        }
        for (int t = 0; t < nc; t++)
            if (threaded[t]) pthread_join(th[t], NULL);
        return nc;
    }
#endif
    fn(ctx, 0, n, slots);
    return 1;
}

/* Chunk contexts + callbacks for the straight (index-preserving) maps. */
typedef struct { const double* id; double* od; const NDUnaryKernel* k; } ndu_hot_ctx;
static bool ndu_hot_chunk(void* c, size_t lo, size_t hi) {
    const ndu_hot_ctx* x = (const ndu_hot_ctx*)c;
    for (size_t j = lo; j < hi; j++)
        if (!x->k->real(x->id[j], &x->od[j])) return false;
    return true;
}

typedef struct { const void* in; void* out; NDType dt; const NDUnaryKernel* k; } ndu_map_ctx;
static bool ndu_cplx_chunk(void* c, size_t lo, size_t hi) {
    const ndu_map_ctx* x = (const ndu_map_ctx*)c;
    for (size_t j = lo; j < hi; j++) {
        double re, im, orr, oii;
        ndt_get(x->in, j, x->dt, &re, &im);
        if (!x->k->cplx(re, im, &orr, &oii)) return false;
        ndt_set(x->out, j, x->dt, orr, oii);
    }
    return true;
}
static bool ndu_realgen_chunk(void* c, size_t lo, size_t hi) {
    const ndu_map_ctx* x = (const ndu_map_ctx*)c;
    for (size_t j = lo; j < hi; j++) {
        double re, im, orr, oii = 0.0;
        ndt_get(x->in, j, x->dt, &re, &im);
        if (x->k->real) { if (!x->k->real(re, &orr)) return false; }
        else if (!x->k->cplx || !x->k->cplx(re, 0.0, &orr, &oii)) return false;
        ndt_set(x->out, j, x->dt, orr, oii);
    }
    return true;
}

typedef struct { const void* in; void* out; NDType dti, dto; const NDUnaryKernel* k; } ndu_proj_ctx;
static bool ndu_proj_chunk(void* c, size_t lo, size_t hi) {
    const ndu_proj_ctx* x = (const ndu_proj_ctx*)c;
    for (size_t j = lo; j < hi; j++) {
        double re, im, orr, oii;
        ndt_get(x->in, j, x->dti, &re, &im);
        if (!x->k->cplx(re, im, &orr, &oii)) return false;
        ndt_set(x->out, j, x->dto, orr, 0.0);
    }
    return true;
}

typedef struct {
    const void* in; void* out; NDType dti, dto; const NDBinaryKernel* k;
    double sre, sim; bool arr_first;
} ndb_ctx;
static bool ndb_map_chunk(void* c, size_t lo, size_t hi) {
    const ndb_ctx* x = (const ndb_ctx*)c;
    /* float64 in and out: raw doubles, no ndt_get/ndt_set indirection. Same
     * lane as ndb2_map_chunk below -- see the note there. Reached whenever the
     * result is real, i.e. the real_closed branch; the complex branch has a
     * complex `dto` and falls to the generic loop. */
    if (x->dti == NDT_FLOAT64 && x->dto == NDT_FLOAT64) {
        const double* in = (const double*)x->in;
        double* o = (double*)x->out;
        double s = x->sre;
        for (size_t j = lo; j < hi; j++) {
            double orr, oii;
            bool ok = x->arr_first ? x->k->cplx(in[j], 0.0, s, 0.0, &orr, &oii)
                                   : x->k->cplx(s, 0.0, in[j], 0.0, &orr, &oii);
            if (!ok) return false;
            o[j] = orr;
        }
        return true;
    }
    for (size_t j = lo; j < hi; j++) {
        double re, im, orr, oii;
        ndt_get(x->in, j, x->dti, &re, &im);
        bool ok = x->arr_first ? x->k->cplx(re, im, x->sre, x->sim, &orr, &oii)
                               : x->k->cplx(x->sre, x->sim, re, im, &orr, &oii);
        if (!ok) return false;
        ndt_set(x->out, j, x->dto, orr, oii);
    }
    return true;
}

/* The same map with TWO array operands rather than one array and a broadcast
 * scalar: element j of the result comes from element j of each input. */
typedef struct {
    const void* in0; const void* in1; void* out;
    NDType dt0, dt1, dto; const NDBinaryKernel* k;
} ndb2_ctx;
static bool ndb2_map_chunk(void* c, size_t lo, size_t hi) {
    const ndb2_ctx* x = (const ndb2_ctx*)c;
    /* float64 throughout: read and write the raw doubles, the binary twin of
     * ndu_hot_chunk. ndt_get/ndt_set are the right default -- they are the one
     * place that knows every dtype -- but they are indirect calls the compiler
     * cannot inline through, and there are three of them per element here
     * against the kernel's one.
     *
     * Measured by building the tree with and without it, 10^6 elements:
     * ArcTan[v, w] 3.26 ms -> 2.49 ms (23%) and Beta[p, p] 5.27 ms -> 4.94 ms
     * (6%). Beta gains less because three lgammas per element dominate what the
     * marshalling costs, where atan2 does not -- so the lane is worth most
     * exactly where the kernel is cheapest, which is the same reason it is most
     * of the win for the unary elementary set (ndu_hot_chunk). */
    if (x->dt0 == NDT_FLOAT64 && x->dt1 == NDT_FLOAT64 && x->dto == NDT_FLOAT64) {
        const double* a = (const double*)x->in0;
        const double* b = (const double*)x->in1;
        double* o = (double*)x->out;
        for (size_t j = lo; j < hi; j++) {
            double orr, oii;
            if (!x->k->cplx(a[j], 0.0, b[j], 0.0, &orr, &oii)) return false;
            o[j] = orr;
        }
        return true;
    }
    for (size_t j = lo; j < hi; j++) {
        double are, aim, bre, bim, orr, oii;
        ndt_get(x->in0, j, x->dt0, &are, &aim);
        ndt_get(x->in1, j, x->dt1, &bre, &bim);
        if (!x->k->cplx(are, aim, bre, bim, &orr, &oii)) return false;
        ndt_set(x->out, j, x->dto, orr, oii);
    }
    return true;
}

/* -------------------- element-wise scalar-function map ------------------- */

/* Escaping-real-input tail: compute the complex kernel into `tmp`, flag whether
 * any element left the real axis. Returns false if the kernel declines (caller
 * degrades). Each chunk sets *any_imag at most once (monotonic 0->1 flag); the
 * concurrent same-value writes are benign and only read after all chunks join. */
typedef struct {
    const void* in; NDType dti; void* tmp; NDType dtc; const NDUnaryKernel* k; int* any_imag;
} ndu_esc_ctx;
static bool ndu_esc_chunk(void* c, size_t lo, size_t hi) {
    const ndu_esc_ctx* x = (const ndu_esc_ctx*)c;
    int ai = 0;
    for (size_t j = lo; j < hi; j++) {
        double re, im, orr, oii;
        ndt_get(x->in, j, x->dti, &re, &im);
        if (!x->k->cplx(re, 0.0, &orr, &oii)) return false;
        if (oii != 0.0) ai = 1;
        ndt_set(x->tmp, j, x->dtc, orr, oii);
    }
    if (ai) *x->any_imag = 1;
    return true;
}
/* Narrow a complex buffer back to a real dtype (drop the zero imaginary parts). */
typedef struct { const void* tmp; NDType dtc; void* out; NDType dto; } ndu_narrow_ctx;
static bool ndu_narrow_chunk(void* c, size_t lo, size_t hi) {
    const ndu_narrow_ctx* x = (const ndu_narrow_ctx*)c;
    for (size_t j = lo; j < hi; j++) {
        double re, im;
        ndt_get(x->tmp, j, x->dtc, &re, &im);
        ndt_set(x->out, j, x->dto, re, im);
    }
    return true;
}

Expr* ndarray_map_unary(const Expr* a, const NDUnaryKernel* k) {
    if (!is_ndarray(a) || !k) return NULL;
    NDType dta = a->data.ndarray.dtype;
    /* A bool array is not a numeric one: Sin[True], Floor[False] and the rest are
     * left unevaluated by Mathematica, not folded to a number. Decline the whole
     * map (including the to_int narrowing arm below) so the caller degrades to
     * ndarray_delist_and_reeval, whose answer is exactly Sin[{True, False}] --
     * i.e. the symbolic List. Every numeric engine in this file takes the same
     * stance; see the NDT_BOOL guards in ndarray_elementwise / _power / _dot. */
    if (dta == NDT_BOOL) return NULL;
    const void* in = a->data.ndarray.data;
    size_t sz = ndarray_size(a);
    int rank = a->data.ndarray.rank;
    const int64_t* dims = a->data.ndarray.dims;

    /* Narrowing (Floor/Ceiling/Round/IntegerPart/Sign/UnitStep): a real or
     * exact-integer input gives an EXACT INTEGER, so the result is NDT_INT64.
     *
     * Deliberately serial and all-or-nothing rather than threaded: a single
     * element that is not representable as an int64 (non-finite, or past the
     * range) abandons the WHOLE array so the ordinary List path can answer with
     * a bignum, and that early exit is worth more than the threading. It is also
     * the same contract the compile engine and Total's int64 accumulate use --
     * abandon, never wrap, never round.
     *
     * Complex input has no integer answer here (Sign[z] is z/|z|), so it falls
     * through to the categories below, which decline for these kernels. */
    if (k->to_int && !ndt_is_complex(dta)) {
        bool from_int = (dta == NDT_INT64);
        if (from_int ? (k->to_int_i != NULL) : (k->to_int_r != NULL)) {
            int64_t* out = (int64_t*)malloc(sizeof(int64_t) * (sz ? sz : 1));
            if (!out) return NULL;
            bool ok = true;
            if (from_int) {
                const int64_t* src = (const int64_t*)in;
                for (size_t i = 0; i < sz && ok; i++)
                    ok = k->to_int_i(src[i], &out[i]);
            } else if (dta == NDT_FLOAT64) {
                const double* src = (const double*)in;
                for (size_t i = 0; i < sz && ok; i++)
                    ok = k->to_int_r(src[i], &out[i]);
            } else {
                for (size_t i = 0; i < sz && ok; i++) {
                    double re, im;
                    ndt_get(in, i, dta, &re, &im);
                    ok = k->to_int_r(re, &out[i]);
                }
            }
            if (!ok) { free(out); return NULL; }   /* degrade to the List path */
            return expr_new_ndarray_like(a, rank, dims, out, NDT_INT64);
        }
        /* No arm for THIS dtype: fall through to the categories below rather
         * than declining the whole call. `to_int` means "prefer the exact
         * integer answer where one exists", not "only integers" -- Abs is
         * exactly that shape: |x| is an exact Integer on an int64 buffer and a
         * projection (to_real) on a Real or complex one, and before this it had
         * to choose one of the two. The six narrowing kernels supply both arms,
         * so nothing they do changes. */
    }

    /* AN int64 BUFFER WITH NO EXACT ARM HAS NO HONEST ANSWER BELOW.
     *
     * Every remaining branch sizes its output buffer from the INPUT dtype and
     * writes through ndt_set, whose NDT_INT64 case is `(int64_t)re` -- so a
     * double result lands in an integer slot and is TRUNCATED. Measured, before
     * this guard:
     *
     *     Sin[NDArray[{1, 2, 3}, DataType -> "int64"]]  ->  {0, 0, 0}
     *     Exp[NDArray[{1, 2, 3}, DataType -> "int64"]]  ->  {2, 7, 20}
     *     Erf, Tanh, Cos, BesselJ, Log[2, .]            ->  likewise
     *
     * silently, for every real_closed kernel without an exact arm -- 56 of the
     * 85 registered. The to_real projections were wrong here too, differently:
     * Re/Im/Arg wrote a FLOAT64 result from an integer buffer, so Re[intArray]
     * gave {1., 2., 3.} where the List gives the exact {1, 2, 3}. That is why
     * Im sits on pack.c's NOT_AWARE list; this guard covers the visible surface
     * the same way.
     *
     * A PACKED list never reached this: src/eval.c's transparency gate
     * materialises an int64 buffer for any head that has not claimed
     * packed_int64_ok, so the List path answered. The VISIBLE NDArray surface
     * has no such gate (`is_packed_list` is false for it, by design), and
     * nothing else stood between an integer buffer and a double kernel.
     *
     * Declining is the same rule the gate applies, for the same reason:
     * materialising is always an option, being wrong is not. The caller
     * degrades to ndarray_delist_and_reeval and the answer is exactly what
     * Sin[{1, 2, 3}] gives, so all three representations now agree.
     *
     * This does not touch the exact arms: Floor/Ceiling/Round/IntegerPart/
     * Sign/UnitStep/Abs (and Mod/Quotient below) return above via to_int_i and
     * keep their integer fast path. */
    if (dta == NDT_INT64) return NULL;

    /* Projection (Abs/Re/Im/Arg): always a real output, even for complex input.
     * The real dtype keeps the input's component width. */
    if (k->to_real) {
        if (!k->cplx) return NULL;   /* degrade */
        NDType dtr = (ndt_comp_size(dta) == sizeof(double)) ? NDT_FLOAT64
                                                            : NDT_FLOAT32;
        void* out = malloc(ndt_elem_size(dtr) * sz);
        ndu_proj_ctx ctx = { in, out, dta, dtr, k };
        if (!nd_parallel_for(sz, ndu_proj_chunk, &ctx)) { free(out); return NULL; }
        return expr_new_ndarray_like(a, rank, dims, out, dtr);
    }

    /* Complex input: complex output, straight map. A real-only kernel (no cplx)
     * declines here, degrading a complex NDArray to the List path. */
    if (ndt_is_complex(dta)) {
        if (!k->cplx) return NULL;
        void* out = malloc(ndt_elem_size(dta) * sz);
        ndu_map_ctx ctx = { in, out, dta, k };
        if (!nd_parallel_for(sz, ndu_cplx_chunk, &ctx)) { free(out); return NULL; }
        return expr_new_ndarray_like(a, rank, dims, out, dta);
    }

    /* Real input, function stays on the real axis: real output. */
    if (k->real_closed) {
        void* out = malloc(ndt_elem_size(dta) * sz);
        if (dta == NDT_FLOAT64 && k->real) {
            /* Hot path: raw double buffers, no widening. */
            ndu_hot_ctx ctx = { (const double*)in, (double*)out, k };
            if (!nd_parallel_for(sz, ndu_hot_chunk, &ctx)) { free(out); return NULL; }
        } else {
            ndu_map_ctx ctx = { in, out, dta, k };
            if (!nd_parallel_for(sz, ndu_realgen_chunk, &ctx)) { free(out); return NULL; }
        }
        return expr_new_ndarray_like(a, rank, dims, out, dta);
    }

    /* Real input, function may leave the real axis (Sqrt/Log/ArcSin/...):
     * compute complex, then narrow to real dtype iff every element is real,
     * else promote (mirrors real_power_escapes in ndarray_scalar_power). */
    if (!k->cplx) return NULL;   /* degrade (sentinel / real-only kernel) */
    NDType dtcplx = ndt_as_complex(dta);
    void* tmp = malloc(ndt_elem_size(dtcplx) * sz);
    int any_imag = 0;
    ndu_esc_ctx ectx = { in, dta, tmp, dtcplx, k, &any_imag };
    if (!nd_parallel_for(sz, ndu_esc_chunk, &ectx)) { free(tmp); return NULL; }
    if (any_imag) return expr_new_ndarray_like(a, rank, dims, tmp, dtcplx);
    /* Narrow back to the original real dtype. */
    void* out = malloc(ndt_elem_size(dta) * sz);
    ndu_narrow_ctx nctx = { tmp, dtcplx, out, dta };
    nd_parallel_for(sz, ndu_narrow_chunk, &nctx);
    free(tmp);
    return expr_new_ndarray_like(a, rank, dims, out, dta);
}

Expr* ndarray_map_binary(const Expr* a0, const Expr* a1, const NDBinaryKernel* k) {
    /* A kernel with no `cplx` is not necessarily a sentinel any more: GCD, LCM
     * and DivisorSigma are defined on the integers and have `to_int_i` only.
     * The test therefore moved down to each branch that actually calls `cplx`,
     * where a sentinel still declines and the caller still degrades to the List
     * path -- unchanged behaviour for every kernel that has one. */
    if (!k) return NULL;
    /* Exactly one operand is an NDArray; the other is a broadcast scalar. */
    const Expr* arr; const Expr* scal; bool arr_first;
    if (is_ndarray(a0) && !is_ndarray(a1))      { arr = a0; scal = a1; arr_first = true; }
    else if (is_ndarray(a1) && !is_ndarray(a0)) { arr = a1; scal = a0; arr_first = false; }
    else return NULL;

    double sre, sim; bool s_cplx;
    if (!scalar_value(scal, &sre, &sim, &s_cplx)) return NULL;

    NDType dta = arr->data.ndarray.dtype;
    if (dta == NDT_BOOL) return NULL;      /* not numeric; see ndarray_map_unary */
    const void* in = arr->data.ndarray.data;
    size_t sz = ndarray_size(arr);
    int rank = arr->data.ndarray.rank;
    const int64_t* dims = arr->data.ndarray.dims;

    /* EXACT INTEGER arms (Mod, Quotient) -- the binary twin of the narrowing
     * branch in ndarray_map_unary, and it comes FIRST for the same reason: an
     * exact answer that exists must not be computed in double and rounded back.
     *
     * Serial and all-or-nothing on purpose: one element that has no exact int64
     * answer (a zero divisor, an INT64_MIN / -1 overflow, a real too large to
     * narrow) abandons the whole array so the List path answers with GMP.
     * Abandon, never wrap. */
    if (k->to_int && !ndt_is_complex(dta) && !s_cplx) {
        int64_t nsi = 0;
        bool scal_exact_int = (scal->type == EXPR_INTEGER);
        if (scal_exact_int) nsi = scal->data.integer;

        bool use_i = (dta == NDT_INT64) && scal_exact_int && k->to_int_i;
        /* The real arm reads elements through ndt_get, which routes an int64
         * through `double` and is exact only to 2^53. That is faithful for a
         * buffer whose elements are ALREADY inexact and a wrong answer for one
         * whose elements are not: Quotient[{10^17, ...}, 1/2] would divide a
         * rounded element where the List path divides in mpq. So an int64
         * buffer takes the exact arm or nothing, and "nothing" means the
         * NDT_INT64 guard below declines and the List path answers. */
        bool use_r = !use_i && k->to_int_r && dta != NDT_INT64;
        if (use_i || use_r) {
            int64_t* out = (int64_t*)malloc(sizeof(int64_t) * (sz ? sz : 1));
            if (!out) return NULL;
            bool ok = true;
            for (size_t i = 0; i < sz && ok; i++) {
                if (use_i) {
                    int64_t e = ndt_get_i(in, i, dta);
                    ok = arr_first ? k->to_int_i(e, nsi, &out[i])
                                   : k->to_int_i(nsi, e, &out[i]);
                } else {
                    double re, im;
                    ndt_get(in, i, dta, &re, &im);
                    ok = arr_first ? k->to_int_r(re, sre, &out[i])
                                   : k->to_int_r(sre, re, &out[i]);
                }
            }
            if (!ok) { free(out); return NULL; }   /* degrade to the List path */
            return expr_new_ndarray_like(arr, rank, dims, out, NDT_INT64);
        }
        /* No arm for THIS dtype/scalar pair -- an int64 buffer with a Real
         * scalar, which is Mod[{1,2,3}, 2.] = {1., 0., 1.}. Falls through to
         * the guard below, which declines so the List path answers with the
         * Real elements the buffer could not have held. */
    }

    /* An int64 buffer with no exact arm above: see the same guard in
     * ndarray_map_unary. Every branch below sizes its output from the INPUT
     * dtype and writes through ndt_set, so a double result would be truncated
     * into an integer slot. Decline instead, and the caller's
     * ndarray_delist_and_reeval answers exactly as the List does. */
    if (dta == NDT_INT64) return NULL;

    if (!k->cplx) return NULL;   /* integer-only kernel, or sentinel: degrade */

    /* Genuinely-complex output: complex input or complex scalar. Map directly. */
    if (ndt_is_complex(dta) || s_cplx) {
        NDType dtc = ndt_as_complex(dta);
        void* out = malloc(ndt_elem_size(dtc) * sz);
        ndb_ctx ctx = { in, out, dta, dtc, k, sre, sim, arr_first };
        if (!nd_parallel_for(sz, ndb_map_chunk, &ctx)) { free(out); return NULL; }
        return expr_new_ndarray_like(arr, rank, dims, out, dtc);
    }

    /* Real inputs, real-closed function: real output. */
    if (k->real_closed) {
        void* out = malloc(ndt_elem_size(dta) * sz);
        ndb_ctx ctx = { in, out, dta, dta, k, sre, 0.0, arr_first };
        if (!nd_parallel_for(sz, ndb_map_chunk, &ctx)) { free(out); return NULL; }
        return expr_new_ndarray_like(arr, rank, dims, out, dta);
    }

    /* Real inputs, may escape (Log[b, x] on negatives): compute complex, narrow
     * to real dtype iff every element is real, else promote. */
    NDType dtcplx = ndt_as_complex(dta);
    void* tmp = malloc(ndt_elem_size(dtcplx) * sz);
    bool any_imag = false;
    for (size_t j = 0; j < sz; j++) {
        double re, im, orr, oii;
        ndt_get(in, j, dta, &re, &im);
        bool ok = arr_first ? k->cplx(re, 0.0, sre, 0.0, &orr, &oii)
                            : k->cplx(sre, 0.0, re, 0.0, &orr, &oii);
        if (!ok) { free(tmp); return NULL; }
        if (oii != 0.0) any_imag = true;
        ndt_set(tmp, j, dtcplx, orr, oii);
    }
    if (any_imag) return expr_new_ndarray_like(arr, rank, dims, tmp, dtcplx);
    void* out = malloc(ndt_elem_size(dta) * sz);
    for (size_t j = 0; j < sz; j++) {
        double re, im;
        ndt_get(tmp, j, dtcplx, &re, &im);
        ndt_set(out, j, dta, re, im);
    }
    free(tmp);
    return expr_new_ndarray_like(arr, rank, dims, out, dta);
}

/* Both operands are arrays. Same branch structure as the scalar-broadcast form
 * above, and deliberately so: the two must agree element for element, since a
 * user can write either and a packed pipeline chooses between them without
 * being asked. The only differences are that the second operand is read through
 * ndt_get at the same index instead of being hoisted out of the loop, and that
 * the result dtype is the PROMOTION of the two inputs rather than the one
 * input's.
 *
 * `ndarray_map_binary` declined this shape outright -- it required exactly one
 * array -- and so did eval.c's dispatch, with the consequence that a
 * two-argument kernel was reachable only when one argument happened to be a
 * scalar. `ArcTan[v, w]`, the numpy `arctan2` of two vectors, and `Beta[p, p]`
 * both fell straight through to the unevaluated call on the visible surface and
 * to one Expr per element on the packed one, despite NDKB_ArcTan and NDKB_Beta
 * having existed all along. Fifteen registered kernels were in the same
 * position, so this is the engine gap rather than fifteen head-level ones. */
Expr* ndarray_map_binary2(const Expr* a0, const Expr* a1, const NDBinaryKernel* k) {
    if (!k) return NULL;               /* `cplx` is tested per branch, as above */
    if (!is_ndarray(a0) || !is_ndarray(a1)) return NULL;

    int rank = a0->data.ndarray.rank;
    if (a1->data.ndarray.rank != rank) return NULL;
    const int64_t* dims = a0->data.ndarray.dims;
    for (int d = 0; d < rank; d++)
        if (a1->data.ndarray.dims[d] != dims[d]) return NULL;

    NDType dt0 = a0->data.ndarray.dtype, dt1 = a1->data.ndarray.dtype;
    if (dt0 == NDT_BOOL || dt1 == NDT_BOOL) return NULL;  /* not numeric; see ndarray_map_unary */
    const void* in0 = a0->data.ndarray.data;
    const void* in1 = a1->data.ndarray.data;
    size_t sz = ndarray_size(a0);
    /* A visible NDArray on either side wins the presentation, which is what
     * nd_present_src2 encodes -- Beta[packed, visible] is a visible answer. */
    const Expr* psrc = nd_present_src2(a0, a1);

    /* EXACT INTEGER arm, as in the scalar form: an exact answer that exists must
     * not be computed in double and rounded back. Both operands have to be
     * int64 for it; a mixed int64/real pair has an inexact answer and falls to
     * the guard below. Serial and all-or-nothing -- one element with no exact
     * int64 answer abandons the whole array so the List path answers in GMP. */
    if (k->to_int && dt0 == NDT_INT64 && dt1 == NDT_INT64 && k->to_int_i) {
        int64_t* out = (int64_t*)malloc(sizeof(int64_t) * (sz ? sz : 1));
        if (!out) return NULL;
        bool ok = true;
        for (size_t i = 0; i < sz && ok; i++)
            ok = k->to_int_i(ndt_get_i(in0, i, dt0), ndt_get_i(in1, i, dt1), &out[i]);
        if (!ok) { free(out); return NULL; }
        return expr_new_ndarray_like(psrc, rank, dims, out, NDT_INT64);
    }

    /* An int64 operand with no exact arm: every branch below reads through
     * ndt_get, which routes an int64 through `double` and is exact only to
     * 2^53. Decline, and ndarray_delist_and_reeval answers as the List does. */
    if (dt0 == NDT_INT64 || dt1 == NDT_INT64) return NULL;

    if (!k->cplx) return NULL;   /* integer-only kernel, or sentinel: degrade */

    NDType dtp = ndt_promote(dt0, dt1);

    /* Genuinely-complex output: either operand complex. */
    if (ndt_is_complex(dtp)) {
        void* out = malloc(ndt_elem_size(dtp) * (sz ? sz : 1));
        if (!out) return NULL;
        ndb2_ctx ctx = { in0, in1, out, dt0, dt1, dtp, k };
        if (!nd_parallel_for(sz, ndb2_map_chunk, &ctx)) { free(out); return NULL; }
        return expr_new_ndarray_like(psrc, rank, dims, out, dtp);
    }

    /* Real inputs, real-closed function: real output. */
    if (k->real_closed) {
        void* out = malloc(ndt_elem_size(dtp) * (sz ? sz : 1));
        if (!out) return NULL;
        ndb2_ctx ctx = { in0, in1, out, dt0, dt1, dtp, k };
        if (!nd_parallel_for(sz, ndb2_map_chunk, &ctx)) { free(out); return NULL; }
        return expr_new_ndarray_like(psrc, rank, dims, out, dtp);
    }

    /* Real inputs, may escape the real axis (Log[b, x] on a negative x):
     * compute complex, narrow back iff every element came out real. */
    NDType dtcplx = ndt_as_complex(dtp);
    void* tmp = malloc(ndt_elem_size(dtcplx) * (sz ? sz : 1));
    if (!tmp) return NULL;
    bool any_imag = false;
    for (size_t j = 0; j < sz; j++) {
        double are, bre, im, orr, oii;
        ndt_get(in0, j, dt0, &are, &im);
        ndt_get(in1, j, dt1, &bre, &im);
        if (!k->cplx(are, 0.0, bre, 0.0, &orr, &oii)) { free(tmp); return NULL; }
        if (oii != 0.0) any_imag = true;
        ndt_set(tmp, j, dtcplx, orr, oii);
    }
    if (any_imag) return expr_new_ndarray_like(psrc, rank, dims, tmp, dtcplx);
    void* out = malloc(ndt_elem_size(dtp) * (sz ? sz : 1));
    if (!out) { free(tmp); return NULL; }
    for (size_t j = 0; j < sz; j++) {
        double re, im;
        ndt_get(tmp, j, dtcplx, &re, &im);
        ndt_set(out, j, dtp, re, im);
    }
    free(tmp);
    return expr_new_ndarray_like(psrc, rank, dims, out, dtp);
}

Expr* ndarray_delist_and_reeval(const Expr* call) {
    if (!call || call->type != EXPR_FUNCTION) return NULL;
    size_t n = call->data.function.arg_count;
    Expr** args = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* a = call->data.function.args[i];
        args[i] = is_ndarray(a) ? ndarray_to_nested_list(a) : expr_copy(a);
    }
    Expr* rebuilt = expr_new_function(expr_copy(call->data.function.head), args, n);
    free(args);
    return eval_and_free(rebuilt);
}

/* An operand whose true answer no machine float buffer can hold: 5/2 *
 * Range[300] is a list of Rationals, Range[300]^(1/2) a list of radicals,
 * Range[300] + I a list of exact Complex, and {True, False} + {True, True} a
 * list of unevaluated Plus[...]. The elementwise fast paths decline those by
 * returning NULL -- and the arithmetic heads' usual degrade for a NULL ("warn
 * that an array met a symbolic operand, leave the call unevaluated") is exactly
 * wrong here, because the List path computes the answer perfectly well.
 *
 * Returns the re-evaluated call when some operand is an int64 OR bool buffer
 * (neither is a machine float, so both must go to the List path), else NULL so
 * the caller keeps its existing degrade. Borrows `call`. */
Expr* ndarray_int64_delist_retry(const Expr* call) {
    if (!call || call->type != EXPR_FUNCTION) return NULL;
    bool any = false;
    for (size_t i = 0; i < call->data.function.arg_count; i++) {
        const Expr* a = call->data.function.args[i];
        if (a && a->type == EXPR_NDARRAY &&
            (a->data.ndarray.dtype == NDT_INT64 ||
             a->data.ndarray.dtype == NDT_BOOL)) {
            any = true;
            break;
        }
    }
    return any ? ndarray_delist_and_reeval(call) : NULL;
}

/* A symbolic operand beside an INVISIBLE packed array: materialise and re-run.
 *
 * `NDArray[{1., 2.}] + x` warning and staying unevaluated is defensible — the
 * user named a purely numeric object. `{1., 2., 3.} + x` doing it is not: that
 * is an ordinary List, the List path threads it to `{1. + x, 2. + x, 3. + x}`,
 * and whether the value happens to be packed is supposed to be invisible.
 * Packing keys on LENGTH, so the same expression answered differently either
 * side of the threshold:
 *
 *     Range[1., 3.]   + x  ->  {1. + x, 2. + x, 3. + x}
 *     Range[1., 300.] + x  ->  unevaluated, with NDArray::sym
 *
 * This is the same shape as ndarray_int64_delist_retry above and sits beside
 * it at each call site: return the re-evaluated call when some argument is a
 * packed LIST (not a visible NDArray) and another is symbolic, else NULL so
 * the caller warns and declines exactly as before.
 *
 * Found on 2026-08-02 by the packed-vs-plain differential, while checking
 * whether Subtract and Divide could be marked packed-aware — they inherit
 * Plus/Times/Power's behaviour here, so the bug had to be fixed at the source
 * rather than propagated to two more heads. */
Expr* ndarray_symbolic_delist_retry(const Expr* call) {
    if (!call || call->type != EXPR_FUNCTION) return NULL;
    bool packed = false, sym = false, visible = false;
    for (size_t i = 0; i < call->data.function.arg_count; i++) {
        const Expr* a = call->data.function.args[i];
        if (!a) continue;
        if (a->type == EXPR_NDARRAY) {
            if (is_packed_list(a)) packed = true;
            else                   visible = true;
        } else if (!expr_is_numeric_like((Expr*)a)) {
            sym = true;
        }
    }
    /* One VISIBLE NDArray anywhere in the call vetoes the retry, even when a
     * packed List is also present.  The alternative was tried and is worse:
     * `NDArray[{1., 2., 3.}] + packedList + x` would then evaluate while
     * `NDArray[{1., 2., 3.}] + plainList + x` still declined, which is the
     * very divergence this function exists to remove, moved rather than
     * fixed.  A visible NDArray dominates a packed one everywhere else
     * (ndarray_result_presentation); it dominates here too. */
    return (packed && sym && !visible) ? ndarray_delist_and_reeval(call) : NULL;
}

/* Format an NDArray's shape as "{d0, d1, ...}" into buf. */
static void format_shape(const Expr* a, char* buf, size_t buflen) {
    size_t off = 0;
    int written = snprintf(buf + off, buflen - off, "{");
    if (written > 0) off += (size_t)written;
    for (int i = 0; i < a->data.ndarray.rank && off < buflen; i++) {
        written = snprintf(buf + off, buflen - off, "%s%lld",
                           i ? ", " : "", (long long)a->data.ndarray.dims[i]);
        if (written > 0) off += (size_t)written;
    }
    if (off < buflen) snprintf(buf + off, buflen - off, "}");
}

/* Emit `msg` to stderr at most once per evaluation. The fixed-point evaluator
 * visits the same node more than once per evaluate(), so dedup on the eval
 * clock (constant within one evaluate(), bumped between REPL inputs) plus a
 * hash of the message: an identical warning prints once per evaluation but
 * still fires on a later input. Respects the arithmetic-warning mute. */
void ndarray_warn_once(const char* msg) {
    if (arith_warnings_muted()) return;
    static uint64_t last_clock = 0;
    static uint64_t last_key = 0;
    uint64_t key = 1469598103934665603ULL; /* FNV-1a over the text */
    for (const char* c = msg; *c; c++)
        key = (key ^ (unsigned char)*c) * 1099511628211ULL;
    uint64_t clock = eval_clock_get();
    if (clock == last_clock && key == last_key) return;
    fputs(msg, stderr);
    last_clock = clock;
    last_key = key;
}

bool ndarray_warn_shape_mismatch(Expr** args, size_t n, const char* verb) {
    if (n < 2) return false;
    for (size_t i = 0; i < n; i++) {
        if (args[i]->type != EXPR_NDARRAY) return false; /* not a pure-NDArray op */
    }
    for (size_t i = 1; i < n; i++) {
        if (!same_shape(args[0], args[i])) {
            char s0[256], si[256], msg[600];
            format_shape(args[0], s0, sizeof s0);
            format_shape(args[i], si, sizeof si);
            snprintf(msg, sizeof msg,
                "NDArray::shape: NDArray objects of shapes %s and %s cannot "
                "be %s elementwise.\n", s0, si, verb);
            ndarray_warn_once(msg);
            return true;
        }
    }
    return false;
}

bool ndarray_warn_symbolic(Expr** args, size_t n, const char* verb) {
    bool has_nd = false, has_sym = false;
    for (size_t i = 0; i < n; i++) {
        if (args[i]->type == EXPR_NDARRAY) has_nd = true;
        else if (!expr_is_numeric_like(args[i])) has_sym = true;
    }
    if (!has_nd || !has_sym) return false;
    char msg[300];
    snprintf(msg, sizeof msg,
        "NDArray::sym: NDArray objects are purely numeric and cannot be %s "
        "with a symbolic expression; the expression is left unevaluated.\n",
        verb);
    ndarray_warn_once(msg);
    return true;
}

/* Structural shape probe independent of leaf numeric-ness: returns the depth
 * of `e` treating any non-List as a depth-0 leaf, writing the shape into dims,
 * or -1 if `e` is internally ragged (siblings of differing shape, or a mix of
 * List and non-List siblings). A bare non-List returns 0. Unlike probe_dims,
 * this does NOT reject symbolic leaves — so NDArray[{1, x}] (rectangular but
 * symbolic) is not ragged, while NDArray[{1, {2}}] is. */
static int shape_probe(const Expr* e, int64_t* dims, int max) {
    if (!head_is(e, SYM_List)) return 0;   /* leaf (any non-List) */
    if (max <= 0) return 0;                /* too deep to record; stop */
    int64_t len = (int64_t)e->data.function.arg_count;
    dims[0] = len;
    if (len == 0) return 1;                /* empty list: shape {0}, not ragged */

    int64_t sub[NDARRAY_MAX_RANK];
    int sub_rank = shape_probe(e->data.function.args[0], sub, max - 1);
    if (sub_rank < 0) return -1;           /* first child already ragged */
    for (int64_t i = 1; i < len; i++) {
        int64_t cur[NDARRAY_MAX_RANK];
        int cur_rank = shape_probe(e->data.function.args[i], cur, max - 1);
        if (cur_rank != sub_rank) return -1;               /* mixed depth */
        for (int j = 0; j < sub_rank; j++)
            if (cur[j] != sub[j]) return -1;               /* unequal dims */
    }
    for (int j = 0; j < sub_rank; j++) dims[j + 1] = sub[j];
    return sub_rank + 1;
}

/* True iff `list` is a List that is internally ragged (non-rectangular). */
static bool ndarray_is_ragged(const Expr* list) {
    if (!head_is(list, SYM_List)) return false;
    int64_t dims[NDARRAY_MAX_RANK];
    return shape_probe(list, dims, NDARRAY_MAX_RANK) < 0;
}

/* Recursively test whether a nested List contains any Complex[...] leaf, so the
 * default (unspecified) dtype can be inferred as complex rather than defaulting
 * to float64 and failing to pack. Recurses through Lists; any non-List node with
 * head Complex counts. */
static bool nested_list_has_complex(const Expr* e) {
    if (head_is(e, SYM_List)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (nested_list_has_complex(e->data.function.args[i])) return true;
        return false;
    }
    return head_is(e, SYM_Complex);
}

/* Strip a trailing `DataType -> "..."` option from res's args (rightmost wins),
 * modelled on options.c's extract_extension_option_full. Sets *dtype (default
 * float64), decrements *argc to the count of non-option args, and returns true
 * unless an explicit DataType option carried an unknown string (in which case
 * *bad is set true and the caller leaves the expression unevaluated). When no
 * DataType option is present, *explicit is set false so the caller can infer the
 * dtype from the data (e.g. complex leaves -> complex64). */
static bool extract_datatype_option(const Expr* res, size_t* argc,
                                    NDType* dtype, bool* bad, bool* explicit_dt) {
    *dtype = NDT_FLOAT64;
    *bad = false;
    size_t n = res->data.function.arg_count;
    bool seen = false;
    while (n > 0) {
        const Expr* opt = res->data.function.args[n - 1];
        if (opt && opt->type == EXPR_FUNCTION && opt->data.function.head
            && opt->data.function.head->type == EXPR_SYMBOL
            && (opt->data.function.head->data.symbol.name == SYM_Rule
                || opt->data.function.head->data.symbol.name == SYM_RuleDelayed)
            && opt->data.function.arg_count == 2) {
            const Expr* lhs = opt->data.function.args[0];
            const Expr* rhs = opt->data.function.args[1];
            if (lhs && lhs->type == EXPR_SYMBOL
                && lhs->data.symbol.name == SYM_DataType) {
                if (!seen) {
                    if (rhs && rhs->type == EXPR_STRING
                        && ndt_from_string(rhs->data.string, dtype)) {
                        /* ok */
                    } else {
                        *bad = true;
                    }
                    seen = true;
                }
                n--;
                continue;
            }
        }
        break;
    }
    *argc = n;
    *explicit_dt = seen;
    return true;
}

Expr* builtin_ndarray(Expr* res) {
    size_t argc;
    NDType dtype;
    bool bad, explicit_dt;
    extract_datatype_option(res, &argc, &dtype, &bad, &explicit_dt);
    if (bad) {
        ndarray_warn_once("NDArray::dtype: Unknown DataType; expected one of "
                          "\"float64\", \"float32\", \"complex64\", "
                          "\"complex32\", \"bool\".\n");
        return NULL;
    }
    if (argc != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    /* Idempotent on an NDArray argument: NDArray[array] is the array itself,
     * presented visibly (a transparent packed-list becomes a visible
     * NDArray[...]).  Without this branch ndarray_from_nested_list below returns
     * NULL for an EXPR_NDARRAY (head_is requires a List), so NDArray[NDArray[m]]
     * stayed unevaluated as a malformed rank-1 value.  Mirrors builtin_tondarray
     * (src/pack.c).  A round-trip through the List form is used only for an
     * explicit DataType re-cast. */
    if (is_ndarray(arg)) {
        bool same_dtype = !explicit_dt || arg->data.ndarray.dtype == dtype;
        if (same_dtype) {
            res->data.function.args[0] = NULL;   /* steal arg; evaluator frees res */
            arg = expr_unshare(arg);             /* private (refcount==1) before mutating */
            arg->data.ndarray.present_as = NDA_HEAD_NDARRAY;
            return arg;
        }
        Expr* list = ndarray_to_nested_list(arg);
        Expr* recast = ndarray_from_nested_list(list, dtype);
        expr_free(list);
        if (recast) recast->data.ndarray.present_as = NDA_HEAD_NDARRAY;
        return recast;   /* NULL leaves the call unevaluated */
    }
    /* No explicit DataType: infer complex64 when the data carries complex
     * leaves, else keep the float64 default. float32/complex32 require an
     * explicit option. */
    if (!explicit_dt && nested_list_has_complex(arg)) dtype = NDT_COMPLEX64;
    Expr* packed = ndarray_from_nested_list(arg, dtype);
    if (packed) return packed;
    /* A ragged (non-rectangular) list can never form an array, so reject it
     * with a warning rather than leaving a valid-looking NDArray[...] around.
     * A non-ragged but non-machine-precision list (e.g. a symbolic entry)
     * stays unevaluated silently — it may become packable after evaluation. */
    if (ndarray_is_ragged(arg)) {
        ndarray_warn_once("NDArray::ragged: The list is not rectangular "
                          "(ragged) and cannot form an NDArray.\n");
    }
    return NULL; /* leave unevaluated */
}

Expr* builtin_datatype(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* arg = res->data.function.args[0];
    if (!is_ndarray(arg)) return NULL;  /* stays symbolic on non-arrays */
    return expr_new_string(ndt_to_string(arg->data.ndarray.dtype));
}

Expr* builtin_ndarrayq(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    return expr_new_symbol(is_ndarray(res->data.function.args[0]) ? SYM_True : SYM_False);
}

/* PackedArrayQ tests the packed-List surface only (is_packed_list), so a visible
 * NDArray[...] gives False -- unlike NDArrayQ, which is True for both surfaces.
 * The head must be packed-aware (AWARE + INT64_OK in src/pack.c) or the gate
 * would materialise a packed argument into a plain List before this runs, making
 * it wrongly answer False. */
Expr* builtin_packedarrayq(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    return expr_new_symbol(is_packed_list(res->data.function.args[0]) ? SYM_True : SYM_False);
}

void ndarray_init(void) {
    symtab_add_builtin("NDArrayQ", builtin_ndarrayq);
    symtab_get_def("NDArrayQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("NDArrayQ",
        "NDArrayQ[expr]\n"
        "\tGives True if expr is an NDArray object, else False.");

    symtab_add_builtin("PackedArrayQ", builtin_packedarrayq);
    symtab_get_def("PackedArrayQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("PackedArrayQ",
        "PackedArrayQ[expr]\n"
        "\tGives True if expr is a packed array -- a List stored as a dense\n"
        "\tmachine-precision buffer -- else False. Provided under Mathematica's\n"
        "\tname; unlike NDArrayQ, a visible NDArray[...] object gives False.");

    symtab_add_builtin("NDArray", builtin_ndarray);
    symtab_get_def("NDArray")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("NDArray",
        "NDArray[nested_list]\n"
        "\tPacks a rectangular, machine-precision (Integer/Real) nested\n"
        "\tlist into a dense N-dimensional array (numpy ndarray style).\n"
        "\tVisibly distinct from List: Head, ListQ, and printing never\n"
        "\ttreat an NDArray as a List. Dimensions gives its shape,\n"
        "\tArrayDepth its rank, Length its leading-axis length. Builtins\n"
        "\tthat recognize NDArray (Dot, Plus, Times) use a fast C-level\n"
        "\tpath; results that would need a non-machine-precision entry\n"
        "\tauto-degrade to an ordinary nested List.\n"
        "NDArray[nested_list, DataType -> \"float32\"]\n"
        "\tPacks at the given element type: \"float64\" (default), \"float32\",\n"
        "\t\"complex64\", \"complex32\", or \"bool\" (a list of True/False;\n"
        "\t\"Boolean\" is accepted too). DataType[a] gives an array's type.\n"
        "\tA ragged (non-rectangular) list is rejected with an\n"
        "\tNDArray::ragged warning; an empty or non-machine-precision\n"
        "\tlist stays unevaluated.");
    /* Default option surfaced by Options[NDArray]. */
    {
        Expr* rule_args[2] = { expr_new_symbol(SYM_DataType),
                               expr_new_string("float64") };
        Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule), rule_args, 2);
        Expr* opts_args[1] = { rule };
        Expr* opts = expr_new_function(expr_new_symbol(SYM_List), opts_args, 1);
        symtab_set_options("NDArray", opts);  /* takes ownership */
    }

    symtab_add_builtin("DataType", builtin_datatype);
    symtab_get_def("DataType")->attributes |= ATTR_PROTECTED;
    /* Docstring lives in info.c (centralized convention). */
}
