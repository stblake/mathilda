/* NDArray — see ndarray.h for the design rationale. */

#include "ndarray.h"
#include "sym_names.h"
#include "attr.h"
#include "symtab.h"
#include "common.h"
#include "arithmetic.h"   /* arith_warnings_muted() */
#include "eval.h"         /* eval_clock_get() — dedup warnings within one eval */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

bool is_ndarray(const Expr* e) {
    return e && e->type == EXPR_NDARRAY;
}

size_t ndarray_size(const Expr* a) {
    size_t n = 1;
    for (int i = 0; i < a->data.ndarray.rank; i++) {
        n *= (size_t)a->data.ndarray.dims[i];
    }
    return n;
}

/* Machine-precision-safe leaf: Integer or Real only. BigInt/Rational/
 * Complex/MPFR/symbolic all risk silent precision loss or aren't numeric
 * at all, so they fail packing here and the caller falls back to List. */
static bool leaf_to_double(const Expr* e, double* out) {
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL) { *out = e->data.real; return true; }
    return false;
}

/* Recursive rank/shape probe, mirroring linalg/util.c's get_tensor_dims but
 * over `int64_t dims[NDARRAY_MAX_RANK]` and rejecting non-machine-precision
 * leaves up front (get_tensor_dims doesn't care about leaf type since the
 * symbolic path handles any leaf). Returns 0 for a non-List leaf, -1 for
 * jagged/non-numeric, else the rank. */
static int probe_dims(const Expr* e, int64_t* dims) {
    if (!head_is(e, SYM_List)) return 0;
    int64_t len = (int64_t)e->data.function.arg_count;
    dims[0] = len;
    if (len == 0) return -1; /* empty: nothing to pack, degrade to List */

    int sub_rank = probe_dims(e->data.function.args[0], dims + 1);
    if (sub_rank == 0) {
        double d;
        if (!leaf_to_double(e->data.function.args[0], &d)) return -1;
    } else if (sub_rank < 0) {
        return -1;
    }
    for (int64_t i = 1; i < len; i++) {
        int64_t cur[NDARRAY_MAX_RANK];
        int cur_rank = probe_dims(e->data.function.args[i], cur);
        if (cur_rank != sub_rank) return -1;
        if (cur_rank == 0) {
            double d;
            if (!leaf_to_double(e->data.function.args[i], &d)) return -1;
        }
        for (int j = 0; j < sub_rank; j++) {
            if (cur[j] != dims[j + 1]) return -1;
        }
    }
    return sub_rank + 1;
}

static void flatten_into(const Expr* e, double* flat, size_t* idx) {
    if (head_is(e, SYM_List)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            flatten_into(e->data.function.args[i], flat, idx);
        }
    } else {
        double d = 0.0;
        leaf_to_double(e, &d); /* already validated by probe_dims */
        flat[(*idx)++] = d;
    }
}

Expr* ndarray_from_nested_list(const Expr* list) {
    if (!head_is(list, SYM_List)) return NULL;
    int64_t dims[NDARRAY_MAX_RANK];
    int rank = probe_dims(list, dims);
    if (rank <= 0 || rank >= NDARRAY_MAX_RANK) return NULL;

    size_t n = 1;
    for (int i = 0; i < rank; i++) n *= (size_t)dims[i];
    double* data = malloc(sizeof(double) * n);
    if (!data) return NULL;
    size_t idx = 0;
    flatten_into(list, data, &idx);

    return expr_new_ndarray(rank, dims, data); /* takes ownership of data */
}

/* Recursive rebuild of one axis-level of the nested List tree, mirroring
 * linalg/dot.c's build_tensor but reading from a flat double* rather than
 * an Expr** buffer. `idx` is the running cursor into `data`. */
static Expr* rebuild_level(const int64_t* dims, int rank, int level,
                            const double* data, size_t* idx) {
    if (level == rank) {
        return expr_new_real(data[(*idx)++]);
    }
    int64_t len = dims[level];
    Expr** args = malloc(sizeof(Expr*) * (size_t)len);
    for (int64_t i = 0; i < len; i++) {
        args[i] = rebuild_level(dims, rank, level + 1, data, idx);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), args, (size_t)len);
    free(args);
    return result;
}

Expr* ndarray_to_nested_list(const Expr* a) {
    size_t idx = 0;
    return rebuild_level(a->data.ndarray.dims, a->data.ndarray.rank, 0,
                          a->data.ndarray.data, &idx);
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

    double* out = malloc(sizeof(double) * (size_t)(Ra * Sb));
    for (int64_t r = 0; r < Ra; r++) {
        for (int64_t s = 0; s < Sb; s++) {
            double sum = 0.0;
            for (int64_t k = 0; k < K; k++) {
                sum += a->data.ndarray.data[r * K + k] * b->data.ndarray.data[k * Sb + s];
            }
            out[r * Sb + s] = sum;
        }
    }

    int rankC = rankA + rankB - 2;
    if (rankC == 0) {
        double scalar = out[0];
        free(out);
        return expr_new_real(scalar);
    }
    int64_t dimsC[2];
    if (rankC == 1) dimsC[0] = (rankA == 2) ? Ra : Sb;
    else { dimsC[0] = Ra; dimsC[1] = Sb; }
    return expr_new_ndarray(rankC, dimsC, out); /* takes ownership of out */
}

static bool same_shape(const Expr* a, const Expr* b) {
    if (a->data.ndarray.rank != b->data.ndarray.rank) return false;
    for (int i = 0; i < a->data.ndarray.rank; i++) {
        if (a->data.ndarray.dims[i] != b->data.ndarray.dims[i]) return false;
    }
    return true;
}

Expr* ndarray_elementwise(Expr** args, size_t n, bool is_plus) {
    if (n == 0) return NULL;
    for (size_t i = 0; i < n; i++) {
        if (args[i]->type != EXPR_NDARRAY) return NULL;
    }
    const Expr* first = args[0];
    for (size_t i = 1; i < n; i++) {
        if (!same_shape(first, args[i])) return NULL;
    }

    size_t sz = ndarray_size(first);
    double* out = malloc(sizeof(double) * sz);
    for (size_t k = 0; k < sz; k++) {
        double acc = is_plus ? 0.0 : 1.0;
        for (size_t i = 0; i < n; i++) {
            double v = args[i]->data.ndarray.data[k];
            acc = is_plus ? acc + v : acc * v;
        }
        out[k] = acc;
    }
    return expr_new_ndarray(first->data.ndarray.rank, first->data.ndarray.dims, out); /* takes ownership of out */
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
static void ndarray_warn_once(const char* msg) {
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

Expr* builtin_ndarray(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    Expr* packed = ndarray_from_nested_list(arg);
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

Expr* builtin_ndarrayq(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    return expr_new_symbol(is_ndarray(res->data.function.args[0]) ? SYM_True : SYM_False);
}

void ndarray_init(void) {
    symtab_add_builtin("NDArrayQ", builtin_ndarrayq);
    symtab_get_def("NDArrayQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("NDArrayQ",
        "NDArrayQ[expr]\n"
        "\tGives True if expr is an NDArray value, else False.");

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
        "\tA ragged (non-rectangular) list is rejected with an\n"
        "\tNDArray::ragged warning; an empty or non-machine-precision\n"
        "\tlist stays unevaluated.");
}
