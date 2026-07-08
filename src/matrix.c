/* Matrix — see matrix.h for the design rationale. */

#include "matrix.h"
#include "sym_names.h"
#include "attr.h"
#include "symtab.h"
#include "common.h"
#include <stdlib.h>
#include <string.h>

bool is_matrix(const Expr* e) {
    return e && e->type == EXPR_MATRIX;
}

size_t matrix_size(const Expr* m) {
    size_t n = 1;
    for (int i = 0; i < m->data.matrix.rank; i++) {
        n *= (size_t)m->data.matrix.dims[i];
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
 * over `int64_t dims[MATRIX_MAX_RANK]` and rejecting non-machine-precision
 * leaves up front (get_tensor_dims doesn't care about leaf type since the
 * symbolic path handles any leaf). Returns 0 for a non-List leaf, -1 for
 * jagged/non-numeric, else the rank. */
#define MATRIX_MAX_RANK 64
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
        int64_t cur[MATRIX_MAX_RANK];
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

Expr* matrix_from_nested_list(const Expr* list) {
    if (!head_is(list, SYM_List)) return NULL;
    int64_t dims[MATRIX_MAX_RANK];
    int rank = probe_dims(list, dims);
    if (rank <= 0 || rank >= MATRIX_MAX_RANK) return NULL;

    size_t n = 1;
    for (int i = 0; i < rank; i++) n *= (size_t)dims[i];
    double* data = malloc(sizeof(double) * n);
    if (!data) return NULL;
    size_t idx = 0;
    flatten_into(list, data, &idx);

    return expr_new_matrix(rank, dims, data); /* takes ownership of data */
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

Expr* matrix_to_nested_list(const Expr* m) {
    size_t idx = 0;
    return rebuild_level(m->data.matrix.dims, m->data.matrix.rank, 0,
                          m->data.matrix.data, &idx);
}

Expr* matrix_dot2(const Expr* a, const Expr* b, bool* shape_error) {
    int rankA = a->data.matrix.rank;
    int rankB = b->data.matrix.rank;
    if (rankA < 1 || rankA > 2 || rankB < 1 || rankB > 2) return NULL;

    const int64_t* dimsA = a->data.matrix.dims;
    const int64_t* dimsB = b->data.matrix.dims;
    int64_t K = dimsA[rankA - 1];
    if (K != dimsB[0]) { *shape_error = true; return NULL; }

    int64_t Ra = (rankA == 2) ? dimsA[0] : 1;
    int64_t Sb = (rankB == 2) ? dimsB[1] : 1;

    double* out = malloc(sizeof(double) * (size_t)(Ra * Sb));
    for (int64_t r = 0; r < Ra; r++) {
        for (int64_t s = 0; s < Sb; s++) {
            double sum = 0.0;
            for (int64_t k = 0; k < K; k++) {
                sum += a->data.matrix.data[r * K + k] * b->data.matrix.data[k * Sb + s];
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
    return expr_new_matrix(rankC, dimsC, out); /* takes ownership of out */
}

static bool same_shape(const Expr* a, const Expr* b) {
    if (a->data.matrix.rank != b->data.matrix.rank) return false;
    for (int i = 0; i < a->data.matrix.rank; i++) {
        if (a->data.matrix.dims[i] != b->data.matrix.dims[i]) return false;
    }
    return true;
}

Expr* matrix_elementwise(Expr** args, size_t n, bool is_plus) {
    if (n == 0) return NULL;
    for (size_t i = 0; i < n; i++) {
        if (args[i]->type != EXPR_MATRIX) return NULL;
    }
    const Expr* first = args[0];
    for (size_t i = 1; i < n; i++) {
        if (!same_shape(first, args[i])) return NULL;
    }

    size_t sz = matrix_size(first);
    double* out = malloc(sizeof(double) * sz);
    for (size_t k = 0; k < sz; k++) {
        double acc = is_plus ? 0.0 : 1.0;
        for (size_t i = 0; i < n; i++) {
            double v = args[i]->data.matrix.data[k];
            acc = is_plus ? acc + v : acc * v;
        }
        out[k] = acc;
    }
    return expr_new_matrix(first->data.matrix.rank, first->data.matrix.dims, out); /* takes ownership of out */
}

Expr* builtin_matrix(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    return matrix_from_nested_list(arg); /* NULL: leave unevaluated */
}

void matrix_init(void) {
    symtab_add_builtin("Matrix", builtin_matrix);
    symtab_get_def("Matrix")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Matrix",
        "Matrix[nested_list]\n"
        "\tPacks a rectangular, machine-precision (Integer/Real) nested\n"
        "\tlist into a dense ndarray. Visibly distinct from List: Head,\n"
        "\tListQ, and printing never treat a Matrix as a List. Builtins\n"
        "\tthat recognize Matrix (Dot, Plus, Times) use a fast C-level\n"
        "\tpath; results that would need a non-machine-precision entry\n"
        "\tauto-degrade to an ordinary nested List.\n"
        "\tGives $Failed-free unevaluated output if the list is jagged,\n"
        "\tempty, or contains a non-machine-precision entry.");
}
