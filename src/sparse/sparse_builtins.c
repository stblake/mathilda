/*
 * sparse_builtins.c -- the SparseArray builtins and their registration:
 * SparseArray, ArrayRules, ArrayQ, ArrayDepth, the s["prop"] properties, and
 * Band (inert). Normal, Part, Length, Dimensions, VectorQ/MatrixQ, N and the
 * Listable heads reach SparseArrays through their own modules' sparse
 * branches; see sparse.h for the representation.
 */

#include "sparse.h"
#include "sparse_internal.h"
#include "sparse_la.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "message.h"
#include "ndarray.h"
#include "print.h"
#include "pack.h"
#include <stdlib.h>
#include <string.h>

Expr* builtin_sparsearray(Expr* res);   /* sparse_build.c */

static bool head_is_sym(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym;
}

static Expr* mk_list(Expr** items, size_t n) {
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, n);
    free(items);
    return out;
}

static Expr** items_alloc(size_t n) { return malloc(sizeof(Expr*) * (n ? n : 1)); }

/* The 1-based position of entry k (row r) as {i1, ..., ir}. */
static Expr* entry_position(const SparseArray* sa, int64_t r, int64_t k) {
    Expr** it = items_alloc((size_t)sa->rank);
    int q = 0;
    if (sa->rank >= 2) it[q++] = expr_new_integer(r + 1);
    for (int j = 0; j < sa->pat->width; j++) it[q++] = expr_new_integer(sp_col(sa->pat, k, j) + 1);
    return mk_list(it, (size_t)sa->rank);
}

/* ---- ArrayRules ---------------------------------------------------------- */

/* {pos1 -> v1, ..., {_, ..., _} -> bg}. */
static Expr* array_rules(const SparseArray* sa) {
    const SparsePattern* p = sa->pat;
    size_t n = (size_t)p->nnz;
    Expr** rules = items_alloc(n + 1);
    size_t q = 0;
    for (int64_t r = 0; r < p->nrows; r++)
        for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1]; k++) {
            Expr* ra[2] = { entry_position(sa, r, k), sp_value_expr(sa, k) };
            rules[q++] = expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
        }
    Expr** blanks = items_alloc((size_t)sa->rank);
    for (int k = 0; k < sa->rank; k++) blanks[k] = expr_new_function(expr_new_symbol(SYM_Blank), NULL, 0);
    Expr* ra[2] = { mk_list(blanks, (size_t)sa->rank), expr_copy(sa->bg) };
    rules[q++] = expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
    return mk_list(rules, q);
}

/* SparseArray[args...] evaluated: the conversion every non-sparse argument of
 * ArrayRules goes through. */
static Expr* to_sparse(Expr* const* args, size_t n) {
    Expr** a = items_alloc(n);
    for (size_t i = 0; i < n; i++) a[i] = expr_copy(args[i]);
    Expr* call = expr_new_function(expr_new_symbol(SYM_SparseArray), a, n);
    free(a);
    Expr* r = evaluate(call);
    expr_free(call);
    if (r && r->type != EXPR_SPARSE) { expr_free(r); return NULL; }
    return r;
}

static Expr* builtin_arrayrules(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t na = res->data.function.arg_count;
    if (na < 1 || na > 2) return NULL;
    Expr* a = res->data.function.args[0];
    Expr* s;
    if (na == 1 && a->type == EXPR_SPARSE) {
        s = expr_copy(a);
    } else if (a->type == EXPR_SPARSE) {
        /* ArrayRules[s, val]: the rules with val as the default */
        Expr* dims = NULL;
        {
            const SparseArray* sa = a->data.sparse;
            Expr** d = items_alloc((size_t)sa->rank);
            for (int k = 0; k < sa->rank; k++) d[k] = expr_new_integer(sa->dims[k]);
            dims = mk_list(d, (size_t)sa->rank);
        }
        Expr* args[3] = { a, dims, res->data.function.args[1] };
        s = to_sparse(args, 3);
        expr_free(dims);
    } else if (head_is_sym(a, SYM_List) || a->type == EXPR_NDARRAY) {
        if (na == 2) {
            /* ArrayRules[list, val] */
            Expr* sp = to_sparse(&a, 1);
            if (!sp) return NULL;
            const SparseArray* sa = sp->data.sparse;
            Expr** d = items_alloc((size_t)sa->rank);
            for (int k = 0; k < sa->rank; k++) d[k] = expr_new_integer(sa->dims[k]);
            Expr* dims = mk_list(d, (size_t)sa->rank);
            Expr* args[3] = { a, dims, res->data.function.args[1] };
            s = to_sparse(args, 3);
            expr_free(dims);
            expr_free(sp);
        } else {
            s = to_sparse(&a, 1);
        }
    } else {
        return NULL;
    }
    if (!s) return NULL;
    Expr* out = array_rules(s->data.sparse);
    expr_free(s);
    return out;
}

/* ---- ArrayQ / ArrayDepth ------------------------------------------------- */

/* Full-depth rectangular shape: rank and whether the structure is an array. */
static bool array_shape(const Expr* e, int* rank, int64_t* dims, int cap) {
    if (e->type == EXPR_SPARSE) {
        const SparseArray* sa = e->data.sparse;
        *rank = sa->rank;
        for (int k = 0; k < sa->rank && k < cap; k++) dims[k] = sa->dims[k];
        return true;
    }
    if (e->type == EXPR_NDARRAY) {
        *rank = e->data.ndarray.rank;
        for (int k = 0; k < *rank && k < cap; k++) dims[k] = e->data.ndarray.dims[k];
        return true;
    }
    if (!head_is_sym(e, SYM_List)) { *rank = 0; return false; }
    size_t n = e->data.function.arg_count;
    if (n == 0) { *rank = 1; dims[0] = 0; return true; }
    int r0 = 0;
    int64_t d0[SP_MAX_RANK];
    bool first_arr = array_shape(e->data.function.args[0], &r0, d0, SP_MAX_RANK);
    if (!first_arr) r0 = 0;
    for (size_t i = 1; i < n; i++) {
        int ri = 0;
        int64_t di[SP_MAX_RANK];
        bool arr = array_shape(e->data.function.args[i], &ri, di, SP_MAX_RANK);
        if (arr != first_arr) return false;
        if (!arr) continue;
        if (ri != r0) return false;
        for (int k = 0; k < ri && k < SP_MAX_RANK; k++) if (di[k] != d0[k]) return false;
    }
    if (1 + r0 > cap || 1 + r0 > SP_MAX_RANK) return false;
    *rank = 1 + r0;
    dims[0] = (int64_t)n;
    for (int k = 0; k < r0; k++) dims[k + 1] = d0[k];
    return true;
}

static bool test_true(const Expr* test, Expr* arg /* adopted */) {
    Expr* a1[1] = { arg };
    Expr* c1 = expr_new_function(expr_copy((Expr*)test), a1, 1);
    Expr* r1 = evaluate(c1);
    expr_free(c1);
    bool ok = r1 && r1->type == EXPR_SYMBOL && r1->data.symbol.name == SYM_True;
    if (r1) expr_free(r1);
    return ok;
}

bool sp_elements_pass(const Expr* s, const Expr* test) {
    const SparseArray* sa = sp_of(s);
    if (!sa) return false;
    /* the background stands for every implicit element, if there is one */
    int64_t total = sp_total(sa->rank, sa->dims);
    bool ok = total == sa->pat->nnz || test_true(test, expr_copy(sa->bg));
    for (int64_t k = 0; ok && k < sa->pat->nnz; k++) ok = test_true(test, sp_value_expr(sa, k));
    return ok;
}

/* Every leaf of a (verified) array passes `test`? */
static bool leaves_pass(const Expr* e, const Expr* test) {
    if (e->type == EXPR_SPARSE) return sp_elements_pass(e, test);
    if (e->type == EXPR_NDARRAY) {
        size_t n = ndarray_size(e);
        for (size_t k = 0; k < n; k++)
            if (!test_true(test, ndarray_element_to_expr(e, k))) return false;
        return true;
    }
    if (head_is_sym(e, SYM_List)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!leaves_pass(e->data.function.args[i], test)) return false;
        return true;
    }
    return test_true(test, expr_copy((Expr*)e));
}

/* ArrayQ[expr], ArrayQ[expr, patt], ArrayQ[expr, patt, test]: a full
 * rectangular array (every leaf at the same depth), whose depth matches patt
 * and whose every element passes test. */
static Expr* builtin_arrayq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t na = res->data.function.arg_count;
    if (na < 1 || na > 3) return NULL;
    const Expr* e = res->data.function.args[0];
    int rank = 0;
    int64_t dims[SP_MAX_RANK];
    bool arr = array_shape(e, &rank, dims, SP_MAX_RANK);
    if (!arr) return expr_new_symbol(SYM_False);
    if (na >= 2) {
        Expr* depth = expr_new_integer(rank);
        Expr* margs[2] = { depth, expr_copy(res->data.function.args[1]) };
        Expr* mq = expr_new_function(expr_new_symbol(SYM_MatchQ), margs, 2);
        Expr* r = evaluate(mq);
        expr_free(mq);
        bool ok = r && r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_True;
        if (r) expr_free(r);
        if (!ok) return expr_new_symbol(SYM_False);
    }
    if (na == 3 && !leaves_pass(e, res->data.function.args[2])) return expr_new_symbol(SYM_False);
    return expr_new_symbol(SYM_True);
}

/* ArrayDepth[expr]: the depth to which expr is a full array (0 for an atom). */
static Expr* builtin_arraydepth(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    const Expr* e = res->data.function.args[0];
    if (e->type == EXPR_SPARSE) return expr_new_integer(e->data.sparse->rank);
    if (e->type == EXPR_NDARRAY) return expr_new_integer(e->data.ndarray.rank);
    /* Descend while every element is a List of one common length. */
    int64_t depth = 0;
    /* The prefix shape of a possibly ragged nest, as Dimensions computes it. */
    Expr* dcall_args[1] = { expr_copy((Expr*)e) };
    Expr* dcall = expr_new_function(expr_new_symbol(intern_symbol("Dimensions")), dcall_args, 1);
    Expr* d = evaluate(dcall);
    expr_free(dcall);
    if (d && head_is_sym(d, SYM_List)) depth = (int64_t)d->data.function.arg_count;
    else if (d && d->type == EXPR_NDARRAY) depth = d->data.ndarray.dims[0];
    if (d) expr_free(d);
    return expr_new_integer(depth);
}

/* ---- properties ---------------------------------------------------------- */

static const char* const PROPERTIES[] = {
    "AdjacencyLists", "Background", "BandWidth", "ColumnIndices", "Density",
    "ExplicitLength", "ExplicitPositions", "ExplicitValues", "ImplicitValue",
    "NonzeroPositions", "NonzeroValues", "Properties", "RowPointers"
};

/* The bulk properties hand back packed buffers when packing accepts them
 * (ndbuild_open: above the threshold, packing enabled), boxed Lists
 * otherwise -- the same value either way. */
static Expr* prop_positions(const SparseArray* sa) {
    const SparsePattern* p = sa->pat;
    {
        int64_t dims[2] = { p->nnz, sa->rank };
        void* buf = NULL;
        Expr* node = p->nnz > 0 ? ndbuild_open(2, dims, NDT_INT64, &buf) : NULL;
        if (node) {
            int64_t* b = buf;
            int64_t q = 0;
            for (int64_t r = 0; r < p->nrows; r++)
                for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1]; k++, q++) {
                    int c = 0;
                    if (sa->rank >= 2) b[q * sa->rank + c++] = r + 1;
                    for (int j = 0; j < p->width; j++) b[q * sa->rank + c++] = sp_col(p, k, j) + 1;
                }
            return node;
        }
    }
    Expr** it = items_alloc((size_t)p->nnz);
    size_t q = 0;
    for (int64_t r = 0; r < p->nrows; r++)
        for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1]; k++) it[q++] = entry_position(sa, r, k);
    return mk_list(it, q);
}

static Expr* prop_values(const SparseArray* sa) {
    int64_t n = sa->pat->nnz;
    if (sa->kind != SPV_EXPR && n > 0) {
        NDType dt = sa->kind == SPV_F64 ? NDT_FLOAT64 : sa->kind == SPV_INT64 ? NDT_INT64 : NDT_COMPLEX64;
        void* buf = NULL;
        Expr* node = ndbuild_open(1, &n, dt, &buf);
        if (node) { memcpy(buf, sa->vals, sp_kind_size(sa->kind) * (size_t)n); return node; }
    }
    Expr** it = items_alloc((size_t)n);
    for (int64_t k = 0; k < n; k++) it[k] = sp_value_expr(sa, k);
    return mk_list(it, (size_t)n);
}

static Expr* sp_property(const SparseArray* sa, const char* name) {
    const SparsePattern* p = sa->pat;
    if (!strcmp(name, "Background") || !strcmp(name, "ImplicitValue")) return expr_copy(sa->bg);
    if (!strcmp(name, "NonzeroValues") || !strcmp(name, "ExplicitValues")) return prop_values(sa);
    if (!strcmp(name, "NonzeroPositions") || !strcmp(name, "ExplicitPositions")) return prop_positions(sa);
    if (!strcmp(name, "ExplicitLength")) return expr_new_integer(p->nnz);
    if (!strcmp(name, "Density")) {
        int64_t t = sp_total(sa->rank, sa->dims);
        return expr_new_real(t > 0 ? (double)p->nnz / (double)t : 0.0);
    }
    if (!strcmp(name, "RowPointers")) {
        int64_t len = p->nrows + 1;
        void* buf = NULL;
        Expr* node = ndbuild_open(1, &len, NDT_INT64, &buf);
        if (node) { memcpy(buf, p->rowptr, sizeof(int64_t) * (size_t)len); return node; }
        Expr** it = items_alloc((size_t)(p->nrows + 1));
        for (int64_t r = 0; r <= p->nrows; r++) it[r] = expr_new_integer(p->rowptr[r]);
        return mk_list(it, (size_t)(p->nrows + 1));
    }
    if (!strcmp(name, "ColumnIndices")) {
        int64_t dims[2] = { p->nnz, p->width };
        void* buf = NULL;
        Expr* node = p->nnz > 0 ? ndbuild_open(2, dims, NDT_INT64, &buf) : NULL;
        if (node) {
            int64_t* b = buf;
            for (int64_t k = 0; k < p->nnz; k++)
                for (int j = 0; j < p->width; j++) b[k * p->width + j] = sp_col(p, k, j) + 1;
            return node;
        }
        Expr** it = items_alloc((size_t)p->nnz);
        for (int64_t k = 0; k < p->nnz; k++) {
            Expr** t = items_alloc((size_t)p->width);
            for (int j = 0; j < p->width; j++) t[j] = expr_new_integer(sp_col(p, k, j) + 1);
            it[k] = mk_list(t, (size_t)p->width);
        }
        return mk_list(it, (size_t)p->nnz);
    }
    if (!strcmp(name, "AdjacencyLists")) {
        if (sa->rank != 2) return NULL;
        Expr** rows = items_alloc((size_t)p->nrows);
        for (int64_t r = 0; r < p->nrows; r++) {
            int64_t lo = p->rowptr[r], hi = p->rowptr[r + 1];
            Expr** c = items_alloc((size_t)(hi - lo));
            for (int64_t k = lo; k < hi; k++) c[k - lo] = expr_new_integer(sp_col(p, k, 0) + 1);
            rows[r] = mk_list(c, (size_t)(hi - lo));
        }
        return mk_list(rows, (size_t)p->nrows);
    }
    if (!strcmp(name, "BandWidth")) {
        /* {lower, upper}: how far the explicit entries reach below and above
         * the main diagonal of a matrix. */
        if (sa->rank != 2) return NULL;
        int64_t lo = 0, up = 0;
        for (int64_t r = 0; r < p->nrows; r++)
            for (int64_t k = p->rowptr[r]; k < p->rowptr[r + 1]; k++) {
                int64_t d = sp_col(p, k, 0) - r;
                if (d > up) up = d;
                if (-d > lo) lo = -d;
            }
        Expr** it = items_alloc(2);
        it[0] = expr_new_integer(lo);
        it[1] = expr_new_integer(up);
        return mk_list(it, 2);
    }
    if (!strcmp(name, "Properties")) {
        size_t n = sizeof PROPERTIES / sizeof PROPERTIES[0];
        Expr** it = items_alloc(n);
        for (size_t i = 0; i < n; i++) it[i] = expr_new_string(PROPERTIES[i]);
        return mk_list(it, n);
    }
    return NULL;
}

Expr* sp_apply(const Expr* s, Expr* const* args, size_t n) {
    const SparseArray* sa = sp_of(s);
    if (!sa || n != 1 || args[0]->type != EXPR_STRING) return NULL;
    Expr* r = sp_property(sa, args[0]->data.string);
    if (!r) {
        char* q = expr_to_string_outputform(args[0]);
        mth_message("SparseArray", "nomthd", "There is no method %s for SparseArray objects.", q);
        free(q);
    }
    return r;
}

/* ---- registration -------------------------------------------------------- */

void sparse_init(void) {
    /* Packed-awareness (and int64 exactness) for these heads is declared in
     * src/pack.c's AWARE / INT64_OK lists, where the audits read it. */
    symtab_add_builtin(SYM_SparseArray, builtin_sparsearray);
    symtab_get_def(SYM_SparseArray)->attributes |= ATTR_PROTECTED;

    symtab_add_builtin(SYM_ArrayRules, builtin_arrayrules);
    symtab_get_def(SYM_ArrayRules)->attributes |= ATTR_PROTECTED;

    symtab_add_builtin(SYM_ArrayQ, builtin_arrayq);
    symtab_get_def(SYM_ArrayQ)->attributes |= ATTR_PROTECTED;

    symtab_add_builtin(SYM_ArrayDepth, builtin_arraydepth);
    symtab_get_def(SYM_ArrayDepth)->attributes |= ATTR_PROTECTED;

    /* Band[...] is an inert specification for SparseArray. */
    symtab_get_def(SYM_Band)->attributes |= ATTR_PROTECTED;

    sparse_gate_init();
    sparse_reduce_init();
    sparse_struct_init();
    sparse_dot_init();
    sparse_func_init();
    sparse_la_init();
}
