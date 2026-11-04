/*
 * src/linalg/la_dense.c -- shared helpers for the dense matrix heads.
 * See la_dense.h for the contract of each function.
 */

#include "la_dense.h"
#include "common.h"
#include "eval.h"
#include "ndarray.h"
#include "print.h"
#include "sym_names.h"
#include "numeric.h"

#include <stdlib.h>
#include <string.h>

LaPrecClass la_prec_class(const Expr* e, long* bits) {
    if (bits) *bits = 53;
    if (is_ndarray(e)) return LA_MACHINE;
    CommonInexactInfo info = common_scan_inexact(e);
    if (!info.has_inexact) return LA_EXACT;
    if (info.min_bits > 53) {
        if (bits) *bits = info.min_bits;
        return LA_MPFR;
    }
    return LA_MACHINE;
}

static bool la_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && head_is(e, SYM_List);
}

bool la_matrix_dims(const Expr* e, int* rows, int* cols) {
    if (!e) return false;
    if (is_ndarray(e)) {
        if (e->data.ndarray.rank != 2) return false;
        *rows = (int)e->data.ndarray.dims[0];
        *cols = (int)e->data.ndarray.dims[1];
        return true;
    }
    if (!la_is_list(e)) return false;
    size_t r = e->data.function.arg_count;
    if (r == 0) return false;
    size_t c = 0;
    for (size_t i = 0; i < r; i++) {
        const Expr* row = e->data.function.args[i];
        if (is_ndarray(row)) {
            if (row->data.ndarray.rank != 1) return false;
            /* a row held as a packed vector: unusual for a List-of-rows,
             * accept its length */
            size_t len = (size_t)row->data.ndarray.dims[0];
            if (i == 0) c = len; else if (len != c) return false;
            continue;
        }
        if (!la_is_list(row)) return false;
        size_t len = row->data.function.arg_count;
        if (i == 0) c = len; else if (len != c) return false;
        for (size_t j = 0; j < len; j++)
            if (la_is_list(row->data.function.args[j])) return false;  /* rank > 2 */
    }
    *rows = (int)r;
    *cols = (int)c;
    return true;
}

bool la_looks_like_matrix(const Expr* e) {
    if (!e) return false;
    if (is_ndarray(e)) return e->data.ndarray.rank >= 2;
    if (!la_is_list(e)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!la_is_list(e->data.function.args[i])
            && !is_ndarray(e->data.function.args[i])) return false;
    return true;
}

Expr* la_boxed_copy(const Expr* e) {
    if (is_ndarray(e)) return ndarray_to_nested_list(e);
    if (la_is_list(e)) {
        size_t n = e->data.function.arg_count;
        bool any_nd = false;
        for (size_t i = 0; i < n; i++)
            if (is_ndarray(e->data.function.args[i])) { any_nd = true; break; }
        if (!any_nd) return expr_copy((Expr*)e);
        Expr** a = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t i = 0; i < n; i++) a[i] = la_boxed_copy(e->data.function.args[i]);
        Expr* out = expr_new_function(expr_new_symbol(SYM_List), a, n);
        free(a);
        return out;
    }
    return expr_copy((Expr*)e);
}

bool la_all_numeric(const Expr* e) {
    if (!e) return false;
    if (is_ndarray(e)) return true;
    if (la_is_list(e)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!la_all_numeric(e->data.function.args[i])) return false;
        return true;
    }
    return expr_is_numeric_quantity(e);
}

bool la_has_complex(const Expr* e) {
    if (!e) return false;
    if (is_ndarray(e)) return ndt_is_complex(e->data.ndarray.dtype);
    if (e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Complex)) return true;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_I) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (la_has_complex(e->data.function.args[i])) return true;
    return false;
}

Expr* la_call1(const char* head, Expr* a) {
    return expr_new_function(expr_new_symbol(head), (Expr*[]){ a }, 1);
}
Expr* la_call2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(head), (Expr*[]){ a, b }, 2);
}
Expr* la_call3(const char* head, Expr* a, Expr* b, Expr* c) {
    return expr_new_function(expr_new_symbol(head), (Expr*[]){ a, b, c }, 3);
}
Expr* la_ev1(const char* head, Expr* a) { return eval_and_free(la_call1(head, a)); }
Expr* la_ev2(const char* head, Expr* a, Expr* b) {
    return eval_and_free(la_call2(head, a, b));
}

Expr* la_empty_list(void) {
    return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
}

Expr* la_matrix_from(Expr** cells, int rows, int cols) {
    Expr** rv = (Expr**)malloc(sizeof(Expr*) * (size_t)(rows > 0 ? rows : 1));
    for (int i = 0; i < rows; i++) {
        Expr** cv = (Expr**)malloc(sizeof(Expr*) * (size_t)(cols > 0 ? cols : 1));
        for (int j = 0; j < cols; j++) cv[j] = cells[(size_t)i * (size_t)cols + (size_t)j];
        rv[i] = expr_new_function(expr_new_symbol(SYM_List), cv, (size_t)cols);
        free(cv);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rv, (size_t)rows);
    free(rv);
    return out;
}

char* la_str(const Expr* e) {
    char* s = expr_to_string((Expr*)e);
    if (!s) {
        s = (char*)malloc(2);
        if (s) { s[0] = '?'; s[1] = '\0'; }
    }
    return s;
}

bool la_parse_take_spec(const Expr* spec, long long* k, bool* upto) {
    if (!spec) return false;
    if (spec->type == EXPR_INTEGER) {
        *k = (long long)spec->data.integer;
        *upto = false;
        return true;
    }
    if (spec->type == EXPR_FUNCTION && head_is(spec, SYM_UpTo)
        && spec->data.function.arg_count == 1
        && spec->data.function.args[0]->type == EXPR_INTEGER
        && spec->data.function.args[0]->data.integer >= 0) {
        *k = (long long)spec->data.function.args[0]->data.integer;
        *upto = true;
        return true;
    }
    return false;
}
