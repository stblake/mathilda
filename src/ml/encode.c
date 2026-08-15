/* encode.c -- the categorical label vocabulary. See encode.h. */
#include <stdlib.h>
#include "expr.h"
#include "sym_names.h"
#include "encode.h"

#ifndef SIZE_MAX
#include <stdint.h>
#endif

void ml_labels_free(MlLabels* l) {
    if (!l) return;
    for (size_t i = 0; i < l->count; i++) expr_free(l->label[i]);
    free(l->label);
    l->label = NULL; l->count = 0;
}

size_t ml_labels_index(const MlLabels* l, const Expr* v) {
    for (size_t i = 0; i < l->count; i++)
        if (expr_eq(l->label[i], v)) return i;
    return (size_t)-1;
}

bool ml_labels_build(Expr** vals, size_t n, MlLabels* out, size_t* index) {
    out->count = 0;
    out->label = malloc(sizeof(Expr*) * (n ? n : 1));   /* at most n distinct */
    if (!out->label) return false;
    for (size_t i = 0; i < n; i++) {
        size_t at = ml_labels_index(out, vals[i]);
        if (at == (size_t)-1) {
            Expr* c = expr_copy(vals[i]);
            if (!c) { ml_labels_free(out); return false; }
            at = out->count;
            out->label[out->count++] = c;
        }
        if (index) index[i] = at;
    }
    return out->count > 0;
}

Expr* ml_labels_to_list(const MlLabels* l) {
    Expr** a = malloc(sizeof(Expr*) * (l->count ? l->count : 1));
    if (!a) return NULL;
    for (size_t i = 0; i < l->count; i++) a[i] = expr_copy(l->label[i]);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), a, l->count);
    free(a);
    return out;
}

bool ml_labels_from_list(Expr* list, MlLabels* out) {
    out->count = 0; out->label = NULL;
    if (!list || list->type != EXPR_FUNCTION) return false;
    if (list->data.function.head->type != EXPR_SYMBOL
        || list->data.function.head->data.symbol.name != SYM_List) return false;
    size_t n = list->data.function.arg_count;
    if (n == 0) return false;
    out->label = malloc(sizeof(Expr*) * n);
    if (!out->label) return false;
    for (size_t i = 0; i < n; i++) {
        out->label[i] = expr_copy(list->data.function.args[i]);
        if (!out->label[i]) { out->count = i; ml_labels_free(out); return false; }
    }
    out->count = n;
    return true;
}
