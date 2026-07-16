#include "list_common.h"
#include "transpose.h"

static int get_array_dimensions(Expr* e, int64_t* dims, const char* head_name) {
    if (!head_is(e, intern_symbol(head_name))) {
        return 0;
    }
    dims[0] = (int64_t)e->data.function.arg_count;
    if (dims[0] == 0) return 1;
    int64_t sub_dims[64];
    int depth = get_array_dimensions(e->data.function.args[0], sub_dims, head_name);
    for (size_t i = 1; i < e->data.function.arg_count; i++) {
        int64_t cur_dims[64];
        if (get_array_dimensions(e->data.function.args[i], cur_dims, head_name) != depth) return 1;
        for (int j = 0; j < depth; j++) if (cur_dims[j] != sub_dims[j]) return 1;
    }
    for (int i = 0; i < depth; i++) dims[i + 1] = sub_dims[i];
    return depth + 1;
}

static Expr* get_element_at(Expr* e, int64_t* indices, size_t depth) {
    Expr* curr = e;
    for (size_t i = 0; i < depth; i++) {
        curr = curr->data.function.args[indices[i]];
    }
    return curr;
}

static Expr* build_transposed(const char* head, int64_t* out_dims, size_t out_depth, int64_t* out_indices_base, int64_t* current_out_indices, 
                             int64_t* in_indices, int64_t* perm, size_t in_depth, Expr* original) {
    if (out_depth == 0) {
        // level k in list is Subscript[n, k]-th level in result.
        // So in_indices[k] = out_indices_base[perm[k] - 1]
        for (size_t k = 0; k < in_depth; k++) {
            in_indices[k] = out_indices_base[perm[k] - 1];
        }
        return expr_copy(get_element_at(original, in_indices, in_depth));
    }

    size_t len = (size_t)out_dims[0];
    Expr** args = malloc(sizeof(Expr*) * len);
    for (size_t i = 0; i < len; i++) {
        current_out_indices[0] = (int64_t)i;
        args[i] = build_transposed(head, out_dims + 1, out_depth - 1, out_indices_base, current_out_indices + 1, in_indices, perm, in_depth, original);
    }
    Expr* res = expr_new_function(expr_new_symbol(head), args, len);
    free(args);
    return res;
}

Expr* builtin_transpose(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* list = res->data.function.args[0];
    if (list->type != EXPR_FUNCTION || list->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* head = list->data.function.head->data.symbol.name;

    int64_t in_dims[64];
    int in_depth = get_array_dimensions(list, in_dims, head);
    if (in_depth < 2) return NULL;

    int64_t* perm = malloc(sizeof(int64_t) * in_depth);
    if (res->data.function.arg_count == 1) {
        perm[0] = 2; perm[1] = 1;
        for (int i = 2; i < in_depth; i++) perm[i] = i + 1;
    } else {
        Expr* spec = res->data.function.args[1];
        if (spec->type != EXPR_FUNCTION || spec->data.function.head->data.symbol.name != SYM_List || 
            spec->data.function.arg_count != (size_t)in_depth) {
            free(perm); return NULL;
        }
        for (int i = 0; i < in_depth; i++) {
            if (spec->data.function.args[i]->type != EXPR_INTEGER) { free(perm); return NULL; }
            perm[i] = spec->data.function.args[i]->data.integer;
        }
    }

    int out_depth = 0;
    for (int i = 0; i < in_depth; i++) if (perm[i] > out_depth) out_depth = (int)perm[i];
    
    int64_t* out_dims = malloc(sizeof(int64_t) * out_depth);
    for (int i = 0; i < out_depth; i++) out_dims[i] = -1;

    for (int i = 0; i < in_depth; i++) {
        int target_idx = (int)perm[i] - 1;
        if (target_idx < 0) { free(perm); free(out_dims); return NULL; }
        if (out_dims[target_idx] == -1 || in_dims[i] < out_dims[target_idx]) {
            out_dims[target_idx] = in_dims[i];
        }
    }

    int64_t* out_indices_base = calloc(out_depth, sizeof(int64_t));
    int64_t* in_indices = malloc(sizeof(int64_t) * in_depth);

    Expr* result = build_transposed(head, out_dims, (size_t)out_depth, out_indices_base, out_indices_base, in_indices, perm, (size_t)in_depth, list);

    free(perm); free(out_dims); free(out_indices_base); free(in_indices);
    return result;
}

/* ConjugateTranspose[m] is equivalent to Conjugate[Transpose[m]];
 * ConjugateTranspose[m, spec] is equivalent to Conjugate[Transpose[m, spec]].
 * For a 1-D vector input (1-arg form) the result keeps the vector shape and
 * the entries are conjugated, matching Mathematica's behaviour. Conjugate is
 * Listable, so wrapping the (possibly nested) transposed array in Conjugate
 * lets the evaluator thread conjugation down to every leaf at fixed point. */
Expr* builtin_conjugate_transpose(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    Expr* m = res->data.function.args[0];

    int64_t dims[64];
    int depth = 0;
    if (m->type == EXPR_FUNCTION && m->data.function.head->type == EXPR_SYMBOL &&
        m->data.function.head->data.symbol.name == SYM_List) {
        depth = get_array_dimensions(m, dims, "List");
    }

    /* Not a (rectangular) List: leave the call unevaluated so that symbolic
     * matrices like ConjugateTranspose[A] stay intact. */
    if (depth == 0) return NULL;

    /* 1-arg form on a 1-D vector: conjugate elementwise, keep the shape. */
    if (argc == 1 && depth == 1) {
        Expr** conj_args = malloc(sizeof(Expr*) * 1);
        conj_args[0] = expr_copy(m);
        Expr* conj = expr_new_function(expr_new_symbol(SYM_Conjugate), conj_args, 1);
        free(conj_args);
        return eval_and_free(conj);
    }

    /* 1-arg form on lower depth (depth < 2) we cannot transpose: leave it. */
    if (argc == 1 && depth < 2) return NULL;

    /* Build Transpose[m] or Transpose[m, spec]. */
    Expr** tr_args = malloc(sizeof(Expr*) * argc);
    tr_args[0] = expr_copy(m);
    if (argc == 2) tr_args[1] = expr_copy(res->data.function.args[1]);
    Expr* transposed_call = expr_new_function(expr_new_symbol(SYM_Transpose), tr_args, argc);
    free(tr_args);

    Expr* transposed = eval_and_free(transposed_call);

    /* If Transpose could not reduce (e.g. invalid spec), surface the
     * unevaluated ConjugateTranspose rather than a spurious
     * Conjugate[Transpose[...]] wrapper. */
    if (head_is(transposed, SYM_Transpose)) {
        expr_free(transposed);
        return NULL;
    }

    Expr** conj_args = malloc(sizeof(Expr*) * 1);
    conj_args[0] = transposed;
    Expr* conj = expr_new_function(expr_new_symbol(SYM_Conjugate), conj_args, 1);
    free(conj_args);
    return eval_and_free(conj);
}
