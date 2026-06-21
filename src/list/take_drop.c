#include "list_common.h"
#include "take_drop.h"

static bool get_seq_spec_indices(Expr* spec, int64_t len, int64_t** out_indices, size_t* out_count) {
    int64_t m = 0, n = 0, s = 1;
    bool is_all = false;
    bool is_none = false;

    if (spec->type == EXPR_SYMBOL) {
        if (spec->data.symbol == SYM_All) {
            is_all = true;
        } else if (spec->data.symbol == SYM_None) {
            is_none = true;
        } else {
            return false;
        }
    } else if (spec->type == EXPR_INTEGER) {
        int64_t k = spec->data.integer;
        if (k >= 0) {
            m = 1;
            n = k;
            if (n > len) return false;
        } else {
            m = len + k + 1;
            n = len;
            if (m < 1) return false;
        }
    } else if (spec->type == EXPR_FUNCTION) {
        const char* head = spec->data.function.head->type == EXPR_SYMBOL ? spec->data.function.head->data.symbol : "";
        if (strcmp(head, "UpTo") == 0 && spec->data.function.arg_count == 1 && spec->data.function.args[0]->type == EXPR_INTEGER) {
            int64_t k = spec->data.function.args[0]->data.integer;
            if (k >= 0) {
                m = 1;
                n = k > len ? len : k;
            } else {
                m = len + k + 1;
                if (m < 1) m = 1;
                n = len;
            }
        } else if (strcmp(head, "List") == 0) {
            size_t count = spec->data.function.arg_count;
            if (count >= 1 && count <= 3) {
                if (spec->data.function.args[0]->type != EXPR_INTEGER) return false;
                m = spec->data.function.args[0]->data.integer;
                m = m < 0 ? len + m + 1 : m;
                
                if (count == 1) {
                    n = m;
                } else {
                    if (spec->data.function.args[1]->type != EXPR_INTEGER) return false;
                    n = spec->data.function.args[1]->data.integer;
                    n = n < 0 ? len + n + 1 : n;
                    
                    if (count == 3) {
                        if (spec->data.function.args[2]->type != EXPR_INTEGER) return false;
                        s = spec->data.function.args[2]->data.integer;
                        if (s == 0) return false;
                    }
                }
            } else {
                return false;
            }
        } else {
            return false;
        }
    } else {
        return false;
    }

    if (is_all) {
        *out_count = (size_t)len;
        if (len > 0) {
            *out_indices = malloc(sizeof(int64_t) * (size_t)len);
            for (int64_t i = 0; i < len; i++) (*out_indices)[i] = i + 1;
        } else {
            *out_indices = NULL;
        }
        return true;
    }
    if (is_none) {
        *out_count = 0;
        *out_indices = NULL;
        return true;
    }

    int64_t count = 0;
    if (s > 0) {
        if (m <= n) count = (n - m) / s + 1;
    } else if (s < 0) {
        if (m >= n) count = (m - n) / (-s) + 1;
    }

    if (count > 0) {
        int64_t* indices = malloc(sizeof(int64_t) * (size_t)count);
        int64_t idx = m;
        for (int64_t i = 0; i < count; i++) {
            if (idx < 1 || idx > len) {
                free(indices);
                return false;
            }
            indices[i] = idx;
            idx += s;
        }
        *out_indices = indices;
        *out_count = (size_t)count;
    } else {
        *out_indices = NULL;
        *out_count = 0;
    }
    return true;
}

static Expr* apply_take_drop(Expr* expr, Expr** specs, size_t nspecs, bool is_take) {
    if (nspecs == 0) return expr_copy(expr);
    if (expr->type != EXPR_FUNCTION) return NULL;

    int64_t len = (int64_t)expr->data.function.arg_count;
    size_t spec_count = 0;
    int64_t* spec_indices = NULL;
    if (!get_seq_spec_indices(specs[0], len, &spec_indices, &spec_count)) return NULL;

    Expr** new_args = NULL;
    size_t new_count = 0;

    if (is_take) {
        new_count = spec_count;
        if (new_count > 0) new_args = malloc(sizeof(Expr*) * new_count);
        for (size_t i = 0; i < new_count; i++) {
            Expr* sub = expr->data.function.args[spec_indices[i] - 1];
            new_args[i] = apply_take_drop(sub, specs + 1, nspecs - 1, is_take);
            if (!new_args[i]) {
                for (size_t j = 0; j < i; j++) expr_free(new_args[j]);
                free(new_args);
                free(spec_indices);
                return NULL;
            }
        }
    } else {
        bool* keep = malloc(sizeof(bool) * (size_t)len);
        for (int64_t i = 0; i < len; i++) keep[i] = true;
        for (size_t i = 0; i < spec_count; i++) keep[spec_indices[i] - 1] = false;

        for (int64_t i = 0; i < len; i++) if (keep[i]) new_count++;
        if (new_count > 0) new_args = malloc(sizeof(Expr*) * new_count);
        
        size_t idx = 0;
        for (int64_t i = 0; i < len; i++) {
            if (keep[i]) {
                Expr* sub = expr->data.function.args[i];
                new_args[idx] = apply_take_drop(sub, specs + 1, nspecs - 1, is_take);
                if (!new_args[idx]) {
                    for (size_t j = 0; j < idx; j++) expr_free(new_args[j]);
                    free(new_args);
                    free(keep);
                    free(spec_indices);
                    return NULL;
                }
                idx++;
            }
        }
        free(keep);
    }

    free(spec_indices);
    Expr* result = expr_new_function(expr_copy(expr->data.function.head), new_args, new_count);
    if (new_args) free(new_args);
    return result;
}

Expr* builtin_take(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    return apply_take_drop(res->data.function.args[0], res->data.function.args + 1, res->data.function.arg_count - 1, true);
}

Expr* builtin_drop(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    return apply_take_drop(res->data.function.args[0], res->data.function.args + 1, res->data.function.arg_count - 1, false);
}
