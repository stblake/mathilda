#include "list.h"
#include <inttypes.h>
#include "common.h"
#include "symtab.h"
#include "eval.h"
#include "iter.h"
#include "core.h"
#include "arithmetic.h"
#include "print.h"
#include "sym_intern.h"
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

Expr* builtin_table(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    
    if (res->data.function.arg_count > 2) {
        Expr** inner_args = malloc(sizeof(Expr*) * 2);
        inner_args[0] = expr_copy(res->data.function.args[0]);
        inner_args[1] = expr_copy(res->data.function.args[res->data.function.arg_count - 1]);
        Expr* inner_table = expr_new_function(expr_new_symbol("Table"), inner_args, 2);
        free(inner_args);
        
        Expr** outer_args = malloc(sizeof(Expr*) * (res->data.function.arg_count - 1));
        outer_args[0] = inner_table;
        for (size_t i = 1; i < res->data.function.arg_count - 1; i++) {
            outer_args[i] = expr_copy(res->data.function.args[i]);
        }
        Expr* outer_table = expr_new_function(expr_new_symbol("Table"), outer_args, res->data.function.arg_count - 1);
        free(outer_args);
        
        Expr* eval_outer = evaluate(outer_table);
        expr_free(outer_table);
        return eval_outer;
    }
    
    Expr* expr = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    /* ---- Parse the iterator spec (shared helper) ---- */
    IterSpec s;
    if (!iter_spec_parse(spec, &s)) return NULL;

    int is_n_times   = (s.kind == ITER_KIND_COUNT);
    int is_list_iter = (s.kind == ITER_KIND_LIST);
    double min_val = 0, max_val = 0, di_val = 0;
    bool is_real = false, is_inf = false;

    /* Table does not iterate to Infinity; allow_inf = false. */
    if (!is_list_iter) {
        if (!iter_spec_resolve_numeric(&s, /*allow_inf=*/false,
                                       &min_val, &max_val, &di_val,
                                       &is_real, &is_inf)) {
            iter_spec_free(&s);
            return NULL;
        }
    }

    /* Convenience aliases into the owned IterSpec (freed via iter_spec_free). */
    Expr* var_sym = s.var;
    Expr* imin_e  = s.imin;
    Expr* imax_e  = s.imax;
    Expr* di_e    = s.di;
    Expr* list_e  = s.list;

    size_t results_cap = 16;
    size_t results_count = 0;
    Expr** results = malloc(sizeof(Expr*) * results_cap);

    Rule* old_own = iter_spec_shadow(var_sym);

    if (is_n_times) {
        int64_t n = imax_e->data.integer;
        for (int64_t i = 0; i < n; i++) {
            Expr* eval_expr = evaluate(expr);
            if (results_count == results_cap) { results_cap *= 2; results = realloc(results, sizeof(Expr*) * results_cap); }
            results[results_count++] = eval_expr;
        }
    } else if (is_list_iter) {
        for (size_t i = 0; i < list_e->data.function.arg_count; i++) {
            symtab_add_own_value(var_sym->data.symbol, var_sym, list_e->data.function.args[i]);
            Expr* eval_expr = evaluate(expr);
            if (results_count == results_cap) { results_cap *= 2; results = realloc(results, sizeof(Expr*) * results_cap); }
            results[results_count++] = eval_expr;
        }
    } else {
        double val = min_val;
        int steps = 0;
        Expr* curr_e = expr_copy(imin_e);
        while ((di_val > 0 && val <= max_val + 1e-14) || (di_val < 0 && val >= max_val - 1e-14)) {
            Expr* i_val = is_real ? expr_new_real(val) : expr_copy(curr_e);
            symtab_add_own_value(var_sym->data.symbol, var_sym, i_val);
            
            Expr* eval_expr = evaluate(expr);
            if (results_count == results_cap) { results_cap *= 2; results = realloc(results, sizeof(Expr*) * results_cap); }
            results[results_count++] = eval_expr;
            
            expr_free(i_val);
            
            Expr* next_args[2] = { expr_copy(curr_e), expr_copy(di_e) };
            Expr* next_expr = expr_new_function(expr_new_symbol("Plus"), next_args, 2);
            Expr* next_e = evaluate(next_expr);
            expr_free(next_expr);
            expr_free(curr_e);
            curr_e = next_e;
            
            val += di_val;
            steps++;
            if (steps > 1000000) break; 
        }
        if (curr_e) expr_free(curr_e);
    }

    iter_spec_restore(var_sym, old_own);
    iter_spec_free(&s);

    Expr* result_list = expr_new_function(expr_new_symbol("List"), results, results_count);
    free(results);
    return result_list;
}

static Expr* array_helper(Expr* f, Expr** n_array, Expr** r_array, size_t dim_count, size_t current_dim, Expr** current_args) {
    if (current_dim == dim_count) {
        Expr** fn_args = malloc(sizeof(Expr*) * dim_count);
        for (size_t i = 0; i < dim_count; i++) fn_args[i] = expr_copy(current_args[i]);
        Expr* fn_expr = expr_new_function(expr_copy(f), fn_args, dim_count);
        free(fn_args);
        Expr* eval_fn = evaluate(fn_expr);
        expr_free(fn_expr);
        return eval_fn;
    }

    Expr* n_expr = n_array[current_dim];
    Expr* r_expr = r_array[current_dim];
    
    int64_t n_val = n_expr->data.integer;
    Expr** results = malloc(sizeof(Expr*) * (size_t)n_val);
    if (n_val == 0) {
        free(results);
        return expr_new_function(expr_new_symbol("List"), NULL, 0);
    }
    
    int is_range = 0;
    Expr* a_expr = NULL;
    Expr* b_expr = NULL;
    
    if (r_expr && r_expr->type == EXPR_FUNCTION && r_expr->data.function.head->type == EXPR_SYMBOL && r_expr->data.function.head->data.symbol == SYM_List && r_expr->data.function.arg_count == 2) {
        is_range = 1;
        a_expr = r_expr->data.function.args[0];
        b_expr = r_expr->data.function.args[1];
    }
    
    Expr* r_base = NULL;
    if (!is_range) {
        r_base = r_expr ? expr_copy(r_expr) : expr_new_integer(1);
    }
    
    for (int64_t i = 0; i < n_val; i++) {
        Expr* arg = NULL;
        if (is_range) {
            if (n_val == 1) {
                arg = expr_copy(a_expr);
            } else {
                Expr* diff = expr_new_function(expr_new_symbol("Plus"), (Expr*[]){expr_copy(b_expr), expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), expr_copy(a_expr)}, 2)}, 2);
                Expr* frac = expr_new_function(expr_new_symbol("Divide"), (Expr*[]){expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(i), diff}, 2), expr_new_integer(n_val - 1)}, 2);
                Expr* sum = expr_new_function(expr_new_symbol("Plus"), (Expr*[]){expr_copy(a_expr), frac}, 2);
                arg = evaluate(sum);
                expr_free(sum);
            }
        } else {
            if (i == 0) {
                arg = expr_copy(r_base);
            } else {
                Expr* sum = expr_new_function(expr_new_symbol("Plus"), (Expr*[]){expr_copy(r_base), expr_new_integer(i)}, 2);
                arg = evaluate(sum);
                expr_free(sum);
            }
        }
        
        current_args[current_dim] = arg;
        results[i] = array_helper(f, n_array, r_array, dim_count, current_dim + 1, current_args);
        expr_free(arg);
    }
    
    if (r_base) expr_free(r_base);
    
    Expr* list_result = expr_new_function(expr_new_symbol("List"), results, (size_t)n_val);
    free(results);
    return list_result;
}

Expr* builtin_array(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2 || res->data.function.arg_count > 3) return NULL;
    
    Expr* f = res->data.function.args[0];
    Expr* n_spec = res->data.function.args[1];
    Expr* r_spec = res->data.function.arg_count == 3 ? res->data.function.args[2] : NULL;
    
    size_t dim_count = 1;
    if (n_spec->type == EXPR_FUNCTION && n_spec->data.function.head->type == EXPR_SYMBOL && n_spec->data.function.head->data.symbol == SYM_List) {
        dim_count = n_spec->data.function.arg_count;
        if (dim_count == 0) return NULL;
    }
    
    Expr** n_array = malloc(sizeof(Expr*) * dim_count);
    Expr** r_array = malloc(sizeof(Expr*) * dim_count);
    
    if (n_spec->type == EXPR_FUNCTION && n_spec->data.function.head->type == EXPR_SYMBOL && n_spec->data.function.head->data.symbol == SYM_List) {
        for(size_t i=0; i<dim_count; i++) n_array[i] = n_spec->data.function.args[i];
    } else {
        n_array[0] = n_spec;
    }
    
    if (r_spec) {
        if (r_spec->type == EXPR_FUNCTION && r_spec->data.function.head->type == EXPR_SYMBOL && r_spec->data.function.head->data.symbol == SYM_List && r_spec->data.function.arg_count == dim_count) {
            for(size_t i=0; i<dim_count; i++) r_array[i] = r_spec->data.function.args[i];
        } else {
            for(size_t i=0; i<dim_count; i++) r_array[i] = r_spec;
        }
    } else {
        for(size_t i=0; i<dim_count; i++) r_array[i] = NULL;
    }
    
    for (size_t i = 0; i < dim_count; i++) {
        if (n_array[i]->type != EXPR_INTEGER || n_array[i]->data.integer < 0) {
            free(n_array);
            free(r_array);
            return NULL;
        }
    }
    
    Expr** current_args = malloc(sizeof(Expr*) * dim_count);
    Expr* result = array_helper(f, n_array, r_array, dim_count, 0, current_args);
    free(current_args);
    
    free(n_array);
    free(r_array);
    
    return result;
}

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

static void flatten_rec(Expr* e, const char* h, int64_t level, Expr*** results, size_t* count, size_t* cap) {
    if (level != 0 && head_is(e, intern_symbol(h))) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            flatten_rec(e->data.function.args[i], h, level == -1 ? -1 : level - 1, results, count, cap);
        }
    } else {
        if (*count == *cap) {
            *cap *= 2;
            *results = realloc(*results, sizeof(Expr*) * (*cap));
        }
        (*results)[(*count)++] = expr_copy(e);
    }
}

Expr* builtin_flatten(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 3) return NULL;

    Expr* list = res->data.function.args[0];
    if (list->type != EXPR_FUNCTION) return expr_copy(list);

    int64_t n = -1; // -1 means infinity
    if (res->data.function.arg_count >= 2) {
        if (res->data.function.args[1]->type != EXPR_INTEGER) return NULL;
        n = res->data.function.args[1]->data.integer;
    }

    const char* h = SYM_List;
    if (res->data.function.arg_count == 3) {
        if (res->data.function.args[2]->type != EXPR_SYMBOL) return NULL;
        h = res->data.function.args[2]->data.symbol;
    }

    size_t cap = 16;
    size_t count = 0;
    Expr** results = malloc(sizeof(Expr*) * cap);

    // Initial call: we flatten children of the head if they also have head h.
    for (size_t i = 0; i < list->data.function.arg_count; i++) {
        flatten_rec(list->data.function.args[i], h, n, &results, &count, &cap);
    }

    Expr* result = expr_new_function(expr_copy(list->data.function.head), results, count);
    free(results);
    return result;
}

static Expr* partition_rec(Expr* list, Expr* n_spec, Expr* d_spec, size_t level_idx) {
    if (list->type != EXPR_FUNCTION) return expr_copy(list);

    // Get n and d for this level
    int64_t n = -1;
    bool n_upto = false;
    Expr* n_e = (n_spec->type == EXPR_FUNCTION && n_spec->data.function.head->data.symbol == SYM_List) ? 
                (level_idx < n_spec->data.function.arg_count ? n_spec->data.function.args[level_idx] : NULL) : 
                (level_idx == 0 ? n_spec : NULL);
    
    if (!n_e) return expr_copy(list);

    if (n_e->type == EXPR_INTEGER) {
        n = n_e->data.integer;
    } else if (n_e->type == EXPR_FUNCTION && n_e->data.function.head->data.symbol == SYM_UpTo && n_e->data.function.arg_count == 1) {
        if (n_e->data.function.args[0]->type == EXPR_INTEGER) {
            n = n_e->data.function.args[0]->data.integer;
            n_upto = true;
        }
    }
    if (n <= 0) return expr_copy(list);

    int64_t d = n;
    if (d_spec) {
        Expr* d_e = (d_spec->type == EXPR_FUNCTION && d_spec->data.function.head->data.symbol == SYM_List) ? 
                    (level_idx < d_spec->data.function.arg_count ? d_spec->data.function.args[level_idx] : NULL) : 
                    (level_idx == 0 ? d_spec : NULL);
        if (d_e && d_e->type == EXPR_INTEGER) {
            d = d_e->data.integer;
        }
    }
    if (d <= 0) return expr_copy(list);

    size_t len = list->data.function.arg_count;
    size_t num_sublists = 0;
    if (n_upto) {
        num_sublists = (len > 0) ? (len - 1) / d + 1 : 0;
    } else {
        if (len >= (size_t)n) {
            num_sublists = (len - (size_t)n) / (size_t)d + 1;
        }
    }

    Expr** sublists = malloc(sizeof(Expr*) * num_sublists);
    for (size_t i = 0; i < num_sublists; i++) {
        size_t start = i * (size_t)d;
        size_t end = start + (size_t)n;
        if (end > len) end = len;
        
        size_t sub_count = end - start;
        Expr** sub_args = malloc(sizeof(Expr*) * sub_count);
        for (size_t j = 0; j < sub_count; j++) {
            sub_args[j] = partition_rec(list->data.function.args[start + j], n_spec, d_spec, level_idx + 1);
        }
        sublists[i] = expr_new_function(expr_copy(list->data.function.head), sub_args, sub_count);
        free(sub_args);
    }

    Expr* result = expr_new_function(expr_copy(list->data.function.head), sublists, num_sublists);
    free(sublists);
    return result;
}

Expr* builtin_partition(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2 || res->data.function.arg_count > 3) return NULL;
    
    Expr* list = res->data.function.args[0];
    Expr* n_spec = res->data.function.args[1];
    Expr* d_spec = (res->data.function.arg_count == 3) ? res->data.function.args[2] : NULL;

    if (list->type != EXPR_FUNCTION) return expr_copy(list);

    return partition_rec(list, n_spec, d_spec, 0);
}

static Expr* rotate_rec(Expr* expr, Expr* n_spec, size_t level_idx) {
    if (expr->type != EXPR_FUNCTION) return expr_copy(expr);

    int64_t n = 0;
    if (n_spec->type == EXPR_INTEGER) {
        if (level_idx == 0) n = n_spec->data.integer;
    } else if (n_spec->type == EXPR_FUNCTION && n_spec->data.function.head->data.symbol == SYM_List) {
        if (level_idx < n_spec->data.function.arg_count) {
            Expr* sub_n = n_spec->data.function.args[level_idx];
            if (sub_n->type == EXPR_INTEGER) n = sub_n->data.integer;
        }
    }

    size_t len = expr->data.function.arg_count;
    Expr** new_args = malloc(sizeof(Expr*) * len);

    if (len > 0) {
        int64_t offset = n % (int64_t)len;
        if (offset < 0) offset += (int64_t)len;

        for (size_t i = 0; i < len; i++) {
            size_t old_idx = (i + (size_t)offset) % len;
            new_args[i] = rotate_rec(expr->data.function.args[old_idx], n_spec, level_idx + 1);
        }
    }

    Expr* result = expr_new_function(expr_copy(expr->data.function.head), new_args, len);
    free(new_args);
    return result;
}

Expr* builtin_rotateleft(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* n_spec = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    
    Expr* default_n = NULL;
    if (!n_spec) {
        default_n = expr_new_integer(1);
        n_spec = default_n;
    }

    Expr* ret = rotate_rec(expr, n_spec, 0);
    if (default_n) expr_free(default_n);
    return ret;
}

Expr* builtin_rotateright(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* n_spec = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    
    Expr* neg_n_spec = NULL;
    if (!n_spec) {
        neg_n_spec = expr_new_integer(-1);
    } else if (n_spec->type == EXPR_INTEGER) {
        neg_n_spec = expr_new_integer(-n_spec->data.integer);
    } else if (n_spec->type == EXPR_FUNCTION && n_spec->data.function.head->data.symbol == SYM_List) {
        Expr** neg_args = malloc(sizeof(Expr*) * n_spec->data.function.arg_count);
        for (size_t i = 0; i < n_spec->data.function.arg_count; i++) {
            if (n_spec->data.function.args[i]->type == EXPR_INTEGER) {
                neg_args[i] = expr_new_integer(-n_spec->data.function.args[i]->data.integer);
            } else {
                neg_args[i] = expr_copy(n_spec->data.function.args[i]);
            }
        }
        neg_n_spec = expr_new_function(expr_new_symbol("List"), neg_args, n_spec->data.function.arg_count);
        free(neg_args);
    } else {
        return NULL;
    }

    Expr* ret = rotate_rec(expr, neg_n_spec, 0);
    expr_free(neg_n_spec);
    return ret;
}

static bool should_reverse_at_level(Expr* level_spec, size_t current_level) {
    if (!level_spec) return current_level == 1;
    if (level_spec->type == EXPR_INTEGER) return (size_t)level_spec->data.integer == current_level;
    if (level_spec->type == EXPR_FUNCTION && level_spec->data.function.head->data.symbol == SYM_List) {
        for (size_t i = 0; i < level_spec->data.function.arg_count; i++) {
            if (level_spec->data.function.args[i]->type == EXPR_INTEGER && 
                (size_t)level_spec->data.function.args[i]->data.integer == current_level) return true;
        }
    }
    return false;
}

static Expr* reverse_rec(Expr* expr, Expr* level_spec, size_t current_level) {
    if (expr->type != EXPR_FUNCTION) return expr_copy(expr);

    size_t len = expr->data.function.arg_count;
    Expr** new_args = malloc(sizeof(Expr*) * len);
    bool do_rev = should_reverse_at_level(level_spec, current_level);

    for (size_t i = 0; i < len; i++) {
        size_t src_idx = do_rev ? (len - 1 - i) : i;
        new_args[i] = reverse_rec(expr->data.function.args[src_idx], level_spec, current_level + 1);
    }

    Expr* result = expr_new_function(expr_copy(expr->data.function.head), new_args, len);
    free(new_args);
    return result;
}

Expr* builtin_reverse(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* level_spec = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;

    return reverse_rec(expr, level_spec, 1);
}

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
    const char* head = list->data.function.head->data.symbol;

    int64_t in_dims[64];
    int in_depth = get_array_dimensions(list, in_dims, head);
    if (in_depth < 2) return NULL;

    int64_t* perm = malloc(sizeof(int64_t) * in_depth);
    if (res->data.function.arg_count == 1) {
        perm[0] = 2; perm[1] = 1;
        for (int i = 2; i < in_depth; i++) perm[i] = i + 1;
    } else {
        Expr* spec = res->data.function.args[1];
        if (spec->type != EXPR_FUNCTION || spec->data.function.head->data.symbol != SYM_List || 
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
        m->data.function.head->data.symbol == SYM_List) {
        depth = get_array_dimensions(m, dims, "List");
    }

    /* Not a (rectangular) List: leave the call unevaluated so that symbolic
     * matrices like ConjugateTranspose[A] stay intact. */
    if (depth == 0) return NULL;

    /* 1-arg form on a 1-D vector: conjugate elementwise, keep the shape. */
    if (argc == 1 && depth == 1) {
        Expr** conj_args = malloc(sizeof(Expr*) * 1);
        conj_args[0] = expr_copy(m);
        Expr* conj = expr_new_function(expr_new_symbol("Conjugate"), conj_args, 1);
        free(conj_args);
        return eval_and_free(conj);
    }

    /* 1-arg form on lower depth (depth < 2) we cannot transpose: leave it. */
    if (argc == 1 && depth < 2) return NULL;

    /* Build Transpose[m] or Transpose[m, spec]. */
    Expr** tr_args = malloc(sizeof(Expr*) * argc);
    tr_args[0] = expr_copy(m);
    if (argc == 2) tr_args[1] = expr_copy(res->data.function.args[1]);
    Expr* transposed_call = expr_new_function(expr_new_symbol("Transpose"), tr_args, argc);
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
    Expr* conj = expr_new_function(expr_new_symbol("Conjugate"), conj_args, 1);
    free(conj_args);
    return eval_and_free(conj);
}

static int compare_expr_ptrs(const void* a, const void* b) {
    Expr* ea = *(Expr**)a;
    Expr* eb = *(Expr**)b;
    return expr_compare(ea, eb);
}

Expr* builtin_union(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;
    
    // Find options
    Expr* same_test = NULL;
    size_t last_arg = res->data.function.arg_count;
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        Expr* arg = res->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol == SYM_Rule &&
            arg->data.function.arg_count == 2 &&
            arg->data.function.args[0]->type == EXPR_SYMBOL &&
            arg->data.function.args[0]->data.symbol == SYM_SameTest) {
            same_test = arg->data.function.args[1];
            if (i < last_arg) last_arg = i;
        }
    }
    
    if (last_arg == 0) return NULL;
    
    // Check if first arg is a function
    Expr* first_list = res->data.function.args[0];
    if (first_list->type != EXPR_FUNCTION) return expr_copy(first_list);
    
    Expr* common_head = first_list->data.function.head;
    
    // Total count of elements
    size_t total_count = 0;
    for (size_t i = 0; i < last_arg; i++) {
        Expr* arg = res->data.function.args[i];
        if (arg->type != EXPR_FUNCTION || !expr_eq(arg->data.function.head, common_head)) {
            // Heads must match
            return NULL;
        }
        total_count += arg->data.function.arg_count;
    }
    
    if (total_count == 0) return expr_copy(first_list);
    
    Expr** all_args = malloc(sizeof(Expr*) * total_count);
    size_t idx = 0;
    for (size_t i = 0; i < last_arg; i++) {
        Expr* arg = res->data.function.args[i];
        for (size_t j = 0; j < arg->data.function.arg_count; j++) {
            all_args[idx++] = expr_copy(arg->data.function.args[j]);
        }
    }
    
    // Sort elements
    qsort(all_args, total_count, sizeof(Expr*), compare_expr_ptrs);
    
    // Remove duplicates
    Expr** unique_args = malloc(sizeof(Expr*) * total_count);
    size_t unique_count = 0;
    
    if (total_count > 0) {
        unique_args[unique_count++] = all_args[0];
        for (size_t i = 1; i < total_count; i++) {
            bool is_dup = false;
            if (same_test == NULL) {
                if (expr_eq(all_args[i], unique_args[unique_count - 1])) {
                    is_dup = true;
                }
            } else {
                Expr* call_args[2] = { expr_copy(all_args[i]), expr_copy(unique_args[unique_count - 1]) };
                Expr* call = expr_new_function(expr_copy(same_test), call_args, 2);
                Expr* eval_res = evaluate(call);
                if (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol == SYM_True) {
                    is_dup = true;
                }
                expr_free(eval_res);
                expr_free(call);
            }
            
            if (is_dup) {
                expr_free(all_args[i]);
            } else {
                unique_args[unique_count++] = all_args[i];
            }
        }
    }
    
    free(all_args);
    
    Expr* result = expr_new_function(expr_copy(common_head), unique_args, unique_count);
    if (unique_args) free(unique_args);
    
    return result;
}

typedef struct HashNode {
    Expr* key;
    size_t index; // Original index or position in unique_elems
    struct HashNode* next;
} HashNode;

typedef struct {
    HashNode** buckets;
    size_t size;
} HashTable;

static HashTable* ht_create(size_t size) {
    HashTable* ht = malloc(sizeof(HashTable));
    ht->size = size;
    ht->buckets = calloc(size, sizeof(HashNode*));
    return ht;
}

static void ht_free(HashTable* ht, bool free_keys) {
    for (size_t i = 0; i < ht->size; i++) {
        HashNode* node = ht->buckets[i];
        while (node) {
            HashNode* next = node->next;
            if (free_keys) expr_free(node->key);
            free(node);
            node = next;
        }
    }
    free(ht->buckets);
    free(ht);
}

static HashNode* ht_find(HashTable* ht, Expr* key) {
    uint64_t h = expr_hash(key);
    size_t bucket = (size_t)(h % ht->size);
    HashNode* node = ht->buckets[bucket];
    while (node) {
        if (expr_eq(node->key, key)) return node;
        node = node->next;
    }
    return NULL;
}

static void ht_insert(HashTable* ht, Expr* key, size_t index) {
    uint64_t h = expr_hash(key);
    size_t bucket = (size_t)(h % ht->size);
    HashNode* node = malloc(sizeof(HashNode));
    node->key = key;
    node->index = index;
    node->next = ht->buckets[bucket];
    ht->buckets[bucket] = node;
}

Expr* builtin_tally(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    
    Expr* list = res->data.function.args[0];
    Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    
    if (list->type != EXPR_FUNCTION) return expr_new_function(expr_new_symbol("List"), NULL, 0);
    
    size_t count = list->data.function.arg_count;
    if (count == 0) return expr_new_function(expr_new_symbol("List"), NULL, 0);
    
    Expr** unique_elems = malloc(sizeof(Expr*) * count);
    int64_t* multiplicities = malloc(sizeof(int64_t) * count);
    size_t unique_count = 0;

    if (test == NULL) {
        HashTable* ht = ht_create(count * 2 + 1);
        for (size_t i = 0; i < count; i++) {
            Expr* elem = list->data.function.args[i];
            HashNode* node = ht_find(ht, elem);
            if (node) {
                multiplicities[node->index]++;
            } else {
                unique_elems[unique_count] = expr_copy(elem);
                multiplicities[unique_count] = 1;
                ht_insert(ht, unique_elems[unique_count], unique_count);
                unique_count++;
            }
        }
        ht_free(ht, false);
    } else {
        // Fallback to O(N^2) for custom test
        for (size_t i = 0; i < count; i++) {
            Expr* elem = list->data.function.args[i];
            int found_idx = -1;
            for (size_t j = 0; j < unique_count; j++) {
                Expr* call_args[2] = { expr_copy(elem), expr_copy(unique_elems[j]) };
                Expr* call = expr_new_function(expr_copy(test), call_args, 2);
                Expr* eval_res = evaluate(call);
                if (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol == SYM_True) {
                    found_idx = (int)j;
                    expr_free(eval_res);
                    expr_free(call);
                    break;
                }
                expr_free(eval_res);
                expr_free(call);
            }
            if (found_idx != -1) {
                multiplicities[found_idx]++;
            } else {
                unique_elems[unique_count] = expr_copy(elem);
                multiplicities[unique_count] = 1;
                unique_count++;
            }
        }
    }
    
    Expr** result_args = malloc(sizeof(Expr*) * unique_count);
    for (size_t i = 0; i < unique_count; i++) {
        Expr** pair_args = malloc(sizeof(Expr*) * 2);
        pair_args[0] = unique_elems[i];
        pair_args[1] = expr_new_integer(multiplicities[i]);
        result_args[i] = expr_new_function(expr_new_symbol("List"), pair_args, 2);
        free(pair_args);
    }
    
    free(unique_elems);
    free(multiplicities);
    
    Expr* result = expr_new_function(expr_new_symbol("List"), result_args, unique_count);
    free(result_args);
    
    return result;
}

Expr* builtin_deleteduplicates(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    
    Expr* list = res->data.function.args[0];
    Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    
    if (list->type != EXPR_FUNCTION) return expr_copy(list);
    
    size_t count = list->data.function.arg_count;
    if (count == 0) return expr_copy(list);
    
    Expr** unique_args = malloc(sizeof(Expr*) * count);
    size_t unique_count = 0;

    if (test == NULL) {
        HashTable* ht = ht_create(count * 2 + 1);
        for (size_t i = 0; i < count; i++) {
            Expr* elem = list->data.function.args[i];
            if (!ht_find(ht, elem)) {
                Expr* copy = expr_copy(elem);
                unique_args[unique_count++] = copy;
                ht_insert(ht, copy, 0);
            }
        }
        ht_free(ht, false);
    } else {
        // Fallback to O(N^2) for custom test
        for (size_t i = 0; i < count; i++) {
            Expr* elem = list->data.function.args[i];
            bool is_duplicate = false;
            for (size_t j = 0; j < unique_count; j++) {
                Expr* call_args[2] = { expr_copy(elem), expr_copy(unique_args[j]) };
                Expr* call = expr_new_function(expr_copy(test), call_args, 2);
                Expr* eval_res = evaluate(call);
                if (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol == SYM_True) {
                    is_duplicate = true;
                    expr_free(eval_res);
                    expr_free(call);
                    break;
                }
                expr_free(eval_res);
                expr_free(call);
            }
            if (!is_duplicate) {
                unique_args[unique_count++] = expr_copy(elem);
            }
        }
    }
    
    Expr* result = expr_new_function(expr_copy(list->data.function.head), unique_args, unique_count);
    if (unique_args) free(unique_args);
    
    return result;
}


Expr* builtin_split(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    
    Expr* list = res->data.function.args[0];
    Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    
    if (list->type != EXPR_FUNCTION) return expr_copy(list);
    
    size_t count = list->data.function.arg_count;
    if (count == 0) return expr_copy(list);
    
    Expr** result_runs = malloc(sizeof(Expr*) * count);
    size_t num_runs = 0;
    
    size_t run_start = 0;
    for (size_t i = 1; i <= count; i++) {
        bool split_here = (i == count);
        if (!split_here) {
            Expr* prev = list->data.function.args[i-1];
            Expr* curr = list->data.function.args[i];
            bool identical = false;
            if (test == NULL) {
                identical = expr_eq(prev, curr);
            } else {
                Expr* call_args[2] = { expr_copy(prev), expr_copy(curr) };
                Expr* call = expr_new_function(expr_copy(test), call_args, 2);
                Expr* eval_res = evaluate(call);
                if (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol == SYM_True) {
                    identical = true;
                }
                expr_free(eval_res);
                expr_free(call);
            }
            if (!identical) split_here = true;
        }
        
        if (split_here) {
            size_t run_len = i - run_start;
            Expr** run_args = malloc(sizeof(Expr*) * run_len);
            for (size_t j = 0; j < run_len; j++) {
                run_args[j] = expr_copy(list->data.function.args[run_start + j]);
            }
            result_runs[num_runs++] = expr_new_function(expr_copy(list->data.function.head), run_args, run_len);
            free(run_args);
            run_start = i;
        }
    }
    
    Expr* result = expr_new_function(expr_copy(list->data.function.head), result_runs, num_runs);
    free(result_runs);
    return result;
}

typedef struct {
    Expr* element;
    int64_t count;
    size_t first_index;
} CommonestItem;

static int compare_commonest_items_desc(const void* a, const void* b) {
    const CommonestItem* item_a = (const CommonestItem*)a;
    const CommonestItem* item_b = (const CommonestItem*)b;
    if (item_a->count != item_b->count) {
        return (item_b->count > item_a->count) ? 1 : -1;
    }
    return (item_a->first_index > item_b->first_index) ? 1 : -1;
}

static int compare_commonest_items_index(const void* a, const void* b) {
    const CommonestItem* item_a = (const CommonestItem*)a;
    const CommonestItem* item_b = (const CommonestItem*)b;
    if (item_a->first_index == item_b->first_index) return 0;
    return (item_a->first_index > item_b->first_index) ? 1 : -1;
}

Expr* builtin_commonest(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    
    Expr* list = res->data.function.args[0];
    if (list->type != EXPR_FUNCTION) return expr_new_function(expr_new_symbol("List"), NULL, 0);
    
    size_t count = list->data.function.arg_count;
    if (count == 0) return expr_new_function(expr_new_symbol("List"), NULL, 0);

    Expr* n_arg = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    int64_t n = -1;
    bool n_upto = false;
    if (n_arg) {
        if (n_arg->type == EXPR_INTEGER) {
            n = n_arg->data.integer;
        } else if (n_arg->type == EXPR_FUNCTION && n_arg->data.function.head->type == EXPR_SYMBOL && 
                   n_arg->data.function.head->data.symbol == SYM_UpTo && n_arg->data.function.arg_count == 1) {
            if (n_arg->data.function.args[0]->type == EXPR_INTEGER) {
                n = n_arg->data.function.args[0]->data.integer;
                n_upto = true;
            } else return NULL;
        } else return NULL;
    }

    // Tally
    Expr** unique_elems = malloc(sizeof(Expr*) * count);
    int64_t* multiplicities = malloc(sizeof(int64_t) * count);
    size_t unique_count = 0;

    HashTable* ht = ht_create(count * 2 + 1);
    for (size_t i = 0; i < count; i++) {
        Expr* elem = list->data.function.args[i];
        HashNode* node = ht_find(ht, elem);
        if (node) {
            multiplicities[node->index]++;
        } else {
            unique_elems[unique_count] = expr_copy(elem);
            multiplicities[unique_count] = 1;
            ht_insert(ht, unique_elems[unique_count], unique_count);
            unique_count++;
        }
    }
    ht_free(ht, false);

    CommonestItem* items = malloc(sizeof(CommonestItem) * unique_count);
    for (size_t i = 0; i < unique_count; i++) {
        items[i].element = unique_elems[i];
        items[i].count = multiplicities[i];
        items[i].first_index = i;
    }
    free(multiplicities);
    free(unique_elems);

    // Sort by count DESC, first_index ASC
    qsort(items, unique_count, sizeof(CommonestItem), compare_commonest_items_desc);

    size_t target_n;
    if (n == -1) {
        // Just the most common ones (highest count)
        int64_t max_count = items[0].count;
        target_n = 0;
        while (target_n < unique_count && items[target_n].count == max_count) {
            target_n++;
        }
    } else {
        if (n < 0) n = 0;
        if ((size_t)n > unique_count) {
            if (!n_upto) {
                printf("Commonest::dstlms: The requested number of elements %" PRId64 " is greater than the number of distinct elements %zu. Only %zu elements will be returned.\n", n, unique_count, unique_count);
            }
            target_n = unique_count;
        } else {
            target_n = (size_t)n;
        }
    }

    // Sort target_n items by first_index ASC to preserve original order
    if (target_n > 0) {
        qsort(items, target_n, sizeof(CommonestItem), compare_commonest_items_index);
    }

    Expr** result_args = malloc(sizeof(Expr*) * target_n);
    for (size_t i = 0; i < target_n; i++) {
        result_args[i] = items[i].element;
    }
    // Free unused elements
    for (size_t i = target_n; i < unique_count; i++) {
        expr_free(items[i].element);
    }
    free(items);

    return expr_new_function(expr_new_symbol("List"), result_args, target_n);
}

void list_init(void) {
    symtab_add_builtin("Table", builtin_table);
    symtab_add_builtin("Range", builtin_range);
    symtab_add_builtin("Array", builtin_array);
    symtab_add_builtin("Take", builtin_take);
    symtab_add_builtin("Drop", builtin_drop);
    symtab_add_builtin("Flatten", builtin_flatten);
    symtab_add_builtin("Partition", builtin_partition);
    symtab_add_builtin("RotateLeft", builtin_rotateleft);
    symtab_add_builtin("RotateRight", builtin_rotateright);
    symtab_add_builtin("Reverse", builtin_reverse);
    symtab_add_builtin("Join", builtin_join);
    symtab_get_def("Join")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Join",
        "Join[list1, list2, ...]\n"
        "\tConcatenates lists or other expressions that share the same head.\n"
        "Join[list1, list2, ..., n]\n"
        "\tJoins the objects at level n in each of the lists.\n"
        "\tHandles ragged arrays by concatenating successive elements at level n.");
    symtab_add_builtin("Transpose", builtin_transpose);
    symtab_add_builtin("ConjugateTranspose", builtin_conjugate_transpose);
    symtab_add_builtin("Tally", builtin_tally);
    symtab_add_builtin("Union", builtin_union);
    symtab_add_builtin("DeleteDuplicates", builtin_deleteduplicates);
    symtab_add_builtin("Split", builtin_split);
    symtab_add_builtin("Total", builtin_total);
    symtab_add_builtin("Accumulate", builtin_accumulate);
    symtab_add_builtin("Differences", builtin_differences);
    symtab_add_builtin("Ratios", builtin_ratios);
    symtab_add_builtin("Commonest", builtin_commonest);
    symtab_add_builtin("Min", builtin_min);
    symtab_add_builtin("Max", builtin_max);
    symtab_add_builtin("ListQ", builtin_listq);
    symtab_add_builtin("VectorQ", builtin_vectorq);
    symtab_add_builtin("MatrixQ", builtin_matrixq);
    symtab_add_builtin("HermitianMatrixQ", builtin_hermitian_matrix_q);
    symtab_add_builtin("SymmetricMatrixQ", builtin_symmetric_matrix_q);
    symtab_add_builtin("SquareMatrixQ", builtin_square_matrix_q);
    symtab_add_builtin("DiagonalMatrixQ", builtin_diagonal_matrix_q);
    symtab_add_builtin("UpperTriangularMatrixQ", builtin_upper_triangular_matrix_q);

    symtab_get_def("Table")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_get_def("Range")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Array")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Take")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Drop")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Flatten")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Partition")->attributes |= ATTR_PROTECTED;
    symtab_get_def("RotateLeft")->attributes |= ATTR_PROTECTED;
    symtab_get_def("RotateRight")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Reverse")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Transpose")->attributes |= ATTR_PROTECTED;
    symtab_get_def("ConjugateTranspose")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Tally")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Union")->attributes |= ATTR_FLAT | ATTR_ONEIDENTITY | ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_get_def("DeleteDuplicates")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Split")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Total")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Accumulate")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Differences")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Ratios")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Commonest")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Min")->attributes |= ATTR_FLAT | ATTR_NUMERICFUNCTION | ATTR_ONEIDENTITY | ATTR_ORDERLESS | ATTR_PROTECTED;
    symtab_get_def("Max")->attributes |= ATTR_FLAT | ATTR_NUMERICFUNCTION | ATTR_ONEIDENTITY | ATTR_ORDERLESS | ATTR_PROTECTED;
    symtab_get_def("ListQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("VectorQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("MatrixQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("HermitianMatrixQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("SymmetricMatrixQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("SquareMatrixQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("DiagonalMatrixQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("UpperTriangularMatrixQ")->attributes |= ATTR_PROTECTED;

    symtab_set_docstring("Total", "Total[list]\n\tgives the total of the elements in list.\nTotal[list, n]\n\ttotals all elements down to level n.\nTotal[list, {n}]\n\ttotals elements at level n.\nTotal[list, {n1, n2}]\n\ttotals elements at levels n1 through n2.");
}

static bool is_overflow(Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol == SYM_Overflow;
}

static bool is_listq(Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol == SYM_List;
}

static bool is_infinity(Expr* e) {
    return e->type == EXPR_SYMBOL && e->data.symbol == SYM_Infinity;
}

static bool is_minus_infinity(Expr* e) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Times &&
        e->data.function.arg_count == 2) {
        Expr* a1 = e->data.function.args[0];
        Expr* a2 = e->data.function.args[1];
        if (a1->type == EXPR_INTEGER && a1->data.integer == -1 && is_infinity(a2)) return true;
        if (a2->type == EXPR_INTEGER && a2->data.integer == -1 && is_infinity(a1)) return true;
    }
    return false;
}

static Expr* make_minus_infinity(void) {
    Expr* args[2] = { expr_new_integer(-1), expr_new_symbol("Infinity") };
    return expr_new_function(expr_new_symbol("Times"), args, 2);
}

static bool is_real_numeric(Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT) return true;
    if (is_rational(e, NULL, NULL)) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    Expr* re, *im;
    if (is_complex(e, &re, &im)) {
        if (im->type == EXPR_INTEGER && im->data.integer == 0) return true;
        if (im->type == EXPR_REAL && im->data.real == 0.0) return true;
        if (im->type == EXPR_BIGINT && mpz_sgn(im->data.bigint) == 0) return true;
#ifdef USE_MPFR
        if (im->type == EXPR_MPFR && mpfr_zero_p(im->data.mpfr)) return true;
#endif
    }
    return false;
}

static int64_t get_depth_for_total(Expr* e) {
    if (e->type != EXPR_FUNCTION) return 1;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol;
        if (h == SYM_Rational || h == SYM_Complex) return 1;
    }
    int64_t max_d = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        int64_t d = get_depth_for_total(e->data.function.args[i]);
        if (d > max_d) max_d = d;
    }
    return 1 + max_d;
}

static Expr* total_at_exactly_level_k(Expr* e, int64_t k) {
    if (k <= 0) return expr_copy(e);
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    if (k == 1) {
        size_t count = e->data.function.arg_count;
        if (count == 0) return expr_new_integer(0);
        Expr** plus_args = malloc(sizeof(Expr*) * count);
        for (size_t i = 0; i < count; i++) plus_args[i] = expr_copy(e->data.function.args[i]);
        Expr* plus_expr = expr_new_function(expr_new_symbol("Plus"), plus_args, count);
        free(plus_args);
        Expr* res = evaluate(plus_expr);
        expr_free(plus_expr);
        return res;
    } else {
        size_t count = e->data.function.arg_count;
        Expr** new_args = malloc(sizeof(Expr*) * count);
        for (size_t i = 0; i < count; i++) {
            new_args[i] = total_at_exactly_level_k(e->data.function.args[i], k - 1);
        }
        Expr* res = expr_new_function(expr_copy(e->data.function.head), new_args, count);
        free(new_args);
        return res;
    }
}

Expr* builtin_total(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    
    Expr* list = res->data.function.args[0];
    Expr* level_spec = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;

    int64_t n1 = 1, n2 = 1;
    int64_t depth = get_depth_for_total(list);

    if (level_spec) {
        if (level_spec->type == EXPR_INTEGER) {
            n1 = 1;
            n2 = level_spec->data.integer;
            if (n2 < 0) n2 = depth + n2;
        } else if (is_infinity(level_spec)) {
            n1 = 1;
            n2 = depth - 1;
        } else if (is_listq(level_spec)) {
            if (level_spec->data.function.arg_count == 1) {
                Expr* arg = level_spec->data.function.args[0];
                if (arg->type == EXPR_INTEGER) {
                    n1 = n2 = arg->data.integer;
                    if (n1 < 0) n1 = n2 = depth + n1;
                } else if (is_infinity(arg)) {
                    n1 = n2 = depth - 1;
                } else return NULL;
            } else if (level_spec->data.function.arg_count == 2) {
                Expr* arg1 = level_spec->data.function.args[0];
                Expr* arg2 = level_spec->data.function.args[1];
                if (arg1->type == EXPR_INTEGER) {
                    n1 = arg1->data.integer;
                    if (n1 < 0) n1 = depth + n1;
                } else return NULL;
                
                if (arg2->type == EXPR_INTEGER) {
                    n2 = arg2->data.integer;
                    if (n2 < 0) n2 = depth + n2;
                } else if (is_infinity(arg2)) {
                    n2 = depth - 1;
                } else return NULL;
            } else return NULL;
        } else return NULL;
    }
    
    if (n1 > n2) return expr_copy(list);
    if (n1 < 1) n1 = 1;
    if (n2 >= depth) n2 = depth - 1;

    Expr* current = expr_copy(list);
    for (int64_t k = n2; k >= n1; k--) {
        Expr* next = total_at_exactly_level_k(current, k);
        expr_free(current);
        current = next;
    }

    return current;
}

/* Accumulate[list] returns the running cumulative totals of list, with the
 * same head and the same length as the input. The intermediate sums are
 * built as Plus[acc, next] and reduced by the evaluator so that integers,
 * arbitrary-precision bignums, machine doubles, lists (matrix columns via
 * Listable Plus), and symbolic expressions all combine correctly.
 *
 * An optional second argument of the form Method -> "CompensatedSummation"
 * triggers Kahan compensated summation when every element reduces to a
 * machine double (EXPR_REAL, EXPR_INTEGER, or EXPR_BIGINT). For other
 * inputs the option is silently ignored and the symbolic accumulation
 * is used. */
static bool accumulate_is_compensated_method(Expr* opt) {
    if (opt->type != EXPR_FUNCTION) return false;
    if (opt->data.function.head->type != EXPR_SYMBOL) return false;
    const char* hd = opt->data.function.head->data.symbol;
    if ((hd != SYM_Rule && hd != SYM_RuleDelayed) ||
        opt->data.function.arg_count != 2) return false;
    Expr* lhs = opt->data.function.args[0];
    Expr* rhs = opt->data.function.args[1];
    if (lhs->type != EXPR_SYMBOL || lhs->data.symbol != SYM_Method) return false;
    if (rhs->type != EXPR_STRING) return false;
    return strcmp(rhs->data.string, "CompensatedSummation") == 0;
}

static bool accumulate_to_double(Expr* e, double* out) {
    if (e->type == EXPR_REAL)    { *out = e->data.real;             return true; }
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer;  return true; }
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint); return true; }
    return false;
}

Expr* builtin_accumulate(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    Expr* lst = res->data.function.args[0];
    if (lst->type != EXPR_FUNCTION) return NULL;

    bool kahan_requested = false;
    if (argc == 2) {
        if (!accumulate_is_compensated_method(res->data.function.args[1])) return NULL;
        kahan_requested = true;
    }

    Expr* head = lst->data.function.head;
    size_t n = lst->data.function.arg_count;

    if (n == 0) return expr_copy(lst);

    /* Kahan compensated summation in double precision when every element
     * is a machine number. */
    if (kahan_requested) {
        bool all_numeric = true;
        for (size_t i = 0; i < n; i++) {
            double tmp;
            if (!accumulate_to_double(lst->data.function.args[i], &tmp)) {
                all_numeric = false;
                break;
            }
        }
        if (all_numeric) {
            Expr** out = malloc(sizeof(Expr*) * n);
            if (!out) return NULL;
            double sum = 0.0;
            double c = 0.0;
            for (size_t i = 0; i < n; i++) {
                double x = 0.0;
                accumulate_to_double(lst->data.function.args[i], &x);
                double y = x - c;
                double t = sum + y;
                c = (t - sum) - y;
                sum = t;
                out[i] = expr_new_real(sum);
            }
            Expr* result = expr_new_function(expr_copy(head), out, n);
            free(out);
            return result;
        }
        /* Mixed/symbolic input: fall through to the standard accumulator. */
    }

    Expr** out = malloc(sizeof(Expr*) * n);
    if (!out) return NULL;

    out[0] = expr_copy(lst->data.function.args[0]);
    for (size_t i = 1; i < n; i++) {
        Expr** plus_args = malloc(sizeof(Expr*) * 2);
        if (!plus_args) {
            for (size_t j = 0; j < i; j++) expr_free(out[j]);
            free(out);
            return NULL;
        }
        plus_args[0] = expr_copy(out[i - 1]);
        plus_args[1] = expr_copy(lst->data.function.args[i]);
        Expr* plus_expr = expr_new_function(expr_new_symbol("Plus"), plus_args, 2);
        free(plus_args);
        out[i] = evaluate(plus_expr);
        expr_free(plus_expr);
    }

    Expr* result = expr_new_function(expr_copy(head), out, n);
    free(out);
    return result;
}

/* ---- Differences ----------------------------------------------------------
 *
 * Differences[list]              successive (first) differences of list.
 * Differences[list, n]           n-fold first differences (length l - n).
 * Differences[list, n, s]        n-fold differences of elements step s apart
 *                                (length l - n|s|).
 * Differences[list, {n1, n2,..}] successive n_k-th differences at level k of a
 *                                nested list / array.
 *
 * The element-wise subtraction list[[i+s]] - list[[i]] is built as
 * Subtract[hi, lo] and reduced by the evaluator, so integers, bignums, machine
 * doubles, symbolic terms, and whole sublists (matrix rows, via the Listable
 * Plus/Times that Subtract rewrites to) all combine correctly. Because row
 * subtraction threads element-wise, Differences[m, n] on a matrix m takes
 * differences of successive rows, matching Differences[m, {n, 0}].
 */

/* Returns the evaluated expression (minuend - subtrahend). The operands are
 * copied; ownership of the result transfers to the caller. */
static Expr* diff_minus(Expr* minuend, Expr* subtrahend) {
    Expr** args = malloc(sizeof(Expr*) * 2);
    args[0] = expr_copy(minuend);
    args[1] = expr_copy(subtrahend);
    Expr* sub = expr_new_function(expr_new_symbol("Subtract"), args, 2);
    free(args);
    return eval_and_free(sub);
}

/* One difference pass with nonzero integer step s over the top level of lst
 * (an EXPR_FUNCTION). The result keeps lst's head. For s > 0 the i-th element
 * is elem[i+s] - elem[i]; for s < 0 it is elem[i] - elem[i+|s|]. A list whose
 * length does not exceed |s| yields the empty list. */
static Expr* diff_once(Expr* lst, int64_t s) {
    Expr* head = lst->data.function.head;
    size_t n = lst->data.function.arg_count;
    int64_t a = (s < 0) ? -s : s;

    if ((int64_t)n <= a) {
        return expr_new_function(expr_copy(head), NULL, 0);
    }

    size_t outn = n - (size_t)a;
    Expr** out = malloc(sizeof(Expr*) * outn);
    for (size_t i = 0; i < outn; i++) {
        Expr* lo = lst->data.function.args[i];
        Expr* hi = lst->data.function.args[i + (size_t)a];
        out[i] = (s > 0) ? diff_minus(hi, lo) : diff_minus(lo, hi);
    }
    Expr* result = expr_new_function(expr_copy(head), out, outn);
    free(out);
    return result;
}

/* Apply diff_once with step s a total of n times. Returns a fresh expression
 * (a copy of lst when n == 0). */
static Expr* diff_n_step(Expr* lst, int64_t n, int64_t s) {
    Expr* cur = expr_copy(lst);
    for (int64_t k = 0; k < n; k++) {
        Expr* nxt = diff_once(cur, s);
        expr_free(cur);
        cur = nxt;
    }
    return cur;
}

/* Multidimensional differences: apply levels[0] step-1 first differences at the
 * top level, then recurse into each element with the remaining levels. */
static Expr* diff_levels(Expr* lst, const int64_t* levels, size_t num) {
    if (num == 0 || lst->type != EXPR_FUNCTION) {
        return expr_copy(lst);
    }

    Expr* cur = diff_n_step(lst, levels[0], 1);
    if (num == 1) {
        return cur;
    }

    size_t m = cur->data.function.arg_count;
    Expr* head = cur->data.function.head;
    Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
    for (size_t i = 0; i < m; i++) {
        out[i] = diff_levels(cur->data.function.args[i], levels + 1, num - 1);
    }
    Expr* result = expr_new_function(expr_copy(head), out, m);
    free(out);
    expr_free(cur);
    return result;
}

Expr* builtin_differences(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return NULL;

    Expr* lst = res->data.function.args[0];
    if (lst->type != EXPR_FUNCTION) return NULL;

    /* Differences[list, {n1, n2, ...}] — per-level differences. */
    if (argc == 2 && is_listq(res->data.function.args[1])) {
        Expr* spec = res->data.function.args[1];
        size_t num = spec->data.function.arg_count;
        if (num == 0) return expr_copy(lst);

        int64_t* levels = malloc(sizeof(int64_t) * num);
        for (size_t i = 0; i < num; i++) {
            Expr* e = spec->data.function.args[i];
            if (e->type != EXPR_INTEGER || e->data.integer < 0) {
                free(levels);
                return NULL;
            }
            levels[i] = e->data.integer;
        }
        Expr* result = diff_levels(lst, levels, num);
        free(levels);
        return result;
    }

    /* Differences[list], Differences[list, n], Differences[list, n, s]. */
    int64_t n = 1, s = 1;
    if (argc >= 2) {
        Expr* an = res->data.function.args[1];
        if (an->type != EXPR_INTEGER || an->data.integer < 0) return NULL;
        n = an->data.integer;
    }
    if (argc == 3) {
        Expr* as = res->data.function.args[2];
        if (as->type != EXPR_INTEGER || as->data.integer == 0) return NULL;
        s = as->data.integer;
    }

    return diff_n_step(lst, n, s);
}

/* ---- Ratios ----------------------------------------------------------------
 *
 * Ratios[list]                   successive (first) ratios of list.
 * Ratios[list, n]                n-fold iterated ratios (length l - n).
 * Ratios[list, {n1, n2, ...}]    successive n_k-th ratios at level k of a
 *                                nested list / array.
 *
 * Ratios is the multiplicative analog of Differences: the element-wise
 * subtraction list[[i+1]] - list[[i]] becomes the division list[[i+1]] /
 * list[[i]], built as Divide[hi, lo] and reduced by the evaluator. Integers,
 * bignums (-> exact Rationals), machine doubles, symbolic terms, and whole
 * sublists (matrix rows, via the Listable Power/Times that Divide rewrites to)
 * all combine correctly. Because row division threads element-wise,
 * Ratios[m, n] on a matrix m takes ratios of successive rows within each
 * column, matching Ratios[m, {n, 0}]. FoldList[Times, x, Ratios[list]] inverts
 * Ratios.
 */

/* Returns the evaluated expression (numerator / denominator). The operands are
 * copied; ownership of the result transfers to the caller. */
static Expr* ratio_divide(Expr* numerator, Expr* denominator) {
    Expr** args = malloc(sizeof(Expr*) * 2);
    args[0] = expr_copy(numerator);
    args[1] = expr_copy(denominator);
    Expr* div = expr_new_function(expr_new_symbol("Divide"), args, 2);
    free(args);
    return eval_and_free(div);
}

/* One ratio pass over the top level of lst (an EXPR_FUNCTION). The result keeps
 * lst's head; the i-th element is elem[i+1] / elem[i]. A list of length <= 1
 * yields the empty list. */
static Expr* ratio_once(Expr* lst) {
    Expr* head = lst->data.function.head;
    size_t n = lst->data.function.arg_count;

    if (n <= 1) {
        return expr_new_function(expr_copy(head), NULL, 0);
    }

    size_t outn = n - 1;
    Expr** out = malloc(sizeof(Expr*) * outn);
    for (size_t i = 0; i < outn; i++) {
        Expr* lo = lst->data.function.args[i];
        Expr* hi = lst->data.function.args[i + 1];
        out[i] = ratio_divide(hi, lo);
    }
    Expr* result = expr_new_function(expr_copy(head), out, outn);
    free(out);
    return result;
}

/* Apply ratio_once a total of n times. Returns a fresh expression (a copy of
 * lst when n == 0). */
static Expr* ratio_n(Expr* lst, int64_t n) {
    Expr* cur = expr_copy(lst);
    for (int64_t k = 0; k < n; k++) {
        Expr* nxt = ratio_once(cur);
        expr_free(cur);
        cur = nxt;
    }
    return cur;
}

/* Multidimensional ratios: apply levels[0] first ratios at the top level, then
 * recurse into each element with the remaining levels. */
static Expr* ratio_levels(Expr* lst, const int64_t* levels, size_t num) {
    if (num == 0 || lst->type != EXPR_FUNCTION) {
        return expr_copy(lst);
    }

    Expr* cur = ratio_n(lst, levels[0]);
    if (num == 1) {
        return cur;
    }

    size_t m = cur->data.function.arg_count;
    Expr* head = cur->data.function.head;
    Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
    for (size_t i = 0; i < m; i++) {
        out[i] = ratio_levels(cur->data.function.args[i], levels + 1, num - 1);
    }
    Expr* result = expr_new_function(expr_copy(head), out, m);
    free(out);
    expr_free(cur);
    return result;
}

Expr* builtin_ratios(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    Expr* lst = res->data.function.args[0];
    if (lst->type != EXPR_FUNCTION) return NULL;

    /* Ratios[list, {n1, n2, ...}] — per-level ratios. */
    if (argc == 2 && is_listq(res->data.function.args[1])) {
        Expr* spec = res->data.function.args[1];
        size_t num = spec->data.function.arg_count;
        if (num == 0) return expr_copy(lst);

        int64_t* levels = malloc(sizeof(int64_t) * num);
        for (size_t i = 0; i < num; i++) {
            Expr* e = spec->data.function.args[i];
            if (e->type != EXPR_INTEGER || e->data.integer < 0) {
                free(levels);
                return NULL;
            }
            levels[i] = e->data.integer;
        }
        Expr* result = ratio_levels(lst, levels, num);
        free(levels);
        return result;
    }

    /* Ratios[list], Ratios[list, n]. */
    int64_t n = 1;
    if (argc == 2) {
        Expr* an = res->data.function.args[1];
        if (an->type != EXPR_INTEGER || an->data.integer < 0) return NULL;
        n = an->data.integer;
    }

    return ratio_n(lst, n);
}

Expr* builtin_listq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    if (is_listq(arg)) {
        return expr_new_symbol("True");
    }
    return expr_new_symbol("False");
}

Expr* builtin_vectorq(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count != 1 && res->data.function.arg_count != 2)) return NULL;
    Expr* arg = res->data.function.args[0];
    if (!is_listq(arg)) return expr_new_symbol("False");

    Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;

    for (size_t i = 0; i < arg->data.function.arg_count; i++) {
        Expr* elem = arg->data.function.args[i];
        if (test == NULL) {
            if (is_listq(elem)) return expr_new_symbol("False");
        } else {
            Expr* call_args[1] = { expr_copy(elem) };
            Expr* call = expr_new_function(expr_copy(test), call_args, 1);
            Expr* eval_res = evaluate(call);
            bool is_true = (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol == SYM_True);
            expr_free(eval_res);
            expr_free(call);
            if (!is_true) return expr_new_symbol("False");
        }
    }
    return expr_new_symbol("True");
}

Expr* builtin_matrixq(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count != 1 && res->data.function.arg_count != 2)) return NULL;
    Expr* arg = res->data.function.args[0];
    if (!is_listq(arg)) return expr_new_symbol("False");

    if (arg->data.function.arg_count == 0) return expr_new_symbol("False");

    Expr* test = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;

    size_t col_count = 0;
    bool first_row = true;

    for (size_t i = 0; i < arg->data.function.arg_count; i++) {
        Expr* row = arg->data.function.args[i];
        if (!is_listq(row)) return expr_new_symbol("False");

        if (first_row) {
            col_count = row->data.function.arg_count;
            first_row = false;
        } else {
            if (row->data.function.arg_count != col_count) return expr_new_symbol("False");
        }

        for (size_t j = 0; j < row->data.function.arg_count; j++) {
            Expr* elem = row->data.function.args[j];
            if (test == NULL) {
                if (is_listq(elem)) return expr_new_symbol("False");
            } else {
                Expr* call_args[1] = { expr_copy(elem) };
                Expr* call = expr_new_function(expr_copy(test), call_args, 1);
                Expr* eval_res = evaluate(call);
                bool is_true = (eval_res->type == EXPR_SYMBOL && eval_res->data.symbol == SYM_True);
                expr_free(eval_res);
                expr_free(call);
                if (!is_true) return expr_new_symbol("False");
            }
        }
    }

    return expr_new_symbol("True");
}

/* --- HermitianMatrixQ ----------------------------------------------------
 *
 * A matrix m is Hermitian iff m == ConjugateTranspose[m], equivalently
 * m[i,j] == Conjugate[m[j,i]] for every pair (i,j).  Diagonal entries
 * must satisfy m[i,i] == Conjugate[m[i,i]], i.e. be self-conjugate
 * (purely real for numeric inputs).
 *
 * The default test is "explicit" (structural) matching: we accept a pair
 * (a, b) as conjugate-mirrored when any of the following holds --
 *
 *   (1) a is Conjugate[c] with c structurally equal to b;
 *   (2) b is Conjugate[c] with c structurally equal to a;
 *   (3) evaluating Conjugate[b] yields a structurally.
 *
 * Branches (1) and (2) catch the symbolic patterns (Conjugate[a], a) and
 * (a, Conjugate[a]) because our Conjugate builtin does NOT fold
 * Conjugate[Conjugate[x]] back to x for symbolic x.  Branch (3) covers
 * the fully numeric case where Conjugate evaluates to a concrete value.
 *
 * Options:
 *   - SameTest -> f : pairs (a, b) are accepted iff f[a, Conjugate[b]]
 *     evaluates to True.  Conjugate[b] is computed once per pair so the
 *     user-supplied predicate sees the actual conjugated entry, mirroring
 *     Mathematica's documented behaviour.
 *   - Tolerance -> t : pairs (a, b) are accepted iff
 *     Abs[a - Conjugate[b]] <= t evaluates to True.
 *
 * SameTest and Tolerance defaults of Automatic fall through to the
 * structural test, which is correct for both symbolic and exact-numeric
 * matrices and degrades to bit-identical comparison for machine reals.
 */

/* True when `n` is an evaluated Conjugate[x] node. */
static bool is_conjugate_node(Expr* n) {
    return n->type == EXPR_FUNCTION &&
           n->data.function.head->type == EXPR_SYMBOL &&
           n->data.function.head->data.symbol == SYM_Conjugate &&
           n->data.function.arg_count == 1;
}

static Expr* eval_conjugate_of(Expr* e) {
    Expr** args = malloc(sizeof(Expr*) * 1);
    args[0] = expr_copy(e);
    Expr* call = expr_new_function(expr_new_symbol("Conjugate"), args, 1);
    free(args);
    return eval_and_free(call);
}

static bool hermitian_pair_structural(Expr* a, Expr* b) {
    /* (1) a == Conjugate[c] structurally and c == b. */
    if (is_conjugate_node(a) && expr_eq(a->data.function.args[0], b)) return true;
    /* (2) b == Conjugate[c] structurally and c == a. */
    if (is_conjugate_node(b) && expr_eq(b->data.function.args[0], a)) return true;
    /* (3) Evaluate Conjugate[b] and compare structurally. */
    Expr* cb = eval_conjugate_of(b);
    bool eq = expr_eq(a, cb);
    expr_free(cb);
    return eq;
}

static bool hermitian_pair_sametest(Expr* a, Expr* b, Expr* test) {
    /* Build test[a, Conjugate[b]] and check for True. */
    Expr* cb = eval_conjugate_of(b);
    Expr** args = malloc(sizeof(Expr*) * 2);
    args[0] = expr_copy(a);
    args[1] = cb;  /* ownership transferred into the call below */
    Expr* call = expr_new_function(expr_copy(test), args, 2);
    free(args);
    Expr* result = eval_and_free(call);
    bool ok = (result->type == EXPR_SYMBOL &&
               result->data.symbol == SYM_True);
    expr_free(result);
    return ok;
}

static bool hermitian_pair_tolerance(Expr* a, Expr* b, Expr* tol) {
    /* Build LessEqual[Abs[a - Conjugate[b]], tol] and check for True. */
    Expr* cb = eval_conjugate_of(b);
    Expr** sub_args = malloc(sizeof(Expr*) * 2);
    sub_args[0] = expr_copy(a);
    sub_args[1] = cb;
    Expr* diff = expr_new_function(expr_new_symbol("Subtract"), sub_args, 2);
    free(sub_args);
    Expr* diff_e = eval_and_free(diff);

    Expr** abs_args = malloc(sizeof(Expr*) * 1);
    abs_args[0] = diff_e;
    Expr* abs_call = expr_new_function(expr_new_symbol("Abs"), abs_args, 1);
    free(abs_args);
    Expr* abs_e = eval_and_free(abs_call);

    Expr** le_args = malloc(sizeof(Expr*) * 2);
    le_args[0] = abs_e;
    le_args[1] = expr_copy(tol);
    Expr* le_call = expr_new_function(expr_new_symbol("LessEqual"), le_args, 2);
    free(le_args);
    Expr* le_e = eval_and_free(le_call);
    bool ok = (le_e->type == EXPR_SYMBOL &&
               le_e->data.symbol == SYM_True);
    expr_free(le_e);
    return ok;
}

Expr* builtin_hermitian_matrix_q(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    Expr* m = res->data.function.args[0];

    /* Parse options.  Each trailing arg must be a Rule with one of the
     * recognised option names; any unrecognised option causes us to
     * leave the call unevaluated so user-typed errors surface. */
    Expr* same_test = NULL;
    Expr* tolerance = NULL;
    for (size_t i = 1; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (!(opt->type == EXPR_FUNCTION &&
              opt->data.function.head->type == EXPR_SYMBOL &&
              opt->data.function.head->data.symbol == SYM_Rule &&
              opt->data.function.arg_count == 2 &&
              opt->data.function.args[0]->type == EXPR_SYMBOL)) {
            return NULL;
        }
        const char* name = opt->data.function.args[0]->data.symbol;
        Expr* val = opt->data.function.args[1];
        if (name == SYM_SameTest) {
            /* Automatic falls through to the structural test. */
            if (!(val->type == EXPR_SYMBOL && val->data.symbol == SYM_Automatic)) {
                same_test = val;
            }
        } else if (name == SYM_Tolerance) {
            if (!(val->type == EXPR_SYMBOL && val->data.symbol == SYM_Automatic)) {
                tolerance = val;
            }
        } else {
            return NULL;
        }
    }

    /* Must be a non-empty square List of Lists with no deeper nesting. */
    if (!is_listq(m)) return expr_new_symbol("False");
    size_t n = m->data.function.arg_count;
    if (n == 0) return expr_new_symbol("False");
    for (size_t i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        if (!is_listq(row)) return expr_new_symbol("False");
        if (row->data.function.arg_count != n) return expr_new_symbol("False");
        for (size_t j = 0; j < n; j++) {
            if (is_listq(row->data.function.args[j])) {
                return expr_new_symbol("False");
            }
        }
    }

    /* Walk the upper triangle (including diagonal) and check each pair
     * against the chosen predicate.  Walking i <= j is sufficient since
     * the pair test is symmetric under (i,j) <-> (j,i): a == Conj[b]
     * iff b == Conj[a]. */
    for (size_t i = 0; i < n; i++) {
        Expr** row_i = m->data.function.args[i]->data.function.args;
        for (size_t j = i; j < n; j++) {
            Expr* a = row_i[j];                                /* m[i,j] */
            Expr* b = m->data.function.args[j]->data.function.args[i]; /* m[j,i] */
            bool ok;
            if (same_test != NULL) {
                ok = hermitian_pair_sametest(a, b, same_test);
            } else if (tolerance != NULL) {
                ok = hermitian_pair_tolerance(a, b, tolerance);
            } else {
                ok = hermitian_pair_structural(a, b);
            }
            if (!ok) return expr_new_symbol("False");
        }
    }

    return expr_new_symbol("True");
}

/* --- SymmetricMatrixQ ----------------------------------------------------
 *
 * A matrix m is symmetric iff m == Transpose[m], i.e. m[i,j] == m[j,i]
 * for every pair (i, j).  This is the complex-symmetric notion -- no
 * conjugation is applied -- so a complex symmetric matrix need not be
 * Hermitian (and a Hermitian complex matrix need not be symmetric).
 *
 * The default test is structural equality via expr_eq.  Options mirror
 * HermitianMatrixQ:
 *   - SameTest -> f : pairs (a, b) accepted iff f[a, b] evaluates to True.
 *   - Tolerance -> t : pairs accepted iff Abs[a - b] <= t evaluates to True.
 *
 * SameTest/Tolerance defaults of Automatic fall through to the structural
 * test, which is correct for both symbolic and exact-numeric matrices and
 * degrades to bit-identical comparison for machine reals.  The diagonal is
 * always trivially symmetric (m[i,i] == m[i,i]) so we walk only the
 * strict upper triangle.
 */

static bool symmetric_pair_sametest(Expr* a, Expr* b, Expr* test) {
    Expr** args = malloc(sizeof(Expr*) * 2);
    args[0] = expr_copy(a);
    args[1] = expr_copy(b);
    Expr* call = expr_new_function(expr_copy(test), args, 2);
    free(args);
    Expr* result = eval_and_free(call);
    bool ok = (result->type == EXPR_SYMBOL &&
               result->data.symbol == SYM_True);
    expr_free(result);
    return ok;
}

static bool symmetric_pair_tolerance(Expr* a, Expr* b, Expr* tol) {
    /* Build LessEqual[Abs[a - b], tol] and check for True. */
    Expr** sub_args = malloc(sizeof(Expr*) * 2);
    sub_args[0] = expr_copy(a);
    sub_args[1] = expr_copy(b);
    Expr* diff = expr_new_function(expr_new_symbol("Subtract"), sub_args, 2);
    free(sub_args);
    Expr* diff_e = eval_and_free(diff);

    Expr** abs_args = malloc(sizeof(Expr*) * 1);
    abs_args[0] = diff_e;
    Expr* abs_call = expr_new_function(expr_new_symbol("Abs"), abs_args, 1);
    free(abs_args);
    Expr* abs_e = eval_and_free(abs_call);

    Expr** le_args = malloc(sizeof(Expr*) * 2);
    le_args[0] = abs_e;
    le_args[1] = expr_copy(tol);
    Expr* le_call = expr_new_function(expr_new_symbol("LessEqual"), le_args, 2);
    free(le_args);
    Expr* le_e = eval_and_free(le_call);
    bool ok = (le_e->type == EXPR_SYMBOL &&
               le_e->data.symbol == SYM_True);
    expr_free(le_e);
    return ok;
}

Expr* builtin_symmetric_matrix_q(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    Expr* m = res->data.function.args[0];

    /* Parse options.  Each trailing arg must be a Rule with one of the
     * recognised option names; any unrecognised option causes us to
     * leave the call unevaluated so user-typed errors surface. */
    Expr* same_test = NULL;
    Expr* tolerance = NULL;
    for (size_t i = 1; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (!(opt->type == EXPR_FUNCTION &&
              opt->data.function.head->type == EXPR_SYMBOL &&
              opt->data.function.head->data.symbol == SYM_Rule &&
              opt->data.function.arg_count == 2 &&
              opt->data.function.args[0]->type == EXPR_SYMBOL)) {
            return NULL;
        }
        const char* name = opt->data.function.args[0]->data.symbol;
        Expr* val = opt->data.function.args[1];
        if (name == SYM_SameTest) {
            if (!(val->type == EXPR_SYMBOL && val->data.symbol == SYM_Automatic)) {
                same_test = val;
            }
        } else if (name == SYM_Tolerance) {
            if (!(val->type == EXPR_SYMBOL && val->data.symbol == SYM_Automatic)) {
                tolerance = val;
            }
        } else {
            return NULL;
        }
    }

    /* Must be a non-empty square List of Lists with no deeper nesting. */
    if (!is_listq(m)) return expr_new_symbol("False");
    size_t n = m->data.function.arg_count;
    if (n == 0) return expr_new_symbol("False");
    for (size_t i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        if (!is_listq(row)) return expr_new_symbol("False");
        if (row->data.function.arg_count != n) return expr_new_symbol("False");
        for (size_t j = 0; j < n; j++) {
            if (is_listq(row->data.function.args[j])) {
                return expr_new_symbol("False");
            }
        }
    }

    /* Walk strict upper triangle (i < j); the diagonal is trivially
     * symmetric under any reasonable equality test that is reflexive. */
    for (size_t i = 0; i < n; i++) {
        Expr** row_i = m->data.function.args[i]->data.function.args;
        for (size_t j = i + 1; j < n; j++) {
            Expr* a = row_i[j];                                /* m[i,j] */
            Expr* b = m->data.function.args[j]->data.function.args[i]; /* m[j,i] */
            bool ok;
            if (same_test != NULL) {
                ok = symmetric_pair_sametest(a, b, same_test);
            } else if (tolerance != NULL) {
                ok = symmetric_pair_tolerance(a, b, tolerance);
            } else {
                ok = expr_eq(a, b);
            }
            if (!ok) return expr_new_symbol("False");
        }
    }

    return expr_new_symbol("True");
}

/* --- SquareMatrixQ -----------------------------------------------------
 *
 * A matrix m is square iff Dimensions[m] == {n, n}, i.e. it has the
 * same number of rows and columns.  Pure shape test -- no element
 * predicate or option is consulted.  Returns False for non-lists,
 * empty lists, lists whose elements are not all lists of the same
 * length, ragged matrices, rectangular matrices, and higher-rank
 * tensors (any entry that is itself a List).  Single-element matrices
 * {{x}} are square for any x (including symbolic x).
 *
 * Exactly one argument is accepted; any other count emits a
 * Mathematica-compatible `SquareMatrixQ::argx` diagnostic to stderr
 * and leaves the call unevaluated, matching the surface behaviour of
 * other wrong-arity builtins (cf. builtin_conjugate).
 */
Expr* builtin_square_matrix_q(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1) {
        size_t n = res->data.function.arg_count;
        fprintf(stderr,
                "SquareMatrixQ::argx: SquareMatrixQ called with %zu "
                "argument%s; 1 argument is expected.\n",
                n, n == 1 ? "" : "s");
        return NULL;
    }

    Expr* m = res->data.function.args[0];

    if (!is_listq(m)) return expr_new_symbol("False");
    size_t n = m->data.function.arg_count;
    if (n == 0) return expr_new_symbol("False");
    for (size_t i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        if (!is_listq(row)) return expr_new_symbol("False");
        if (row->data.function.arg_count != n) return expr_new_symbol("False");
        for (size_t j = 0; j < n; j++) {
            if (is_listq(row->data.function.args[j])) {
                return expr_new_symbol("False");
            }
        }
    }
    return expr_new_symbol("True");
}

/* --- DiagonalMatrixQ -----------------------------------------------------
 *
 * `DiagonalMatrixQ[m]` returns True iff every off-(main-diagonal) entry of
 * the matrix `m` is structurally zero; False otherwise.
 *
 * `DiagonalMatrixQ[m, k]` checks the k-th diagonal: positive k selects a
 * superdiagonal above the main diagonal, negative k selects a subdiagonal
 * below it. The matrix is "k-diagonal" iff every entry m[i,j] with
 * j - i != k is zero. Rectangular matrices are supported -- only shape and
 * the entry-zero predicate matter.
 *
 * `DiagonalMatrixQ[m, ..., Tolerance -> t]` relaxes the zero test so that
 * an entry e is considered zero iff `Abs[e] <= t` evaluates to True. This
 * mirrors the behaviour of HermitianMatrixQ / SymmetricMatrixQ.
 *
 * Diagnostics:
 *   - 0 args -> `DiagonalMatrixQ::argt` to stderr, call unevaluated.
 *   - >= 3 positional args (or a non-Rule arg where an option is expected)
 *     -> `DiagonalMatrixQ::nonopt` to stderr, call unevaluated.
 *
 * Shape rejections that return False (rather than unevaluated): non-list
 * input, the empty list `{}`, ragged rows, and any matrix whose entries
 * are themselves Lists (rank >= 3 tensor). `{{}, {}, ...}` (n x 0) is
 * vacuously diagonal and returns True.
 *
 * For symbolic entries the structural zero test is exact: only the
 * literal numeric zeros (Integer 0, Real 0.0, BigInt 0) count as zero.
 * Pure symbols, function calls, and nonzero literals are non-zero, so the
 * predicate is conservative -- it never proves a matrix is diagonal that
 * is not provably so under structural reasoning alone.
 */

static bool diag_entry_is_exact_zero(Expr* e) {
    if (e == NULL) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_REAL) return e->data.real == 0.0;
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint) == 0;
    return false;
}

static bool diag_entry_under_tolerance(Expr* e, Expr* tol) {
    /* Build LessEqual[Abs[e], tol] and require evaluation to True. */
    Expr** abs_args = malloc(sizeof(Expr*) * 1);
    abs_args[0] = expr_copy(e);
    Expr* abs_call = expr_new_function(expr_new_symbol("Abs"), abs_args, 1);
    free(abs_args);
    Expr* abs_e = eval_and_free(abs_call);

    Expr** le_args = malloc(sizeof(Expr*) * 2);
    le_args[0] = abs_e;
    le_args[1] = expr_copy(tol);
    Expr* le_call = expr_new_function(expr_new_symbol("LessEqual"), le_args, 2);
    free(le_args);
    Expr* le_e = eval_and_free(le_call);
    bool ok = (le_e->type == EXPR_SYMBOL && le_e->data.symbol == SYM_True);
    expr_free(le_e);
    return ok;
}

/* Print a Mathematica-compatible argt diagnostic and return NULL so the
 * evaluator leaves the call unevaluated, matching builtin_square_matrix_q
 * / builtin_conjugate's surface behaviour. */
static Expr* diag_emit_argt(size_t argc) {
    fprintf(stderr,
            "DiagonalMatrixQ::argt: DiagonalMatrixQ called with %zu "
            "argument%s; 1 or 2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/* Print a Mathematica-compatible nonopt diagnostic.  We mirror Wolfram's
 * surface text: "Options expected (instead of <expr>) beyond position 2".
 * `bad` is the offending non-Rule expression that broke option parsing. */
static Expr* diag_emit_nonopt(Expr* bad, Expr* res) {
    char* bad_str = expr_to_string(bad);
    char* call_str = expr_to_string(res);
    fprintf(stderr,
            "DiagonalMatrixQ::nonopt: Options expected (instead of %s) "
            "beyond position 2 in %s. An option must be a rule or a list "
            "of rules.\n",
            bad_str ? bad_str : "?", call_str ? call_str : "?");
    free(bad_str);
    free(call_str);
    return NULL;
}

Expr* builtin_diagonal_matrix_q(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return diag_emit_argt(0);

    Expr* m = res->data.function.args[0];

    /* Argument layout:
     *   args[0] = matrix m
     *   args[1] = optional integer k (default 0) OR first option Rule.
     *   args[2..] = options.
     */
    int64_t k = 0;
    size_t opt_start = 1;
    if (argc >= 2) {
        Expr* a1 = res->data.function.args[1];
        bool is_rule = (a1->type == EXPR_FUNCTION &&
                        a1->data.function.head->type == EXPR_SYMBOL &&
                        a1->data.function.head->data.symbol == SYM_Rule);
        if (is_rule) {
            /* k defaults to 0; options start at position 1. */
            opt_start = 1;
        } else {
            /* Position 1 must be an integer k. */
            if (a1->type != EXPR_INTEGER) {
                /* Non-integer, non-Rule -> treat as a bad option starting
                 * at position 2 (1-indexed). */
                return diag_emit_nonopt(a1, res);
            }
            k = a1->data.integer;
            opt_start = 2;
        }
    }

    /* Parse options. Only Tolerance is recognised; anything else is a
     * nonopt diagnostic (unknown-option rules left over would otherwise
     * silently no-op). */
    Expr* tolerance = NULL;
    Expr* last_bad = NULL;  /* report the LAST non-option encountered */
    for (size_t i = opt_start; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        bool is_rule = (opt->type == EXPR_FUNCTION &&
                        opt->data.function.head->type == EXPR_SYMBOL &&
                        opt->data.function.head->data.symbol == SYM_Rule &&
                        opt->data.function.arg_count == 2 &&
                        opt->data.function.args[0]->type == EXPR_SYMBOL);
        if (!is_rule) {
            last_bad = opt;
            continue;
        }
        const char* name = opt->data.function.args[0]->data.symbol;
        Expr* val = opt->data.function.args[1];
        if (name == SYM_Tolerance) {
            if (!(val->type == EXPR_SYMBOL && val->data.symbol == SYM_Automatic)) {
                tolerance = val;
            }
        } else {
            /* Unknown option name -> treat as a bad option. */
            last_bad = opt;
        }
    }
    if (last_bad != NULL) {
        return diag_emit_nonopt(last_bad, res);
    }

    /* Validate matrix shape: must be a List of Lists, all rows the same
     * length, no entries that are themselves Lists.  The empty list `{}`
     * (a vector with zero entries) is not a matrix and returns False;
     * `{{}, {}, ...}` (n x 0) is accepted as a vacuous matrix. */
    if (!is_listq(m)) return expr_new_symbol("False");
    size_t nrows = m->data.function.arg_count;
    if (nrows == 0) return expr_new_symbol("False");
    size_t ncols = 0;
    for (size_t i = 0; i < nrows; i++) {
        Expr* row = m->data.function.args[i];
        if (!is_listq(row)) return expr_new_symbol("False");
        size_t this_ncols = row->data.function.arg_count;
        if (i == 0) {
            ncols = this_ncols;
        } else if (this_ncols != ncols) {
            return expr_new_symbol("False");
        }
        for (size_t j = 0; j < this_ncols; j++) {
            if (is_listq(row->data.function.args[j])) {
                return expr_new_symbol("False");
            }
        }
    }

    /* Walk every cell; an entry off the k-th diagonal (j - i != k) must
     * be zero under the chosen predicate.  Entries on the diagonal are
     * unconstrained -- they may be any value, including symbolic. */
    for (size_t i = 0; i < nrows; i++) {
        Expr** row_i = m->data.function.args[i]->data.function.args;
        for (size_t j = 0; j < ncols; j++) {
            int64_t off = (int64_t)j - (int64_t)i;
            if (off == k) continue;  /* on the k-th diagonal -- ignore */
            Expr* entry = row_i[j];
            bool zero;
            if (tolerance != NULL) {
                zero = diag_entry_under_tolerance(entry, tolerance);
            } else {
                zero = diag_entry_is_exact_zero(entry);
            }
            if (!zero) return expr_new_symbol("False");
        }
    }

    return expr_new_symbol("True");
}

/* --- UpperTriangularMatrixQ --------------------------------------------
 *
 * `UpperTriangularMatrixQ[m]` returns True iff every entry of `m` strictly
 * below the main diagonal is zero; False otherwise.
 *
 * `UpperTriangularMatrixQ[m, k]` shifts the cut-off to the k-th diagonal:
 * every entry m[i,j] with `j - i < k` must be zero.  Equivalently, all
 * nonzero entries are confined to the region on or above the k-th
 * diagonal.  Positive k selects a superdiagonal above the main diagonal
 * (so the test becomes stricter); negative k selects a subdiagonal below
 * it (the test becomes more permissive).  Rectangular matrices are
 * supported.
 *
 * `UpperTriangularMatrixQ[m, ..., Tolerance -> t]` relaxes the zero test
 * so that an entry e is considered zero iff `Abs[e] <= t` evaluates to
 * True.  Mirrors the surface API of DiagonalMatrixQ / HermitianMatrixQ /
 * SymmetricMatrixQ.
 *
 * Diagnostics (Mathematica-compatible, to stderr):
 *   - 0 args -> `UpperTriangularMatrixQ::argt`, call left unevaluated.
 *   - >= 3 positional args (or a non-Rule arg where an option is
 *     expected) -> `UpperTriangularMatrixQ::nonopt`, call left
 *     unevaluated.
 *
 * Shape rejections that return False (rather than unevaluated): non-list
 * input, the empty list `{}`, ragged rows, and any matrix whose entries
 * are themselves Lists (rank >= 3 tensor).  `{{}, {}, ...}` (n x 0) is
 * vacuously upper-triangular and returns True.
 *
 * For symbolic entries the structural zero test is exact: only literal
 * numeric zeros (Integer 0, Real 0.0, BigInt 0) count as zero.  Symbolic
 * sub-diagonal entries cause the matrix to be rejected, so the predicate
 * is conservative.
 */

static Expr* utri_emit_argt(size_t argc) {
    fprintf(stderr,
            "UpperTriangularMatrixQ::argt: UpperTriangularMatrixQ called "
            "with %zu argument%s; 1 or 2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

static Expr* utri_emit_nonopt(Expr* bad, Expr* res) {
    char* bad_str = expr_to_string(bad);
    char* call_str = expr_to_string(res);
    fprintf(stderr,
            "UpperTriangularMatrixQ::nonopt: Options expected (instead of "
            "%s) beyond position 2 in %s. An option must be a rule or a "
            "list of rules.\n",
            bad_str ? bad_str : "?", call_str ? call_str : "?");
    free(bad_str);
    free(call_str);
    return NULL;
}

Expr* builtin_upper_triangular_matrix_q(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return utri_emit_argt(0);

    Expr* m = res->data.function.args[0];

    /* Argument layout, mirroring DiagonalMatrixQ:
     *   args[0]    = matrix m
     *   args[1]    = optional integer k (default 0) OR first option Rule.
     *   args[2..]  = options.
     */
    int64_t k = 0;
    size_t opt_start = 1;
    if (argc >= 2) {
        Expr* a1 = res->data.function.args[1];
        bool is_rule = (a1->type == EXPR_FUNCTION &&
                        a1->data.function.head->type == EXPR_SYMBOL &&
                        a1->data.function.head->data.symbol == SYM_Rule);
        if (is_rule) {
            opt_start = 1;
        } else {
            if (a1->type != EXPR_INTEGER) {
                return utri_emit_nonopt(a1, res);
            }
            k = a1->data.integer;
            opt_start = 2;
        }
    }

    /* Parse options.  Only Tolerance is recognised. */
    Expr* tolerance = NULL;
    Expr* last_bad = NULL;
    for (size_t i = opt_start; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        bool is_rule = (opt->type == EXPR_FUNCTION &&
                        opt->data.function.head->type == EXPR_SYMBOL &&
                        opt->data.function.head->data.symbol == SYM_Rule &&
                        opt->data.function.arg_count == 2 &&
                        opt->data.function.args[0]->type == EXPR_SYMBOL);
        if (!is_rule) {
            last_bad = opt;
            continue;
        }
        const char* name = opt->data.function.args[0]->data.symbol;
        Expr* val = opt->data.function.args[1];
        if (name == SYM_Tolerance) {
            if (!(val->type == EXPR_SYMBOL && val->data.symbol == SYM_Automatic)) {
                tolerance = val;
            }
        } else {
            last_bad = opt;
        }
    }
    if (last_bad != NULL) {
        return utri_emit_nonopt(last_bad, res);
    }

    /* Validate matrix shape -- same gate as DiagonalMatrixQ. */
    if (!is_listq(m)) return expr_new_symbol("False");
    size_t nrows = m->data.function.arg_count;
    if (nrows == 0) return expr_new_symbol("False");
    size_t ncols = 0;
    for (size_t i = 0; i < nrows; i++) {
        Expr* row = m->data.function.args[i];
        if (!is_listq(row)) return expr_new_symbol("False");
        size_t this_ncols = row->data.function.arg_count;
        if (i == 0) {
            ncols = this_ncols;
        } else if (this_ncols != ncols) {
            return expr_new_symbol("False");
        }
        for (size_t j = 0; j < this_ncols; j++) {
            if (is_listq(row->data.function.args[j])) {
                return expr_new_symbol("False");
            }
        }
    }

    /* Every entry strictly below the k-th diagonal (j - i < k) must be
     * zero under the chosen predicate.  Entries on or above the k-th
     * diagonal are unconstrained. */
    for (size_t i = 0; i < nrows; i++) {
        Expr** row_i = m->data.function.args[i]->data.function.args;
        for (size_t j = 0; j < ncols; j++) {
            int64_t off = (int64_t)j - (int64_t)i;
            if (off >= k) continue;  /* on or above k-th diagonal -- ignore */
            Expr* entry = row_i[j];
            bool zero;
            if (tolerance != NULL) {
                zero = diag_entry_under_tolerance(entry, tolerance);
            } else {
                zero = diag_entry_is_exact_zero(entry);
            }
            if (!zero) return expr_new_symbol("False");
        }
    }

    return expr_new_symbol("True");
}
Expr* builtin_min(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 0) return expr_new_symbol("Infinity");
    
    // Check for List arguments to flatten
    bool has_list = false;
    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && 
            arg->data.function.head->data.symbol == SYM_List) {
            has_list = true;
            break;
        }
    }
    
    if (has_list) {
        size_t new_count = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* arg = res->data.function.args[i];
            if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && 
                arg->data.function.head->data.symbol == SYM_List) {
                new_count += arg->data.function.arg_count;
            } else {
                new_count++;
            }
        }
        
        Expr** new_args = malloc(sizeof(Expr*) * new_count);
        size_t k = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* arg = res->data.function.args[i];
            if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && 
                arg->data.function.head->data.symbol == SYM_List) {
                for (size_t j = 0; j < arg->data.function.arg_count; j++) {
                    new_args[k++] = expr_copy(arg->data.function.args[j]);
                }
            } else {
                new_args[k++] = expr_copy(arg);
            }
        }
        Expr* ret = expr_new_function(expr_copy(res->data.function.head), new_args, new_count);
        free(new_args);
        return ret;
    }
    
    // Check for Overflow[] and -Infinity
    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (is_overflow(arg)) return expr_new_function(expr_new_symbol("Overflow"), NULL, 0);
        if (is_minus_infinity(arg)) return expr_copy(arg);
    }
    
    // Combine numbers and remove duplicates
    size_t unique_count = 0;
    Expr** unique_args = malloc(sizeof(Expr*) * n);
    Expr* min_num = NULL;
    bool needs_simplification = false;
    
    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (is_real_numeric(arg)) {
            if (min_num == NULL) {
                min_num = expr_copy(arg);
            } else {
                if (expr_compare(arg, min_num) < 0) {
                    expr_free(min_num);
                    min_num = expr_copy(arg);
                }
                needs_simplification = true;
            }
            continue;
        }
        
        if (is_infinity(arg)) {
            if (n > 1) needs_simplification = true;
            continue;
        }
        
        if (unique_count > 0 && expr_eq(arg, unique_args[unique_count - 1])) {
            needs_simplification = true;
            continue;
        }
        unique_args[unique_count++] = expr_copy(arg);
    }
    
    if (needs_simplification) {
        size_t final_count = unique_count + (min_num ? 1 : 0);
        if (final_count == 0) {
            if (min_num) expr_free(min_num);
            for (size_t i = 0; i < unique_count; i++) expr_free(unique_args[i]);
            free(unique_args);
            return expr_new_symbol("Infinity");
        }
        if (final_count == 1) {
            Expr* single = min_num ? min_num : unique_args[0];
            free(unique_args);
            return single;
        }
        Expr** final_args = malloc(sizeof(Expr*) * final_count);
        size_t k = 0;
        if (min_num) final_args[k++] = min_num;
        for (size_t i = 0; i < unique_count; i++) final_args[k++] = unique_args[i];

        Expr* ret = expr_new_function(expr_copy(res->data.function.head), final_args, final_count);
        free(final_args);
        free(unique_args);
        return ret;
    }
    
    if (min_num) expr_free(min_num);
    for (size_t i = 0; i < unique_count; i++) expr_free(unique_args[i]);
    free(unique_args);
    return NULL;
}

Expr* builtin_max(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 0) return make_minus_infinity();
    
    // Check for List arguments to flatten
    bool has_list = false;
    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && 
            arg->data.function.head->data.symbol == SYM_List) {
            has_list = true;
            break;
        }
    }
    
    if (has_list) {
        size_t new_count = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* arg = res->data.function.args[i];
            if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && 
                arg->data.function.head->data.symbol == SYM_List) {
                new_count += arg->data.function.arg_count;
            } else {
                new_count++;
            }
        }
        
        Expr** new_args = malloc(sizeof(Expr*) * new_count);
        size_t k = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* arg = res->data.function.args[i];
            if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && 
                arg->data.function.head->data.symbol == SYM_List) {
                for (size_t j = 0; j < arg->data.function.arg_count; j++) {
                    new_args[k++] = expr_copy(arg->data.function.args[j]);
                }
            } else {
                new_args[k++] = expr_copy(arg);
            }
        }
        Expr* ret = expr_new_function(expr_copy(res->data.function.head), new_args, new_count);
        free(new_args);
        return ret;
    }
    
    // Check for Overflow[] and Infinity
    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (is_overflow(arg)) return expr_new_function(expr_new_symbol("Overflow"), NULL, 0);
        if (is_infinity(arg)) return expr_copy(arg);
    }
    
    // Combine numbers and remove duplicates
    size_t unique_count = 0;
    Expr** unique_args = malloc(sizeof(Expr*) * n);
    Expr* max_num = NULL;
    bool needs_simplification = false;
    
    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (is_real_numeric(arg)) {
            if (max_num == NULL) {
                max_num = expr_copy(arg);
            } else {
                if (expr_compare(arg, max_num) > 0) {
                    expr_free(max_num);
                    max_num = expr_copy(arg);
                }
                needs_simplification = true;
            }
            continue;
        }
        
        if (is_minus_infinity(arg)) {
            if (n > 1) needs_simplification = true;
            continue;
        }
        
        if (unique_count > 0 && expr_eq(arg, unique_args[unique_count - 1])) {
            needs_simplification = true;
            continue;
        }
        unique_args[unique_count++] = expr_copy(arg);
    }
    
    if (needs_simplification) {
        size_t final_count = unique_count + (max_num ? 1 : 0);
        if (final_count == 0) {
            if (max_num) expr_free(max_num);
            for (size_t i = 0; i < unique_count; i++) expr_free(unique_args[i]);
            free(unique_args);
            return make_minus_infinity();
        }
        if (final_count == 1) {
            Expr* single = max_num ? max_num : unique_args[0];
            free(unique_args);
            return single;
        }
        Expr** final_args = malloc(sizeof(Expr*) * final_count);
        size_t k = 0;
        if (max_num) final_args[k++] = max_num;
        for (size_t i = 0; i < unique_count; i++) final_args[k++] = unique_args[i];

        Expr* ret = expr_new_function(expr_copy(res->data.function.head), final_args, final_count);
        free(final_args);
        free(unique_args);
        return ret;
    }
    
    if (max_num) expr_free(max_num);
    for (size_t i = 0; i < unique_count; i++) expr_free(unique_args[i]);
    free(unique_args);
    return NULL;
}

Expr* builtin_range(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 3) return NULL;
    
    size_t len = res->data.function.arg_count;
    Expr* imin_e = NULL;
    Expr* imax_e = NULL;
    Expr* di_e = NULL;
    
    if (len == 1) {
        imin_e = expr_new_integer(1);
        imax_e = expr_copy(res->data.function.args[0]);
        di_e = expr_new_integer(1);
    } else if (len == 2) {
        imin_e = expr_copy(res->data.function.args[0]);
        imax_e = expr_copy(res->data.function.args[1]);
        di_e = expr_new_integer(1);
    } else if (len == 3) {
        imin_e = expr_copy(res->data.function.args[0]);
        imax_e = expr_copy(res->data.function.args[1]);
        di_e = expr_copy(res->data.function.args[2]);
    }
    
    bool is_real = false;
    double min_val = 0, max_val = 0, di_val = 0;
    int64_t n, d;
    
    if (imin_e->type == EXPR_REAL || imax_e->type == EXPR_REAL || di_e->type == EXPR_REAL) is_real = true;
    
    if (imin_e->type == EXPR_INTEGER) min_val = (double)imin_e->data.integer;
    else if (imin_e->type == EXPR_REAL) min_val = imin_e->data.real;
    else if (is_rational(imin_e, &n, &d)) min_val = (double)n / d;
    else goto L_fail_range;
    
    if (imax_e->type == EXPR_INTEGER) max_val = (double)imax_e->data.integer;
    else if (imax_e->type == EXPR_REAL) max_val = imax_e->data.real;
    else if (is_rational(imax_e, &n, &d)) max_val = (double)n / d;
    else goto L_fail_range;
    
    if (di_e->type == EXPR_INTEGER) di_val = (double)di_e->data.integer;
    else if (di_e->type == EXPR_REAL) di_val = di_e->data.real;
    else if (is_rational(di_e, &n, &d)) di_val = (double)n / d;
    else goto L_fail_range;
    
    if (di_val == 0) goto L_fail_range;
    if ((di_val > 0 && min_val > max_val) || (di_val < 0 && min_val < max_val)) {
        expr_free(imin_e); expr_free(imax_e); expr_free(di_e);
        return expr_new_function(expr_new_symbol("List"), NULL, 0);
    }
    
    size_t results_cap = 16;
    size_t results_count = 0;
    Expr** results = malloc(sizeof(Expr*) * results_cap);
    
    double val = min_val;
    int steps = 0;
    Expr* curr_e = expr_copy(imin_e);
    
    while ((di_val > 0 && val <= max_val + 1e-14) || (di_val < 0 && val >= max_val - 1e-14)) {
        Expr* i_val = is_real ? expr_new_real(val) : expr_copy(curr_e);
        
        if (results_count == results_cap) { results_cap *= 2; results = realloc(results, sizeof(Expr*) * results_cap); }
        results[results_count++] = i_val;
        
        Expr* next_args[2] = { expr_copy(curr_e), expr_copy(di_e) };
        Expr* next_expr = expr_new_function(expr_new_symbol("Plus"), next_args, 2);
        Expr* next_e = evaluate(next_expr);
        expr_free(next_expr);
        expr_free(curr_e);
        curr_e = next_e;
        
        val += di_val;
        steps++;
        if (steps > 1000000) break; 
    }
    
    if (curr_e) expr_free(curr_e);
    expr_free(imin_e);
    expr_free(imax_e);
    expr_free(di_e);
    
    Expr* result_list = expr_new_function(expr_new_symbol("List"), results, results_count);
    free(results);
    return result_list;

L_fail_range:
    if (imin_e) expr_free(imin_e);
    if (imax_e) expr_free(imax_e);
    if (di_e) expr_free(di_e);
    return NULL;
}

/* ------------------- Join ------------------- */

/*
 * join_at_level: recursively join expressions at a given depth.
 * level 1 = concatenate the top-level arguments (basic join).
 * level 2 = descend one level into each list, joining corresponding
 *           sub-elements across the input lists. Ragged arrays are
 *           handled by concatenating successive elements at that level.
 *
 * lists: array of n_lists expressions to join
 * n_lists: number of lists
 * level: remaining depth to descend (>= 1)
 *
 * Returns a new Expr* on success, or NULL if types don't match.
 */
static Expr* join_at_level(Expr** lists, size_t n_lists, int level) {
    if (n_lists == 0) return NULL;

    Expr* first = lists[0];
    if (first->type != EXPR_FUNCTION) return NULL;

    Expr* head = first->data.function.head;

    /* Verify all lists share the same head */
    for (size_t i = 1; i < n_lists; i++) {
        if (lists[i]->type != EXPR_FUNCTION ||
            !expr_eq(lists[i]->data.function.head, head))
            return NULL;
    }

    if (level == 1) {
        /* Base case: concatenate all arguments */
        size_t total = 0;
        for (size_t i = 0; i < n_lists; i++)
            total += lists[i]->data.function.arg_count;

        Expr** new_args = malloc(sizeof(Expr*) * (total > 0 ? total : 1));
        size_t curr = 0;
        for (size_t i = 0; i < n_lists; i++) {
            Expr* li = lists[i];
            for (size_t j = 0; j < li->data.function.arg_count; j++)
                new_args[curr++] = expr_copy(li->data.function.args[j]);
        }
        Expr* result = expr_new_function(expr_copy(head), new_args, total);
        free(new_args);
        return result;
    }

    /* Recursive case: level > 1.
     * Find the maximum number of sub-elements across all lists.
     * For each position k, gather the sub-elements from each list
     * that have a k-th element and recursively join them at level-1.
     * Lists that don't have a k-th element simply contribute nothing. */
    size_t max_len = 0;
    for (size_t i = 0; i < n_lists; i++) {
        if (lists[i]->data.function.arg_count > max_len)
            max_len = lists[i]->data.function.arg_count;
    }

    Expr** result_args = malloc(sizeof(Expr*) * (max_len > 0 ? max_len : 1));
    size_t result_count = 0;

    /* Temporary buffer for gathering sub-elements at position k */
    Expr** subs = malloc(sizeof(Expr*) * n_lists);

    for (size_t k = 0; k < max_len; k++) {
        size_t n_subs = 0;
        for (size_t i = 0; i < n_lists; i++) {
            if (k < lists[i]->data.function.arg_count)
                subs[n_subs++] = lists[i]->data.function.args[k];
        }

        if (n_subs == 0) continue;

        if (n_subs == 1) {
            /* Only one list contributes at this position: copy as-is */
            result_args[result_count++] = expr_copy(subs[0]);
        } else {
            Expr* joined = join_at_level(subs, n_subs, level - 1);
            if (!joined) {
                /* Cleanup on failure */
                for (size_t j = 0; j < result_count; j++)
                    expr_free(result_args[j]);
                free(result_args);
                free(subs);
                return NULL;
            }
            result_args[result_count++] = joined;
        }
    }

    free(subs);
    Expr* result = expr_new_function(expr_copy(head), result_args, result_count);
    free(result_args);
    return result;
}

Expr* builtin_join(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1)
        return NULL;

    size_t n_args = res->data.function.arg_count;
    int level = 1;
    size_t n_lists = n_args;

    /* Check if the last argument is an integer (level specification) */
    Expr* last = res->data.function.args[n_args - 1];
    if (n_args >= 2 && last->type == EXPR_INTEGER) {
        level = (int)last->data.integer;
        if (level < 1) return NULL;
        n_lists = n_args - 1;
    }

    if (n_lists < 1) return NULL;

    Expr* result = join_at_level(res->data.function.args, n_lists, level);
    if (!result) return NULL;

    return result;
}
