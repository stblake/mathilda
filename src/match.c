#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "match.h"
#include <gmp.h>
#include "part.h"
#include "eval.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include "default_helper.h"
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static bool is_atom(const Expr* e) {
    if (!e) return true;
    if (e->type != EXPR_FUNCTION) return true;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* head_name = e->data.function.head->data.symbol;
        if (head_name == SYM_Complex || head_name == SYM_Rational) {
            return true;
        }
    }
    return false;
}

MatchEnv* env_new(void) {
    MatchEnv* env = malloc(sizeof(MatchEnv));
    env->count = 0;
    env->capacity = 8;
    env->callback = NULL;
    env->callback_data = NULL;
    env->symbols = malloc(sizeof(char*) * env->capacity);
    env->values = malloc(sizeof(Expr*) * env->capacity);
    return env;
}

void env_free(MatchEnv* env) {
    if (!env) return;
    for (size_t i = 0; i < env->count; i++) {
        free(env->symbols[i]);
        expr_free(env->values[i]);
    }
    free(env->symbols);
    free(env->values);
    free(env);
}

void env_set(MatchEnv* env, const char* symbol, Expr* value) {
    for (size_t i = 0; i < env->count; i++) {
        if (strcmp(env->symbols[i], symbol) == 0) {
            expr_free(env->values[i]);
            env->values[i] = expr_copy(value);
            return;
        }
    }
    if (env->count >= env->capacity) {
        env->capacity *= 2;
        env->symbols = realloc(env->symbols, sizeof(char*) * env->capacity);
        env->values = realloc(env->values, sizeof(Expr*) * env->capacity);
    }
    env->symbols[env->count] = strdup(symbol);
    env->values[env->count] = expr_copy(value);
    env->count++;
}

static void env_rollback(MatchEnv* env, size_t saved_count) {
    while (env->count > saved_count) {
        env->count--;
        free(env->symbols[env->count]);
        expr_free(env->values[env->count]);
    }
}

Expr* env_get(MatchEnv* env, const char* symbol) {
    for (size_t i = 0; i < env->count; i++) {
        if (strcmp(env->symbols[i], symbol) == 0) {
            return env->values[i];
        }
    }
    return NULL;
}

static bool is_pattern(Expr* e, Expr** sym_out, Expr** pat_out) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        if (e->data.function.head->data.symbol == SYM_Pattern) {
            if (e->data.function.arg_count == 2) {
                if (sym_out) *sym_out = e->data.function.args[0];
                if (pat_out) *pat_out = e->data.function.args[1];
                return true;
            }
        }
    }
    return false;
}

static bool is_blank(Expr* e, Expr** head_out) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        if (e->data.function.head->data.symbol == SYM_Blank) {
            if (e->data.function.arg_count == 0) {
                if (head_out) *head_out = NULL;
                return true;
            } else if (e->data.function.arg_count == 1) {
                if (head_out) *head_out = e->data.function.args[0];
                return true;
            }
        }
    }
    return false;
}

static bool is_sequence_blank(Expr* e, Expr** head_out, int* min_len) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        const char* head = e->data.function.head->data.symbol;
        if (head == SYM_BlankSequence) {
            if (min_len) *min_len = 1;
            if (head_out) *head_out = (e->data.function.arg_count == 1) ? e->data.function.args[0] : NULL;
            return true;
        }
        if (head == SYM_BlankNullSequence) {
            if (min_len) *min_len = 0;
            if (head_out) *head_out = (e->data.function.arg_count == 1) ? e->data.function.args[0] : NULL;
            return true;
        }
    }
    return false;
}

static bool is_repeated(Expr* e, Expr** rep_pat, int* min_len, int* max_len) {
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* head = e->data.function.head->data.symbol;
    bool is_rep = (head == SYM_Repeated);
    bool is_rep_null = (head == SYM_RepeatedNull);
    if (!is_rep && !is_rep_null) return false;

    *min_len = is_rep ? 1 : 0;
    *max_len = -1; // -1 means infinity
    
    if (e->data.function.arg_count >= 1) {
        *rep_pat = e->data.function.args[0];
    } else {
        return false;
    }
    
    if (e->data.function.arg_count >= 2) {
        Expr* spec = e->data.function.args[1];
        if (spec->type == EXPR_INTEGER) {
            *max_len = (int)spec->data.integer;
        } else if (spec->type == EXPR_FUNCTION && spec->data.function.head->type == EXPR_SYMBOL && spec->data.function.head->data.symbol == SYM_List) {
            if (spec->data.function.arg_count == 1 && spec->data.function.args[0]->type == EXPR_INTEGER) {
                *min_len = (int)spec->data.function.args[0]->data.integer;
                *max_len = *min_len;
            } else if (spec->data.function.arg_count == 2) {
                if (spec->data.function.args[0]->type == EXPR_INTEGER) {
                    *min_len = (int)spec->data.function.args[0]->data.integer;
                }
                if (spec->data.function.args[1]->type == EXPR_INTEGER) {
                    *max_len = (int)spec->data.function.args[1]->data.integer;
                } else if (spec->data.function.args[1]->type == EXPR_SYMBOL && spec->data.function.args[1]->data.symbol == SYM_Infinity) {
                    *max_len = -1;
                }
            }
        }
    }
    return true;
}

typedef struct ParentMatch {
    Expr** exprs;
    size_t n_exprs;
    Expr** pats;
    size_t n_pats;
    Expr* condition;
    Expr* pat_head;
    size_t total_pats;
    struct ParentMatch* parent;
} ParentMatch;

static bool match_internal(Expr* expr, Expr* pattern, MatchEnv* env, ParentMatch* parent);
static bool match_args_internal(Expr** exprs, size_t n_exprs, Expr** pats, size_t n_pats, MatchEnv* env, Expr* condition, Expr* pat_head, size_t total_pats, ParentMatch* parent);

#include "part.h" // for expr_head

/* Helper to get head without allocating */
static Expr* get_expr_head_borrowed(Expr* e) {
    if (is_atom(e)) return NULL; // atoms return symbol head usually, but match handles it
    if (e->type == EXPR_FUNCTION) return e->data.function.head;
    return NULL;
}

static bool next_combination(int* comb, int n, int k) {
    if (k == 0) return false;
    int i = k - 1;
    while (i >= 0 && comb[i] == n - k + i) {
        i--;
    }
    if (i < 0) return false;
    comb[i]++;
    for (int j = i + 1; j < k; j++) {
        comb[j] = comb[j - 1] + 1;
    }
    return true;
}

static void extract_subset(Expr** exprs, size_t n_exprs, int* comb, int k, Expr** subset, Expr** remainder) {
    int c_idx = 0;
    int r_idx = 0;
    for (int i = 0; i < (int)n_exprs; i++) {
        if (c_idx < k && comb[c_idx] == i) {
            subset[c_idx++] = exprs[i];
        } else {
            remainder[r_idx++] = exprs[i];
        }
    }
}

// Call parent if exists, otherwise return true (or check top condition)

static ParentMatch top_level_sentinel = {0};

static bool call_parent(MatchEnv* env, ParentMatch* parent) {
    if (parent == &top_level_sentinel) {
        if (env->callback) {
            return env->callback(env, env->callback_data);
        }
        return true;
    }
    if (parent) {
        return match_args_internal(parent->exprs, parent->n_exprs, parent->pats, parent->n_pats, env, parent->condition, parent->pat_head, parent->total_pats, parent->parent);
    }
    return true;
}


static bool match_internal(Expr* expr, Expr* pattern, MatchEnv* env, ParentMatch* parent) {
    if (!pattern) return false;
    if (!expr) return false;

    /* HoldPattern[p] is equivalent to p for matching purposes. The
     * HoldAll attribute on HoldPattern prevents the inner expression from
     * being evaluated when the rule LHS is formed, so patterns that would
     * otherwise simplify (e.g. _+_ becoming 2 Blank[]) are preserved. */
    if (pattern->type == EXPR_FUNCTION && pattern->data.function.head->type == EXPR_SYMBOL &&
        pattern->data.function.head->data.symbol == SYM_HoldPattern &&
        pattern->data.function.arg_count == 1) {
        return match_internal(expr, pattern->data.function.args[0], env, parent);
    }

    // Handle Except
    if (pattern->type == EXPR_FUNCTION && pattern->data.function.head->type == EXPR_SYMBOL &&
        pattern->data.function.head->data.symbol == SYM_Except) {
        if (pattern->data.function.arg_count == 1) {
            size_t saved_env_count = env->count;
            if (match_internal(expr, pattern->data.function.args[0], env, NULL)) {
                env_rollback(env, saved_env_count);
                return false;
            }
            return call_parent(env, parent);
        }
        if (pattern->data.function.arg_count == 2) {
            size_t saved_env_count = env->count;
            if (!match_internal(expr, pattern->data.function.args[1], env, NULL)) {
                env_rollback(env, saved_env_count);
                return false;
            }
            size_t inner_saved = env->count;
            if (match_internal(expr, pattern->data.function.args[0], env, NULL)) {
                env_rollback(env, saved_env_count);
                return false;
            }
            env_rollback(env, inner_saved);
            return call_parent(env, parent);
        }
    }

    // Handle Condition
    if (pattern->type == EXPR_FUNCTION && pattern->data.function.head->type == EXPR_SYMBOL &&
        pattern->data.function.head->data.symbol == SYM_Condition) {
        if (pattern->data.function.arg_count != 2) return false;
        
        Expr* inner_pat = pattern->data.function.args[0];
        Expr* cond = pattern->data.function.args[1];
        
        if (inner_pat->type == EXPR_FUNCTION && expr->type == EXPR_FUNCTION) {
            size_t saved_env_count = env->count;
            if (!match_internal(expr->data.function.head, inner_pat->data.function.head, env, NULL)) {
                env_rollback(env, saved_env_count);
                return false;
            }
            bool res = match_args_internal(expr->data.function.args, expr->data.function.arg_count,
                                 inner_pat->data.function.args, inner_pat->data.function.arg_count, env, cond, inner_pat->data.function.head, inner_pat->data.function.arg_count, parent);
            if (!res) env_rollback(env, saved_env_count);
            return res;
        }

        size_t saved_env_count = env->count;
        if (!match_internal(expr, inner_pat, env, NULL)) {
            env_rollback(env, saved_env_count);
            return false;
        }
        Expr* expanded_test = replace_bindings(cond, env);
        Expr* result = evaluate(expanded_test);
        expr_free(expanded_test);
        bool success = (result->type == EXPR_SYMBOL && result->data.symbol == SYM_True);
        expr_free(result);
        if (!success) {
            env_rollback(env, saved_env_count);
            return false;
        }
        return call_parent(env, parent);
    }

    Expr* b_head = NULL;
    if (is_blank(pattern, &b_head)) {
        bool head_ok = false;
        if (!b_head) head_ok = true;
        else {
            Expr* h = get_expr_head_borrowed(expr);
            if (h) head_ok = expr_eq(h, b_head);
            else if (b_head->type == EXPR_SYMBOL) {
                const char* hn = b_head->data.symbol;
                if (expr->type == EXPR_INTEGER && hn == SYM_Integer) head_ok = true;
                else if (expr->type == EXPR_REAL && hn == SYM_Real) head_ok = true;
                else if (expr->type == EXPR_SYMBOL && hn == SYM_Symbol) head_ok = true;
                else if (expr->type == EXPR_STRING && hn == SYM_String) head_ok = true;
            }
        }
        if (head_ok) return call_parent(env, parent);
        return false;
    }

    int min_len = 0;
    if (is_sequence_blank(pattern, &b_head, &min_len)) {
        if (min_len > 1) return false;
        bool head_ok = false;
        if (!b_head) head_ok = true;
        else {
            Expr* h = get_expr_head_borrowed(expr);
            if (h) head_ok = expr_eq(h, b_head);
            else if (b_head->type == EXPR_SYMBOL) {
                const char* hn = b_head->data.symbol;
                if (expr->type == EXPR_INTEGER && hn == SYM_Integer) head_ok = true;
                else if (expr->type == EXPR_REAL && hn == SYM_Real) head_ok = true;
                else if (expr->type == EXPR_SYMBOL && hn == SYM_Symbol) head_ok = true;
                else if (expr->type == EXPR_STRING && hn == SYM_String) head_ok = true;
            }
        }
        if (head_ok) return call_parent(env, parent);
        return false;
    }

    Expr* p_sym = NULL;
    Expr* p_pat = NULL;
    if (is_pattern(pattern, &p_sym, &p_pat)) {
        if (p_sym->type == EXPR_SYMBOL) {
            size_t saved_env_count = env->count;
            // Match p_pat first
            if (!match_internal(expr, p_pat, env, NULL)) {
                return false;
            }
            // Then handle binding and parent
            Expr* existing = env_get(env, p_sym->data.symbol);
            if (existing) {
                if (!expr_eq(expr, existing)) {
                    env_rollback(env, saved_env_count);
                    return false;
                }
                // Match with parent
                if (call_parent(env, parent)) return true;
                env_rollback(env, saved_env_count);
                return false;
            } else {
                env_set(env, p_sym->data.symbol, expr);
                if (call_parent(env, parent)) return true;
                env_rollback(env, saved_env_count);
                return false;
            }
        }
    }

    // Handle Alternatives
    if (pattern->type == EXPR_FUNCTION && pattern->data.function.head->type == EXPR_SYMBOL &&
        pattern->data.function.head->data.symbol == SYM_Alternatives) {
        size_t saved_env_count = env->count;
        for (size_t i = 0; i < pattern->data.function.arg_count; i++) {
            if (match_internal(expr, pattern->data.function.args[i], env, parent)) return true;
            env_rollback(env, saved_env_count);
        }
        return false;
    }

    // Handle PatternTest
    if (pattern->type == EXPR_FUNCTION && pattern->data.function.head->type == EXPR_SYMBOL &&
        pattern->data.function.head->data.symbol == SYM_PatternTest) {
        if (pattern->data.function.arg_count != 2) return false;
        size_t saved_env_count = env->count;
        if (!match_internal(expr, pattern->data.function.args[0], env, NULL)) {
            return false;
        }
        
        Expr* test_func = pattern->data.function.args[1];
        Expr* call_args[1] = { expr_copy(expr) };
        Expr* test_call = expr_new_function(expr_copy(test_func), call_args, 1);
        Expr* result = evaluate(test_call);
        expr_free(test_call);
        
        bool success = (result->type == EXPR_SYMBOL && result->data.symbol == SYM_True);
        expr_free(result);
        if (!success) {
            env_rollback(env, saved_env_count);
            return false;
        }
        if (call_parent(env, parent)) return true;
        env_rollback(env, saved_env_count);
        return false;
    }

    if (expr->type == pattern->type) {
        switch (expr->type) {
            case EXPR_INTEGER:
                if (expr->data.integer == pattern->data.integer) return call_parent(env, parent);
                return false;
            case EXPR_REAL:
                if (expr->data.real == pattern->data.real) return call_parent(env, parent);
                return false;
            case EXPR_SYMBOL: {
                Expr* bound = env_get(env, pattern->data.symbol);
                if (bound) {
                    if (expr_eq(expr, bound)) return call_parent(env, parent);
                    return false;
                }
                if (strcmp(expr->data.symbol, pattern->data.symbol) == 0) return call_parent(env, parent);
                return false;
            }
            case EXPR_STRING:
                if (strcmp(expr->data.string, pattern->data.string) == 0) return call_parent(env, parent);
                return false;
            case EXPR_BIGINT:
                if (mpz_cmp(expr->data.bigint, pattern->data.bigint) == 0) return call_parent(env, parent);
                return false;
#ifdef USE_MPFR
            case EXPR_MPFR:
                /* Match identical-precision MPFR values bit-for-bit. */
                if (mpfr_get_prec(expr->data.mpfr) == mpfr_get_prec(pattern->data.mpfr)
                    && mpfr_equal_p(expr->data.mpfr, pattern->data.mpfr)) {
                    return call_parent(env, parent);
                }
                return false;
#endif
            case EXPR_FUNCTION: {
                size_t saved_env_count = env->count;
                if (match_internal(expr->data.function.head, pattern->data.function.head, env, NULL)) {
                    if (match_args_internal(expr->data.function.args, expr->data.function.arg_count,
                                   pattern->data.function.args, pattern->data.function.arg_count, env, NULL, pattern->data.function.head, pattern->data.function.arg_count, parent)) {
                        return true;
                    }
                }
                env_rollback(env, saved_env_count);
                break;
            }
        }
    }

    // Try OneIdentity if pattern is a function
    if (pattern->type == EXPR_FUNCTION && pattern->data.function.head->type == EXPR_SYMBOL) {
        SymbolDef* def = symtab_get_def(pattern->data.function.head->data.symbol);
        if (def && (def->attributes & ATTR_ONEIDENTITY)) {
            size_t saved_env_count = env->count;
            Expr* args[1] = { expr };
            if (match_args_internal(args, 1, pattern->data.function.args, pattern->data.function.arg_count, env, NULL, pattern->data.function.head, pattern->data.function.arg_count, parent)) {
                return true;
            }
            env_rollback(env, saved_env_count);
        }
    }

    return false;
}

bool match(Expr* expr, Expr* pattern, MatchEnv* env) {
    return match_internal(expr, pattern, env, &top_level_sentinel);
}

static bool match_args_internal(Expr** exprs, size_t n_exprs, Expr** pats, size_t n_pats, MatchEnv* env, Expr* condition, Expr* pat_head, size_t total_pats, ParentMatch* parent) {
    
    if (n_pats == 0) {
        if (n_exprs != 0) return false;
        
        if (condition) {
            Expr* expanded_test = replace_bindings(condition, env);
            Expr* result = evaluate(expanded_test);
            expr_free(expanded_test);
            bool success = (result->type == EXPR_SYMBOL && result->data.symbol == SYM_True);
            expr_free(result);
            if (!success) return false;
        }
        
        return call_parent(env, parent);
    }

    Expr* p = pats[0];
    
    bool is_optional = false;
    bool is_shortest = false;
    bool is_longest = false;
    Expr* opt_container = p;
    Expr* bind_sym = NULL;
    Expr* current_p = p;

    while (current_p->type == EXPR_FUNCTION && current_p->data.function.head->type == EXPR_SYMBOL && current_p->data.function.arg_count >= 1) {
        const char* head = current_p->data.function.head->data.symbol;
        if (head == SYM_Pattern && current_p->data.function.arg_count == 2) {
            bind_sym = current_p->data.function.args[0];
            current_p = current_p->data.function.args[1];
        } else if (head == SYM_Shortest) {
            is_shortest = true;
            current_p = current_p->data.function.args[0];
        } else if (head == SYM_Longest) {
            is_longest = true;
            is_shortest = false;
            current_p = current_p->data.function.args[0];
        } else if (head == SYM_Optional) {
            is_optional = true;
            opt_container = current_p;
            current_p = current_p->data.function.args[0];
        } else {
            break;
        }
    }
    
    Expr* opt_pat = current_p;
    if (is_optional && is_shortest) {
        size_t saved_env_count = env->count;
        Expr* def_val = NULL;
        if (opt_container->data.function.arg_count == 2) {
            def_val = expr_copy(opt_container->data.function.args[1]);
        } else {
            def_val = get_default_value(pat_head, total_pats - n_pats + 1, total_pats);
        }

        if (def_val) {
            Expr* p_sym = bind_sym;
            Expr* inner_p = opt_pat;
            if (!p_sym) is_pattern(opt_pat, &p_sym, &inner_p);
            bool consistent = true;
            if (p_sym && p_sym->type == EXPR_SYMBOL) {
                Expr* existing = env_get(env, p_sym->data.symbol);
                if (existing) {
                    if (!expr_eq(existing, def_val)) consistent = false;
                } else {
                    env_set(env, p_sym->data.symbol, def_val);
                }
            }
            expr_free(def_val);


            if (consistent && match_args_internal(exprs, n_exprs, pats + 1, n_pats - 1, env, condition, pat_head, total_pats, parent)) {
                return true;
            }
            env_rollback(env, saved_env_count);
        }
    }

    Expr* inner_p = opt_pat;
    Expr* p_sym = bind_sym;
    int min_len = 0;
    int max_len = -1;
    Expr* b_head = NULL;
    Expr* rep_pat = NULL;

    bool is_seq = false;
    bool is_rep = false;
    if (!p_sym && is_pattern(opt_pat, &p_sym, &inner_p)) {
        is_seq = is_sequence_blank(inner_p, &b_head, &min_len);
        if (!is_seq) is_rep = is_repeated(inner_p, &rep_pat, &min_len, &max_len);
    } else {
        is_seq = is_sequence_blank(inner_p, &b_head, &min_len);
        if (!is_seq) is_rep = is_repeated(inner_p, &rep_pat, &min_len, &max_len);
    }


    bool is_orderless = false;
    bool is_flat = false;
    if (pat_head && pat_head->type == EXPR_SYMBOL) {
        SymbolDef* def = symtab_get_def(pat_head->data.symbol);
        if (def) {
            if (def->attributes & ATTR_FLAT) is_flat = true;
            if (def->attributes & ATTR_ORDERLESS) is_orderless = true;
        }
    }

    size_t min_k, max_k;
    if (is_seq) {
        min_k = min_len;
        max_k = n_exprs;
    } else if (is_rep) {
        min_k = min_len;
        max_k = (max_len == -1 || max_len > (int)n_exprs) ? n_exprs : (size_t)max_len;
    } else if (is_flat) {
        min_k = 1;
        max_k = n_exprs;
    } else {
        min_k = 1;
        max_k = (n_exprs > 0) ? 1 : 0;
    }

    if (min_k <= max_k) {
        for (size_t step = 0; step <= max_k - min_k; step++) {
            size_t k = is_longest ? (max_k - step) : (min_k + step);
            if (k > n_exprs) continue;
            
            int* comb = NULL;
            if (k > 0) {
                comb = malloc(k * sizeof(int));
                for (int i = 0; i < (int)k; i++) comb[i] = i;
            }

        do {
            Expr** subset = NULL;
            Expr** remainder = NULL;
            if (k > 0) subset = malloc(k * sizeof(Expr*));
            if (n_exprs - k > 0) remainder = malloc((n_exprs - k) * sizeof(Expr*));
            
            extract_subset(exprs, n_exprs, comb, (int)k, subset, remainder);

            size_t saved_env = env->count;

            if (is_seq || is_rep) {
                bool type_ok = true;
                if (is_seq) {
                    if (b_head) {
                        for (size_t i = 0; i < k; i++) {
                            Expr* h = get_expr_head_borrowed(subset[i]);
                            bool ok = false;
                            if (h) ok = expr_eq(h, b_head);
                            else if (b_head->type == EXPR_SYMBOL) {
                                const char* hn = b_head->data.symbol;
                                if (subset[i]->type == EXPR_INTEGER && hn == SYM_Integer) ok = true;
                                else if (subset[i]->type == EXPR_REAL && hn == SYM_Real) ok = true;
                                else if (subset[i]->type == EXPR_SYMBOL && hn == SYM_Symbol) ok = true;
                                else if (subset[i]->type == EXPR_STRING && hn == SYM_String) ok = true;
                            }
                            if (!ok) { type_ok = false; break; }
                        }
                    }
                } else if (is_rep) {
                    for (size_t i = 0; i < k; i++) {
                        if (!match_internal(subset[i], rep_pat, env, NULL)) { type_ok = false; break; }
                    }
                }

                if (type_ok) {
                    if (p_sym) {
                        Expr* seq_val = expr_new_function(expr_new_symbol("Sequence"), NULL, k);
                        for (size_t i = 0; i < k; i++) seq_val->data.function.args[i] = expr_copy(subset[i]);
                        Expr* existing = env_get(env, p_sym->data.symbol);
                        if (existing) {
                            if (expr_eq(seq_val, existing)) {
                                if (match_args_internal(remainder, n_exprs - k, pats + 1, n_pats - 1, env, condition, pat_head, total_pats, parent)) {
                                    expr_free(seq_val); { if (subset) free(subset); if (remainder) free(remainder); if (comb) free(comb); return true; }
                                }
                            }
                        } else {
                            env_set(env, p_sym->data.symbol, seq_val);
                            if (match_args_internal(remainder, n_exprs - k, pats + 1, n_pats - 1, env, condition, pat_head, total_pats, parent)) {
                                expr_free(seq_val); { if (subset) free(subset); if (remainder) free(remainder); if (comb) free(comb); return true; }
                            }
                        }
                        expr_free(seq_val);
                    } else {
                        if (match_args_internal(remainder, n_exprs - k, pats + 1, n_pats - 1, env, condition, pat_head, total_pats, parent)) {
                            { if (subset) free(subset); if (remainder) free(remainder); if (comb) free(comb); return true; }
                        }
                    }
                }
            } else {
                Expr* matched_val = NULL;
                if (k == 1) {
                    matched_val = expr_copy(subset[0]);
                } else if (k > 1) {
                    matched_val = expr_new_function(expr_copy(pat_head), NULL, k);
                    for (size_t i = 0; i < k; i++) {
                        matched_val->data.function.args[i] = expr_copy(subset[i]);
                    }
                }

                if (matched_val) {
                    size_t saved_env_inner = env->count;
                    if (match_internal(matched_val, inner_p, env, NULL)) {
                        if (p_sym) {
                            Expr* existing = env_get(env, p_sym->data.symbol);
                            if (existing) {
                                if (expr_eq(matched_val, existing)) {
                                    if (match_args_internal(remainder, n_exprs - k, pats + 1, n_pats - 1, env, condition, pat_head, total_pats, parent)) {
                                        expr_free(matched_val); { if (subset) free(subset); if (remainder) free(remainder); if (comb) free(comb); return true; }
                                    }
                                }
                            } else {
                                env_set(env, p_sym->data.symbol, matched_val);
                                if (match_args_internal(remainder, n_exprs - k, pats + 1, n_pats - 1, env, condition, pat_head, total_pats, parent)) {
                                    expr_free(matched_val); { if (subset) free(subset); if (remainder) free(remainder); if (comb) free(comb); return true; }
                                }
                            }
                        } else {
                            if (match_args_internal(remainder, n_exprs - k, pats + 1, n_pats - 1, env, condition, pat_head, total_pats, parent)) {
                                expr_free(matched_val); { if (subset) free(subset); if (remainder) free(remainder); if (comb) free(comb); return true; }
                            }
                        }
                    }
                    env_rollback(env, saved_env_inner);
                    expr_free(matched_val);
                }

            }

            env_rollback(env, saved_env);
            if (subset) free(subset);
            if (remainder) free(remainder);

        } while (is_orderless && next_combination(comb, n_exprs, (int)k));
        
        if (comb) free(comb);
        }
    }

    if (is_optional && !is_shortest) {
        size_t saved_env_count = env->count;
        Expr* def_val = NULL;
        if (opt_container->data.function.arg_count == 2) {
            def_val = expr_copy(opt_container->data.function.args[1]);
        } else {
            def_val = get_default_value(pat_head, total_pats - n_pats + 1, total_pats);
        }

        if (def_val) {
            Expr* p_sym = bind_sym;
            Expr* inner_p = opt_pat;
            if (!p_sym) is_pattern(opt_pat, &p_sym, &inner_p);
            bool consistent = true;
            if (p_sym && p_sym->type == EXPR_SYMBOL) {
                Expr* existing = env_get(env, p_sym->data.symbol);
                if (existing) {
                    /* Don't overwrite an existing binding for this name; the
                     * Optional fallback can only commit to the default value if
                     * it is consistent with what is already bound. Otherwise the
                     * pattern would produce inconsistent bindings such as
                     * `a_. Sin[x_]^2 + a_. Cos[x_]^2` matching `Sin[x]^2 -
                     * Cos[x]^2` with one a=1 and another a=-1. */
                    if (!expr_eq(existing, def_val)) consistent = false;
                } else {
                    env_set(env, p_sym->data.symbol, def_val);
                }
            }
            expr_free(def_val);


            if (consistent && match_args_internal(exprs, n_exprs, pats + 1, n_pats - 1, env, condition, pat_head, total_pats, parent)) {
                return true;
            }
            env_rollback(env, saved_env_count);
        }
    }

    return false;
}

Expr* replace_bindings(Expr* expr, MatchEnv* env) {
    if (!expr) return NULL;
    
    if (expr->type == EXPR_SYMBOL) {
        Expr* bound = env_get(env, expr->data.symbol);
        if (bound) {
            return expr_copy(bound);
        }
    }
    
    if (expr->type == EXPR_FUNCTION) {
        Expr* new_head = replace_bindings(expr->data.function.head, env);
        
        bool skip_flattening = false;
        if (new_head->type == EXPR_SYMBOL && (
            new_head->data.symbol == SYM_Set ||
            new_head->data.symbol == SYM_SetDelayed ||
            new_head->data.symbol == SYM_Rule ||
            new_head->data.symbol == SYM_RuleDelayed)) {
            skip_flattening = true;
        }

        // Temporarily store replaced args to check for Sequence expansion
        Expr** temp_args = malloc(sizeof(Expr*) * expr->data.function.arg_count);
        size_t new_count = 0;
        for (size_t i = 0; i < expr->data.function.arg_count; i++) {
            temp_args[i] = replace_bindings(expr->data.function.args[i], env);
            if (!skip_flattening && temp_args[i]->type == EXPR_FUNCTION && 
                temp_args[i]->data.function.head->type == EXPR_SYMBOL &&
                temp_args[i]->data.function.head->data.symbol == SYM_Sequence) {
                new_count += temp_args[i]->data.function.arg_count;
            } else {
                new_count++;
            }
        }
        
        Expr** final_args = malloc(sizeof(Expr*) * new_count);
        size_t idx = 0;
        for (size_t i = 0; i < expr->data.function.arg_count; i++) {
            if (!skip_flattening && temp_args[i]->type == EXPR_FUNCTION && 
                temp_args[i]->data.function.head->type == EXPR_SYMBOL &&
                temp_args[i]->data.function.head->data.symbol == SYM_Sequence) {
                for (size_t j = 0; j < temp_args[i]->data.function.arg_count; j++) {
                    final_args[idx++] = expr_copy(temp_args[i]->data.function.args[j]);
                }
                expr_free(temp_args[i]);
            } else {
                final_args[idx++] = temp_args[i];
            }
        }
        free(temp_args);
        
        Expr* result = expr_new_function(new_head, final_args, new_count);
        free(final_args);
        return result;
    }
    
    return expr_copy(expr);
}

Expr* builtin_matchq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    
    Expr* expr = evaluate(res->data.function.args[0]);
    Expr* pattern = res->data.function.args[1];
    
    MatchEnv* env = env_new();
    bool result = match(expr, pattern, env);
    env_free(env);
    expr_free(expr);
    
    return expr_new_symbol(result ? "True" : "False");
}
