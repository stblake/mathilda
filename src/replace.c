#include "replace.h"
#include "symtab.h"
#include "eval.h"
#include "match.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

typedef struct {
    Expr* pos_spec;
    Expr* replacement;
    bool delayed;
} RawRule;

static void add_raw_rule(RawRule** rules, size_t* cap, size_t* count, Expr* path, Expr* replacement, bool delayed) {
    if (*count >= *cap) {
        *cap = (*cap == 0) ? 4 : (*cap * 2);
        *rules = realloc(*rules, sizeof(RawRule) * (*cap));
    }
    (*rules)[*count].pos_spec = expr_copy(path);
    (*rules)[*count].replacement = expr_copy(replacement);
    (*rules)[*count].delayed = delayed;
    (*count)++;
}

static void parse_single_rule(Expr* rule_expr, RawRule** rules, size_t* cap, size_t* count) {
    if (rule_expr->type != EXPR_FUNCTION) return;
    bool delayed = false;
    if (rule_expr->data.function.head->data.symbol == SYM_Rule) delayed = false;
    else if (rule_expr->data.function.head->data.symbol == SYM_RuleDelayed) delayed = true;
    else return;

    if (rule_expr->data.function.arg_count != 2) return;
    Expr* lhs = rule_expr->data.function.args[0];
    Expr* rhs = rule_expr->data.function.args[1];

    if (lhs->type == EXPR_FUNCTION && lhs->data.function.head->data.symbol == SYM_List) {
        if (lhs->data.function.arg_count > 0 && 
            lhs->data.function.args[0]->type == EXPR_FUNCTION &&
            lhs->data.function.args[0]->data.function.head->data.symbol == SYM_List) {
            for (size_t i = 0; i < lhs->data.function.arg_count; i++) {
                add_raw_rule(rules, cap, count, lhs->data.function.args[i], rhs, delayed);
            }
        } else {
            add_raw_rule(rules, cap, count, lhs, rhs, delayed);
        }
    } else {
        Expr* lhs_copy = expr_copy(lhs);
        Expr* wrapped = expr_new_function(expr_new_symbol(SYM_List), &lhs_copy, 1);
        add_raw_rule(rules, cap, count, wrapped, rhs, delayed);
        expr_free(wrapped);
    }
}

typedef struct {
    int64_t pos_index;
    int64_t neg_index;
} PathElement;

/*
 * normalize_pos_spec:
 * Convert negative integer indices in a rule's position spec to their positive
 * equivalents using the known path lengths. Returns a new expression that must
 * be freed, or NULL if no normalization was needed.
 */
static Expr* normalize_pos_spec(Expr* pos_spec, PathElement* path, size_t depth) {
    if (pos_spec->type != EXPR_FUNCTION || depth == 0) return NULL;
    if (pos_spec->data.function.arg_count != depth) return NULL;

    bool needs_norm = false;
    for (size_t d = 0; d < depth; d++) {
        Expr* elem = pos_spec->data.function.args[d];
        if (elem->type == EXPR_INTEGER && elem->data.integer < 0) {
            needs_norm = true;
            break;
        }
    }
    if (!needs_norm) return NULL;

    Expr** new_args = malloc(sizeof(Expr*) * depth);
    for (size_t d = 0; d < depth; d++) {
        Expr* elem = pos_spec->data.function.args[d];
        if (elem->type == EXPR_INTEGER && elem->data.integer < 0) {
            /* Recover length from pos_index and neg_index:
             * neg_index = pos_index - length - 1, so length = pos_index - neg_index - 1 */
            int64_t length = path[d].pos_index - path[d].neg_index - 1;
            int64_t positive = elem->data.integer + length + 1;
            new_args[d] = expr_new_integer(positive);
        } else {
            new_args[d] = expr_copy(elem);
        }
    }
    Expr* result = expr_new_function(expr_copy(pos_spec->data.function.head), new_args, depth);
    free(new_args);
    return result;
}

static bool check_match(RawRule* rules, size_t num_rules, PathElement* path, size_t depth,
                        Expr** out_replacement, bool* out_delayed, MatchEnv** out_env) {
    for (int64_t r = (int64_t)num_rules - 1; r >= 0; r--) {
        /* Build the current position as a List of positive indices */
        Expr** alias_args = NULL;
        if (depth > 0) alias_args = malloc(sizeof(Expr*) * depth);
        for (size_t d = 0; d < depth; d++) {
            alias_args[d] = expr_new_integer(path[d].pos_index);
        }
        Expr* alias_expr = expr_new_function(expr_new_symbol(SYM_List), alias_args, depth);

        /* Normalize any negative indices in the rule's pos_spec to positive */
        Expr* norm_spec = normalize_pos_spec(rules[r].pos_spec, path, depth);
        Expr* spec_to_match = norm_spec ? norm_spec : rules[r].pos_spec;

        MatchEnv* env = env_new();
        if (match(alias_expr, spec_to_match, env)) {
            *out_replacement = rules[r].replacement;
            *out_delayed = rules[r].delayed;
            *out_env = env;
            expr_free(alias_expr);
            if (alias_args) free(alias_args);
            if (norm_spec) expr_free(norm_spec);
            return true;
        }
        env_free(env);
        expr_free(alias_expr);
        if (alias_args) free(alias_args);
        if (norm_spec) expr_free(norm_spec);
    }
    return false;
}

static Expr* replace_part_rec(Expr* expr, RawRule* rules, size_t num_rules, PathElement* current_path, size_t depth) {
    Expr* replacement_expr = NULL;
    bool delayed = false;
    MatchEnv* env = NULL;

    if (check_match(rules, num_rules, current_path, depth, &replacement_expr, &delayed, &env)) {
        Expr* res = replace_bindings(replacement_expr, env);
        if (delayed) {
            Expr* eval_res = evaluate(res);
            expr_free(res);
            res = eval_res;
        }
        env_free(env);
        return res;
    }

    if (expr->type == EXPR_FUNCTION) {
        size_t len = expr->data.function.arg_count;
        Expr** new_args = NULL;
        if (len > 0) new_args = malloc(sizeof(Expr*) * len);

        PathElement* next_path = malloc(sizeof(PathElement) * (depth + 1));
        if (depth > 0) memcpy(next_path, current_path, sizeof(PathElement) * depth);

        /* Check head (index 0) -- only replace if a rule explicitly targets index 0
         * with a literal integer. Default Heads -> False means patterns like _ or
         * Except[...] should not match the head position. */
        next_path[depth].pos_index = 0;
        next_path[depth].neg_index = 0;
        Expr* new_head = NULL;
        bool head_replaced = false;
        for (size_t r = num_rules; r > 0; r--) {
            Expr* spec = rules[r - 1].pos_spec;
            /* Only consider this rule for the head if it has a literal 0 at this depth */
            if (spec->type == EXPR_FUNCTION && spec->data.function.arg_count == depth + 1) {
                Expr* idx = spec->data.function.args[depth];
                if (idx->type == EXPR_INTEGER && idx->data.integer == 0) {
                    new_head = replace_part_rec(expr->data.function.head, rules, num_rules, next_path, depth + 1);
                    head_replaced = true;
                    break;
                }
            }
        }
        if (!head_replaced) {
            new_head = expr_copy(expr->data.function.head);
        }

        for (size_t i = 0; i < len; i++) {
            next_path[depth].pos_index = i + 1;
            next_path[depth].neg_index = (int64_t)(i + 1) - (int64_t)len - 1;
            new_args[i] = replace_part_rec(expr->data.function.args[i], rules, num_rules, next_path, depth + 1);
        }
        
        free(next_path);
        
        Expr* result = expr_new_function(new_head, new_args, len);
        if (new_args) free(new_args);
        return result;
    }
    
    return expr_copy(expr);
}

Expr* builtin_replace_part(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count == 1) {
        // Operator form ReplacePart[rules]
        return NULL; // Stays unevaluated
    } else if (res->data.function.arg_count == 2) {
        Expr* expr = res->data.function.args[0];
        Expr* rules_expr = res->data.function.args[1];

        size_t cap = 4;
        size_t num_rules = 0;
        RawRule* rules = malloc(sizeof(RawRule) * cap);

        if (rules_expr->type == EXPR_FUNCTION && rules_expr->data.function.head->data.symbol == SYM_List) {
            for (size_t i = 0; i < rules_expr->data.function.arg_count; i++) {
                parse_single_rule(rules_expr->data.function.args[i], &rules, &cap, &num_rules);
            }
        } else {
            parse_single_rule(rules_expr, &rules, &cap, &num_rules);
        }

        if (num_rules == 0) {
            free(rules);
            return expr_copy(expr);
        }

        Expr* result = replace_part_rec(expr, rules, num_rules, NULL, 0);

        for (size_t i = 0; i < num_rules; i++) {
            expr_free(rules[i].pos_spec);
            expr_free(rules[i].replacement);
        }
        free(rules);

        return result;
    }
    return NULL;
}

typedef struct {
    Expr* pattern;
    Expr* replacement;
    bool delayed;
} ReplaceRule;

static int64_t get_expr_depth_replace(Expr* e, bool heads) {
    if (e->type != EXPR_FUNCTION) return 1;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol;
        if (h == SYM_Rational || h == SYM_Complex) return 1;
    }
    int64_t max_d = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        int64_t d = get_expr_depth_replace(e->data.function.args[i], heads);
        if (d > max_d) max_d = d;
    }
    if (heads) {
        int64_t d_head = get_expr_depth_replace(e->data.function.head, heads);
        if (d_head > max_d) max_d = d_head;
    }
    return max_d + 1;
}

static Expr* do_replace_at_level(Expr* e, int64_t current_level, int64_t min_l, int64_t max_l, bool heads, ReplaceRule* rules, size_t num_rules) {
    if (max_l >= 0 && current_level > max_l) return expr_copy(e);

    int64_t d = get_expr_depth_replace(e, heads);

    bool match_level = true;
    if (min_l >= 0) {
        if (current_level < min_l) match_level = false;
    } else {
        if (d > -min_l) match_level = false;
    }

    if (max_l >= 0) {
        if (current_level > max_l) match_level = false;
    } else {
        if (d < -max_l) match_level = false;
    }

    Expr* new_e = NULL;
    if (e->type == EXPR_FUNCTION) {
        Expr* new_head = e->data.function.head;
        if (heads) {
            new_head = do_replace_at_level(e->data.function.head, current_level + 1, min_l, max_l, heads, rules, num_rules);
        } else {
            new_head = expr_copy(new_head);
        }

        size_t count = e->data.function.arg_count;
        Expr** new_args = NULL;
        if (count > 0) new_args = malloc(sizeof(Expr*) * count);
        for (size_t i = 0; i < count; i++) {
            new_args[i] = do_replace_at_level(e->data.function.args[i], current_level + 1, min_l, max_l, heads, rules, num_rules);
        }

        new_e = expr_new_function(new_head, new_args, count);
        if (new_args) free(new_args);
    } else {
        new_e = expr_copy(e);
    }

    if (match_level) {
        for (size_t i = 0; i < num_rules; i++) {
            MatchEnv* env = env_new();
            if (match(new_e, rules[i].pattern, env)) {
                Expr* repl = replace_bindings(rules[i].replacement, env);
                env_free(env);
                expr_free(new_e);
                return repl;
            }
            env_free(env);
        }
    }

    return new_e;
}

static bool is_rule(Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol;
    return (h == SYM_Rule || h == SYM_RuleDelayed);
}

static Expr* apply_replace_nested(Expr* expr, Expr* rules_expr, int64_t min_l, int64_t max_l, bool heads) {
    if (rules_expr->type == EXPR_FUNCTION && rules_expr->data.function.head->data.symbol == SYM_List) {
        if (rules_expr->data.function.arg_count > 0 && !is_rule(rules_expr->data.function.args[0])) {
            size_t count = rules_expr->data.function.arg_count;
            Expr** args = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                args[i] = apply_replace_nested(expr, rules_expr->data.function.args[i], min_l, max_l, heads);
            }
            Expr* res = expr_new_function(expr_new_symbol(SYM_List), args, count);
            free(args);
            return res;
        }
    }
    
    size_t cap = 4;
    size_t num_rules = 0;
    ReplaceRule* rules = malloc(sizeof(ReplaceRule) * cap);
    
    if (rules_expr->type == EXPR_FUNCTION && rules_expr->data.function.head->data.symbol == SYM_List) {
        for (size_t i = 0; i < rules_expr->data.function.arg_count; i++) {
            Expr* r = rules_expr->data.function.args[i];
            if (is_rule(r) && r->data.function.arg_count == 2) {
                if (num_rules == cap) { cap *= 2; rules = realloc(rules, sizeof(ReplaceRule) * cap); }
                rules[num_rules].pattern = r->data.function.args[0];
                rules[num_rules].replacement = r->data.function.args[1];
                rules[num_rules].delayed = r->data.function.head->data.symbol == SYM_RuleDelayed;
                num_rules++;
            }
        }
    } else if (is_rule(rules_expr) && rules_expr->data.function.arg_count == 2) {
        rules[num_rules].pattern = rules_expr->data.function.args[0];
        rules[num_rules].replacement = rules_expr->data.function.args[1];
        rules[num_rules].delayed = rules_expr->data.function.head->data.symbol == SYM_RuleDelayed;
        num_rules++;
    }
    
    Expr* result;
    if (num_rules > 0) {
        result = do_replace_at_level(expr, 0, min_l, max_l, heads, rules, num_rules);
    } else {
        result = expr_copy(expr);
    }
    free(rules);
    return result;
}

Expr* builtin_replace(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* rules_expr = res->data.function.args[1];
    
    int64_t min_l = 0, max_l = 0;
    bool heads = false;
    
    if (res->data.function.arg_count >= 3) {
        Expr* ls = res->data.function.args[2];
        if (ls->type == EXPR_INTEGER) {
            min_l = 1; max_l = ls->data.integer;
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol == SYM_All) {
            min_l = 0; max_l = 1000000;
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol == SYM_Infinity) {
            min_l = 1; max_l = 1000000;
        } else if (ls->type == EXPR_FUNCTION && ls->data.function.head->data.symbol == SYM_List) {
            if (ls->data.function.arg_count == 1 && ls->data.function.args[0]->type == EXPR_INTEGER) {
                min_l = max_l = ls->data.function.args[0]->data.integer;
            } else if (ls->data.function.arg_count == 2) {
                if (ls->data.function.args[0]->type == EXPR_INTEGER) min_l = ls->data.function.args[0]->data.integer;
                if (ls->data.function.args[1]->type == EXPR_INTEGER) max_l = ls->data.function.args[1]->data.integer;
                else if (ls->data.function.args[1]->type == EXPR_SYMBOL && ls->data.function.args[1]->data.symbol == SYM_Infinity) max_l = 1000000;
            }
        }
    }
    
    for (size_t i = 2; i < res->data.function.arg_count; i++) {
        Expr* opt = res->data.function.args[i];
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->data.symbol == SYM_Rule && opt->data.function.arg_count == 2) {
            if (opt->data.function.args[0]->type == EXPR_SYMBOL && opt->data.function.args[0]->data.symbol == SYM_Heads) {
                if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol == SYM_True) heads = true;
                else if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol == SYM_False) heads = false;
            }
        }
    }
    
    return apply_replace_nested(expr, rules_expr, min_l, max_l, heads);
}

static Expr* do_replace_all(Expr* e, ReplaceRule* rules, size_t num_rules) {
    for (size_t i = 0; i < num_rules; i++) {
        MatchEnv* env = env_new();
        if (match(e, rules[i].pattern, env)) {
            Expr* repl = replace_bindings(rules[i].replacement, env);
            env_free(env);
            return repl;
        }
        env_free(env);
    }

    if (e->type == EXPR_FUNCTION) {
        Expr* new_head = do_replace_all(e->data.function.head, rules, num_rules);
        size_t count = e->data.function.arg_count;
        Expr** new_args = NULL;
        if (count > 0) new_args = malloc(sizeof(Expr*) * count);
        for (size_t i = 0; i < count; i++) {
            new_args[i] = do_replace_all(e->data.function.args[i], rules, num_rules);
        }
        Expr* new_e = expr_new_function(new_head, new_args, count);
        if (new_args) free(new_args);
        return new_e;
    }

    return expr_copy(e);
}

static Expr* apply_replace_all_nested(Expr* expr, Expr* rules_expr) {
    if (rules_expr->type == EXPR_FUNCTION && rules_expr->data.function.head->data.symbol == SYM_List) {
        if (rules_expr->data.function.arg_count > 0 && !is_rule(rules_expr->data.function.args[0])) {
            size_t count = rules_expr->data.function.arg_count;
            Expr** args = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                args[i] = apply_replace_all_nested(expr, rules_expr->data.function.args[i]);
            }
            Expr* res = expr_new_function(expr_new_symbol(SYM_List), args, count);
            free(args);
            return res;
        }
    }
    
    size_t cap = 4;
    size_t num_rules = 0;
    ReplaceRule* rules = malloc(sizeof(ReplaceRule) * cap);
    
    if (rules_expr->type == EXPR_FUNCTION && rules_expr->data.function.head->data.symbol == SYM_List) {
        for (size_t i = 0; i < rules_expr->data.function.arg_count; i++) {
            Expr* r = rules_expr->data.function.args[i];
            if (is_rule(r) && r->data.function.arg_count == 2) {
                if (num_rules == cap) { cap *= 2; rules = realloc(rules, sizeof(ReplaceRule) * cap); }
                rules[num_rules].pattern = r->data.function.args[0];
                rules[num_rules].replacement = r->data.function.args[1];
                rules[num_rules].delayed = r->data.function.head->data.symbol == SYM_RuleDelayed;
                num_rules++;
            }
        }
    } else if (is_rule(rules_expr) && rules_expr->data.function.arg_count == 2) {
        rules[num_rules].pattern = rules_expr->data.function.args[0];
        rules[num_rules].replacement = rules_expr->data.function.args[1];
        rules[num_rules].delayed = rules_expr->data.function.head->data.symbol == SYM_RuleDelayed;
        num_rules++;
    }
    
    Expr* result;
    if (num_rules > 0) {
        result = do_replace_all(expr, rules, num_rules);
    } else {
        result = expr_copy(expr);
    }
    free(rules);
    return result;
}

Expr* builtin_replace_all(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* rules_expr = res->data.function.args[1];
    
    return apply_replace_all_nested(expr, rules_expr);
}

static Expr* apply_replace_repeated_nested(Expr* expr, Expr* rules_expr) {
    if (rules_expr->type == EXPR_FUNCTION && rules_expr->data.function.head->data.symbol == SYM_List) {
        if (rules_expr->data.function.arg_count > 0 && !is_rule(rules_expr->data.function.args[0])) {
            size_t count = rules_expr->data.function.arg_count;
            Expr** args = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                args[i] = apply_replace_repeated_nested(expr, rules_expr->data.function.args[i]);
            }
            Expr* res = expr_new_function(expr_new_symbol(SYM_List), args, count);
            free(args);
            return res;
        }
    }

    Expr* current = expr_copy(expr);
    int iterations = 0;
    while (iterations < 65536) {
        Expr* next_raw = apply_replace_all_nested(current, rules_expr);
        Expr* next = eval_and_free(next_raw);
        if (expr_eq(current, next)) {
            expr_free(next);
            return current;
        }
        expr_free(current);
        current = next;
        iterations++;
    }
    
    return current;
}

Expr* builtin_replace_repeated(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* rules_expr = res->data.function.args[1];
    
    return apply_replace_repeated_nested(expr, rules_expr);
}


typedef struct {
    Expr* replacement;
    bool delayed;
    Expr** results;
    size_t count;
    size_t cap;
    int64_t limit;
} ReplaceListState;

static bool replacelist_callback(struct MatchEnv* env, void* user_data) {
    ReplaceListState* state = (ReplaceListState*)user_data;
    
    if (state->limit >= 0 && (int64_t)state->count >= state->limit) return true;
    
    Expr* repl = replace_bindings(state->replacement, env);
    if (state->delayed) {
        repl = eval_and_free(repl);
    }
    
    if (state->count == state->cap) {
        state->cap *= 2;
        state->results = realloc(state->results, sizeof(Expr*) * state->cap);
    }
    state->results[state->count++] = repl;
    
    return false;
}

Expr* builtin_replacelist(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;

    Expr* expr = res->data.function.args[0];
    Expr* rules_expr = res->data.function.args[1];
    
    int64_t limit = -1;
    if (argc == 3) {
        Expr* l_expr = res->data.function.args[2];
        if (l_expr->type == EXPR_INTEGER) {
            limit = l_expr->data.integer;
        }
    }
    
    size_t cap = 16;
    size_t num_rules = 0;
    ReplaceRule* rules = malloc(sizeof(ReplaceRule) * cap);
    
    if (rules_expr->type == EXPR_FUNCTION && rules_expr->data.function.head->data.symbol == SYM_List) {
        for (size_t i = 0; i < rules_expr->data.function.arg_count; i++) {
            Expr* r = rules_expr->data.function.args[i];
            if (is_rule(r) && r->data.function.arg_count == 2) {
                if (num_rules == cap) { cap *= 2; rules = realloc(rules, sizeof(ReplaceRule) * cap); }
                rules[num_rules].pattern = r->data.function.args[0];
                rules[num_rules].replacement = r->data.function.args[1];
                rules[num_rules].delayed = r->data.function.head->data.symbol == SYM_RuleDelayed;
                num_rules++;
            }
        }
    } else if (is_rule(rules_expr) && rules_expr->data.function.arg_count == 2) {
        rules[num_rules].pattern = rules_expr->data.function.args[0];
        rules[num_rules].replacement = rules_expr->data.function.args[1];
        rules[num_rules].delayed = rules_expr->data.function.head->data.symbol == SYM_RuleDelayed;
        num_rules++;
    }
    
    ReplaceListState state;
    state.results = malloc(sizeof(Expr*) * 16);
    state.count = 0;
    state.cap = 16;
    state.limit = limit;
    
    for (size_t i = 0; i < num_rules; i++) {
        state.replacement = rules[i].replacement;
        state.delayed = rules[i].delayed;
        
        MatchEnv* env = env_new();
        env->callback = replacelist_callback;
        env->callback_data = &state;
        
        match(expr, rules[i].pattern, env);
        env_free(env);
        
        if (state.limit >= 0 && (int64_t)state.count >= state.limit) {
            break;
        }
    }
    
    free(rules);
    
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), state.results, state.count);
    free(state.results);
    return list;
}

/* ------------------- ReplaceAt ------------------- */

/*
 * Parse rules_expr (a Rule, RuleDelayed, or List of such) into a freshly
 * allocated array of ReplaceRule. The pattern/replacement pointers borrow
 * from rules_expr; the array itself must be freed by the caller. Sets
 * *out_count to the number of valid rules collected.
 */
static ReplaceRule* parse_replace_rules(Expr* rules_expr, size_t* out_count) {
    size_t cap = 4;
    size_t count = 0;
    ReplaceRule* rules = malloc(sizeof(ReplaceRule) * cap);

    if (rules_expr->type == EXPR_FUNCTION &&
        rules_expr->data.function.head->type == EXPR_SYMBOL &&
        rules_expr->data.function.head->data.symbol == SYM_List) {
        for (size_t i = 0; i < rules_expr->data.function.arg_count; i++) {
            Expr* r = rules_expr->data.function.args[i];
            if (is_rule(r) && r->data.function.arg_count == 2) {
                if (count == cap) { cap *= 2; rules = realloc(rules, sizeof(ReplaceRule) * cap); }
                rules[count].pattern = r->data.function.args[0];
                rules[count].replacement = r->data.function.args[1];
                rules[count].delayed = r->data.function.head->data.symbol == SYM_RuleDelayed;
                count++;
            }
        }
    } else if (is_rule(rules_expr) && rules_expr->data.function.arg_count == 2) {
        rules[count].pattern = rules_expr->data.function.args[0];
        rules[count].replacement = rules_expr->data.function.args[1];
        rules[count].delayed = rules_expr->data.function.head->data.symbol == SYM_RuleDelayed;
        count++;
    }
    *out_count = count;
    return rules;
}

/*
 * Try each rule in order against target. The first rule that matches wins.
 * Returns a freshly-owned expression: either the replacement (with bindings
 * substituted, evaluated for delayed rules) or expr_copy(target) when no
 * rule matches.
 */
static Expr* replaceat_apply_rules(ReplaceRule* rules, size_t num_rules, Expr* target) {
    for (size_t i = 0; i < num_rules; i++) {
        MatchEnv* env = env_new();
        if (match(target, rules[i].pattern, env)) {
            Expr* repl = replace_bindings(rules[i].replacement, env);
            if (rules[i].delayed) {
                repl = eval_and_free(repl);
            }
            env_free(env);
            return repl;
        }
        env_free(env);
    }
    return expr_copy(target);
}

/*
 * Walk path inside expr, applying rules at the leaf indicated by path.
 *
 * A path element may be:
 *   - integer k (>0 from start, <0 from end, 0 targets the head)
 *   - the symbol All (recurse into every child at this level)
 *   - Span[a, b] or Span[a, b, step]
 *
 * Returns a freshly-allocated expression. Out-of-range integers are silently
 * ignored, matching Mathematica's permissive behaviour.
 */
static Expr* replaceat_at_path(ReplaceRule* rules, size_t num_rules, Expr* expr, Expr** path, size_t plen) {
    if (plen == 0) {
        return replaceat_apply_rules(rules, num_rules, expr);
    }
    if (expr->type != EXPR_FUNCTION) {
        return expr_copy(expr);
    }

    Expr* idx = path[0];
    size_t len = expr->data.function.arg_count;

    Expr** new_args = NULL;
    if (len > 0) {
        new_args = malloc(sizeof(Expr*) * len);
        for (size_t i = 0; i < len; i++) {
            new_args[i] = expr_copy(expr->data.function.args[i]);
        }
    }
    Expr* new_head = expr_copy(expr->data.function.head);

    if (idx->type == EXPR_INTEGER) {
        int64_t k = idx->data.integer;
        if (k == 0) {
            Expr* r = replaceat_at_path(rules, num_rules, expr->data.function.head, path + 1, plen - 1);
            expr_free(new_head);
            new_head = r;
        } else {
            if (k < 0) k = (int64_t)len + k + 1;
            if (k >= 1 && k <= (int64_t)len) {
                Expr* r = replaceat_at_path(rules, num_rules, expr->data.function.args[k - 1], path + 1, plen - 1);
                expr_free(new_args[k - 1]);
                new_args[k - 1] = r;
            }
        }
    } else if (idx->type == EXPR_SYMBOL && idx->data.symbol == SYM_All) {
        for (size_t i = 0; i < len; i++) {
            Expr* r = replaceat_at_path(rules, num_rules, expr->data.function.args[i], path + 1, plen - 1);
            expr_free(new_args[i]);
            new_args[i] = r;
        }
    } else if (idx->type == EXPR_FUNCTION &&
               idx->data.function.head->type == EXPR_SYMBOL &&
               idx->data.function.head->data.symbol == SYM_Span) {
        int64_t start = 1, end = (int64_t)len, step = 1;
        size_t span_argc = idx->data.function.arg_count;
        if (span_argc >= 1) {
            Expr* a1 = idx->data.function.args[0];
            if (a1->type == EXPR_INTEGER) {
                start = a1->data.integer;
                if (start < 0) start = (int64_t)len + start + 1;
            } else if (a1->type == EXPR_SYMBOL && a1->data.symbol == SYM_All) {
                start = 1;
            }
        }
        if (span_argc >= 2) {
            Expr* a2 = idx->data.function.args[1];
            if (a2->type == EXPR_INTEGER) {
                end = a2->data.integer;
                if (end < 0) end = (int64_t)len + end + 1;
            } else if (a2->type == EXPR_SYMBOL && a2->data.symbol == SYM_All) {
                end = (int64_t)len;
            }
        }
        if (span_argc >= 3) {
            Expr* a3 = idx->data.function.args[2];
            if (a3->type == EXPR_INTEGER) step = a3->data.integer;
        }

        if (step > 0) {
            for (int64_t i = start; i <= end && i >= 1 && i <= (int64_t)len; i += step) {
                Expr* r = replaceat_at_path(rules, num_rules, expr->data.function.args[i - 1], path + 1, plen - 1);
                expr_free(new_args[i - 1]);
                new_args[i - 1] = r;
            }
        } else if (step < 0) {
            for (int64_t i = start; i >= end && i >= 1 && i <= (int64_t)len; i += step) {
                Expr* r = replaceat_at_path(rules, num_rules, expr->data.function.args[i - 1], path + 1, plen - 1);
                expr_free(new_args[i - 1]);
                new_args[i - 1] = r;
            }
        }
    }

    Expr* result = expr_new_function(new_head, new_args, len);
    if (new_args) free(new_args);
    return result;
}

/* True iff e is a List expression. */
static bool replaceat_is_list(Expr* e) {
    return e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol == SYM_List;
}

Expr* builtin_replace_at(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;

    Expr* expr = res->data.function.args[0];
    Expr* rules_expr = res->data.function.args[1];
    Expr* pos = res->data.function.args[2];

    size_t num_rules = 0;
    ReplaceRule* rules = parse_replace_rules(rules_expr, &num_rules);

    if (num_rules == 0) {
        free(rules);
        return expr_copy(expr);
    }

    /* Disambiguate single-path vs. multiple-paths:
     *   - multiple-paths iff pos is a non-empty List whose first element is itself a List
     *   - otherwise, single path (possibly wrapped in a List of indices) */
    bool multi = false;
    if (replaceat_is_list(pos) && pos->data.function.arg_count > 0 &&
        replaceat_is_list(pos->data.function.args[0])) {
        multi = true;
    }

    Expr* result;
    if (!multi) {
        Expr** path;
        size_t plen;
        if (replaceat_is_list(pos)) {
            plen = pos->data.function.arg_count;
            path = pos->data.function.args;
        } else {
            plen = 1;
            path = &pos;
        }
        result = replaceat_at_path(rules, num_rules, expr, path, plen);
    } else {
        /* Apply paths sequentially. Repeated positions apply rules repeatedly. */
        Expr* current = expr_copy(expr);
        for (size_t i = 0; i < pos->data.function.arg_count; i++) {
            Expr* sub = pos->data.function.args[i];
            Expr** path;
            size_t plen;
            if (replaceat_is_list(sub)) {
                plen = sub->data.function.arg_count;
                path = sub->data.function.args;
            } else {
                plen = 1;
                path = &sub;
            }
            Expr* next = replaceat_at_path(rules, num_rules, current, path, plen);
            expr_free(current);
            current = next;
        }
        result = current;
    }

    free(rules);
    return result;
}

void replace_init(void) {
    symtab_add_builtin("ReplacePart", builtin_replace_part);
    symtab_add_builtin("Replace", builtin_replace);
    symtab_get_def("Replace")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("ReplaceAll", builtin_replace_all);
    symtab_get_def("ReplaceAll")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("ReplaceRepeated", builtin_replace_repeated);
    symtab_get_def("ReplaceRepeated")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("ReplaceList", builtin_replacelist);
    symtab_get_def("ReplaceList")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("ReplaceAt", builtin_replace_at);
    symtab_get_def("ReplaceAt")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ReplaceAt",
        "ReplaceAt[expr, rules, n]\n"
        "\ttransforms expr by replacing the n-th element using rules.\n"
        "ReplaceAt[expr, rules, {i, j, ...}]\n"
        "\treplaces the part of expr at position {i, j, ...}.\n"
        "ReplaceAt[expr, rules, {{i1, j1, ...}, {i2, j2, ...}, ...}]\n"
        "\treplaces parts at several positions.\n"
        "\n"
        "Rules may be a single Rule/RuleDelayed or a list of them; rules are tried\n"
        "in order and the first match wins. Negative indices count from the end;\n"
        "0 targets the head. All and Span specifications are supported. Repeated\n"
        "positions cause rules to be applied repeatedly to that part.");
}
