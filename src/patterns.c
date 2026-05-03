#include "patterns.h"
#include "symtab.h"
#include "eval.h"
#include "match.h"
#include "core.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

static int64_t get_expr_depth_patterns(Expr* e, bool heads) {
    if (e->type != EXPR_FUNCTION) return 1;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol;
        if (h == SYM_Rational || h == SYM_Complex) return 1;
    }
    int64_t max_d = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        int64_t d = get_expr_depth_patterns(e->data.function.args[i], heads);
        if (d > max_d) max_d = d;
    }
    if (heads) {
        int64_t d_head = get_expr_depth_patterns(e->data.function.head, heads);
        if (d_head > max_d) max_d = d_head;
    }
    return max_d + 1;
}

static void do_cases_at_level(Expr* e, int64_t current_level, int64_t min_l, int64_t max_l, bool heads, Expr* pattern, Expr* replacement, bool delayed, Expr*** results, size_t* count, size_t* cap, int64_t max_results) {
    if (max_results >= 0 && (int64_t)(*count) >= max_results) return;

    if (e->type == EXPR_FUNCTION) {
        if (heads) {
            do_cases_at_level(e->data.function.head, current_level + 1, min_l, max_l, heads, pattern, replacement, delayed, results, count, cap, max_results);
        }
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (max_results >= 0 && (int64_t)(*count) >= max_results) break;
            do_cases_at_level(e->data.function.args[i], current_level + 1, min_l, max_l, heads, pattern, replacement, delayed, results, count, cap, max_results);
        }
    }

    if (max_results >= 0 && (int64_t)(*count) >= max_results) return;
    if (min_l >= 0 && max_l >= 0 && current_level > max_l) return;

    int64_t d = get_expr_depth_patterns(e, heads);

    bool match_level = true;
    if (min_l >= 0) {
        if (current_level < min_l || current_level > max_l) match_level = false;
    } else {
        if (min_l < 0 && max_l == min_l && d != -min_l) match_level = false;
        else if (min_l < 0 && max_l < 0 && (d < -max_l || d > -min_l)) match_level = false;
    }

    if (match_level) {
        MatchEnv* env = env_new();
        if (match(e, pattern, env)) {
            Expr* res;
            if (replacement) {
                Expr* repl_bound = replace_bindings(replacement, env);
                if (delayed) {
                    res = eval_and_free(repl_bound);
                } else {
                    res = repl_bound;
                }
            } else {
                res = expr_copy(e);
            }
            if (*count >= *cap) {
                *cap = (*cap == 0) ? 16 : (*cap * 2);
                *results = realloc(*results, sizeof(Expr*) * (*cap));
            }
            (*results)[(*count)++] = res;
        }
        env_free(env);
    }
}

Expr* builtin_cases(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    
    if (argc == 1) {
        Expr* slot_args[1] = { expr_new_integer(1) };
        Expr* slot = expr_new_function(expr_new_symbol("Slot"), slot_args, 1);
        Expr* inner_args[2] = { slot, expr_copy(res->data.function.args[0]) };
        Expr* inner_cases = expr_new_function(expr_new_symbol("Cases"), inner_args, 2);
        Expr* func_args[1] = { inner_cases };
        return expr_new_function(expr_new_symbol("Function"), func_args, 1);
    }
    
    if (argc < 2) return NULL;

    Expr* expr = res->data.function.args[0];
    Expr* patt_arg = res->data.function.args[1];

    int64_t min_l = 1, max_l = 1;
    bool heads = false;

    if (argc >= 3) {
        Expr* ls = res->data.function.args[2];
        if (ls->type == EXPR_INTEGER) {
            if (ls->data.integer < 0) {
                min_l = ls->data.integer; max_l = ls->data.integer;
            } else {
                min_l = 1; max_l = ls->data.integer;
            }
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol == SYM_All) {
            min_l = 1; max_l = 1000000;
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

    for (size_t i = 2; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->data.symbol == SYM_Rule && opt->data.function.arg_count == 2) {
            if (opt->data.function.args[0]->type == EXPR_SYMBOL && opt->data.function.args[0]->data.symbol == SYM_Heads) {
                if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol == SYM_True) heads = true;
                else if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol == SYM_False) heads = false;
            }
        }
    }

    int64_t max_results = -1;
    if (argc >= 4) {
        Expr* n_expr = res->data.function.args[3];
        if (n_expr->type == EXPR_INTEGER && n_expr->data.integer >= 0) {
            max_results = n_expr->data.integer;
        }
    }

    Expr* pattern = patt_arg;
    Expr* replacement = NULL;
    bool delayed = false;

    if (patt_arg->type == EXPR_FUNCTION && (patt_arg->data.function.head->data.symbol == SYM_Rule || patt_arg->data.function.head->data.symbol == SYM_RuleDelayed) && patt_arg->data.function.arg_count == 2) {
        pattern = patt_arg->data.function.args[0];
        replacement = patt_arg->data.function.args[1];
        delayed = (patt_arg->data.function.head->data.symbol == SYM_RuleDelayed);
    } else if (patt_arg->type == EXPR_FUNCTION && patt_arg->data.function.head->data.symbol == SYM_HoldPattern && patt_arg->data.function.arg_count == 1) {
        Expr* hp_arg = patt_arg->data.function.args[0];
        if (hp_arg->type == EXPR_FUNCTION && (hp_arg->data.function.head->data.symbol == SYM_Rule || hp_arg->data.function.head->data.symbol == SYM_RuleDelayed) && hp_arg->data.function.arg_count == 2) {
            pattern = hp_arg; // Just match the rule itself
        }
    }

    size_t count = 0;
    size_t cap = 16;
    Expr** results = malloc(sizeof(Expr*) * cap);

    do_cases_at_level(expr, 0, min_l, max_l, heads, pattern, replacement, delayed, &results, &count, &cap, max_results);

    Expr* list = expr_new_function(expr_new_symbol("List"), results, count);
    free(results);
    return list;
}

static void do_position_at_level(Expr* e, int64_t current_level, int64_t min_l, int64_t max_l, bool heads, Expr* pattern, Expr*** results, size_t* count, size_t* cap, int64_t max_results, int64_t* current_path, size_t path_len) {
    if (max_results >= 0 && (int64_t)(*count) >= max_results) return;

    if (e->type == EXPR_FUNCTION) {
        int64_t* next_path = malloc(sizeof(int64_t) * (path_len + 1));
        if (path_len > 0) memcpy(next_path, current_path, sizeof(int64_t) * path_len);
        
        if (heads) {
            next_path[path_len] = 0;
            do_position_at_level(e->data.function.head, current_level + 1, min_l, max_l, heads, pattern, results, count, cap, max_results, next_path, path_len + 1);
        }
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (max_results >= 0 && (int64_t)(*count) >= max_results) break;
            next_path[path_len] = i + 1;
            do_position_at_level(e->data.function.args[i], current_level + 1, min_l, max_l, heads, pattern, results, count, cap, max_results, next_path, path_len + 1);
        }
        free(next_path);
    }

    if (max_results >= 0 && (int64_t)(*count) >= max_results) return;
    if (min_l >= 0 && max_l >= 0 && current_level > max_l) return;

    int64_t d = get_expr_depth_patterns(e, heads);

    bool match_level = true;
    if (min_l >= 0) {
        if (current_level < min_l || current_level > max_l) match_level = false;
    } else {
        if (min_l < 0 && max_l == min_l && d != -min_l) match_level = false;
        else if (min_l < 0 && max_l < 0 && (d < -max_l || d > -min_l)) match_level = false;
    }

    if (match_level) {
        MatchEnv* env = env_new();
        if (match(e, pattern, env)) {
            // Add path to results
            Expr** path_exprs = malloc(sizeof(Expr*) * path_len);
            for (size_t i = 0; i < path_len; i++) {
                path_exprs[i] = expr_new_integer(current_path[i]);
            }
            Expr* pos_expr = expr_new_function(expr_new_symbol("List"), path_exprs, path_len);
            free(path_exprs);
            
            if (*count >= *cap) {
                *cap = (*cap == 0) ? 16 : (*cap * 2);
                *results = realloc(*results, sizeof(Expr*) * (*cap));
            }
            (*results)[(*count)++] = pos_expr;
        }
        env_free(env);
    }
}

Expr* builtin_position(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    
    if (argc == 1) {
        Expr* slot_args[1] = { expr_new_integer(1) };
        Expr* slot = expr_new_function(expr_new_symbol("Slot"), slot_args, 1);
        Expr* inner_args[2] = { slot, expr_copy(res->data.function.args[0]) };
        Expr* inner_pos = expr_new_function(expr_new_symbol("Position"), inner_args, 2);
        Expr* func_args[1] = { inner_pos };
        return expr_new_function(expr_new_symbol("Function"), func_args, 1);
    }
    
    if (argc < 2) return NULL;

    Expr* expr = res->data.function.args[0];
    Expr* pattern = res->data.function.args[1];

    int64_t min_l = 0, max_l = 1000000;
    bool heads = true;

    if (argc >= 3) {
        Expr* ls = res->data.function.args[2];
        if (ls->type == EXPR_INTEGER) {
            if (ls->data.integer < 0) {
                min_l = ls->data.integer; max_l = ls->data.integer;
            } else {
                min_l = 1; max_l = ls->data.integer;
            }
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol == SYM_All) {
            min_l = 1; max_l = 1000000;
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

    for (size_t i = 2; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->data.symbol == SYM_Rule && opt->data.function.arg_count == 2) {
            if (opt->data.function.args[0]->type == EXPR_SYMBOL && opt->data.function.args[0]->data.symbol == SYM_Heads) {
                if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol == SYM_True) heads = true;
                else if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol == SYM_False) heads = false;
            }
        }
    }

    int64_t max_results = -1;
    if (argc >= 4) {
        Expr* n_expr = res->data.function.args[3];
        if (n_expr->type == EXPR_INTEGER && n_expr->data.integer >= 0) {
            max_results = n_expr->data.integer;
        }
    }

    size_t count = 0;
    size_t cap = 16;
    Expr** results = malloc(sizeof(Expr*) * cap);

    do_position_at_level(expr, 0, min_l, max_l, heads, pattern, &results, &count, &cap, max_results, NULL, 0);

    Expr* list = expr_new_function(expr_new_symbol("List"), results, count);
    free(results);
    return list;
}

static void do_count_at_level(Expr* e, int64_t current_level, int64_t min_l, int64_t max_l, bool heads, Expr* pattern, size_t* count) {
    if (e->type == EXPR_FUNCTION) {
        if (heads) {
            do_count_at_level(e->data.function.head, current_level + 1, min_l, max_l, heads, pattern, count);
        }
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            do_count_at_level(e->data.function.args[i], current_level + 1, min_l, max_l, heads, pattern, count);
        }
    }

    if (min_l >= 0 && max_l >= 0 && current_level > max_l) return;

    int64_t d = get_expr_depth_patterns(e, heads);

    bool match_level = true;
    if (min_l >= 0) {
        if (current_level < min_l || current_level > max_l) match_level = false;
    } else {
        if (min_l < 0 && max_l == min_l && d != -min_l) match_level = false;
        else if (min_l < 0 && max_l < 0 && (d < -max_l || d > -min_l)) match_level = false;
    }

    if (match_level) {
        MatchEnv* env = env_new();
        if (match(e, pattern, env)) {
            (*count)++;
        }
        env_free(env);
    }
}

Expr* builtin_count(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    
    if (argc == 1) {
        Expr* slot_args[1] = { expr_new_integer(1) };
        Expr* slot = expr_new_function(expr_new_symbol("Slot"), slot_args, 1);
        Expr* inner_args[2] = { slot, expr_copy(res->data.function.args[0]) };
        Expr* inner_count = expr_new_function(expr_new_symbol("Count"), inner_args, 2);
        Expr* func_args[1] = { inner_count };
        return expr_new_function(expr_new_symbol("Function"), func_args, 1);
    }
    
    if (argc < 2) return NULL;

    Expr* expr = res->data.function.args[0];
    Expr* pattern = res->data.function.args[1];

    int64_t min_l = 1, max_l = 1;
    bool heads = false;

    if (argc >= 3) {
        Expr* ls = res->data.function.args[2];
        if (ls->type == EXPR_INTEGER) {
            if (ls->data.integer < 0) {
                min_l = ls->data.integer; max_l = ls->data.integer;
            } else {
                min_l = 1; max_l = ls->data.integer;
            }
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol == SYM_All) {
            min_l = 1; max_l = 1000000;
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

    for (size_t i = 2; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->data.symbol == SYM_Rule && opt->data.function.arg_count == 2) {
            if (opt->data.function.args[0]->type == EXPR_SYMBOL && opt->data.function.args[0]->data.symbol == SYM_Heads) {
                if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol == SYM_True) heads = true;
                else if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol == SYM_False) heads = false;
            }
        }
    }

    size_t count = 0;
    do_count_at_level(expr, 0, min_l, max_l, heads, pattern, &count);

    return expr_new_integer(count);
}


static bool do_member_at_level(Expr* e, int64_t current_level, int64_t min_l, int64_t max_l, bool heads, Expr* pattern) {
    if (min_l >= 0 && max_l >= 0 && current_level > max_l) return false;

    int64_t d = get_expr_depth_patterns(e, heads);

    bool match_level = true;
    if (min_l >= 0) {
        if (current_level < min_l || current_level > max_l) match_level = false;
    } else {
        if (min_l < 0 && max_l == min_l && d != -min_l) match_level = false;
        else if (min_l < 0 && max_l < 0 && (d < -max_l || d > -min_l)) match_level = false;
    }

    if (match_level) {
        MatchEnv* env = env_new();
        if (match(e, pattern, env)) {
            env_free(env);
            return true;
        }
        env_free(env);
    }

    if (e->type == EXPR_FUNCTION) {
        if (heads) {
            if (do_member_at_level(e->data.function.head, current_level + 1, min_l, max_l, heads, pattern)) return true;
        }
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (do_member_at_level(e->data.function.args[i], current_level + 1, min_l, max_l, heads, pattern)) return true;
        }
    }
    return false;
}

Expr* builtin_memberq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    
    if (argc == 1) {
        Expr* slot_args[1] = { expr_new_integer(1) };
        Expr* slot = expr_new_function(expr_new_symbol("Slot"), slot_args, 1);
        Expr* inner_args[2] = { slot, expr_copy(res->data.function.args[0]) };
        Expr* inner_memberq = expr_new_function(expr_new_symbol("MemberQ"), inner_args, 2);
        Expr* func_args[1] = { inner_memberq };
        return expr_new_function(expr_new_symbol("Function"), func_args, 1);
    }
    
    if (argc < 2) return NULL;

    Expr* expr = res->data.function.args[0];
    Expr* pattern = res->data.function.args[1];

    int64_t min_l = 1, max_l = 1;
    bool heads = false;

    if (argc >= 3) {
        Expr* ls = res->data.function.args[2];
        if (ls->type == EXPR_INTEGER) {
            if (ls->data.integer < 0) {
                min_l = ls->data.integer; max_l = ls->data.integer;
            } else {
                min_l = 1; max_l = ls->data.integer;
            }
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol == SYM_All) {
            min_l = 1; max_l = 1000000;
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

    for (size_t i = 2; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->data.symbol == SYM_Rule && opt->data.function.arg_count == 2) {
            if (opt->data.function.args[0]->type == EXPR_SYMBOL && opt->data.function.args[0]->data.symbol == SYM_Heads) {
                if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol == SYM_True) heads = true;
                else if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol == SYM_False) heads = false;
            }
        }
    }

    if (do_member_at_level(expr, 0, min_l, max_l, heads, pattern)) {
        return expr_new_symbol("True");
    } else {
        return expr_new_symbol("False");
    }
}

void patterns_init(void) {
    symtab_add_builtin("Cases", builtin_cases);
    symtab_get_def("Cases")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Position", builtin_position);
    symtab_get_def("Position")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Count", builtin_count);
    symtab_get_def("Count")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("MemberQ", builtin_memberq);
    symtab_get_def("MemberQ")->attributes |= ATTR_PROTECTED;
}