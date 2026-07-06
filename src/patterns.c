#include "patterns.h"
#include "symtab.h"
#include "eval.h"
#include "match.h"
#include "core.h"
#include "sym_names.h"
#include "assoc.h"
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

/* FirstCase[expr, patt] / FirstCase[expr, patt, default] — the first element
 * matching patt (or the first matching value, for an association), else the
 * default or Missing["NotFound"]. Reuses Cases (which already handles the
 * pattern and association-value threading) and takes the first match. */
Expr* builtin_first_case(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2 && argc != 3) return NULL;

    Expr* cases_args[2] = { expr_copy(res->data.function.args[0]),
                            expr_copy(res->data.function.args[1]) };
    Expr* cases = expr_new_function(expr_new_symbol(SYM_Cases), cases_args, 2);
    Expr* matches = evaluate(cases);
    expr_free(cases);

    Expr* result = NULL;
    if (matches && matches->type == EXPR_FUNCTION &&
        matches->data.function.head->data.symbol == SYM_List &&
        matches->data.function.arg_count >= 1) {
        result = expr_copy(matches->data.function.args[0]);
    }
    if (matches) expr_free(matches);
    if (result) return result;

    if (argc == 3) return expr_copy(res->data.function.args[2]);
    Expr* margs[1] = { expr_new_string("NotFound") };
    return expr_new_function(expr_new_symbol(SYM_Missing), margs, 1);
}

/* DeleteMissing[expr] — remove all Missing[...] elements. Equivalent to
 * DeleteCases[expr, _Missing], so it inherits list and association-value
 * handling (over an association it drops entries whose value is Missing[...]). */
Expr* builtin_delete_missing(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    /* pattern _Missing == Blank[Missing] */
    Expr* blank_args[1] = { expr_new_symbol(SYM_Missing) };
    Expr* pattern = expr_new_function(expr_new_symbol(SYM_Blank), blank_args, 1);
    Expr* dc_args[2] = { expr_copy(res->data.function.args[0]), pattern };
    Expr* dc = expr_new_function(expr_new_symbol(SYM_DeleteCases), dc_args, 2);
    Expr* result = evaluate(dc);
    expr_free(dc);
    return result;
}

Expr* builtin_cases(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    
    if (argc == 1) {
        Expr* slot_args[1] = { expr_new_integer(1) };
        Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot), slot_args, 1);
        Expr* inner_args[2] = { slot, expr_copy(res->data.function.args[0]) };
        Expr* inner_cases = expr_new_function(expr_new_symbol(SYM_Cases), inner_args, 2);
        Expr* func_args[1] = { inner_cases };
        return expr_new_function(expr_new_symbol(SYM_Function), func_args, 1);
    }

    /* Cases[assoc, patt, ...] collects matching values (Cases[Values[assoc], ...]). */
    if (argc >= 2 && is_association(res->data.function.args[0])) {
        Expr* r = assoc_apply_over_values(res); if (r) return r;
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

    Expr* list = expr_new_function(expr_new_symbol(SYM_List), results, count);
    free(results);
    return list;
}

/*
 * do_delete_cases_at_level
 *
 * Returns a freshly-allocated transformed copy of `e` with descendants matching
 * `pattern` (within levelspec [min_l, max_l]) removed. Traversal is depth-first
 * post-order ("leaves before roots"), mirroring `do_cases_at_level`.
 *
 * Out-parameter `*delete_me` is set to true when the parent should drop `e`
 * from its argument list. Mirroring Mathematica's semantics, the match test is
 * applied to the ORIGINAL `e` (not the transformed copy), so a parent that
 * matches the pattern is dropped even if some of its children were already
 * removed; the work done on the discarded children is harmless but charged
 * against the `n` budget like any other match.
 *
 * Heads -> True (heads = true) lets the head of a function be tested too. A
 * matching head behaves like FlattenAt: the function call is replaced by
 * Sequence[args...], which the surrounding loop splices into the parent.
 *
 * `*count_remaining` follows the same convention as Cases:
 *   -1 -> unlimited deletions, > 0 -> deletions left, 0 -> budget exhausted.
 * It is decremented only for positive budgets.
 */
static Expr* do_delete_cases_at_level(Expr* e, int64_t current_level, int64_t min_l, int64_t max_l, bool heads, Expr* pattern, int64_t* count_remaining, bool* delete_me) {
    *delete_me = false;
    Expr* result;

    if (e->type == EXPR_FUNCTION) {
        size_t orig_count = e->data.function.arg_count;
        size_t cap = orig_count > 0 ? orig_count : 1;
        Expr** new_args = malloc(sizeof(Expr*) * cap);
        size_t new_count = 0;

        Expr* new_head;
        bool head_delete = false;
        if (heads) {
            new_head = do_delete_cases_at_level(e->data.function.head, current_level + 1, min_l, max_l, heads, pattern, count_remaining, &head_delete);
        } else {
            new_head = expr_copy(e->data.function.head);
        }

        for (size_t i = 0; i < orig_count; i++) {
            bool arg_delete = false;
            Expr* new_arg = do_delete_cases_at_level(e->data.function.args[i], current_level + 1, min_l, max_l, heads, pattern, count_remaining, &arg_delete);
            if (arg_delete) {
                expr_free(new_arg);
                continue;
            }
            /* Splice Sequence[...] inline so head-deletions flatten outwards. */
            if (new_arg->type == EXPR_FUNCTION &&
                new_arg->data.function.head->type == EXPR_SYMBOL &&
                new_arg->data.function.head->data.symbol == SYM_Sequence) {
                size_t seq_count = new_arg->data.function.arg_count;
                while (new_count + seq_count > cap) {
                    cap = cap * 2 + 1;
                    new_args = realloc(new_args, sizeof(Expr*) * cap);
                }
                for (size_t j = 0; j < seq_count; j++) {
                    new_args[new_count++] = new_arg->data.function.args[j];
                    new_arg->data.function.args[j] = NULL;
                }
                expr_free(new_arg);
            } else {
                if (new_count >= cap) {
                    cap = cap * 2 + 1;
                    new_args = realloc(new_args, sizeof(Expr*) * cap);
                }
                new_args[new_count++] = new_arg;
            }
        }

        if (head_delete) {
            expr_free(new_head);
            result = expr_new_function(expr_new_symbol(SYM_Sequence), new_args, new_count);
        } else {
            result = expr_new_function(new_head, new_args, new_count);
        }
        free(new_args);
    } else {
        result = expr_copy(e);
    }

    /* Decide whether `e` should be removed from its parent. */
    if (*count_remaining != 0) {
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
                *delete_me = true;
                if (*count_remaining > 0) (*count_remaining)--;
            }
            env_free(env);
        }
    }

    return result;
}

Expr* builtin_delete_cases(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 1) {
        Expr* slot_args[1] = { expr_new_integer(1) };
        Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot), slot_args, 1);
        Expr* inner_args[2] = { slot, expr_copy(res->data.function.args[0]) };
        Expr* inner_dc = expr_new_function(expr_new_symbol(SYM_DeleteCases), inner_args, 2);
        Expr* func_args[1] = { inner_dc };
        return expr_new_function(expr_new_symbol(SYM_Function), func_args, 1);
    }

    /* DeleteCases[assoc, patt] drops entries whose value matches patt, keeping
     * the result an association. */
    if (argc == 2 && is_association(res->data.function.args[0]))
        return assoc_delete_cases(res->data.function.args[0], res->data.function.args[1]);

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
        } else if (ls->type == EXPR_FUNCTION && ls->data.function.head->type == EXPR_SYMBOL && ls->data.function.head->data.symbol == SYM_List) {
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
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->type == EXPR_SYMBOL && opt->data.function.head->data.symbol == SYM_Rule && opt->data.function.arg_count == 2) {
            if (opt->data.function.args[0]->type == EXPR_SYMBOL && opt->data.function.args[0]->data.symbol == SYM_Heads) {
                if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol == SYM_True) heads = true;
                else if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol == SYM_False) heads = false;
            }
        }
    }

    int64_t count_remaining = -1;
    if (argc >= 4) {
        Expr* n_expr = res->data.function.args[3];
        if (n_expr->type == EXPR_INTEGER && n_expr->data.integer >= 0) {
            count_remaining = n_expr->data.integer;
        }
    }

    bool dummy = false;
    Expr* result = do_delete_cases_at_level(expr, 0, min_l, max_l, heads, pattern, &count_remaining, &dummy);

    return result;
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
            Expr* pos_expr = expr_new_function(expr_new_symbol(SYM_List), path_exprs, path_len);
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
        Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot), slot_args, 1);
        Expr* inner_args[2] = { slot, expr_copy(res->data.function.args[0]) };
        Expr* inner_pos = expr_new_function(expr_new_symbol(SYM_Position), inner_args, 2);
        Expr* func_args[1] = { inner_pos };
        return expr_new_function(expr_new_symbol(SYM_Function), func_args, 1);
    }

    if (argc < 2) return NULL;

    /* Position[assoc, patt] reports positions inside the *values* as
     * {Key[k], subpos...} (Wolfram semantics). Delegate to Position over
     * Values[assoc] (which handles nested descent and head positions), then
     * remap each result's leading value-index to Key[key]. */
    if (argc == 2 && is_association(res->data.function.args[0])) {
        Expr* assoc = res->data.function.args[0];
        size_t na = assoc->data.function.arg_count;
        Expr* vals = assoc_values_list(assoc);                     /* List of values */
        Expr* pargs[2] = { vals, expr_copy(res->data.function.args[1]) };
        Expr* pcall = expr_new_function(expr_new_symbol(SYM_Position), pargs, 2);
        Expr* raw = evaluate(pcall);                               /* {{i, sub...}, ...} */
        expr_free(pcall);
        if (!(raw && raw->type == EXPR_FUNCTION &&
              raw->data.function.head->data.symbol == SYM_List)) return raw;

        size_t nr = raw->data.function.arg_count;
        Expr** out = malloc(sizeof(Expr*) * (nr ? nr : 1));
        size_t nout = 0;
        for (size_t i = 0; i < nr; i++) {
            Expr* pl = raw->data.function.args[i];
            /* Position lists look like {idx, sub...}; idx is the 1-based value
             * index. idx == 0 is the synthetic Values-list head — not a real
             * entry — so drop it. */
            if (!(pl->type == EXPR_FUNCTION && pl->data.function.head->data.symbol == SYM_List &&
                  pl->data.function.arg_count >= 1 &&
                  pl->data.function.args[0]->type == EXPR_INTEGER)) continue;
            int64_t vi = pl->data.function.args[0]->data.integer;
            if (vi < 1 || vi > (int64_t)na) continue;
            Expr* rule = assoc->data.function.args[vi - 1];
            if (!(rule->type == EXPR_FUNCTION && rule->data.function.arg_count == 2)) continue;
            size_t plen = pl->data.function.arg_count;
            Expr** np = malloc(sizeof(Expr*) * plen);
            Expr* kargs[1] = { expr_copy(rule->data.function.args[0]) };
            np[0] = expr_new_function(expr_new_symbol(SYM_Key), kargs, 1);
            for (size_t j = 1; j < plen; j++) np[j] = expr_copy(pl->data.function.args[j]);
            out[nout++] = expr_new_function(expr_new_symbol(SYM_List), np, plen);
            free(np);
        }
        expr_free(raw);
        Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, nout);
        free(out);
        return result;
    }

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

    Expr* list = expr_new_function(expr_new_symbol(SYM_List), results, count);
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
        Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot), slot_args, 1);
        Expr* inner_args[2] = { slot, expr_copy(res->data.function.args[0]) };
        Expr* inner_count = expr_new_function(expr_new_symbol(SYM_Count), inner_args, 2);
        Expr* func_args[1] = { inner_count };
        return expr_new_function(expr_new_symbol(SYM_Function), func_args, 1);
    }

    /* Count[assoc, patt] counts matching values (Count[Values[assoc], patt]). */
    if (argc >= 2 && is_association(res->data.function.args[0])) {
        Expr* r = assoc_apply_over_values(res); if (r) return r;
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
        Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot), slot_args, 1);
        Expr* inner_args[2] = { slot, expr_copy(res->data.function.args[0]) };
        Expr* inner_memberq = expr_new_function(expr_new_symbol(SYM_MemberQ), inner_args, 2);
        Expr* func_args[1] = { inner_memberq };
        return expr_new_function(expr_new_symbol(SYM_Function), func_args, 1);
    }

    /* MemberQ[assoc, form] tests the association's values. */
    if (argc >= 2 && is_association(res->data.function.args[0])) {
        Expr* r = assoc_apply_over_values(res); if (r) return r;
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
        return expr_new_symbol(SYM_True);
    } else {
        return expr_new_symbol(SYM_False);
    }
}

void patterns_init(void) {
    symtab_add_builtin("Cases", builtin_cases);
    symtab_get_def("Cases")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("FirstCase", builtin_first_case);
    symtab_get_def("FirstCase")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FirstCase",
        "FirstCase[expr, patt]\n\tGives the first element of expr matching patt,\n"
        "\tor Missing[\"NotFound\"]. FirstCase[expr, patt, default] uses default.\n"
        "\tOver an association, matches values and returns the first match.");
    symtab_add_builtin("DeleteCases", builtin_delete_cases);
    symtab_get_def("DeleteCases")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("DeleteMissing", builtin_delete_missing);
    symtab_get_def("DeleteMissing")->attributes |= ATTR_PROTECTED;
    /* KeyValuePattern is an inert pattern head (handled by the matcher, no
     * builtin); mark it Protected and give it a docstring. */
    symtab_get_def("KeyValuePattern")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("KeyValuePattern",
        "KeyValuePattern[{k1 -> p1, ...}]\n\tA pattern matching an association (or\n"
        "\tlist of rules) that contains keys matching k1, ... with values matching\n"
        "\tp1, .... Value patterns may bind (e.g. KeyValuePattern[{\"a\" -> v_}]).\n"
        "\tKeyValuePattern[k -> p] is the single-key form.");
    symtab_set_docstring("DeleteMissing",
        "DeleteMissing[expr]\n\tRemoves all Missing[...] elements (equivalent to\n"
        "\tDeleteCases[expr, _Missing]). Over an association, drops entries whose\n"
        "\tvalue is Missing[...].");
    symtab_add_builtin("Position", builtin_position);
    symtab_get_def("Position")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Count", builtin_count);
    symtab_get_def("Count")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("MemberQ", builtin_memberq);
    symtab_get_def("MemberQ")->attributes |= ATTR_PROTECTED;
}