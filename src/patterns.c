#include "patterns.h"
#include "symtab.h"
#include "eval.h"
#include "match.h"
#include "core.h"
#include "sym_names.h"
#include "assoc.h"
#include "assoc_struct.h"  /* associations are atomic: walkers see values only */
#include "ndarray.h"   /* is_ndarray — see patterns_delist_visible */
#include <stdlib.h>
#include <string.h>

/* A VISIBLE NDArray[...] is an ATOM to everything in this file.
 *
 * The matcher walks `data.function.args`, and an EXPR_NDARRAY has none — so
 * every head here searched an expression with no elements and answered, with
 * complete confidence, that it found nothing:
 *
 *     MemberQ[NDArray[Range[1., 300.]], 5.]   ->  False   (List: True)
 *     Count[NDArray[Range[300]], 5]           ->  0       (List: 1)
 *     Position[NDArray[Range[300]], 5]        ->  {}      (List: {{5}})
 *     Cases[NDArray[Range[300]], 5]           ->  {}      (List: {5})
 *
 * That is worse than declining: a wrong answer that looks like a right one.
 *
 * The PACKED representation was never affected, because these heads are not on
 * pack.c's AWARE list and the transparency gate materialises their arguments
 * before they run. This is the same asymmetry that let every real kernel
 * truncate a visible int64 buffer (ndarray_map_unary): the gate protects one
 * representation and the other has to be handled where it arrives.
 *
 * Materialising is the right answer rather than a buffer-aware matcher: a
 * pattern here is an arbitrary expression, so there is no fast path to reach —
 * only a correct answer to produce. By the time a builtin runs, a packed
 * argument has already been materialised by the gate, so this fires only for a
 * genuinely visible array and costs nothing otherwise.
 *
 * ONLY args[0], the expression being searched -- never the PATTERN. The two
 * arguments are not the same kind of thing: materialising the searched
 * expression preserves the answer by the transparency contract, while
 * materialising a pattern would change what it matches, since a visible
 * NDArray[...] and the plain List of the same values are distinct expressions.
 * `MemberQ[list, NDArray[{1., 2.}]]` asks for elements equal to that array and
 * must keep asking for exactly that. The single-argument OPERATOR forms
 * (`Cases[patt]`, `MemberQ[patt]`) are skipped for the same reason: their one
 * argument is the pattern.
 *
 * Rebuilt by hand rather than through ndarray_delist_and_reeval, which
 * materialises EVERY array argument -- the pattern included, which is precisely
 * what must not happen here.
 *
 * Returns the re-evaluated call, or NULL when args[0] is not an array. */
static Expr* patterns_delist_visible(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n < 2) return NULL;                                 /* operator form */
    if (!is_ndarray(res->data.function.args[0])) return NULL;

    Expr** args = malloc(sizeof(Expr*) * n);
    if (!args) return NULL;
    args[0] = ndarray_to_nested_list(res->data.function.args[0]);
    for (size_t i = 1; i < n; i++)
        args[i] = expr_copy(res->data.function.args[i]);
    Expr* rebuilt = expr_new_function(expr_copy(res->data.function.head), args, n);
    free(args);
    return eval_and_free(rebuilt);
}

static int64_t get_expr_depth_patterns(Expr* e, bool heads) {
    if (e->type != EXPR_FUNCTION) return 1;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_Rational || h == SYM_Complex) return 1;
    }
    bool assoc = assoc_is_wellformed(e);     /* parts = values */
    int64_t max_d = assoc ? 1 : 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        int64_t d = get_expr_depth_patterns(struct_part(e, i, assoc), heads);
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
        bool assoc = assoc_is_wellformed(e);
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (max_results >= 0 && (int64_t)(*count) >= max_results) break;
            do_cases_at_level(struct_part(e, i, assoc), current_level + 1, min_l, max_l, heads, pattern, replacement, delayed, results, count, cap, max_results);
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
    /* A visible NDArray is an atom to the matcher; materialise it first.
     * See patterns_delist_visible. */
    { Expr* nd_ = patterns_delist_visible(res); if (nd_) return nd_; }
    size_t argc = res->data.function.arg_count;
    if (argc != 2 && argc != 3) return NULL;

    Expr* cases_args[2] = { expr_copy(res->data.function.args[0]),
                            expr_copy(res->data.function.args[1]) };
    Expr* cases = expr_new_function(expr_new_symbol(SYM_Cases), cases_args, 2);
    Expr* matches = evaluate(cases);
    expr_free(cases);

    Expr* result = NULL;
    if (matches && matches->type == EXPR_FUNCTION &&
        matches->data.function.head->data.symbol.name == SYM_List &&
        matches->data.function.arg_count >= 1) {
        result = expr_copy(matches->data.function.args[0]);
    }
    if (matches) expr_free(matches);
    if (result) return result;

    if (argc == 3) return expr_copy(res->data.function.args[2]);
    Expr* margs[1] = { expr_new_string("NotFound") };
    return expr_new_function(expr_new_symbol(SYM_Missing), margs, 1);
}

/* FirstPosition[expr, pattern] gives the position (a list of indices) of the
 * first element of expr matching pattern in depth-first order, or
 * Missing["NotFound"] if none is found.
 * FirstPosition[expr, pattern, default] returns default instead of the Missing.
 * default is held (FirstPosition is HoldRest) and only evaluated when returned.
 * FirstPosition[expr, pattern, default, levelspec] searches only the specified
 * levels; a Heads -> True|False option is honoured (default True).
 *
 * This delegates to Position and takes the first result, so it inherits
 * Position's levelspec parsing, Heads handling, association value -> Key[...]
 * remapping, depth-first ordering, and first-n-match early stop unchanged. */
Expr* builtin_first_position(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    /* A visible NDArray is an atom to the matcher; materialise it first.
     * See patterns_delist_visible. */
    { Expr* nd_ = patterns_delist_visible(res); if (nd_) return nd_; }
    size_t argc = res->data.function.arg_count;

    /* Split the (held) trailing Heads -> True|False option from the positional
     * arguments expr, pattern, default, levelspec (default before levelspec). */
    Expr* heads_opt = NULL;              /* borrowed pointer into res */
    Expr* pos_args[4];
    size_t np = 0;
    for (size_t i = 0; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (a->type == EXPR_FUNCTION && a->data.function.arg_count == 2 &&
            a->data.function.head->type == EXPR_SYMBOL &&
            a->data.function.head->data.symbol.name == SYM_Rule &&
            a->data.function.args[0]->type == EXPR_SYMBOL &&
            a->data.function.args[0]->data.symbol.name == SYM_Heads) {
            heads_opt = a;
        } else if (np < 4) {
            pos_args[np++] = a;
        } else {
            return NULL;                 /* too many positional arguments */
        }
    }
    if (np < 2) return NULL;             /* need at least expr and pattern */

    Expr* expr    = pos_args[0];
    Expr* pattern = pos_args[1];
    Expr* deflt   = (np >= 3) ? pos_args[2] : NULL;   /* held; only if returned */
    Expr* lspec   = (np >= 4) ? pos_args[3] : NULL;

    /* Build the delegated Position call. */
    Expr* pcall;
    if (is_association(expr) && lspec == NULL) {
        /* Only the 2-arg Position form performs the value -> Key[k] remapping. */
        Expr* pa[2] = { expr_copy(expr), expr_copy(pattern) };
        pcall = expr_new_function(expr_new_symbol(SYM_Position), pa, 2);
    } else {
        /* levelspec: the supplied one, or the {0, Infinity} default (which
         * includes level 0 == the whole expression, giving position {}). */
        Expr* ls;
        if (lspec) {
            ls = expr_copy(lspec);
        } else {
            Expr* largs[2] = { expr_new_integer(0), expr_new_symbol(SYM_Infinity) };
            ls = expr_new_function(expr_new_symbol(SYM_List), largs, 2);
        }
        Expr* pa[5];
        size_t k = 0;
        pa[k++] = expr_copy(expr);
        pa[k++] = expr_copy(pattern);
        pa[k++] = ls;
        pa[k++] = expr_new_integer(1);   /* n = 1: stop at the first match */
        if (heads_opt) pa[k++] = expr_copy(heads_opt);
        pcall = expr_new_function(expr_new_symbol(SYM_Position), pa, k);
    }

    Expr* raw = evaluate(pcall);
    expr_free(pcall);

    Expr* result = NULL;
    if (raw && raw->type == EXPR_FUNCTION &&
        raw->data.function.head->data.symbol.name == SYM_List &&
        raw->data.function.arg_count >= 1) {
        result = expr_copy(raw->data.function.args[0]);
    }
    if (raw) expr_free(raw);
    if (result) return result;

    /* No match: the (held) default, evaluated by the fixed-point loop when it is
     * returned, else Missing["NotFound"]. */
    if (deflt) return expr_copy(deflt);
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
    /* A visible NDArray is an atom to the matcher; materialise it first.
     * See patterns_delist_visible. */
    { Expr* nd_ = patterns_delist_visible(res); if (nd_) return nd_; }
    size_t argc = res->data.function.arg_count;
    
    if (argc == 1) {
        Expr* slot_args[1] = { expr_new_integer(1) };
        Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot), slot_args, 1);
        Expr* inner_args[2] = { slot, expr_copy(res->data.function.args[0]) };
        Expr* inner_cases = expr_new_function(expr_new_symbol(SYM_Cases), inner_args, 2);
        Expr* func_args[1] = { inner_cases };
        return expr_new_function(expr_new_symbol(SYM_Function), func_args, 1);
    }

    /* Cases[assoc, patt] collects matching values (Cases[Values[assoc], patt]).
     * With a level spec or option the association-aware walker below runs, so
     * that e.g. Cases[<|a -> 1|>, _, {0}] is {<|a -> 1|>} and not {{1}}. */
    if (argc == 2 && assoc_is_wellformed(res->data.function.args[0])) {
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
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol.name == SYM_All) {
            min_l = 1; max_l = 1000000;
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol.name == SYM_Infinity) {
            min_l = 1; max_l = 1000000;
        } else if (ls->type == EXPR_FUNCTION && ls->data.function.head->data.symbol.name == SYM_List) {
            if (ls->data.function.arg_count == 1 && ls->data.function.args[0]->type == EXPR_INTEGER) {
                min_l = max_l = ls->data.function.args[0]->data.integer;
            } else if (ls->data.function.arg_count == 2) {
                if (ls->data.function.args[0]->type == EXPR_INTEGER) min_l = ls->data.function.args[0]->data.integer;
                if (ls->data.function.args[1]->type == EXPR_INTEGER) max_l = ls->data.function.args[1]->data.integer;
                else if (ls->data.function.args[1]->type == EXPR_SYMBOL && ls->data.function.args[1]->data.symbol.name == SYM_Infinity) max_l = 1000000;
            }
        }
    }

    for (size_t i = 2; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->data.symbol.name == SYM_Rule && opt->data.function.arg_count == 2) {
            if (opt->data.function.args[0]->type == EXPR_SYMBOL && opt->data.function.args[0]->data.symbol.name == SYM_Heads) {
                if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol.name == SYM_True) heads = true;
                else if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol.name == SYM_False) heads = false;
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

    if (patt_arg->type == EXPR_FUNCTION && (patt_arg->data.function.head->data.symbol.name == SYM_Rule || patt_arg->data.function.head->data.symbol.name == SYM_RuleDelayed) && patt_arg->data.function.arg_count == 2) {
        pattern = patt_arg->data.function.args[0];
        replacement = patt_arg->data.function.args[1];
        delayed = (patt_arg->data.function.head->data.symbol.name == SYM_RuleDelayed);
    } else if (patt_arg->type == EXPR_FUNCTION && patt_arg->data.function.head->data.symbol.name == SYM_HoldPattern && patt_arg->data.function.arg_count == 1) {
        Expr* hp_arg = patt_arg->data.function.args[0];
        if (hp_arg->type == EXPR_FUNCTION && (hp_arg->data.function.head->data.symbol.name == SYM_Rule || hp_arg->data.function.head->data.symbol.name == SYM_RuleDelayed) && hp_arg->data.function.arg_count == 2) {
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

    if (assoc_is_wellformed(e)) {
        /* An atomic association: its parts are its values, and deleting a
         * value drops the whole entry, so the result stays well formed:
         * DeleteCases[{<|a -> 1, b -> 2|>}, 1, Infinity] is {<|b -> 2|>}.
         * Keys are never tested (assoc_struct.h). */
        size_t n = e->data.function.arg_count;
        Expr** kept = malloc(sizeof(Expr*) * (n ? n : 1));
        size_t nk = 0;
        bool head_delete = false;
        Expr* new_head = heads
            ? do_delete_cases_at_level(e->data.function.head, current_level + 1, min_l, max_l, heads, pattern, count_remaining, &head_delete)
            : expr_copy(e->data.function.head);
        for (size_t i = 0; i < n; i++) {
            Expr* entry = e->data.function.args[i];
            bool del = false;
            Expr* nv = do_delete_cases_at_level(struct_part(e, i, true), current_level + 1, min_l, max_l, heads, pattern, count_remaining, &del);
            if (del) { expr_free(nv); continue; }
            Expr* rargs[2] = { expr_copy(entry->data.function.args[0]), nv };
            kept[nk++] = expr_new_function(expr_copy(entry->data.function.head), rargs, 2);
        }
        if (head_delete) {
            expr_free(new_head);
            new_head = expr_new_symbol(SYM_Sequence);
        }
        result = expr_new_function(new_head, kept, nk);
        free(kept);
    } else if (e->type == EXPR_FUNCTION) {
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
                new_arg->data.function.head->data.symbol.name == SYM_Sequence) {
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
    /* A visible NDArray is an atom to the matcher; materialise it first.
     * See patterns_delist_visible. */
    { Expr* nd_ = patterns_delist_visible(res); if (nd_) return nd_; }
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
    if (argc == 2 && assoc_is_wellformed(res->data.function.args[0]))
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
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol.name == SYM_All) {
            min_l = 1; max_l = 1000000;
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol.name == SYM_Infinity) {
            min_l = 1; max_l = 1000000;
        } else if (ls->type == EXPR_FUNCTION && ls->data.function.head->type == EXPR_SYMBOL && ls->data.function.head->data.symbol.name == SYM_List) {
            if (ls->data.function.arg_count == 1 && ls->data.function.args[0]->type == EXPR_INTEGER) {
                min_l = max_l = ls->data.function.args[0]->data.integer;
            } else if (ls->data.function.arg_count == 2) {
                if (ls->data.function.args[0]->type == EXPR_INTEGER) min_l = ls->data.function.args[0]->data.integer;
                if (ls->data.function.args[1]->type == EXPR_INTEGER) max_l = ls->data.function.args[1]->data.integer;
                else if (ls->data.function.args[1]->type == EXPR_SYMBOL && ls->data.function.args[1]->data.symbol.name == SYM_Infinity) max_l = 1000000;
            }
        }
    }

    for (size_t i = 2; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->type == EXPR_SYMBOL && opt->data.function.head->data.symbol.name == SYM_Rule && opt->data.function.arg_count == 2) {
            if (opt->data.function.args[0]->type == EXPR_SYMBOL && opt->data.function.args[0]->data.symbol.name == SYM_Heads) {
                if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol.name == SYM_True) heads = true;
                else if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol.name == SYM_False) heads = false;
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

/* `current_path` holds the integer position components; `current_keys` runs
 * parallel to it and is non-NULL at a component that lies inside an atomic
 * association, where the component is Key[current_keys[i]] (a borrowed key)
 * instead of an integer: Position[{<|a -> 1|>}, 1] is {{1, Key[a]}}. */
static void do_position_at_level(Expr* e, int64_t current_level, int64_t min_l, int64_t max_l, bool heads, Expr* pattern, Expr*** results, size_t* count, size_t* cap, int64_t max_results, int64_t* current_path, Expr** current_keys, size_t path_len) {
    if (max_results >= 0 && (int64_t)(*count) >= max_results) return;

    if (e->type == EXPR_FUNCTION) {
        int64_t* next_path = malloc(sizeof(int64_t) * (path_len + 1));
        Expr** next_keys = malloc(sizeof(Expr*) * (path_len + 1));
        if (path_len > 0) {
            memcpy(next_path, current_path, sizeof(int64_t) * path_len);
            memcpy(next_keys, current_keys, sizeof(Expr*) * path_len);
        }
        bool assoc = assoc_is_wellformed(e);   /* parts = values, at Key[k] */

        if (heads) {
            next_path[path_len] = 0;
            next_keys[path_len] = NULL;
            do_position_at_level(e->data.function.head, current_level + 1, min_l, max_l, heads, pattern, results, count, cap, max_results, next_path, next_keys, path_len + 1);
        }
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (max_results >= 0 && (int64_t)(*count) >= max_results) break;
            next_path[path_len] = i + 1;
            next_keys[path_len] = assoc ? e->data.function.args[i]->data.function.args[0] : NULL;
            do_position_at_level(struct_part(e, i, assoc), current_level + 1, min_l, max_l, heads, pattern, results, count, cap, max_results, next_path, next_keys, path_len + 1);
        }
        free(next_path);
        free(next_keys);
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
            Expr** path_exprs = malloc(sizeof(Expr*) * (path_len ? path_len : 1));
            for (size_t i = 0; i < path_len; i++) {
                if (current_keys[i]) {
                    Expr* k = expr_copy(current_keys[i]);
                    path_exprs[i] = expr_new_function(expr_new_symbol(SYM_Key), &k, 1);
                } else {
                    path_exprs[i] = expr_new_integer(current_path[i]);
                }
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
    /* A visible NDArray is an atom to the matcher; materialise it first.
     * See patterns_delist_visible. */
    { Expr* nd_ = patterns_delist_visible(res); if (nd_) return nd_; }
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

    /* Position inside an association reports {Key[k], subpos...} for its
     * values (and {0} for its head, {} for itself) -- handled by the
     * association-aware walker, at any nesting depth (assoc_struct.h). */

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
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol.name == SYM_All) {
            min_l = 1; max_l = 1000000;
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol.name == SYM_Infinity) {
            min_l = 1; max_l = 1000000;
        } else if (ls->type == EXPR_FUNCTION && ls->data.function.head->data.symbol.name == SYM_List) {
            if (ls->data.function.arg_count == 1 && ls->data.function.args[0]->type == EXPR_INTEGER) {
                min_l = max_l = ls->data.function.args[0]->data.integer;
            } else if (ls->data.function.arg_count == 2) {
                if (ls->data.function.args[0]->type == EXPR_INTEGER) min_l = ls->data.function.args[0]->data.integer;
                if (ls->data.function.args[1]->type == EXPR_INTEGER) max_l = ls->data.function.args[1]->data.integer;
                else if (ls->data.function.args[1]->type == EXPR_SYMBOL && ls->data.function.args[1]->data.symbol.name == SYM_Infinity) max_l = 1000000;
            }
        }
    }

    for (size_t i = 2; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->data.symbol.name == SYM_Rule && opt->data.function.arg_count == 2) {
            if (opt->data.function.args[0]->type == EXPR_SYMBOL && opt->data.function.args[0]->data.symbol.name == SYM_Heads) {
                if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol.name == SYM_True) heads = true;
                else if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol.name == SYM_False) heads = false;
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

    do_position_at_level(expr, 0, min_l, max_l, heads, pattern, &results, &count, &cap, max_results, NULL, NULL, 0);

    Expr* list = expr_new_function(expr_new_symbol(SYM_List), results, count);
    free(results);
    return list;
}

static void do_count_at_level(Expr* e, int64_t current_level, int64_t min_l, int64_t max_l, bool heads, Expr* pattern, size_t* count) {
    if (e->type == EXPR_FUNCTION) {
        if (heads) {
            do_count_at_level(e->data.function.head, current_level + 1, min_l, max_l, heads, pattern, count);
        }
        bool assoc = assoc_is_wellformed(e);   /* parts = values */
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            do_count_at_level(struct_part(e, i, assoc), current_level + 1, min_l, max_l, heads, pattern, count);
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
    /* A visible NDArray is an atom to the matcher; materialise it first.
     * See patterns_delist_visible. */
    { Expr* nd_ = patterns_delist_visible(res); if (nd_) return nd_; }
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
    if (argc == 2 && assoc_is_wellformed(res->data.function.args[0])) {
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
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol.name == SYM_All) {
            min_l = 1; max_l = 1000000;
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol.name == SYM_Infinity) {
            min_l = 1; max_l = 1000000;
        } else if (ls->type == EXPR_FUNCTION && ls->data.function.head->data.symbol.name == SYM_List) {
            if (ls->data.function.arg_count == 1 && ls->data.function.args[0]->type == EXPR_INTEGER) {
                min_l = max_l = ls->data.function.args[0]->data.integer;
            } else if (ls->data.function.arg_count == 2) {
                if (ls->data.function.args[0]->type == EXPR_INTEGER) min_l = ls->data.function.args[0]->data.integer;
                if (ls->data.function.args[1]->type == EXPR_INTEGER) max_l = ls->data.function.args[1]->data.integer;
                else if (ls->data.function.args[1]->type == EXPR_SYMBOL && ls->data.function.args[1]->data.symbol.name == SYM_Infinity) max_l = 1000000;
            }
        }
    }

    for (size_t i = 2; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->data.symbol.name == SYM_Rule && opt->data.function.arg_count == 2) {
            if (opt->data.function.args[0]->type == EXPR_SYMBOL && opt->data.function.args[0]->data.symbol.name == SYM_Heads) {
                if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol.name == SYM_True) heads = true;
                else if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol.name == SYM_False) heads = false;
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
        bool assoc = assoc_is_wellformed(e);   /* parts = values */
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (do_member_at_level(struct_part(e, i, assoc), current_level + 1, min_l, max_l, heads, pattern)) return true;
        }
    }
    return false;
}

Expr* builtin_memberq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    /* A visible NDArray is an atom to the matcher; materialise it first.
     * See patterns_delist_visible. */
    { Expr* nd_ = patterns_delist_visible(res); if (nd_) return nd_; }
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
    if (argc == 2 && assoc_is_wellformed(res->data.function.args[0])) {
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
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol.name == SYM_All) {
            min_l = 1; max_l = 1000000;
        } else if (ls->type == EXPR_SYMBOL && ls->data.symbol.name == SYM_Infinity) {
            min_l = 1; max_l = 1000000;
        } else if (ls->type == EXPR_FUNCTION && ls->data.function.head->data.symbol.name == SYM_List) {
            if (ls->data.function.arg_count == 1 && ls->data.function.args[0]->type == EXPR_INTEGER) {
                min_l = max_l = ls->data.function.args[0]->data.integer;
            } else if (ls->data.function.arg_count == 2) {
                if (ls->data.function.args[0]->type == EXPR_INTEGER) min_l = ls->data.function.args[0]->data.integer;
                if (ls->data.function.args[1]->type == EXPR_INTEGER) max_l = ls->data.function.args[1]->data.integer;
                else if (ls->data.function.args[1]->type == EXPR_SYMBOL && ls->data.function.args[1]->data.symbol.name == SYM_Infinity) max_l = 1000000;
            }
        }
    }

    for (size_t i = 2; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->data.symbol.name == SYM_Rule && opt->data.function.arg_count == 2) {
            if (opt->data.function.args[0]->type == EXPR_SYMBOL && opt->data.function.args[0]->data.symbol.name == SYM_Heads) {
                if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol.name == SYM_True) heads = true;
                else if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol.name == SYM_False) heads = false;
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
    symtab_add_builtin("FirstPosition", builtin_first_position);
    symtab_get_def("FirstPosition")->attributes |= ATTR_HOLDREST | ATTR_PROTECTED;
    symtab_add_builtin("Count", builtin_count);
    symtab_get_def("Count")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("MemberQ", builtin_memberq);
    symtab_get_def("MemberQ")->attributes |= ATTR_PROTECTED;
}