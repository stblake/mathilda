/*
 * eval.c
 *
 * This file implements the core evaluation engine of PicoCAS.
 * It follows the "infinite evaluation" semantics of the Mathematica:
 * expressions are repeatedly transformed until they no longer change.
 *
 * The main entry point is evaluate(), which calls evaluate_step() in a loop.
 */

#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "purefunc.h"
#include "print.h"
#include "deriv.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/*
 * The maximum number of evaluation steps to prevent infinite recursion
 * in cases of circular definitions.
 */
#define MAX_ITERATIONS 4096

/*
 * $RecursionLimit guard. The REPL is single-threaded so a static counter
 * suffices. Each call to evaluate() bumps eval_recursion_depth on entry and
 * decrements on exit; if the depth would exceed eval_recursion_limit we
 * return the expression wrapped in Hold[] (so it cannot re-enter the
 * evaluator) and emit a $RecursionLimit::reclim message.
 *
 * The default of 512 leaves comfortable headroom under the typical 8 MB
 * thread stack while still catching pathological recursion (e.g. a self-
 * referential rule like f[x_] := f[x] + 1).
 */
#define DEFAULT_RECURSION_LIMIT 512
static int  eval_recursion_depth = 0;
static int  eval_recursion_limit = DEFAULT_RECURSION_LIMIT;
/* Once the recursion limit is hit anywhere in the evaluation tree, this
 * sticky flag tells *every* enclosing evaluate() loop to stop iterating
 * and return its current value. Without it, the wrap-in-Hold result
 * causes outer fixed-point loops to perceive endless "progress" (each
 * iteration adds another Hold wrapper) and chew through all 4096
 * outer iterations at every level of the unwind. The flag is cleared
 * at the top of each top-level evaluate() call. */
static bool eval_overflow = false;

int  eval_get_recursion_limit(void) { return eval_recursion_limit; }
int  eval_get_recursion_depth(void) { return eval_recursion_depth; }
void eval_set_recursion_limit(int n) {
    eval_recursion_limit = (n > 0) ? n : DEFAULT_RECURSION_LIMIT;
}

/*
 * eval_compare_expr_ptrs:
 * Helper for sorting expression arguments when a head has the Orderless attribute.
 * Uses the canonical expr_compare to ensure a stable, deterministic order.
 */
int eval_compare_expr_ptrs(const void* a, const void* b) {
    Expr* ea = *(Expr**)a;
    Expr* eb = *(Expr**)b;
    return expr_compare(ea, eb);
}

/*
 * flatten_args:
 * Implements the Flat (associative) attribute.
 * If a function has the same head as some of its arguments, those arguments
 * are "unwrapped" and their elements are promoted to be direct arguments 
 * of the parent function.
 * Example: f[a, f[b, c], d] -> f[a, b, c, d]
 */
void eval_flatten_args(Expr* e, const char* head_name) {
    size_t new_count = 0;
    bool needs_flattening = false;
    
    /* First pass: calculate the total number of arguments after flattening */
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && 
            strcmp(arg->data.function.head->data.symbol, head_name) == 0) {
            new_count += arg->data.function.arg_count;
            needs_flattening = true;
        } else {
            new_count++;
        }
    }
    
    /* If no nested occurrences of the head were found, we are done */
    if (!needs_flattening) return;
    
    /* Second pass: allocate new argument array and copy elements */
    Expr** new_args = malloc(sizeof(Expr*) * new_count);
    size_t idx = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && 
            strcmp(arg->data.function.head->data.symbol, head_name) == 0) {
            /* Splat nested arguments into the new array */
            for (size_t j = 0; j < arg->data.function.arg_count; j++) {
                new_args[idx++] = expr_copy(arg->data.function.args[j]);
            }
            /* Free the intermediate nested function node */
            expr_free(arg); 
        } else {
            new_args[idx++] = arg;
        }
    }
    
    /* Replace old argument array with the flattened one */
    free(e->data.function.args);
    e->data.function.args = new_args;
    e->data.function.arg_count = new_count;
}

/*
 * has_list_arg:
 * Helper for the Listable attribute.
 * Checks if any argument of the function is an explicit List[...].
 */
static bool has_list_arg(Expr* e) {
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && strcmp(arg->data.function.head->data.symbol, "List") == 0) {
            return true;
        }
    }
    return false;
}

/*
 * apply_listable:
 * Implements automatic threading of functions over lists.
 * If a function is Listable and contains list arguments, it maps itself over them.
 * Example: f[{a, b}, c] -> {f[a, c], f[b, c]}
 */
static Expr* apply_listable(Expr* e) {
    /* Determine the required length of the result list */
    size_t list_len = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && strcmp(arg->data.function.head->data.symbol, "List") == 0) {
            list_len = arg->data.function.arg_count;
            break;
        }
    }
    
    if (list_len == 0) return NULL;
    
    /* Construct a new List containing the threaded evaluations */
    Expr** new_list_args = malloc(sizeof(Expr*) * list_len);
    for (size_t j = 0; j < list_len; j++) {
        Expr** new_func_args = malloc(sizeof(Expr*) * e->data.function.arg_count);
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            Expr* arg = e->data.function.args[i];
            if (arg->type == EXPR_FUNCTION && strcmp(arg->data.function.head->data.symbol, "List") == 0) {
                /* All list arguments must have identical lengths */
                if (arg->data.function.arg_count != list_len) {
                    char* s = expr_to_string(e);
                    printf("Thread::tdlen: Objects of unequal length in %s cannot be combined.\n", s);
                    free(s);
                    for (size_t k = 0; k < i; k++) expr_free(new_func_args[k]);
                    free(new_func_args);
                    for (size_t k = 0; k < j; k++) expr_free(new_list_args[k]);
                    free(new_list_args);
                    return NULL;
                }
                new_func_args[i] = expr_copy(arg->data.function.args[j]);
            } else {
                /* Non-list arguments are repeated for every element */
                new_func_args[i] = expr_copy(arg);
            }
        }
        /* Recursively evaluate each threaded call */
        Expr* tmp = expr_new_function(expr_copy(e->data.function.head), new_func_args, e->data.function.arg_count);
        free(new_func_args);
        new_list_args[j] = evaluate(tmp);
        expr_free(tmp);
    }
    
    Expr* final_res = expr_new_function(expr_new_symbol("List"), new_list_args, list_len);
    free(new_list_args);
    return final_res;
}

/*
 * assignment_target_symbol:
 * Returns the name of the symbol whose OwnValue/DownValue/Part would be
 * written if lhs were used as the LHS of a Set or SetDelayed, or NULL if
 * lhs does not name a specific symbol. Wrappers that do not change the
 * ultimate target (Condition[pat, test], HoldPattern[pat], Part[x, ...])
 * are unwrapped so that Protected can be detected on the underlying head.
 */
static const char* assignment_target_symbol(Expr* lhs) {
    if (!lhs) return NULL;
    if (lhs->type == EXPR_SYMBOL) return lhs->data.symbol;
    if (lhs->type == EXPR_FUNCTION &&
        lhs->data.function.head->type == EXPR_SYMBOL &&
        lhs->data.function.arg_count >= 1) {
        const char* h = lhs->data.function.head->data.symbol;
        if (strcmp(h, "Condition") == 0 ||
            strcmp(h, "HoldPattern") == 0 ||
            strcmp(h, "Part") == 0) {
            return assignment_target_symbol(lhs->data.function.args[0]);
        }
        return h;
    }
    return NULL;
}

/*
 * apply_assignment:
 * Helper to handle the 'Set' (=) and 'SetDelayed' (:=) primitives.
 * Supports recursive list destructuring.
 * Example: {x, y} = {1, 2}
 *
 * Returns true if the caller should respond as if the assignment
 * succeeded (return the RHS for Set, Null for SetDelayed). Attempts to
 * assign to a Protected symbol emit a Set::wrsym message, leave state
 * unchanged, and still return true so the caller yields the RHS, matching
 * Mathematica semantics.
 */
static bool apply_assignment(Expr* lhs, Expr* rhs, bool is_delayed) {
    /* Block writes to Protected symbols. List destructuring is recursed
     * into below and each child runs through apply_assignment again, so
     * per-element protection checks happen naturally -- we only skip the
     * outer check when lhs itself is a List. */
    bool lhs_is_list = (lhs->type == EXPR_FUNCTION &&
                        lhs->data.function.head->type == EXPR_SYMBOL &&
                        strcmp(lhs->data.function.head->data.symbol, "List") == 0);
    if (!lhs_is_list) {
        const char* target = assignment_target_symbol(lhs);
        if (target && (get_attributes(target) & ATTR_PROTECTED)) {
            fprintf(stderr, "%s::wrsym: Symbol %s is Protected.\n",
                    is_delayed ? "SetDelayed" : "Set", target);
            return true;
        }
    }

    if (lhs->type == EXPR_SYMBOL) {
        /* Standard symbol assignment */
        symtab_add_own_value(lhs->data.symbol, lhs, rhs);
        return true;
    } else if (lhs->type == EXPR_FUNCTION) {
        if (lhs->data.function.head->type == EXPR_SYMBOL && 
            strcmp(lhs->data.function.head->data.symbol, "List") == 0 &&
            rhs->type == EXPR_FUNCTION &&
            rhs->data.function.head->type == EXPR_SYMBOL &&
            strcmp(rhs->data.function.head->data.symbol, "List") == 0) {
            
            /* List destructuring: match lengths and recurse. Any child that
             * cannot be assigned (e.g. LHS element is a literal number) fails
             * the whole destructuring so the caller does not misreport success. */
            if (lhs->data.function.arg_count != rhs->data.function.arg_count) {
                return false;
            }

            bool all_ok = true;
            for (size_t i = 0; i < lhs->data.function.arg_count; i++) {
                if (!apply_assignment(lhs->data.function.args[i], rhs->data.function.args[i], is_delayed)) {
                    all_ok = false;
                }
            }
            return all_ok;
        } else if (lhs->data.function.head->type == EXPR_SYMBOL && strcmp(lhs->data.function.head->data.symbol, "Part") == 0) {
            Expr* expr_part_assign(Expr* lhs, Expr* rhs); // Forward declare or include part.h
            Expr* assigned = expr_part_assign(lhs, rhs);
            if (assigned) {
                expr_free(assigned);
                return true;
            }
            return false;
        } else if (lhs->data.function.head->type == EXPR_SYMBOL) {
            /* Pattern-based assignment (DownValues) */
            /* We use the entire lhs as the pattern, and its head as the key */
            const char* symbol_name = lhs->data.function.head->data.symbol;

            Expr* actual_pattern = lhs;
            Expr* actual_rhs = rhs;

            /* If the RHS is Condition[body, test], move the condition to the LHS.
             * This makes f[x_] := body /; test equivalent to f[x_] /; test := body.
             * This is standard Mathematica semantics. */
            if (is_delayed && rhs->type == EXPR_FUNCTION &&
                rhs->data.function.head->type == EXPR_SYMBOL &&
                strcmp(rhs->data.function.head->data.symbol, "Condition") == 0 &&
                rhs->data.function.arg_count == 2) {
                /* Build Condition[lhs, test] as the new pattern */
                Expr** cond_args = malloc(sizeof(Expr*) * 2);
                cond_args[0] = expr_copy(lhs);
                cond_args[1] = expr_copy(rhs->data.function.args[1]);
                actual_pattern = expr_new_function(expr_new_symbol("Condition"), cond_args, 2);
                free(cond_args);
                /* The actual replacement is just the body */
                actual_rhs = rhs->data.function.args[0];
            }

            if (strcmp(symbol_name, "Condition") == 0 && actual_pattern->data.function.arg_count == 2) {
                Expr* inner_lhs = actual_pattern->data.function.args[0];
                if (inner_lhs->type == EXPR_FUNCTION && inner_lhs->data.function.head->type == EXPR_SYMBOL) {
                    symbol_name = inner_lhs->data.function.head->data.symbol;
                } else if (inner_lhs->type == EXPR_SYMBOL) {
                    symbol_name = inner_lhs->data.symbol;
                }
            }
            symtab_add_down_value(symbol_name, actual_pattern, actual_rhs);
            if (actual_pattern != lhs) expr_free(actual_pattern);
            return true;
        }
    }
    return false;
}

/*
 * flatten_sequences:
 * Flattens any Sequence[...] heads found in the arguments of e.
 */
static void flatten_sequences(Expr* e) {
    if (e->type != EXPR_FUNCTION) return;
    
    size_t new_count = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
            strcmp(arg->data.function.head->data.symbol, "Sequence") == 0) {
            new_count += arg->data.function.arg_count;
        } else {
            new_count++;
        }
    }
    
    if (new_count == e->data.function.arg_count) return;
    
    Expr** new_args = malloc(sizeof(Expr*) * new_count);
    size_t k = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
            strcmp(arg->data.function.head->data.symbol, "Sequence") == 0) {
            for (size_t j = 0; j < arg->data.function.arg_count; j++) {
                new_args[k++] = expr_copy(arg->data.function.args[j]);
            }
            expr_free(arg);
        } else {
            new_args[k++] = arg;
        }
    }
    
    free(e->data.function.args);
    e->data.function.args = new_args;
    e->data.function.arg_count = new_count;
}

/*
 * evaluate_step:
 * Performs exactly one level of evaluation transformation.
 */
Expr* evaluate_step(Expr* e) {
    if (!e) return NULL;

    switch (e->type) {        /* Atomics evaluate to themselves */
        case EXPR_INTEGER:
        case EXPR_REAL:
        case EXPR_STRING:
        case EXPR_BIGINT:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return expr_copy(e);
            
        case EXPR_SYMBOL: {
            /* Check for immediate assignments (OwnValues) like x = 5 */
            Expr* own = apply_own_values(e);
            if (own) return own;
            return expr_copy(e);
        }
            
        case EXPR_FUNCTION: {
            /* 1. Evaluate the head recursively (e.g. f[x][y]) */
            Expr* head = evaluate(e->data.function.head);
            
            uint32_t attrs = ATTR_NONE;
            if (head->type == EXPR_SYMBOL) {
                attrs = get_attributes(head->data.symbol);
            } else if (head->type == EXPR_FUNCTION) {
                attrs = pure_function_attributes(head);
            }
            
            /* 2. Handle 'Hold' attributes.
             *
             * HoldAllComplete is like HoldAll but additionally suppresses
             * Sequence flattening, Unevaluated stripping, Flat flattening,
             * and (eventually) UpValues lookup. Inside a HoldAllComplete
             * head, Evaluate[expr] does NOT force evaluation. Built-ins and
             * DownValues attached to the head still apply -- this is what
             * lets Hold-style heads do useful work (e.g. Length[Hold[a,b,c]]
             * via the down code on Hold). */
            bool hold_all_complete = (attrs & ATTR_HOLDALLCOMPLETE) != 0;

            /* Evaluate arguments unless suppressed by HoldFirst, HoldRest,
             * HoldAll, or HoldAllComplete. */
            Expr** new_args = malloc(sizeof(Expr*) * e->data.function.arg_count);
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                bool hold = hold_all_complete;
                if (i == 0 && (attrs & ATTR_HOLDFIRST)) hold = true;
                if (i > 0 && (attrs & ATTR_HOLDREST)) hold = true;

                if (hold) {
                    /* Check for Evaluate[expr] - overrides HoldFirst/HoldRest/HoldAll
                     * but NOT HoldAllComplete. */
                    Expr* arg = e->data.function.args[i];
                    if (!hold_all_complete &&
                        arg->type == EXPR_FUNCTION &&
                        arg->data.function.head->type == EXPR_SYMBOL &&
                        strcmp(arg->data.function.head->data.symbol, "Evaluate") == 0 &&
                        arg->data.function.arg_count == 1) {
                        new_args[i] = evaluate(arg->data.function.args[0]);
                    } else {
                        new_args[i] = expr_copy(e->data.function.args[i]);
                    }
                } else {
                    new_args[i] = evaluate(e->data.function.args[i]);
                }
            }
            
            Expr* res = expr_new_function(head, new_args, e->data.function.arg_count);
            free(new_args);
            
    /* 2.5 Flatten Sequences - must happen before attributes.
     * Suppressed under HoldAllComplete and inside Set/SetDelayed/Rule/RuleDelayed
     * (whose pattern syntax preserves Sequence as a literal head). */
    if (hold_all_complete) {
        /* HoldAllComplete leaves Sequence intact */
    } else if (head->type == EXPR_SYMBOL &&
        (strcmp(head->data.symbol, "Set") == 0 || strcmp(head->data.symbol, "SetDelayed") == 0 ||
         strcmp(head->data.symbol, "Rule") == 0 || strcmp(head->data.symbol, "RuleDelayed") == 0)) {
        // Do not flatten sequences in assignments or rules
    } else {
        flatten_sequences(res);
    }

            /* 2.6 Strip Unevaluated wrappers in non-held positions.
             * f[Unevaluated[expr]] passes expr (unevaluated) to f, with the
             * wrapper removed. This runs AFTER flatten_sequences so that any
             * Sequence[...] directly inside Unevaluated is preserved (e.g.
             * Length[Unevaluated[Sequence[a,b]]] gives 2 because Sequence is
             * not flattened into Length's argument list).
             * The wrapper is kept in held positions (so Hold[Unevaluated[1+2]]
             * stays as Hold[Unevaluated[1+2]]) and for HoldAllComplete heads
             * (already handled by the early return above). */
            for (size_t i = 0; i < res->data.function.arg_count; i++) {
                bool held = hold_all_complete;
                if (i == 0 && (attrs & ATTR_HOLDFIRST)) held = true;
                if (i > 0 && (attrs & ATTR_HOLDREST)) held = true;
                if (held) continue;

                Expr* arg = res->data.function.args[i];
                if (arg->type == EXPR_FUNCTION &&
                    arg->data.function.head->type == EXPR_SYMBOL &&
                    strcmp(arg->data.function.head->data.symbol, "Unevaluated") == 0 &&
                    arg->data.function.arg_count == 1) {
                    Expr* stripped = expr_copy(arg->data.function.args[0]);
                    expr_free(arg);
                    res->data.function.args[i] = stripped;
                }
            }

            /* 3. Apply structural and semantic attributes.
             * Order follows Withoff §3.1: Flat → Sequence (already done above) →
             * Listable → Orderless. Flat must run before Listable so that lists
             * exposed by flattening get threaded, e.g. Plus[Plus[a,{1,2}],3]
             * → Plus[a,{1,2},3] → {Plus[a,1,3], Plus[a,2,3]}. */

            /* Flat: associative flattening (requires symbolic head, suppressed by HoldAllComplete) */
            if (head->type == EXPR_SYMBOL && (attrs & ATTR_FLAT) && !hold_all_complete) {
                eval_flatten_args(res, head->data.symbol);
            }

            /* Listable: automatic threading */
            if ((attrs & ATTR_LISTABLE) && has_list_arg(res)) {
                Expr* list_res = apply_listable(res);
                if (list_res) {
                    expr_free(res);
                    return list_res;
                }
            }

            if (head->type == EXPR_SYMBOL) {
                const char* head_name = head->data.symbol;

                /* Orderless: commutative sorting */
                if (attrs & ATTR_ORDERLESS) {
                    qsort(res->data.function.args, res->data.function.arg_count, sizeof(Expr*), eval_compare_expr_ptrs);
                }

                /* 4. Apply user-defined DownValues FIRST (Withoff §3.1).
                 * In Mathematica's evaluation pipeline, user-defined
                 * DownValues take precedence over internal "down code"
                 * (built-in implementations). This lets a user override
                 * a built-in for non-Protected symbols, while Protected
                 * symbols (which is most builtin-bearing heads) are
                 * unaffected because apply_assignment refuses to install
                 * DownValues on a Protected target. */
                Expr* down = apply_down_values(res);
                if (down) {
                    expr_free(res);
                    return down;
                }

                /* 5. Call C-level Built-in Functions (internal "down code") */
                SymbolDef* def = symtab_get_def(head_name);
                if (def && def->builtin_func) {
                    Expr* ret = def->builtin_func(res);
                    if (ret) {
                        expr_free(res);
                        return ret;
                    }
                }

                /* 6. Special primitives (Set, SetDelayed) */
                if ((strcmp(head_name, "Set") == 0 || strcmp(head_name, "SetDelayed") == 0) && res->data.function.arg_count == 2) {
                    Expr* lhs = res->data.function.args[0];
                    Expr* rhs = res->data.function.args[1];
                    int is_delayed = (strcmp(head_name, "SetDelayed") == 0);
                    
                    /* For Set and SetDelayed, we evaluate the arguments of the LHS to find the actual target */
                    /* e.g. f[x] = 1 where x=c should define f[c]=1 */
                    /* Patterns must also be evaluated to canonical form to match evaluated inputs. */
                    Expr* target_lhs = lhs;
                    bool free_target = false;
                    if (lhs->type == EXPR_FUNCTION) {
                        /* Only evaluate arguments, not the head, to avoid matching existing rules */
                        Expr** eval_args = malloc(sizeof(Expr*) * lhs->data.function.arg_count);
                        bool is_part = (lhs->data.function.head->type == EXPR_SYMBOL && strcmp(lhs->data.function.head->data.symbol, "Part") == 0);
                        /* List destructuring: {a, b, ...} = {...}. Each element that is
                         * a Symbol is a binding target and must NOT be evaluated (otherwise
                         * prior OwnValues clobber the targets -- e.g. {a,b}={1,2} then
                         * {a,b}={3,4} would try to assign to the values 1,2 instead of a,b).
                         * Non-symbol elements (e.g. a[x] in {a[x], b[y]} = ...) still need
                         * their inner arguments evaluated so the target pattern is correct. */
                        bool is_list = (lhs->data.function.head->type == EXPR_SYMBOL && strcmp(lhs->data.function.head->data.symbol, "List") == 0);

                        uint32_t lhs_attrs = ATTR_NONE;
                        if (lhs->data.function.head->type == EXPR_SYMBOL) {
                            lhs_attrs = get_attributes(lhs->data.function.head->data.symbol);
                        }

                        for (size_t i = 0; i < lhs->data.function.arg_count; i++) {
                            bool hold = false;
                            if ((lhs_attrs & ATTR_HOLDALLCOMPLETE) == ATTR_HOLDALLCOMPLETE) hold = true;
                            else if ((lhs_attrs & ATTR_HOLDALL) == ATTR_HOLDALL) hold = true;
                            else if (i == 0 && (lhs_attrs & ATTR_HOLDFIRST)) hold = true;
                            else if (i > 0 && (lhs_attrs & ATTR_HOLDREST)) hold = true;

                            if (is_part && i == 0) hold = true; // Hold the first argument of Part

                            /* In a List-LHS, hold any element that is itself a symbol or
                             * a nested List (binding targets / nested destructuring). */
                            if (is_list) {
                                Expr* child = lhs->data.function.args[i];
                                if (child->type == EXPR_SYMBOL) {
                                    hold = true;
                                } else if (child->type == EXPR_FUNCTION &&
                                           child->data.function.head->type == EXPR_SYMBOL &&
                                           strcmp(child->data.function.head->data.symbol, "List") == 0) {
                                    hold = true;
                                }
                            }

                            if (hold) {
                                eval_args[i] = expr_copy(lhs->data.function.args[i]);
                            } else {
                                eval_args[i] = evaluate(lhs->data.function.args[i]);
                            }
                        }
                        target_lhs = expr_new_function(expr_copy(lhs->data.function.head), eval_args, lhs->data.function.arg_count);
                        free(eval_args);
                        free_target = true;
                    }

                    if (apply_assignment(target_lhs, rhs, is_delayed)) {
                        Expr* ret = is_delayed ? expr_new_symbol("Null") : evaluate(rhs);
                        if (free_target) expr_free(target_lhs);
                        expr_free(res);
                        return ret;
                    }
                    if (free_target) expr_free(target_lhs);
                }
                
                /* OneIdentity is intentionally NOT rewritten at evaluation
                 * time. In Mathematica it is purely a pattern-matching
                 * attribute: it lets f[x_, y_:def] match a literal `a`.
                 * The 1-arg collapse f[x] -> x is the responsibility of
                 * each head's builtin (Plus, Times, Power, GCD, LCM, And,
                 * Or, Dot all handle the n==1 case explicitly), so a
                 * user-defined OneIdentity head like
                 *   SetAttributes[g, OneIdentity]
                 * leaves g[x] as g[x] rather than rewriting to x.
                 * The pattern-matching half lives in src/match.c (search
                 * for ATTR_ONEIDENTITY). */
            } else if (head->type == EXPR_FUNCTION && head->data.function.head->type == EXPR_SYMBOL &&
                       strcmp(head->data.function.head->data.symbol, "Function") == 0) {

                /* 7. Apply Pure Function */
                Expr* applied = apply_pure_function(head, res->data.function.args, res->data.function.arg_count);
                if (applied) {
                    expr_free(res);
                    return applied;
                }
            } else if (head->type == EXPR_FUNCTION && head->data.function.head->type == EXPR_SYMBOL &&
                       strcmp(head->data.function.head->data.symbol, "Derivative") == 0 &&
                       res->data.function.arg_count == 1 &&
                       res->data.function.args[0]->type == EXPR_FUNCTION &&
                       res->data.function.args[0]->data.function.head->type == EXPR_SYMBOL &&
                       strcmp(res->data.function.args[0]->data.function.head->data.symbol, "Function") == 0) {
                /* 7b. Derivative[n1,...,nm][Function[...]] reduces to a new
                 * Function whose body has been differentiated. Without this
                 * step, f'[x] for a pure-function f would remain stuck as
                 * Derivative[1][Function[...]][x]. */
                Expr* reduced = derivative_of_pure_function(head, res->data.function.args[0]);
                if (reduced) {
                    expr_free(res);
                    return reduced;
                }
            }
            
            return res;
        }
    }
    return expr_copy(e);
}

/*
 * evaluate:
 * The primary evaluator loop. Repeatedly applies evaluate_step until the 
 * expression reaches a fixed point (no further changes) or the iteration 
 * limit is reached.
 */
Expr* evaluate(Expr* e) {
    if (!e) return NULL;

    bool is_top_level = (eval_recursion_depth == 0);
    if (is_top_level) eval_overflow = false;

    /* Guard the C stack: when nested evaluate() calls would exceed the
     * recursion limit, wrap the input in Hold[] so it stops re-entering
     * the evaluator, set the sticky overflow flag so all enclosing
     * fixed-point loops bail, and emit a message exactly once per
     * top-level evaluation. */
    if (eval_recursion_depth >= eval_recursion_limit) {
        if (!eval_overflow) {
            fprintf(stderr,
                    "$RecursionLimit::reclim: Recursion depth of %d exceeded.\n",
                    eval_recursion_limit);
        }
        eval_overflow = true;
        Expr** wrap = malloc(sizeof(Expr*));
        wrap[0] = expr_copy(e);
        Expr* held = expr_new_function(expr_new_symbol("Hold"), wrap, 1);
        free(wrap);
        return held;
    }

    eval_recursion_depth++;

    Expr* current = expr_copy(e);
    Expr* next = NULL;
    int iterations = 0;

    while (iterations < MAX_ITERATIONS) {
        next = evaluate_step(current);

        /* Fixed point check: if step produced an identical expression, we are done */
        if (expr_eq(current, next)) {
            expr_free(next);
            eval_recursion_depth--;
            return current;
        }

        /* Prepare for the next iteration */
        expr_free(current);
        current = next;
        iterations++;

        /* If a deeper call hit the recursion limit, the rewrites
         * above are no longer making real progress -- bail out so the
         * unwind doesn't burn $IterationLimit at every level. */
        if (eval_overflow) break;
    }

    if (iterations >= MAX_ITERATIONS) {
        fprintf(stderr, "$IterationLimit exceeded\n");
    }
    eval_recursion_depth--;
    return current;
}
