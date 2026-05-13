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
#include "sym_names.h"
#include "sym_intern.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

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
 * The default of 1024 matches modern Mathematica and leaves comfortable
 * headroom under the typical 8 MB thread stack while still catching
 * pathological recursion (e.g. a self-referential rule like
 * f[x_] := f[x] + 1).
 *
 * Minimum enforced at 20 (Mathematica-compatible) so users cannot brick
 * the evaluator by setting a value below the depth its own bookkeeping
 * needs.
 */
#define DEFAULT_RECURSION_LIMIT 1024
#define MIN_RECURSION_LIMIT     20
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
    eval_recursion_limit = (n >= MIN_RECURSION_LIMIT) ? n : DEFAULT_RECURSION_LIMIT;
}

/* See eval.h.  Used after a siglongjmp out of evaluate(): the matching
 * decrements never ran, so we restore the depth counter to the value it
 * had before the timed call entered the evaluator.  We also clear
 * eval_overflow so a future evaluate() is not falsely poisoned by the
 * aborted call. */
void eval_reset_recursion_depth(int n) {
    if (n < 0) n = 0;
    eval_recursion_depth = n;
    eval_overflow = false;
}

/* M3 phase-3 evaluation clock. Starts at 1 so a freshly-allocated Expr
 * (last_evaluated_at == 0) is never mistaken for "already evaluated".
 * Bumped by symtab.c (add_rule, symtab_clear_symbol) and attr.c
 * (set_attributes, add/remove_single_attribute). 64 bits is enough to
 * absorb roughly 2^64 mutations, far beyond any practical session. */
static uint64_t g_eval_clock = 1;
uint64_t eval_clock_get(void) { return g_eval_clock; }
void     eval_clock_bump(void) { g_eval_clock++; }

/*
 * eval_classify_return:
 * See the contract in eval.h. Pointer-equality on interned symbols is
 * the dispatch primitive: every Expr_Symbol's `data.symbol` field is the
 * canonical interned pointer (sym_intern.c), so checks like
 * `head->data.symbol == SYM_Return` and `target->data.symbol ==
 * boundary_head` are O(1) and never strcmp.
 *
 * Care is taken to keep this side-effect free: no eval_clock_bump,
 * no expr_free, no allocation when the answer is NONE/PROPAGATE. The
 * single allocation on CONSUME is either a fresh `Null` symbol (for the
 * 0-arg form) or an expr_copy of args[0] (for 1-arg / 2-arg forms).
 * That copy is necessary because the caller will free `e` after
 * yielding the value.
 */
EvalReturnAction eval_classify_return(Expr* e,
                                      const char* boundary_head,
                                      Expr** out_value) {
    if (out_value) *out_value = NULL;
    if (!e) return EVAL_RETURN_NONE;
    if (e->type != EXPR_FUNCTION) return EVAL_RETURN_NONE;
    if (e->data.function.head->type != EXPR_SYMBOL) return EVAL_RETURN_NONE;
    if (e->data.function.head->data.symbol != SYM_Return) return EVAL_RETURN_NONE;

    size_t argc = e->data.function.arg_count;

    /* Return[]: yield Null at the nearest boundary. */
    if (argc == 0) {
        if (out_value) *out_value = expr_new_symbol("Null");
        return EVAL_RETURN_CONSUME;
    }

    /* Return[expr]: yield expr at the nearest boundary, regardless of
     * boundary_head. */
    if (argc == 1) {
        if (out_value) *out_value = expr_copy(e->data.function.args[0]);
        return EVAL_RETURN_CONSUME;
    }

    /* Return[expr, h, ...]: target the nearest boundary whose head is h.
     * Extra arguments are ignored (Mathematica accepts up to 2 args; we
     * accept the same and leave further args to a caller-side message
     * if the user supplies them).
     *
     * The target must be an EXPR_SYMBOL for the comparison to be
     * meaningful. If it isn't, fall back to PROPAGATE so that the
     * marker isn't accidentally consumed by an arbitrary boundary. */
    Expr* target = e->data.function.args[1];
    if (boundary_head &&
        target->type == EXPR_SYMBOL &&
        target->data.symbol == boundary_head) {
        if (out_value) *out_value = expr_copy(e->data.function.args[0]);
        return EVAL_RETURN_CONSUME;
    }
    return EVAL_RETURN_PROPAGATE;
}

/*
 * eval_init:
 * Registers the user-visible $RecursionLimit symbol with its default value
 * as an OwnValue so that the user can read or assign to it from the REPL.
 * The C-side state is kept in sync via the hook in apply_assignment.
 *
 * Must be called after symtab_init().
 */
void eval_init(void) {
    Expr* sym = expr_new_symbol("$RecursionLimit");
    Expr* val = expr_new_integer(eval_recursion_limit);
    symtab_add_own_value("$RecursionLimit", sym, val);
    expr_free(sym);
    expr_free(val);

    symtab_set_docstring("$RecursionLimit",
        "$RecursionLimit\n"
        "\tgives the maximum length of the evaluation stack -- the maximum\n"
        "\tnumber of nested invocations of the evaluator that can occur.\n"
        "\n"
        "Assigning a positive integer N (>= 20) updates the limit; smaller\n"
        "values are rejected with a $RecursionLimit::limset message.");
}

/*
 * sync_recursion_limit_from_value:
 * Inspect a candidate value (typically the RHS of $RecursionLimit = ...)
 * and, if it is a positive integer >= MIN_RECURSION_LIMIT, push it into
 * the C-level limit. Otherwise emit a $RecursionLimit::limset message and
 * leave the C state untouched. Bigints are clamped to INT_MAX.
 */
static void sync_recursion_limit_from_value(Expr* value) {
    long n = -1;
    if (value->type == EXPR_INTEGER) {
        n = (long)value->data.integer;
    } else if (value->type == EXPR_BIGINT) {
        /* Anything large enough not to fit in a long is far beyond any
         * useful recursion limit; treat it as "huge and acceptable". */
        if (mpz_fits_slong_p(value->data.bigint)) {
            n = mpz_get_si(value->data.bigint);
        } else if (mpz_sgn(value->data.bigint) > 0) {
            n = (long)INT_MAX;
        }
    }

    if (n < MIN_RECURSION_LIMIT) {
        fprintf(stderr,
                "$RecursionLimit::limset: Cannot set $RecursionLimit to a value below %d.\n",
                MIN_RECURSION_LIMIT);
        /* Restore the OwnValue to the current C-side limit so the symbol
         * does not lie about the active value. */
        Expr* sym  = expr_new_symbol("$RecursionLimit");
        Expr* curr = expr_new_integer(eval_recursion_limit);
        symtab_add_own_value("$RecursionLimit", sym, curr);
        expr_free(sym);
        expr_free(curr);
        return;
    }
    if (n > INT_MAX) n = INT_MAX;
    eval_set_recursion_limit((int)n);
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
/* Returns true iff the call actually flattened nested same-head children
 * (i.e. produced a structurally different argument list). When false,
 * `e` is byte-for-byte unchanged and the §3.4 fixed-point detector can
 * count this step as a no-op. */
bool eval_flatten_args(Expr* e, const char* head_name) {
    /* Callers may hand us a C-string literal (e.g. internal_call_impl
     * passes "Plus") rather than the interned canonical pointer. Funnel
     * through the interner so the per-arg head check below is a pointer
     * compare against the same pointer that lives on every interned
     * EXPR_SYMBOL. */
    head_name = intern_symbol(head_name);

    size_t new_count = 0;
    bool needs_flattening = false;

    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol == head_name) {
            new_count += arg->data.function.arg_count;
            needs_flattening = true;
        } else {
            new_count++;
        }
    }

    /* If no nested occurrences of the head were found, we are done */
    if (!needs_flattening) return false;

    /* Second pass: allocate new argument array and copy elements */
    Expr** new_args = malloc(sizeof(Expr*) * new_count);
    size_t idx = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol == head_name) {
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
    return true;
}

/*
 * has_list_arg:
 * Helper for the Listable attribute.
 * Checks if any argument of the function is an explicit List[...].
 */
static bool has_list_arg(Expr* e) {
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_FUNCTION &&
            arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol == SYM_List) {
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
        if (arg->type == EXPR_FUNCTION &&
            arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol == SYM_List) {
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
            if (arg->type == EXPR_FUNCTION &&
            arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol == SYM_List) {
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
        if (h == SYM_Condition || h == SYM_HoldPattern || h == SYM_Part) {
            return assignment_target_symbol(lhs->data.function.args[0]);
        }
        return h;
    }
    return NULL;
}

/*
 * is_assignable_lhs:
 * Validate that `lhs` shaped against `rhs` is a structurally legal target
 * for a Set/SetDelayed, including all destructured sub-elements. Used as a
 * pre-flight check on List destructuring so a malformed element (e.g. a
 * literal integer on the LHS) cannot land partial assignments before being
 * detected by the in-loop failure path.
 */
static bool is_assignable_lhs(Expr* lhs, Expr* rhs) {
    if (!lhs) return false;
    if (lhs->type == EXPR_SYMBOL) return true;
    if (lhs->type != EXPR_FUNCTION) return false;
    if (lhs->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = lhs->data.function.head->data.symbol;
    if (h != SYM_List) return true; /* downvalue / part / etc. handled in apply_assignment */

    if (rhs->type != EXPR_FUNCTION ||
        rhs->data.function.head->type != EXPR_SYMBOL ||
        rhs->data.function.head->data.symbol != SYM_List) {
        return false;
    }
    if (lhs->data.function.arg_count != rhs->data.function.arg_count) return false;
    for (size_t i = 0; i < lhs->data.function.arg_count; i++) {
        if (!is_assignable_lhs(lhs->data.function.args[i], rhs->data.function.args[i])) {
            return false;
        }
    }
    return true;
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
                        lhs->data.function.head->data.symbol == SYM_List);
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

        /* Special system variables: keep their C-side mirror state in sync.
         * Set has HoldFirst, so for `$RecursionLimit = expr` the rhs is
         * already evaluated; for SetDelayed, we evaluate a copy here so the
         * limit reflects the value the user expects to see when they read
         * the symbol back. If validation fails, the OwnValue is rolled back
         * to the current C-side limit. */
        if (strcmp(lhs->data.symbol, "$RecursionLimit") == 0) {
            Expr* probe = is_delayed ? evaluate(expr_copy(rhs)) : expr_copy(rhs);
            sync_recursion_limit_from_value(probe);
            expr_free(probe);
        }
        return true;
    } else if (lhs->type == EXPR_FUNCTION) {
        if (lhs->data.function.head->type == EXPR_SYMBOL &&
            lhs->data.function.head->data.symbol == SYM_List &&
            rhs->type == EXPR_FUNCTION &&
            rhs->data.function.head->type == EXPR_SYMBOL &&
            rhs->data.function.head->data.symbol == SYM_List) {
            
            /* List destructuring: match lengths and recurse. Pre-flight every
             * element so a malformed child (e.g. a literal integer on the LHS)
             * fails the whole destructuring before any sibling is assigned --
             * partial assignments would otherwise leak past the failure. */
            if (lhs->data.function.arg_count != rhs->data.function.arg_count) {
                return false;
            }
            for (size_t i = 0; i < lhs->data.function.arg_count; i++) {
                if (!is_assignable_lhs(lhs->data.function.args[i], rhs->data.function.args[i])) {
                    return false;
                }
            }

            bool all_ok = true;
            for (size_t i = 0; i < lhs->data.function.arg_count; i++) {
                if (!apply_assignment(lhs->data.function.args[i], rhs->data.function.args[i], is_delayed)) {
                    all_ok = false;
                }
            }
            return all_ok;
        } else if (lhs->data.function.head->type == EXPR_SYMBOL && lhs->data.function.head->data.symbol == SYM_Part) {
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
                rhs->data.function.head->data.symbol == SYM_Condition &&
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

            if (symbol_name == SYM_Condition && actual_pattern->data.function.arg_count == 2) {
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
 * Returns true iff the args list was actually rewritten.
 */
static bool flatten_sequences(Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;

    size_t new_count = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol == SYM_Sequence) {
            new_count += arg->data.function.arg_count;
        } else {
            new_count++;
        }
    }

    if (new_count == e->data.function.arg_count) return false;
    
    Expr** new_args = malloc(sizeof(Expr*) * new_count);
    size_t k = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* arg = e->data.function.args[i];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol == SYM_Sequence) {
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
    return true;
}

/*
 * evaluate_step:
 * Performs exactly one level of evaluation transformation.
 *
 * `changed` is an out-parameter set to true iff a real rewrite fired
 * during this step (M3 §3.4 — eager early-exit fixed-point loop). The
 * outer loop in evaluate() uses it to skip the O(tree) expr_eq compare.
 * See the contract in eval.h. NULL is permitted for callers that do not
 * care about the signal.
 */
Expr* evaluate_step(Expr* e, bool* changed) {
    /* Local sink so we can write through `*changed` unconditionally. */
    bool sink = false;
    if (!changed) changed = &sink;
    *changed = false;

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
            if (own) { *changed = true; return own; }
            return expr_copy(e);
        }
            
        case EXPR_FUNCTION: {
            /* 1. Evaluate the head recursively (e.g. f[x][y]).
             * Refcount sharing means evaluate() returns the same pointer
             * when nothing rewrote — pointer-inequality is a sound
             * "head changed" signal, and we only fall through to the
             * structural compare if needed. */
            Expr* orig_head = e->data.function.head;
            Expr* head = evaluate(orig_head);
            if (head != orig_head) *changed = true;

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
             * HoldAll, or HoldAllComplete.
             *
             * Pointer-identity check against the original arg signals
             * "did sub-evaluation change anything"; with refcount sharing
             * (M3 phase-2) evaluate() returns the same pointer for
             * already-stable inputs. Stripping an Evaluate[] wrapper is
             * itself a rewrite even if the wrapped expression evaluates
             * to itself — flag it explicitly. */
            Expr** new_args = malloc(sizeof(Expr*) * e->data.function.arg_count);
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                bool hold = hold_all_complete;
                if (i == 0 && (attrs & ATTR_HOLDFIRST)) hold = true;
                if (i > 0 && (attrs & ATTR_HOLDREST)) hold = true;

                Expr* orig_arg = e->data.function.args[i];
                if (hold) {
                    /* Check for Evaluate[expr] - overrides HoldFirst/HoldRest/HoldAll
                     * but NOT HoldAllComplete. */
                    if (!hold_all_complete &&
                        orig_arg->type == EXPR_FUNCTION &&
                        orig_arg->data.function.head->type == EXPR_SYMBOL &&
                        orig_arg->data.function.head->data.symbol == SYM_Evaluate &&
                        orig_arg->data.function.arg_count == 1) {
                        new_args[i] = evaluate(orig_arg->data.function.args[0]);
                        *changed = true; /* Evaluate[] wrapper stripped */
                    } else {
                        new_args[i] = expr_copy(orig_arg);
                    }
                } else {
                    new_args[i] = evaluate(orig_arg);
                    if (new_args[i] != orig_arg) *changed = true;
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
        (head->data.symbol == SYM_Set || head->data.symbol == SYM_SetDelayed ||
         head->data.symbol == SYM_Rule || head->data.symbol == SYM_RuleDelayed)) {
        // Do not flatten sequences in assignments or rules
    } else {
        if (flatten_sequences(res)) *changed = true;
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
                    arg->data.function.head->data.symbol == SYM_Unevaluated &&
                    arg->data.function.arg_count == 1) {
                    Expr* stripped = expr_copy(arg->data.function.args[0]);
                    expr_free(arg);
                    res->data.function.args[i] = stripped;
                    *changed = true; /* Unevaluated wrapper removed */
                }
            }

            /* 3. Apply structural and semantic attributes.
             * Order follows Withoff §3.1: Flat → Sequence (already done above) →
             * Listable → Orderless. Flat must run before Listable so that lists
             * exposed by flattening get threaded, e.g. Plus[Plus[a,{1,2}],3]
             * → Plus[a,{1,2},3] → {Plus[a,1,3], Plus[a,2,3]}. */

            /* Flat: associative flattening (requires symbolic head, suppressed by HoldAllComplete) */
            if (head->type == EXPR_SYMBOL && (attrs & ATTR_FLAT) && !hold_all_complete) {
                if (eval_flatten_args(res, head->data.symbol)) *changed = true;
            }

            /* Listable: automatic threading */
            if ((attrs & ATTR_LISTABLE) && has_list_arg(res)) {
                Expr* list_res = apply_listable(res);
                if (list_res) {
                    expr_free(res);
                    *changed = true; /* List threading reshaped the call */
                    return list_res;
                }
            }

            if (head->type == EXPR_SYMBOL) {
                const char* head_name = head->data.symbol;

                /* Orderless: commutative sorting. Pre-check whether the
                 * args are already in canonical order so the §3.4 detector
                 * can skip a no-op qsort on stable expressions
                 * (Plus[a,b,c] re-evaluating, etc.). */
                if (attrs & ATTR_ORDERLESS) {
                    bool already_sorted = true;
                    for (size_t i = 1; i < res->data.function.arg_count; i++) {
                        if (eval_compare_expr_ptrs(&res->data.function.args[i - 1],
                                                   &res->data.function.args[i]) > 0) {
                            already_sorted = false;
                            break;
                        }
                    }
                    if (!already_sorted) {
                        qsort(res->data.function.args, res->data.function.arg_count, sizeof(Expr*), eval_compare_expr_ptrs);
                        *changed = true;
                    }
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
                    *changed = true; /* DownValue rule fired */
                    return down;
                }

                /* 5. Call C-level Built-in Functions (internal "down code") */
                SymbolDef* def = symtab_get_def(head_name);
                if (def && def->builtin_func) {
                    Expr* ret = def->builtin_func(res);
                    if (ret) {
                        expr_free(res);
                        *changed = true; /* Built-in produced a rewrite */
                        return ret;
                    }
                }

                /* 6. Special primitives (Set, SetDelayed) */
                if ((head_name == SYM_Set || head_name == SYM_SetDelayed) && res->data.function.arg_count == 2) {
                    Expr* lhs = res->data.function.args[0];
                    Expr* rhs = res->data.function.args[1];
                    int is_delayed = (head_name == SYM_SetDelayed);
                    
                    /* For Set and SetDelayed, we evaluate the arguments of the LHS to find the actual target */
                    /* e.g. f[x] = 1 where x=c should define f[c]=1 */
                    /* Patterns must also be evaluated to canonical form to match evaluated inputs. */
                    Expr* target_lhs = lhs;
                    bool free_target = false;
                    if (lhs->type == EXPR_FUNCTION) {
                        /* Only evaluate arguments, not the head, to avoid matching existing rules */
                        Expr** eval_args = malloc(sizeof(Expr*) * lhs->data.function.arg_count);
                        bool is_part = (lhs->data.function.head->type == EXPR_SYMBOL && lhs->data.function.head->data.symbol == SYM_Part);
                        /* List destructuring: {a, b, ...} = {...}. Each element that is
                         * a Symbol is a binding target and must NOT be evaluated (otherwise
                         * prior OwnValues clobber the targets -- e.g. {a,b}={1,2} then
                         * {a,b}={3,4} would try to assign to the values 1,2 instead of a,b).
                         * Non-symbol elements (e.g. a[x] in {a[x], b[y]} = ...) still need
                         * their inner arguments evaluated so the target pattern is correct. */
                        bool is_list = (lhs->data.function.head->type == EXPR_SYMBOL && lhs->data.function.head->data.symbol == SYM_List);

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
                                           child->data.function.head->data.symbol == SYM_List) {
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
                        *changed = true; /* Set/SetDelayed installed a rule */
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
                       head->data.function.head->data.symbol == SYM_Function) {

                /* 7. Apply Pure Function */
                Expr* applied = apply_pure_function(head, res->data.function.args, res->data.function.arg_count);
                if (applied) {
                    expr_free(res);
                    *changed = true; /* Pure Function applied */
                    return applied;
                }
            } else if (head->type == EXPR_FUNCTION && head->data.function.head->type == EXPR_SYMBOL &&
                       head->data.function.head->data.symbol == SYM_Derivative &&
                       res->data.function.arg_count == 1 &&
                       res->data.function.args[0]->type == EXPR_FUNCTION &&
                       res->data.function.args[0]->data.function.head->type == EXPR_SYMBOL &&
                       res->data.function.args[0]->data.function.head->data.symbol == SYM_Function) {
                /* 7b. Derivative[n1,...,nm][Function[...]] reduces to a new
                 * Function whose body has been differentiated. Without this
                 * step, f'[x] for a pure-function f would remain stuck as
                 * Derivative[1][Function[...]][x]. */
                Expr* reduced = derivative_of_pure_function(head, res->data.function.args[0]);
                if (reduced) {
                    expr_free(res);
                    *changed = true; /* Derivative-of-Function reduced */
                    return reduced;
                }
            } else if (head->type == EXPR_FUNCTION && head->data.function.head->type == EXPR_SYMBOL &&
                       head->data.function.head->data.symbol == SYM_Composition &&
                       head->data.function.arg_count >= 1) {
                /* 7c. Composition[f1, ..., fn][args...] -> f1[f2[...[fn[args...]]]].
                 * The innermost call carries all the user-supplied arguments;
                 * each outer fk wraps the previous result as a single argument. */
                size_t nf = head->data.function.arg_count;
                size_t na = res->data.function.arg_count;
                Expr** call_args = malloc(sizeof(Expr*) * na);
                for (size_t i = 0; i < na; i++) {
                    call_args[i] = expr_copy(res->data.function.args[i]);
                }
                Expr* inner = expr_new_function(
                    expr_copy(head->data.function.args[nf - 1]),
                    call_args, na);
                free(call_args);
                for (size_t k = nf - 1; k > 0; k--) {
                    Expr* one[1] = { inner };
                    inner = expr_new_function(
                        expr_copy(head->data.function.args[k - 1]),
                        one, 1);
                }
                expr_free(res);
                *changed = true; /* Composition unrolled */
                return inner;
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

    /* M3 phase-3 timestamp early-exit. If this expression has been fully
     * evaluated under the current symbol-table state (clock unchanged
     * since), return an inc-ref'd view immediately and skip the entire
     * fixed-point loop and all of evaluate_step. This lifts the cost of
     * a re-evaluation from O(tree size) to O(1). Atoms and bare symbols
     * already short-circuit cheaply inside evaluate_step (atom returns
     * expr_copy, symbol-no-OwnValue returns expr_copy), so we limit the
     * pre-check to FUNCTION nodes -- both to avoid an extra branch on
     * the common atom path and because atoms are never expensive to
     * "re-evaluate" anyway. */
    if (e->type == EXPR_FUNCTION && e->last_evaluated_at == g_eval_clock) {
        return expr_copy(e);
    }

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
        /* TimeConstrained's cooperative wall-clock backstop.  No-op
         * unless we're inside an active TimeConstrained call; on
         * hosts where ITIMER_PROF/SIGPROF are reliable the signal
         * normally fires first and this check stays a cheap inactive
         * read.  On hosts where they aren't (WSL 1), this is what
         * actually enforces the deadline.  When tripped, the call
         * siglongjmp's out of this loop straight to TimeConstrained's
         * sigsetjmp; no further cleanup runs here, exactly matching
         * the signal-handler path. */
        tc_check_deadline();

        bool step_changed = false;
        next = evaluate_step(current, &step_changed);

        /* M3 phase-4 (§3.4): eager early exit. evaluate_step signals via
         * the `step_changed` out-parameter whether any rewrite fired
         * during the step (head re-evaluation, arg evaluation, Sequence
         * flatten, Unevaluated strip, Flat, Listable, Orderless,
         * DownValue, built-in, special primitive, pure Function,
         * Derivative-of-Function, Composition-unfold). When nothing
         * fired, the result is structurally identical to `current` and
         * we are at a fixed point — skip the O(tree) expr_eq compare.
         *
         * Some built-ins (Plus, Times, ...) unconditionally rebuild
         * their output even when no terms combined; those trip the
         * change flag without producing a structural difference. Use
         * expr_eq as a fallback in the changed-true branch so those
         * "false positives" still converge in one iteration, matching
         * the old semantics. Cost: identical to the pre-§3.4 path on
         * the slow case; the win is the cheap boolean fast-path on the
         * common case where nothing fires (atoms, bare symbols, fully
         * reduced functions). */
        bool is_fixed_point = !step_changed || expr_eq(current, next);
        if (is_fixed_point) {
            expr_free(next);
            /* M3 phase-3: stamp the fully-evaluated result with the
             * current clock so a subsequent evaluate(current) hits the
             * early-exit above. We deliberately stamp ONLY on a clean
             * fixed-point exit; the iteration-cap and recursion-overflow
             * paths below leave the timestamp untouched, so a later
             * evaluator gets a fresh chance to make progress. The
             * write is benign metadata, so it is safe even when
             * `current` is shared (refcount > 1). */
            if (!eval_overflow) {
                current->last_evaluated_at = g_eval_clock;
            }
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
