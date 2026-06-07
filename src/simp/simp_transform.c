/*
 * simp_transform.c -- user-supplied TransformationFunctions for Simplify.
 *
 * Mathematica's Simplify and FullSimplify take a TransformationFunctions
 * option giving the list of functions to apply to try to transform parts of
 * an expression. The default (Automatic) uses the built-in transformation
 * collection; an explicit list of functions {f1, f2, ...} uses only those;
 * {Automatic, f1, ...} uses the built-ins together with the user functions.
 *
 * builtin_simplify (simp_builtins.c) parses the option, decides whether the
 * built-in pipeline runs, and -- when user functions are present -- hands the
 * best-so-far expression to simp_apply_transformations below.
 *
 * The driver is a bounded global-best fixed point: each round it applies every
 * user function to the whole expression and to every subexpression, keeps the
 * candidate of strictly lowest complexity (scored with the same
 * score_with_func used everywhere else in the search), and repeats until no
 * function lowers the score or a budget is exhausted. Scoring the *whole*
 * expression after each subexpression rewrite -- rather than the rewritten part
 * in isolation -- matches Mathematica's "simplest overall form wins" semantics
 * and lets a local rewrite that enables a global collapse be accepted.
 */

#include "simp.h"
#include "simp_internal.h"
#include "eval.h"
#include "expr.h"

#include <stdlib.h>

/* Budgets. A user can supply arbitrary functions, so we bound the work:
 *   - MAX_ROUNDS   fixed-point iterations (each round restarts on improvement);
 *   - NODE_CAP     subexpressions scanned per round (DFS pre-order, capped);
 *   - EVAL_BUDGET  total function applications across the whole call.
 * The caps are generous for the small expressions Simplify handles while
 * guaranteeing termination on pathological transformation functions. */
#define TF_MAX_ROUNDS  32
#define TF_NODE_CAP    400
#define TF_EVAL_BUDGET 20000

/* Number of nodes in DFS pre-order (root + all arguments, recursively). The
 * head of a function is not treated as a separate transformable node -- it is
 * almost always a bare symbol and applying f to it is never useful. */
static size_t tf_node_count(const Expr* e) {
    if (!e) return 0;
    size_t n = 1;
    if (e->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            n += tf_node_count(e->data.function.args[i]);
    }
    return n;
}

/* Borrowed pointer to the DFS pre-order node at index `target` (0 = root).
 * `*idx` must start at 0. Returns NULL if target is out of range. */
static const Expr* tf_node_at(const Expr* root, size_t target, size_t* idx) {
    if (*idx == target) return root;
    (*idx)++;
    if (root->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < root->data.function.arg_count; i++) {
            const Expr* r = tf_node_at(root->data.function.args[i], target, idx);
            if (r) return r;
        }
    }
    return NULL;
}

/* Deep copy of `root` with the DFS pre-order node at index `target` replaced
 * by a copy of `repl`. `*idx` must start at 0. The pre-order numbering matches
 * tf_node_at for every node up to and including `target`, so the replacement
 * lands on exactly the node tf_node_at would return. */
static Expr* tf_copy_replace_at(const Expr* root, size_t target,
                                const Expr* repl, size_t* idx) {
    if (*idx == target) {
        (*idx)++;
        return expr_copy((Expr*)repl);
    }
    (*idx)++;
    if (root->type == EXPR_FUNCTION) {
        size_t n = root->data.function.arg_count;
        Expr** args = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
        for (size_t i = 0; i < n; i++)
            args[i] = tf_copy_replace_at(root->data.function.args[i],
                                         target, repl, idx);
        Expr* out = expr_new_function(expr_copy((Expr*)root->data.function.head),
                                      args, n);
        free(args);
        return out;
    }
    return expr_copy((Expr*)root);
}

/* Apply the user-supplied transformation functions to `expr` and return the
 * lowest-complexity result found. `funcs` are borrowed (not freed). The
 * returned expression is freshly allocated and owned by the caller; it is
 * never NULL (worst case it is a copy of the input). */
Expr* simp_apply_transformations(const Expr* expr, Expr* const* funcs,
                                 size_t nfuncs, const Expr* complexity_func) {
    Expr* best = expr_copy((Expr*)expr);
    if (nfuncs == 0) return best;

    size_t best_score = score_with_func(best, complexity_func);
    size_t evals = 0;

    for (int round = 0; round < TF_MAX_ROUNDS; round++) {
        bool improved = false;
        size_t n = tf_node_count(best);
        if (n > TF_NODE_CAP) n = TF_NODE_CAP;

        for (size_t node = 0; node < n && !improved; node++) {
            size_t walk = 0;
            const Expr* sub = tf_node_at(best, node, &walk);
            if (!sub) break;

            for (size_t fi = 0; fi < nfuncs; fi++) {
                if (evals >= TF_EVAL_BUDGET) { round = TF_MAX_ROUNDS; break; }
                evals++;

                /* Build and evaluate funcs[fi][sub]. */
                Expr** call_args = (Expr**)calloc(1, sizeof(Expr*));
                call_args[0] = expr_copy((Expr*)sub);
                Expr* call = expr_new_function(expr_copy(funcs[fi]),
                                               call_args, 1);
                free(call_args);
                Expr* val = eval_and_free(call);
                if (!val) continue;
                if (expr_eq(val, sub)) { expr_free(val); continue; }

                /* Splice the transformed part back into the whole, then
                 * re-evaluate so the surrounding context can fold/canonicalise
                 * around it before we score the global form. */
                size_t place = 0;
                Expr* cand = tf_copy_replace_at(best, node, val, &place);
                expr_free(val);
                Expr* cand_eval = eval_and_free(cand);

                size_t s = score_with_func(cand_eval, complexity_func);
                if (s < best_score) {
                    expr_free(best);
                    best = cand_eval;
                    best_score = s;
                    improved = true;
                    break; /* restart the scan from the new best */
                }
                expr_free(cand_eval);
            }
        }

        if (!improved) break;
    }

    return best;
}
