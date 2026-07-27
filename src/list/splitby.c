/* SplitBy[list, f] — split a list into runs of *consecutive* elements that
 * share the same value of f[element].
 *
 * This is the key-function counterpart of Split (src/list/split.c): Split
 * compares adjacent elements directly (or via a two-argument test), whereas
 * SplitBy compares the evaluated keys f[e]. Only adjacent elements are ever
 * grouped, which is what distinguishes SplitBy from GatherBy (src/assoc.c) —
 * GatherBy collects *all* elements sharing a key, no matter where they sit.
 *
 *   SplitBy[{1, 3, 2, 4, 5}, EvenQ]     -> {{1, 3}, {2, 4}, {5}}
 *   SplitBy[{1, 2, 3, 4, 5, 6}, EvenQ]  -> {{1}, {2}, {3}, {4}, {5}, {6}}
 *   SplitBy[{1, 1, 2, 2, 3}, Identity]  -> {{1, 1}, {2, 2}, {3}}
 *
 * The list form SplitBy[list, {f1, f2, ...}] splits by f1, then splits each
 * resulting run by f2, and so on, nesting one level deeper per function:
 *
 *   SplitBy[{1, 3, 2, 4}, {EvenQ}]      -> {{1, 3}, {2, 4}}
 *
 * Cost: f is evaluated exactly once per element per level, i.e. O(n) calls per
 * function in the key spec, plus O(n) structural copying. Keys are compared
 * with expr_eq — the same structural equality the rest of the kernel uses — so
 * two adjacent elements whose keys stay unevaluated but identical still group
 * together. Only one key is held live at a time (the previous element's), so
 * peak overhead beyond the result itself is a single key plus the per-level
 * run vector. */

#include "list_common.h"
#include "splitby.h"

/* Evaluate f[e] and return the resulting key. Caller owns the result. */
static Expr* splitby_key(Expr* f, Expr* e) {
    Expr* args[1] = { expr_copy(e) };
    Expr* call = expr_new_function(expr_copy(f), args, 1);
    Expr* key = evaluate(call);
    expr_free(call);
    return key;
}

/* Build one run: elems[0..n) copied under a fresh `head`. */
static Expr* splitby_make_run(Expr** elems, size_t n, Expr* head) {
    Expr** args = n ? malloc(sizeof(Expr*) * n) : NULL;
    for (size_t i = 0; i < n; i++) args[i] = expr_copy(elems[i]);
    Expr* run = expr_new_function(expr_copy(head), args, n);
    free(args);
    return run;
}

/* Split elems[0..n) into runs by fs[0], recursing into fs[1..nfs) so each run
 * is itself split one level deeper. Returns the runs collected under `head`;
 * n == 0 yields an empty expression, so SplitBy[{}, f] is {}. */
static Expr* splitby_level(Expr** elems, size_t n, Expr* head,
                           Expr* const* fs, size_t nfs) {
    Expr** runs = malloc(sizeof(Expr*) * (n ? n : 1));
    size_t num_runs = 0;
    size_t run_start = 0;

    /* Key of elems[i-1], carried across iterations so f runs once per element. */
    Expr* prev_key = n ? splitby_key(fs[0], elems[0]) : NULL;

    for (size_t i = 1; i <= n; i++) {
        bool cut = (i == n);          /* the final run always closes at the end */
        Expr* key = NULL;
        if (!cut) {
            key = splitby_key(fs[0], elems[i]);
            cut = !expr_eq(key, prev_key);
        }
        if (key) { expr_free(prev_key); prev_key = key; }

        if (cut) {
            size_t len = i - run_start;
            runs[num_runs++] = (nfs > 1)
                ? splitby_level(elems + run_start, len, head, fs + 1, nfs - 1)
                : splitby_make_run(elems + run_start, len, head);
            run_start = i;
        }
    }
    expr_free(prev_key);

    Expr* result = expr_new_function(expr_copy(head), runs, num_runs);
    free(runs);
    return result;
}

Expr* builtin_splitby(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;

    Expr* list  = res->data.function.args[0];
    Expr* fspec = res->data.function.args[1];

    /* Atoms have no elements to split; leave the call unevaluated. */
    if (list->type != EXPR_FUNCTION) return NULL;

    /* fspec is either a single key function or a List of them. An empty list
     * specifies no split at all, which we leave unevaluated rather than
     * guessing at an interpretation. */
    Expr* single[1];
    Expr* const* fs;
    size_t nfs;
    if (head_is(fspec, SYM_List)) {
        nfs = fspec->data.function.arg_count;
        if (nfs == 0) return NULL;
        fs = fspec->data.function.args;
    } else {
        single[0] = fspec;
        fs = single;
        nfs = 1;
    }

    return splitby_level(list->data.function.args, list->data.function.arg_count,
                         list->data.function.head, fs, nfs);
}
