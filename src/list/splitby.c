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
 * with SameQ (expr_sameq: structural identity, machine reals within
 * Internal`$SameQTolerance), so 0.1 + 0.2 and 0.3 are one key, and two
 * adjacent elements whose keys stay unevaluated but identical still group
 * together. Only one key is held live at a time (the previous element's), so
 * peak overhead beyond the result itself is a single key plus the per-level
 * run vector. */

#include "eval.h"
#include "list_common.h"
#include "splitby.h"
#include "linkedlist_seq.h"
#include "comparisons.h"   /* expr_sameq */

/* Evaluate f[e] and return the resulting key. Caller owns the result. */
static Expr* splitby_key(Expr* f, Expr* e) {
    Expr* args[1] = { expr_copy(e) };
    Expr* call = expr_new_function(expr_copy(f), args, 1);
    Expr* key = eval_callback(call);
    expr_free(call);
    return key;
}

/* Split elements [lo, hi) of q into runs by f_k (element k of fq, a List or
 * LinkedList of functions -- or the single function `f1` when fq is NULL),
 * recursing with k + 1 so each run is split one level deeper. Elements are
 * read on demand; a run is copied span-wise (llq_range). An empty range
 * yields an empty expression, so SplitBy[{}, f] is {}. */
static Expr* splitby_level(LLSeq* q, size_t lo, size_t hi, LLSeq* fq, Expr* f1, size_t k, size_t nfs) {
    size_t n = hi - lo;
    Expr** runs = malloc(sizeof(Expr*) * (n ? n : 1));
    size_t num_runs = 0;
    size_t run_start = lo;
    Expr* f = expr_copy(fq ? llq_peek(fq, k) : f1);

    /* Key of element i-1, carried across iterations so f runs once per element. */
    Expr* prev_key = n ? splitby_key(f, llq_peek(q, lo)) : NULL;

    for (size_t i = lo + 1; i <= hi; i++) {
        bool cut = (i == hi);          /* the final run always closes at the end */
        Expr* key = NULL;
        if (!cut) {
            key = splitby_key(f, llq_peek(q, i));
            cut = !expr_sameq(key, prev_key);   /* SameQ, tolerance included */
        }
        if (key) { expr_free(prev_key); prev_key = key; }

        if (cut) {
            runs[num_runs++] = (k + 1 < nfs)
                ? splitby_level(q, run_start, i, fq, f1, k + 1, nfs)
                : llq_range(q, run_start, i);
            run_start = i;
        }
    }
    expr_free(prev_key);
    expr_free(f);

    Expr* result = llq_result(q, runs, num_runs);
    free(runs);
    return result;
}

Expr* builtin_splitby(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;

    Expr* list  = res->data.function.args[0];
    Expr* fspec = res->data.function.args[1];

    /* Atoms have no elements to split; leave the call unevaluated. */
    if (list->type != EXPR_FUNCTION && list->type != EXPR_LINKEDLIST) return NULL;

    /* fspec is either a single key function or a List (or LinkedList) of
     * them. An empty list specifies no split at all, which we leave
     * unevaluated rather than guessing at an interpretation. */
    LLSeq fq;
    bool multi = head_is(fspec, SYM_List) || fspec->type == EXPR_LINKEDLIST;
    size_t nfs = 1;
    if (multi) {
        llq_open(fspec, &fq);
        nfs = fq.n;
        if (nfs == 0) { llq_close(&fq); return NULL; }
    }

    /* A LinkedList splits from a snapshot of its slots into LinkedList runs. */
    LLSeq q;
    llq_open(list, &q);
    Expr* r = splitby_level(&q, 0, q.n, multi ? &fq : NULL, fspec, 0, nfs);
    llq_close(&q);
    if (multi) llq_close(&fq);
    return r;
}
