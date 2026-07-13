/*
 * tests/bench_eval.c
 *
 * Baseline benchmark for the core evaluator (evaluate() / evaluate_step()).
 *
 * Phase 0 of EVAL_SYMTAB_IMPROVEMENTS.md. The upcoming symbol-object rework
 * (precomputed base attributes, unified intern+symtab table, EXPR_SYMBOL that
 * caches its SymbolDef*) targets the per-function-node dispatch tax paid on
 * every evaluation pass:
 *     get_attributes()      -- linear strcmp over a ~140-entry table
 *     symtab_get_def() x3   -- attributes, DownValues, builtin dispatch
 * This harness measures workloads where that tax is a meaningful fraction of
 * evaluate() time, so Phase 1 (base attrs) and Phase 3 (cached def) can be
 * validated against a committed baseline on the *same* machine.
 *
 * Design mirrors bench_assoc.c: median-of-trials wall time, plus a
 * machine-normalized cost (each workload's time divided by a plain
 * Total[Range[...]] measured on the same run) so a checked-in baseline gives a
 * portable "it got much slower" tripwire. Absolute microseconds are printed for
 * eyeballing but never gate; the normalized ratio does.
 *
 * FIRST RUN: every baseline_norm is 0.0, so every row prints "(record as
 * baseline)" and the run passes. Copy the printed `norm` values into
 * baseline_norm below to arm the gate, and re-record after an *intended*
 * evaluator performance change (that is the whole point -- Phases 1/3 should
 * move these down).
 *
 * Wired into ctest as `bench_eval`. Runtime is a second or two.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"

#define N_TRIALS 5

static double now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1.0e9 + (double)ts.tv_nsec;
}

static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

/* Evaluate a source string once for its side effect (e.g. a symbol
 * definition), discarding the result. */
static void eval_discard(const char* src) {
    Expr* parsed = parse_expression(src);
    if (!parsed) { fprintf(stderr, "parse failed: %s\n", src); exit(2); }
    Expr* res = evaluate(parsed);
    expr_free(res);
    expr_free(parsed);
}

/* Median wall time (microseconds) of parse+evaluate of `expr` over N_TRIALS,
 * after one untimed warm-up so first-touch symbol interning / table growth does
 * not skew the median. Each trial re-parses to get a fresh (uncached) tree, so
 * we measure a full evaluation rather than an eval-clock cache hit. */
static double median_us(const char* expr) {
    { Expr* p = parse_expression(expr);
      if (!p) { fprintf(stderr, "parse failed: %s\n", expr); exit(2); }
      Expr* r = evaluate(p); expr_free(r); expr_free(p); }
    double s[N_TRIALS];
    for (int t = 0; t < N_TRIALS; t++) {
        Expr* p = parse_expression(expr);
        double t0 = now_ns();
        Expr* r = evaluate(p);
        double t1 = now_ns();
        expr_free(r);
        expr_free(p);
        s[t] = (t1 - t0) / 1000.0;
    }
    qsort(s, N_TRIALS, sizeof(double), cmp_double);
    return s[N_TRIALS / 2];
}

/* Each workload is chosen to stress the evaluator dispatch path with many small
 * function nodes and/or many evaluation passes, rather than a single heavy
 * builtin body (which would swamp the dispatch signal we care about). */
typedef struct { const char* label; const char* expr; double baseline_norm; } Bench;

/* Fail a workload only if it is this many times slower (relative to the
 * calibration) than its recorded baseline -- loose enough to absorb machine and
 * measurement variation, tight enough to catch a real regression. */
#define SLOWDOWN_MAX 2.5

/* baseline_norm recorded on an Apple M-series build, USE_ECM=OFF. Stable to a
 * few percent across runs; the SLOWDOWN_MAX margin swamps that. Re-record after
 * each intended speedup so the gate keeps catching accidental regressions on the
 * new, lower floor.
 *
 * History (norm on the pure-dispatch row, the primary dispatch signal):
 *   2026-07-13  Phase 0 (before rework)                  0.60
 *   2026-07-13  Phase 1 (base-attr fold in get_attributes) 0.40  (-33%)
 * Values below are the post-Phase-1 floor. The dispatch-light rows (orderless
 * sort, replacerepeated fold) are matcher/compare-bound and did not move; their
 * baselines are unchanged. */
static Bench BENCHES[] = {
    /* 8000 structurally DISTINCT undefined nodes hh[k, k] (k differs each row),
     * so no eval-clock cache hit or refcount DAG-sharing collapses them. No rule
     * fires, so every node pays pure dispatch: head eval + get_attributes +
     * apply_down_values miss + builtin-lookup miss. The cleanest dispatch signal
     * in the set. (A Nest[f[u,u]] tree would instead share args by refcount and
     * collapse to a ~13-node DAG -- not a dispatch test.) */
    { "pure dispatch (8k distinct)", "Length[Table[hh[k, k], {k, 8000}]]", 0.42 },

    /* 4000 sequential DownValue rewrites g[k] -> k+1. Exercises
     * apply_down_values + get_attributes on the user symbol g every step. */
    { "downvalue rewrite (Nest 4k)", "Nest[g, 0, 4000]", 0.16 },

    /* Listable threading of Cos over a 4000-element list: 4000 threaded
     * sub-evaluations, each dispatching Cos[k] (stays symbolic). */
    { "listable thread (Cos 4k)", "Length[Cos[Range[4000]]]", 0.13 },

    /* Build a 1200-term Plus and let Orderless sort it: stresses eval_sort_args
     * + expr_compare in addition to per-term dispatch. */
    { "orderless sort (Plus 1200)", "Length[Plus @@ Table[c[k] x, {k, 1200}]]", 0.92 },

    /* Pairwise fold of a 400-element list to a fixed point via ReplaceRepeated:
     * ~399 rewrites, each a full matcher + re-evaluation pass (O(n^2), so kept
     * small). Stresses the pattern matcher, not just dispatch. */
    { "replacerepeated fold (400)", "First[Range[400] //. {a_, b_, r___} :> {a + b, r}]", 3.00 },

    /* 3000 sequential Plus re-evaluations (each Nest step re-evaluates the
     * running sum). Exercises Flat/Orderless combination on every pass. */
    { "nested plus collapse (3k)", "Nest[# + a &, x, 3000] /. a -> 0", 0.37 },
};
#define N_BENCH ((int)(sizeof(BENCHES) / sizeof(BENCHES[0])))

int main(void) {
    symtab_init();
    core_init();

    /* Untimed setup: the one workload that needs a definition. */
    eval_discard("g[x_] := x + 1");

    printf("Evaluator baseline benchmark (median of %d trials)\n", N_TRIALS);
    printf("  norm = workload_us / calibration_us (machine-independent)\n");
    printf("  gate: fail if norm > %.1fx its recorded baseline\n\n", SLOWDOWN_MAX);

    /* Calibration: a plain arithmetic reduction over a list. Cheap, dispatch-
     * light, and its wall time tracks raw machine speed, so dividing by it
     * cancels most machine/measurement variance. */
    double calib_us = median_us("Total[Range[40000]]");
    if (calib_us <= 0.0) calib_us = 1e-6;
    printf("calibration  Total[Range[40000]]  = %.1f us\n\n", calib_us);

    printf("%-30s %12s %10s %10s %8s\n",
           "workload", "median(us)", "norm", "baseline", "x base");
    printf("--------------------------------------------------------------------------------\n");

    int slow = 0;
    for (int i = 0; i < N_BENCH; i++) {
        double us = median_us(BENCHES[i].expr);
        double norm = us / calib_us;
        if (BENCHES[i].baseline_norm <= 0.0) {
            printf("%-30s %12.1f %10.2f %10s %8s   (record as baseline)\n",
                   BENCHES[i].label, us, norm, "-", "-");
            continue;
        }
        double x = norm / BENCHES[i].baseline_norm;
        int bad = (x > SLOWDOWN_MAX);
        printf("%-30s %12.1f %10.2f %10.2f %8.2f%s\n",
               BENCHES[i].label, us, norm, BENCHES[i].baseline_norm, x,
               bad ? "  <== MUCH SLOWER" : "");
        if (bad) slow++;
    }

    printf("--------------------------------------------------------------------------------\n");
    if (slow) {
        printf("FAIL: %d workload(s) more than %.1fx slower than baseline\n",
               slow, SLOWDOWN_MAX);
        return 1;
    }
    printf("PASS: no workload exceeded %.1fx its baseline cost\n", SLOWDOWN_MAX);
    return 0;
}
