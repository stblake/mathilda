/*
 * tests/bench_assoc.c
 *
 * Performance regression gate for the Association data structure.
 *
 * The thing we actually want to catch is *algorithmic* regression -- an
 * association operation that was O(n) silently becoming O(n^2) (a per-element
 * linear scan creeping back in where a hash lookup belonged). Absolute wall
 * times are machine-dependent and make flaky CI gates, but the **doubling
 * ratio** t(2n)/t(n) is machine-independent: it is ~2 for an O(n) op and ~4 for
 * an O(n^2) one. So this harness times each op at n and 2n and fails if the
 * ratio exceeds RATIO_MAX -- comfortably above 2 (noise headroom) yet well below
 * 4 (an accidental quadratic).
 *
 * It also prints the absolute median ns/element so a human can eyeball whether
 * constant factors are drifting up over time; those numbers are informational
 * and never fail the run.
 *
 * Wired into ctest as `bench_assoc`. Runtime is a few hundred ms.
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

#define N_TRIALS   5
#define N_SMALL    20000
#define N_LARGE    40000     /* must be 2 * N_SMALL for the ratio to mean 2x */
/* An O(n) op doubles to ~2.0; allow generous noise headroom but stay well under
 * the 4.0 an accidental O(n^2) would produce. */
#define RATIO_MAX  3.3

static double now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1.0e9 + (double)ts.tv_nsec;
}

static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

/* Evaluate a source string once for its side effect (e.g. a symbol assignment),
 * discarding the result. */
static void eval_discard(const char* src) {
    Expr* parsed = parse_expression(src);
    if (!parsed) { fprintf(stderr, "parse failed: %s\n", src); exit(2); }
    Expr* res = evaluate(parsed);
    expr_free(res);
    expr_free(parsed);
}

/* Build the per-size fixtures as symbols suffixed by n:
 *   rulesN   = {1 -> 0, 2 -> 1, ..., n -> Mod[n,100]}   (list of rules)
 *   valsN    = {Mod[1,100], ...}                        (values, ~100 distinct)
 *   keysN    = {1, ..., n}
 *   assocN   = <|k -> Mod[k,100]|>
 *   assoc2N  = <|k -> k|> over an overlapping key range   (for Merge/KeyUnion)
 * Built once, outside every timing block, so we measure the op and not setup. */
static void setup(int n) {
    char buf[256];
    snprintf(buf, sizeof(buf), "rules%d = Table[k -> Mod[k, 100], {k, %d}]", n, n);   eval_discard(buf);
    snprintf(buf, sizeof(buf), "vals%d = Table[Mod[k, 100], {k, %d}]", n, n);         eval_discard(buf);
    snprintf(buf, sizeof(buf), "keys%d = Range[%d]", n, n);                            eval_discard(buf);
    snprintf(buf, sizeof(buf), "assoc%d = Association @@ rules%d", n, n);             eval_discard(buf);
    snprintf(buf, sizeof(buf), "assoc2%d = Association @@ Table[k -> k, {k, %d, %d}]",
             n, n / 2, n + n / 2);                                                     eval_discard(buf);
}

/* Median wall time (microseconds) of evaluating `expr` N_TRIALS times, after one
 * untimed warm-up so symbol-table / cache priming does not skew the median. */
static double median_us(const char* expr) {
    { Expr* p = parse_expression(expr); Expr* r = evaluate(p); expr_free(r); expr_free(p); }
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

/* Each op is a printf format with one %d for the size suffix. */
typedef struct { const char* label; const char* fmt; } Op;

static const Op OPS[] = {
    { "Association @@ rules", "Association @@ rules%d" },
    { "Counts",               "Counts[vals%d]" },
    { "CountsBy",             "CountsBy[keys%d, EvenQ]" },
    { "GroupBy",              "GroupBy[keys%d, Mod[#, 100] &]" },
    { "Merge (Total)",        "Merge[{assoc%d, assoc2%d}, Total]" },
    { "KeyUnion",             "KeyUnion[{assoc%d, assoc2%d}]" },
    { "Lookup (bulk keys)",   "Lookup[assoc%d, keys%d]" },
    { "Map over values",      "Map[# + 1 &, assoc%d]" },
    { "KeySort",              "KeySort[assoc%d]" },
};
#define N_OPS ((int)(sizeof(OPS) / sizeof(OPS[0])))

/* Fill `out` with the op expression for size n (handles the 1- and 2-%d fmts). */
static void format_op(char* out, size_t cap, const char* fmt, int n) {
    /* All fmts use the same n for every %d, so a plain vsnprintf-style call with
     * the size repeated a few times is safe (extra args are ignored). */
    snprintf(out, cap, fmt, n, n, n);
}

int main(void) {
    symtab_init();
    core_init();

    printf("Association performance gate (median of %d trials)\n", N_TRIALS);
    printf("  scaling: t(%d)/t(%d) must stay < %.1f  (O(n) ~ 2.0, O(n^2) ~ 4.0)\n\n",
           N_LARGE, N_SMALL, RATIO_MAX);

    setup(N_SMALL);
    setup(N_LARGE);

    printf("%-22s %10s %10s %8s %6s\n", "operation", "n=" "20k(us)", "n=40k(us)", "ns/elem", "ratio");
    printf("--------------------------------------------------------------------\n");

    int failures = 0;
    char expr[256];
    for (int i = 0; i < N_OPS; i++) {
        format_op(expr, sizeof(expr), OPS[i].fmt, N_SMALL);
        double us_small = median_us(expr);
        format_op(expr, sizeof(expr), OPS[i].fmt, N_LARGE);
        double us_large = median_us(expr);

        double ratio = (us_small > 0.0) ? us_large / us_small : 0.0;
        double ns_per_elem = (us_large * 1000.0) / (double)N_LARGE;
        /* KeySort is O(n log n): its doubling ratio is a touch above 2, still
         * well under RATIO_MAX. Everything else is hash-backed O(n). */
        int bad = (ratio > RATIO_MAX);
        printf("%-22s %10.1f %10.1f %8.1f %6.2f%s\n",
               OPS[i].label, us_small, us_large, ns_per_elem, ratio,
               bad ? "  <== REGRESSION" : "");
        if (bad) failures++;
    }

    printf("--------------------------------------------------------------------\n");
    if (failures) {
        printf("FAIL: %d operation(s) scaled worse than O(n) (ratio > %.1f)\n",
               failures, RATIO_MAX);
        return 1;
    }
    printf("PASS: all operations scaled linearly (ratio < %.1f)\n", RATIO_MAX);
    return 0;
}
