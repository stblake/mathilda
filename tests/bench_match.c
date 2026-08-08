/*
 * tests/bench_match.c
 *
 * Performance regression gate for the pattern matcher (src/match.c).
 *
 * The failure mode that matters is an *algorithmic* regression -- a per-element
 * matcher path that was O(n) silently becoming O(n^2) (an accidental rescan in
 * the traversal, or a matcher fast-path that stopped firing). Absolute wall
 * times are machine-dependent and make flaky gates, but the **doubling ratio**
 * t(2n)/t(n) is machine-independent: ~2 for an O(n) op, ~4 for an O(n^2) one.
 * Each op is timed at n and 2n and fails if the ratio exceeds RATIO_MAX --
 * comfortably above 2 (noise headroom) yet well below 4 (accidental quadratic).
 *
 * The fixtures are lists of f[k] (function nodes), NOT Range[n]: a packed
 * integer array would send Cases/Count/Position/Replace down a buffer fast path
 * that never touches the matcher, which is precisely the code under test. f[k]
 * cannot pack, so every element is matched structurally.
 *
 * The ns/element column is informational (constant-factor drift for a human);
 * only the doubling ratio gates. Wired into ctest as `bench_match`.
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
#define N_SMALL    10000
#define N_LARGE    20000     /* must be 2 * N_SMALL for the ratio to mean 2x */
/* O(n) doubles to ~2.0; generous noise headroom, well under the 4.0 an
 * accidental O(n^2) would produce. */
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

static void eval_discard(const char* src) {
    Expr* parsed = parse_expression(src);
    if (!parsed) { fprintf(stderr, "parse failed: %s\n", src); exit(2); }
    Expr* res = evaluate(parsed);
    expr_free(res);
    expr_free(parsed);
}

/* fs<n> = {f[1], f[2], ..., f[n]} -- function nodes, so the list cannot pack and
 * every Cases/Count/Replace element is matched by the matcher, not a buffer. */
static void setup(int n) {
    char buf[128];
    snprintf(buf, sizeof(buf), "fs%d = Table[f[k], {k, %d}]", n, n);   eval_discard(buf);
}

/* Median wall time (us) of evaluating `expr` N_TRIALS times, after one untimed
 * warm-up so caching does not skew the median. */
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

/* Each op is a traversal that matches one pattern per element -- genuinely O(n)
 * and dominated by the matcher. The one %d is the size suffix. */
typedef struct { const char* label; const char* fmt; } Op;

static Op OPS[] = {
    { "Cases f[_]",          "Length[Cases[fs%d, f[_]]]" },
    { "Count f[_?OddQ]",     "Count[fs%d, f[_?OddQ]]" },
    { "Position f[_]",       "Length[Position[fs%d, f[_]]]" },
    { "DeleteCases f[_?Ev]", "Length[DeleteCases[fs%d, f[_?EvenQ]]]" },
    { "ReplaceAll f[x_]",    "Length[fs%d /. f[x_] :> x]" },
    { "MemberQ f[_?Prime]",  "MemberQ[fs%d, f[_?PrimeQ]]" },
};
#define N_OPS ((int)(sizeof(OPS) / sizeof(OPS[0])))

int main(void) {
    symtab_init();
    core_init();

    /* A DownValue whose dispatch we also scale-test below. */
    eval_discard("g[x_] := x");

    printf("Pattern-matcher performance gate (median of %d trials)\n", N_TRIALS);
    printf("  scaling: t(%d)/t(%d) must stay < %.1f  (O(n) ~ 2.0, O(n^2) ~ 4.0)\n\n",
           N_LARGE, N_SMALL, RATIO_MAX);

    setup(N_SMALL);
    setup(N_LARGE);

    printf("%-22s %10s %10s %8s %6s\n", "operation", "n=10k(us)", "n=20k(us)", "ns/elem", "ratio");
    printf("--------------------------------------------------------------------\n");

    int failures = 0;
    char expr[128];
    for (int i = 0; i < N_OPS; i++) {
        snprintf(expr, sizeof(expr), OPS[i].fmt, N_SMALL);
        double us_small = median_us(expr);
        snprintf(expr, sizeof(expr), OPS[i].fmt, N_LARGE);
        double us_large = median_us(expr);

        double ratio = (us_small > 0.0) ? us_large / us_small : 0.0;
        double ns_per_elem = (us_large * 1000.0) / (double)N_LARGE;
        int bad = (ratio > RATIO_MAX);
        printf("%-22s %10.1f %10.1f %8.1f %6.2f%s\n",
               OPS[i].label, us_small, us_large, ns_per_elem, ratio,
               bad ? "  <== REGRESSION" : "");
        if (bad) failures++;
    }

    /* DownValue dispatch: Map[g, fs<n>] fires g[x_]:=x per element. Separate row
     * because it exercises apply_down_values dispatch, not the Cases/Replace
     * level-walker. */
    {
        char e[64];
        snprintf(e, sizeof(e), "Length[Map[g, fs%d]]", N_SMALL);
        double us_small = median_us(e);
        snprintf(e, sizeof(e), "Length[Map[g, fs%d]]", N_LARGE);
        double us_large = median_us(e);
        double ratio = (us_small > 0.0) ? us_large / us_small : 0.0;
        double ns_per_elem = (us_large * 1000.0) / (double)N_LARGE;
        int bad = (ratio > RATIO_MAX);
        printf("%-22s %10.1f %10.1f %8.1f %6.2f%s\n",
               "DownValue g[x_] Map", us_small, us_large, ns_per_elem, ratio,
               bad ? "  <== REGRESSION" : "");
        if (bad) failures++;
    }

    printf("--------------------------------------------------------------------\n");
    if (failures) {
        printf("FAIL: %d matcher op(s) scaled worse than O(n) (ratio > %.1f)\n",
               failures, RATIO_MAX);
        return 1;
    }
    printf("PASS: all matcher ops scaled linearly (ratio < %.1f)\n", RATIO_MAX);
    return 0;
}
