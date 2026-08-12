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
#include "assoc.h"    /* assoc_lookup_value — the single-key lookup primitive */

#define N_TRIALS   5
#define N_SMALL    20000
#define N_LARGE    40000     /* must be 2 * N_SMALL for the ratio to mean 2x */
/* An O(n) op doubles to ~2.0; allow generous noise headroom but stay well under
 * the 4.0 an accidental O(n^2) would produce. */
#define RATIO_MAX  3.3

/* Single-key gate: a FIXED number of single-key probes over an association of
 * size n vs 2n.  With the persistent key index each probe is O(1), so the total
 * is independent of n and the doubling ratio is ~1.0; the pre-index O(n) linear
 * scan made each probe O(n), so a regression back to scanning shows up as ratio
 * ~2.0.  The generous RATIO_MAX above would wave that through, so this half uses
 * its own tight threshold — the machine-independent proof that single-key
 * Lookup / KeyExistsQ / Part are genuinely O(1). */
#define SINGLE_KEY_REPS  1000000
#define SINGLE_KEY_RATIO_MAX 1.6

/* Loop-invariant gate: a FIXED number of Do-loop iterations, each looking a key
 * up in a large loop-invariant association, at size n vs 2n.  Do rebinds its
 * iterator every step, bumping the eval clock; pre-fix that invalidated the
 * association's cached fixed point so each iteration re-canonicalised it O(n)
 * and the whole loop was O(reps*n) -- ratio ~2.  The GROUND fixed-point flag
 * lets the loop-invariant value survive iterator churn, making the loop O(reps)
 * -- ratio ~1.  This is the direct regression guard for the eval-clock-churn
 * fix (distinct from the Map gate, which never binds a named iterator). */
#define DO_LOOP_REPS      20000
#define DO_LOOP_RATIO_MAX 1.6

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
    /* The calibration list. Deliberately Rationals, which no machine buffer can
     * hold, so a change in how ordinary Lists are STORED cannot move the divisor.
     * Total[keys<n>] was the calibration until automatic packed arrays made it
     * ~114x faster and inflated every normalized cost by the same factor, failing
     * all nine ops at once with nothing actually slower. Same landmine, and same
     * fix, as tests/bench_eval.c's. */
    snprintf(buf, sizeof(buf), "calib%d = Table[k/3, {k, %d}]", n, n);                 eval_discard(buf);
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

/* Each op has a scaling format (one %d for the size suffix) and an absolute-cost
 * baseline, expressed in "calibration units" -- its per-element time divided by
 * the per-element time of a plain list sum measured on the *same* run. That ratio
 * is machine-independent (a slow box scales both numerator and denominator), so
 * comparing it to a checked-in baseline gives a portable "much slower" tripwire.
 * baseline_norm is the value recorded on the reference machine; SLOWDOWN_MAX is
 * how many times worse than that we tolerate before failing. */
typedef struct { const char* label; const char* fmt; double baseline_norm; } Op;

/* Fail an op only if it is this many times slower (relative to the calibration)
 * than its recorded baseline -- loose enough to absorb machine/measurement
 * variation, tight enough to catch a real "it got much slower" regression. */
#define SLOWDOWN_MAX 2.5

/* baseline_norm recorded on an Apple M-series build, USE_ECM=OFF, 2026-07-30 --
 * median of three quiet runs. RE-RECORDED from the 2026-07-06 set because the
 * CALIBRATION changed, not because any operation did: it was Total[Range[n]],
 * which automatic packed arrays made ~114x faster, so every normalized cost
 * inflated by the same factor and all nine operations failed at once with
 * nothing actually slower. The divisor is now a Rational list, which no machine
 * buffer can hold. Run-to-run spread is ~15% on the association-heavy rows;
 * SLOWDOWN_MAX swamps it. */
static Op OPS[] = {
    { "Association @@ rules", "Association @@ rules%d", 5.5 },
    { "Counts",               "Counts[vals%d]", 0.5 },
    { "CountsBy",             "CountsBy[keys%d, EvenQ]", 1.4 },
    { "GroupBy",              "GroupBy[keys%d, Mod[#, 100] &]", 4.6 },
    { "Merge (Total)",        "Merge[{assoc%d, assoc2%d}, Total]", 24.0 },
    { "KeyUnion",             "KeyUnion[{assoc%d, assoc2%d}]", 18.6 },
    { "Lookup (bulk keys)",   "Lookup[assoc%d, keys%d]", 3.0 },
    { "Map over values",      "Map[# + 1 &, assoc%d]", 14.8 },
    { "KeySort",              "KeySort[assoc%d]", 3.9 },
};
#define N_OPS ((int)(sizeof(OPS) / sizeof(OPS[0])))

/* Evaluate `src`, returning its OWNED result (caller frees). */
static Expr* eval_to_expr(const char* src) {
    Expr* parsed = parse_expression(src);
    if (!parsed) { fprintf(stderr, "parse failed: %s\n", src); exit(2); }
    Expr* res = evaluate(parsed);
    expr_free(parsed);
    return res;
}

/* Median ns per single-key lookup on an n-entry association, measured on the
 * lookup PRIMITIVE (assoc_lookup_value) directly.
 *
 * This deliberately bypasses evaluate(): a Lookup[a, k] inside a Do/Table loop
 * re-evaluates the whole association O(n) every iteration (iterator binding
 * bumps the eval clock, invalidating the value's timestamp), so a loop measures
 * that re-evaluation, not the lookup.  The primitive is what "O(1) lookup"
 * actually means, and the index it lazily builds is exactly what a real caller
 * hits once the association is a stable value.  First call warms (builds the
 * index); the timed calls are the steady state. */
static double single_key_ns(int n, int reps, int hit) {
    char buf[128];
    snprintf(buf, sizeof(buf), "Association @@ Table[k -> Mod[k, 100], {k, %d}]", n);
    Expr* a = eval_to_expr(buf);
    Expr* key = expr_new_integer(hit ? (n / 2) : (-1));   /* present : absent */

    (void)assoc_lookup_value(a, key);                     /* warm: build index */

    double best = 1e30;
    volatile uintptr_t sink = 0;                          /* defeat DCE, portably */
    for (int t = 0; t < N_TRIALS; t++) {
        double t0 = now_ns();
        for (int i = 0; i < reps; i++)
            sink ^= (uintptr_t)assoc_lookup_value(a, key);
        double t1 = now_ns();
        double per = (t1 - t0) / (double)reps;
        if (per < best) best = per;
    }
    (void)sink;
    expr_free(key);
    expr_free(a);
    return best;
}

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
    double us_large_of[N_OPS];   /* remembered for the absolute-cost check below */
    char expr[256];
    for (int i = 0; i < N_OPS; i++) {
        format_op(expr, sizeof(expr), OPS[i].fmt, N_SMALL);
        double us_small = median_us(expr);
        format_op(expr, sizeof(expr), OPS[i].fmt, N_LARGE);
        double us_large = median_us(expr);
        us_large_of[i] = us_large;

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
    printf("PASS: all operations scaled linearly (ratio < %.1f)\n\n", RATIO_MAX);

    /* ---- Single-key O(1) gate -------------------------------------------------
     * The one that motivated the persistent index. Measured on the lookup
     * PRIMITIVE (assoc_lookup_value), size n vs 2n: O(1) => ratio ~1, the
     * pre-index O(n) scan => ratio ~2. See single_key_ns for why this bypasses
     * evaluate(). */
    printf("Single-key lookup primitive: %d probes on an association of size n vs 2n\n",
           SINGLE_KEY_REPS);
    printf("  ratio must stay < %.1f  (O(1) ~ 1.0; an O(n) scan regression ~ 2.0)\n\n",
           SINGLE_KEY_RATIO_MAX);
    printf("%-14s %12s %12s %6s\n", "case", "n=20k(ns)", "n=40k(ns)", "ratio");
    printf("--------------------------------------------------------------------\n");
    int sk_fail = 0;
    const struct { const char* label; int hit; } sk_cases[] = { {"hit", 1}, {"miss", 0} };
    for (int i = 0; i < 2; i++) {
        double ns_small = single_key_ns(N_SMALL, SINGLE_KEY_REPS, sk_cases[i].hit);
        double ns_large = single_key_ns(N_LARGE, SINGLE_KEY_REPS, sk_cases[i].hit);
        double ratio = (ns_small > 0.0) ? ns_large / ns_small : 0.0;
        int bad = (ratio > SINGLE_KEY_RATIO_MAX);
        printf("%-14s %12.2f %12.2f %6.2f%s\n",
               sk_cases[i].label, ns_small, ns_large, ratio,
               bad ? "  <== O(n) SCAN?" : "");
        if (bad) sk_fail++;
    }
    printf("--------------------------------------------------------------------\n");
    if (sk_fail) {
        printf("FAIL: %d single-key case(s) scaled with n (ratio > %.1f) -- lookup is not O(1)\n",
               sk_fail, SINGLE_KEY_RATIO_MAX);
        return 1;
    }
    printf("PASS: single-key lookup is O(1) (ratio < %.1f)\n", SINGLE_KEY_RATIO_MAX);

    /* End-to-end interpreter check: Total[Map[Lookup[a, #]&, keys]] over an
     * association of size n vs 2n. Map does not bind a named iterator, so the
     * eval clock is stable across its elements; evaluate()'s in-loop timestamp
     * short-circuit then keeps `a` from being re-canonicalised each element, and
     * the lazily-built key index makes every probe O(1). Ratio ~1 proves the
     * whole interpreter path (not just the C primitive) is O(1); a regression in
     * either the index or the short-circuit shows up as ratio ~2. */
    {
        char e[96];
        snprintf(e, sizeof(e), "Total[Map[Lookup[assoc%d, #] &, keys%d]]", N_SMALL, N_SMALL);
        double us_small = median_us(e);
        snprintf(e, sizeof(e), "Total[Map[Lookup[assoc%d, #] &, keys%d]]", N_LARGE, N_LARGE);
        /* keys<N_LARGE> is Range[N_LARGE]; the map is over 2x as many keys, so
         * normalise the ratio by the key count to isolate per-probe scaling. */
        double us_large = median_us(e);
        double ratio = (us_small > 0.0) ? (us_large / 2.0) / us_small : 0.0;
        /* A looser bound than the C primitive: this path also pays evaluator and
         * list-building overhead with more run-to-run variance. It still cleanly
         * separates O(1) (~1.2) from an O(n) re-canonicalisation regression (~2). */
        const double interp_max = 1.8;
        int bad = (ratio > interp_max);
        printf("Interpreter Map[Lookup]: n=%d %.0fus, n=%d %.0fus (per-key ratio %.2f)%s\n",
               N_SMALL, us_small, N_LARGE, us_large, ratio, bad ? "  <== NOT O(1)" : "");
        if (bad) {
            printf("FAIL: interpreter repeated Lookup is not O(1) (per-key ratio > %.1f)\n",
                   interp_max);
            return 1;
        }
    }
    printf("PASS: interpreter repeated Lookup is O(1)\n\n");

    /* ---- Loop-invariant O(1) gate (the eval-clock-churn fix) -----------------
     * The Map gate above never binds a named iterator, so the eval clock is
     * stable across its elements.  This gate uses a Do loop, which DOES rebind
     * its iterator every step and bumps the clock -- the exact churn that made a
     * loop-invariant association re-canonicalise O(n) per iteration before the
     * GROUND fixed-point flag.  A fixed rep count timed at n vs 2n gives ratio
     * ~1 when the value survives the churn (O(reps)) and ~2 if it does not
     * (O(reps*n)).  `assoc<n>` is the global set up above; the key varies
     * (exercises the index) while the association stays invariant. */
    {
        char e[160];
        snprintf(e, sizeof(e),
            "Module[{s = 0}, Do[s = s + Lookup[assoc%d, Mod[i, 100] + 1, 0], {i, 1, %d}]; s]",
            N_SMALL, DO_LOOP_REPS);
        double us_small = median_us(e);
        snprintf(e, sizeof(e),
            "Module[{s = 0}, Do[s = s + Lookup[assoc%d, Mod[i, 100] + 1, 0], {i, 1, %d}]; s]",
            N_LARGE, DO_LOOP_REPS);
        double us_large = median_us(e);
        double ratio = (us_small > 0.0) ? us_large / us_small : 0.0;
        int bad = (ratio > DO_LOOP_RATIO_MAX);
        printf("Do-loop Lookup (loop-invariant assoc): n=%d %.0fus, n=%d %.0fus (ratio %.2f)%s\n",
               N_SMALL, us_small, N_LARGE, us_large, ratio, bad ? "  <== NOT O(1)" : "");
        if (bad) {
            printf("FAIL: Do-loop over a loop-invariant association is not O(1) "
                   "(ratio > %.1f) -- eval-clock churn re-canonicalises it each iteration\n",
                   DO_LOOP_RATIO_MAX);
            return 1;
        }
    }
    printf("PASS: Do-loop over a loop-invariant association is O(1)\n\n");

    /* ---- Absolute-cost check (machine-normalized) ----------------------------
     * Calibrate against a plain list sum of the same size, then express each op's
     * per-element cost as a multiple of that. This "cost in calibration units" is
     * machine-independent, so we can gate on it: fail if an op is > SLOWDOWN_MAX x
     * its recorded baseline. */
    char calib_expr[64];
    snprintf(calib_expr, sizeof(calib_expr), "Total[calib%d]", N_LARGE);
    double calib_us = median_us(calib_expr);
    if (calib_us <= 0.0) calib_us = 1e-6;    /* guard against a zero divide */

    printf("Absolute cost vs. calibration (Total[Rational list] of the same size)\n");
    printf("  fail if an op is > %.1fx its recorded baseline cost\n\n", SLOWDOWN_MAX);
    printf("%-22s %10s %10s %8s\n", "operation", "norm", "baseline", "x base");
    printf("--------------------------------------------------------------------\n");

    int slow = 0;
    for (int i = 0; i < N_OPS; i++) {
        double norm = us_large_of[i] / calib_us;   /* cost in calibration units */
        if (OPS[i].baseline_norm <= 0.0) {
            /* Baseline not yet recorded: report the measured value, do not gate. */
            printf("%-22s %10.2f %10s %8s   (record as baseline)\n",
                   OPS[i].label, norm, "-", "-");
            continue;
        }
        double x = norm / OPS[i].baseline_norm;
        int bad = (x > SLOWDOWN_MAX);
        printf("%-22s %10.2f %10.2f %8.2f%s\n",
               OPS[i].label, norm, OPS[i].baseline_norm, x,
               bad ? "  <== MUCH SLOWER" : "");
        if (bad) slow++;
    }

    printf("--------------------------------------------------------------------\n");
    if (slow) {
        printf("FAIL: %d operation(s) more than %.1fx slower than baseline\n",
               slow, SLOWDOWN_MAX);
        return 1;
    }
    printf("PASS: no operation exceeded %.1fx its baseline cost\n", SLOWDOWN_MAX);
    return 0;
}
