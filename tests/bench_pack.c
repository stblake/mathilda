/*
 * tests/bench_pack.c
 *
 * Performance gate for the automatic packed-array surface
 * (docs/design/packed_arrays.md).
 *
 * WHY THIS EXISTS. The packing gate's rule is "materialise for any head that has
 * not opted in", so a missing opt-in is CORRECT and SILENT: the right answer, on
 * the slow path, with every test still green. That failure mode has now happened
 * four times, and each one cost between 30x and 658x:
 *
 *     the 26 linear-algebra heads   LinearSolve of a 90x90 real system
 *                                   stack-overflowed (the same test selected
 *                                   both the fast path and the algorithm)
 *     Nest / NestList / FixedPoint  118x -- Fold was aware, Nest was not
 *     the structural family         30x-237x each (RotateLeft, Join, Partition,
 *                                   Differences, Riffle, PadLeft, PadRight)
 *     user DownValues               158x -- `jac[u_] := ...` materialised on
 *                                   every call
 *
 * None was found by reading code or by running the suite. They were found by
 * timing one operation at a time against Mathematica. This file is the standing
 * version of that: every row below is an operation whose slow path is at least
 * an order of magnitude worse, so a regression cannot hide inside the 2.5x
 * threshold.
 *
 * Design mirrors bench_eval.c and bench_assoc.c: median-of-trials wall time plus
 * a machine-normalized cost (workload / calibration measured on the same run),
 * so a checked-in baseline is a portable "it got much slower" tripwire. Absolute
 * microseconds are printed for eyeballing but never gate.
 *
 * FIRST RUN: every baseline_norm is 0.0, so every row prints
 * "(record as baseline)" and the run passes. Copy the printed `norm` values into
 * baseline_norm to arm the gate, and re-record only after an INTENDED speedup.
 *
 * Wired into ctest as `bench_pack`.
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
#include "pack.h"

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

static void eval_discard(const char* src) {
    Expr* parsed = parse_expression(src);
    if (!parsed) { fprintf(stderr, "parse failed: %s\n", src); exit(2); }
    Expr* res = evaluate(parsed);
    expr_free(res);
    expr_free(parsed);
}

/* Median wall time (us) of parse+evaluate over N_TRIALS, after one untimed
 * warm-up. Re-parses per trial so each is a full evaluation, not a cache hit. */
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

typedef struct { const char* label; const char* expr; double baseline_norm; } Bench;

#define SLOWDOWN_MAX 2.5

/* baseline_norm: median of three runs, 2026-07-30, Intel i9-9880H, USE_ECM=ON.
 * Run-to-run spread is a few percent on most rows and up to ~30% on Det and
 * Inverse, all comfortably inside the 2.5x threshold.
 *
 * The first armed run of this file paid for itself: `Differences int64` read
 * 89.3 ms against 0.83 ms for the same-size float64 buffer. The int64 arm of
 * ndstruct_differences tested `!ci_sub_i64(...)`, but ci_sub_i64 is
 * __builtin_sub_overflow and returns TRUE on OVERFLOW -- so the loop abandoned
 * on the first SUCCESSFUL subtraction and every call fell back to
 * delist_repack. Correct answer, correct dtype, 1193x too slow, and invisible
 * to every value test. That is exactly the failure mode this file exists for.
 *
 * Each row's SLOW path (the one this gate exists to catch) is noted, so that a
 * failure is self-diagnosing: if a row trips, the named opt-in is what to check
 * in src/pack.c's AWARE / INT64_OK lists. */
static Bench BENCHES[] = {
    /* ---- the structural family: native buffer walks in src/ndstruct.c ----
     * Slow path for all of these is ndstruct_delist_repack -- one Expr per
     * element, the generic List implementation, then a re-sniff and re-pack. */
    { "RotateLeft matrix (200^2)",   "Length[RotateLeft[pm, {1, 0}]]",        0.0150 },
    { "RotateLeft vector (10^5)",    "Length[RotateLeft[pv, 3]]",             0.0348 },
    { "Reverse matrix (200^2)",      "Length[Reverse[pm]]",                   0.0089 },
    { "Join two vectors (10^5)",     "Length[Join[pv, pv]]",                  0.0342 },
    { "Partition (10^5, 2)",         "Length[Partition[pv, 2]]",              0.0177 },
    { "Differences (10^5)",          "Length[Differences[pv]]",               0.4717 },
    { "Riffle Real sep (10^5)",      "Length[Riffle[pv, 0.]]",                0.3850 },
    { "PadRight Real fill (10^5)",   "Length[PadRight[pv, 100001, 0.]]",      0.0177 },

    /* ---- the iteration family (Nest et al.) ----
     * Slow path: the gate materialises the state on EVERY iteration, so the cost
     * is O(iterations x elements) instead of O(elements). Was 118x. */
    { "Nest over packed (20x)",      "Length[Nest[Reverse, pm, 20]]",         0.1593 },
    { "FoldList over packed (10^4)", "Length[FoldList[Plus, 0., fv]]",        2.2774 },

    /* ---- user DownValues that bind opaquely (dv_binds_opaquely, src/eval.c) --
     * Slow path: a user symbol has no packed_aware bit, so `jac[u_]` used to
     * materialise on every call. Was 158x on the composite stencil below. */
    { "user f[v_] over packed",      "Length[usr[pv]]",                       0.0238 },
    { "Jacobi stencil (200^2, 20)",  "Length[Nest[jac, pm, 20]]",             3.1153 },

    /* ---- linear algebra (the 26 heads; ndla_* in src/linalg/ndlinalg.c) ----
     * Slow path is not merely slow: with the buffer materialised first,
     * linalg_call_has_ndarray reads false and a machine-real solve takes the
     * exact fraction-free path, which recurses. This row is a crash canary. */
    { "LinearSolve (120^2)",         "Length[LinearSolve[sm, sv]]",           0.1041 },
    { "Det (120^2)",                 "Head[Det[sm]]",                         0.1305 },
    { "Inverse (120^2)",             "Length[Inverse[sm]]",                   0.3974 },

    /* ---- reductions and order statistics already on the buffer ----
     * Regression guards for paths that are currently fast. */
    { "Total (10^5)",                "Total[pv]",                             0.0231 },
    { "Sort (10^5)",                 "Length[Sort[pv]]",                      1.7303 },
    { "Accumulate (10^5)",           "Length[Accumulate[pv]]",                0.5034 },

    /* ---- exact int64 buffers ----
     * Slow path: any head not on INT64_OK materialises an integer buffer. The
     * answers must stay exact Integers either way, so only the timing tells. */
    { "Total int64 (10^5)",          "Head[Total[iv]]",                       0.0225 },
    { "RotateLeft int64 (10^5)",     "Length[RotateLeft[iv, 3]]",             0.0348 },
    { "Differences int64 (10^5)",    "Length[Differences[iv]]",               0.0489 },
};
#define N_BENCH ((int)(sizeof(BENCHES) / sizeof(BENCHES[0])))

int main(void) {
    symtab_init();
    core_init();

    /* Untimed setup. Every value here is above PACK_MIN_ELEMENTS, so it is a
     * packed buffer; the point of each row is that it STAYS one. */
    eval_discard("pv = Range[1., 100000.]");
    eval_discard("iv = Range[100000]");
    /* FoldList applies its function per element through the evaluator, so it
     * is inherently O(n) evaluations rather than a buffer walk; kept small so
     * it does not dominate this benchmark's runtime. */
    eval_discard("fv = Range[1., 10000.]");
    eval_discard("pm = Table[N[i + j], {i, 200}, {j, 200}]");
    eval_discard("sm = Table[N[Mod[7 i + 13 j, 101] + If[i == j, 200, 0]], "
                 "{i, 120}, {j, 120}]");
    eval_discard("sv = Table[N[Mod[11 k, 97]], {k, 120}]");
    eval_discard("usr[v_] := Total[v] + Length[v]");
    eval_discard("jac[u_] := (RotateLeft[u, {1, 0}] + RotateRight[u, {1, 0}] + "
                 "RotateLeft[u, {0, 1}] + RotateRight[u, {0, 1}])/4.");

    /* Assert the inputs really are packed. Without this the whole file could
     * pass while measuring nothing: if packing were switched off, every row
     * would simply be slow *and consistent*, and the normalized baselines would
     * be re-recorded against the slow path. */
    /* `sv` is deliberately NOT in this list. At 120 elements it is below
     * PACK_MIN_ELEMENTS, so LinearSolve[sm, sv] arrives with the matrix packed
     * (14400 elements) and the right-hand side not -- which is exactly the
     * configuration that stack-overflowed before nd_load_rhs learned to accept a
     * plain numeric List. Making it pack would remove the canary. */
    static const char* const MUST_PACK[] = { "pv", "iv", "fv", "pm", "sm" };
    for (size_t i = 0; i < sizeof(MUST_PACK) / sizeof(MUST_PACK[0]); i++) {
        char buf[64];
        snprintf(buf, sizeof buf, "NDArrayQ[%s]", MUST_PACK[i]);
        Expr* p = parse_expression(buf);
        Expr* r = evaluate(p);
        int ok = r && r->type == EXPR_SYMBOL && strcmp(r->data.symbol.name, "True") == 0;
        expr_free(r); expr_free(p);
        if (!ok) {
            fprintf(stderr,
                    "FAIL: %s is not packed -- this benchmark would measure the "
                    "unpacked path and silently re-baseline against it.\n"
                    "(MATHILDA_NO_PACK set, or $AutoArrayPacking False?)\n",
                    MUST_PACK[i]);
            return 2;
        }
    }

    printf("Packed-array performance gate (median of %d trials)\n", N_TRIALS);
    printf("  norm = workload_us / calibration_us (machine-independent)\n");
    printf("  gate: fail if norm > %.1fx its recorded baseline\n\n", SLOWDOWN_MAX);

    /* Calibration: a list of exact RATIONALS, and it must stay that way.
     *
     * The same trap bench_eval.c documents applies here with extra force: a
     * calibration that packing can touch would move with the very thing this
     * file measures, and a large packing speedup would divide every norm at once
     * and report a catastrophic regression on the day of an improvement. An
     * exact Rational leaf can never live in a machine buffer, so pack_sniff
     * rejects the FIRST element and the packing attempt is O(1) -- immune by
     * construction, not by luck. */
    double calib_us = median_us("Length[Table[k/3, {k, 1200}]]");
    if (calib_us <= 0.0) calib_us = 1e-6;
    printf("calibration  Length[Table[k/3, {k, 1200}]]  = %.1f us\n\n", calib_us);

    printf("%-30s %12s %10s %10s %8s\n",
           "workload", "median(us)", "norm", "baseline", "x base");
    printf("--------------------------------------------------------------------------------\n");

    int slow = 0;
    for (int i = 0; i < N_BENCH; i++) {
        double us = median_us(BENCHES[i].expr);
        double norm = us / calib_us;
        if (BENCHES[i].baseline_norm <= 0.0) {
            printf("%-30s %12.1f %10.4f %10s %8s   (record as baseline)\n",
                   BENCHES[i].label, us, norm, "-", "-");
            continue;
        }
        double x = norm / BENCHES[i].baseline_norm;
        int bad = (x > SLOWDOWN_MAX);
        printf("%-30s %12.1f %10.4f %10.4f %8.2f%s\n",
               BENCHES[i].label, us, norm, BENCHES[i].baseline_norm, x,
               bad ? "  <== MUCH SLOWER" : "");
        if (bad) slow++;
    }

    printf("--------------------------------------------------------------------------------\n");
    if (slow) {
        printf("FAIL: %d workload(s) more than %.1fx slower than baseline\n",
               slow, SLOWDOWN_MAX);
        printf("A packed fast path has probably stopped being reached. Check the\n"
               "AWARE / INT64_OK lists in src/pack.c against the head named above.\n");
        return 1;
    }
    printf("PASS: no workload exceeded %.1fx its baseline cost\n", SLOWDOWN_MAX);
    return 0;
}
