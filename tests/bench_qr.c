/*
 * tests/bench_qr.c
 *
 * Phase-5 microbench for QRDecomposition.  Prints per-(shape x precision
 * x pivoting) wall-clock timings.  Not wired into the CI / ctest as a
 * pass/fail target -- the plan calls for ratio checks against a checked-in
 * `tests/data/qr_mathematica.json` baseline once that file lands.  Until
 * then this binary is a spot-check tool: run it after a perf-relevant
 * change to confirm the LAPACK / MPFR fast paths are still winning.
 *
 * The plan's success bar from S7 of tasks/qr_perf_plan.md is reproduced
 * inline below.  Eyeball against the printed timings.
 *
 *   | Case                | Mma t   | target                |
 *   |---------------------|---------|-----------------------|
 *   | symbolic int n=4    |  696 us |    <= 1 ms            |
 *   | machine real n=20   |    8 us |    <= 100 us          |
 *   | machine real n=100  |    ~ms  |    <= 5 ms            |
 *   | machine cplx n=20   |   18 us |    <= 200 us          |
 *   | MPFR-30 n=20        |  2.7 ms |    <= 10 ms           |
 *   | MPFR-100 n=20       |  2.3 ms |    <= 30 ms           |
 *
 * Each case runs N_TRIALS evaluations of QRDecomposition[m] on the same
 * pre-parsed input and prints the median (so a cold-cache outlier doesn't
 * distort the report).
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"

#define N_TRIALS 7

/* Get wall-clock time in nanoseconds. */
static double now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1.0e9 + (double)ts.tv_nsec;
}

static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

/* Pre-evaluate `m_src` once, then time N_TRIALS QR factorisations on
 * the cached matrix.  This mirrors the convention used in published
 * Mathematica benchmarks (m is built once outside the timing block;
 * QRDecomposition[m] is what gets timed).  Without this, the report
 * mostly measures rational-to-MPFR conversion of the literal matrix. */
static double bench_qr(const char* label, const char* m_src,
                       bool pivoting, double target_us) {
    /* Pre-evaluate m once.  Symbol-store it so the QR call below can
     * reference it via a short token. */
    {
        char setup[64 + 1];
        snprintf(setup, sizeof(setup), "benchMlocal =");
        size_t setlen = strlen(m_src) + 64;
        char* full = (char*)malloc(setlen);
        snprintf(full, setlen, "benchMlocal = %s", m_src);
        Expr* parsed = parse_expression(full);
        Expr* res = evaluate(parsed);
        expr_free(res);
        expr_free(parsed);
        free(full);
    }

    const char* qr_expr = pivoting
        ? "QRDecomposition[benchMlocal, Pivoting -> True]"
        : "QRDecomposition[benchMlocal]";

    /* Warm up once so the symbol-table / printer / etc caches are
     * primed -- otherwise the first trial reads ~10x slower than the
     * rest and skews the median on small cases. */
    {
        Expr* parsed = parse_expression(qr_expr);
        Expr* res = evaluate(parsed);
        expr_free(res);
        expr_free(parsed);
    }

    double samples[N_TRIALS];
    for (int t = 0; t < N_TRIALS; t++) {
        Expr* parsed = parse_expression(qr_expr);
        double t0 = now_ns();
        Expr* res = evaluate(parsed);
        double t1 = now_ns();
        expr_free(res);
        expr_free(parsed);
        samples[t] = (t1 - t0) / 1000.0;   /* nanoseconds -> microseconds */
    }
    qsort(samples, N_TRIALS, sizeof(double), cmp_double);
    double median = samples[N_TRIALS / 2];
    const char* tag = median <= target_us ? "  OK " : " SLOW";
    printf("%s  %-40s  median = %10.2f us   (target %.0f us)\n",
           tag, label, median, target_us);
    return median;
}

/* Build a "dense pseudo-random" n x n matrix literal at the given
 * precision wrapper.
 *   prec  = 0   -> exact integer matrix (no fractions, so symbolic MGS
 *                  doesn't pay rational-arithmetic blow-up tax).
 *   prec  = 53  -> machine (N[..., 53]); entries are 5/100-style
 *                  fractions wrapped in N[] so they become Reals.
 *   prec  > 53  -> MPFR; entries are exact fractions wrapped in
 *                  N[..., digits] so they become EXPR_MPFR.
 *   is_complex makes every other entry have an I component. */
static char* build_matrix(int n, long prec, bool is_complex) {
    /* Allocate generously: each entry is ~24 chars worst case. */
    size_t cap = (size_t)n * (size_t)n * 32 + 256;
    char* buf = (char*)malloc(cap);
    char* p = buf;
    char* end = buf + cap;
    if (prec >= 53) {
        p += snprintf(p, (size_t)(end - p), "N[");
    }
    p += snprintf(p, (size_t)(end - p), "{");
    for (int i = 0; i < n; i++) {
        if (i) p += snprintf(p, (size_t)(end - p), ", ");
        p += snprintf(p, (size_t)(end - p), "{");
        for (int j = 0; j < n; j++) {
            if (j) p += snprintf(p, (size_t)(end - p), ", ");
            /* Diagonally-dominant pattern so the matrix stays well
             * conditioned.  Mirrors the 10x10 / 50x50 fixtures used by
             * the machine-kernel tests.  For prec == 0 (symbolic
             * exact) we drop the /100 to keep entries integer -- the
             * MGS pipeline pays an order-of-magnitude rational-
             * arithmetic tax when entries carry a common factor of
             * 1/100 across every cell. */
            int v_off = (i * 7 + j * 3) % 11;
            if (i == j) {
                p += snprintf(p, (size_t)(end - p), "%d", 5 * n);
            } else if (prec == 0) {
                p += snprintf(p, (size_t)(end - p), "%d", v_off);
            } else {
                p += snprintf(p, (size_t)(end - p), "%d/100", v_off);
            }
            if (is_complex && ((i + j) % 2)) {
                p += snprintf(p, (size_t)(end - p), " + %d/100 I",
                              (v_off + 3) % 11);
            }
        }
        p += snprintf(p, (size_t)(end - p), "}");
    }
    p += snprintf(p, (size_t)(end - p), "}");
    if (prec > 53) {
        int digits = (int)(prec * 0.30103);
        p += snprintf(p, (size_t)(end - p), ", %d]", digits);
    } else if (prec == 53) {
        p += snprintf(p, (size_t)(end - p), "]");
    }
    return buf;
}

int main(void) {
    symtab_init();
    core_init();

    printf("QRDecomposition microbench (median of %d trials, matrix\n"
           "pre-evaluated outside the timing block).\n\n",
           N_TRIALS);
    printf("kernel/case                                median = ...    (target)\n");
    printf("--------------------------------------------------------------------\n");

    /* Symbolic n=4 (exact int).  Mma: 696 us; target: <= 1000 us. */
    {
        char* m = build_matrix(4, 0, false);
        bench_qr("symbolic int n=4", m, false, 1000.0);
        free(m);
    }

    /* Machine real n=20.  Mma: 8 us; target: <= 100 us. */
    {
        char* m = build_matrix(20, 53, false);
        bench_qr("machine real n=20", m, false, 100.0);
        free(m);
    }

    /* Machine real n=100.  Mma: ~ms; target: <= 5 ms = 5000 us. */
    {
        char* m = build_matrix(100, 53, false);
        bench_qr("machine real n=100", m, false, 5000.0);
        free(m);
    }

    /* Machine complex n=20.  Mma: 18 us; target: <= 200 us. */
    {
        char* m = build_matrix(20, 53, true);
        bench_qr("machine complex n=20", m, false, 200.0);
        free(m);
    }

    /* Machine real n=20 with pivoting. */
    {
        char* m = build_matrix(20, 53, false);
        bench_qr("machine real n=20 + pivot", m, true, 200.0);
        free(m);
    }

    /* MPFR-30 n=20.  Mma: 2.7 ms; target: <= 10 ms = 10000 us. */
    {
        char* m = build_matrix(20, 100, false);   /* 100 bits ~= 30 digits */
        bench_qr("MPFR-30 real n=20", m, false, 10000.0);
        free(m);
    }

    /* MPFR-100 n=20.  Mma: 2.3 ms; target: <= 30 ms = 30000 us. */
    {
        char* m = build_matrix(20, 333, false);   /* 333 bits ~= 100 digits */
        bench_qr("MPFR-100 real n=20", m, false, 30000.0);
        free(m);
    }

    /* MPFR-30 n=10 + pivot. */
    {
        char* m = build_matrix(10, 100, false);
        bench_qr("MPFR-30 real n=10 + pivot", m, true, 5000.0);
        free(m);
    }

    printf("\nBench complete.\n");
    return 0;
}
