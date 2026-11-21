/*
 * bench_calibrate.h -- the machine-speed unit the bench_* gates divide by.
 *
 * Each gate reports norm = workload_us / calibration_us and fails when a norm
 * exceeds SLOWDOWN_MAX x its recorded baseline.  The calibration exists only to
 * cancel the speed of the MACHINE (and its momentary load); it must not move
 * when MATHILDA changes.  Every calibration spelled in the Wolfram language
 * eventually did:
 *
 *   Total[Range[40000]]           automatic packed arrays made it ~20x faster;
 *   Table[k/3, {k, n}]            v0.6319 (terminal Rationals) made it 1.65x
 *                                 faster (baselines rescaled 2026-11-07), and
 *                                 the v0.6653-57 scoping work 2.9x more --
 *
 * and each speedup of the unit divides every norm at once, so the day of an
 * improvement read as a regression of every row.  This unit is plain C, so no
 * change to the evaluator can touch it.  It mixes the two costs the gates
 * measure: small-node allocation with dependent pointer chasing (what an Expr
 * tree costs) and a large buffer allocated, written, read and freed (what a
 * packed-array result costs, page faults included -- the part that differs most
 * between the macOS and glibc allocators).
 *
 * Baselines are still per PLATFORM (bench_*.c carry a Darwin and a Linux
 * column): one C loop cannot capture how differently two OSes price fresh
 * memory, and the Linux runner pays ~7x macOS for a 10^5-element RotateLeft.
 */
#ifndef MATHILDA_BENCH_CALIBRATE_H
#define MATHILDA_BENCH_CALIBRATE_H

#include <stdlib.h>
#include <stdint.h>
#include <time.h>

typedef struct BenchCalNode { struct BenchCalNode* next; uint64_t v; double pad[4]; } BenchCalNode;

static volatile uint64_t bench_cal_sink;

static double bench_cal_now_us(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1.0e6 + (double)ts.tv_nsec / 1.0e3;
}

/* One pass: ~4000 small nodes linked in a scrambled order and walked, then a
 * 10^5-double buffer allocated, filled, summed and freed.  A few hundred us. */
static void bench_cal_pass(void) {
    enum { NODES = 4000, BUF = 100000 };
    BenchCalNode* nodes[NODES];
    for (int i = 0; i < NODES; i++) {
        nodes[i] = (BenchCalNode*)malloc(sizeof(BenchCalNode));
        nodes[i]->v = (uint64_t)i * UINT64_C(2654435761);
    }
    /* link in a fixed pseudo-random permutation (stride coprime to NODES) */
    for (int i = 0; i < NODES; i++)
        nodes[(i * 1597) % NODES]->next = nodes[((i + 1) * 1597) % NODES];
    uint64_t acc = 0;
    BenchCalNode* p = nodes[0];
    for (int i = 0; i < 8 * NODES; i++) { acc = acc * 31 + p->v; p = p->next; }
    for (int i = 0; i < NODES; i++) free(nodes[i]);

    double* buf = (double*)malloc(BUF * sizeof(double));
    for (int i = 0; i < BUF; i++) buf[i] = (double)(i ^ (int)(acc & 1023));
    double s = 0.0;
    for (int i = 0; i < BUF; i++) s += buf[i];
    free(buf);
    bench_cal_sink = acc + (uint64_t)s;
}

static int bench_cal_cmp(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

/* Median over `trials` of `reps` passes, in us.  One untimed warm-up pass. */
static double bench_calibrate_us(int trials, int reps) {
    double s[64];
    if (trials > 64) trials = 64;
    if (trials < 1) trials = 1;
    bench_cal_pass();
    for (int t = 0; t < trials; t++) {
        double t0 = bench_cal_now_us();
        for (int r = 0; r < reps; r++) bench_cal_pass();
        s[t] = bench_cal_now_us() - t0;
    }
    qsort(s, (size_t)trials, sizeof(double), bench_cal_cmp);
    double m = s[trials / 2];
    return m > 0.0 ? m : 1e-6;
}

/* Pick the recorded baseline for the platform this binary was built for.  A
 * zero means "not recorded on this platform yet": the row prints its norm so it
 * can be recorded, and does not gate. */
#if defined(__linux__)
#  define BENCH_PLATFORM "Linux"
#  define BENCH_BASELINE(darwin, linux_) (linux_)
#elif defined(__APPLE__)
#  define BENCH_PLATFORM "Darwin"
#  define BENCH_BASELINE(darwin, linux_) (darwin)
#else
#  define BENCH_PLATFORM "other (ungated)"
#  define BENCH_BASELINE(darwin, linux_) (0.0)
#endif

#endif /* MATHILDA_BENCH_CALIBRATE_H */
