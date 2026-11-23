/*
 * tests/bench_ds_updates.c
 *
 * Update-loop baseline matrix for the data-structure program.
 *
 * Each row runs an interpreter loop that performs n single-element updates on
 * one structure, timed at n and 2n.  The doubling ratio t(2n)/t(n) is
 * machine-independent: ~2 means O(1) amortised per update, ~4 means each update
 * copies the whole structure (O(n) per update, O(n^2) total).
 *
 * Default mode REPORTS only (exit 0) so the known-quadratic rows do not break
 * ctest while the baseline is being recorded.  With DS_BENCH_STRICT=1 every row
 * with a `want_linear` flag must stay under RATIO_MAX, which is how each
 * optimisation phase locks in its win.
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

#define N_TRIALS  3
#define RATIO_MAX 3.0       /* linear ~2.0, accidental quadratic ~4.0 */

typedef struct { const char* label; const char* fmt; int want_linear; int n; } Row;
/* Per-row n: rows whose per-update cost is a memcpy (packed buffers) look linear
 * at n~1e3 and only show their O(n)-per-update quadratic at n~1e4+, so they
 * are measured larger. Slow pointer-copying rows stay small to bound runtime. */

static const Row ROWS[] = {
    /* label                    loop (one %d = n)                                         linear? */
    { "assoc AssociateTo",      "Module[{a = <||>}, Do[AssociateTo[a, i -> i], {i, %d}]; Length[a]]", 1, 25000 },
    { "assoc AppendTo",         "Module[{a = <||>}, Do[AppendTo[a, i -> i], {i, %d}]; Length[a]]", 1, 25000 },
    { "assoc a[[Key[k]]]=v",    "Module[{a = <|1 -> 1|>}, Do[a[[Key[i]]] = i, {i, 2, %d}]; Length[a]]", 1, 25000 },
    { "assoc pure Append",      "Module[{a = <||>}, Do[a = Append[a, i -> i], {i, %d}]; Length[a]]", 0, 1500 },
    { "list pure Append",       "Module[{l = {}}, Do[l = Append[l, i], {i, %d}]; Length[l]]", 0, 6000 },
    { "list AppendTo",          "Module[{l = {}}, Do[AppendTo[l, i], {i, %d}]; Length[l]]", 1, 25000 },
    { "list AppendTo strings",   "Module[{l = {\"a\"}}, Do[AppendTo[l, \"b\"], {i, %d}]; Length[l]]", 1, 25000 },
    { "packed Part-set",        "Module[{p = Range[%d]}, Do[p[[i]] = 0, {i, %d}]; Total[p]]", 1, 25000 },
    { "LinkedList AppendTo",    "Module[{ll = LinkedList[]}, Do[AppendTo[ll, i], {i, %d}]; Length[ll]]", 1, 25000 },
    { "LinkedList pure Append", "Module[{ll = LinkedList[]}, Do[ll = Append[ll, i], {i, %d}]; Length[ll]]", 0, 25000 },
    { "assoc a[[Key[k]]] read", "Module[{a = Association[Table[i -> i, {i, %d}]], s = 0}, Do[s += a[[Key[i]]], {i, %d}]; s]", 1, 25000 },
    { "assoc bulk build",       "Length[Association[Table[i -> i, {i, %d}]]]", 1, 25000 },
};
#define N_ROWS ((int)(sizeof(ROWS) / sizeof(ROWS[0])))

static double now_us(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1.0e6 + (double)ts.tv_nsec / 1.0e3;
}
static int cmp_d(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}
static double median_us(const char* src) {
    double s[N_TRIALS];
    for (int t = -1; t < N_TRIALS; t++) {           /* t==-1 is an untimed warm-up */
        Expr* p = parse_expression(src);
        if (!p) { fprintf(stderr, "parse failed: %s\n", src); exit(2); }
        double t0 = now_us();
        Expr* r = evaluate(p);
        double t1 = now_us();
        expr_free(r); expr_free(p);
        if (t >= 0) s[t] = t1 - t0;
    }
    qsort(s, N_TRIALS, sizeof(double), cmp_d);
    return s[N_TRIALS / 2];
}

int main(void) {
    symtab_init();
    core_init();
    const char* strict_env = getenv("DS_BENCH_STRICT");
    int strict = strict_env && strict_env[0] == '1';

    printf("Data-structure update-loop matrix (per-row n vs 2n, median of %d)  mode=%s\n",
           N_TRIALS, strict ? "STRICT" : "report");
    printf("%-24s %7s %12s %12s %9s %6s  %s\n", "row", "n", "n(us)", "2n(us)", "ns/upd", "ratio", "status");
    printf("--------------------------------------------------------------------------------\n");
    int fail = 0;
    char e[512];
    for (int i = 0; i < N_ROWS; i++) {
        snprintf(e, sizeof(e), ROWS[i].fmt, ROWS[i].n, ROWS[i].n);
        double a = median_us(e);
        snprintf(e, sizeof(e), ROWS[i].fmt, 2 * ROWS[i].n, 2 * ROWS[i].n);
        double b = median_us(e);
        double ratio = a > 0.0 ? b / a : 0.0;
        int quad = ratio > RATIO_MAX;
        const char* st = quad ? (ROWS[i].want_linear ? "REGRESSION" : "quadratic (known)")
                              : "linear";
        printf("%-24s %7d %12.0f %12.0f %9.0f %6.2f  %s\n",
               ROWS[i].label, ROWS[i].n, a, b, b * 1000.0 / (2.0 * ROWS[i].n), ratio, st);
        if (quad && ROWS[i].want_linear) fail++;
    }
    printf("--------------------------------------------------------------------------------\n");
    if (fail) { printf("%s: %d row(s) that must be linear are quadratic\n", strict ? "FAIL" : "NOTE", fail); }
    else printf("PASS: every row marked linear is linear\n");
    return (strict && fail) ? 1 : 0;
}
