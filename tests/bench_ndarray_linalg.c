/* Scaling / stress gate for the NDArray linear-algebra fast paths.
 *
 * Machine-independent doubling-ratio test (same idea as bench_assoc.c): time
 * each op at size n and 2n and check t(2n)/t(n). The dense factorisation fast
 * paths are O(n^3), so the ratio is ~8; we fail only if it blows past RATIO_MAX,
 * which catches the catastrophic failure mode this whole feature exists to
 * prevent — Det silently falling back to the O(n!) symbolic Laplace path.
 *
 * Also a plain stress pass: build and solve/invert/factor large NDArrays and
 * confirm they complete and return a numeric result (no crash, no hang).
 *
 * Exits nonzero on any gate failure so it can wire into CI.
 */
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

extern struct Expr* parse_expression(const char*);
extern struct Expr* evaluate(struct Expr*);
extern void expr_free(struct Expr*);

#define N_TRIALS   5
#define RATIO_MAX  16.0     /* O(n^3) ~ 8; O(n^4) ~ 16; O(n!) is astronomically larger */

static double now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1.0e9 + (double)ts.tv_nsec;
}

static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

/* A diagonally-dominant n x n real NDArray literal (nonsingular, well-conditioned)
 * built with Table so we don't emit a giant literal string. */
static void matrix_expr(char* buf, size_t buflen, const char* op, int n) {
    /* Diagonal scales with n so the matrix stays strictly diagonally dominant
     * (hence nonsingular and well-conditioned) at every size. */
    snprintf(buf, buflen,
        "%s[NDArray[Table[N[If[i == j, 100*%d + i, Mod[i*7 + j*3, 11]]], {i, %d}, {j, %d}]]]",
        op, n, n, n);
}

/* median wall time (us) to evaluate `expr` once. */
static double median_us(const char* expr) {
    { struct Expr* p = parse_expression(expr); struct Expr* r = evaluate(p); expr_free(r); expr_free(p); }
    double s[N_TRIALS];
    for (int t = 0; t < N_TRIALS; t++) {
        struct Expr* p = parse_expression(expr);
        double t0 = now_ns();
        struct Expr* r = evaluate(p);
        double t1 = now_ns();
        expr_free(r); expr_free(p);
        s[t] = (t1 - t0) / 1000.0;
    }
    qsort(s, N_TRIALS, sizeof(double), cmp_double);
    return s[N_TRIALS / 2];
}

static int check_scaling(const char* op, int n) {
    char e1[512], e2[512];
    matrix_expr(e1, sizeof(e1), op, n);
    matrix_expr(e2, sizeof(e2), op, 2 * n);
    double t1 = median_us(e1);
    double t2 = median_us(e2);
    double ratio = (t1 > 0.0) ? (t2 / t1) : 0.0;
    int bad = (ratio > RATIO_MAX);
    printf("  %-14s n=%3d: %9.1f us   2n=%3d: %9.1f us   ratio=%5.2f%s\n",
           op, n, t1, 2 * n, t2, ratio, bad ? "   <== REGRESSION" : "");
    return bad ? 1 : 0;
}

/* Stress: a big op must complete and return a numeric (non-symbolic) result. */
static int stress(const char* op, int n, const char* wrap) {
    char inner[512], full[640];
    matrix_expr(inner, sizeof(inner), op, n);
    snprintf(full, sizeof(full), "%s", inner);
    if (wrap) snprintf(full, sizeof(full), wrap, inner);
    struct Expr* p = parse_expression(full);
    struct Expr* r = evaluate(p);
    /* Result must be a concrete numeric object (NDArray / Real / Integer /
     * List), i.e. the call did not stay unevaluated. `Blank[NDArray]` does not
     * match the NDArray leaf node, so gate on Head instead. */
    char chk[768];
    snprintf(chk, sizeof(chk),
             "MemberQ[{NDArray, Real, Integer, List}, Head[%s]]", full);
    struct Expr* pc = parse_expression(chk);
    struct Expr* rc = evaluate(pc);
    char* s = expr_to_string(rc);
    int ok = (strcmp(s, "True") == 0);
    printf("  stress %-14s n=%d: %s\n", op, n, ok ? "ok" : "FAILED");
    free(s); expr_free(rc); expr_free(pc); expr_free(r); expr_free(p);
    return ok ? 0 : 1;
}

int main(void) {
    symtab_init();
    core_init();

    int fails = 0;
    printf("NDArray linalg scaling gate (RATIO_MAX=%.1f):\n", RATIO_MAX);
    fails += check_scaling("Det", 16);
    fails += check_scaling("Inverse", 16);   /* same LU kernel as LinearSolve */
    fails += check_scaling("MatrixRank", 16);

    printf("NDArray linalg stress (large sizes complete & stay numeric):\n");
    fails += stress("Det", 96, NULL);
    fails += stress("Inverse", 96, NULL);
    fails += stress("MatrixRank", 96, NULL);

    if (fails) {
        printf("BENCH FAILED: %d gate(s) tripped.\n", fails);
        return 1;
    }
    printf("All NDArray linalg scaling/stress gates passed.\n");
    return 0;
}
