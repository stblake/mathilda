/* bench_intrischnorman.c — wall-clock benchmark for Integrate`RischNorman.
 *
 * Runs a curated corpus of integrands (drawn from test_intrischnorman.c
 * phases 4–7) through Integrate`RischNorman[f, x] and prints per-case
 * wall-clock time, sorted by slowest at the end.
 *
 * Usage:
 *   ./bench_intrischnorman               # one warm-up + 1 measured run
 *   ./bench_intrischnorman --reps N      # N measured runs per case
 *   PMINT_PROFILE=1 ./bench_intrischnorman  # show per-phase times
 *                                          # (when phase timers are built in)
 *
 * Exit status is 0 even if a case returns the call unevaluated — this
 * is a perf tool, not a correctness check.  See test_intrischnorman.c
 * for the correctness predicate.
 */

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* --------------------------------------------------------------------- */
/* Monotonic wall-clock helper.                                          */
/* --------------------------------------------------------------------- */

static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* --------------------------------------------------------------------- */
/* Corpus.  Curated across phases 4–7 of test_intrischnorman.c.          */
/* --------------------------------------------------------------------- */

typedef struct {
    const char* label;        /* short tag, e.g. "p4-exp_x" */
    const char* integrand;    /* picocas-syntax integrand   */
} BenchCase;

static const BenchCase corpus[] = {
    /* Phase 4 — Q-rational base cases. */
    { "p4-exp_x",            "Exp[x]" },
    { "p4-exp_2x",           "Exp[2 x]" },
    { "p4-x_exp_x",          "x Exp[x]" },
    { "p4-x2_exp_x",         "x^2 Exp[x]" },
    { "p4-log_x",            "Log[x]" },
    { "p4-x_log_x",          "x Log[x]" },
    { "p4-sin_exp",          "Sin[x] Exp[x]" },
    { "p4-cos_exp",          "Cos[x] Exp[x]" },

    /* Phase 5 — log-bearing + Weierstrass trig. */
    { "p5-inv_x",            "1/x" },
    { "p5-sin_x",            "Sin[x]" },
    { "p5-cos_x",            "Cos[x]" },
    { "p5-tan_x",            "Tan[x]" },
    { "p5-cot_x",            "Cot[x]" },
    { "p5-inv_1pexp",        "1/(1 + Exp[x])" },
    { "p5-exp_over_1pexp",   "Exp[x]/(1 + Exp[x])" },
    { "p5-inv_x_logx",       "1/(x Log[x])" },
    { "p5-log_sq",           "Log[x]^2" },
    { "p5-1plogxlogx",       "(1 + Log[x])/(x Log[x])" },

    /* Phase 6 — assorted, including trig powers. */
    { "p6-exp_ax",           "Exp[a x]" },
    { "p6-x3_exp_x",         "x^3 Exp[x]" },
    { "p6-log_cube",         "Log[x]^3" },
    { "p6-sinh_x",           "Sinh[x]" },
    { "p6-cosh_x",           "Cosh[x]" },
    { "p6-tanh_x",           "Tanh[x]" },
    { "p6-sin_sq",           "Sin[x]^2" },
    { "p6-cos_sq",           "Cos[x]^2" },
    { "p6-tan_sq",           "Tan[x]^2" },
    { "p6-inv_x2p1",         "1/(1+x^2)" },

    /* Phase 7 — pmint paper / textbook hard cases. */
    { "p7-pmint_paper_1",
      "(x^7 - 24 x^4 - 4 x^2 + 8 x - 8)/(x^8 + 6 x^6 + 12 x^4 + 8 x^2)" },
    { "p7-pmint_paper_2",
      "(-4 x^2 - 4 x^3 - x^4)/((-1 + x^2) (1 + x + x^2)^2)" },
    { "p7-x3_over_xp1",      "x^3/(1 + x)" },
    { "p7-rat_log_complex",
      "(x^4 - 3 x^2 + 6)/(x^6 - 5 x^4 + 5 x^2 + 4)" },
    { "p7-log_neg2",         "(-1 + Log[x])/Log[x]^2" },
    { "p7-exp_over_log_shifted",
      "(E^x (-1 - Log[-1 + x] + x Log[-1 + x]))/((-1 + x) Log[-1 + x]^2)" },
    { "p7-log_diff_sq",      "(1 - Log[x])/(x^2 - Log[x]^2)" },
    { "p7-sin_over_x_plus_log_cos",
      "Sin[x]/x + Log[x] Cos[x]" },
    { "p7-exp_inv_x",        "((x + 1)/x^4) Exp[1/x]" },
    { "p7-exp_x_plus_inv_log",
      "(1 - 1/(x Log[x]^2)) Exp[1/Log[x] + x]" },
    { "p7-exp_rat",
      "((-1 - x - x^2 + x^3)/(1 - 2 x^2 + x^4)) Exp[x]" },
    { "p7-exp_inv_log",
      "(Exp[1/Log[x]] (Log[x]^2 - 1))/Log[x]^2" },
    { "p7-exp_quad_denom",
      "(E^x (-1 + 1289 x + 278 x^2 + 55 x^3 + 56 x^4))/(-392 + x - 56 x^2)^2" },
    { "p7-exp_x2_diff",
      "((4 x^2 + 4 x - 1) (Exp[x^2] + 1) (Exp[x^2] - 1))/(x + 1)^2" },
    { "p7-x3_exp_x2",        "x^3 Exp[x^2]" },
    { "p7-exp_inv_x_x3",     "Exp[1/x]/x^3" },
    { "p7-x2mx_exp3x",       "(x^2 - x) Exp[3 x]" },
    { "p7-trig_x2_poly",
      "(x^3 + 2 x) Cos[x^2] + x Sin[x^2]" },
    { "p7-x3_sin_x2p1",      "x^3 Sin[x^2 + 1]" },
    { "p7-arctan_over_x2",   "ArcTan[x]/x^2" },
    { "p7-tan_div",          "(x - Tan[x])/Tan[x]^2 + Tan[x]" },
    { "p7-cot_half",         "(1 - Cos[x] + Sin[x])/(-1 + Cos[x])^2" },
    { "p7-exp_sin_2x",       "Exp[x] Sin[2 x]" },
    { "p7-sin_sin_2x",       "Sin[x] Sin[2 x]" },
};

static const size_t corpus_n = sizeof(corpus) / sizeof(corpus[0]);

/* --------------------------------------------------------------------- */
/* One case.                                                              */
/* --------------------------------------------------------------------- */

typedef struct {
    const char* label;
    double seconds;   /* fastest of N runs (less noisy than mean)        */
    int succeeded;    /* 1 if returned != unevaluated head call          */
} BenchResult;

static int is_unevaluated_call(const Expr* res) {
    /* Heuristic: a leaf head returning `Integrate`RischNorman[...]` is
     * unevaluated.  We check the printed string. */
    if (!res) return 1;
    char* s = expr_to_string((Expr*)res);
    int u = (s && strstr(s, "Integrate`RischNorman[") == s);
    free(s);
    return u;
}

static double time_one(const char* integrand, int* succeeded) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
             "Integrate`RischNorman[%s, x]", integrand);

    Expr* parsed = parse_expression(buf);
    if (!parsed) {
        if (succeeded) *succeeded = 0;
        return -1.0;
    }

    double t0 = now_seconds();
    Expr* res = evaluate(parsed);
    double t1 = now_seconds();

    if (succeeded) *succeeded = !is_unevaluated_call(res);

    expr_free(parsed);
    expr_free(res);
    return t1 - t0;
}

/* --------------------------------------------------------------------- */
/* Driver.                                                                */
/* --------------------------------------------------------------------- */

static int cmp_result_desc(const void* a, const void* b) {
    const BenchResult* ra = (const BenchResult*)a;
    const BenchResult* rb = (const BenchResult*)b;
    if (ra->seconds < rb->seconds) return 1;
    if (ra->seconds > rb->seconds) return -1;
    return 0;
}

int main(int argc, char** argv) {
    int reps = 1;
    const char* only = NULL;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--reps") == 0 && i + 1 < argc) {
            reps = atoi(argv[++i]);
            if (reps < 1) reps = 1;
        } else if (strcmp(argv[i], "--only") == 0 && i + 1 < argc) {
            only = argv[++i];
        }
    }

    symtab_init();
    core_init();

    /* Warm up: parse / evaluate something tiny so static caches inside
     * the evaluator / symbol table are populated before the first
     * measured call. */
    {
        Expr* e = parse_expression("Integrate`RischNorman[1, x]");
        Expr* r = evaluate(e);
        expr_free(e); expr_free(r);
    }

    BenchResult* results = (BenchResult*)calloc(corpus_n, sizeof(BenchResult));

    fprintf(stdout,
            "%-34s  %-8s  %-7s  %s\n",
            "label", "ok?", "reps", "best (s)");
    fprintf(stdout,
            "----------------------------------  --------  -------  --------\n");
    fflush(stdout);

    double total = 0.0;
    for (size_t i = 0; i < corpus_n; i++) {
        if (only && strcmp(only, corpus[i].label) != 0) {
            results[i].label = corpus[i].label;
            results[i].seconds = 0.0;
            continue;
        }
        double best = 1e18;
        int ok = 0;
        for (int r = 0; r < reps; r++) {
            int this_ok = 0;
            double t = time_one(corpus[i].integrand, &this_ok);
            if (t >= 0.0 && t < best) best = t;
            ok = ok || this_ok;
        }
        results[i].label = corpus[i].label;
        results[i].seconds = best;
        results[i].succeeded = ok;
        total += best;

        fprintf(stdout, "%-34s  %-8s  %-7d  %8.4f\n",
                corpus[i].label,
                ok ? "ok" : "FAIL",
                reps,
                best);
        fflush(stdout);
    }

    fprintf(stdout,
            "\nTotal across %zu cases (best-of-%d): %.4f s\n",
            corpus_n, reps, total);

    /* Top-N slowest. */
    qsort(results, corpus_n, sizeof(BenchResult), cmp_result_desc);
    fprintf(stdout, "\nSlowest 10:\n");
    size_t top = corpus_n < 10 ? corpus_n : 10;
    for (size_t i = 0; i < top; i++) {
        fprintf(stdout, "  %-34s  %8.4f s  %s\n",
                results[i].label,
                results[i].seconds,
                results[i].succeeded ? "" : "(unevaluated)");
    }

    free(results);
    return 0;
}
