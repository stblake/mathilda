/* bench_separable_decomposition.c -- performance gate for
 * SeparablePolynomialDecomposition (src/internal/simp/SeparablePolynomialDecomposition.m).
 *
 * Two machine-independent gates (nonzero exit on failure), modelled on
 * bench_assoc.c:
 *
 *   1. Absolute cost, normalised.  The decomposition of the dense rank-1
 *      Expand[(1 + a + b + c)^6 (1 + x + y + z)^6] (7056 terms) is timed
 *      against a calibration measured in the same process: the Expand that
 *      built it.  That ratio is machine-independent.  Fails above 2.5x
 *      its baseline.  This is what the per-(row, column) CoefficientRules
 *      in the first version of the flattening would trip: 88x on the matrix,
 *      3.6x on the whole call at 27225 terms.
 *   2. Scaling.  Random rank-16 and rank-32 sums of separable products in
 *      {a, b, c} x {x, y, z} (1405 and 2928 terms) give the doubling ratio
 *      t(2n)/t(n); measured 2.4-2.5.  Fails above 3.3 (an O(n^2) step makes
 *      it ~4).
 *
 * Baselines: Apple M-series, 2026-11-05.  Each time is the best of 3 runs.
 * Re-record only after an INTENDED performance change.  Nightly tier
 * (bench_*): wall-clock gates are noise on a loaded shared runner.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const double BASELINE_NORM    = 7.2;   /* t(decompose) / t(Expand) */
static const double NORM_TOLERANCE   = 2.5;
static const double MAX_DOUBLING     = 3.3;

/* Evaluate `src` and return it as a double (it must evaluate to a number). */
static double eval_number(const char* src) {
    Expr* e = parse_expression(src);
    if (!e) { fprintf(stderr, "bench: does not parse: %s\n", src); exit(2); }
    Expr* r = evaluate(e);
    expr_free(e);
    double v = -1.0;
    if (r && r->type == EXPR_REAL) v = r->data.real;
    else if (r && r->type == EXPR_INTEGER) v = (double)r->data.integer;
    if (v < 0.0) {
        char* s = r ? expr_to_string(r) : NULL;
        fprintf(stderr, "bench: not a number: %s -> %s\n", src, s ? s : "(null)");
        free(s);
        exit(2);
    }
    expr_free(r);
    return v;
}

static void run(const char* src) {
    Expr* e = parse_expression(src);
    if (!e) { fprintf(stderr, "bench: does not parse: %s\n", src); exit(2); }
    Expr* r = evaluate(e);
    expr_free(e);
    expr_free(r);
}

int main(void) {
    symtab_init();
    core_init();
    int fails = 0;

    run("spdBest[f_] := Min[Table[First@AbsoluteTiming[f[]], {3}]]; Null");
    run("spdPr = (1 + a + b + c)^6 (1 + x + y + z)^6; spdP = Expand[spdPr]; Null");
    run("SeedRandom[3]; spdRP[vars_, deg_, nt_] := Total@Table[RandomInteger[{-9, 9}] "
        "Times @@ (vars^RandomInteger[{0, deg}, Length[vars]]), {k, nt}]; "
        "spdQ16 = Expand[Total@Table[spdRP[{a, b, c}, 8, 10] spdRP[{x, y, z}, 8, 10], {16}]]; "
        "spdQ32 = Expand[Total@Table[spdRP[{a, b, c}, 8, 10] spdRP[{x, y, z}, 8, 10], {32}]]; Null");

    /* correctness first: a fast wrong answer is not a pass */
    double ok = eval_number(
        "Boole[Expand[SeparablePolynomialDecomposition[spdP, {a, b, c}, {x, y, z}] - spdP] === 0 && "
        "Expand[SeparablePolynomialDecomposition[spdQ32, {a, b, c}, {x, y, z}] - spdQ32] === 0]");
    if (ok != 1.0) { fprintf(stderr, "FAIL: decomposition is not exact\n"); return 1; }

    double cal = eval_number("spdBest[Expand[spdPr] &]");
    double t   = eval_number("spdBest[SeparablePolynomialDecomposition[spdP, {a, b, c}, {x, y, z}] &]");
    double norm = t / cal;
    printf("dense rank 1, 7056 terms: %.4f s, calibration %.4f s, norm %.2f "
           "(baseline %.2f, limit %.2f)\n",
           t, cal, norm, BASELINE_NORM, BASELINE_NORM * NORM_TOLERANCE);
    if (norm > BASELINE_NORM * NORM_TOLERANCE) {
        fprintf(stderr, "FAIL: normalised cost %.2f exceeds %.2f\n",
                norm, BASELINE_NORM * NORM_TOLERANCE);
        fails++;
    }

    double t16 = eval_number("spdBest[SeparablePolynomialDecomposition[spdQ16, {a, b, c}, {x, y, z}] &]");
    double t32 = eval_number("spdBest[SeparablePolynomialDecomposition[spdQ32, {a, b, c}, {x, y, z}] &]");
    double n16 = eval_number("Length[spdQ16]");
    double n32 = eval_number("Length[spdQ32]");
    double ratio = t32 / t16;
    printf("rank 16 -> 32 (%.0f -> %.0f terms): %.4f s -> %.4f s, doubling ratio %.2f (limit %.2f)\n",
           n16, n32, t16, t32, ratio, MAX_DOUBLING);
    if (ratio > MAX_DOUBLING) {
        fprintf(stderr, "FAIL: doubling ratio %.2f exceeds %.2f\n", ratio, MAX_DOUBLING);
        fails++;
    }

    if (fails) return 1;
    printf("SeparablePolynomialDecomposition performance gates passed.\n");
    return 0;
}
