/* Tests for interval arithmetic — Interval[{min, max}] and its union form.
 *
 * Covers canonicalization (ordering, merging, disjoint unions), the arithmetic
 * kernels threaded through Plus/Times/Power/Divide/Subtract, elementary-function
 * threading with critical-point range analysis (Sin/Cos/Tan/Exp/Log/Sqrt),
 * Min/Max endpoint extraction, and interval-vs-scalar comparisons.
 *
 * Results are compared on the printed form (interval endpoints are exact where
 * the inputs are exact), matching the Wolfram Language reference outputs.
 */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char* eval_str(const char* input) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    expr_free(e);
    return s;
}

/* Assert that `input` evaluates to exactly `expected` (printed form). */
static bool check(const char* input, const char* expected) {
    char* s = eval_str(input);
    bool ok = (strcmp(s, expected) == 0);
    if (!ok) fprintf(stderr, "  FAIL: %s  =>  %s   (expected %s)\n", input, s, expected);
    free(s);
    return ok;
}

int main(void) {
    symtab_init();
    core_init();
    int failures = 0;

    /* --- canonicalization --- */
    if (!check("Interval[{6, 1}]", "Interval[{1, 6}]")) failures++;
    if (!check("Interval[3]", "Interval[{3, 3}]")) failures++;
    if (!check("Interval[{1, 3}, {2, 4}]", "Interval[{1, 4}]")) failures++;   /* overlap merges */
    if (!check("Interval[{1, 2}, {2, 3}]", "Interval[{1, 3}]")) failures++;   /* touching merges */
    if (!check("Interval[{1, 2}, {5, 6}]", "Interval[{1, 2}, {5, 6}]")) failures++; /* disjoint kept */
    if (!check("Interval[{a, b}]", "Interval[{a, b}]")) failures++;           /* symbolic preserved */

    /* --- arithmetic (exact endpoints stay exact) --- */
    if (!check("Interval[{1, 6}] + Interval[{0, 2}]", "Interval[{1, 8}]")) failures++;
    if (!check("Interval[{1, 2}] + 3", "Interval[{4, 5}]")) failures++;
    if (!check("Interval[{5, 8}] - Interval[{1, 2}]", "Interval[{3, 7}]")) failures++;
    if (!check("-Interval[{1, 3}]", "Interval[{-3, -1}]")) failures++;
    if (!check("2 Interval[{1, 3}]", "Interval[{2, 6}]")) failures++;
    if (!check("Interval[{-1, 2}] Interval[{3, 4}]", "Interval[{-4, 8}]")) failures++;
    if (!check("Interval[{-2, 5}]^2", "Interval[{0, 25}]")) failures++;       /* even power straddles 0 */
    if (!check("Interval[{-3, 2}]^3", "Interval[{-27, 8}]")) failures++;      /* odd power */
    if (!check("1/Interval[{-2, 5}]", "Interval[{-Infinity, -1/2}, {1/5, Infinity}]")) failures++;
    if (!check("1/Interval[{-2, -1}]", "Interval[{-1, -1/2}]")) failures++;
    if (!check("Interval[{2, 3}]^(-1)", "Interval[{1/3, 1/2}]")) failures++;

    /* --- elementary functions --- */
    if (!check("Sin[Interval[{2, 7}]]", "Interval[{-1, Sin[2]}]")) failures++;  /* trough pinned, max symbolic */
    if (!check("Sin[Interval[{2, 10}]]", "Interval[{-1, 1}]")) failures++;      /* crest + trough */
    if (!check("Cos[Interval[{0, 1}]]", "Interval[{Cos[1], 1}]")) failures++;   /* crest at 0 */
    if (!check("Exp[Interval[{0, 1}]]", "Interval[{1, E}]")) failures++;        /* monotone, exact */
    if (!check("Log[Interval[{1, E}]]", "Interval[{0, 1}]")) failures++;
    if (!check("Sqrt[Interval[{4, 9}]]", "Interval[{2, 3}]")) failures++;
    if (!check("Interval[{4, 9}]^(1/2)", "Interval[{2, 3}]")) failures++;
    if (!check("Tan[Interval[{0, 1}]]", "Interval[{0, Tan[1]}]")) failures++;
    if (!check("Abs[Interval[{-2, 1}]]", "Interval[{0, 2}]")) failures++;    /* straddles 0: min 0 */
    if (!check("Abs[Interval[{-5, -2}]]", "Interval[{2, 5}]")) failures++;   /* all negative */
    if (!check("Abs[Interval[{2, 5}]]", "Interval[{2, 5}]")) failures++;     /* all positive */

    /* --- reciprocal trig / hyperbolic (built as 1/base) --- */
    if (!check("Sec[Interval[{0, 1}]]", "Interval[{1, Sec[1]}]")) failures++;
    if (!check("Sec[Interval[{1, 2}]]",
               "Interval[{-Infinity, Sec[2]}, {Sec[1], Infinity}]")) failures++;   /* pole at Pi/2 */
    if (!check("Sech[Interval[{0, 1}]]", "Interval[{Sech[1], 1}]")) failures++;
    if (!check("Csch[Interval[{1, 2}]]", "Interval[{Csch[2], Csch[1]}]")) failures++;
    if (!check("Coth[Interval[{1, 2}]]", "Interval[{Coth[2], Coth[1]}]")) failures++;

    /* --- Sign / Floor / Ceiling (monotone non-decreasing) --- */
    if (!check("Sign[Interval[{-2, 3}]]", "Interval[{-1, 1}]")) failures++;
    if (!check("Sign[Interval[{2, 5}]]", "Interval[{1, 1}]")) failures++;
    if (!check("Floor[Interval[{1.2, 3.8}]]", "Interval[{1, 3}]")) failures++;
    if (!check("Ceiling[Interval[{1.2, 3.8}]]", "Interval[{2, 4}]")) failures++;

    /* --- special functions (monotone on a region) --- */
    if (!check("Erf[Interval[{0, 1}]]", "Interval[{0, Erf[1]}]")) failures++;
    if (!check("Erfc[Interval[{0, 1}]]", "Interval[{Erfc[1], 1}]")) failures++;       /* decreasing */
    if (!check("Gamma[Interval[{2, 3}]]", "Interval[{1, 2}]")) failures++;            /* increasing branch */
    if (!check("Gamma[Interval[{1/2, 1}]]", "Interval[{1, Sqrt[Pi]}]")) failures++;   /* decreasing branch */
    if (!check("Gamma[Interval[{1, 2}]]", "Gamma[Interval[{1, 2}]]")) failures++;     /* straddles min: symbolic */
    if (!check("LogGamma[Interval[{2, 3}]]", "Interval[{0, Log[2]}]")) failures++;
    if (!check("Zeta[Interval[{2, 3}]]", "Interval[{Zeta[3], 1/6 Pi^2}]")) failures++; /* decreasing on (1,inf) */
    if (!check("Zeta[Interval[{1, 2}]]", "Zeta[Interval[{1, 2}]]")) failures++;        /* pole at 1: symbolic */
    if (!check("PolyGamma[0, Interval[{1, 2}]]", "Interval[{-EulerGamma, 1 - EulerGamma}]")) failures++;
    if (!check("PolyGamma[1, Interval[{1, 2}]]", "Interval[{-1 + 1/6 Pi^2, 1/6 Pi^2}]")) failures++; /* decreasing */

    /* --- inverse-reciprocal functions (built as inverse_base(1/x)) --- */
    if (!check("ArcCot[Interval[{1, 2}]]", "Interval[{ArcTan[1/2], 1/4 Pi}]")) failures++;
    if (!check("ArcCot[Interval[{-1, 1}]]",
               "Interval[{-1/2 Pi, -1/4 Pi}, {1/4 Pi, 1/2 Pi}]")) failures++;   /* discontinuity at 0 */
    if (!check("ArcSec[Interval[{2, 3}]]", "Interval[{1/3 Pi, ArcCos[1/3]}]")) failures++;
    if (!check("ArcCsc[Interval[{2, 3}]]", "Interval[{ArcSin[1/3], 1/6 Pi}]")) failures++;
    if (!check("ArcCoth[Interval[{2, 3}]]", "Interval[{ArcTanh[1/3], ArcTanh[1/2]}]")) failures++;
    if (!check("ArcSech[Interval[{1/2, 1}]]", "Interval[{0, ArcCosh[2]}]")) failures++;
    if (!check("ArcCsch[Interval[{1, 2}]]", "Interval[{ArcSinh[1/2], ArcSinh[1]}]")) failures++;
    if (!check("ArcCsch[Interval[{-2, 3}]]",
               "Interval[{-Infinity, -ArcSinh[1/2]}, {ArcSinh[1/3], Infinity}]")) failures++;

    /* --- domain guard: out-of-domain intervals stay symbolic (not complex) --- */
    if (!check("ArcSec[Interval[{1/5, 1/2}]]", "ArcSec[Interval[{1/5, 1/2}]]")) failures++;
    if (!check("ArcSin[Interval[{2, 3}]]", "ArcSin[Interval[{2, 3}]]")) failures++;
    if (!check("Log[Interval[{-2, -1}]]", "Log[Interval[{-2, -1}]]")) failures++;

    /* --- Min / Max --- */
    if (!check("Min[Interval[{2, 7}]]", "2")) failures++;
    if (!check("Max[Interval[{2, 7}]]", "7")) failures++;

    /* --- comparisons (disjoint decide, straddling stays symbolic) --- */
    if (!check("Interval[{5, 8}] > Pi", "True")) failures++;
    if (!check("Interval[{5, 8}] < Pi", "False")) failures++;
    if (!check("Interval[{1, 2}] == 5", "False")) failures++;
    if (!check("Interval[{1, 4}] > Pi", "Interval[{1, 4}] > Pi")) failures++;   /* Pi inside: symbolic */

    /* --- N threads over endpoints --- */
    if (!check("N[Interval[{1, 2}]]", "Interval[{1.0, 2.0}]")) failures++;

    /* --- companion set operations --- */
    if (!check("IntervalUnion[Interval[{1, 3}], Interval[{5, 7}]]", "Interval[{1, 3}, {5, 7}]")) failures++;
    if (!check("IntervalUnion[Interval[{1, 3}], Interval[{2, 5}]]", "Interval[{1, 5}]")) failures++;
    if (!check("IntervalIntersection[Interval[{1, 5}], Interval[{3, 8}]]", "Interval[{3, 5}]")) failures++;
    if (!check("IntervalIntersection[Interval[{1, 2}], Interval[{5, 6}]]", "Interval[]")) failures++;
    if (!check("IntervalMemberQ[Interval[{1, 5}], 3]", "True")) failures++;
    if (!check("IntervalMemberQ[Interval[{1, 5}], 7]", "False")) failures++;
    if (!check("IntervalMemberQ[Interval[{1, 5}], Pi]", "True")) failures++;
    if (!check("IntervalMemberQ[Interval[{1, 5}], Interval[{2, 3}]]", "True")) failures++;

    /* --- interval vs interval comparisons --- */
    if (!check("Interval[{1, 2}] < Interval[{5, 6}]", "True")) failures++;
    if (!check("Interval[{1, 2}] > Interval[{5, 6}]", "False")) failures++;
    if (!check("Interval[{1, 4}] < Interval[{3, 6}]", "Interval[{1, 4}] < Interval[{3, 6}]")) failures++;

    /* --- more elementary functions --- */
    if (!check("Sinh[Interval[{0, 1}]]", "Interval[{0, Sinh[1]}]")) failures++;
    if (!check("Cosh[Interval[{-1, 2}]]", "Interval[{1, Cosh[2]}]")) failures++;      /* even, min 1 at 0 */
    if (!check("Tanh[Interval[{0, 1}]]", "Interval[{0, Tanh[1]}]")) failures++;
    if (!check("ArcTan[Interval[{0, 1}]]", "Interval[{0, 1/4 Pi}]")) failures++;
    if (!check("ArcSin[Interval[{0, 1}]]", "Interval[{0, 1/2 Pi}]")) failures++;
    if (!check("ArcCos[Interval[{0, 1}]]", "Interval[{0, 1/2 Pi}]")) failures++;       /* decreasing */

    /* --- general power --- */
    if (!check("2^Interval[{1, 3}]", "Interval[{2, 8}]")) failures++;                  /* scalar^interval */
    if (!check("(1/2)^Interval[{1, 3}]", "Interval[{1/8, 1/2}]")) failures++;          /* base < 1: decreasing */
    if (!check("Interval[{2, 3}]^Interval[{1, 2}]", "Interval[{2, 9}]")) failures++;   /* interval^interval */

    if (failures == 0) {
        printf("All interval tests passed.\n");
        return 0;
    }
    fprintf(stderr, "%d interval test(s) FAILED.\n", failures);
    return 1;
}
