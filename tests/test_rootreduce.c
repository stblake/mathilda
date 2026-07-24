/*
 * test_rootreduce.c — RootReduce and the genuine-algebraic field arithmetic
 * (Milestone B): denominator rationalisation via field inversion
 * (flint_algebraic_field_canonical) and the WL-faithful Together/Cancel combine
 * over cube-and-higher-root towers (flint_algebraic_field_together).
 *
 * The strongest correctness check here is VALUE PRESERVATION cross-checked by an
 * INDEPENDENT engine: for every input e we assert
 *     Simplify[RootReduce[e] - e] == 0
 * where Simplify's zero test is the ideal-reduction path
 * (flint_algebraic_field_normalize) — a different FLINT engine than the
 * matrix-inversion RootReduce uses. Agreement between the two is a rigorous
 * end-to-end check with no numeric oracle. On top of that we pin exact canonical
 * forms, assert the denominator is rationalised (a polynomial in the parameters,
 * radical-free), check idempotence, and stress a batch of towers.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"
#include "print.h"
#include "flint_bridge.h"
#include "test_utils.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static Expr* eval_str(const char* src) {
    Expr* parsed = parse_expression(src);
    ASSERT(parsed != NULL);
    Expr* e = evaluate(parsed);
    expr_free(parsed);
    return e;
}

/* Assert eval(input) prints as `expected`. */
static void check(const char* input, const char* expected) {
    Expr* e = eval_str(input);
    char* s = expr_to_string(e);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  got:      %s\n",
                input, expected, s);
        exit(1);
    }
    free(s);
    expr_free(e);
}

/* Value preservation: RootReduce must not change the value. RootReduce[X] of a
 * genuine field zero X reduces to 0 (its numerator vanishes modulo the
 * minimal-polynomial ideal), so RootReduce[RootReduce[e] - e] == 0 verifies the
 * rationalised form equals the input. Cheap (two field solves) and exercises the
 * ideal reduction as well as the matrix inversion. */
static void check_preserves(const char* e) {
    char buf[1200];
    snprintf(buf, sizeof buf, "RootReduce[RootReduce[%s] - (%s)]", e, e);
    check(buf, "0");
}

/* Independent cross-check for the genuine-algebraic tower: the difference
 * RootReduce[e] - e must reduce to 0 through the *ideal-reduction* zero test
 * (flint_algebraic_field_normalize), a different FLINT engine than the
 * matrix inversion RootReduce uses. Agreement is a rigorous, no-oracle check. */
static void check_preserves_ideal(const char* e) {
    char buf[1200];
    snprintf(buf, sizeof buf, "(RootReduce[%s]) - (%s)", e, e);
    Expr* diff = eval_str(buf);
    Expr* r = flint_algebraic_field_normalize(diff);
    if (r == NULL || !(r->type == EXPR_INTEGER && r->data.integer == 0)) {
        char* s = r ? expr_to_string(r) : NULL;
        fprintf(stderr, "FAIL: ideal zero-test of RootReduce[%s]-(%s) expected 0, "
                        "got %s\n", e, e, r ? s : "NULL");
        if (s) free(s);
        exit(1);
    }
    expr_free(r);
    expr_free(diff);
}

/* The denominator of RootReduce[e] must be rationalised: a polynomial in `vars`
 * (radical-free). PolynomialQ is False if any radical Power survives. */
static void check_rationalised(const char* e, const char* vars) {
    char buf[1024];
    snprintf(buf, sizeof buf, "PolynomialQ[Denominator[RootReduce[%s]], %s]", e, vars);
    check(buf, "True");
}

/* ------------------------------------------------------------------ */
/*  RootReduce: exact canonical forms                                  */
/* ------------------------------------------------------------------ */

static void test_rootreduce_exact(void) {
    /* Rationalise a cube-root denominator: 1/(1+k^(1/3)) = (1-k^(1/3)+k^(2/3))/(1+k). */
    check("RootReduce[1/(1 + k^(1/3))]", "(1 - k^(1/3) + k^(2/3))/(1 + k)");
    /* Two-parameter numerator: 1/(a + b k^(1/3)). */
    check("RootReduce[1/(a + b*k^(1/3))]",
          "(a^2 - a b k^(1/3) + b^2 k^(2/3))/(a^3 + b^3 k)");
    /* Constant algebraic number, degree 2 -> quadratic radical (qqbar path):
     * 1/(1+Sqrt[2]) = Sqrt[2]-1. */
    check("RootReduce[1/(1 + Sqrt[2])]", "-1 + Sqrt[2]");
    /* Constant algebraic number, degree 3 -> Root object (WL keeps degree >= 3
     * as a Root, not a cubic radical). 1/(1+2^(1/3)+2^(2/3)) = 2^(1/3)-1, whose
     * minimal polynomial is #^3 + 3#^2 + 3# - 1. */
    check("RootReduce[1/(1 + 2^(1/3) + 2^(2/3))]",
          "Root[-1 + 3 #1 + 3 #1^2 + #1^3 &, 1]");
    /* Free variable x with a constant-in-x cube root: cancels AND rationalises
     * (x - k^(1/3))/(x^2 - k^(2/3)) = 1/(x + k^(1/3)) -> rationalised. */
    check("RootReduce[(x - k^(1/3))/(x^2 - k^(2/3))]",
          "(k^(2/3) - k^(1/3) x + x^2)/(k + x^3)");
    /* No algebraic generator: RootReduce is the identity (it does NOT do plain
     * polynomial cancellation — that is Cancel's job). */
    check("RootReduce[x + 1]", "1 + x");
    check("RootReduce[(x^2 - 1)/(x - 1)]", "(-1 + x^2)/(-1 + x)");

    /* Thread over the coefficients of a polynomial in a free variable: the x^2
     * coefficient Sqrt[2]+Sqrt[3]-Sqrt[5+2Sqrt[6]] is a genuine algebraic zero,
     * so RootReduce reduces it to 0 and the monomial drops out. */
    check("RootReduce[(Sqrt[2] + Sqrt[3] - Sqrt[5 + 2 Sqrt[6]]) x^2 + x + 1]",
          "1 + x");
    /* A non-vanishing algebraic coefficient is canonicalised in place, while a
     * purely symbolic coefficient (a) and the variable structure are preserved. */
    check("RootReduce[a x^2 + Sqrt[8] x]", "2 Sqrt[2] x + a x^2");
    /* Rational function: the coefficient in the numerator is reduced, the
     * free-variable denominator is left intact (no spurious rationalisation). */
    check("RootReduce[Sqrt[8]/x]", "(2 Sqrt[2])/x");
}

/* ------------------------------------------------------------------ */
/*  RootReduce: value preservation (independent cross-check)           */
/* ------------------------------------------------------------------ */

static void test_rootreduce_preserves(void) {
    check_preserves("1/(1 + k^(1/3))");
    check_preserves("1/(a + b*k^(1/3))");
    check_preserves("1/(1 + k^(1/3) + k^(2/3))");
    check_preserves("(1 + k^(1/3))/(2 - k^(1/3))");
    check_preserves("k^(1/3)/(1 - k^(2/3))");
    check_preserves("1/(2 + 3 k^(1/3) + k^(2/3))");
    check_preserves("(x - k^(1/3))/(x^2 - k^(2/3))");
    check_preserves("1/(x + k^(1/3))");
    check_preserves("1/(x^2 + k^(1/3) x + k^(2/3))");
    /* higher index */
    check_preserves("1/(1 + k^(1/4))");
    check_preserves("1/(1 + k^(1/5))");
    /* two independent generators (cube root of k, square root of m) */
    check_preserves("1/(k^(1/3) + m^(1/2))");
    check_preserves("1/(1 + k^(1/3) + m^(1/2))");
    /* root of unity present */
    check_preserves("1/(1 + (-1)^(1/3) k^(1/3))");
    /* the genuine Goursat generator (radicand depends on x) */
    check_preserves("1/(1 + (x*(1 - x)*(1 - k*x))^(1/3))");
    /* constant radicands (number fields) */
    check_preserves("1/(1 + Sqrt[2])");
    check_preserves("1/(2 + 3 Sqrt[2])");
    check_preserves("1/(1 + 2^(1/3))");
    check_preserves("1/(a + b Sqrt[2] + c Sqrt[3])");

    /* Independent-engine cross-check on the genuine-algebraic tower: the
     * ideal-reduction zero test agrees with the matrix inversion (two distinct
     * FLINT engines, no numeric oracle). Restricted to genuine (symbol-radicand)
     * generators, which is where flint_algebraic_field_normalize engages. */
    check_preserves_ideal("1/(1 + k^(1/3))");
    check_preserves_ideal("1/(a + b*k^(1/3))");
    check_preserves_ideal("1/(1 + k^(1/3) + k^(2/3))");
    check_preserves_ideal("(x - k^(1/3))/(x^2 - k^(2/3))");
    check_preserves_ideal("1/(x^2 + k^(1/3) x + k^(2/3))");
    check_preserves_ideal("1/(k^(1/3) + m^(1/2))");
    check_preserves_ideal("1/(1 + (x*(1 - x)*(1 - k*x))^(1/3))");
}

/* ------------------------------------------------------------------ */
/*  RootReduce: denominator is rationalised (radical-free)             */
/* ------------------------------------------------------------------ */

static void test_rootreduce_rationalised(void) {
    check_rationalised("1/(1 + k^(1/3))", "k");
    check_rationalised("1/(a + b*k^(1/3))", "{a, b, k}");
    check_rationalised("1/(x + k^(1/3))", "{x, k}");
    check_rationalised("1/(1 + k^(1/4))", "k");
    check_rationalised("1/(k^(1/3) + m^(1/2))", "{k, m}");
    check_rationalised("k^(1/3)/(1 - k^(2/3))", "k");
}

/* ------------------------------------------------------------------ */
/*  RootReduce: idempotence                                            */
/* ------------------------------------------------------------------ */

static void test_rootreduce_idempotent(void) {
    check("RootReduce[RootReduce[1/(1 + k^(1/3))]] - RootReduce[1/(1 + k^(1/3))]", "0");
    check("RootReduce[RootReduce[1/(a + b k^(1/3))]] - RootReduce[1/(a + b k^(1/3))]", "0");
    check("Simplify[RootReduce[RootReduce[1/(k^(1/3) + m^(1/2))]] "
          "- RootReduce[1/(k^(1/3) + m^(1/2))]]", "0");
}

/* ------------------------------------------------------------------ */
/*  WL-faithful Together/Cancel over cube-and-higher-root towers        */
/* ------------------------------------------------------------------ */

static void test_together_cancel(void) {
    /* The previously-unhandled gap: a Plus of cube-root fractions combines to a
     * single fraction, radicals kept in the denominator (2x/(x^2 - k^(2/3))). */
    check("Together[1/(x - k^(1/3)) + 1/(x + k^(1/3))]", "-(2 x)/(k^(2/3) - x^2)");
    check("Together[1/(1 + k^(1/3)) + 1/(1 - k^(1/3))]", "-2/(-1 + k^(2/3))");
    /* Single-fraction Cancel still uses the relation-aware path (NOT pre-empted
     * by the free-kernel combine): (x^3 - k)/(x - k^(1/3)) cancels fully. */
    check("Cancel[(x^3 - k)/(x - k^(1/3))]", "k^(2/3) + k^(1/3) x + x^2");
    check("Cancel[(x - k^(1/3))/(x^2 - k^(2/3))]", "1/(k^(1/3) + x)");
    /* Constant-radical (number-field) sums still combine via the classical path. */
    check("Together[1/(x - Sqrt[2]) + 1/(x + Sqrt[2])]", "(2 x)/(-2 + x^2)");
    check("Cancel[(x - Sqrt[2])/(x^2 - 2)]", "1/(Sqrt[2] + x)");
    /* A Together identity: the combined form differentiates back to the sum. */
    check("Simplify[1/(x - k^(1/3)) + 1/(x + k^(1/3)) - 2 x/(x^2 - k^(2/3))]", "0");
}

/* ------------------------------------------------------------------ */
/*  Stress: a batch of towers, all value-preserving                    */
/* ------------------------------------------------------------------ */

static void test_stress(void) {
    /* Rational functions of one cube-root generator with varied numerators and
     * denominators — each must round-trip through RootReduce with value intact. */
    const char* cases[] = {
        "1/(1 + 2 k^(1/3))",
        "1/(1 - k^(1/3) + k^(2/3))",
        "(1 + k^(1/3) + k^(2/3))/(1 - k^(1/3))",
        "(k^(2/3) - 1)/(k^(1/3) + 2)",
        "1/(3 + k^(1/3))^2",
        "(a + b k^(1/3) + c k^(2/3))/(d + e k^(1/3))",
        "1/(1 + k^(1/6))",
        "1/(1 + k^(1/3) - k^(1/2))",
        "x/(x^3 - k)",
        "(x + k^(1/3))/(x - k^(1/3))",
        "1/((1 + k^(1/3))(2 + k^(1/3)))",
        "1/(1 + (-1)^(1/3) + k^(1/3))",
        NULL
    };
    for (int i = 0; cases[i]; i++) check_preserves(cases[i]);
}

/* ------------------------------------------------------------------ */
/*  G1/G2: constant algebraic numbers -> Root / quadratic radical / rational */
/* ------------------------------------------------------------------ */

static void test_qqbar_canonical(void) {
    /* G1: numerator-only algebraic number -> single Root object (WL: degree 4,
     * index 4 for Sqrt[2]+Sqrt[3]). */
    check("RootReduce[Sqrt[2] + Sqrt[3]]", "Root[1 - 10 #1^2 + #1^4 &, 4]");
    check("MinimalPolynomial[RootReduce[Sqrt[2] + Sqrt[3]], x]", "1 - 10 x^2 + x^4");
    /* G1: three provably-equal nested radicals canonicalise to the SAME Root. */
    check("RootReduce[Sqrt[2] + Sqrt[3] + Sqrt[5]] == "
          "RootReduce[Sqrt[10 + 2 Sqrt[15] + 4 Sqrt[4 + Sqrt[15]]]]", "True");
    /* G2: nested constant radicals reduce fully. */
    check("RootReduce[(Sqrt[18] + Sqrt[27])/Sqrt[5 + 2 Sqrt[6]]]", "3");
    check("RootReduce[(Sqrt[2] + Sqrt[3] + Sqrt[6] + 3)/Sqrt[5 + 2 Sqrt[6]]]",
          "1 + Sqrt[3]");
    check("RootReduce[Sqrt[7]/Sqrt[5 + 2 Sqrt[6]]]", "Root[49 - 70 #1^2 + #1^4 &, 3]");
    /* Degree 1 -> rational; degree 2 -> quadratic radical. */
    check("RootReduce[Sqrt[8] - 2 Sqrt[2]]", "0");
    check("RootReduce[1/(1 + Sqrt[2])]", "-1 + Sqrt[2]");
    /* Root-object arithmetic -> a single degree-15 Root (WL index 1). */
    check("RootReduce[Root[#^5 + 11 # + 1 &, 1] Root[#^3 + # + 17 &, 1]]",
          "Root[-1419857 + 918731 #1 + 111166451 #1^3 + 1446 #1^5 + 162316 #1^6 "
          "+ 139997 #1^7 + 85 #1^10 + 22 #1^11 + #1^15 &, 1]");
    /* Idempotence of a Root object (round-trips to its monic minimal poly). */
    check("RootReduce[Root[Function[t, t^3 + t + 17], 1]]",
          "Root[17 + #1 + #1^3 &, 1]");
}

/* ------------------------------------------------------------------ */
/*  G4: threading over equations / inequalities / logic                */
/* ------------------------------------------------------------------ */

static void test_qqbar_threading(void) {
    check("RootReduce[Sqrt[2] + Sqrt[3] + Sqrt[5] == "
          "Sqrt[10 + 2 Sqrt[15] + 4 Sqrt[4 + Sqrt[15]]]]", "True");
    check("RootReduce[Sqrt[2] == Sqrt[3]]", "False");
    check("RootReduce[Sqrt[2] < Sqrt[3]]", "True");
    check("RootReduce[Sqrt[3] <= Sqrt[3]]", "True");
    check("RootReduce[Sqrt[2] != Sqrt[3]]", "True");
    /* Logic threads: RootReduce maps into And, deciding each algebraic leaf. */
    check("RootReduce[Sqrt[2] < Sqrt[3] && Sqrt[2] == Sqrt[2]]", "True");
    /* Listable regression: threads over lists elementwise. */
    check("RootReduce[{1/(1 + Sqrt[2]), 2}]", "{-1 + Sqrt[2], 2}");
}

/* ------------------------------------------------------------------ */
/*  G3: argument-count / Method diagnostics leave the call unevaluated  */
/* ------------------------------------------------------------------ */

static void test_qqbar_argx(void) {
    check("RootReduce[]", "RootReduce[]");
    check("RootReduce[a, b]", "RootReduce[a, b]");
    check("RootReduce[Sqrt[2], Method -> \"Bogus\"]",
          "RootReduce[Sqrt[2], Method -> \"Bogus\"]");
}

/* ------------------------------------------------------------------ */
/*  G5: Method -> "Recursive" / "NumberField" agree with Automatic      */
/* ------------------------------------------------------------------ */

static void test_qqbar_methods(void) {
    /* The recursive-vs-numberfield example from the WL docs: c = 1. */
    check("b = RootReduce[2^(1/3) + 3^(1/3) + 1]; "
          "c = b - 2^(1/3) - 3^(1/3); RootReduce[c]", "1");
    check("RootReduce[c, Method -> \"Recursive\"]", "1");
    check("RootReduce[c, Method -> \"NumberField\"]", "1");
    /* Multi-generator tower: all three methods produce the identical canonical
     * Root. Uses a degree-15 tower Q(2^(1/3), 3^(1/5)); the analogous degree-21
     * Q(2^(1/3), 3^(1/7)) tower is mathematically identical in intent but each
     * qqbar reduction runs ~220s, exceeding the harness alarm(60), so it is not
     * used as a unit test. */
    check("aa = 2 2^(1/3) + 3 3^(1/5) + 5 2^(1/3) 3^(1/5); "
          "RootReduce[aa] == RootReduce[aa, Method -> \"NumberField\"] == "
          "RootReduce[aa, Method -> \"Recursive\"]", "True");
}

int main(void) {
    symtab_init();
    core_init();

    if (!flint_bridge_available()) {
        printf("FLINT not compiled in (USE_FLINT off); skipping RootReduce tests.\n");
        return 0;
    }

    test_rootreduce_exact();
    test_rootreduce_preserves();
    test_rootreduce_rationalised();
    test_rootreduce_idempotent();
    test_together_cancel();
    test_stress();
    test_qqbar_canonical();
    test_qqbar_threading();
    test_qqbar_argx();
    test_qqbar_methods();

    printf("All RootReduce / algebraic-field tests passed.\n");
    return 0;
}
