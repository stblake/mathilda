/* test_algebraicnumber.c — AlgebraicNumber and ToNumberField.
 *
 * Covers canonicalisation (algebraic-integer reduction, power-basis reduction,
 * over-length folding, rational collapse), ToNumberField in every argument
 * form, number-field arithmetic (+,*,/,^), N to high precision, and the
 * Re/Im/Abs/Round/Less/Equal/NumericQ operations. Correctness without a numeric
 * oracle is cross-checked with the independent RootReduce zero test
 * (RootReduce[result - original] == 0).
 *
 * Requires FLINT (the qqbar engine); the whole suite SKIPs cleanly when off.
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

/* Assert the short-form printed value of `input` equals `expected`. */
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

/* Assert the printed value of `input` starts with `prefix` (arbitrary-precision
 * numeric outputs). */
static void check_prefix(const char* input, const char* prefix) {
    Expr* e = eval_str(input);
    char* s = expr_to_string(e);
    if (strncmp(s, prefix, strlen(prefix)) != 0) {
        fprintf(stderr, "FAIL(prefix): %s\n  expected prefix: %s\n  got: %s\n",
                input, prefix, s);
        exit(1);
    }
    free(s);
    expr_free(e);
}

/* Value-preservation: assert RootReduce[expr] prints as "0". */
static void check_zero(const char* expr) {
    char buf[1024];
    snprintf(buf, sizeof buf, "RootReduce[%s]", expr);
    Expr* e = eval_str(buf);
    char* s = expr_to_string(e);
    if (strcmp(s, "0") != 0) {
        fprintf(stderr, "FAIL(nonzero): RootReduce[%s] = %s (expected 0)\n", expr, s);
        exit(1);
    }
    free(s);
    expr_free(e);
}

/* ---------------- Canonicalisation ----------------------------------------- */

static void test_canonicalisation(void) {
    /* Already canonical (Root generator, coeff length = degree). */
    check("AlgebraicNumber[Root[#^3+#+1&,3],{1,2,1}]",
          "AlgebraicNumber[Root[1 + #1 + #1^3 &, 3], {1, 2, 1}]");
    /* Rational generator collapses to the rational value. */
    check("AlgebraicNumber[3,{1,2}]", "7");
    check("AlgebraicNumber[5/2,{3,2}]", "8");
    /* Non-algebraic-integer generator rescales (lc = 2): (1+I)/2 -> 1+I. */
    check("AlgebraicNumber[(1+I)/2,{1,3}]", "AlgebraicNumber[1 + I, {1, 3/2}]");
    /* Nested radical generator becomes a Root object. */
    check("AlgebraicNumber[Sqrt[Sqrt[2]+1],{1,2,1,2}]",
          "AlgebraicNumber[Root[-1 - 2 #1^2 + #1^4 &, 2], {1, 2, 1, 2}]");
    /* Root generator with lc = 5: rescales to a monic minimal polynomial. */
    check("AlgebraicNumber[Root[5#^5+11#+1&,1],{1,1,2}]",
          "AlgebraicNumber[Root[625 + 1375 #1 + #1^5 &, 1], {1, 1/5, 2/25, 0, 0}]");
    /* Nested AlgebraicNumber generator. */
    check("AlgebraicNumber[AlgebraicNumber[Root[-3+#1^3&,1],{1,2,1}],{1,1,2}]",
          "AlgebraicNumber[Root[-16 - 15 #1 - 3 #1^2 + #1^3 &, 1], {1, 1, 2}]");
    /* Sum of radicals -> degree-4 Root, coeffs padded. */
    check("AlgebraicNumber[Sqrt[2]+Sqrt[5],{1,1/2}]",
          "AlgebraicNumber[Root[9 - 14 #1^2 + #1^4 &, 4], {1, 1/2, 0, 0}]");
    /* Short coefficient list padded to the degree. */
    check("AlgebraicNumber[3^(1/5),{1,2,1}]",
          "AlgebraicNumber[Root[-3 + #1^5 &, 1], {1, 2, 1, 0, 0}]");
    /* Over-length coefficient list folds modulo the minimal polynomial. */
    check("AlgebraicNumber[3^(1/5),{1,2,1,3,3,1}]",
          "AlgebraicNumber[Root[-3 + #1^5 &, 1], {4, 2, 1, 3, 3}]");
    /* Empty coefficient list is 0; every higher coeff zero collapses. */
    check("AlgebraicNumber[Sqrt[2],{}]", "0");
    check("AlgebraicNumber[Sqrt[2],{5,0}]", "5");
    /* Idempotence: re-evaluating a canonical object is a fixed point. */
    check("AlgebraicNumber[Root[#^3+#+1&,3],{1,2,1}] === "
          "AlgebraicNumber[Root[#^3+#+1&,3],{1,2,1}]", "True");
}

/* ---------------- Value preservation --------------------------------------- */

static void test_value_preservation(void) {
    check_zero("AlgebraicNumber[(1+I)/2,{1,3}] - (1 + 3 (1+I)/2)");
    check_zero("AlgebraicNumber[Root[5#^5+11#+1&,1],{1,1,2}] - "
               "(1 + Root[5#^5+11#+1&,1] + 2 Root[5#^5+11#+1&,1]^2)");
    check_zero("AlgebraicNumber[Sqrt[2]+Sqrt[5],{1,1/2}] - (1 + (Sqrt[2]+Sqrt[5])/2)");
    check_zero("AlgebraicNumber[3^(1/5),{1,2,1,3,3,1}] - "
               "(1 + 2 3^(1/5) + 3^(2/5) + 3 3^(3/5) + 3 3^(4/5) + 3)");
}

/* ---------------- ToNumberField -------------------------------------------- */

static void test_to_number_field(void) {
    check("ToNumberField[Sqrt[2], 2^(1/4)]",
          "AlgebraicNumber[Root[-2 + #1^4 &, 2], {0, 0, 1, 0}]");
    check("ToNumberField[2, 1/2]", "2");                  /* field is Q */
    /* a not in Q(theta): unevaluated. */
    check("ToNumberField[Sqrt[3], Sqrt[2]]", "ToNumberField[Sqrt[3], Sqrt[2]]");
    /* a in Q(theta), value-preserving. */
    check_zero("ToNumberField[Sqrt[2], 2^(1/4)] - Sqrt[2]");
    check_zero("ToNumberField[Root[-25-24#1-12#1^2+4#1^3&,1], Root[-3+2#1^3&,1]] "
               "- Root[-25-24#1-12#1^2+4#1^3&,1]");
    /* ToNumberField[x]: explicit AlgebraicNumber, value-preserving. */
    check_zero("ToNumberField[Sqrt[Sqrt[2]+Sqrt[3]]] - Sqrt[Sqrt[2]+Sqrt[3]]");
    /* Common field: primitive element of Q(Sqrt[2], I) is Sqrt[2]+I. */
    check("ToNumberField[{Sqrt[2],I},All][[1,1]]", "Root[9 - 2 #1^2 + #1^4 &, 4]");
    /* Common field, All: each element value-preserving, shared generator. */
    check_zero("ToNumberField[{Sqrt[3],(1+I Sqrt[3])/2}][[1]] - Sqrt[3]");
    check_zero("ToNumberField[{Sqrt[3],(1+I Sqrt[3])/2}][[2]] - (1+I Sqrt[3])/2");
    check_zero("ToNumberField[{AlgebraicNumber[Root[1-10#1^2+#1^4&,4],{0,-9/2,0,1/2}],"
               "Sqrt[5]},All][[1]] - Sqrt[2]");
    check_zero("ToNumberField[{AlgebraicNumber[Root[1-10#1^2+#1^4&,4],{0,-9/2,0,1/2}],"
               "Sqrt[5]},All][[2]] - Sqrt[5]");
    /* ToNumberField[{a..}, theta]: express each in Q(theta). */
    check_zero("ToNumberField[{1, Sqrt[2]}, Sqrt[2]][[2]] - Sqrt[2]");
}

/* ---------------- Field arithmetic ----------------------------------------- */

static void test_arithmetic(void) {
    check("1 + AlgebraicNumber[Root[#^3+#+1&,3],{1,2,1}]^2",
          "AlgebraicNumber[Root[1 + #1 + #1^3 &, 3], {-2, -1, 5}]");
    check("AlgebraicNumber[Sqrt[2],{1,1/2}]+AlgebraicNumber[Sqrt[2],{1,2}]",
          "AlgebraicNumber[Sqrt[2], {2, 5/2}]");
    check("AlgebraicNumber[Sqrt[2],{1,1/2}]*AlgebraicNumber[Sqrt[2],{1,2}]",
          "AlgebraicNumber[Sqrt[2], {3, 5/2}]");
    check("1/AlgebraicNumber[Sqrt[2],{1,1/2}]", "AlgebraicNumber[Sqrt[2], {2, -1}]");
    check("AlgebraicNumber[Sqrt[2],{1,1/2}]^3", "AlgebraicNumber[Sqrt[2], {5/2, 7/4}]");
    check("3*AlgebraicNumber[Sqrt[2],{1,2}]", "AlgebraicNumber[Sqrt[2], {3, 6}]");
    check("AlgebraicNumber[Sqrt[2],{1,2}]^0", "1");
    /* Cancellation collapses to a rational. */
    check("AlgebraicNumber[Sqrt[2],{1,1/2}] - AlgebraicNumber[Sqrt[2],{1,1/2}]", "0");
    /* Distinct generators do not combine. */
    check("AlgebraicNumber[Sqrt[2],{1,1}]+AlgebraicNumber[Sqrt[3],{1,1}]",
          "AlgebraicNumber[Sqrt[2], {1, 1}] + AlgebraicNumber[Sqrt[3], {1, 1}]");
    /* Arithmetic is value-preserving. */
    check_zero("(AlgebraicNumber[Sqrt[2],{1,1/2}]*AlgebraicNumber[Sqrt[2],{1,2}]) "
               "- (1+Sqrt[2]/2)(1+2 Sqrt[2])");
    check_zero("(1/AlgebraicNumber[Sqrt[2],{1,1/2}]) - 1/(1+Sqrt[2]/2)");
}

/* ---------------- N / operations ------------------------------------------- */

static void test_numeric_and_operations(void) {
    /* N of AlgebraicNumber[Sqrt[2] I, {1,-1}] = 1 - I Sqrt[2]. */
    check_prefix("N[AlgebraicNumber[Sqrt[2] I,{1,-1}]]", "1.0 - 1.41421");
    check_prefix("N[AlgebraicNumber[Sqrt[2] I,{1,-1}],50]",
                 "1.0 - 1.4142135623730950488016887242096980785696718753769");

    /* Real algebraic number a ~ 1.18: comparisons, Round, Re, Im, Abs. */
    const char* a = "AlgebraicNumber[Root[-1+#1+#1^3&,1],{0,-1,4}]";
    char buf[512];
    snprintf(buf, sizeof buf, "1 < %s", a);              check(buf, "True");
    snprintf(buf, sizeof buf, "%s < 1", a);              check(buf, "False");
    snprintf(buf, sizeof buf, "Round[%s]", a);           check(buf, "1");
    snprintf(buf, sizeof buf, "Re[%s]", a);
    check(buf, "AlgebraicNumber[Root[-1 + #1 + #1^3 &, 1], {0, -1, 4}]");
    snprintf(buf, sizeof buf, "Im[%s]", a);              check(buf, "0");
    snprintf(buf, sizeof buf, "Abs[%s]", a);
    check(buf, "AlgebraicNumber[Root[-1 + #1 + #1^3 &, 1], {0, -1, 4}]");
    /* Abs of a negative real algebraic number negates. */
    check("Abs[AlgebraicNumber[Sqrt[2],{0,-1}]]", "AlgebraicNumber[Sqrt[2], {0, 1}]");

    /* NumericQ / exact (in)equality. */
    check("NumericQ[AlgebraicNumber[Sqrt[2],{1,2}]]", "True");
    check("AlgebraicNumber[I,{0,1}] == AlgebraicNumber[I,{2,1}]", "False");
    check("AlgebraicNumber[I,{0,1}] == AlgebraicNumber[I,{0,1}]", "True");
    check("AlgebraicNumber[I,{0,1}] != AlgebraicNumber[I,{2,1}]", "True");

    /* RootReduce converts an AlgebraicNumber into a Root object. */
    check("RootReduce[AlgebraicNumber[Root[#^3+#+1&,3],{1,2,1}]]",
          "Root[-1 + 10 #1 - #1^2 + #1^3 &, 3]");
    /* E^(I Pi r) is recognised as a root of unity. */
    check("RootReduce[E^(Pi I/4)]", "Root[1 + #1^4 &, 4]");
    check_prefix("N[ToNumberField[E^(Pi I/4), I AlgebraicNumber[Sqrt[2],{1,2}]]]",
                 "0.707107 + 0.707107");
}

/* ---------------- Declines -------------------------------------------------- */

static void test_declines(void) {
    /* Symbolic generator: unevaluated. */
    check("AlgebraicNumber[x, {1,2}]", "AlgebraicNumber[x, {1, 2}]");
    /* Transcendental generator: unevaluated. */
    check("AlgebraicNumber[Pi, {1,2}]", "AlgebraicNumber[Pi, {1, 2}]");
}

int main(void) {
    symtab_init();
    core_init();
    if (!flint_bridge_available()) {
        printf("FLINT not compiled in (USE_FLINT off); skipping AlgebraicNumber tests.\n");
        return 0;
    }
    test_canonicalisation();
    test_value_preservation();
    test_to_number_field();
    test_arithmetic();
    test_numeric_and_operations();
    test_declines();
    printf("test_algebraicnumber: all passed\n");
    return 0;
}
