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
    /* A14: a degree-5 Root generator whose membership relation needs > 64-bit
     * precision. Q(Root[2869 x^5 + ...]) == Q(Root[x^5-x-1]); before the
     * ToNumberField-path escalation this declined (unevaluated). Value-preserving
     * both directions (MATHILDA_DIVERGENCES A14). */
    check_zero("ToNumberField[Root[-1 + 15 #1 - 80 #1^2 + 160 #1^3 + 2869 #1^5 &, 1], "
               "Root[-1 - #1 + #1^5 &, 1]] "
               "- Root[-1 + 15 #1 - 80 #1^2 + 160 #1^3 + 2869 #1^5 &, 1]");
    check_zero("ToNumberField[Root[-1 - #1 + #1^5 &, 1], "
               "Root[-1 + 15 #1 - 80 #1^2 + 160 #1^3 + 2869 #1^5 &, 1]] "
               "- Root[-1 - #1 + #1^5 &, 1]");
    /* Degree-8 nested-radical compositum: 64-bit field-membership precision is
     * too low, so ToNumberField escalates precision (degree-gated, only for
     * degree > 6).  Before the escalation this declined (unevaluated); the guard
     * is that it now constructs a single theta and preserves value.  (Regression
     * guard for the P8/P4 degree-16 precision fix.) */
    check("Head[ToNumberField[{Sqrt[5], Sqrt[1 + Sqrt[5]], Sqrt[1 - Sqrt[5]]}]]", "List");
    check("Abs[N[ToNumberField[{Sqrt[5], Sqrt[1 + Sqrt[5]], Sqrt[1 - Sqrt[5]]}][[2]] "
          "- Sqrt[1 + Sqrt[5]], 40]] < 10^-30", "True");
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

/* An AlgebraicNumber is a CONSTANT of its number field, not a polynomial
 * variable: Variables[] must not surface it, and the extension machinery must
 * not mine the field's label theta out of it as a tower generator (which made
 * PolynomialGCD[.., Extension -> Automatic] answer 1 on a pair with a common
 * factor). See MATHILDA_DIVERGENCES A26. */
static void test_not_a_polynomial_variable(void) {
    check("Variables[AlgebraicNumber[Sqrt[2], {0, 1}] + x]", "{x}");
    check("Variables[AlgebraicNumber[Sqrt[2], {0, 1}] x^2 + 3]", "{x}");
    check("PolynomialGCD[Expand[(x + AlgebraicNumber[Sqrt[2], {0, 1}]) (x + 1)], "
          "Expand[(x + AlgebraicNumber[Sqrt[2], {0, 1}]) (x + 2)]]",
          "AlgebraicNumber[Sqrt[2], {0, 1}] + x");
    check("PolynomialGCD[Expand[(x + AlgebraicNumber[Sqrt[2], {0, 1}]) (x + 1)], "
          "Expand[(x + AlgebraicNumber[Sqrt[2], {0, 1}]) (x + 2)], Extension -> Automatic]",
          "AlgebraicNumber[Sqrt[2], {0, 1}] + x");
}

/* Assert that PolynomialGCD[f, g] really is the gcd, given a known answer `d`.
 *
 * Asserting the PRINTED form would be asserting an associate: a gcd over a
 * field is defined only up to a constant of that field, and the paths here
 * disagree on which one they pick (flint_field_gcd is monic, the classical PRS
 * is not).  So assert the two facts that actually matter — the result divides
 * both operands, and the known answer divides the result.  Together those pin
 * it to `d` up to a field constant, which is the whole of the specification.
 *
 * Division is exact polynomial division in x rather than Cancel, because Cancel
 * deliberately declines on AlgebraicNumber coefficients under its default
 * Extension -> None (as Mathematica's does). */
static void check_gcd(const char* f, const char* g, const char* d) {
    /* Sized from the operands rather than fixed: f and g each appear TWICE in the
     * template, and the large-coefficient cases carry a ~900-digit literal, so a
     * fixed buffer would truncate the assertion into something that no longer
     * tests what it names. */
    size_t need = 2 * (strlen(f) + strlen(g)) + strlen(d) + 512;
    char* buf = malloc(need);
    ASSERT(buf != NULL);
    snprintf(buf, need,
        "Module[{gg = PolynomialGCD[%s, %s], dv},"
        " dv = Function[{p, q}, Expand[p - PolynomialQuotient[p, q, x] q] === 0];"
        " dv[%s, gg] && dv[%s, gg] && dv[gg, %s]]", f, g, f, g, d);
    Expr* e = eval_str(buf);
    char* s = expr_to_string(e);
    if (strcmp(s, "True") != 0) {
        snprintf(buf, need, "PolynomialGCD[%s, %s]", f, g);
        Expr* got = eval_str(buf);
        char* gs = expr_to_string(got);
        fprintf(stderr, "FAIL(gcd): PolynomialGCD[%s, %s]\n"
                        "  expected an associate of: %s\n  got: %s\n", f, g, d, gs);
        free(gs); expr_free(got);
        exit(1);
    }
    free(s);
    expr_free(e);
    free(buf);
}

/* MULTIVARIATE gcd over a number field — MATHILDA_DIVERGENCES A26.
 *
 * Every engine other than flint_field_gcd is univariate (gr_poly over an antic
 * nf_t), so these fell through to the classical pseudo-remainder PRS, whose
 * content is INTEGER content.  A coefficient in K contributes content 1, both
 * operands stay non-primitive over K, and the first pseudo-remainder
 * lc(B)*A - lc(A)*B vanishes identically whenever the two share a factor and
 * agree in degree -- so the loop returned the SECOND OPERAND, which is not a
 * common divisor at all.  There was no multivariate test here, which is exactly
 * how it went unnoticed. */
static void test_multivariate_number_field_gcd(void) {
    /* the two A26 repros: both answered (x + a y)(x + 2) before */
    check_gcd("Expand[(x + AlgebraicNumber[Sqrt[2], {0, 1}] y) (x + 1)]",
              "Expand[(x + AlgebraicNumber[Sqrt[2], {0, 1}] y) (x + 2)]",
              "x + AlgebraicNumber[Sqrt[2], {0, 1}] y");
    check_gcd("x + AlgebraicNumber[Sqrt[2], {0, 1}] y",
              "Expand[(x + AlgebraicNumber[Sqrt[2], {0, 1}] y) (x + 2)]",
              "x + AlgebraicNumber[Sqrt[2], {0, 1}] y");
    /* operand order must not matter (it returned whichever was second) */
    check_gcd("Expand[(x + AlgebraicNumber[Sqrt[2], {0, 1}] y) (x + 2)]",
              "Expand[(x + AlgebraicNumber[Sqrt[2], {0, 1}] y) (x + 1)]",
              "x + AlgebraicNumber[Sqrt[2], {0, 1}] y");
    /* the same in the radical spelling, with and without the option */
    check_gcd("Expand[(x + Sqrt[2] y) (x + 1)]", "Expand[(x + Sqrt[2] y) (x + 2)]",
              "x + Sqrt[2] y");

    /* Cofactors that themselves carry algebraic constants.  This is the shape
     * the Phase D tower path gets wrong in the other direction -- it computes
     * the Q[gamma, x, y]-GCD rather than the Q(gamma)[x, y]-GCD and answered 1. */
    check_gcd("Expand[(x^3 + Sqrt[2] x y + y^2 + 1) (x^2 + Sqrt[2] y + 3)]",
              "Expand[(x^3 + Sqrt[2] x y + y^2 + 1) (x^2 + 2 Sqrt[2] y - 1)]",
              "x^3 + Sqrt[2] x y + y^2 + 1");
    check_gcd("Expand[(x^3 + AlgebraicNumber[Sqrt[2], {0, 1}] x y + y^2 + 1) "
              "(x^2 + AlgebraicNumber[Sqrt[2], {0, 1}] y + 3)]",
              "Expand[(x^3 + AlgebraicNumber[Sqrt[2], {0, 1}] x y + y^2 + 1) "
              "(x^2 + 2 AlgebraicNumber[Sqrt[2], {0, 1}] y - 1)]",
              "x^3 + AlgebraicNumber[Sqrt[2], {0, 1}] x y + y^2 + 1");

    /* a compositum (no prime keeps its minpoly irreducible -- the residue field
     * has to be SPLIT and the components CRT'd back), a degree-3 and a degree-4
     * field, a Root generator, Q(i), and three variables */
    check_gcd("Expand[(x + Sqrt[2] y + Sqrt[3]) (x + 1)]",
              "Expand[(x + Sqrt[2] y + Sqrt[3]) (x + 2)]", "x + Sqrt[2] y + Sqrt[3]");
    check_gcd("Expand[(x + 2^(1/3) y) (x + 1)]", "Expand[(x + 2^(1/3) y) (x + 2)]",
              "x + 2^(1/3) y");
    check_gcd("Expand[(x + 2^(1/4) y) (x + 1)]", "Expand[(x + 2^(1/4) y) (x + 2)]",
              "x + 2^(1/4) y");
    check_gcd("Expand[(x + Root[#^3 - # - 1 &, 1] y) (x + 1)]",
              "Expand[(x + Root[#^3 - # - 1 &, 1] y) (x + 2)]",
              "x + Root[#^3 - # - 1 &, 1] y");
    check_gcd("Expand[(x + I y) (x + 1)]", "Expand[(x + I y) (x + 2)]", "x + I y");
    check_gcd("Expand[(x + Sqrt[2] y + z) (x + 1)]",
              "Expand[(x + Sqrt[2] y + z) (x + 2)]", "x + Sqrt[2] y + z");
    /* a NESTED radical: its base is itself a Plus, which the constant/parametric
     * pre-filter has to accept (the head symbol `Plus` is not a free variable) */
    check_gcd("Expand[(x + Sqrt[1 + Sqrt[2]] y) (x + 1)]",
              "Expand[(x + Sqrt[1 + Sqrt[2]] y) (x + 2)]", "x + Sqrt[1 + Sqrt[2]] y");
    check_gcd("Expand[(x + 2^(1/6) y) (x + 1)]", "Expand[(x + 2^(1/6) y) (x + 2)]",
              "x + 2^(1/6) y");

    /* A PARAMETRIC radical is a rational function field, not a number field: it
     * must keep going to flint_parametric_sqrt_gcd rather than being swallowed
     * here (and must not be caught by the classical path's 1-guard either). */
    check_gcd("Expand[(x + Sqrt[k] y) (x + 1)]", "Expand[(x + Sqrt[k] y) (x + 2)]",
              "x + Sqrt[k] y");

    /* Coprime operands must still answer 1 -- the certificate has to reject a
     * spurious common factor, not just confirm a real one. */
    check("PolynomialGCD[Expand[(x + Sqrt[2] y) (x + 1)], Expand[(x + Sqrt[3] y) (x + 2)]]", "1");
    check("PolynomialGCD[Expand[(x + AlgebraicNumber[Sqrt[2], {0, 1}] y) (x + 1)], "
          "Expand[(x + AlgebraicNumber[Sqrt[2], {0, 3}] y) (x + 2)]]", "1");

    /* The radical spelling must come back in radicals, not in qqbar's Root
     * spelling of the primitive element (the answer is mapped through the
     * product basis of the caller's own atoms to keep this true). */
    check("PolynomialGCD[Expand[(x + 2^(1/3) y) (x + 1)], Expand[(x + 2^(1/3) y) (x + 2)]]",
          "x + 2^(1/3) y");
    check("PolynomialGCD[Expand[(x + Sqrt[2] y + Sqrt[3]) (x + 1)], "
          "Expand[(x + Sqrt[2] y + Sqrt[3]) (x + 2)]]", "Sqrt[3] + x + Sqrt[2] y");
}

/* The defects the v0.231 stress pass turned up.  Each of these ANSWERED -- with 1,
 * or with an unreadable spelling -- rather than failing visibly, which is why none
 * was caught by the v0.230 tests. */
static void test_field_gcd_stress_regressions(void) {
    /* 1. DEGREE 6 declined for every radical generator, so the gcd was lost and
     * the caller's post-check answered 1.  The cause was upstream in the qqbar
     * compositum: it picked a primitive element by TRIAL membership, and
     * qqbar_express_in_field_esc refuses to escalate past 64 bits of working
     * precision when the generator has degree <= 6.  Degrees 2-5 resolve inside
     * 64 bits and 7+ are allowed to escalate, so 6 alone failed -- for every c in
     * the search, hence the whole field.  The primitive element is now chosen by
     * degree, which needs no membership test at all. */
    check_gcd("Expand[(x^3 + 2^(1/6) x y + y^2 + 1) (x^2 + 2^(1/6) y + 3)]",
              "Expand[(x^3 + 2^(1/6) x y + y^2 + 1) (x^2 + 2 2^(1/6) y - 1)]",
              "x^3 + 2^(1/6) x y + y^2 + 1");
    check_gcd("Expand[(x^3 + 3^(1/6) x y + y^2 + 1) (x^2 + 3^(1/6) y + 3)]",
              "Expand[(x^3 + 3^(1/6) x y + y^2 + 1) (x^2 + 2 3^(1/6) y - 1)]",
              "x^3 + 3^(1/6) x y + y^2 + 1");
    /* the same degree 6 reached as a genuine compositum of two generators */
    check_gcd("Expand[(x^2 + (Sqrt[2] + 2^(1/3)) y + 1) (x + y + 1)]",
              "Expand[(x^2 + (Sqrt[2] + 2^(1/3)) y + 1) (x - y + 2)]",
              "x^2 + (Sqrt[2] + 2^(1/3)) y + 1");
    /* and ToNumberField itself, which returned unevaluated for this atom set */
    check("Head[ToNumberField[{2^(1/6), 2^(1/3)}]]", "List");

    /* 2. The COMPOSITUM render-back was abandoned whenever the operands mentioned
     * more atoms than a basis needs.  Expand folds Sqrt[2] Sqrt[3] into Sqrt[6],
     * so a Q(sqrt2, sqrt3) problem has atoms {sqrt2, sqrt3, sqrt6} whose product
     * basis would hold 8 members against [K:Q] = 4; the overshoot was read as
     * "these atoms do not span" and the answer came back as a degree-4 Root per
     * coefficient.  Atoms are now selected greedily until the product reaches n. */
    check("PolynomialGCD[Expand[(x^3 + Sqrt[2] x y + y^2 + 1) (x^2 + Sqrt[3] y + 3)], "
          "Expand[(x^3 + Sqrt[2] x y + y^2 + 1) (x^2 + Sqrt[2] y - 1)]]",
          "1 + x^3 + Sqrt[2] x y + y^2");

    /* 3. The COEFFICIENT CEILING.  64 primes of 29 bits is 1856 bits of modulus,
     * and rational reconstruction needs about twice the coefficient size, so the
     * engine declined -- and the caller answered 1 -- above ~831 bits.  Measured,
     * 10^250 passed and 10^300 did not.  Both must now hold, and well beyond. */
    check_gcd("Expand[(x^3 + 10^301 Sqrt[2] x y + y^2 + 1) (x^2 + Sqrt[2] y + 3)]",
              "Expand[(x^3 + 10^301 Sqrt[2] x y + y^2 + 1) (x^2 + 2 Sqrt[2] y - 1)]",
              "x^3 + 10^301 Sqrt[2] x y + y^2 + 1");
    check_gcd("Expand[(x^3 + 10^903 Sqrt[2] x y + y^2 + 1) (x^2 + Sqrt[2] y + 3)]",
              "Expand[(x^3 + 10^903 Sqrt[2] x y + y^2 + 1) (x^2 + 2 Sqrt[2] y - 1)]",
              "x^3 + 10^903 Sqrt[2] x y + y^2 + 1");
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
    test_not_a_polynomial_variable();
    test_multivariate_number_field_gcd();
    test_field_gcd_stress_regressions();
    printf("test_algebraicnumber: all passed\n");
    return 0;
}
