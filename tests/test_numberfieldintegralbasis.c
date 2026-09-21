/* test_numberfieldintegralbasis.c — NumberFieldIntegralBasis and AlgebraicIntegerQ.
 *
 * NumberFieldIntegralBasis[a] returns a Z-basis of the ring of integers O_K of
 * Q(a).  A returned basis is not unique, so correctness is checked by
 * representation-independent invariants rather than exact element values:
 *   - the basis has [K:Q] elements and begins with the rational integer 1;
 *   - every element is an algebraic integer (AlgebraicIntegerQ), as is every
 *     integer combination of them;
 *   - for canonical fields the elements match known values, cross-checked
 *     oracle-free with RootReduce[element - expected] == 0;
 *   - the maximal (non-monogenic) order is found (Q(Sqrt[5]) -> golden ratio,
 *     not Z[Sqrt[5]]; Q(Sqrt[2]+Sqrt[3]) -> a 4-element basis).
 * AlgebraicIntegerQ is tested directly across integers, rationals, radicals,
 * roots of unity, Root objects, and non-algebraic / non-scalar arguments.
 *
 * Requires FLINT (the qqbar + number-field engine); the suite SKIPs cleanly off.
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

/* Assert the printed value of `input` starts with `prefix` (numeric N output). */
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

/* Value-preservation without a numeric oracle: assert RootReduce[expr] == "0". */
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

/* ---------------- AlgebraicIntegerQ ---------------------------------------- */

static void test_algebraic_integer_q(void) {
    /* Rational integers are algebraic integers; non-integer rationals are not. */
    check("AlgebraicIntegerQ[0]", "True");
    check("AlgebraicIntegerQ[3]", "True");
    check("AlgebraicIntegerQ[-5]", "True");
    check("AlgebraicIntegerQ[1/2]", "False");
    check("AlgebraicIntegerQ[2/3]", "False");

    /* Radicals and roots of unity. */
    check("AlgebraicIntegerQ[Sqrt[2]]", "True");
    check("AlgebraicIntegerQ[2^(1/3)]", "True");
    check("AlgebraicIntegerQ[I]", "True");
    check("AlgebraicIntegerQ[E^(I Pi/4)]", "True");
    check("AlgebraicIntegerQ[(1 + Sqrt[5])/2]", "True");   /* golden ratio */

    /* Algebraic but not an algebraic integer (non-monic minimal polynomial). */
    check("AlgebraicIntegerQ[2^(1/3)/2]", "False");
    check("AlgebraicIntegerQ[Sqrt[5]/2]", "False");

    /* Root objects: monic defining polynomial -> yes, otherwise no. */
    check("AlgebraicIntegerQ[Root[#^3 + # + 1 &, 1]]", "True");
    check("AlgebraicIntegerQ[Root[3 #^3 + 2 &, 1]]", "False");

    /* Non-algebraic, symbolic, and non-scalar arguments are all False (the head
     * is not Listable, so a list is not an algebraic integer). */
    check("AlgebraicIntegerQ[Pi]", "False");
    check("AlgebraicIntegerQ[x]", "False");
    check("AlgebraicIntegerQ[{1, 2}]", "False");

    /* Wrong arity stays unevaluated. */
    check("AlgebraicIntegerQ[1, 2]", "AlgebraicIntegerQ[1, 2]");
}

/* ---------------- NumberFieldIntegralBasis: shape & Q --------------------- */

static void test_shape(void) {
    /* Q(a) = Q: the integral basis is {1}. */
    check("NumberFieldIntegralBasis[5]", "{1}");
    check("NumberFieldIntegralBasis[1/3]", "{1}");
    check("NumberFieldIntegralBasis[-2]", "{1}");

    /* Cardinality = field degree, and the basis begins with 1. */
    check("Length[NumberFieldIntegralBasis[2^(1/3)]]", "3");
    check("Length[NumberFieldIntegralBasis[Sqrt[2] + Sqrt[3]]]", "4");
    check("Length[NumberFieldIntegralBasis[E^(I Pi/4)]]", "4");
    check("First[NumberFieldIntegralBasis[2^(1/3)]]", "1");
    check("First[NumberFieldIntegralBasis[Sqrt[2] + Sqrt[3]]]", "1");
    check("First[NumberFieldIntegralBasis[E^(I Pi/4)]]", "1");
    check("First[NumberFieldIntegralBasis[Sqrt[2 + Sqrt[3]]]]", "1");
}

/* ---------------- Canonical fields: numeric + oracle-free ----------------- */

static void test_canonical_fields(void) {
    /* Q(2^(1/3)) is monogenic: the power basis {1, 2^(1/3), 2^(2/3)}. */
    check_prefix("N[NumberFieldIntegralBasis[2^(1/3)]]", "{1.0, 1.25992, 1.5874");
    check_zero("NumberFieldIntegralBasis[2^(1/3)][[2]] - 2^(1/3)");
    check_zero("NumberFieldIntegralBasis[2^(1/3)][[3]] - 2^(2/3)");

    /* Q(i) = Z[i]. */
    check_prefix("N[NumberFieldIntegralBasis[I]]", "{1.0, 0.0 + 1.0*I}");
    check_zero("NumberFieldIntegralBasis[I][[2]] - I");

    /* The cyclotomic field Q(zeta_8) = Z[zeta_8]: {1, zeta, zeta^2, zeta^3}. */
    check_prefix("N[NumberFieldIntegralBasis[E^(I Pi/4)]]", "{1.0, 0.707107");
    check("AlgebraicIntegerQ /@ NumberFieldIntegralBasis[E^(I Pi/4)]",
          "{True, True, True, True}");
    check_zero("NumberFieldIntegralBasis[E^(I Pi/4)][[2]] - E^(I Pi/4)");
    check_zero("NumberFieldIntegralBasis[E^(I Pi/4)][[3]] - I");
}

/* ---------------- Maximal (non-monogenic) orders -------------------------- */

static void test_maximal_orders(void) {
    /* Q(Sqrt[5]): the maximal order is Z[(1+Sqrt[5])/2] (the golden ratio),
     * NOT Z[Sqrt[5]] -- so the second basis element is ~1.618, not ~2.236. */
    check_prefix("N[NumberFieldIntegralBasis[Sqrt[5]]]", "{1.0, 1.61803}");
    check_zero("NumberFieldIntegralBasis[Sqrt[5]][[2]] - (1 + Sqrt[5])/2");

    /* Q(Sqrt[2]+Sqrt[3]): the equation order Z[theta] is not maximal; Round 2
     * enlarges it to O_K.  All four basis elements are algebraic integers. */
    check("AlgebraicIntegerQ /@ NumberFieldIntegralBasis[Sqrt[2] + Sqrt[3]]",
          "{True, True, True, True}");
    /* A complex-embedded cubic from a non-algebraic-integer generator
     * (Root[2+3#^3&,2] has minimal polynomial 3x^3+2): the generator is
     * rescaled to an algebraic integer before the basis is built. */
    check("Length[NumberFieldIntegralBasis[Root[2 + 3 #^3 &, 2]]]", "3");
    check("AlgebraicIntegerQ /@ NumberFieldIntegralBasis[Root[2 + 3 #^3 &, 2]]",
          "{True, True, True}");

    /* AlgebraicNumber input generating Q(10^(1/3)). */
    check("Length[NumberFieldIntegralBasis[AlgebraicNumber[Root[-10 + #^3 &, 1], {1, 0, 1/3}]]]", "3");
    check("AlgebraicIntegerQ /@ NumberFieldIntegralBasis[AlgebraicNumber[Root[-10 + #^3 &, 1], {1, 0, 1/3}]]",
          "{True, True, True}");

    /* Every integer combination of an integral basis is an algebraic integer. */
    check("AlgebraicIntegerQ[NumberFieldIntegralBasis[Sqrt[2 + Sqrt[3]]] . {-1, 2, 2, 3}]",
          "True");
    check("AlgebraicIntegerQ[NumberFieldIntegralBasis[Sqrt[2] + Sqrt[3]] . {5, -3, 2, 7}]",
          "True");
}

/* ---------------- Listable threading -------------------------------------- */

static void test_threading(void) {
    /* Listable: one basis per generator (each is its own field). */
    check("Length[NumberFieldIntegralBasis[{Sqrt[2], Sqrt[5]}]]", "2");
    check("Map[Length, NumberFieldIntegralBasis[{Sqrt[2], Sqrt[5]}]]", "{2, 2}");
    check_prefix("N[NumberFieldIntegralBasis[{Sqrt[2], Sqrt[5]}]]",
                 "{{1.0, 1.41421}, {1.0, 1.61803}}");
}

/* ---------------- ToNumberField primitive element ------------------------- */

static void test_primitive_element(void) {
    /* Q(Sqrt[2], I) via a primitive element (Sqrt[2] + I): degree-4 basis, all
     * algebraic integers. */
    check("Length[NumberFieldIntegralBasis[ToNumberField[{Sqrt[2], I}, All][[1, 1]]]]", "4");
    check("AlgebraicIntegerQ /@ NumberFieldIntegralBasis[ToNumberField[{Sqrt[2], I}, All][[1, 1]]]",
          "{True, True, True, True}");
}

/* ---------------- Declines ------------------------------------------------- */

static void test_declines(void) {
    /* Non-algebraic / symbolic arguments leave the head unevaluated (a message
     * is emitted to stderr). */
    check("NumberFieldIntegralBasis[x]", "NumberFieldIntegralBasis[x]");
    check("NumberFieldIntegralBasis[Pi]", "NumberFieldIntegralBasis[Pi]");
    /* Wrong arity stays unevaluated. */
    check("NumberFieldIntegralBasis[1, 2]", "NumberFieldIntegralBasis[1, 2]");
}

int main(void) {
    symtab_init();
    core_init();
    if (!flint_bridge_available()) {
        printf("FLINT not compiled in (USE_FLINT off); skipping "
               "NumberFieldIntegralBasis tests.\n");
        return 0;
    }
    test_algebraic_integer_q();
    test_shape();
    test_canonical_fields();
    test_maximal_orders();
    test_threading();
    test_primitive_element();
    test_declines();
    printf("test_numberfieldintegralbasis: all passed\n");
    return 0;
}
