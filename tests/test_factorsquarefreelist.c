/* Tests for FactorSquareFreeList[poly] and the Extension option on
 * FactorSquareFree.
 *
 * FactorSquareFreeList is a thin wrapper over FactorSquareFree (as FactorList
 * is over Factor): it returns the square-free factors with their multiplicities
 * as {factor, exponent} pairs, with a leading {c, 1} numerical factor.  The
 * Extension option (default Extension -> None, factoring over Q) is handled in
 * FactorSquareFree and forwarded verbatim; this file exercises both surfaces.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>

/* Shared driver: parse, evaluate, compare the printed form. */
static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* res_str = expr_to_string(res);
    if (strcmp(res_str, expected) != 0) {
        printf("FactorSquareFreeList test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- The spec example and basic square-free decompositions ---- */
void test_fsfl_basic() {
    /* The defining example: repeated linear factor (-1+x)^2 plus a square-free
     * cubic; note the cubic is NOT further factored (that is FactorList's job). */
    check("FactorSquareFreeList[x^5 - x^3 - x^2 + 1]",
          "{{1, 1}, {-1 + x, 2}, {1 + 2 x + 2 x^2 + x^3, 1}}");
    /* Already square-free (multiplicity 1), irreducible-or-not is irrelevant. */
    check("FactorSquareFreeList[x^2+1]", "{{1, 1}, {1 + x^2, 1}}");
    /* x^2-1 is square-free, so it is returned whole, NOT split into (x-1)(x+1). */
    check("FactorSquareFreeList[x^2-1]", "{{1, 1}, {-1 + x^2, 1}}");
    /* A different doubled factor. */
    check("FactorSquareFreeList[x^4-9x^3+29x^2-39x+18]",
          "{{1, 1}, {-3 + x, 2}, {2 - 3 x + x^2, 1}}");
    /* Pure square-free power: 1 - x^3 is square-free. */
    check("FactorSquareFreeList[(1-x^3)^2]", "{{1, 1}, {1 - x^3, 2}}");
}

/* ---- The leading numerical factor {c, 1} ---- */
void test_fsfl_numeric_factor() {
    /* Content 2 is pulled into the leading pair. */
    check("FactorSquareFreeList[2x^3+2x^2-2x-2]",
          "{{2, 1}, {-1 + x, 1}, {1 + x, 2}}");
}

/* ---- Listable: threads element-wise over a list argument ---- */
void test_fsfl_listable() {
    check("FactorSquareFreeList[{x^2-1, x^2+1}]",
          "{{{1, 1}, {-1 + x^2, 1}}, {{1, 1}, {1 + x^2, 1}}}");
}

/* ---- Round-trip invariant: Times @@ Power @@@ FactorSquareFreeList[p] == p ----
 * Order- and normalization-independent. */
void test_fsfl_roundtrip() {
    check("Expand[(Times@@Power@@@FactorSquareFreeList[x^5-x^3-x^2+1]) - (x^5-x^3-x^2+1)]", "0");
    check("Expand[(Times@@Power@@@FactorSquareFreeList[x^4-9x^3+29x^2-39x+18]) - "
          "(x^4-9x^3+29x^2-39x+18)]", "0");
    check("Expand[(Times@@Power@@@FactorSquareFreeList[2x^3+2x^2-2x-2]) - "
          "(2x^3+2x^2-2x-2)]", "0");
}

/* ---- Extension option on FactorSquareFree ---- */
void test_fsf_extension() {
    /* Extension -> None equals the plain form. */
    check("FactorSquareFree[x^5-x^3-x^2+1, Extension->None]",
          "(-1 + x)^2 (1 + 2 x + 2 x^2 + x^3)");
    /* Algebraic-coefficient input: x^2 - 2 Sqrt[2] x + 2 = (x - Sqrt[2])^2.
     * This genuinely needs the extension -- over Q it is not a rational poly. */
    check("FactorSquareFree[x^2 - 2 Sqrt[2] x + 2, Extension->Sqrt[2]]",
          "(-Sqrt[2] + x)^2");
    /* Square-free structure is field-independent for rational inputs, so an
     * extension leaves the answer unchanged. */
    check("FactorSquareFree[x^4-2, Extension->Sqrt[2]]", "-2 + x^4");
    check("Expand[FactorSquareFree[x^4-4x^2+4, Extension->Sqrt[2]] - "
          "FactorSquareFree[x^4-4x^2+4]]", "0");
    /* Real-coefficient input over Q(I): the (x-I)(x+I) split recombines on
     * expansion, giving the field-independent square-free part (1+x^2)^2. */
    check("FactorSquareFree[x^4+2x^2+1, Extension->I]", "(1 + x^2)^2");
}

/* ---- Extension option forwarded through FactorSquareFreeList ---- */
void test_fsfl_extension() {
    check("FactorSquareFreeList[x^2 - 2 Sqrt[2] x + 2, Extension->Sqrt[2]]",
          "{{1, 1}, {-Sqrt[2] + x, 2}}");
    check("FactorSquareFreeList[x^4+2x^2+1, Extension->I]",
          "{{1, 1}, {1 + x^2, 2}}");
    /* Round-trip over the extension. */
    check("Expand[(Times@@Power@@@FactorSquareFreeList[x^2 - 2 Sqrt[2] x + 2, "
          "Extension->Sqrt[2]]) - (x^2 - 2 Sqrt[2] x + 2)]", "0");
}

/* ---- Attributes ---- */
void test_fsfl_attributes() {
    check("Attributes[FactorSquareFreeList]", "{Listable, Protected}");
    check("MemberQ[Attributes[FactorSquareFreeList], Listable]", "True");
    check("MemberQ[Attributes[FactorSquareFreeList], Protected]", "True");
}

/* ---- Argument / option errors: left unevaluated ---- */
void test_fsfl_arity() {
    check("FactorSquareFreeList[]", "FactorSquareFreeList[]");
    /* A non-option trailing argument declines all the way up. */
    check("FactorSquareFreeList[x, y]", "FactorSquareFreeList[x, y]");
    check("FactorSquareFreeList[x^2-1, Foo->1]", "FactorSquareFreeList[-1 + x^2, Foo -> 1]");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_fsfl_basic);
    TEST(test_fsfl_numeric_factor);
    TEST(test_fsfl_listable);
    TEST(test_fsfl_roundtrip);
    TEST(test_fsf_extension);
    TEST(test_fsfl_extension);
    TEST(test_fsfl_attributes);
    TEST(test_fsfl_arity);

    printf("All FactorSquareFreeList tests passed!\n");
    return 0;
}
