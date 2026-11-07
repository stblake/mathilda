/* Tests for CubeRoot[x] -- the real cube root, Surd[x, 3] under its own name.
 * Semantics follow the reference: exact roots, -2^(1/3) for CubeRoot[-2],
 * Surd[x, 3] for symbolic x, CubeRoot::preal (unevaluated) for a non-real
 * number; machine values from the shared sf_surd_r kernel on every surface. */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0) {
        printf("CubeRoot test failed: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
        ASSERT(0);
    }
    free(s);
    expr_free(e);
    expr_free(r);
}

static void test_exact(void) {
    check("CubeRoot[8]", "2");
    check("CubeRoot[-8]", "-2");
    check("CubeRoot[-27/8]", "-3/2");
    check("CubeRoot[2]", "2^(1/3)");
    check("CubeRoot[-2]", "-2^(1/3)");
    check("CubeRoot[16]", "2 2^(1/3)");
    check("CubeRoot[0]", "0");
    check("CubeRoot[x]", "Surd[x, 3]");
    check("CubeRoot[x^3]", "Surd[x^3, 3]");
}

static void test_numeric(void) {
    check("CubeRoot[8.]", "2.0");
    check("CubeRoot[-2.] == -2.^(1/3)", "True");
    check("CubeRoot[-2.] == Surd[-2., 3]", "True");
    check("Abs[CubeRoot[N[-2, 30]] + 1.25992104989487316476721060727822835057`35] < 10^-29", "True");
    check("Precision[CubeRoot[N[-2, 30]]]", "30.0");
}

static void test_nonreal(void) {
    check("Quiet[CubeRoot[I]]", "CubeRoot[I]");
    check("Check[CubeRoot[1. + I], fired]", "fired");
}

static void test_surfaces(void) {
    check("ys = {-8., 27., 2., -3.5, 0.001, 1.*^300}; cg = Compile[{{v, _Real, 1}}, CubeRoot[v]]; "
          "CubeRoot /@ ys == CubeRoot[ys] == CubeRoot[Developer`ToPackedArray[ys]] == "
          "Normal[CubeRoot[NDArray[ys]]] == cg[ys] == (Compile[{{x, _Real}}, CubeRoot[x]] /@ ys)", "True");
    check("Developer`PackedArrayQ[CubeRoot[Developer`ToPackedArray[{-8., 27.}]]]", "True");
    check("CompileDiagnostics[{{x, _Real}}, CubeRoot[x]][[1]]", "\"Compiled\" -> True");
    check("CompileDiagnostics[{{v, _Real, 1}}, CubeRoot[v]][[1]]", "\"Compiled\" -> True");
}

static void test_attributes(void) {
    check("Attributes[CubeRoot]", "{Listable, NumericFunction, Protected}");
    check("CubeRoot[{8, -1, 27}]", "{2, -1, 3}");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_exact);
    TEST(test_numeric);
    TEST(test_nonreal);
    TEST(test_surfaces);
    TEST(test_attributes);
    printf("All CubeRoot tests passed!\n");
    return 0;
}
