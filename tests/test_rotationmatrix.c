/* Tests for RotationMatrix (2-D and 3-D).  Expected forms follow the
 * reference where it gives a canonical one ({{Cos[t], -Sin[t]}, ...}, the
 * rational RotationMatrix[{{1,2,3},{3,2,1}}], the messages); the exact radical
 * forms are checked structurally (orthogonal, det 1, maps the vector). */

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
        printf("RotationMatrix test failed: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
        ASSERT(0);
    }
    free(s);
    expr_free(e);
    expr_free(r);
}

static void test_2d(void) {
    check("RotationMatrix[t]", "{{Cos[t], -Sin[t]}, {Sin[t], Cos[t]}}");
    check("RotationMatrix[Pi/2]", "{{0, -1}, {1, 0}}");
    check("RotationMatrix[Pi/3] . {1, 0}", "{1/2, 1/2 Sqrt[3]}");
    check("RotationMatrix[0.5] == {{Cos[0.5], -Sin[0.5]}, {Sin[0.5], Cos[0.5]}}", "True");
    check("RotationMatrix[{{1, 0}, {0, 1}}]", "{{0, -1}, {1, 0}}");
    check("RotationMatrix[{{1, 1}, {0, 1}}] . {1, 1}", "{0, Sqrt[2]}");
    check("RotationMatrix[{{1., 0.}, {-1., 0.}}]", "{{-1.0, 0.0}, {0.0, -1.0}}");
    check("RotationMatrix[t, {{0, 1}, {1, 0}}] == RotationMatrix[-t]", "True");
}

static void test_3d(void) {
    check("RotationMatrix[t, {0, 0, 1}]", "{{Cos[t], -Sin[t], 0}, {Sin[t], Cos[t], 0}, {0, 0, 1}}");
    check("RotationMatrix[Pi/2, {0, 0, 2}]", "{{0, -1, 0}, {1, 0, 0}, {0, 0, 1}}");
    check("RotationMatrix[{{1, 0, 0}, {0, 0, 1}}]", "{{0, 0, -1}, {0, 1, 0}, {1, 0, 0}}");
    check("RotationMatrix[{{1, 0, 0}, {0, 1, 0}}]", "{{0, -1, 0}, {1, 0, 0}, {0, 0, 1}}");
    /* the unnormalised Rodrigues form keeps this all-rational */
    check("RotationMatrix[{{1, 2, 3}, {3, 2, 1}}]",
          "{{16/21, 4/21, 13/21}, {-8/21, 19/21, 4/21}, {-11/21, -8/21, 16/21}}");
    check("RotationMatrix[{{1, 2, 3}, {3, 2, 1}}] . {1, 2, 3}", "{3, 2, 1}");
    check("RotationMatrix[{{1, 0, 0}, {2, 0, 0}}]", "{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}");
    check("m = RotationMatrix[Pi/2, {1, 1, 1}]; Simplify[{m . Transpose[m], Det[m], m . {1, 1, 1}}]",
          "{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, 1, {1, 1, 1}}");
    check("Simplify[RotationMatrix[t, {a, b, c}] . {a, b, c}]", "{a, b, c}");
    check("Simplify[Det[RotationMatrix[t, {1, 2, 2}]]]", "1");
    check("Simplify[RotationMatrix[t, {{1, 0, 0}, {0, 1, 0}}] - RotationMatrix[t, {0, 0, 1}]]",
          "{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}");
}

static void test_messages(void) {
    check("Check[RotationMatrix[{1, 0}], fired]", "fired");
    check("Quiet[RotationMatrix[{1, 0}]]", "RotationMatrix[{1, 0}]");
    check("Check[RotationMatrix[t, {1, 0}], fired]", "fired");
    check("Check[RotationMatrix[{{1, 0, 0}, {-1, 0, 0}}], fired]", "fired");
    check("Quiet[RotationMatrix[{{1, 0, 0}, {-1, 0, 0}}]]", "RotationMatrix[{{1, 0, 0}, {-1, 0, 0}}]");
}

static void test_attributes(void) {
    check("Attributes[RotationMatrix]", "{Protected}");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_2d);
    TEST(test_3d);
    TEST(test_messages);
    TEST(test_attributes);
    printf("All RotationMatrix tests passed!\n");
    return 0;
}
