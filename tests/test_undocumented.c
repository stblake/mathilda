/*
 * test_undocumented.c -- Indexed, $InitialDirectory / $RootDirectory and the
 * Internal` utilities (Internal`InheritedBlock, ExceptionFreeQ,
 * PartitionRagged, the integer / real-valued predicates).  GitHub issue #127.
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

static void run(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string(res);
    if (strcmp(s, expected) != 0)
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
    ASSERT_STR_EQ(s, expected);
    free(s);
    expr_free(e);
    expr_free(res);
}

static void test_indexed(void) {
    run("Indexed[{a, b, c}, 2]", "b");
    run("Indexed[{a, b, c}, -1]", "c");
    run("Indexed[{{1, 2}, {3, 4}}, {2, 1}]", "3");
    run("Indexed[x, 1]", "Indexed[x, 1]");
    run("Indexed[x, {1, 2}]", "Indexed[x, {1, 2}]");
    run("Indexed[{a, b}, 5]", "Indexed[{a, b}, 5]");     /* out of range: inert */
    run("Indexed[{{1, 2}}, {1, 3}]", "Indexed[{{1, 2}}, {1, 3}]");
    run("N[Indexed[x, 1/2]]", "Indexed[x, 1/2]");       /* NHoldRest */
    run("Attributes[Indexed]", "{NHoldRest, Protected}");
    run("Indexed[x, 1] + Indexed[x, 1]", "2 Indexed[x, 1]");
}

static void test_directories(void) {
    run("$RootDirectory", "/");
    run("$InitialDirectory === Directory[]", "True");
    run("StringQ[$InitialDirectory]", "True");
    run("MemberQ[Attributes[$InitialDirectory], Protected]", "True");
}

static void test_exceptionfreeq(void) {
    run("Internal`ExceptionFreeQ[{1, x, f[y], 1.5}]", "True");
    run("Internal`ExceptionFreeQ[{1, Infinity}]", "False");
    run("Internal`ExceptionFreeQ[-Infinity]", "False");
    run("Internal`ExceptionFreeQ[f[ComplexInfinity]]", "False");
    run("Internal`ExceptionFreeQ[Indeterminate]", "False");
    run("Internal`ExceptionFreeQ[Overflow[]]", "False");
}

static void test_predicates(void) {
    run("Internal`PositiveIntegerQ /@ {1, 0, -2, 10^30, 1.0, x}",
        "{True, False, False, True, False, False}");
    run("Internal`NonNegativeIntegerQ /@ {1, 0, -2, 10^30, 1.0, x}",
        "{True, True, False, True, False, False}");
    run("Internal`RealValuedNumberQ /@ {1, 1/2, 1.5, I, 1 + 2 I, Pi, x}",
        "{True, True, True, False, False, False, False}");
    run("Internal`RealValuedNumericQ /@ {1, Pi, Sqrt[2], Log[3], I, Sqrt[-2], x}",
        "{True, True, True, True, False, False, False}");
}

static void test_partitionragged(void) {
    run("Internal`PartitionRagged[{a, b, c, d, e, f}, {1, 2, 3}]", "{{a}, {b, c}, {d, e, f}}");
    run("Internal`PartitionRagged[{1, 2, 3, 4}, {0, 4, 0}]", "{{}, {1, 2, 3, 4}, {}}");
    run("Internal`PartitionRagged[{1, 2, 3, 4}, {1, 2}]",
        "Internal`PartitionRagged[{1, 2, 3, 4}, {1, 2}]");
    run("Internal`PartitionRagged[{1, 2}, {3}]", "Internal`PartitionRagged[{1, 2}, {3}]");
}

static void test_inheritedblock(void) {
    run("ibg[x_] := x^2; iby = 7; Null", "Null");
    /* The body starts from the current definitions and may change them ... */
    run("Internal`InheritedBlock[{ibg, iby}, ibg[x_] := x^3; iby = iby + 1; {ibg[2], iby}]",
        "{8, 8}");
    /* ... and every change is undone on exit. */
    run("{ibg[2], iby, Length[DownValues[ibg]]}", "{4, 7, 1}");
    /* Block, by contrast, starts from an unset symbol. */
    run("Block[{iby}, iby]", "iby");
    run("Internal`InheritedBlock[{iby}, iby]", "7");
    run("Internal`InheritedBlock[{iby = 2}, iby]", "2");
    run("iby", "7");
    /* Attributes are inherited and restored. */
    run("SetAttributes[ibh, Orderless]; ibh[x_, y_] := {x, y}; Null", "Null");
    run("Internal`InheritedBlock[{ibh}, ClearAttributes[ibh, Orderless]; ibh[b, a]]", "{b, a}");
    run("{ibh[b, a], Attributes[ibh]}", "{{a, b}, {Orderless}}");
    /* System symbols: temporary rules on a Protected builtin. */
    run("Internal`InheritedBlock[{Plus}, Unprotect[Plus]; Plus[ia, ib] = 0; ia + ib]", "0");
    run("{ia + ib, MemberQ[Attributes[Plus], Protected]}", "{ia + ib, True}");
    run("Internal`InheritedBlock[{iby}, Return[5]; 6]", "5");
    run("Attributes[Internal`InheritedBlock]", "{HoldAll, Protected}");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_indexed);
    TEST(test_directories);
    TEST(test_exceptionfreeq);
    TEST(test_predicates);
    TEST(test_partitionragged);
    TEST(test_inheritedblock);
    printf("All undocumented tests passed!\n");
    return 0;
}
