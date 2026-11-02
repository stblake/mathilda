/* test_numeric_literals_log.c
 *
 * Machine-real literals, InputForm of machine reals, Log[b, x] at machine
 * precision, ComplexExpand's grouping, and HeavisideTheta/DiracDelta
 * attributes (FeynCalc's FCCompareNumbers and FCLoopFindSectors, 7 tests).
 * Expected values computed by the reference kernel.
 *
 *  1. A literal of up to 17 significant digits is a machine real; 18 or more
 *     make it arbitrary precision (the cut was 15.95 implied digits).
 *  2. InputForm of a machine real is its shortest round-trip text, with no
 *     precision mark.
 *  3. Log[b, x] for a machine x and a numeric base is log(x)/log(b) in
 *     doubles, bit for bit as the reference.
 *  4. ComplexExpand returns re + I im with I outside the imaginary part.
 *  5. HeavisideTheta and DiracDelta are Listable, Orderless, Protected.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>

static void test_literal_precision(void) {
    assert_eval_eq("Precision[1.2887955878011237]", "MachinePrecision", 0);
    assert_eval_eq("Precision[9.9999999999999999]", "MachinePrecision", 0);
    assert_eval_eq("Precision[0.12345678901234567]", "MachinePrecision", 0);
    assert_eval_eq("MachineNumberQ[1.00000000000000001]", "False", 0);
    assert_eval_eq("ToString[-0.30102999566398114, InputForm]", "\"-0.30102999566398114\"", 0);
}

static void test_log_base_machine(void) {
    assert_eval_eq("ToString[Log[10, 0.008196721311475417], InputForm]", "\"-2.086359830674748\"", 0);
    assert_eval_eq("ToString[Log[10, 3.513194478131193*^-9], InputForm]", "\"-8.454297808522389\"", 0);
    assert_eval_eq("ToString[Log[2, 0.008196721311475417], InputForm]", "\"-6.930737337562886\"", 0);
    assert_eval_eq("ToString[Log[3, 0.008196721311475417], InputForm]", "\"-4.37280840045705\"", 0);
}

static void test_complexexpand_grouping(void) {
    assert_eval_eq("ComplexExpand[(1 + 2 I) x + (3 + 4 I) y]", "x + 3 y + I (2 x + 4 y)", 0);
    assert_eval_eq("ComplexExpand[(1 + I) x]", "(1 + I) x", 0);
    assert_eval_eq("ComplexExpand[Exp[I x]]", "Cos[x] + I Sin[x]", 0);
}

static void test_heaviside_attributes(void) {
    /* the rules and attributes live in src/internal/distributions.m */
    /* searched upward, so the suite runs from tests/ and from CI's tests/build/ */
    assert_eval_eq("Get[SelectFirst[{\"src\", \"../src\", \"../../src\", \"../../../src\"}, "
                   "FileExistsQ[# <> \"/internal/distributions.m\"] &] <> "
                   "\"/internal/distributions.m\"]; 0", "0", 0);
    assert_eval_eq("Attributes[HeavisideTheta]", "{Listable, Orderless, Protected}", 0);
    assert_eval_eq("Attributes[DiracDelta]", "{Listable, Orderless, Protected}", 0);
    assert_eval_eq("HeavisideTheta[{1, 0, 1} - 1/2]", "{1, 0, 1}", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_literal_precision);
    TEST(test_log_base_machine);
    TEST(test_complexexpand_grouping);
    TEST(test_heaviside_attributes);
    printf("All numeric literal / Log / ComplexExpand tests passed.\n");
    return 0;
}
