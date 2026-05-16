#include "patterns.h"
#include "print.h"
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>

void test_cases() {
    assert_eval_eq("Cases[{1, 1, f[a], 2, 3, y, f[8], 9, f[10]}, _Integer]", "{1, 1, 2, 3, 9}", 0);
    assert_eval_eq("Cases[{1, 1, f[a], 2, 3, y, f[8], 9, f[10]}, f[x_] -> x]", "{a, 8, 10}", 0);
    assert_eval_eq("Cases[{{1, 2}, {2}, {3, 4, 1}, {5, 4}, {3, 3}}, {_, _}]", "{{1, 2}, {5, 4}, {3, 3}}", 0);
    assert_eval_eq("Cases[{{1, 4, a, 0}, {b, 3, 2, 2}, {c, c, 5, 5}}, _Integer, 2]", "{1, 4, 0, 3, 2, 2, 5, 5}", 0);
    assert_eval_eq("Cases[{f[{a, b}], f[{a}], g[{a}], f[{a, b, c, d}]}, f[x_] :> Length[x]]", "{2, 1, 4}", 0);
}

void test_position() {
    assert_eval_eq("Position[{a,b,a,a,b,c,b},b]", "{{2}, {5}, {7}}", 0);
    assert_eval_eq("Position[{{a,a,b},{b,a,a},{a,b,a}},b]", "{{1, 3}, {2, 1}, {3, 2}}", 0);
    assert_eval_eq("Position[{1+x^2,5,x^4,a+(1+x^2)^2},x^_]", "{{1, 2}, {3}, {4, 2, 1, 2}}", 0);
    assert_eval_eq("Position[{1+x^2,5,x^4,a+(1+x^2)^2},x^_,2]", "{{1, 2}, {3}}", 0);
    assert_eval_eq("Position[{a,b,a,a,b,c,b,a,b},b,1,2]", "{{2}, {5}}", 0);
    assert_eval_eq("Position[f[g[h[x]]],_,Infinity]", "{{0}, {1, 0}, {1, 1, 0}, {1, 1, 1}, {1, 1}, {1}}", 0);
    assert_eval_eq("Position[f[a,b,b,a],b]", "{{2}, {3}}", 0);
    assert_eval_eq("Position[x^2+y^2,Power]", "{{1, 0}, {2, 0}}", 0);
    assert_eval_eq("Position[x^2+y^2,Power,Heads->False]", "{}", 0);
}

void test_count() {
    assert_eval_eq("Count[{a,b,a,a,b,c,b},b]", "3", 0);
    assert_eval_eq("Count[{a,2,a,a,1,c,b,3,3},_Integer]", "4", 0);
    assert_eval_eq("Count[{a,b,a,a,b,c,b,a,a},Except[b]]", "6", 0);
    assert_eval_eq("Count[{{a,a,b},b,{a,b,a}},b,2]", "3", 0);
    assert_eval_eq("Count[{{a,a,b},b,{a,b,a}},b,{2}]", "2", 0);
    assert_eval_eq("Count[x^3+1.5x^2+Pi x +7,_?NumberQ,-1]", "4", 0);
    assert_eval_eq("Count[5,_?NumberQ,-1]", "1", 0);
    assert_eval_eq("Count[5,_?NumberQ,{0,-1}]", "0", 0);
}

void test_delete_cases() {
    /* Documented examples from the spec. */
    assert_eval_eq("DeleteCases[{1, 1, x, 2, 3, y, 9, y}, _Integer]", "{x, y, y}", 0);
    assert_eval_eq("DeleteCases[_Integer][{1, 1, x, 2, 3, y, 9, y}]", "{x, y, y}", 0);
    assert_eval_eq("DeleteCases[{1, f[2, 3], 4}, f, {2}, Heads -> True]", "{1, 2, 3, 4}", 0);

    /* Default levelspec is {1}: only level-1 elements are tested. */
    assert_eval_eq("DeleteCases[{a, b, b, a, b, c, b}, b]", "{a, a, c}", 0);
    assert_eval_eq("DeleteCases[{1, 2, 3, x, y, z}, _Symbol]", "{1, 2, 3}", 0);
    assert_eval_eq("DeleteCases[{1, 2, 3, x, y, z}, _Integer]", "{x, y, z}", 0);

    /* Nothing matches -- unchanged copy. */
    assert_eval_eq("DeleteCases[{1, 2, 3}, _Real]", "{1, 2, 3}", 0);
    /* Empty list is preserved. */
    assert_eval_eq("DeleteCases[{}, _Integer]", "{}", 0);

    /* Pattern matches against compound elements at level 1. */
    assert_eval_eq("DeleteCases[{f[1], g[2], f[3], h[4]}, f[_]]", "{g[2], h[4]}", 0);
    assert_eval_eq("DeleteCases[{{1, 2}, {2, 3}, {3, 4}}, {_, 3}]", "{{1, 2}, {3, 4}}", 0);
    assert_eval_eq("DeleteCases[{{1, 2}, {2}, {3, 4, 1}, {5, 4}, {3, 3}}, {_, _}]", "{{2}, {3, 4, 1}}", 0);

    /* DeleteCases works on arbitrary heads, not just List. */
    assert_eval_eq("DeleteCases[a + b + c, b]", "a + c", 0);
    assert_eval_eq("DeleteCases[f[1, 2, 3, 4], 2]", "f[1, 3, 4]", 0);

    /* levelspec: positive integer N = levels {1..N}. */
    assert_eval_eq("DeleteCases[{{1, 4, a, 0}, {b, 3, 2, 2}, {c, c, 5, 5}}, _Integer, 2]",
                   "{{a}, {b}, {c, c}}", 0);
    /* levelspec: explicit {n}. */
    assert_eval_eq("DeleteCases[{{a, a, b}, b, {a, b, a}}, b, {2}]", "{{a, a}, b, {a, a}}", 0);
    /* levelspec: {2} -- only level 2, top-level b is preserved. */
    assert_eval_eq("DeleteCases[{a, b, {a, b, c}}, b, {2}]", "{a, b, {a, c}}", 0);
    /* levelspec: Infinity / All -- all levels from 1 down. */
    assert_eval_eq("DeleteCases[{a, b, {a, b, {a, b}}}, b, Infinity]",
                   "{a, {a, {a}}}", 0);

    /* Count limit n. With levelspec {1}, only the first n level-1 matches are removed. */
    assert_eval_eq("DeleteCases[{1, 2, 3, 4, 5}, _Integer, {1}, 2]", "{3, 4, 5}", 0);
    assert_eval_eq("DeleteCases[{1, 2, 3, 4, 5}, _Integer, {1}, 0]", "{1, 2, 3, 4, 5}", 0);
    /* Limit larger than the number of matches: deletes all matches, no error. */
    assert_eval_eq("DeleteCases[{1, 2, x}, _Integer, {1}, 99]", "{x}", 0);

    /* Heads -> True deletes heads, behaving like FlattenAt. */
    assert_eval_eq("DeleteCases[{1, f[2, 3], 4}, f, {2}, Heads -> True]", "{1, 2, 3, 4}", 0);
    /* Heads -> False (default) leaves head matches alone. */
    assert_eval_eq("DeleteCases[{1, f[2, 3], 4}, f, {2}, Heads -> False]", "{1, f[2, 3], 4}", 0);
    assert_eval_eq("DeleteCases[{1, f[2, 3], 4}, f, {2}]", "{1, f[2, 3], 4}", 0);

    /* Test pattern semantics: matches based on the ORIGINAL expression
     * (Mathematica behaviour). An outer node that matches the pattern is
     * deleted even after its children have already been transformed. */
    assert_eval_eq("DeleteCases[{f[1, 2], f[3, f[4]]}, f[__], Infinity]", "{}", 0);
    assert_eval_eq("DeleteCases[{f[1, 2], f[3, f[4]]}, f[__], Infinity, 1]",
                   "{f[3, f[4]]}", 0);

    /* Operator form preserved across applications. */
    assert_eval_eq("DeleteCases[_Symbol][{1, x, 2, y, 3}]", "{1, 2, 3}", 0);

    /* PatternTest and Condition guards. */
    assert_eval_eq("DeleteCases[{1, 2, 3, 4, 5, 6}, _?EvenQ]", "{1, 3, 5}", 0);
    assert_eval_eq("DeleteCases[{1, 2, 3, 4, 5}, x_ /; x > 2]", "{1, 2}", 0);

    /* Negative levelspec: level -1 = atomic leaves. */
    assert_eval_eq("DeleteCases[{1, {2, 3}, 4}, _Integer, {-1}]", "{{}}", 0);
}

void test_memberq() {
    assert_eval_eq("MemberQ[{1,3,4,1,2},2]", "True", 0);
    assert_eval_eq("MemberQ[{1,3,4,1,5},2]", "False", 0);
    assert_eval_eq("MemberQ[{x^2,y^2,x^3},x^_]", "True", 0);
    assert_eval_eq("MemberQ[{{1,1,3,0},{2,1,2,2}},0,2]", "True", 0);
    assert_eval_eq("MemberQ[{{1,1,3,0},{2,1,2,2}},0]", "False", 0);
    assert_eval_eq("MemberQ[Table[Mod[2^i,7],{i,10}],1]", "True", 0);
    assert_eval_eq("MemberQ[f[a,b,c],b]", "True", 0);
    assert_eval_eq("MemberQ[f[a,b,c],d]", "False", 0);
    assert_eval_eq("MemberQ[1]", "Function[MemberQ[Slot[1], 1]]", 1);
}

int main() {
    symtab_init();
    core_init();
    
    TEST(test_cases);
    TEST(test_delete_cases);
    TEST(test_position);
    TEST(test_count);
    TEST(test_memberq);
    
    printf("All patterns tests passed!\n");
    symtab_clear();
    return 0;
}