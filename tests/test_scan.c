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

/* Evaluate `input` and compare its infix printed form to `expected`.
 *
 * Scan itself returns Null and discards f's results, so most tests observe a
 * side effect: they thread a list (or counter) through an assignment inside a
 * CompoundExpression and read it back. To pin down the *order* and *set* of the
 * parts a level spec selects, the side effect accumulates Depth[#] of each
 * visited part: for nested single-element lists the depth uniquely identifies
 * each part, so the collected integer sequence encodes visitation order (and
 * post-order = deepest-in-range first). */
static void run(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    ASSERT_MSG(strcmp(s, expected) == 0, "Scan %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* ---------- Default level {1} ---------- */

static void test_scan_default_level1(void) {
    run("r={};Scan[(AppendTo[r,#])&,{a,b,c}];r", "{a, b, c}");
    /* Only the level-1 element (the single inner list {{{{a}}}}) is visited; the
     * collected list wraps it, adding one brace layer. */
    run("r={};Scan[(AppendTo[r,#])&,{{{{{a}}}}}];r", "{{{{{a}}}}}");
}

static void test_scan_returns_null(void) {
    run("Scan[gg,{a,b,c}]", "Null");
    run("Scan[Print,{a,b}]", "Null");
}

/* ---------- Positive depth / single level / ranges ----------
 * Depths of {{{{{a}}}}}: L1=5, L2=4, L3=3, L4=2, L5=1 (a). */

static void test_scan_positive_depth(void) {
    /* level 2 -> {1,2}: post-order visits L2 (depth 4) then L1 (depth 5). */
    run("r={};Scan[(AppendTo[r,Depth[#]])&,{{{{{a}}}}},2];r", "{4, 5}");
    /* level 3 -> {1,3}: L3(3), L2(4), L1(5). */
    run("r={};Scan[(AppendTo[r,Depth[#]])&,{{{{{a}}}}},3];r", "{3, 4, 5}");
}

static void test_scan_single_level(void) {
    /* {2} -> level 2 only. */
    run("r={};Scan[(AppendTo[r,Depth[#]])&,{{{{{a}}}}},{2}];r", "{4}");
}

static void test_scan_range_with_zero(void) {
    /* {0,2} -> levels 0,1,2: L2(4), L1(5), L0(6). */
    run("r={};Scan[(AppendTo[r,Depth[#]])&,{{{{{a}}}}},{0,2}];r", "{4, 5, 6}");
}

/* ---------- Infinity ----------
 * Depths of {{{a}}}: L1=3, L2=2, L3=1 (a); whole = 4. */

static void test_scan_infinity(void) {
    /* {1,Infinity}: a(1), {a}(2), {{a}}(3). */
    run("r={};Scan[(AppendTo[r,Depth[#]])&,{{{a}}},Infinity];r", "{1, 2, 3}");
    /* {0,Infinity} also includes the whole expression (depth 4). */
    run("r={};Scan[(AppendTo[r,Depth[#]])&,{{{a}}},{0,Infinity}];r", "{1, 2, 3, 4}");
}

/* ---------- Negative levels ----------
 * Depths of {{{{a}}}}: L1=4, L2=3, L3=2, L4=1 (a); whole = 5. */

static void test_scan_negative_all(void) {
    /* -1 -> {1,-1}: every part. */
    run("r={};Scan[(AppendTo[r,Depth[#]])&,{{{{a}}}},-1];r", "{1, 2, 3, 4}");
}

static void test_scan_negative_2(void) {
    /* -2 -> {1,-2}: excludes the depth-1 leaf a. */
    run("r={};Scan[(AppendTo[r,Depth[#]])&,{{{{a}}}},-2];r", "{2, 3, 4}");
}

static void test_scan_negative_single(void) {
    /* {-2} -> depth 2 only (the part {a}). */
    run("r={};Scan[(AppendTo[r,Depth[#]])&,{{{{a}}}},{-2}];r", "{2}");
}

/* ---------- Mixed positive/negative ----------
 * {2,-3}: level >= 2 AND depth >= 3. */

static void test_scan_mixed_levels_list(void) {
    /* {{{{{a}}}}}: selects L2(depth4) and L3(depth3); post-order -> {3,4}. */
    run("r={};Scan[(AppendTo[r,Depth[#]])&,{{{{{a}}}}},{2,-3}];r", "{3, 4}");
}

static void test_scan_mixed_levels_heads(void) {
    /* Non-List heads: h0[h1[h2[h3[h4[a]]]]] -> h3[h4[a]], h2[h3[h4[a]]]. */
    run("r={};Scan[(AppendTo[r,#])&,h0[h1[h2[h3[h4[a]]]]],{2,-3}];r",
        "{h3[h4[a]], h2[h3[h4[a]]]}");
}

/* ---------- Non-List heads ---------- */

static void test_scan_plus_head(void) {
    run("r={};Scan[(AppendTo[r,#])&,a+b+c];r", "{a, b, c}");
}

static void test_scan_power_level2(void) {
    /* x^2 + y^2 at level 2 visits the leaves of each Power: x,2,y,2. */
    run("r={};Scan[(AppendTo[r,#])&,x^2+y^2,{2}];r", "{x, 2, y, 2}");
}

/* ---------- Heads option ---------- */

static void test_scan_heads_true(void) {
    /* Heads->True visits the head (List) before the elements. */
    run("r={};Scan[(AppendTo[r,#])&,{a,b},Heads->True];r", "{List, a, b}");
}

static void test_scan_heads_default_false(void) {
    run("r={};Scan[(AppendTo[r,#])&,{a,b}];r", "{a, b}");
}

/* ---------- Throw / Catch ---------- */

static void test_scan_throw(void) {
    run("Catch[Scan[If[#>5,Throw[#]]&,{2,4,6,8}]]", "6");
    run("Catch[Scan[If[#>2,Throw[#]]&,{1,2,3,4}]]", "3");
}

/* ---------- Return ---------- */

static void test_scan_return_boundary(void) {
    /* Return[ret, Scan] targets the Scan boundary explicitly. */
    run("Scan[If[#>2,Return[#,Scan]]&,{1,2,3,4}]", "3");
}

static void test_scan_return_symbol_f(void) {
    /* A symbol f whose DownValue yields a raw Return[...] (no Function boundary
     * to consume it) exits Scan with that value. */
    run("Clear[fret];fret[x_]:=If[x>2,Return[x]];Scan[fret,{1,2,3,4}]", "3");
}

/* ---------- Association (values scanned) ---------- */

static void test_scan_association_sum(void) {
    run("sc=0;Scan[(sc=sc+#)&,<|\"a\"->1,\"b\"->2,\"c\"->3|>];sc", "6");
}

static void test_scan_association_null(void) {
    run("Scan[Identity,<|\"a\"->1|>]", "Null");
}

/* ---------- NDArray fast path ---------- */

static void test_scan_ndarray_rank1(void) {
    /* Rank-1 default level: elements box as Real; sum is 10. */
    run("s=0;Scan[(s=s+#)&,NDArray[Range[4]]];s==10", "True");
}

static void test_scan_ndarray_rank2_rows(void) {
    /* Rank-2 default level: each level-1 part is a sub-NDArray row. */
    run("r={};Scan[(AppendTo[r,Head[#]])&,NDArray[Table[i+j,{i,2},{j,3}]]];r",
        "{NDArray, NDArray}");
}

static void test_scan_ndarray_throw(void) {
    run("Catch[Scan[If[#>3.5,Throw[#]]&,NDArray[Range[6]]]]==4", "True");
}

static void test_scan_ndarray_deep_level_fallback(void) {
    /* A non-default level spec over an NDArray materializes to a nested list and
     * scans generically, matching the plain-List result. Round removes the
     * Real-vs-Integer formatting difference. */
    run("t1=0;Scan[(t1=t1+#)&,NDArray[Table[i+j,{i,2},{j,3}]],Infinity];"
        "Round[t1]=={26,28,30}", "True");
}

/* ---------- Atoms / arity ---------- */

static void test_scan_atom(void) {
    run("Scan[Print,x]", "Null");
}

static void test_scan_arity(void) {
    /* Scan[] leaves itself unevaluated (and prints Scan::argb to stderr). */
    run("Scan[]", "Scan[]");
    /* Scan[f] (single arg) is left unevaluated (no operator form). */
    run("Scan[gg]", "Scan[gg]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_scan_default_level1);
    TEST(test_scan_returns_null);
    TEST(test_scan_positive_depth);
    TEST(test_scan_single_level);
    TEST(test_scan_range_with_zero);
    TEST(test_scan_infinity);
    TEST(test_scan_negative_all);
    TEST(test_scan_negative_2);
    TEST(test_scan_negative_single);
    TEST(test_scan_mixed_levels_list);
    TEST(test_scan_mixed_levels_heads);
    TEST(test_scan_plus_head);
    TEST(test_scan_power_level2);
    TEST(test_scan_heads_true);
    TEST(test_scan_heads_default_false);
    TEST(test_scan_throw);
    TEST(test_scan_return_boundary);
    TEST(test_scan_return_symbol_f);
    TEST(test_scan_association_sum);
    TEST(test_scan_association_null);
    TEST(test_scan_ndarray_rank1);
    TEST(test_scan_ndarray_rank2_rows);
    TEST(test_scan_ndarray_throw);
    TEST(test_scan_ndarray_deep_level_fallback);
    TEST(test_scan_atom);
    TEST(test_scan_arity);

    printf("All Scan tests passed!\n");
    return 0;
}
