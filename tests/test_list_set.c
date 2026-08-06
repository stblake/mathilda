/*
 * test_list_set.c -- unit tests for list-destructuring assignment
 * ({a,b,...} = {...}), with a particular focus on the regression where
 * reassigning to symbols that already held OwnValues silently no-op'd
 * because Set evaluated the LHS elements (binding targets) to their
 * current values before destructuring.
 *
 * Reported symptom:
 *     {a,b,c,d} = {1,1,1,1}; {a++, ++b, c--, --d}     (* -> {1,2,1,0} *)
 *     {a,b,c,d} = {1,1,1,1}; {a++, ++b, c--, --d}     (* was {2,3,0,-1} *)
 * The second line returned the same mutation applied to the stale post-
 * increment state because the reassignment never took effect.
 */

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

/* Evaluate an input string and discard the result. Handy for setup lines. */
static void eval_and_discard(const char* src) {
    Expr* parsed = parse_expression(src);
    assert(parsed != NULL);
    Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    expr_free(evaluated);
}

/* Evaluate input and return its printed form; caller must free. */
static char* eval_to_string(const char* src) {
    Expr* parsed = parse_expression(src);
    assert(parsed != NULL);
    Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    char* s = expr_to_string(evaluated);
    expr_free(evaluated);
    return s;
}

static void clear_symbols(const char* const* names, size_t n) {
    for (size_t i = 0; i < n; i++) symtab_clear_symbol(names[i]);
}

/* --- Tests ---------------------------------------------------------- */

/* Baseline: first-time destructuring still works. */
static void test_destructure_fresh(void) {
    const char* syms[] = {"a", "b", "c", "d"};
    clear_symbols(syms, 4);

    eval_and_discard("{a, b, c, d} = {1, 1, 1, 1}");
    assert_eval_eq("{a, b, c, d}", "{1, 1, 1, 1}", 0);
}

/* Reassigning symbols that already have OwnValues must actually update
 * them. Before the fix, Set would evaluate a,b,c,d on the LHS to their
 * existing values and apply_assignment would silently no-op. */
static void test_reassignment_updates_targets(void) {
    const char* syms[] = {"a", "b", "c", "d"};
    clear_symbols(syms, 4);

    eval_and_discard("{a, b, c, d} = {1, 2, 3, 4}");
    eval_and_discard("{a, b, c, d} = {10, 20, 30, 40}");
    assert_eval_eq("{a, b, c, d}", "{10, 20, 30, 40}", 0);
}

/* The exact scenario from the user's bug report. */
static void test_reassign_then_increment_decrement(void) {
    const char* syms[] = {"a", "b", "c", "d"};
    clear_symbols(syms, 4);

    /* First round: fresh state, should produce {1, 2, 1, 0}. */
    eval_and_discard("{a, b, c, d} = {1, 1, 1, 1}");
    assert_eval_eq("{a++, ++b, c--, --d}", "{1, 2, 1, 0}", 0);

    /* After first round, values are {2, 2, 0, 0}. Reassign and repeat -- the
     * reported bug returned {2, 3, 0, -1} because the reassignment no-op'd. */
    eval_and_discard("{a, b, c, d} = {1, 1, 1, 1}");
    assert_eval_eq("{a++, ++b, c--, --d}", "{1, 2, 1, 0}", 0);
    assert_eval_eq("{a, b, c, d}", "{2, 2, 0, 0}", 0);
}

/* Simultaneous swap: {x, y} = {y, x} requires that the RHS is evaluated
 * in the old environment, and the LHS bindings take effect. */
static void test_simultaneous_swap(void) {
    const char* syms[] = {"x", "y"};
    clear_symbols(syms, 2);

    eval_and_discard("x = 1");
    eval_and_discard("y = 2");
    eval_and_discard("{x, y} = {y, x}");
    assert_eval_eq("{x, y}", "{2, 1}", 0);

    /* And swapping back must also work. */
    eval_and_discard("{x, y} = {y, x}");
    assert_eval_eq("{x, y}", "{1, 2}", 0);
}

/* Nested destructuring with reassignment: inner Lists are nested binding
 * groups, not expressions to evaluate. */
static void test_nested_destructuring(void) {
    const char* syms[] = {"a", "b", "c"};
    clear_symbols(syms, 3);

    eval_and_discard("{{a, b}, c} = {{10, 20}, 30}");
    assert_eval_eq("{a, b, c}", "{10, 20, 30}", 0);

    /* Reassign nested targets. */
    eval_and_discard("{{a, b}, c} = {{100, 200}, 300}");
    assert_eval_eq("{a, b, c}", "{100, 200, 300}", 0);
}

/* Non-symbol elements in a List LHS are function-shaped targets (DownValue
 * creation). Their inner arguments must still be evaluated so that, e.g.,
 * {a[x], b[y]} = {p, q} (with x=5, y=7) produces a[5]=p and b[7]=q. */
static void test_downvalue_destructuring(void) {
    const char* syms[] = {"a", "b", "x", "y"};
    clear_symbols(syms, 4);

    eval_and_discard("x = 5");
    eval_and_discard("y = 7");
    eval_and_discard("{a[x], b[y]} = {100, 200}");

    /* Values of x and y untouched. */
    assert_eval_eq("x", "5", 0);
    assert_eval_eq("y", "7", 0);

    /* DownValues created against the *evaluated* x and y. */
    assert_eval_eq("a[5]", "100", 0);
    assert_eval_eq("b[7]", "200", 0);
}

/* Pattern-bearing destructuring: {a[p_], b[q_]} = {1, 2} should create
 * generic DownValues (any a[anything] -> 1, any b[anything] -> 2). */
static void test_pattern_downvalue_destructuring(void) {
    const char* syms[] = {"a", "b"};
    clear_symbols(syms, 2);

    eval_and_discard("{a[p_], b[q_]} = {1, 2}");
    assert_eval_eq("a[42]", "1", 0);
    assert_eval_eq("b[99]", "2", 0);
}

/* Length mismatch leaves targets untouched. */
static void test_length_mismatch(void) {
    const char* syms[] = {"a", "b"};
    clear_symbols(syms, 2);

    eval_and_discard("a = 7");
    eval_and_discard("b = 8");

    /* Attempt a mismatched destructuring -- should NOT alter a or b. */
    Expr* parsed = parse_expression("{a, b} = {1, 2, 3}");
    Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    expr_free(evaluated);

    assert_eval_eq("{a, b}", "{7, 8}", 0);
}

/* Literal-integer LHS element: {1, a} = {1, 2} -- the prior code would
 * silently return {1, 2} and leave a unbound. Post-fix, the destructuring
 * reports failure, so the Set expression is returned unevaluated (or, at
 * minimum, 'a' does NOT get set to 2 via the literal path).
 *
 * We assert only the observable invariant: that `a` was not assigned to.
 */
static void test_literal_lhs_element_does_not_pretend_success(void) {
    const char* syms[] = {"a"};
    clear_symbols(syms, 1);

    Expr* parsed = parse_expression("{1, a} = {1, 2}");
    Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    expr_free(evaluated);

    /* Because the literal 1-vs-1 child recursion fails, the whole
     * destructuring is now treated as a failure and 'a' is untouched.
     * (Mathematica's own behavior here is to emit an error and leave a
     * unbound; we accept "a stays symbolic" as the correct outcome.) */
    char* s = eval_to_string("a");
    ASSERT_STR_EQ(s, "a");
    free(s);
}

/* Destructuring must not leak or corrupt when assigning in a loop (smoke
 * test for memory behavior under repeated reassignment). */
static void test_repeated_reassignment_smoke(void) {
    const char* syms[] = {"a", "b", "c"};
    clear_symbols(syms, 3);

    for (int i = 0; i < 50; i++) {
        char buf[64];
        snprintf(buf, sizeof(buf), "{a, b, c} = {%d, %d, %d}", i, i + 1, i + 2);
        eval_and_discard(buf);
    }
    assert_eval_eq("{a, b, c}", "{49, 50, 51}", 0);
}

/* ------------------------------------------------------------------------
 * Packed-array RHS destructuring (regression for the automatic-packing bug).
 *
 * `Set` is a packed-aware head so that a whole-value binding (x = Range[10^6])
 * keeps its argument packed rather than materialising it. The transparency gate
 * therefore leaves a packed argument intact for Set -- correct for x = big, but
 * destructuring is the one assignment path that reads the RHS *structure*. Before
 * the fix, {xc, yc} = {Range[m], Range[n]} evaluated the RHS to a packed rank-2
 * array (a single EXPR_NDARRAY, not a List node), the destructuring branch's
 * `rhs is a List` test failed, and the assignment silently bound NOTHING -- xc
 * and yc stayed symbolic. apply_assignment now normalises a packed RHS to a List
 * of its top-level slices, keeping the slices packed.
 * --------------------------------------------------------------------------- */

/* The reported case: a list of Range[] rows packs to a rank-2 buffer, and
 * destructuring must still bind each row. */
static void test_packed_rhs_destructure(void) {
    const char* syms[] = {"xc", "yc"};
    clear_symbols(syms, 2);

    /* Sanity: the RHS really does pack (else the test proves nothing). */
    assert_eval_eq("NDArrayQ[{Range[4], Range[4]}]", "True", 0);

    eval_and_discard("{xc, yc} = {Range[4], Range[4]}");
    assert_eval_eq("xc", "{1, 2, 3, 4}", 0);
    assert_eval_eq("yc", "{1, 2, 3, 4}", 0);
    assert_eval_eq("Head[xc]", "List", 0);
}

/* A rank-1 packed RHS (Range[4] packs on its own) destructures element-wise,
 * and each bound element must stay an EXACT Integer, not a coerced Real. */
static void test_packed_rank1_rhs_exact(void) {
    const char* syms[] = {"a", "b", "c", "d"};
    clear_symbols(syms, 4);

    assert_eval_eq("NDArrayQ[Range[4]]", "True", 0);

    eval_and_discard("{a, b, c, d} = Range[4]");
    assert_eval_eq("{a, b, c, d}", "{1, 2, 3, 4}", 0);
    assert_eval_eq("Head[a]", "Integer", 0);
    assert_eval_eq("a === 1", "True", 0);
}

/* Packing must SURVIVE the destructuring: a symbol bound to a whole row gets a
 * packed vector, not a materialised List of one Expr per element. */
static void test_packed_rhs_rows_stay_packed(void) {
    const char* syms[] = {"p", "q"};
    clear_symbols(syms, 2);

    eval_and_discard("{p, q} = {Range[10], Range[10]}");
    assert_eval_eq("NDArrayQ[p]", "True", 0);
    assert_eval_eq("NDArrayQ[q]", "True", 0);
    assert_eval_eq("p === Range[10]", "True", 0);
    assert_eval_eq("Length[p]", "10", 0);
}

/* A float64 packed RHS keeps Real leaves. */
static void test_packed_float_rhs(void) {
    const char* syms[] = {"fa", "fb"};
    clear_symbols(syms, 2);

    eval_and_discard("{fa, fb} = {Range[4.0], Range[4.0]}");
    assert_eval_eq("fa === Range[4.0]", "True", 0);
    assert_eval_eq("Head[fa[[1]]]", "Real", 0);
}

/* Fully-nested destructuring against a packed rank-2 buffer reaches scalar
 * leaves, exact integers throughout. */
static void test_packed_nested_destructure(void) {
    const char* syms[] = {"w", "x", "y", "z"};
    clear_symbols(syms, 4);

    eval_and_discard("mat = Table[10 i + j, {i, 2}, {j, 2}]");
    assert_eval_eq("NDArrayQ[mat]", "True", 0);

    eval_and_discard("{{w, x}, {y, z}} = mat");
    assert_eval_eq("{w, x, y, z}", "{11, 12, 21, 22}", 0);
    assert_eval_eq("Head[w]", "Integer", 0);
}

/* Mixed LHS: a bare symbol binds a whole (packed) row while its sibling nests
 * down to scalars. */
static void test_packed_mixed_lhs(void) {
    const char* syms[] = {"a", "b", "c"};
    clear_symbols(syms, 3);

    eval_and_discard("mat = Table[10 i + j, {i, 2}, {j, 2}]");
    eval_and_discard("{a, {b, c}} = mat");
    assert_eval_eq("a", "{11, 12}", 0);
    assert_eval_eq("NDArrayQ[a]", "True", 0);
    assert_eval_eq("{b, c}", "{21, 22}", 0);
}

/* Length mismatch against a packed RHS must leave the targets untouched (no
 * partial assignment), exactly as for a plain-List RHS. */
static void test_packed_length_mismatch(void) {
    const char* syms[] = {"a", "b"};
    clear_symbols(syms, 2);

    eval_and_discard("a = 7");
    eval_and_discard("b = 8");
    eval_and_discard("{a, b} = Range[4]");   /* 2 targets, 4 elements */
    assert_eval_eq("{a, b}", "{7, 8}", 0);
}

/* A literal on the LHS fails the pre-flight even with a packed RHS, so nothing
 * is bound (no partial assignment past the failure). */
static void test_packed_literal_lhs(void) {
    const char* syms[] = {"a", "c", "d"};
    clear_symbols(syms, 3);

    eval_and_discard("{a, 7, c, d} = Range[4]");
    assert_eval_eq("a", "a", 0);
    assert_eval_eq("c", "c", 0);
}

/* A nested-List LHS deeper than the packed array's rank is not an assignable
 * shape and must be rejected cleanly (targets untouched). */
static void test_packed_lhs_deeper_than_rank(void) {
    const char* syms[] = {"a", "b", "c", "d"};
    clear_symbols(syms, 4);

    eval_and_discard("{{a, b}, {c, d}} = Range[4]");   /* rank-1 RHS */
    assert_eval_eq("a", "a", 0);
    assert_eval_eq("d", "d", 0);
}

/* The invariant that mattered all along: destructuring gives the same bound
 * values whether or not the RHS happened to pack. */
static void test_packed_vs_unpacked_equivalence(void) {
    const char* syms[] = {"oa", "ob", "pa", "pb"};
    clear_symbols(syms, 4);

    eval_and_discard("$AutoArrayPacking = False");
    eval_and_discard("{oa, ob} = {Range[4], Range[4]}");
    eval_and_discard("$AutoArrayPacking = True");
    eval_and_discard("{pa, pb} = {Range[4], Range[4]}");
    assert_eval_eq("oa === pa && ob === pb", "True", 0);
}

/* ------------------------------------------------------------------------
 * Scalar-threading over a list LHS (Wolfram Set semantics).
 *
 * {a, b} = c with a non-List c binds a = c and b = c and returns c. Before the
 * fix it returned c without binding anything AND installed a garbage
 * DownValues[List] = {c -> c}, corrupting the List head globally.
 * --------------------------------------------------------------------------- */

/* {a, b} = scalar threads the scalar to both, returns the scalar. */
static void test_thread_scalar(void) {
    const char* syms[] = {"a", "b"};
    clear_symbols(syms, 2);

    assert_eval_eq("{a, b} = 5", "5", 0);
    assert_eval_eq("{a, b}", "{5, 5}", 0);
}

/* A non-atomic, non-List RHS threads whole (each target gets the expression). */
static void test_thread_nonatomic(void) {
    const char* syms[] = {"c", "d"};
    clear_symbols(syms, 2);

    eval_and_discard("{c, d} = g[1, 2]");
    assert_eval_eq("c", "g[1, 2]", 0);
    assert_eval_eq("d", "g[1, 2]", 0);
}

/* Threading recurses into a nested-List target: {p, {q, s}, t} = 9. */
static void test_thread_nested(void) {
    const char* syms[] = {"p", "q", "s", "t"};
    clear_symbols(syms, 4);

    eval_and_discard("{p, {q, s}, t} = 9");
    assert_eval_eq("{p, q, s, t}", "{9, 9, 9, 9}", 0);
}

/* The regression that mattered: threading must NOT install a DownValue on the
 * List head. */
static void test_thread_no_list_corruption(void) {
    const char* syms[] = {"a", "b"};
    clear_symbols(syms, 2);
    symtab_clear_symbol("List");   /* start from a clean List */

    eval_and_discard("{a, b} = 5");
    assert_eval_eq("DownValues[List]", "{}", 0);
    /* And an ordinary List still evaluates to itself, not to the threaded RHS. */
    assert_eval_eq("{10, 20, 30}", "{10, 20, 30}", 0);
}

/* A literal on the LHS fails the pre-flight even in the threaded case, so
 * nothing is bound. */
static void test_thread_literal_lhs(void) {
    const char* syms[] = {"a"};
    clear_symbols(syms, 1);

    eval_and_discard("{a, 7} = 5");
    assert_eval_eq("a", "a", 0);
}

int main(void) {
    symtab_init();
    core_init();

    /* The packed-RHS tests need automatic packing on regardless of the
     * MATHILDA_NO_PACK env var a CI run might set. */
    eval_and_discard("$AutoArrayPacking = True");

    TEST(test_destructure_fresh);
    TEST(test_reassignment_updates_targets);
    TEST(test_reassign_then_increment_decrement);
    TEST(test_simultaneous_swap);
    TEST(test_nested_destructuring);
    TEST(test_downvalue_destructuring);
    TEST(test_pattern_downvalue_destructuring);
    TEST(test_length_mismatch);
    TEST(test_literal_lhs_element_does_not_pretend_success);
    TEST(test_repeated_reassignment_smoke);

    /* Packed-array RHS destructuring (automatic-packing regression). */
    TEST(test_packed_rhs_destructure);
    TEST(test_packed_rank1_rhs_exact);
    TEST(test_packed_rhs_rows_stay_packed);
    TEST(test_packed_float_rhs);
    TEST(test_packed_nested_destructure);
    TEST(test_packed_mixed_lhs);
    TEST(test_packed_length_mismatch);
    TEST(test_packed_literal_lhs);
    TEST(test_packed_lhs_deeper_than_rank);
    TEST(test_packed_vs_unpacked_equivalence);

    /* Scalar threading over a list LHS (Wolfram Set semantics). */
    TEST(test_thread_scalar);
    TEST(test_thread_nonatomic);
    TEST(test_thread_nested);
    TEST(test_thread_no_list_corruption);
    TEST(test_thread_literal_lhs);

    printf("All tests passed!\n");
    return 0;
}
