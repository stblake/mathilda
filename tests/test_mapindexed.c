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

/* MapIndexed[f, expr] / MapIndexed[f, expr, levelspec] — f applied to the
 * selected parts, each paired with its position.
 *
 * Every expected value below is Wolfram's own documented output for the same
 * input, so this file doubles as a conformance record. Two properties do most
 * of the work and are asserted directly rather than left implicit:
 *
 *   - the position is the one Part/Extract take, so Extract[expr, #2] === #1
 *     at every visited node (test_roundtrip_*);
 *   - MapIndexed "constructs a complete new expression and then evaluates it",
 *     which is observable through a held head (test_hold_defers_evaluation) —
 *     Map, which evaluates f per node, gets that case wrong.
 *
 * Bare machine reals are never asserted on directly (their printed form is
 * formatting-sensitive); NDArray cases go through Head, a Part equality, or a
 * position, following tests/test_map_ndarray.c. */

/* Evaluate `input` and compare its infix printed form to `expected`. */
static void run(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "MapIndexed %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* Same, against FullForm — used where the surface syntax would hide the tree
 * (Rule vs RuleDelayed, Rational, Association). */
static void run_full(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string_fullform(r);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "MapIndexed %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* ---------- Default level {1} ---------- */

static void test_default_level_list(void) {
    run("MapIndexed[f, {a, b, c, d}]",
        "{f[a, {1}], f[b, {2}], f[c, {3}], f[d, {4}]}");
    run("MapIndexed[f, {10, 20, 30}]",
        "{f[10, {1}], f[20, {2}], f[30, {3}]}");
}

static void test_default_level_is_shallow(void) {
    /* Only level 1: the nested lists are handed to f whole. */
    run("MapIndexed[f, {{a, b}, {c, d, e}}]",
        "{f[{a, b}, {1}], f[{c, d, e}, {2}]}");
    run("MapIndexed[f, {{{{{a}}}}}]", "{f[{{{{a}}}}, {1}]}");
}

static void test_empty_list(void) {
    run("MapIndexed[f, {}]", "{}");
    run("MapIndexed[f, g[]]", "g[]");
}

static void test_non_list_head(void) {
    /* The head is preserved; only the arguments are mapped. */
    run("MapIndexed[f, g[a, b, c]]", "g[f[a, {1}], f[b, {2}], f[c, {3}]]");
    /* A compound head is left alone without Heads->True. */
    run("MapIndexed[f, p[x][a, b]]", "p[x][f[a, {1}], f[b, {2}]]");
}

static void test_expression_is_itself_a_rule(void) {
    /* a -> b in the *expression* slot is mapped, not mistaken for an option:
     * the option guard only ever inspects argument 3. */
    run("MapIndexed[f, a -> b]", "f[a, {1}] -> f[b, {2}]");
}

static void test_orderless_head(void) {
    /* Plus re-canonicalises on re-evaluation; the elements are already in
     * canonical order here, so the result is stable. */
    run("MapIndexed[f, {a + b, c}]", "{f[a + b, {1}], f[c, {2}]}");
}

static void test_hold_defers_evaluation(void) {
    /* MapIndexed builds the whole new expression and lets the evaluator reduce
     * it, so Hold suppresses the arithmetic inside. Evaluating f per node (as
     * Map does) would give Hold[f[2, {1}], f[4, {2}]]. */
    run("MapIndexed[f, Hold[1 + 1, 2 + 2]]",
        "Hold[f[1 + 1, {1}], f[2 + 2, {2}]]");
}

/* ---------- Slot access ---------- */

static void test_slots(void) {
    run("MapIndexed[#1 &, {a, b, c}]", "{a, b, c}");
    run("MapIndexed[#2 &, {a, b, c}]", "{{1}, {2}, {3}}");
    run("MapIndexed[First[#2] + f[#1] &, {a, b, c, d}]",
        "{1 + f[a], 2 + f[b], 3 + f[c], 4 + f[d]}");
}

/* ---------- Positive level specs ---------- */

static void test_level_n_maps_one_through_n(void) {
    /* n means levels 1..n, and the traversal is bottom-up: the level-1 node f
     * sees has already had its own parts wrapped. */
    run("MapIndexed[f, {{{{{a}}}}}, 2]", "{f[{f[{{{a}}}, {1, 1}]}, {1}]}");
    run("MapIndexed[f, {{{{{a}}}}}, 3]",
        "{f[{f[{f[{{a}}, {1, 1, 1}]}, {1, 1}]}, {1}]}");
}

static void test_level_brace_n_is_that_level_only(void) {
    run("MapIndexed[f, {{{{{a}}}}}, {2}]", "{{f[{{{a}}}, {1, 1}]}}");
    run("MapIndexed[f, {{a, b}, {c, d, e}}, {2}]",
        "{{f[a, {1, 1}], f[b, {1, 2}]}, "
        "{f[c, {2, 1}], f[d, {2, 2}], f[e, {2, 3}]}}");
}

static void test_level_range(void) {
    run("MapIndexed[f, {{a}, b}, {1, 2}]",
        "{f[{f[a, {1, 1}]}, {1}], f[b, {2}]}");
    run("MapIndexed[f, {{{a}}}, {2, 3}]",
        "{{f[{f[a, {1, 1, 1}]}, {1, 1}]}}");
}

static void test_level_zero_has_empty_position(void) {
    /* Level 0 is the whole expression; its position is {}. */
    run("MapIndexed[f, {a, b}, {0, 1}]", "f[{f[a, {1}], f[b, {2}]}, {}]");
    run("MapIndexed[f, {a, b}, {0}]", "f[{a, b}, {}]");
    run("MapIndexed[f, a, {0}]", "f[a, {}]");
    run("MapIndexed[f, g[], {0}]", "f[g[], {}]");
}

static void test_level_infinity(void) {
    run("MapIndexed[f, {{a, b}, {c, d, {e}}}, Infinity]",
        "{f[{f[a, {1, 1}], f[b, {1, 2}]}, {1}], "
        "f[{f[c, {2, 1}], f[d, {2, 2}], f[{f[e, {2, 3, 1}]}, {2, 3}]}, {2}]}");
    /* Infinity in the second slot of a two-element spec, which the lenient
     * parser Map still uses cannot express. */
    run("MapIndexed[f, {a, b}, {0, Infinity}]",
        "f[{f[a, {1}], f[b, {2}]}, {}]");
    run("MapIndexed[f, {a, b}, {2, Infinity}]", "{a, b}");
}

/* ---------- Negative and mixed level specs ---------- */

static void test_negative_level_brace(void) {
    /* {-1}: the leaves only. */
    run("MapIndexed[f, {{a, b}, {c, d, {e}}}, {-1}]",
        "{{f[a, {1, 1}], f[b, {1, 2}]}, "
        "{f[c, {2, 1}], f[d, {2, 2}], {f[e, {2, 3, 1}]}}}");
    /* {-2}: exactly the parts of depth 2. */
    run("MapIndexed[f, {{{{a}}}}, {-2}]", "{{{f[{a}, {1, 1, 1}]}}}");
}

static void test_negative_level_n(void) {
    /* -n means levels 1 through -n, i.e. everything of depth >= n. */
    run("MapIndexed[f, {{{{a}}}}, -1]",
        "{f[{f[{f[{f[a, {1, 1, 1, 1}]}, {1, 1, 1}]}, {1, 1}]}, {1}]}");
    run("MapIndexed[f, {{{{a}}}}, -2]",
        "{f[{f[{f[{a}, {1, 1, 1}]}, {1, 1}]}, {1}]}");
    run("MapIndexed[f, {{{{a}}}}, -3]", "{f[{f[{{a}}, {1, 1}]}, {1}]}");
}

static void test_negative_range(void) {
    /* {-3, -2}: depth 2 through 3. */
    run("MapIndexed[f, {{{{a}}}}, {-3, -2}]",
        "{{f[{f[{a}, {1, 1, 1}]}, {1, 1}]}}");
    /* {-1, -3} is empty: nothing has depth <= 1 and >= 3. */
    run("MapIndexed[f, {{{{a}}}}, {-1, -3}]", "{{{{a}}}}");
}

static void test_mixed_bounds(void) {
    /* {2, -3}: at level 2 or deeper AND of depth 3 or more. Different heads at
     * each level make the selected nodes unambiguous. */
    run("MapIndexed[f, h0[h1[h2[h3[h4[a]]]]], {2, -3}]",
        "h0[h1[f[h2[f[h3[h4[a]], {1, 1, 1}]], {1, 1}]]]");
    /* A negative lower bound with a positive upper one. */
    run("MapIndexed[f, {a, b}, {-3, 2}]", "f[{f[a, {1}], f[b, {2}]}, {}]");
}

/* ---------- Empty level ranges are the identity ---------- */

static void test_empty_range_is_identity(void) {
    run("MapIndexed[f, {a, b}, 0]", "{a, b}");
    run("MapIndexed[f, {a, b}, {3, 1}]", "{a, b}");
    run("MapIndexed[f, {a, b}, {5}]", "{a, b}");
}

/* ---------- Heads -> True ---------- */

static void test_heads_true(void) {
    /* A head is part 0 of its node. */
    run("MapIndexed[f, {a, b, c}, Heads -> True]",
        "f[List, {0}][f[a, {1}], f[b, {2}], f[c, {3}]]");
    run("MapIndexed[f, g[], Heads -> True]", "f[g, {0}][]");
}

static void test_heads_true_all_levels(void) {
    run("MapIndexed[f, p[x][a, b, c], Infinity, Heads -> True]",
        "f[f[p, {0, 0}][f[x, {0, 1}]], {0}][f[a, {1}], f[b, {2}], f[c, {3}]]");
}

static void test_heads_option_forms(void) {
    /* Explicit False is the default. */
    run("MapIndexed[f, {a, b, c}, Heads -> False]",
        "{f[a, {1}], f[b, {2}], f[c, {3}]}");
    /* The option may follow an explicit level spec, or take its slot. */
    run("MapIndexed[f, {a, b}, {1}, Heads -> True]",
        "f[List, {0}][f[a, {1}], f[b, {2}]]");
    run("MapIndexed[f, {a, b}, Heads -> True]",
        "f[List, {0}][f[a, {1}], f[b, {2}]]");
    /* A delayed rule is a legal option too. */
    run("MapIndexed[f, {a, b}, Heads :> True]",
        "f[List, {0}][f[a, {1}], f[b, {2}]]");
}

static void test_heads_true_negative_level(void) {
    /* With Heads->True the head counts towards depth: Depth[p[q[r]][a]] is 2,
     * but Depth[p[q[r]][a], Heads->True] is 4. {-2} therefore selects q[r]
     * (depth 2) at position {0, 1}. A heads-blind depth would select the root
     * instead. */
    run("Depth[p[q[r]][a], Heads -> True]", "4");
    run("MapIndexed[f, p[q[r]][a], {-2}, Heads -> True]",
        "p[f[q[r], {0, 1}]][a]");
}

/* ---------- Atoms ---------- */

static void test_atoms_pass_through(void) {
    /* An atom has no level-1 parts, so the default spec leaves it alone —
     * matching Map[f, a] -> a. */
    run("MapIndexed[f, a]", "a");
    run("MapIndexed[f, 5]", "5");
    run("MapIndexed[f, \"abc\"]", "\"abc\"");
}

static void test_rational_and_complex_are_atomic(void) {
    /* Rational and Complex are atomic here, as they are for Depth and Level. */
    run_full("MapIndexed[f, 1/2]", "Rational[1, 2]");
    run_full("MapIndexed[f, 2 + 3 I]", "Complex[2, 3]");
    run_full("MapIndexed[f, 1/2, Infinity]", "Rational[1, 2]");
}

/* ---------- Associations ---------- */

static void test_association_default_level(void) {
    /* The parts of an association are its values, positioned by Key[k]; the
     * keys are preserved. Key types are not restricted. */
    run("MapIndexed[f, <|\"a\" -> 1, a -> 2, 1 -> 1|>]",
        "<|\"a\" -> f[1, {Key[\"a\"]}], a -> f[2, {Key[a]}], "
        "1 -> f[1, {Key[1]}]|>");
    run("MapIndexed[f, <|\"a\" -> 10, \"b\" -> 20|>]",
        "<|\"a\" -> f[10, {Key[\"a\"]}], \"b\" -> f[20, {Key[\"b\"]}]|>");
    run("MapIndexed[#1 &, <|\"a\" -> 10, \"b\" -> 20|>]",
        "<|\"a\" -> 10, \"b\" -> 20|>");
    run("MapIndexed[#2 &, <|\"a\" -> 10, \"b\" -> 20|>]",
        "<|\"a\" -> {Key[\"a\"]}, \"b\" -> {Key[\"b\"]}|>");
}

static void test_association_empty_and_single(void) {
    run("MapIndexed[f, <||>]", "<||>");
    run("MapIndexed[f, <|a -> 1|>]", "<|a -> f[1, {Key[a]}]|>");
}

static void test_association_nested(void) {
    /* A Rule wrapper is not a level of its own, so a value two associations
     * deep is at level 2 with a two-component position. */
    run("MapIndexed[h, <|a -> <|b -> c, p -> <|q -> r|>|>, d -> {e}|>, {2}]",
        "<|a -> <|b -> h[c, {Key[a], Key[b]}], "
        "p -> h[<|q -> r|>, {Key[a], Key[p]}]|>, "
        "d -> {h[e, {Key[d], 1}]}|>");
    /* An association value that is a list mixes Key and integer components. */
    run("MapIndexed[f, <|a -> {b, c}|>, {2}]",
        "<|a -> {f[b, {Key[a], 1}], f[c, {Key[a], 2}]}|>");
}

static void test_association_head_is_never_mapped(void) {
    run("MapIndexed[f, <|a -> 1|>, Heads -> True]",
        "<|a -> f[1, {Key[a]}]|>");
}

static void test_association_ruledelayed_preserved(void) {
    /* <|a :> 1|> survives evaluation as Association[RuleDelayed[a, 1]] and
     * must stay delayed after mapping. */
    run_full("MapIndexed[f, <|a :> 1|>]",
             "Association[RuleDelayed[a, f[1, List[Key[a]]]]]");
}

static void test_malformed_association_does_not_crash(void) {
    /* is_association() only checks the head, and a malformed Association[a, b]
     * survives evaluation, so its entries are not necessarily rules. Reading
     * args[1] of the symbol `a` used to type-pun a SymbolDef* and segfault. */
    run_full("MapIndexed[f, Association[a, b]]",
             "Association[f[a, List[1]], f[b, List[2]]]");
    /* The same bug was reachable through Map's assoc_map_values. */
    run_full("Map[f, Association[a, b]]", "Association[a, b]");
}

/* ---------- NDArray ---------- */

static void test_ndarray_default_level(void) {
    /* Rank 1: each leading-axis part is a scalar, positioned {i}. */
    run("MapIndexed[#2 &, NDArray[Range[3]]]", "{{1}, {2}, {3}}");
    run("MapIndexed[Head[#1] &, NDArray[Range[3]]]", "{Real, Real, Real}");
    /* Rank 2: each part is a sub-NDArray row. */
    run("MapIndexed[#2 &, NDArray[Table[i + j, {i, 2}, {j, 2}]]]",
        "{{1}, {2}}");
    run("MapIndexed[Head[#1] &, NDArray[Table[i + j, {i, 2}, {j, 2}]]]",
        "{NDArray, NDArray}");
}

static void test_ndarray_level_spec(void) {
    /* A non-default spec materializes the array and maps generically, so the
     * scalars of a rank-2 array become reachable at level 2. */
    run("MapIndexed[#2 &, NDArray[Table[i + j, {i, 2}, {j, 2}]], {2}]",
        "{{{1, 1}, {1, 2}}, {{2, 1}, {2, 2}}}");
    /* An empty range keeps the array packed rather than unpacking it. */
    run("Head[MapIndexed[f, NDArray[Range[3]], 0]]", "NDArray");
}

static void test_ndarray_nested_is_atomic(void) {
    /* An NDArray inside a list is an atom: the traversal does not descend. */
    run("MapIndexed[#2 &, {NDArray[Range[2]], x}, Infinity]", "{{1}, {2}}");
}

/* ---------- Arity and malformed arguments ---------- */

static void test_arity_zero_reports_argb(void) {
    /* MapIndexed::argb goes to stderr and the call is left unevaluated. */
    run("MapIndexed[]", "MapIndexed[]");
}

static void test_one_argument_is_left_alone(void) {
    run("MapIndexed[f]", "MapIndexed[f]");
}

static void test_bad_level_spec_left_unevaluated(void) {
    run("MapIndexed[f, x, foo]", "MapIndexed[f, x, foo]");
    run("MapIndexed[f, {a, b}, 1.5]", "MapIndexed[f, {a, b}, 1.5]");
    run("MapIndexed[f, {a, b}, {1, 2, 3}]", "MapIndexed[f, {a, b}, {1, 2, 3}]");
}

static void test_non_option_trailing_argument(void) {
    /* MapIndexed::nonopt, then unevaluated — never silently ignored. */
    run("MapIndexed[f, x, 2, 3]", "MapIndexed[f, x, 2, 3]");
}

static void test_unknown_option(void) {
    /* MapIndexed::optx, then unevaluated. */
    run("MapIndexed[f, {a, b}, x -> y]", "MapIndexed[f, {a, b}, x -> y]");
}

/* ---------- Options[] ---------- */

static void test_options(void) {
    run("Options[MapIndexed]", "{Heads -> False}");
}

/* ---------- Round trip: Extract[expr, #2] === #1 ---------- */

static void test_roundtrip_list(void) {
    run("MapIndexed[Extract[{{a, b}, {c, d, {e}}}, #2] === #1 &, "
        "{{a, b}, {c, d, {e}}}, {-1}]",
        "{{True, True}, {True, True, {True}}}");
}

static void test_roundtrip_heads(void) {
    /* Head positions carry a 0 component, which Extract also understands. */
    run("MapIndexed[Extract[p[x][a, b], #2] === #1 &, p[x][a, b], {-1}, "
        "Heads -> True]",
        "True[True][True, True]");
}

static void test_roundtrip_association(void) {
    run("MapIndexed[Extract[<|\"a\" -> 1, b -> <|c -> 2|>|>, #2] === #1 &, "
        "<|\"a\" -> 1, b -> <|c -> 2|>|>, {-1}]",
        "<|\"a\" -> True, b -> <|c -> True|>|>");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_default_level_list);
    TEST(test_default_level_is_shallow);
    TEST(test_empty_list);
    TEST(test_non_list_head);
    TEST(test_expression_is_itself_a_rule);
    TEST(test_orderless_head);
    TEST(test_hold_defers_evaluation);

    TEST(test_slots);

    TEST(test_level_n_maps_one_through_n);
    TEST(test_level_brace_n_is_that_level_only);
    TEST(test_level_range);
    TEST(test_level_zero_has_empty_position);
    TEST(test_level_infinity);

    TEST(test_negative_level_brace);
    TEST(test_negative_level_n);
    TEST(test_negative_range);
    TEST(test_mixed_bounds);
    TEST(test_empty_range_is_identity);

    TEST(test_heads_true);
    TEST(test_heads_true_all_levels);
    TEST(test_heads_option_forms);
    TEST(test_heads_true_negative_level);

    TEST(test_atoms_pass_through);
    TEST(test_rational_and_complex_are_atomic);

    TEST(test_association_default_level);
    TEST(test_association_empty_and_single);
    TEST(test_association_nested);
    TEST(test_association_head_is_never_mapped);
    TEST(test_association_ruledelayed_preserved);
    TEST(test_malformed_association_does_not_crash);

    TEST(test_ndarray_default_level);
    TEST(test_ndarray_level_spec);
    TEST(test_ndarray_nested_is_atomic);

    TEST(test_arity_zero_reports_argb);
    TEST(test_one_argument_is_left_alone);
    TEST(test_bad_level_spec_left_unevaluated);
    TEST(test_non_option_trailing_argument);
    TEST(test_unknown_option);

    TEST(test_options);

    TEST(test_roundtrip_list);
    TEST(test_roundtrip_heads);
    TEST(test_roundtrip_association);

    printf("All MapIndexed tests passed!\n");
    return 0;
}
