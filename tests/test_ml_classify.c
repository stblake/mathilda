/* test_ml_classify.c -- Classify, ClassifierFunction, and the label vocabulary.
 *
 * Two of these assertions are ABSOLUTE rather than comparative, which is what makes them
 * worth more than an accuracy figure:
 *
 *   - a k = 1 nearest-neighbour classifier must reproduce EVERY training label exactly,
 *     because each training point is its own nearest neighbour. Not "high accuracy" --
 *     exactly right, every row;
 *   - class probabilities must sum to 1. A classifier whose probabilities were each
 *     individually plausible but did not sum to 1 would pass any per-class check.
 *
 * The rest establish that a class really is an arbitrary expression -- strings, symbols
 * and integers all tested -- which is the point of the vocabulary this family added.
 */
#include <stdio.h>
#include "test_utils.h"

extern void symtab_init(void);
extern void core_init(void);

#define TR "{{0., 0.} -> \"red\", {1., 0.} -> \"red\", {0., 1.} -> \"red\", " \
           "{10., 10.} -> \"blue\", {11., 10.} -> \"blue\", {10., 11.} -> \"blue\"}"

static void test_k1_reproduces_every_training_label(void) {
    /* The exact, absolute check. Each training point is its own nearest neighbour, so a
     * k = 1 classifier that gets any row wrong has a bug in its distance, its vocabulary
     * or its vote -- there is no "close enough" here. */
    assert_eval_eq("Module[{tr = " TR ", c},"
                   " c = Classify[tr];"
                   " And @@ Table[c[First[tr[[i]]]] === Last[tr[[i]]], {i, Length[tr]}]]",
                   "True", 0);
    /* And it generalises the obvious way on either side of the gap. */
    assert_eval_eq("Classify[" TR "][{0.5, 0.5}]", "\"red\"", 0);
    assert_eval_eq("Classify[" TR "][{10.5, 10.5}]", "\"blue\"", 0);
}

static void test_class_probabilities_sum_to_one(void) {
    /* The second absolute check. Asserted at k = 5 as well, where the vote is genuinely
     * split (0.6/0.4) rather than degenerate -- a 1.0/0.0 split would sum to 1 even if the
     * normalisation were wrong. */
    assert_eval_eq("Chop[Total[Map[Last, Classify[" TR "][{0.5, 0.5}, \"Probabilities\"]]] - 1.]",
                   "0", 0);
    assert_eval_eq("Chop[Total[Map[Last, Classify[" TR ", "
                   "Method -> {\"NearestNeighbors\", \"NeighborsNumber\" -> 5}]"
                   "[{5., 5.}, \"Probabilities\"]]] - 1.]", "0", 0);
    /* The split at k = 5 must be non-degenerate, or the row above proves little. */
    assert_eval_eq("Module[{p = Classify[" TR ", Method -> {\"NearestNeighbors\", "
                   "\"NeighborsNumber\" -> 5}][{5., 5.}, \"Probabilities\"]},"
                   " 0. < Last[First[p]] < 1.]", "True", 0);
}

static void test_a_class_can_be_any_expression(void) {
    /* The whole point of the label vocabulary: earlier families took numeric responses,
     * and a class is not a number. Strings, symbols and integers must all work, and a
     * symbol must not be confused with the string of the same name. */
    assert_eval_eq("Classify[{1. -> a, 2. -> a, 9. -> b, 10. -> b}][\"Classes\"]",
                   "{a, b}", 0);
    assert_eval_eq("Classify[{1. -> a, 2. -> a, 9. -> b, 10. -> b}][1.5]", "a", 0);
    assert_eval_eq("Classify[{1. -> 0, 2. -> 0, 9. -> 1, 10. -> 1}][\"Classes\"]",
                   "{0, 1}", 0);
    assert_eval_eq("Classify[{1. -> 0, 2. -> 0, 9. -> 1, 10. -> 1}][9.5]", "1", 0);
    /* Comparison is structural, so "a" and a are two classes rather than one. */
    assert_eval_eq("Length[Classify[{1. -> \"a\", 2. -> a, 9. -> \"b\"}][\"Classes\"]]",
                   "3", 0);
}

static void test_class_order_is_first_appearance_and_deterministic(void) {
    /* Not sorted -- first appearance in the training data. What matters is that it is
     * deterministic and readable next to the data; reversing the data reverses the
     * vocabulary, which demonstrates the rule rather than an accident. */
    assert_eval_eq("Classify[" TR "][\"Classes\"]", "{\"red\", \"blue\"}", 0);
    assert_eval_eq("Classify[Reverse[" TR "]][\"Classes\"]", "{\"blue\", \"red\"}", 0);
    /* The prediction must not depend on that ordering. */
    assert_eval_eq("Classify[" TR "][{0.5, 0.5}] === "
                   "Classify[Reverse[" TR "]][{0.5, 0.5}]", "True", 0);
}

static void test_three_classes_and_properties(void) {
    assert_eval_eq("Classify[{{0.,0.} -> \"a\", {5.,5.} -> \"b\", "
                   "{10.,10.} -> \"c\"}][{4.9, 4.9}]", "\"b\"", 0);
    assert_eval_eq("Classify[" TR "][\"FeatureCount\"]", "2", 0);
    assert_eval_eq("Classify[" TR "][\"NeighborCount\"]", "1", 0);
    assert_eval_eq("Classify[" TR "][\"Method\"]", "\"NearestNeighbors\"", 0);
    /* A fitted classifier elides, like every other fitted model. */
    assert_eval_eq("Classify[" TR "]", "ClassifierFunction[\"NearestNeighbors\", <>]", 0);
}

static void test_malformed_input_declines(void) {
    /* Only the rule form is accepted, and that is not arbitrary: a matrix with the class
     * in its last column cannot work, because a class need not be a number. */
    assert_eval_eq("Head[Classify[{{1., 2.}, {3., 4.}}]]", "Classify", 0);
    assert_eval_eq("Head[Classify[{}]]", "Classify", 0);
    assert_eval_eq("Head[Classify[x]]", "Classify", 0);
    assert_eval_eq("Head[Classify[" TR ", Method -> \"RandomForest\"]]", "Classify", 0);
    assert_eval_eq("Head[Classify[" TR ", Method -> {\"NearestNeighbors\", "
                   "\"NeighborsNumber\" -> 0}]]", "Classify", 0);
    /* A wrongly-shaped query, and an unknown property, leave the application
     * unevaluated -- so no class comes back. */
    assert_eval_eq("StringQ[Classify[" TR "][{1.}]]", "False", 0);
    assert_eval_eq("StringQ[Classify[" TR "][\"NoSuchProperty\"]]", "False", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_k1_reproduces_every_training_label);
    TEST(test_class_probabilities_sum_to_one);
    TEST(test_a_class_can_be_any_expression);
    TEST(test_class_order_is_first_appearance_and_deterministic);
    TEST(test_three_classes_and_properties);
    TEST(test_malformed_input_declines);

    printf("All ml Classify tests passed.\n");
    return 0;
}
